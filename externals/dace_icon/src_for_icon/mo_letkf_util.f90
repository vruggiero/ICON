!
!+ Low level utility functions for the LETKF
!
MODULE mo_letkf_util
!
! Description:
!   Low level utility functions for the LETKF
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_28        2014/02/26 Andreas Rhodin
!  new module with low level utility functions for the LETKF
! V1_31        2014-08-21 Andreas Rhodin
!  move time-correlation random pattern generator from mo_2d_random here
! V1_37        2014-12-23 Andreas Rhodin
!  define 'gaspari_cohn' as elemental function
! V1_42        2015-06-08 Andreas Rhodin
!  fix for VarEnKF (vertical localisation)
! V1_43        2015-08-19 Andreas Rhodin
!  avoid double implementation of apply_H_member0/1
! V1_45        2015-12-15 Harald Anlauf
!  gaspari_cohn: fix logic, rewrite using Horner scheme
! V1_46        2016-02-05 Valerio Cardinali
!  height dependent horizontal localisation (COMET, italian met. center)
! V1_48        2016-10-06 Andreas Rhodin
!  gaspari_cohn: return result >=0 in any case (fix rounding errors)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2014
!==============================================================================

  !=============
  ! Modules used
  !=============
  use mo_kind,       only: wp             ! working precision kind
  use mo_exception,  only: finish         ! abort in case of error
  use mo_atm_grid,   only: t_grid,       &! grid meta data type
                           construct,    &! grid constructor routine
                           print,        &! print grid meta data
                           phi2phirot,   &! latitude  -> rotated system
                           rla2rlarot     ! longitude -> rotated system
  use mo_atm_state,  only: t_atm          ! atmospheric state derived type
  use mo_mpi_dace,   only: dace,         &! MPI group info
                           p_min, p_max, &! generic MPI routines
                           p_sum          !
  use mo_time,       only: t_time,       &! date & time derived type
                           hours          ! derive hours from t_time
  use mo_obs_set,    only: t_obs_set      ! observation data type
! use mo_algorithms, only: init_splinex, &! set coefficients (2nd deriv.) for spline interpol.
!                          splint         ! spline interpol. with optimization for subseq.calls
  use mo_fdbk_tables,only: OT_RAD         ! Radiances          report id.
! use mo_physics,    only: r2d            ! 180 / pi
! use mo_wmo_tables, only: WMO6_ROTLL,      &! rotated lat/lon grid  type
!                          WMO6_LATLON,     &! lat/lon         grid  type
!                          DWD6_ICON,       &! ICON            grid  type
!                          DWD6_ICOSAHEDRON  ! GME             grid  type

  implicit none

  !================
  ! Public entities
  !================
  private
  public :: set_lv           ! set vertical localisation length scale
  public :: set_lh           ! set horizontal localisation length scale
  public :: set_w_f          ! set coarse weighting function vectors
  public :: gaspari_cohn     ! Gaspari&Cohn function
  public :: t_rndm_tcor      ! temporal correlation meta data
  public :: construct        ! set up temporal correlation meta data
  public :: flat_copy        ! flat copy of atmospheric state

  !-----------------
  ! Type definitions
  !-----------------
  type t_rndm_tcor
    real(wp) :: scale ! requested correlation time scale
    integer  :: n     ! number of time slices required
    real(wp) :: w (5) ! weights for temporal interpolation
    integer  :: ht(5) ! seed for random number generator
  end type t_rndm_tcor

  !-----------
  ! Interfaces
  !-----------
  interface construct
    module procedure construct_rndm_tcor
  end interface construct

!==============================================================================
contains
!==============================================================================
  subroutine set_lv (nzr, lv_surf, lv_top, r_0, flr_d, pr, flr, enkf_mean,&
                     fr, fp, f0, b, d)
  !---------------------------------------------------------------------
  ! set vertical localisation length scale constant or linear in log(p)
  ! and pressure levels on appropriate vertical grid
  ! special values for nzr (number of height levels):
  !   nzr > 1 : number of height levels prescribed
  !   nzr = 1 : no vertical localisation
  !   nzr = 0 : no vertical coarse grid (use native grid)
  !   nzr < 0 : number of height levels to be determined (<= loc.length)
  !---------------------------------------------------------------------
  integer          ,intent(inout) :: nzr       ! number of height levels
  real(wp)         ,intent(in)    :: lv_surf   ! localisation length at surface
  real(wp)         ,intent(in)    :: lv_top    ! localisation length at top
  real(wp)         ,intent(out)   :: r_0       ! log(p) reference
  real(wp)         ,intent(out)   :: flr_d     ! localisation length increment
  real(wp)         ,pointer       :: pr (:)    ! p-levels on reduced grid (Pa)
  real(wp)         ,pointer       :: flr(:)    ! loc.length on reduced grid
  type(t_atm),optional,intent(in) :: enkf_mean ! ensemble mean (for pressure)
  real(wp),optional,pointer       :: fr(:)     ! levels in transformed coord.
  real(wp),optional,intent(out)   :: fp        ! transformation parameter
  real(wp),optional,intent(out)   :: f0        ! transformation parameter
  real(wp),optional,intent(out)   :: b         ! transformation parameter
  real(wp),optional,intent(out)   :: d         ! transformation parameter

  !--------------------------------------------------------------------
  ! Method:
  !
  ! The localisation length scale 'flr' (in log(p)) is required to
  ! vary linearly in r=log(p):
  ! flr = flr_d * (r - r0)
  ! with appropriate constants flr_d and r0 so that:
  ! lv_top  = flr (log(pmin)) = flr_d * (log(pmin) - r0)
  ! lv_surf = flr (log(pmax)) = flr_d * (log(pmax) - r0) .
  !
  ! For a given number of grid-points (nzr) of the 'coarse' grid to be
  ! used in the LETKF calculations optimal grid-points are sought:
  ! Thus we look for a further coordinate transformation f(r) with
  ! dr/df ~ flr or f ~ log(r) . In this coordinate we chose an
  ! equidistant grid spacing f = fb + (i-1) * df for grid point
  ! indices i and appropriate constants fb and df.
  !
  ! Final output variables flr_d and flr are multiplied by
  ! sqrt(10./3.) to be used as parameters in the Gaspari& Cohn
  ! localisation function.
  !
  ! For usage in the VarEnKF a suitable number of grid-points nzr (set
  ! <0 on program entry in this case) is chosen so that the actual
  ! grid spacing in log(p) is equal to or slightly less than the
  ! localisation length scale flr. For this application the parameters
  ! fp and f0 are returned to calculate f(r) = fp * log(r0-r) + f0 , as
  ! well as the values of fr at the chosen gridpoints. These
  ! parameters are scaled by sqrt(10./3.) as well.
  !---------------------------------------------------------------------


    !----------------
    ! local variables
    !----------------
    integer  :: i               ! loop index
    real(wp) :: delta_p         ! delta p in log(p)
    real(wp) :: pmin, pmin_lim  ! min. pressure (fg,limit)
    real(wp) :: pmax, pmax_lim  ! max. pressure (fg,limit)
    real(wp) :: s1,s2,r,fb,ft   ! aux. variables
    real(wp) :: df,rt,rb,f      ! aux. variables
    integer  :: nzg             ! number of levels of native grid
    logical  :: lnative         ! use native grid

    real(wp) ,parameter ::   sq103 = sqrt (10._wp/3._wp)
    real(wp) ,parameter :: edsq103 = 1._wp / sq103

    !----------------------------------
    ! derive pressure bounds (min, max)
    !----------------------------------
    pmax     = 104000._wp     ! typical values for current ICON
    pmin     =      1.34_wp   ! used for verification (ICON state not available)
    pmax_lim = 105000._wp               ! max press. on diagn. levels
    pmin_lim =   4000._wp               ! min press. on diagn. levels
    nzg      = 1
    if (present(enkf_mean)) then
      nzg    =         enkf_mean% ub(3)
      pmax   = maxval (enkf_mean% pf(:,:,nzg,:))
      pmin   = minval (enkf_mean% pf(:,:, 1 ,:))
    endif
    pmax     = p_max  (pmax)
    pmax     = max    (pmax,pmax_lim)   ! check if true pmax is larger
    pmin     = p_min  (pmin)
    pmin     = min    (pmin,pmin_lim)   ! check if true pmin is smaller

    !--------------------------------------------
    ! set number of levels if native grid is used
    !--------------------------------------------
    lnative = (nzr == 0)
    if (lnative) nzr = nzg

    !------------------------------------------------------
    ! for native grid: derive mean pressure of model levels
    !------------------------------------------------------
    if (lnative) then
      allocate (pr (nzr))
      if (.not.present (enkf_mean))                                  &
        call finish ('set_lv','native grid and enkf_mean not present')
      if (.not.associated (enkf_mean% pf))                    &
        call finish ('set_lv','native grid and pf not present')
      pr = 0._wp
      do i = 1, nzr
        pr (i) = sum (enkf_mean% pf(:,:,i,:))
      end do
      pr = log (p_sum (pr) / enkf_mean% grid% nxny)
    endif

    if (nzr/=1) then      ! only if vertical localisation is used
      !---------------------------------------------
      ! derive vertical coordinate (linear in log p)
      !---------------------------------------------
      if(dace% lpio) then
        write(6,*)
        write(6,*) ' set_lv: set coarse grid pressure levels'
        write(6,*)
        write(6,*) '   pmin, pmax, lv_surf, lv_top = ' &
                     , pmin, pmax, lv_surf, lv_top
        write(6,*)
      end if
      !---------------------------------------------------
      ! set # levels (if not prescribed) and log(p) values
      !---------------------------------------------------
      if (lv_surf == lv_top) then
        !------------------------------
        ! for constant spacing in ln(p)
        !------------------------------
        if (lv_surf <= 0._wp) call finish ("set_lv", "lv_surf = lv_top <= 0")
        if (nzr < 1) nzr = ceiling (log (pmax/pmin) / lv_surf) + 1
        delta_p  = (log(pmax) - log(pmin)) / max (nzr-1, 1)
        !------------------------------------------------
        ! allocate arrays for precalculated length scales
        ! and pressure values on coarse grid
        !------------------------------------------------
        allocate (flr(nzr))
        if (.not. lnative) then
          allocate (pr (nzr))
          pr(1) = log(pmin)
          do i = 1, nzr - 1
            pr(i+1) = log(pmin) + i*delta_p
          end do
          if (present(fr)) then
            allocate (fr(nzr))
!NEC$ ivdep
            do i = 1, nzr
              fr(i) = sq103 * sqrt(2._wp) * pr(nzr+1-i)
            end do
          end if
        endif

        flr   = lv_surf
        r_0   = log (pmax)
        flr_d = 0._wp
        if (present(d)) d = -sq103 * sqrt(2._wp) * delta_p
        if (present(b)) b =  sq103 * sqrt(2._wp) * log(pmax)

      else
        !----------------------------
        ! for linear spacing in ln(p)
        !----------------------------
        if (lv_surf <= 0._wp) call finish ("set_lv", "lv_surf <= 0")
        rb    = log (pmax)
        rt    = log (pmin)
        flr_d = (lv_top - lv_surf) / (rt - rb)
        r_0   = rb - (lv_surf / flr_d)
        s1    = 1._wp / (rb - r_0)
        s2    = -(rb-r_0)/lv_surf
        fb    = s2 * log ((rb-r_0)/s1)
        ft    = s2 * log ((rt-r_0)/s1)
        if (nzr < 1) nzr = ceiling (ft - fb) + 1
        df    = (ft - fb) / (nzr-1)
        if(dace% lpio) then
          write(6,*) '   r_0, s1, s2, df = ', r_0, s1, s2, df
          write(6,*)
        end if
        !------------------------------------------------
        ! allocate arrays for precalculated length scales
        ! and pressure values on coarse grid
        !------------------------------------------------
        allocate (flr(nzr))
        if (.not. lnative) allocate (pr (nzr))
!NEC$ ivdep
        do i = 1, nzr
          f = fb + (i-1) * df
          if (lnative) then
            r = pr (nzr+1-i)
          else
            r = s1*exp(f/s2)+r_0
            pr  (nzr+1-i) = r
          endif
          flr (nzr+1-i) = flr_d * (r - r_0)
        end do
        if (present (fp)) fp =   edsq103 * s2
        if (present (f0)) f0 = - edsq103 * s2 * log(-s1)
        if (present (b )) b  =   edsq103 * fb
        if (present (d )) d  =   edsq103 * df
        if (present (fr)) then
          allocate (fr(nzr))
          do i = 1, nzr
            fr (i) = edsq103 * (fb + (i-1) * df)
          end do
        endif
      end if
      !-----------------------------------------------------------
      ! return pressure values and modified length scale, printout
      !-----------------------------------------------------------
      pr = exp(pr)
      if(dace% lpio) then
        write(6,*) '   p, log(p), lv(p)'
        do i = 1, nzr
          write(6,*) i, pr(i), log(pr(i)), flr(i)
        end do
        write(6,*)
      end if
      !--------------------------------------------------
      ! scale parameters for use in Gaspari&Cohn function
      !--------------------------------------------------
      flr   = sq103 * flr
      flr_d = sq103 * flr_d
    end if

  end subroutine set_lv

!==============================================================================

  subroutine set_lh (nzr, lh_surf, lh_top, h_0, flh_d, pr, flh, enkf_mean)
  !---------------------------------------------------------------------
  ! set horizontal localisation length scale constant or linear in log(p)
  ! special values for nzr (number of height levels):
  !   nzr > 1 : number of height levels prescribed (coarse grid)
  !   nzr = 1 : number of height levels = enkf_mean% ub(3) (regular grid)
  !---------------------------------------------------------------------
  integer          ,intent(in) :: nzr       ! number of height levels
  real(wp)         ,intent(in) :: lh_surf   ! hor. localisation length at surface
  real(wp)         ,intent(in) :: lh_top    ! hor. localisation length at top
  real(wp)         ,intent(out):: h_0       ! log(p) reference
  real(wp)         ,intent(out):: flh_d     ! localisation length increment
  real(wp)         ,pointer    :: pr (:)    ! p-levels on reduced grid (Pa)
  real(wp)         ,pointer    :: flh(:)    ! hor. loc.length on reduced grid
  type(t_atm)      ,intent(in) :: enkf_mean ! ensemble mean (for pressure)

  !--------------------------------------------------------------------
  ! Method:
  !
  ! The localisation length scale 'flh' (in km) is required to
  ! vary linearly in r=log(p):
  ! flh = flh_d * (r - r0)
  ! with appropriate constants flh_d and h0 so that:
  ! lh_top  = flh (log(pmin)) = flh_d * (log(pmin) - r0)
  ! lh_surf = flh (log(pmax)) = flh_d * (log(pmax) - r0) .
  !
  ! For a given number of grid-points (nzr) of the 'coarse' grid to be used in
  ! the LETKF calculations has previously been determined by subroutine
  ! set_lv.
  !
  ! Final output variables flh_d and flh are multiplied by
  ! sqrt(10./3.) to be used as parameters in the Gaspari& Cohn
  ! localisation function.
  !---------------------------------------------------------------------


    !----------------
    ! local variables
    !----------------
    integer  :: i               ! loop index
    real(wp) :: pmin, pmin_lim  ! min. pressure (fg,limit)
    real(wp) :: pmax, pmax_lim  ! max. pressure (fg,limit)
    real(wp) :: r               ! log (pressure)
    real(wp) :: rt,rb           ! aux. variables
!   real(wp) :: s1,s2,fb,ft     ! aux. variables
!   real(wp) :: df,f            ! aux. variables
    integer  :: nz              ! number of vertical levels

    real(wp) ,parameter ::   sq103 = sqrt (10._wp/3._wp)
!   real(wp) ,parameter :: edsq103 = 1._wp / sq103

    nullify (flh)

    if (nzr >= 1) then      ! coarse/reduced grid used if nzr>1, else no localisation
      !-------------------------------
      ! derive number of levels
      !-------------------------------
      nz       =         enkf_mean% ub(3)
      !---------------------------------------------
      ! derive vertical coordinate (linear in log p)
      !---------------------------------------------
      pmax_lim = 105000._wp               ! max press. on diagn. levels
      pmin_lim = 4000._wp                 ! min press. on diagn. levels

      pmax     = maxval (enkf_mean% pf(:,:,nz,:))
      pmin     = minval (enkf_mean% pf(:,:, 1,:))

      pmax     = p_max  (pmax)
      pmax     = max    (pmax,pmax_lim)   ! check if true pmax is larger
      pmin     = p_min  (pmin)
      pmin     = min    (pmin,pmin_lim)   ! check if true pmin is smaller

      !---------------------------------------------------
      ! set # levels (if not prescribed) and log(p) values
      !---------------------------------------------------
      if (lh_surf == lh_top) then
        !------------------------------------------------
        ! allocate arrays for precalculated length scales
        ! and pressure values on coarse grid
        !------------------------------------------------
        allocate (flh(nzr))
        flh   = sq103 * lh_surf
        flh_d = 0._wp
!       h_0   = 0._wp
      else
        !----------------------------
        ! for linear spacing in ln(p)
        !----------------------------
        if(dace% lpio) then
          write(6,*)
          write(6,*) ' set_lh: set horizontal localisation length scale'
          write(6,*)
          write(6,*) '   pmin, pmax, lh_surf, lh_top = ' &
                       , pmin, pmax, lh_surf, lh_top
          write(6,*)
        end if
        rb    = log (pmax)
        rt    = log (pmin)
        flh_d = (lh_top - lh_surf) / (rt - rb)
        h_0   = rb - (lh_surf / flh_d)
!       s1    = 1._wp / (rb - h_0)
!       s2    = -(rb-h_0)/lh_surf
!       fb    = s2 * log ((rb-h_0)/s1)
!       ft    = s2 * log ((rt-h_0)/s1)
!       df    = (ft - fb) / (nzr-1)
        !---------------------------------------------------------------
        ! allocate arrays for precalculated length scales on coarse grid
        !---------------------------------------------------------------
        if (associated (pr)) then
          allocate (flh(nzr))
!NEC$ ivdep
          do i = 1, nzr
!           f = fb + (i-1) * df
!           r = s1*exp(f/s2)+h_0
            r = log (pr (nzr+1-i))
            flh (nzr+1-i) = flh_d * (r - h_0)
          end do
          !---------------------------------------------------
          ! printout pressure values and modified length scale
          !---------------------------------------------------
          if (dace% lpio) then
            write(6,*) '   p, log(p), lh(p)'
            do i = 1, nzr
              write(6,*) i, pr(i), log(pr(i)), flh(i)
            end do
            write(6,*)
          end if
          !--------------------------------------------------
          ! scale parameters for use in Gaspari&Cohn function
          !--------------------------------------------------
          flh   = sq103 * flh
          flh_d = sq103 * flh_d
        end if
      end if
    end if

  end subroutine set_lh

!==============================================================================

  subroutine set_w_f (nzr, obs, pr, enkf_mean, nz, v_coarse)
  !------------------------------------------------
  ! Interpolate observation weight from COSMO grid
  ! to coarse grid used for vertical interpolation.
  !------------------------------------------------
  integer             ,intent(inout) :: nzr       ! number of height levels
  type(t_obs_set)     ,intent(inout) :: obs       ! observation data
  real(wp)            ,pointer       :: pr (:)    ! p-levels on reduced grid (Pa)
  type(t_atm),optional,intent(in)    :: enkf_mean ! ensemble mean (for pressure)
  integer             ,intent(in)    :: nz        ! number of height levels read
  logical             ,intent(in)    :: v_coarse  ! use coarse grid

    !----------------
    ! local variables
    !----------------
    integer  :: l,ll            ! loop index
!   real(wp) :: y2_tmp (nz)     ! second derivative of function
!   real(wp) :: p               ! actual pressure value

    if ( nz /= nzr .and..not. v_coarse)               &
      call finish('set_w_f','nz/=nzr and not v_coarse')
    if (nzr > 1) then      ! coarse/reduced grid is used
      !------------------------------------------------------------------------------
      ! calculation of coarse weighting function vector for each radiance observation
      !------------------------------------------------------------------------------
      do l = 1, obs% o(1)% n_spot
        if (obs% o(1)% spot(l)% hd% obstype == OT_RAD) then
          do ll = obs% o(1)% spot(l)% o% i + 1, obs% o(1)% spot(l)% o% i + obs% o(1)% spot(l)% o% n
!c          allocate (obs% o(1)% body(ll)% w_f_coarse(nzr))
!c          if (v_coarse) then
!c            !----------------------------------------------------
!c            ! Jakobian was derived on different grid: interpolate
!c            !----------------------------------------------------
!c            call init_splinex (obs% o(1)% spot(l)% p_wf(1:nz), obs% o(1)% body(ll)% w_f(1:nz), y2_tmp)
!c            do i = 1,nzr
!c              if      (pr(i) < obs% o(1)% spot(l)% p_wf(1))  then
!c                obs% o(1)% body(ll)% w_f_coarse(i) = 0._wp
!c              else if (pr(i) > obs% o(1)% spot(l)% p_wf(nz)) then
!c                obs% o(1)% body(ll)% w_f_coarse(i) = obs% o(1)% body(ll)% w_f(nz)
!c              else
!c                call splint (obs% o(1)% spot(l) % p_wf(1:nz),  &
!c                             obs% o(1)% body(ll)% w_f (1:nz),  &
!c                             y2_tmp, pr(i),                    &
!c                             obs% o(1)% body(ll)% w_f_coarse(i))
!c              endif
!c            enddo
!c          else
!c            !-----------------------------------------------
!c            ! Jakobian was derived on same native COSMO grid
!c            !-----------------------------------------------
!c            obs% o(1)% body(ll)% w_f_coarse = obs% o(1)% body(ll)% w_f(1:nz)
!c          endif
          enddo
        endif
      enddo
    end if

  end subroutine set_w_f

!==================================================================================

  elemental function gaspari_cohn (dis, f_l) result (obs_w)
  !-----------------------------------------------------
  !calculates the distance-depending observation weights
  !-----------------------------------------------------
  real(wp)             :: obs_w ! weight
  real(wp) ,intent(in) :: dis   ! distance
  real(wp) ,intent(in) :: f_l   ! length scale times sqrt(10./3.)
    real(wp)           :: x

    if (dis <= 0._wp .or. f_l <= 0._wp) then
      obs_w = 1._wp
    else
      x = dis / f_l
      if (x <= 1._wp) then
!       obs_w = -0.25_wp*x**5 + 0.5_wp*x**4 &
!             + (5._wp/8._wp)*x**3          &
!             - (5._wp/3._wp)*x**2 + 1._wp
        obs_w = (((-0.25_wp * x + 0.5_wp) * x + 5._wp/8._wp) * x  &
              -  5._wp/3._wp) *x*x + 1._wp

      else if (x < 2._wp) then
!       obs_w = (1._wp/12._wp)*x**5 &
!             - 0.5_wp*x**4 + (5._wp/8._wp)*x**3     &
!             + (5._wp/3._wp)*x**2 - 5._wp*x + 4._wp &
!             - (2._wp/3._wp)/x
        obs_w = max ((((((1._wp/12._wp) * x - 0.5_wp) * x + 5._wp/8._wp) * x  &
                        +  5._wp/3._wp) * x  - 5._wp) * x + 4._wp             &
                        - (2._wp/3._wp) / x, 0._wp)
      else
        obs_w = 0._wp
      end if
    end if
  end function gaspari_cohn

!==============================================================================

  subroutine construct_rndm_tcor (ct, scale, time, verbose)
  type (t_rndm_tcor) ,intent(out) :: ct      ! time correlation meta data
  real(wp)           ,intent(in)  :: scale   ! correlation time scale (h)
  type (t_time)      ,intent(in)  :: time    ! actual time
  logical  ,optional ,intent(in)  :: verbose ! steering of printout
  !--------------------------------------
  ! set up temporal correlation meta data
  !--------------------------------------

    real(wp)              :: mo          ! offset of first axiliary time
    real(wp)              :: h           ! time of destination grid (hours)
    integer               :: nm          ! number of time steps (<=5)
    real(wp)              :: w (5)       ! weights for eech time step
    integer               :: ht(5)       ! hour (for unique seed generation)
    real(wp)              :: dh(5)       ! time difference to destination (h)
    real(wp)              :: t_loc       ! scaled correlation time scale
    integer               :: ts          ! time for random pattern generator
    integer               :: j           ! loop index

    h  = hours (time)
    if (scale /= 0) then
      !--------------------------------------------------------
      ! temporal correlation scale >= 1h
      ! interpolate correlation patterns in between time slices
      !--------------------------------------------------------
      mo = mod (h, scale)  ! offset
      ts = int (h/ scale)  ! time slice for random numbers
      if (mo == 0._wp) then
        !-----------------------------------------------------
        ! destination time corresponds with central time slice
        !-----------------------------------------------------
        nm = 5
        mo = -2 * scale
        ts = ts - 2
      else
        !--------------------------------------------
        ! destination time in between two time slices
        !--------------------------------------------
        nm = 4
        mo = - scale - mo
        ts = ts - 1
      endif
      !----------------------------------------------
      ! calculate weights for temporal interpolations
      !----------------------------------------------
      t_loc = sqrt (10._wp/3._wp) * scale / sqrt(2._wp)
      do j = 1, nm
        dh(j) = (j-1) * scale + mo
        w (j) = gaspari_cohn (abs(dh(j)), t_loc)
        ht(j) = (j-1) + ts
      end do
      w(1:nm) = w(1:nm) / sqrt(sum(w(1:nm)**2))
    else
      !------------------------------------------------
      ! no temporal correlations: set no.of slices to 1
      !------------------------------------------------
      mo    = 0._wp
      nm    = 1
      w (1) = 1._wp
      ht(1) = nint (h)
      dh(1) = 0._wp
      t_loc = 0._wp
    endif
    !--------------------------------------------
    ! print out temporal interpolation parameters
    !--------------------------------------------
    if (present (verbose)) then
      if (verbose .and. dace% lpio) then
        write(6,'()')
        write(6,'(a)') '  temporal correlation parameters:'
        write(6,'()')
        write(6,'(a,f5.1 ,a)') '   scale  = ',scale ,' ! time scale (h)'
        write(6,'(a,f5.1 ,a)') '   t_loc  = ',t_loc ,' ! normalised scale'
        write(6,'(a,i3,2x,a)') '   nm     = ',nm    ," ! number of 'modes'"
        write(6,'(a,f5.1 ,a)') '   mo     = ',mo    ,' ! offset'
        write(6,'()')
        do j = 1, nm
          write(6,'(a,i1,2f8.2,i10)')  '   w dh h = ',j,w(j), dh(j), ht(j)
        end do
        write(6,'()')
      endif
    endif
    !-------------------------------
    ! store in derived type variable
    !-------------------------------
    ct% scale = scale
    ct% n     = nm
    ct% w     = w
    ct% ht    = ht

  end subroutine construct_rndm_tcor

!==============================================================================

  subroutine flat_copy (y, x)
  type (t_atm) ,intent(out) :: y
  type (t_atm) ,intent(in)  :: x
  !-------------------------------
  ! flat copy of atmospheric state
  !-------------------------------
    y = x
  end subroutine flat_copy

!=============================================================================
end module mo_letkf_util
