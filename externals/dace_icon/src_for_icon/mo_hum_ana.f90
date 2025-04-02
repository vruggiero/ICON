!
!+ handling of stratospheric humidity
!
MODULE mo_hum_ana
!
! Description:
!   Special handling of stratospheric humidity:
!   Relative humidity increment zero for p < hum_ana_top.
!   For p < tropopause level,
!     set constant specific humidity          (strat_q_mode=1),
!     take specific humidity from first guess (strat_q_mode=2),
!     relaxation to first guess/param.climate (strat_q_mode=3),
!     relaxation to analysis/param.climate    (strat_q_mode=4),
!     relaxation to first guess/climatology   (strat_q_mode=5).
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  Changes for COSMO vertical coordinate
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Harald Anlauf
!  Namelist /HUM_ANA/: rh0_qcl0 (rel.humidity below which we set qcl to 0)
! V1_26        2013/06/27 Harald Anlauf
!  mo_hum_ana: adapt tropopause diagnostics to ICON
! V1_27        2013-11-08 Harald Anlauf
!  implement strat_q_mode=2 to keep specific humidity in stratosphere
! V1_35        2014-11-07 Harald Anlauf
!  set_strat_q: rh_over_ice(flag) for choice of q_rh or q_rhi for limiting
!               stratospheric humidity
!  WMO tropopause diagnostics, revised treatment of stratospheric qv
!               (strat_q_mode=3)
! V1_42        2015-06-08 Harald Anlauf
!  further refinement of q climatology
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! ??????????????  ???  ?????????  original OI sources
! Andreas Rhodin  DWD  2004-2007  port to 3DVAR/PSAS
! Harald Anlauf   DWD  2012
!------------------------------------------------------------------------------

  !-------------
  ! modules used
  !-------------
  use mo_kind,       only: wp, sp          ! working precision kind parameter
  use mo_exception,  only: finish          ! abort on error
  use mo_atm_state,  only: t_atm           ! atmospheric state data type
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_bcast         ! broadcast routine
  use mo_namelist,   only: position_nml,  &! routine to position nml group
                           nnml,          &! namelist fortran unit number
                           POSITIONED      ! position_nml: OK    return flag
  use mo_physics,    only: q_rh,          &! spec.humidity <- rel.humidity
                           q_rhi,         &! spec.humidity <- rel.humidity over ice
                           rd    => R,    &! gas constant of dry air
                           g     => gacc, &! gravity acceleration
                           tmelt => t0c,  &! melting temperature
                           RdRv  => RDRD, &! R_dry / R_vapor (~0.622)
                           r2d,           &! conversion factor radians->degrees
                           d2r             ! conversion factor degrees->radians
  use mo_wmo_tables, only: WMO3_HYBRIDB    ! level type: hybrid
  use mo_time,       only: operator(-),   &! calculate time difference
                           seconds,       &! convert to seconds
                           t_time,        &! derive type
                           iyyddd,        &! year + day-of-year
                           days_year       ! days in current year
  use mo_run_params, only: ana_time,      &! analysis time
                           data,          &! data path
                           path_file       ! add a path to a file name
  !---------------------
  ! NetCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_open,             &
                           nf90_close,            &
                           nf90_get_att,          &
                           nf90_get_var,          &
                           nf90_inq_dimid,        &
                           nf90_inq_varid,        &
                           nf90_inquire_dimension,&
                           nf90_strerror,         &
                           NF90_FILL_REAL,        &
                           NF90_GLOBAL,           &
                           NF90_NOWRITE,          &
                           NF90_NOERR

  implicit none
  !----------------
  ! public entities
  !----------------
  private
  public :: hum_ana_incr     ! set upper level analysis increment to zero
  public :: set_strat_q      ! set stratospheric spec.hum. according to mode
  public :: read_nml_hum_ana ! read namelist /HUM_ANA/
  public :: e_hum_top        ! nominal FG error above hum_ana_top
  public :: hum_ana_top      ! top of humidity analysis [Pa]
  public :: ln_hum_ana_top   ! logarithm of hum_ana_top
  public :: rh0_qcl0         ! rel.hum. below which we set qcl to 0
  public :: strat_q          ! stratospheric humidity
  public :: strat_q_mode     ! stratospheric humidity analysis mode
  public :: strat_q_ens      ! relax stratospheric hum. of ensemble members

  !=================
  ! Module variables
  !=================
  !-------------------
  ! Namelist variables
  !-------------------
  real(wp) :: hum_ana_top    = 25000._wp ! top of humidity analysis
  real(wp) :: e_hum_top      = 0.01_wp   ! nominal FG error above hum_ana_top
  real(wp) :: ln_hum_ana_top = 10.12663110385033780125_wp ! log(hum_ana_top)
  real(wp) :: rh0_qcl0       = 0.9_wp    ! rel.humidity threshold for qcl=0
  real(wp) :: strat_q        = 4.0e-6_wp ! upper stratospheric humidity
  real(wp) :: strat_q_low    = 2.5e-6_wp ! lower stratospheric humidity
  real(wp) :: strat_q_low_tr = 2.5e-6_wp ! lower stratospheric humidity, tropics
  integer  :: strat_q_mode   = 1         ! 0=keep ana, 1=const. value, 2=use fg
                                         ! 3=relaxation to fg/climate profile
                                         ! 4=relaxation to ana/climate profile
                                         ! 5=relaxation to fg/ext.climatology
  logical  :: strat_q_ens    = .false.   ! relax stratospheric hum. of ensemble
  logical  :: rh_over_ice    = .false.   ! use q_rhi instead of q_rhw
  logical  :: wmo_lrt_diag   = .false.   ! use WMO lapse rate tropopause diag.
  real(wp) :: asurf_min      =  2000._wp ! min. TP height above surface [m]
  real(wp) :: pt_min         =  7000._wp ! min. TP pressure             [Pa]
  real(wp) :: pt_fallback    = 30000._wp ! TP height when search fails  [Pa]
  real(wp) :: tau_strat      =    30._wp ! time constant for q relaxation [d]
  real(wp) :: tau_strat_tp   =    0.5_wp ! relaxation time near tropopause [d]
  real(wp) :: tau_strat_scal = 1 / 3._wp ! scaling factor for temp.gradient term
  real(wp) :: pmax_qclim     =100000._wp ! max. pressure to use clim. QV [Pa]
  real(wp) :: qsat_ice_scal  =    1.0_wp ! scaling factor for saturation over ice
  real(wp) :: zlo_clq_tr     = 30000._wp ! lower ref. height of q clim., tropics
  real(wp) :: zlo_clq_pol    = 30000._wp ! lower ref. height of q clim., poles
  real(wp) :: zhi_clq_tr     = 60000._wp ! upper ref. height of q clim., tropics
  real(wp) :: zhi_clq_pol    = 60000._wp ! upper ref. height of q clim., poles
  real(wp) :: exclude_lat(2) =  -999._wp ! latitude range to exclude from relax.
  real(wp) :: exclude_p  (2) =  -999._wp ! pressure range to exclude     [Pa]
  character(128) :: file     = ""        ! humidity climatology file

  namelist /HUM_ANA/ hum_ana_top, e_hum_top, rh0_qcl0,      &
                     strat_q, strat_q_low, strat_q_mode,    &
                     rh_over_ice, wmo_lrt_diag, asurf_min,  &
                     pt_min, pt_fallback, strat_q_low_tr,   &
                     tau_strat, zlo_clq_pol, zlo_clq_tr,    &
                     zhi_clq_pol, zhi_clq_tr, strat_q_ens,  &
                     tau_strat_tp, tau_strat_scal,          &
                     pmax_qclim, qsat_ice_scal,             &
                     exclude_lat, exclude_p, file

  !------------------------------------------------------
  ! Derived type definition:
  ! This structure is used to hold an externally provided
  ! humidity climatology in the stratosphere
  !------------------------------------------------------
  type t_hum_clim
     integer               :: nlat = 0      ! number of latitudes
     integer               :: nlev = 0      ! levels
     real(wp)              :: dlat = 0._wp  ! latitude spacing
     real(wp), allocatable :: lat(:)        ! latitude [deg]    (nlat)
     real(wp), allocatable :: lnp(:)        ! level (log(p/Pa)) (nlev)
     real(wp), allocatable :: qv (:,:)      ! climatological qv (nlev,nlat)
  end type t_hum_clim

  type(t_hum_clim), save   :: hum_clim      ! Humidity climatology

  !-----------
  ! Interfaces
  !-----------
  interface set_strat_q
    module procedure set_strat_q_0  ! scalar version
    module procedure set_strat_q_1  ! vector version
    module procedure set_strat_q_2  ! vector/vector version
  end interface

  !----------------------------------------------------------------------------
  ! Beware: these parameters are implemented with wrong units (km instead of m)
  ! Keep these as default (bug compatibility)
  real(wp) ,parameter :: z30 = 30._wp  ! Reference heights of ...
  real(wp) ,parameter :: z60 = 60._wp  ! ... climatological humidity profile
  !----------------------------------------------------------------------------

contains
!------------------------------------------------------------------------------
  subroutine hum_ana_incr (aninc, bg)
  type (t_atm) ,intent(inout) :: aninc ! analysis increment
  type (t_atm) ,intent(in)    :: bg    ! background (for pressure)
  !-----------------------------------------------------
  ! modify analysis increment:
  ! relative humidity increment zero for p < hum_ana_top
  !-----------------------------------------------------
    integer :: i,j,k,l
    if (0._wp >= hum_ana_top) return
    if (.not. associated (bg% pf)) &
      call finish ('hum_ana_incr','pf is not set')
    if (.not. associated (aninc% rh)) &
      call finish ('hum_ana_incr','rh is missing')
    do       l = aninc% lb(4), aninc% ub(4)
      do     j = aninc% lb(2), aninc% ub(2)
        do   i = aninc% lb(1), aninc% ub(1)
          do k = aninc% lb(3), aninc% ub(3)      ! From model top to bottom
            if (bg% pf (i,j,k,l) >= hum_ana_top) exit
            aninc% rh (i,j,k,l) = 0._wp
          end do
        end do
      end do
    end do
  end subroutine hum_ana_incr
!------------------------------------------------------------------------------
!  subroutine set_strat_q (atm)
!  type (t_atm) ,intent(inout) :: atm
!  !----------------------------------------------------
!  ! set stratospheric humidity:
!  ! constant for p < tropopause level.
!  ! tropopause level tropo_p is defined by the extremum
!  ! of d2T/dp2 between 7000. and 50000. hPa
!  !----------------------------------------------------
!    integer  :: i,j,k,l
!    real(wp) :: d2t (atm% ub(3)), d2tmax, rh
!    integer  :: ktr, ke
!    if (strat_q < 0._wp) return
!    ke = atm% ub(3)
!    do       l = atm% lb(4), atm% ub(4)
!      do     j = atm% lb(2), atm% ub(2)
!        do   i = atm% lb(1), atm% ub(1)
!          !---------------------------------------
!          ! determine the height of the tropopause
!          !---------------------------------------
!          call init_splinex (log(atm% pf(i,j,:,l)), atm% t (i,j,:,l), d2t)
!          ktr    = 1
!          d2tmax = 0._wp
!          do k = 2, ke
!            if (atm% pf(i,j,k,l) > 50000._wp) exit
!            if (atm% pf(i,j,k,l) <  7000._wp) cycle
!            if (d2t(k)         > d2tmax           .and. &
!                atm%t(i,j,k,l) < atm%t(i,j,k-1,l) .and. &
!                atm%t(i,j,k,l) < atm%t(i,j,k+1,l)       ) then
!              d2tmax = d2t(k)
!              ktr    = k
!            endif
!          end do
!          !------------------------------------------------------
!          ! above tropopause level set constant specific humidity
!          !------------------------------------------------------
!          do k = 1, ktr-1
!            rh = rh_q (atm% q(i,j,k,l), atm% t(i,j,k,l), atm% pf(i,j,k,l))
!            atm% q(i,j,k,l) = strat_q
!            !------------------------------------------------
!            ! restrict stratospheric relative humidity to 90%
!            !------------------------------------------------
!            if (rh > 0.9_wp) atm% q(i,j,k,l) = atm% q(i,j,k,l) * 0.9_wp / rh
!          end do
!        end do
!      end do
!    end do
!  end subroutine set_strat_q
!------------------------------------------------------------------------------
  subroutine set_strat_q_1 (atm, ref)
    type (t_atm) ,intent(inout) :: atm (:)
    type (t_atm) ,intent(in)    :: ref     ! Reference state for spec.humidity
    optional                    :: ref
    integer :: i
    do i = 1, size (atm)
       call set_strat_q_0 (atm(i), ref)
    end do
  end subroutine set_strat_q_1
!------------------------------------------------------------------------------
  subroutine set_strat_q_2 (atm, ref)
    type (t_atm) ,intent(inout) :: atm (:)
    type (t_atm) ,intent(in)    :: ref (:) ! Reference states for spec.humidity
    integer :: i
    do i = 1, size (atm)
       call set_strat_q_0 (atm(i), ref(i))
    end do
  end subroutine set_strat_q_2
!------------------------------------------------------------------------------
  subroutine set_strat_q_0 (atm, ref)
  type (t_atm) ,intent(inout) :: atm
  type (t_atm) ,intent(in)    :: ref     ! Reference state for spec.humidity
  optional                    :: ref
  !----------------------------------------
  ! set stratospheric humidity:
  ! This routine is taken over from the OI.
  !----------------------------------------
    integer           :: i,j,jl,l        ! grid point indices
    real(wp)          :: zpt2min         ! lowest pressure level for tropopause
    real(wp)          :: zpt2max         ! high.  pressure level for tropopause
    real(wp)          :: zpf(atm%ub(3)+1)! pressure at full levels
    real(wp) ,pointer :: tana (:)        ! temperature
    real(wp) ,pointer :: an_qv(:)        ! specific humidity
    real(wp) ,pointer :: aa(:), bb(:)    ! hybrid pressure level coefficients
    real(wp)          :: zpsana          ! extract analysed surface pressure
    real(wp)          :: qstrat          ! prescribed stratospheric spec. hum.
    integer           :: nlevm           ! number of model levels
    integer           :: ilst, ilnd      ! do loop bounds
    real(wp)          :: zdz             ! height differences betw. full levels
!   real(wp)          :: zdz1(atm%ub(3)) ! height differences betw. full levels
    real(wp)          :: zdt1(atm%ub(3)) ! derivative of t with respect to z
    integer           :: iltrop
    integer           :: iltroph
!   real(wp)          :: zttrop
!   real(wp)          :: zttroph
    real(wp)          :: zdt12
    real(wp)          :: cenlat, cenlon  ! latitude, longitude (degrees)
    real(wp)          :: zqsat           ! saturation specific humidity
    real(wp)          :: zqm
    real(wp)          :: ref_q(atm%ub(3))! Reference specific humidity
    logical           :: l_q_ref         ! Use stratospheric q from reference
    logical           :: excl_reg        ! Exclude region from relaxation?
    logical           :: excl_lat        ! Latitude in exclusion region?
    logical           :: excl_p          ! Pressure value in exclusion region?
    !-------------------------------------------------------------------------
    ! Work variables for tropopause diagnostics based on proper WMO definition
    !-------------------------------------------------------------------------
    integer, parameter:: mtp = 2         ! Max. number of accepted tropopauses
    integer           :: ntp             ! Tropopause count
    integer           :: ltp(mtp)        ! Model levels of tropopause(s)
    logical           :: scan            ! Scan for tropopause?
    logical           :: lthr            ! Lapse rate threshold exceeded?
    logical           :: crit            ! Lapse rate criterion always satisfied?
    integer           :: kl              ! Search index
    real(wp)          :: tt              ! (Tropopause) level temperature
    real(wp)          :: zt              ! (Tropopause) level height
    real(wp)          :: zfmin           ! Min. tropopause height
    real(wp)          :: zf  (atm%ub(3)) ! height of full level [m]
    real(wp)          :: dtdz            ! lapse rate dT/dz [K/m]
    integer           :: jtop            ! Model level above last scanned
    !---------------------------------------------------------
    ! Work variables for revised treatment of stratospheric QV
    !---------------------------------------------------------
    integer           :: ktp             ! Level of accepted thermal tropopause
    integer           :: k2, k6          ! Levels approx. 2km and 6km above "
    integer           :: k               ! Auxiliary index
    integer           :: nl              ! No. levels in climatology
    real(wp)          :: alpha           ! Relaxation rate
    real(wp)          :: w               ! Relaxation weight
    real(wp)          :: dt              ! Forecast time
    real(wp)          :: qv6             ! Lower stratospheric humidity
    real(wp)          :: lnp             ! log(p)
    !-------------------------------------------
    ! Climatological specific humidity profiles:
    !-------------------------------------------
    real(wp)          :: cl_q_pol(atm%ub(3)) ! climatology at the poles
    real(wp)          :: cl_q_tr (atm%ub(3)) ! climatology in the tropics
    real(wp)          :: cl_q    (atm%ub(3)) ! local climatology profile
    real(wp), allocatable :: cl_q_tmp(:)     ! Interpolated climatology profile

    l_q_ref = .false.
    select case (strat_q_mode)
    case (0)
       return                            ! Do nothing
    case (1)
       if (strat_q < 0._wp) return       ! immediately return if strat_q < 0
    case (2)
       if (.not. present (ref)) &
            call finish ('set_strat_q','strat_q_mode=2, but ref is missing')
       if (.not. associated (ref% q)) &
            call finish ('set_strat_q','ref%q not associated')
       l_q_ref = .true.
    case (3,4,5)
       if (.not. wmo_lrt_diag) &
            call finish ('set_strat_q','strat_q_mode>=3 requires wmo_lrt_diag=T')
       if (.not. rh_over_ice) &
            call finish ('set_strat_q','strat_q_mode>=3 requires rh_over_ice=T')
       if (.not. present (ref)) &
            call finish ('set_strat_q','strat_q_mode>=3, but ref is missing')
       if (.not. associated (ref% q)) &
            call finish ('set_strat_q','ref%q not associated')
       if (strat_q_mode == 5) then
          if (hum_clim% nlat == 0 .or. hum_clim% nlev == 0) &
            call finish ('set_strat_q','strat_q_mode==5, but no climatology')
          allocate (cl_q_tmp(hum_clim% nlev))
       end if
       l_q_ref = .true.
    case default
       call finish ('set_strat_q','invalid strat_q_mode')
    end select
    if (atm% grid% levtyp == WMO3_HYBRIDB) then
       aa => atm% grid% ak
       bb => atm% grid% bk
       if (.not. associated (atm% ps)) &
            call finish ('set_strat_q_0','ps is not present !')
    else
       !----------------------------------------
       ! COSMO / ICON height-based hybrid levels
       !----------------------------------------
       nullify (aa, bb)
       if (.not.associated (atm% pf)) &
            call finish('set_strat_q_0','pf is not present !')
    end if
    if (wmo_lrt_diag) then
       if (.not. associated (atm% geof)) &
            call finish('set_strat_q_0','geof is not allocated!')
    end if
    !---------------------------
    ! Check for exclusion region
    !---------------------------
    excl_reg = (exclude_lat(1) < exclude_lat(2) &
         .and.  exclude_lat(1) <     90._wp     &
         .and.  exclude_lat(2) >    -90._wp     &
         .and.  exclude_p  (1) < exclude_p  (2) &
         .and.  exclude_p  (1) < 100000._wp     &
         .and.  exclude_p  (2) >      0._wp     )
    !-----------------------
    ! Derive relaxation time
    !-----------------------
    if (present (ref)) then
       dt = seconds (atm% time - ref% ref_time)
!      if (dace% lpio) then
!         print *, "### dt[s] =", dt
!         print *, "### atm% time =", (atm% time)
!         print *, "### ref% ref_time =", (ref% ref_time)
!      end if
       dt = min (dt, 30*86400._wp)      ! Limit dt to 30 days
       if (dt <= 0._wp) dt = 10800._wp  ! Default: 3 hours.
    else
       dt = 0._wp
    end if
    ref_q = -HUGE (0._wp)               ! Set to undef.
    !--------------------------------
    ! loop over horizontal gridpoints
    !--------------------------------
    do       l = atm% lb(4), atm% ub(4)
      do     j = atm% lb(2), atm% ub(2)
        do   i = atm% lb(1), atm% ub(1)
          !----------------------------------------------------
          ! pick up quantities from atmospheric state data type
          !----------------------------------------------------
          nlevm   =  atm% ub(3)                      ! number of model levels
          tana    => atm% t         (i,j,:,l)        ! temperature
          an_qv   => atm% q         (i,j,:,l)        ! specific humidity
          cenlat  =  atm% grid% rlat(i,j,1,l) * r2d  ! latitude  (degrees)
!         cenlon  =  atm% grid% rlon(i,j,1,l) * r2d  ! longitude (degrees)
          !----------------------------
          ! set pressure at full levels
          !----------------------------
          select case (atm% grid% levtyp)
          case (WMO3_HYBRIDB)
            zpsana = atm% ps        (i,j,1,l)        ! surface pressure
            !------------------------
            ! GME / HRM hybrid levels
            !------------------------
            DO 435 jl = 1, nlevm+1 ! nlvmp1
              zpf (jl) = 0.5_wp * (aa (max (1, jl - 1) ) + aa (jl) ) &
                       + 0.5_wp * (bb (max (1, jl - 1) ) + bb (jl) ) * zpsana
  435       END DO
            !-------------------------------------------
            !*  highest level uses unmodified aa (or bb)
            !-------------------------------------------
!           zpf (2) = 0.5_wp * (0._wp + aa (2) ) + (0._wp + bb (2) ) * zpsana ! ???
          case default
            !----------------------------------------
            ! COSMO / ICON height-based hybrid levels
            !----------------------------------------
            zpf (2:) = atm% pf (i,j,:,l)
            zpf (1)  = zpf (2)
            zpsana   = zpf (nlevm+1)          ! proxy for surface pressure
          end select
          !================================================
          !*   7.8   modify stratospheric specific humidity
          !
          !* assign min and max pressure of tropopause,
          !  stratospheric mixing ratio and max rel hum.
          !=================================================
!         zpt2min =  max (zpf (2), 7000._wp)  ! pressure level range
          zpt2min =  max (zpf (2),   pt_min)  ! pressure level range
          zpt2max =  min (zpsana, 50000._wp)  !   for tropopause

          IF (abs (zpf (1) - zpf (2) ) .lt. 0.000001_wp) THEN
            ilst = 3
            ilnd = nlevm - 2
          ELSE
            ilst = 2
            ilnd = nlevm - 1
          END IF
          !===============================================
          !  7.81  first derivative of t with respect to z
          !===============================================
          DO 781 jl = ilst, nlevm
            !---------------------------------------------
            ! height differences between full model levels
            !---------------------------------------------
            zdz = rd / g * 0.5_wp * (tana (jl - 1) + tana (jl) ) &
                * log (zpf (jl + 1) / zpf (jl) )
            zdt1 (jl) = (tana (jl - 1) - tana (jl) ) / zdz
!           zdz1 (jl) = zdz
  781     END DO
          iltrop  = 0
          iltroph = 0
!         zttrop  = 0._wp
!         zttroph = 0._wp

          !--------------------------------------
          ! WMO Lapse Rate Tropopause diagnostics
          !--------------------------------------
          if (wmo_lrt_diag) then
             zf(:) = atm% geof(i,j,1:nlevm,l) / g
             zfmin = zf(nlevm) + asurf_min
             ntp   = 0
             ltp   = 0
             scan  = .true.
             crit  = .true.
scan_lrt:    do jl = ilnd, ilst, -1
                if (zf (jl)   < zfmin  ) cycle scan_lrt
                if (zpf(jl+1) > zpt2max) cycle scan_lrt
                if (zpf(jl+1) < zpt2min) exit  scan_lrt

                tt   =  tana(jl)
                zt   =  zf  (jl)
                dtdz = (tana(jl-1) - tt) / (zf(jl-1) - zt)
                lthr =  dtdz > -2.e-3_wp    ! Lapse rate above threshold?
                crit = (crit .and. lthr)    ! Remember any drop below threshold
                !------------------------------------------------
                ! Scan for lapse rate tropopause:
                ! Lapse rate > -2 K/km for all levels within 2km,
                ! but not for the level below.
                !------------------------------------------------
                if (scan) then
                   if (.not. lthr)                                cycle scan_lrt
                   if ((tt-tana(jl+1))/(zt-zf(jl+1)) > -2.e-3_wp) cycle scan_lrt
                   do kl = jl-2, ilst-1, -1
                      if (zf(kl) > zt + 2000._wp) exit
                      if ((tana(kl)-tt)/(zf(kl)-zt)  < -2.e-3_wp) cycle scan_lrt
                   end do
                   ntp      = ntp + 1    ! Tropopause count
                   ltp(ntp) = jl
                   if (ntp >= mtp)       exit  scan_lrt
                   scan     = .false.

                !----------------------------------------------
                ! Check for region separating tropopauses:
                ! Lapse rate < -3K/km for all levels within 1km
                !----------------------------------------------
                else
                   if (dtdz > -3.e-3_wp) cycle scan_lrt
                   do kl = jl-2, ilst-1, -1
                      if (zf(kl) > zt + 1000._wp) exit
                      if ((tana(kl)-tt)/(zf(kl)-zt)  > -3.e-3_wp) cycle scan_lrt
                   end do
                   scan = .true.         ! Scan for next tropopause
                end if
             end do scan_lrt
             jtop = jl
             !----------------------------------------------------
             ! Fallback if no tropopause found:
             ! - use pt_fallback if dT/dz > -2 K/km in full range,
             ! - use pt_min otherwise
             !----------------------------------------------------
             if (ntp == 0) then
                if (crit) then
                   ltp(1) = minloc (abs (zpf(2:ilnd) - pt_fallback), 1)
                else
                   ltp(1) = minloc (abs (zpf(2:ilnd) - pt_min     ), 1)
                end if
             end if
             !-------------------------------------------------------
             ! If secondary tropopause found, use it if below 200 hPa
             !-------------------------------------------------------
             iltrop  = ltp(1)
             iltroph = ltp(2)
             if (ntp == 2 .and. iltroph > 0) then
                if (zpf(iltroph+1) > 20000._wp) iltrop = iltroph
             end if

#if 0  /* DEBUG */
cenlon  =  atm% grid% rlon(i,j,1,l) * r2d  ! longitude (degrees)
if(ntp==2) then
kl=ltp(1)
if(iltrop == iltroph) then
write(3000+dace% pe,*) real(cenlat),real(cenlon),ltp,real(zpf(kl+1)/100),real(zpf(ltp(2)+1)/100),&
     iltrop == iltroph !, real(zf(iltrop))
write(3000+dace% pe,*) "dtdz:",(real((tana(jl-1)-tana(jl))/(zf(jl-1)-zf(jl))*1000),jl=kl,ltp(2)-1,-1)
write(3000+dace% pe,*) "DtDz:",(real((tana(jl-1)-tana(kl))/(zf(jl-1)-zf(kl))*1000),jl=kl,ltp(2)-1,-1)
write(3000+dace% pe,*) "  Dz:",(real((zf(jl-1)-zf(kl))),jl=kl,ltp(2)-1,-1)
write(3000+dace% pe,*)
else
write(2000+dace% pe,*) real(cenlat),real(cenlon),ltp,real(zpf(kl+1)/100),real(zpf(ltp(2)+1)/100),&
     iltrop == iltroph !, real(zf(iltrop))
end if
else if (ntp == 0) then
write(9000+dace% pe,*) real(cenlat),real(cenlon),ntp,ltp(1),real(zpf(iltrop+1)/100), &
     "p_above =", real(zpf(jtop+1)/100)
write(9000+dace% pe,*) "dtdz:",(real((tana(jl-1)-tana(jl))/(zf(jl-1)-zf(jl))*1000),jl=ltp(1),jtop,-1)
write(9000+dace% pe,*)
else if (ntp == 1) then
write(1000+dace% pe,*) real(cenlat),real(cenlon),ntp,ltp(1),real(zpf(iltrop+1)/100)!,real(zf(iltrop))
end if
#endif

          else   ! .not. wmo_lrt_diag
          !==============================================================
          !  7.82  the tropopause is  defined as the lowest level above
          !        which the temperature lapse rate is greater than
          !        -2 K / km and with temperature less than -45 degree C.
          !        for temperatures at all levels > -45 C, take only
          !        lapse rate condition
          !        find primary and secondary tropopauses
          !==============================================================
          DO 782 jl = ilnd, ilst, - 1
            IF (zpf(jl+1) .lt. zpt2min .or. zpf(jl+1) .gt. zpt2max) CYCLE
            IF (tana (jl) .lt. tmelt - 45._wp) THEN
              IF (       ( zpf (jl + 1) .LT. 20000._wp )                &
                   .AND. ( iltrop          .EQ. 0         )) THEN
                iltrop = jl
              END IF
              IF (zdt1 (jl) .gt. - 2.e-3_wp) THEN
                iltrop = jl
!               zttrop = tana (jl) - tmelt
                GOTO 7821
              END IF
            ELSE
              IF (        ( zpf (jl) .LT. 35000._wp )        &
                    .AND. ( zpf (jl) .GT. 7000._wp  ))              THEN
                IF (tana (jl) .lt.tmelt - 40._wp) THEN
                  IF (zdt1 (jl) .gt. - 2.e-3_wp) THEN
                    IF (iltroph.eq.0) THEN
                      iltroph = jl
!                     zttroph = tana (jl) - tmelt
                    END IF
                  END IF
                END IF
              END IF
            END IF

  782       END DO
 7821       CONTINUE
            !----------------------------------------------------------
            ! if two tropopauses are found , take the lower one of both
            !----------------------------------------------------------
            IF (iltrop .ne. 0 .and. iltroph .ne. 0) THEN

!             IF (iltroph - iltrop.gt.2) THEN
!               print *,'mype ',myproc,' thgrpev.f two tropopauses found ',&
!                'lev1 ',iltrop ,zttrop ,' [c] ',                          &
!                'lev2 ',iltroph,zttroph,' [c] ',                          &
!                'irow ',irow,' jlon=',jlon
!              END IF

              !---------------------------------------------
              ! height differences between full model levels
              !---------------------------------------------
              IF (iltrop.lt.iltroph) THEN
                zdz = rd / g * 0.5_wp * (tana (iltrop) + tana (iltroph) ) &
                    * log (zpf (iltrop) / zpf (iltroph) )
                zdt12 = (tana (iltrop) - tana (iltroph) ) / zdz

                IF (zdt12.gt. - 2.e-3_wp) THEN
                  IF (cenlat .lt. - 35._wp .or. cenlat .gt. 35._wp) THEN

!                   IF (iltroph - iltrop.gt.2) THEN
!                     print *,'mype ',myproc,' thgrpev.f two tropopauses ',&
!                       ' extratropics, take lower one',                   &
!                       'lev1 ',iltrop ,zttrop ,' [c] ',                   &
!                       'lev2 ',iltroph,zttroph,' [c] ',                   &
!                       'irow ',irow,' jlon=',jlon
!                   END IF

                    !------------------------------------------
                    ! take the first one(lower one in z-system)
                    !------------------------------------------
                    iltrop = iltroph

                  ELSE

                    IF (iltroph - iltrop.gt.2) THEN
!!!                   PRINT * , 'mype ', dace% pe,                         &
!!!                    ' thgrpev.f two tropopauses ',                      &
!!!                    ' tropics, take lower one', 'lev1 ', iltrop,        &
!!!                    zttrop, ' [c] ', 'lev2 ', iltroph, zttroph, ' [c] ',&
!!!                    'lat ', cenlat, ' lon=', cenlon
                    END IF
                    !------------------------------------------
                    ! take the first one(lower one in z-system)
                    !------------------------------------------
                    iltrop = iltroph
                  END IF
                END IF
              END IF

            END IF
            !------------------------------------------------------------
            ! if no tropopause found yet, take the tropopause
            ! for the lapse rate condition
            ! otherwise take the default for the middle of the atmsophere
            !------------------------------------------------------------
            IF (iltrop.eq.0) THEN

!             print *,'mype ',myproc,' thgrpev.f 7.82 no standard ', &
!               'tropopause found: take auxiliary tropopause ',      &
!               'irow ',irow,' jlon=',jlon,                          &
!               'iltrop ',iltrop,' iltroph',iltroph

              IF (iltroph.ne.0) THEN
                iltrop = iltroph
              ELSE
                iltrop = nlevm / 2
!!!             PRINT * , 'mype ', dace% pe,' thgrpev.f 7.82 neither     ', &
!!!              'tropopause found: take nlevm/2', 'lat ', cenlat, ' lon=', &
!!!              cenlon,'iltrop ', iltrop, ' iltroph', iltroph
              END IF
            END IF

          end if  ! wmo_lrt_diag

          !---------------------------------------
          ! Set stratospheric humidity in analysis
          !---------------------------------------
          if (l_q_ref) then
            if (strat_q_mode == 4) then
              ref_q(1:iltrop) = an_qv(     1:iltrop  )  ! use QV analysis as reference
            else
              ref_q(1:iltrop) = ref% q(i,j,1:iltrop,l)  ! use FG as reference humidity
            endif
          end if

          select case (strat_q_mode)
          case (3,4,5)
             ktp  = iltrop
             k2   = minloc (abs (zf(:) - (zf(ktp) + 2000._wp)), 1)
             k6   = minloc (abs (zf(:) - (zf(ktp) + 6000._wp)), 1)
             k2   = min (k2, ktp-1)      ! Safeguard z(k2) > z(ktp)
             k6   = min (k6, k2 -1)      ! Safeguard z(k6) > z(k2)
             !---------------------------------
             ! Derive "climatological humidity"
             !---------------------------------
             if (strat_q_mode == 3 .OR. strat_q_mode == 4) then
                w = cos (min (abs (cenlat), 45._wp) * 2*d2r) ** 2  ! weight(lat)
                do jl = 1, ktp
                   zt = zf(jl)
                   !------------
                   ! Near poles:
                   !------------
                   if      (zt < zlo_clq_pol) then
                      cl_q_pol(jl) = strat_q_low
                   else if (zt > zhi_clq_pol) then
                      cl_q_pol(jl) = strat_q
                   else
                      cl_q_pol(jl) = ( (zt - zlo_clq_pol) * strat_q   +  &
                                       (zhi_clq_pol - zt) * strat_q_low  &
                                     ) / (zhi_clq_pol - zlo_clq_pol)
                   end if
                   !---------
                   ! Tropics:
                   !---------
                   if      (zt < zlo_clq_tr) then
                      cl_q_tr(jl) = strat_q_low_tr
                   else if (zt > zhi_clq_tr) then
                      cl_q_tr(jl) = strat_q
                   else
                      cl_q_tr(jl) = ( (zt - zlo_clq_tr) * strat_q      +  &
                                      (zhi_clq_tr - zt) * strat_q_low_tr  &
                                    ) / (zhi_clq_tr - zlo_clq_tr)
                   end if
                   !---------------------
                   ! Derive local profile
                   !---------------------
                   cl_q(jl) = (1._wp-w) * cl_q_pol(jl) + w * cl_q_tr(jl)
                end do
             else ! strat_q_mode == 5
                !--------------------------------
                ! Local profile (pressure levels)
                !--------------------------------
                w = (cenlat - hum_clim% lat(1)) / hum_clim% dlat + 1._wp
                k = int (w)
                w = w - k
                if (k < 1) then
                   k = 1
                   w = 0._wp
                else if (k >= hum_clim% nlat) then
                   k = hum_clim% nlat - 1
                   w = 1._wp
                end if
                cl_q_tmp(:) = hum_clim% qv(:,k  ) * (1._wp - w) &
                            + hum_clim% qv(:,k+1) *          w
                !-----------------------------
                ! Local profile (model levels)
                !-----------------------------
                nl = hum_clim% nlev
                do jl = 1, ktp
                   lnp = log (zpf(jl+1))
                   if      (lnp <= hum_clim% lnp(1) ) then
                      k = 1
                      w = 0._wp
                   else if (lnp >= hum_clim% lnp(nl)) then
                      k = nl - 1
                      w = 1._wp
                   else
                      do k = 1, nl-1
                         if (lnp < hum_clim% lnp(k+1)) then
                            w = (          lnp      - hum_clim% lnp(k)) &
                              / (hum_clim% lnp(k+1) - hum_clim% lnp(k))
                            exit
                         end if
                      end do
                      if (k >= nl) call finish ("set_strat_q","fatal: k>=nl")
                   end if
                   cl_q(jl) = (1._wp - w) * cl_q_tmp(k) + w * cl_q_tmp(k+1)
                end do
             end if
             !----------------------------------
             ! Check against saturation over ice
             !----------------------------------
             do jl = 1, ktp
                zqsat    = q_rhi (1._wp, tana(jl), zpf(jl+1))
                cl_q(jl) = min (cl_q(jl), zqsat * qsat_ice_scal)
             end do
             !--------------------------
             ! Relaxation for z > tp+6km
             !--------------------------
             alpha = 1 / (tau_strat * 86400._wp)
             w     = exp (-alpha * dt)
             excl_lat  = excl_reg .and. cenlat >= exclude_lat(1) &
                                  .and. cenlat <= exclude_lat(2)
             do jl = 1, k6
                excl_p = excl_lat .and. zpf(jl+1) >= exclude_p(1) &
                                  .and. zpf(jl+1) <= exclude_p(2)
                if (.not. excl_p) then
                   ref_q(jl) = ref_q(jl) * w + cl_q(jl) * (1._wp-w)
                end if
             end do
             !-------------------------------
             ! Relaxation for tp < z < tp+6km
             !-------------------------------
             dtdz  = (tana(k2) - tana(ktp)) / (zf(k2) - zf(ktp))
             alpha = 1 / (tau_strat_tp * 86400._wp)
             alpha = min (alpha,                                            &
                          (max (dtdz + 2.e-3_wp, 0._wp))**2 * tau_strat_scal)
             w     = exp (-alpha * dt)
             qv6   = ref_q(k6)
             do jl = k6+1, ktp
                ref_q(jl) = ref_q(jl) * w + min (ref_q(jl), qv6) * (1._wp-w)
             end do
             !---------------------------------------------
             ! Use qv profile *above* diagnosed tropopause,
             ! keep analyzed qv at tropopause level.
             !---------------------------------------------
             do jl = 1, ktp-1
                !-------------------------------------------
                ! blending between reference and climatology
                ! between pmax_qclim and 0.9*pmax_qclim
                !-------------------------------------------
                w = min (1._wp, max (0._wp, (pmax_qclim - zpf(jl+1)) &
                                            / (0.1_wp*pmax_qclim)   ))
                an_qv(jl) = w * ref_q(jl) + (1._wp - w) * an_qv(jl)
             end do

          case default ! strat_q_mode /= 3,4,5

            !========================================================
            !* 7.83  ascribe constant stratospheric specific humidity
            !========================================================
            DO 783 jl = iltrop, 1, - 1
!             loice = .false.
!             IF (loice) THEN
!               zesat = svp (c1es, c3ies, c4ies, tana (jl), tmelt)
!             ELSE IF (.not.loice) THEN
!               zesat = svp (c1es, c3les, c4les, tana (jl), tmelt)
!             END IF
!             zqsat = rd / rv * zesat / (zpf (jl + 1) + zepsm1 * zesat)

!             zqsat = q_e (es_t (tana (jl)), zpf (jl))

              !-------------------------------------------------------------
              ! derive limits using relative humidity over ice or over water
              !-------------------------------------------------------------
              if (rh_over_ice) then
                 zqsat = q_rhi (1._wp, tana(jl), zpf(jl+1))
              else
                 zqsat = q_rh  (1._wp, tana(jl), zpf(jl+1))
              end if

              if (l_q_ref) then
                 qstrat = ref_q (jl)    ! spec.hum. from reference model layer
              else
                 qstrat = strat_q       ! spec.hum. from namelist
              end if

              IF (jl .eq. iltrop) THEN
                an_qv (jl) = min (an_qv (jl), 0.97_wp * zqsat)
              ELSE IF (jl .eq. iltrop - 1) THEN
                zqm = 0.5_wp * (an_qv(jl) + qstrat)
                an_qv (jl) = min (an_qv (jl), zqm, 0.7_wp * zqsat, qstrat)
              ELSE
                an_qv (jl) = min (qstrat, 0.5_wp * zqsat)
              END IF
  783       END DO

          end select ! strat_q_mode

        end do
      end do
    end do
  end subroutine set_strat_q_0
!------------------------------------------------------------------------------
  subroutine read_nml_hum_ana
  !------------------------
  ! read namelist /HUM_ANA/
  !------------------------
    integer  :: ierr
    !-------------
    ! set defaults
    !-------------
    hum_ana_top    = 25000._wp ! top of humidity analysis
    e_hum_top      = 0.01_wp   ! nominal FG error above hum_ana_top
    rh0_qcl0       = 0.9_wp    ! rel.hum. below which we set qcl to 0
    strat_q        = 4.0e-6_wp ! upper stratospheric humidity
    strat_q_low    = 2.5e-6_wp ! lower stratospheric humidity
    strat_q_low_tr = 2.5e-6_wp ! lower stratospheric humidity, tropics
    strat_q_mode   = 1         ! 0=keep ana, 1=const. value, 2=use fg
    strat_q_ens    = .false.   ! relax stratospheric hum. of ensemble
    rh_over_ice    = .false.   ! use q_rhi instead of q_rhw
    wmo_lrt_diag   = .false.   ! use WMO lapse rate tropopause diag.
    asurf_min      =  2000._wp ! min. TP height above surface [m]
    pt_min         =  7000._wp ! min. TP pressure             [Pa]
    pt_fallback    = 30000._wp ! TP height when search fails  [Pa]
    tau_strat      =    30._wp ! time constant for q relaxation [d]
    tau_strat_tp   =    0.5_wp ! relaxation time near tropopause [d]
    tau_strat_scal = 1 / 3._wp ! scaling factor for temp.gradient term
    pmax_qclim     =100000._wp ! max. pressure to use clim. QV [Pa]
    qsat_ice_scal  =    1.0_wp ! scaling factor for saturation over ice
    exclude_lat    =  -999._wp ! latitude range to exclude from relax.
    exclude_p      =  -999._wp ! pressure range to exclude     [Pa]
    file           = ""        ! humidity climatology file
    !-------------------------------------------------
    ! Read comments on bug compatibility further up...
    !-------------------------------------------------
    zlo_clq_tr     = z30       ! lower ref. height of q clim., tropics
    zlo_clq_pol    = z30       ! lower ref. height of q clim., poles
    zhi_clq_tr     = z60       ! upper ref. height of q clim., tropics
    zhi_clq_pol    = z60       ! upper ref. height of q clim., poles
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write(6,'(a/)') repeat('-',79)
      call position_nml ('HUM_ANA', status=ierr)
      select case (ierr)
      case (POSITIONED)
        write(6,'(a)')        ' Namelist /HUM_ANA/:'
#if defined(__ibm__)
        read (nnml ,nml=HUM_ANA, iostat=ierr)
        if (ierr/=0) &
          call finish ('read_nml_hum_ana','ERROR in namelist /HUM_ANA/')
#else
        read (nnml ,nml=HUM_ANA)
#endif
      case default
        write(6,'(a)')        ' Namelist /HUM_ANA/ is not present'
        write(6,'(a)')        ' using default parameters'
      end select
      write(6,'(a)')
      write(6,'(a,f12.2,a)')  ' hum_ana_top   =', hum_ana_top/100, ' hPa'
      write(6,'(a,es12.2)')   ' e_hum_top     =', e_hum_top
      write(6,'(a,f12.2,a)')  ' rh0_qcl0      =', rh0_qcl0
      write(6,'(a)')
      write(6,'(a,i12)')      ' strat_q_mode  =', strat_q_mode
      write(6,'(a,l12)')      ' strat_q_ens   =', strat_q_ens
      write(6,'(a,es12.2)')   ' strat_q       =', strat_q
      write(6,'(a,es12.2)')   ' strat_q_low   =', strat_q_low
      write(6,'(a,es12.2)')   ' strat_q_low_tr=', strat_q_low_tr
      write(6,'(a,f12.2,a)')  ' tau_strat     =', tau_strat,   ' d'
      write(6,'(a,f12.2,a)')  ' tau_strat_tp  =', tau_strat_tp,' d'
      write(6,'(a,es12.2)')   ' tau_strat_scal=', tau_strat_scal
      write(6,'(a,f12.2,a)')  ' pmax_qclim    =', pmax_qclim/100, ' hPa'
      write(6,'(a,f12.2)')    ' qsat_ice_scal =', qsat_ice_scal
      write(6,'(a,f12.2,a)')  ' zhi_clq_tr    =', zhi_clq_tr,  ' m'
      write(6,'(a,f12.2,a)')  ' zlo_clq_tr    =', zlo_clq_tr,  ' m'
      write(6,'(a,f12.2,a)')  ' zhi_clq_pol   =', zhi_clq_pol, ' m'
      write(6,'(a,f12.2,a)')  ' zlo_clq_pol   =', zlo_clq_pol, ' m'
      write(6,'(a,l12)')      ' rh_over_ice   =', rh_over_ice
      write(6,'(a,l12)')      ' wmo_lrt_diag  =', wmo_lrt_diag
      write(6,'(a,f12.2,a)')  ' asurf_min     =', asurf_min, ' m'
      write(6,'(a,f12.2,a)')  ' pt_min        =', pt_min     /100, ' hPa'
      write(6,'(a,f12.2,a)')  ' pt_fallback   =', pt_fallback/100, ' hPa'
      write(6,'(a,2f12.2)')   ' exclude_lat   =', exclude_lat
      write(6,'(a,2f12.2,a)') ' exclude_p     =', exclude_p / 100, ' hPa'
      write(6,'(a,a)')        ' file          = ',trim (file)
      write(6,'(a)')
    end if
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (hum_ana_top   ,dace% pio)
    call p_bcast (e_hum_top     ,dace% pio)
    call p_bcast (rh0_qcl0      ,dace% pio)
    call p_bcast (strat_q_mode  ,dace% pio)
    call p_bcast (strat_q_ens   ,dace% pio)
    call p_bcast (strat_q       ,dace% pio)
    call p_bcast (strat_q_low   ,dace% pio)
    call p_bcast (strat_q_low_tr,dace% pio)
    call p_bcast (pmax_qclim    ,dace% pio)
    call p_bcast (rh_over_ice   ,dace% pio)
    call p_bcast (wmo_lrt_diag  ,dace% pio)
    call p_bcast (asurf_min     ,dace% pio)
    call p_bcast (pt_min        ,dace% pio)
    call p_bcast (pt_fallback   ,dace% pio)
    call p_bcast (tau_strat     ,dace% pio)
    call p_bcast (tau_strat_tp  ,dace% pio)
    call p_bcast (tau_strat_scal,dace% pio)
    call p_bcast (qsat_ice_scal ,dace% pio)
    call p_bcast (zlo_clq_tr    ,dace% pio)
    call p_bcast (zlo_clq_pol   ,dace% pio)
    call p_bcast (zhi_clq_tr    ,dace% pio)
    call p_bcast (zhi_clq_pol   ,dace% pio)
    call p_bcast (exclude_lat   ,dace% pio)
    call p_bcast (exclude_p     ,dace% pio)
    call p_bcast (file          ,dace% pio)
    ln_hum_ana_top = -huge(ln_hum_ana_top)
    if (hum_ana_top>0._wp) ln_hum_ana_top = log (hum_ana_top)
    if (strat_q_mode < 0 .or. strat_q_mode > 5) &
      call finish ('read_nml_hum_ana','invalid strat_q_mode')
    if (strat_q_mode >= 3 .and. .not. wmo_lrt_diag) &
      call finish ('read_nml_hum_ana','strat_q_mode=3/4 requires wmo_lrt_diag=T')
    if (strat_q_mode >= 3 .and. .not. rh_over_ice) &
      call finish ('read_nml_hum_ana','strat_q_mode=3/4 requires rh_over_ice=T')
    if (strat_q_ens       .and. strat_q_mode < 3) &
      call finish ('read_nml_hum_ana','strat_q_ens=T requires strat_q_mode=3 or 4')
    if (exclude_lat(1) > exclude_lat(2)) &
      call finish ('read_nml_hum_ana','exclude_lat(1) > exclude_lat(2)')
    if (exclude_p  (1) > exclude_p  (2)) &
      call finish ('read_nml_hum_ana','exclude_p(1) > exclude_p(2)')
    if (strat_q_mode == 5 .and. file == "") &
      call finish ('read_nml_hum_ana','strat_q_mode=5 but no climatology file')

    if (strat_q_mode == 5) then
       call read_hum_climate (path_file (data, file), date=ana_time)
    end if
  end subroutine read_nml_hum_ana
!------------------------------------------------------------------------------
  subroutine read_hum_climate (file, date)
    character(*), intent(in) :: file            ! Climatology filename
    type(t_time), intent(in) :: date            ! Analysis date
    !--------------------------
    ! Read humidity climatology
    !--------------------------
    character(32) :: rname = ""   ! netcdf routine called
    character(32) :: title, version, created
!   character(99) :: history
    character(16) :: h2o_units, lat_units, lev_units
    character(32) :: long_name
    integer       :: ncid
    integer       :: dimid_lat, dimid_lev, dimid_time
    integer       :: varid_lat, varid_lev, varid_time, varid_h2o
    integer       :: nlat, nlev, ntime
    integer       :: doy, ydays
    integer       :: m, m1, m2    ! Time indices (month)
    real(wp)      :: w            ! Interpolation weight
    integer,  allocatable :: month(:)
    real(wp), allocatable :: lat(:), lev(:)
    real(sp), allocatable :: h2o(:,:,:)

    !----------------------------------
    ! Read climatology on I/O processor
    !----------------------------------
    if (dace% lpio) then
       print *, "opening file: ", trim (file)
       rname = "nf90_open"
       call chk (nf90_open (file, NF90_NOWRITE, ncid))

       rname = "nf90_get_att(global)"
       call chk (nf90_get_att (ncid, NF90_GLOBAL, 'title',   title))
       print *, "title  : ", trim (title)
       call chk (nf90_get_att (ncid, NF90_GLOBAL, 'version', version))
       print *, "version: ", trim (version)
       call chk (nf90_get_att (ncid, NF90_GLOBAL, 'created', created))
       print *, "created: ", trim (created)
!      call chk (nf90_get_att (ncid, NF90_GLOBAL, 'history', history))
!      print *, "history: ", trim (history)

       rname = "nf90_inq_dimid"
       call chk (nf90_inq_dimid         (ncid, "lat",   dimid_lat))
       call chk (nf90_inq_dimid         (ncid, "level", dimid_lev))
       call chk (nf90_inq_dimid         (ncid, "month", dimid_time))

       rname = "nf90_inquire_dimension"
       call chk (nf90_inquire_dimension (ncid,          dimid_lat,  len=nlat))
       call chk (nf90_inquire_dimension (ncid,          dimid_lev,  len=nlev))
       call chk (nf90_inquire_dimension (ncid,          dimid_time, len=ntime))

       allocate (lat  (nlat))
       allocate (lev  (nlev))
       allocate (month(ntime))
       allocate (h2o  (nlat,nlev,ntime))

       rname = "nf90_inq_varid"
       call chk (nf90_inq_varid (ncid, 'lat',          varid_lat))
       call chk (nf90_inq_varid (ncid, 'level',        varid_lev))
       call chk (nf90_inq_varid (ncid, 'month',        varid_time))
       call chk (nf90_inq_varid (ncid, 'h2o_seasonal', varid_h2o))

       rname = "nf90_get_var"
       call chk (nf90_get_var (ncid, varid_lat,  lat))
       call chk (nf90_get_var (ncid, varid_lev,  lev))
       call chk (nf90_get_var (ncid, varid_time, month))
       call chk (nf90_get_var (ncid, varid_h2o,  h2o))

       rname = "nf90_get_att"
       call chk (nf90_get_att (ncid, varid_lat, 'units',     lat_units))
       call chk (nf90_get_att (ncid, varid_lev, 'units',     lev_units))
       call chk (nf90_get_att (ncid, varid_h2o, 'units',     h2o_units))
       call chk (nf90_get_att (ncid, varid_h2o, 'long_name', long_name))
       print *, "loading: ", trim (long_name)
       print *

       rname = "nf90_close"
       call chk (nf90_close   (ncid))
    end if

    call p_bcast (nlat,  dace% pio)
    call p_bcast (nlev,  dace% pio)
    call p_bcast (ntime, dace% pio)

    if (.not. dace% lpio) then
       allocate (lat  (nlat))
       allocate (lev  (nlev))
       allocate (month(ntime))
       allocate (h2o  (nlat,nlev,ntime))
    end if

    call p_bcast (lat,       dace% pio)
    call p_bcast (lev,       dace% pio)
    call p_bcast (month,     dace% pio)
    call p_bcast (h2o,       dace% pio)
    call p_bcast (lat_units, dace% pio)
    call p_bcast (lev_units, dace% pio)
    call p_bcast (h2o_units, dace% pio)

    !-------------------------------------
    ! Check sanity of data (units, values)
    !-------------------------------------
    if (lat_units /= "degrees_north") &
         call finish ("read_hum_climate","lat units must be 'degrees_north'")
    if (any (lat(1:nlat-1) >= lat(2:nlat))) &
         call finish ("read_hum_climate","lat not increasing")

    if (lev_units /= "hPa") &
         call finish ("read_hum_climate","level units must be 'hPa'")
    if (any (lev <= 0._wp)) &
         call finish ("read_hum_climate","bad levels")
    if (.not. (all (lev(1:nlev-1) < lev(2:nlev))   &
         .or.  all (lev(1:nlev-1) > lev(2:nlev)))) &
         call finish ("read_hum_climate","level not monotonous")
    if (maxval (lev) < 300._wp) &
         call finish ("read_hum_climate","max. level value must reach 300 hPa")

    if (ntime /= 12) &
         call finish ("read_hum_climate","ntime /= 12")
    if (any (month(1:12) /= [ (m,m=1,12) ])) &
         call finish ("read_hum_climate","unexpected ordering of months")

    if (h2o_units /= "ppmv") &
         call finish ("read_hum_climate",                     &
                      "unsupported unit: " // trim (h2o_units))
    if (any (h2o == NF90_FILL_REAL)) &
         call finish ("read_hum_climate","cannot handle missing values")
    if (any (h2o < 1._sp) .or. any (h2o > 1000._sp)) &
         call finish ("read_hum_climate","check range of h2o_seasonal")

    !-----------------------
    ! Temporal interpolation
    !-----------------------
    doy   = mod (iyyddd (date), 1000)   ! Day of year
    ydays = days_year   (date)          ! Days in current year
    w     = 12._wp * (doy + 15) / ydays ! Temporal interpolation weight
    m1    = int (w)                     ! Jan 16 -> m1=1, m2=2
    m2    = m1 + 1
    w     = w - m1
    if (m1 < 1 .or. m2 > 12) then
       m1 = 12
       m2 = 1
    end if

    hum_clim% nlat = nlat
    hum_clim% nlev = nlev
    allocate (hum_clim% lat(nlat))
    allocate (hum_clim% lnp(nlev))
    allocate (hum_clim% qv (nlev,nlat))
    hum_clim% lnp(:)   = log (lev(:) * 100._wp)   ! hPa -> log(Pa)
    hum_clim% lat(:)   =      lat(:)
    hum_clim% dlat     =     (lat(nlat) - lat(1)) / (nlat-1)
    !-----------------------------------------------
    ! Interpolate, convert ppmv to specific humidity
    !-----------------------------------------------
    hum_clim% qv (:,:) = transpose (( h2o(:,:,m1) * (1._wp - w)   &
                                    + h2o(:,:,m2) *          w  ) &
                                    * (RdRv * 1.e-6_wp)           )
    !---------------------------------------------
    ! Flip levels such that pressure is increasing
    !---------------------------------------------
    if (lev(1) > lev(nlev)) then
       hum_clim% lnp(1:nlev  ) = hum_clim% lnp(nlev:1:-1  )
       hum_clim% qv (1:nlev,:) = hum_clim% qv (nlev:1:-1,:)
    end if

    if (dace% lpio) then
      write(6,'(a,2(f8.2,4x,a))') ' levels : ',minval (lev),          " ...", &
                                               maxval (lev),          " hPa"
      write(6,'(a,2(es12.2,a))')  ' qv     : ',minval (hum_clim% qv), " ...", &
                                               maxval (hum_clim% qv), " kg/kg"
      write(6,'(a)')
    end if

  contains

    subroutine chk (status, ierr)
      integer, intent(in)            :: status
      integer, intent(out) ,optional :: ierr
      !-----------------------------------------------------------
      ! check error status after each NetCDF call, abort on error.
      !-----------------------------------------------------------
      if (present (ierr)) then
         ierr = status
      else
         if (status /= NF90_NOERR) then
            call finish ("read_hum_climate", &
                         trim (rname) // " : " // trim (nf90_strerror (status)))
         end if
      endif
    end subroutine chk

  end subroutine read_hum_climate
!------------------------------------------------------------------------------
end module mo_hum_ana
