!
!+ Post multiplication routines
!
!------------------------------------------------------------------------------
! This compiler directive must appear before the MODULE statement !
! The IBM seems to be unable to allocate spill space
! (=stack space for temporary storage of register contents) dynamically.
! Default according to docs: SPILLSIZE(512), but finish_ana needs much more.
#if defined (__ibm__) && (__SIZE_PTR==64)
@PROCESS SPILLSIZE(2048)
#endif
!------------------------------------------------------------------------------
MODULE mo_postmult
!
! Description:
!   Post multiplication routines.
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
! V1_4         2009/03/26 Hendrik Reich
!  Changes for COSMO: rotated coordinates, C-grid, vertical coordinate
! V1_5         2009/05/25 Harald Anlauf
!  use vectorized tq_tvgh_vec instead of elemental tq_tvgh
! V1_7         2009/08/24 Andreas Rhodin
!  bug fix (CHECK!) for COSMO ensemble generation
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  finish_ana: new optional parameter lanagp (derive tv from geoh)
! V1_13        2011/11/01 Andreas Rhodin
!  bugfix for LETKF ensemble generation
! V1_20        2012-06-18 Andreas Rhodin
!  rename mo_profiles to mo_profile (make module name and file name consistent)
! V1_22        2013-02-13 Harald Anlauf
!  finish_ana: simple hydrostatic balancing of pressure increments
! V1_26        2013/06/27 Harald Anlauf
!  Modify workaround for bug for sxf90 rev.451
! V1_27        2013-11-08 Harald Anlauf
!  for ICON: re-diagnosis of ph, implement strat_q_mode=2
! V1_28        2014/02/26 Harald Anlauf
!  finish_ana: improve stability of interpolation of analysis for ICON
!  Fix for rev.9619 (due to attempt to take water loading into account)
! V1_35        2014-11-07 Harald Anlauf
!  WMO tropopause diagnostics, revised treatment of stratospheric qv
!    in analysis (strat_q_mode=3)
! V1_37        2014-12-23 Harald Anlauf
!  finish_ana: option to keep surface geop. fixed for ICON analysis
! V1_42        2015-06-08 Harald Anlauf
!  finish_ana: improve fix_geosp=T by relaxation near model top
! V1_43        2015-08-19 Harald Anlauf
!  namelist /POSTMULT/, finish_ana: T-inversion preserving interpolation
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_48        2016-10-06 Harald Anlauf
!  wind adjustment for analysis increment
!  implement 'lanagp' for ICON EnVar (A.Rhodin)
! V1_50        2017-01-09 Andreas Rhodin
!  for COSMO: handle vertical coordinates like ICON
! V1_51        2017-02-24 Valerio Cardinali
!  adjust_wind: print out some more diagnostics on 'dpsdt' (for COMET)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006-2008  post_mult moved from module mo_psas
! Harald Anlauf   DWD  2007-2008  work around IBM and NEC/SX8 problems
!=========================================================================

  !-------------
  ! modules used
  !-------------
  use mo_kind,        only: wp                 ! kind parameter
  use mo_mpi_dace,    only: dace,             &! MPI group info
                            p_barrier,        &! generic MPI barrier routine
                            p_min,            &! generic min-function
                            p_max,            &! generic max-function
                            p_sum,            &! generic MPI sum
                            p_bcast            ! generic broadcast routine
  use mo_exception,   only: finish,           &! exit in case of error
                            message            ! write message
  use mo_cpu_time,    only: stop_time          ! determine cpu and wall time
  use mo_namelist,    only: position_nml,     &! position namelist
                            nnml,             &! namelist Fortran unit number
                            POSITIONED         ! ok    code from position_nml
  use mo_t_obs,       only: t_obs,            &! observation data type
                            t_spot,           &!   component of t_obs
                            OBS_TV, OBS_RH,   &! observation type id.s
                            OBS_H,            &! observation type id.s
                            OBS_U, OBS_V       ! observation type id.s
  use mo_obs_set,     only: t_obs_set          ! observation data type
  use mo_dec_matrix,  only: t_vector,         &! decomposed vector
                            t_vector_segm      !   component of t_vector
  use mo_t_col,       only: t_cols,           &! data type to hold columns
                            t_col,            &! hold one column of the model
                            alloc_cols,       &! allocate components of t_cols
                            dealloc_cols,     &! deallocate comp.    of t_cols
                            assignment(=),    &! type(t_cols) = real(wp)
                            p_send,           &! generic MPI send routine
                            p_recv,           &! generic MPI receive routine
                            alltoall_cols,    &! redistribute columns
                            set_xuv,          &! set cartes.coord. from lon/lat
                            COL_T, COL_UV,    &! identifier for temp., hor.wind,
                            COL_RH,           &!   rel.humidity
                            COL_P              !   pressure
  use mo_fg_cov,      only: t_rowcol,         &! type for precalc. quantities
                            x_cut,            &! cutoff length scale
                            L_max,            &! max.hor.lgth.scale(unitsphere)
                            init_spot,        &! precalc. quantities for spot
                            get_corr,         &! get correl. for pair of spots
                            construct,        &! construct var.of type t_rowcol
                            destruct           ! destruct var. of type t_rowcol
  use mo_physics,     only: gacc,             &! gravity acceleration
                            Rd => R,          &! gas constant of dry air
                            RDDRM1,           &! Rv/Rd - 1._wp
                            pi,               &! 3.1415...
                            rh_q,             &! derive rel.hum. from spec.h.
 !                          tv_t_q,           &! derive Tv from T,q
                            t_tv_q             ! derive temperature from Tv
  use mo_t_bg_err_op, only: covm               ! NMC fitted covariances
  use mo_bg_err_op,   only: repr_2dh,         &! hor.covariance represent. flag
                            n_av_pole          ! average anal.incr. at poles
  use mo_bg_err_2d,   only: apply_B_mi         ! apply_B_mi
  use mo_atm_state,   only: t_atm,            &! atmospheric state derived type
                            operator(+),      &! add atmospheric fields
                            operator(-),      &! subtract atmospheric fields
                            operator(*),      &! scale atmospheric fields
                            assignment(=),    &! generic assignment
                            set_geo,          &! set geopotential
                            set_p,            &! set pressure
                            set_pf,           &! set pressure at full levels
!                           set_ph,           &! set pressure at half levels
                            set_ps,           &! compute surface pressure
                            set_tv_geo,       &! virtual temperat. from height
                            select_params,    &! gather fields from states
                            vert_intp,        &! vertical (spline)interpolation
                            p_sum_atm,        &! sum over processors
                            construct,        &! construct atmospheric state
                            destruct,         &! free atmospheric state
                            allocate,         &! allocate atmospheric fields
                            deallocate,       &! deallocate atmospheric fields
                            to_grads,         &! write t_atm%m(:) to grads file
                            print,            &! printout
                            A_to_C_grid,      &! transform from A-grid to C-grid
                            average_pole_wind,&! average poles from neighbours
                            average_pole_scalar,&! average poles from neighbours
                            uvtrafo,          &! transform of (u,v) geo<-->rot
                            mean_stdev         ! calculate mean and stdev
  use mo_wmo_tables,  only: WMO3_ISOBARIC,    &! level type
                            WMO6_GAUSSIAN,    &! grid  type
                            WMO6_LATLON,      &! latitude/longitude grid
                            WMO6_ROTLL,       &! dto., rotated
                            DWD6_ICOSAHEDRON, &! GME  triangular
                            DWD6_ICON          ! ICON triangular
  use mo_atm_grid,    only: t_grid,           &! atmospheric grid derived type
                            init_ctl,         &! preset .ctl file from grid
                            to_grads,         &! write t_grid to grads file
                            VCT_P_ISO,        &!       isobaric p coordinate
                            VCT_P_HYB,        &! GME/HRM hybrid p coordinate
                            VCT_Z_HYB,        &!         hybrid z coordinate
                            VCT_Z_GEN,        &!    generalised z coordinate
                            MO_ICON, MO_GME,  &! ICON / GME / IFS  model
                            MO_IFS, MO_COSMO
  use mo_profile,     only: write_profiles     ! write atmospheric profiles
  use mo_hum_ana,     only: hum_ana_incr,     &! set analysis increment to zero
                            set_strat_q,      &! set stratospheric spec.hum.
                            rh0_qcl0           ! rel.hum. threshold for qcl=0
  use mo_run_params,  only: aux,              &! path for additional output
                            path_file,        &! add a path to a file name
                            pass_fields        ! pass some fields to the output
  use mo_grads,       only: t_ctl,            &! grads  .ctl file data type
                            write_ctl,        &! write  .ctl file
                            destruct           ! clean up .ctl data type
  use mo_gmt,         only: atm_to_gmt         ! write GMT plotfile
  use mo_cntrlvar,    only: gh_rh,            &! generalized humid. <- rel.hum.
!                           tq_tvgh,          &! t,spc.hum. from vrt.t,gen.hum.
                            tq_tvgh_vec,      &! t,spc.hum. from vrt.t,gen.hum.
                            rh0                ! <=0 : no generalized humidity
  use meteo_utilities,only: calrho             ! compute air density
  use mo_grib,        only: write_grib         ! write atmospheric state
  use mo_math_oper,   only: math_oper_init,   &! Initialized math operators
                            math_oper_cleanup,&! Clean up module
                            hor_div,          &! Divergence of vector field
                            hor_div_t,        &! Dto., transposed version
                            hor_curl,         &! Curl of vector field
                            hor_curl_t         ! Dto., transposed version
  use mo_dace_string, only: tolower,          &! Lowercase of string
                            char3              ! conversion: integer -> char(3)
  use mo_atm_util,    only: prep_extana,      &! prepare external analysis
                            read_extana        ! parameters from external analysis
  use mo_synop,       only: version            ! observation operator version
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: post_mult         ! perform post multiplication
  public :: setup_post        ! set up model columns (t_cols) for post multipl.
  public :: finish_post       ! finish post multiplication: copy t_cols to t_atm
  public :: finish_ana        ! finish analysis
  public :: read_nml_postmult ! read namelist /POSTMULT/
  public :: adjust_wind       ! adjust wind
  public :: adjust_wind_ens   ! adjust wind (ensemble version)
  public :: wind_adjust       ! wind adjustment mode

  !-----------------
  ! Module variables
  !-----------------
  real(wp), parameter :: h0         =  36000._wp ! Start level of height relaxation
  real(wp), parameter :: dtdz_min   = -8.5e-3_wp ! min. vertical gradient [K/m]
  real(wp), parameter :: dtdz_max   = -1.5e-3_wp ! max. vertical gradient [K/m]

  !===================
  ! namelist variables
  !===================
  logical             :: lkeep_tinv = .false.    ! Preserve temp.inversions
  logical             :: lfix_geosp = .false.    ! Debugging (do not use!)
  real(wp)            :: zpbl1      =   500._wp  ! AGL heights used in
  real(wp)            :: zpbl2      =  1000._wp  ! lapse rate computations
  !-------------------------
  ! Analysis wind adjustment
  !-------------------------
  integer, protected  :: wind_adjust = 0         ! Wind adjustment mode:
                                                 ! 0 = off
                                                 ! 1 = diagnose dpsdt of states
                                                 ! 2 = derive adjustment(u,v)
                                                 ! 3 = apply  adjustment(u,v)
  integer             :: diag_level  = 1         ! Level of diagnostics
  logical             :: adjust_det  = .true.    ! Adjust deterministic analysis
  logical             :: adjust_ens  = .false.   ! Adjust ensemble      members
  integer             :: adjust_mode = 0         ! dpsdt source for adjustment:
                                                 ! 0 = dpsdt from LETKF fg
                                                 ! 1 = dpsdt from LETKF ana
  real(wp)            :: w_spread_an = 1._wp     ! Rel. weight ana/fg spread
  real(wp)            :: sqz_trans(2)= 0._wp     ! Mod. for squeezed levels
  character(len=8)    :: solver      = "cgls"    ! Solver type (name)
  integer             :: min_iter    = 1         ! Min. no. iterations
  integer             :: max_iter    = 20        ! Max. no. iterations
  real(wp)            :: tol(2)      = 0._wp     ! Convergence criteria
  character(len=128)  :: diag_base   = "diag"    ! Basename of diagnostics file

contains
!==============================================================================

  subroutine post_mult  (a_m, cbg, obs, z_i, e_fi, lhskep, lnewpl, pert_B)
  type(t_cols)    ,intent(out)   :: a_m(:) ! analysis increment (model space)
  type(t_cols)    ,intent(in)    :: cbg(:) ! reference state
  type(t_obs_set) ,intent(inout) :: obs    ! observation data type array
  type(t_vector)  ,intent(in)    :: z_i    ! analysis increment (intp.space)
  type(t_vector)  ,intent(in)    :: e_fi   ! background error   (intp.space)
  integer         ,intent(in)    :: lhskep ! keep lhs in post-multiplication
  logical         ,intent(in)    :: lnewpl ! analysis increment on new p-levs
  real(wp)        ,intent(in)    :: pert_B ! random perturbation factor

    if (covm% valid >= 4) then
      call apply_B_mi   (a_m, cbg, obs,    z_i,         lnewpl, e_fi, pert_B)
    else
      call post_mult_ex (a_m, cbg, obs% o, z_i, lhskep, lnewpl)
    endif

  end subroutine post_mult

!------------------------------------------------------------------------------

  subroutine post_mult_ex  (a_m, cbg, obs, z_i, lhskep, lnewpl)
  type(t_cols)        ,intent(out) :: a_m(:) ! analysis increment (model space)
  type(t_cols)        ,intent(in)  :: cbg(:) ! reference state
  type(t_obs) ,target ,intent(in)  :: obs(:) ! observation data type array
  type(t_vector)      ,intent(in)  :: z_i    ! analysis increment (int. space)
  integer             ,intent(in)  :: lhskep ! keep lhs in post-multiplication
  logical             ,intent(in)  :: lnewpl ! analysis increment on new p-levs

    !----------------
    ! local variables
    !----------------
    integer          ,parameter :: nm = 4 ! number of multilevel fields to ass.
    integer          ,parameter :: ns = 1 ! number of singlelevel fields
    integer                     :: l      ! background box index
    integer                     :: nb     ! number of background boxes
    integer                     :: k      ! index
    integer                     :: ic     ! grid column index
    integer                     :: js     ! obsercation spot index
    integer                     :: is     ! lhs spot index
    integer                     :: jb     ! observationboxes/vectorsegment indx
    type(t_vector_segm),pointer :: zi     ! vector segment pointer
    type(t_obs)        ,pointer :: ob     ! observation box pointer
    type(t_spot)       ,pointer :: s      ! observation spot pointer
    type(t_col)        ,pointer :: am     ! analysis   column pointer
    type(t_col)        ,pointer :: bg     ! background column pointer
    real(wp)                    :: ds     ! secant: grid column - obsercation
    integer                     :: t_int (cbg(1)%ke*nm+ns) ! obs.ids.
    real(wp)                    :: lev   (cbg(1)%ke*nm+ns) ! ln(p)
    real(wp)                    :: elev  (cbg(1)%ke   +1 ) ! temporary
    real(wp)                    :: ac    (cbg(1)%ke*nm+ns) ! analysis increm.
    integer                     :: ke     ! number of levels
    integer                     :: nvc    ! number of variables / column
    integer                     :: nonz   ! number of nonzero elements/matrix
    type(t_rowcol)              :: lhs    ! description of left hand side
    type(t_rowcol)              :: rhs    ! desgription of right hand side
!   integer                     :: nr     ! number of rows
    real(wp)                    :: dh     ! geopotential increment
    integer                     :: m      ! outer loop index variable
    integer                     :: m1, m2 ! outer loop index bounds
    integer                     :: lhsk   ! local copy of lhskep
    !++++++++++++++++
    ! temporary flags
    !++++++++++++++++
    lhsk = lhskep

    !--------------------------------------------
    ! outer loop: m = 1 : set pressure level only
    !             m = 2 : set all variables
    !--------------------------------------------
    m2 = 2;     m1 = 2
    if (lnewpl) m1 = 1
    do m = m1, m2
      !------------------------------
      ! set identifier for parameters
      !------------------------------
      ke = cbg(1)% ke
      t_int (1)     = OBS_H
      if (m==1) then
        nvc = 1
      else
        nvc = nm * ke + ns
        t_int (2::nm) = OBS_TV
        t_int (3::nm) = OBS_RH
        t_int (4::nm) = OBS_U
        t_int (5::nm) = OBS_V
      endif
      !-------------
      ! allocate lhs
      !-------------
      if(lhsk == 0) then
        call construct (lhs,   &! row or column to set
                        nvc,   &! number of control parameters
                        1)      ! number of spots
      endif
      !------------------------
      ! loop over grid segments
      !------------------------
      nb = size(cbg)
      if (dace% lpio) write(6,*) 'post_mult: nb',nb
      do l = 1, nb
        !+++++++++++++++++++++++++++++++++++++++++++++++
        ! some communication in this part of the program
        !+++++++++++++++++++++++++++++++++++++++++++++++
        call p_barrier
        if (dace% lpio) write(6,*) 'post_mult: gridbox,ncol',l,cbg(l)% ncol
        call p_barrier
        !--------------------------------------
        ! allocate result variable, set to zero
        !--------------------------------------
        if (cbg(l)% ncol == 0) cycle
        if (m==m1) then
          call alloc_cols (a_m(l), tmp=cbg(l), ids=COL_T+COL_UV+COL_RH)
          a_m(l) = 0._wp
        endif
        !-------------
        ! allocate lhs
        !-------------
        if(lhsk /= 0) then
          call construct (lhs,              &! row or column to set
                          a_m(l)% ncol*nvc, &! number of control parameters
                          a_m(l)% ncol)      ! number of spots
        endif
        !----------------------------------------------------------------
        ! loop over observation boxes / vector segments (right hand side)
        !----------------------------------------------------------------
        do jb = 1, z_i% n_s
          zi  => z_i% s(jb)
          ob  => obs   (jb)
          !----------
          ! setup rhs
          !----------
          call construct (rhs,        &!
                          ob% n_int,  &! number of columns
                          ob% n_spot)  ! number of spots in a column
          do js = 1, ob% n_spot
            s => ob% spot(js)
            !------------------
            ! set rhs structure
            !------------------
            call init_spot (rhs,             &!
                            js,              &! spot index within the box
                            s% col% c,       &! coordinates
                            ob% t_int(s%i%i+1:s%i%i+s%i%n) ,&!
                            ob% lev  (s%i%i+1:s%i%i+s%i%n))
          end do
          !----------------------------------------
          ! loop over grid columns (left hand side)
          !----------------------------------------
          do ic = 1, a_m(l)% ncol
            am => a_m(l)% col(ic)
            bg => cbg(l)% col(ic)
            !------------------
            ! setup lhs: levels
            !------------------
            if (lhsk == 0) then
              is = 1
              lhs% spots% ih = 0
              lhs% in        = 0
            else
              is = ic
            endif
            if (lhsk == 0 .or. jb == 1) then
              !--------------------
              ! set pressure levels
              !--------------------
              lev   (1)               = log (bg% s% ps)
              select case (cbg(l)% vctype)
              case (VCT_P_ISO, VCT_P_HYB)
                do k=1,ke+1
                  elev (k) = a_m(l)% ak(k) + a_m(l)% bk(k) * bg% s% ps
                end do
                do k=1,ke
                  lev (nm*(k-1)+2:nm*k+1) = log(0.5_wp * (elev (k) + elev (k+1)))
                end do
              case (VCT_Z_HYB, VCT_Z_GEN)
                do k=1,ke
                  lev (nm*(k-1)+2:nm*k+1) = log (bg% p (k))
                end do
              case default
                call finish ('post_mult_ex','invalid vctype')
              end select
              !----------------
              ! setup lhs: spot
              !----------------
              call init_spot (lhs,        &!
                              is,         &! spot index within the box
                              am% c,      &! coordinates
                              t_int(:nvc),&!
                              lev  (:nvc))
            endif
            !----------------------------------
            ! loop over spots (right hand side)
            !----------------------------------
            dh = 0._wp
            do js = 1, ob% n_spot
              s => ob% spot (js)
!             nr = s%i% n
              !-------------------------------------------------------
              ! calculate distance column-spot (secant of unit sphere)
              !-------------------------------------------------------
              ds = sqrt( sum( (am%c% x - s% col% c% x)**2 ) )

              if (ds > x_cut*l_max) cycle
              !----------------------------------
              ! set matrix for interpolated space
              !----------------------------------
              nonz = 0
              ac (:nvc) = 0._wp
              call get_corr (lhs, rhs, is, js, .TRUE., .false., nonz, &
                             l=ac(:nvc), r=zi% x(s%i%i+1:s%i%i+s%i%n) )
              !----------------------
              ! matrix multiplication
              !----------------------
              if (nonz>0) then
                if (m==m1) then
                  dh      =     dh    + ac (1)
                endif
                if (m==m2) then
                  am%  t  = am% t     + ac (2:nm*ke+1:nm)
                  am%  rh = am% rh    + ac (3:nm*ke+1:nm)
                  am%  u  = am% u     + ac (4:nm*ke+1:nm)
                  am%  v  = am% v     + ac (5:nm*ke+1:nm)
                endif
              endif
              !------------------
              ! end do's, cleanup
              !------------------
            end do ! js = 1, ob% n_spot   ! loop over observation spots
            !--------------------------------------------------
            ! convert geopotential increment to pressure change
            !--------------------------------------------------
            am% s% ps = am% s% ps + dh * (gacc / Rd) * (bg% s% ps / bg% s% t2m)
          end do   ! ic = 1, a_m(l)% ncol ! loop over grid columns
          call destruct (rhs)
        end do     ! jb = 1, z_i% n_s     ! loop over observation boxes
        if (lhsk /= 0) call destruct (lhs)

        !-------------------------------------------------
        ! special handling for pressure if 'lnewpl' is set
        !-------------------------------------------------
        if (lnewpl) then
          do ic = 1, a_m(l)% ncol ! loop over grid columns
            am => a_m(l)% col(ic)
            bg => cbg(l)% col(ic)
            if (m==m1) bg% s% ps = bg% s% ps + am% s% ps
          end do
        end if

      end do       ! l  = 1, nb           ! loop over grid segments
      if (lhsk == 0) call destruct (lhs)

    end do ! m = m1, m2
  end subroutine post_mult_ex
!==============================================================================
  subroutine setup_post (cols, bg, pes_post)
  type (t_cols) ,pointer     :: cols(:)  ! model columns
  type (t_atm)  ,intent(in)  :: bg       ! background atmospheric state
  integer       ,intent(in)  :: pes_post ! stride for Pes in post-multiplic.

    type (t_grid) ,pointer :: g                  ! grid pointer
    logical                :: ldec               ! analysis distrib. over PE's
    integer                :: npe                ! number of source PE's
!   integer                :: gridtype           !
    integer                :: l
    integer                :: pe
    type (t_cols)          :: col

    g => bg% grid
    ldec     = g% ldec
    npe      = g% dc% npe
!   gridtype = g% gridtype
    !--------------------------------------------------------
    ! allocate columns (currently only surface pressure, t2m)
    !--------------------------------------------------------
    allocate (cols (npe))
    if (.not.ldec) then
      call spread_bg (cols(1), bg, dace% pe, pes_post)
    else
      do pe = 0, npe-1
        l = pe+1
        if (dace% pe < pe) then
          call spread_bg (col, bg, pe, pes_post)
          call p_send (col     ,pe)
          call p_recv (cols(l),pe)
          call dealloc_cols (col)
        else if (dace% pe == pe) then
          call spread_bg (cols(l), bg, pe, pes_post)
        else if (dace% pe > pe) then
          call spread_bg (col, bg, pe, pes_post)
          call p_recv (cols(l),pe)
          call p_send (col     ,pe)
          call dealloc_cols (col)
        endif
      end do
    endif

  end subroutine setup_post

!------------------------------------------------------------------------------

  subroutine spread_bg (cols, bg, pe, pes_post)
  type (t_cols) ,intent(out) :: cols     ! model columns
  type (t_atm)  ,intent(in)  :: bg       ! background atmospheric state
  integer       ,intent(in)  :: pe       ! pe to send data to
  integer       ,intent(in)  :: pes_post ! stride for Pes in post-multiplic.
  !-------------------------------------------------
  ! Scatter background over processors (round robin)
  !-------------------------------------------------
    !----------------
    ! local variables
    !----------------
    real(wp) ,parameter    :: rdp = 180._wp / pi ! conversion factor
    integer                :: ncol               ! number of grid columns
    integer                :: i,j,id,l,m         ! indices
    type (t_grid) ,pointer :: g                  ! grid pointer
    type (t_col)  ,pointer :: c                  ! grid column pointer
    integer                :: nx,ny,nd           ! local grid size
    integer                :: nprocs             ! actual number of PEs used
    integer                :: pes                ! pe / pes_post
    !------------------------------------------------------------
    ! determine number of processors holding the background state
    !------------------------------------------------------------
    nprocs = dace% npe / pes_post
    pes = pe / pes_post
    g => bg% grid
    nx       = g% ub(1) - g% lb(1) + 1
    ny       = g% ub(2) - g% lb(2) + 1
    nd       = g% ub(4) - g% lb(4) + 1
    if (mod(pe,pes_post)==0) then
      ncol = ( (nx * ny * nd) + (nprocs-1) - pes ) / nprocs
    else
      ncol = 0
    endif
    !--------------------------------------------------------
    ! allocate columns (currently only surface pressure, t2m)
    !--------------------------------------------------------
    !-----------------------------
    ! if PE is used spread columns
    !-----------------------------
    select case (g% vct)
    case (VCT_P_ISO)
      !---------------------
      ! pressure coordinates
      !---------------------
      call alloc_cols (cols, ncol, g% nz, ntr=0)
      cols% ak   = g% akf (:g% nz+1)
      cols% bk   = 0._wp
    case (VCT_P_HYB)
      !------------------------------------
      ! GME/HRM hybrid pressure coordinates
      !------------------------------------
      call alloc_cols (cols, ncol, g% nz, ntr=0)
      cols% ak   = g% ak  (:g% nz+1)
      cols% bk   = g% bk  (:g% nz+1)
    case default
      !--------------------------
      ! COSMO hybrid z-coordinate
      !--------------------------
      if (.not.associated(bg% pf))                     &
        call finish('spread_bg','bg% pf is not present')
      call alloc_cols (cols, ncol, g% nz, ids=COL_P)
    end select
    cols% ke      = g%     nz
    cols% levtyp  = g%     levtyp
    cols% vctype  = g% vct

    if (ncol /= 0) then
      m = 0
      l = 0
      do id = g%lb(4), g%ub(4)
        do j = g%lb(2), g%ub(2)
          do i = g%lb(1), g%ub(1)
            if (mod(m,nprocs)==pes) then
              l = l + 1
              c => cols% col(l)
              c% c% dlat = rdp * g% rlat (i,j,1,id)
              c% c% dlon = rdp * g% rlon (i,j,1,id)
              c%    i    = i
              c%    j    = j
              c%    l    = id
              c%    pe   = dace% pe
              call set_xuv (c)
              c%   s% ps  = 0._wp
              c%   s% t2m = 0._wp    ! put t(nz) into t2m !!
              if(associated(bg% t )) c% s% t2m = bg% t  (i,j,g%nz,id)
              if(associated(bg% ps)) c% s% ps  = bg% ps (i,j,   1,id)
              if(c%s% t2m == 0._wp ) c% s% t2m =    300._wp
              if(c%s% ps  == 0._wp ) c% s% ps  = 100000._wp
              if(g%vct/=VCT_P_HYB  ) c%    p   = bg% pf (i,j,   :,id)
!            cols% col(i)% x   = tmp% col(i)% x
!            cols% col(i)% du  = tmp% col(i)% du
!            cols% col(i)% dv  = tmp% col(i)% dv
!            cols% col(i)% lat = tmp% col(i)% lat
!            cols% col(i)% lon = tmp% col(i)% lon
!            cols% col(i)% i   = tmp% col(i)% i
!            cols% col(i)% j   = tmp% col(i)% j
!            cols% col(i)% l   = tmp% col(i)% l
!            cols% col(i)% pe  = tmp% col(i)% pe

            endif
            m = m + 1
          end do
        end do
      end do
      if (l/=ncol) call finish('spread_bg','l/=ncol')
    endif
  end subroutine spread_bg

!==============================================================================

  subroutine finish_post (ai, cols, bg)
  type (t_atm)  ,intent(out) :: ai      ! analysis increment
  type (t_cols) ,intent(in)  :: cols(:) ! model columns
  type (t_atm)  ,intent(in)  :: bg      ! background atmospheric state
  !------------------------------------------------------------------------
  ! finish post multiplication: copy data from derived type t_cols to t_atm
  !------------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    integer                :: l
    integer                :: pe
    logical                :: ldec     ! analysis distrib. over PE's
    type (t_cols)          :: col
    type (t_cols)          :: temp(dace% npe)

    ldec     = bg% grid% ldec
    !-------------------------------------
    ! allocate analysis increment variable
    !-------------------------------------
    if (covm% valid >= 4 .and. repr_2dh == 2) then
      select case (bg% grid% vct)
      case (VCT_P_ISO)
        !--------------------
        ! pressure coordinate
        !--------------------
        call construct (ai, template = bg, alloc = 't rh u v ps geof')
      case (VCT_P_HYB)
        !-----------------------------------
        ! GME/HRM hybrid pressure coordinate
        !-----------------------------------
        call construct (ai, template = bg, alloc = 't rh u v ps geoh')
      case default
        !-------------------------------
        ! COSMO/ICON hybrid z-coordinate
        !-------------------------------
        call construct (ai, template = bg, alloc = 't rh u v ps geof')
      end select
    else
      call construct (ai, template = bg, alloc = 't rh u v ps')
    endif
    ai = 0._wp
    if (.not.ldec) then
      call col2atm (ai, cols(1))
      !---------------------------------------
      ! quick and dirty spread over processors
      !---------------------------------------
      call p_sum_atm (ai)
#ifndef NO_OPT_MPI_ATOA
    else if (size (cols) == dace% npe) then
      !-------------------
      ! Optimized alltoall
      !-------------------
      call alltoall_cols (cols, temp)
      do l = 1, dace% npe
        pe = l-1
        if (dace% pe == pe) then
          call col2atm (ai,  cols(l))
        else ! if (pe /= dace% pe) then
          call col2atm (ai,  temp(l))
          call dealloc_cols (temp(l))
        end if
      end do
#endif
    else
      do l = 1, size(cols)
        pe = l-1
        if (dace% pe < pe) then
          call p_send (cols(l),pe)
          call p_recv (col    ,pe)
          call col2atm (ai, col)
          call dealloc_cols (col)
        else if (dace% pe == pe) then
          call col2atm (ai, cols(l))
        else if (dace% pe > pe) then
          call p_recv (col    ,pe)
          call p_send (cols(l),pe)
          call col2atm (ai, col)
          call dealloc_cols (col)
        endif
      end do
    endif
  end subroutine finish_post

!------------------------------------------------------------------------------

  subroutine col2atm (ai, cols)
  type (t_atm)  ,intent(inout) :: ai   ! analysis increment
  type (t_cols) ,intent(in)    :: cols ! model columns
    !----------------
    ! local variables
    !----------------
    integer                :: i,j,id,l ! indices
    type (t_col)  ,pointer :: c        ! grid column pointer
    !----------
    ! copy data
    !----------
    do l = 1, cols% ncol
      c => cols% col(l); i = c%i; j = c%j; id = c%l
      ai%       t     (i,j,:,id) = c%    t
      ai%       rh    (i,j,:,id) = c%    rh
      ai%       u     (i,j,:,id) = c%    u
      ai%       v     (i,j,:,id) = c%    v
!     ai% grid% geosp (i,j,:,id) = c% s% geosp
      ai%       ps    (i,j,:,id) = c% s% ps
      if (associated(ai% geoh)) &
      ai%       geoh  (i,j,:,id) = c%    geoh
      if (associated(ai% geof)) &
      ai%       geof  (i,j,:,id) = c%    geoh (2:)
    end do
  end subroutine col2atm

!==============================================================================

  subroutine finish_ana (ana, ani, fc, lvint, lanagp, lerror_tq, &
                         lplot, verb, a2c, lextana               )
  !-------------------------------------------------------------------------
  ! finish analysis:
  !
  !   1) average poles from neighbours (on icosahedron)
  !   2) take vertical coordinates from forecast
  !   3) optionally write profiles for diagnostics
  !   4) deallocate fields not used any more
  !   5) special handling for stratospheric humidity
  !   6) optionally write plot file for diagnostics
  !   7) set analysis time (ana and fc)
  !   8) print out analysis increment
  !   9) restore C-grid
  !  10) allocate analysis state
  !  11) temporarily derive vt, gh in forecast; save t, rh
  !  12) analysis = analysis increment + forecast
  !  13) force negative RH values to zero
  !  14) derive spec.humid. and temperature from gen.humid. and virtual temp.
  !  15) special handling for stratospheric humidity
  !  16) restore temperature, relative humidity
  !  17) copy some quantities from the forecast
  !  18) zero increments for quantities not yet analysed; no qcl for rh < 0.9
  !  19) restore temperature and rel.humidity in forecast
  !  20) re-derive ps from p (COSMO)
  !--------------------------------------------------------------------------
  type(t_atm) ,intent(out)          :: ana       ! analysis
  type(t_atm) ,intent(inout)        :: ani       ! analysis increment
  type(t_atm) ,intent(inout)        :: fc        ! forecast
  logical     ,intent(in)           :: lvint     ! interpolate bg to p-levs.
  logical     ,intent(in)           :: lanagp    ! analyse gp instead of t.
  integer     ,intent(in)           :: lerror_tq ! treat errors in tq_tvgh
  integer     ,intent(in)           :: lplot     ! flag: write diagnostics
  logical     ,intent(in)           :: verb      ! verbose (call stop_time?)
  logical     ,intent(in) ,optional :: a2c       ! interpol./rotate to C grid
! logical     ,intent(in) ,optional :: keep_tinv ! keep near-surface temp.inv.?
  logical     ,intent(in) ,optional :: lextana   ! use external analysis?
    !----------------
    ! local variables
    !----------------
    character(len=128)    :: filename     ! temporary file name
    type (t_ctl)          :: ctl          ! grads ctl file data (model)
    real(wp) ,allocatable :: to (:,:,:,:) ! temporary: old temperature
    real(wp) ,allocatable :: rh (:,:,:,:) ! temporary: old rel.humidity
    real(wp) ,allocatable :: tv (:,:,:,:) ! temporary: virtual temp.
!   real(wp) ,allocatable :: x  (:,:,:,:) ! temporary: water loading
    integer  ,allocatable :: er (:,:,:,:) ! error code from tq_tvgh
    real(wp)              :: minrh        ! minimum rh value analysed
    integer               :: iel(3)       ! indices of error
    integer               :: ieb(3)       ! indices of error
    integer               :: i            ! index variable
    integer               :: k            ! level index
    logical               :: ac           ! copy of a2c
    integer               :: ie,je,ke,nd  ! dimensions
    type(t_grid) ,pointer :: g            ! grid pointer
    real(wp) ,allocatable :: dp_(:,:,:)   ! temporary for hydrostatic balancing
    real(wp) ,allocatable :: dph(:,:,:)   ! balanced increment of ln(ph)
    real(wp) ,allocatable :: dz_(:,:,:)   ! temporary for integration
    real(wp) ,allocatable :: fac(:,:,:)   ! scaling factor for geop. adjustment
    integer               :: k0           ! Start level of height relaxation
    real(wp), allocatable :: tinv(:,:,:,:)! near-surface temperature inversion
    logical               :: licon        ! branch for ICON (or COSMO) model
    logical               :: l_extana     ! local copy of lextana
    type (t_atm) ,pointer :: extana       ! external analysis

    if(.not.verb) call stop_time ('finish_ana')
    ac = .true.; if (present(a2c)) ac = a2c
    ie = fc% grid% shape (1)
    je = fc% grid% shape (2)
    ke = fc% grid% shape (3)
    nd = fc% grid% shape (4)
    licon = fc% grid% model == MO_ICON &  ! handle ICON ..
       .or. fc% grid% model == MO_COSMO   ! .. and COSMO in same way
    l_extana = .false.; if (present(lextana)) l_extana = lextana
    !--------------------------------------------------
    ! 1) average poles from neighbours (on icosahedron)
    !--------------------------------------------------
    if(verb) call stop_time ('(1) average poles')
    call average_pole_wind   (ani          ,n_av_pole)
    call average_pole_scalar (ani ,ani% ps ,n_av_pole)
    call average_pole_scalar (ani ,ani% t  ,n_av_pole)
    call average_pole_scalar (ani ,ani% rh ,n_av_pole)
    !-------------------------------------------
    ! 2) take vertical coordinates from forecast
    !-------------------------------------------
    if(verb) call stop_time ('(2) vertical coordinates from fc')
    select case (fc% grid% vct)
    case (VCT_P_HYB)
      !------------------------------------
      ! GME/HRM hybrid pressure coordinates
      !------------------------------------
      call allocate (ani, 'pf')
      call allocate (ani, 'ph')
      call allocate (ani, 'psr')
      ani% pf  = fc% pf
      ani% ph  = fc% ph
      ani% psr = fc% ps
    case default
      if (licon) then
        !--------------------------
        ! ICON hybrid z-coordinates
        !--------------------------
!       call set_geo  (ani, geof=.true.)
        !--------------------------------------------------------------
        ! Derive pressure increment consistent with hydrostatic balance
        !--------------------------------------------------------------
        if (.not. associated (fc% tv)) then
           call finish ('finish_ana','fc%tv not associated for ICON')
        end if
!       if (.not. associated (ani% tv)) then
!          call finish ('finish_ana','ani%tv not associated for ICON')
!       end if
        call allocate (ani, 'pf')
        call allocate (ani, 'ph')
        g => fc% grid
        allocate (dph(ie,je,nd))
        allocate (dp_(ie,je,nd))
        dph(:,:,:) = ani% ps(:,:,1,:) / fc% ps(:,:,1,:)   ! Increment of ln(ps)
        ani% ph(:,:,ke+1,:) = ani% ps(:,:,1,:)
        do k = ke,1,-1
          dp_(:,:,:) = gacc/Rd * (g% hhl(:,:,k,:) - g% hhl(:,:,k+1,:)) &
                               * ani% t (:,:,k,:) / fc% tv(:,:,k,:)**2
!!!                              this ^ should actually be the tv increment!
          ani% pf(:,:,k,:) = fc% pf(:,:,k,:) * (dph(:,:,:) + 0.5_wp*dp_(:,:,:))
          dph(:,:,:)       =                    dph(:,:,:) +        dp_(:,:,:)
          ani% ph(:,:,k,:) = fc% ph(:,:,k,:) *  dph(:,:,:)
        end do
        deallocate (dph, dp_)
      else
        !---------------------------
        ! COSMO hybrid z-coordinates
        !---------------------------
        call allocate (ani, 'pf')
        call calrho  ( fc % t              ,&
                       fc % pp             ,&
                       fc % q              ,&
                       fc % qcl + fc % qci ,&
                       fc % qr  + fc % qs  ,&
                       fc % grid% p0       ,&
                       ani% pf             ,&! -> rho
                       ie                  ,&
                       je                  ,&
                       ke                  ,&
                       Rd                  ,&
                       RDDRM1               )
!       ani% pf = gacc * ani% pf * ani% geof
        ani% pf = Rd   * ani% pf             ! p = R rho
        call deallocate (ani, 'geof')
!       call   allocate (ani, 'geoh')
!       ani% geoh = fc% geoh
      endif
    end select
    !---------------------------------------------
    ! 3) optionally write profiles for diagnostics
    !---------------------------------------------
    if (lplot /= 0) then
      if(verb) call stop_time ('(3) write profiles for diagnostics')
      call write_profiles (ani, 'analysis_incr',&
                           'analysis increment requested')
    endif
    !---------------------------------------
    ! 4) deallocate fields not used any more
    !---------------------------------------
    if (.not. licon) then
                       call deallocate (ani, 'ph'  )
      if (.not.lanagp) call deallocate (ani, 'geoh')
    end if
    !-----------------------------------------------
    ! 5) special handling for stratospheric humidity
    !-----------------------------------------------
    if(verb) call stop_time ('(5) stratospheric humidity')
    call hum_ana_incr (ani, fc)
    !----------------------------------------------
    ! 6) optionally write plot file for diagnostics
    !----------------------------------------------
    select case (ani% grid% vct)
    case (VCT_P_HYB)
      call deallocate (ani, 'pf')
    end select
    if (lplot /= 0) then
      if(verb) call stop_time ('(6) plot file for diagnostics')
      select case (fc% grid% gridtype)
      case (WMO6_GAUSSIAN,WMO6_LATLON)
        if (dace% lpio) then
          filename = path_file (aux, '^psas.grads')
          call init_ctl  (ctl, filename, fc% grid)
          call to_grads  (ctl, fc% grid)
          call to_grads  (ctl, fc)
          call to_grads  (ctl, ani, t=2)
          call write_ctl (ctl)
          call destruct  (ctl)
        endif
      case (DWD6_ICOSAHEDRON)
        call atm_to_gmt (ani, 'plot/anainc')
      end select
    endif
    !----------------------------------
    ! 7) set analysis time (ani and fc)
    !----------------------------------
    ani% time     = fc% time    !
    ani% ref_time = fc% time    ! set reference time to analysis time
!   fc % ref_time = fc% time    ! Note: fc%ref_time is needed by set_strat_q
    !--------------------------------
    ! 8) print out analysis increment
    !--------------------------------
    if(verb) call stop_time ('(8) print analysis increment')
    if(dace% lpio) write(6,'(/a//3x,a/)') repeat('-',79), 'analysis increment'
    call print (ani, verbose=.false.)
!   call print (ani, verbose=verb)
    !------------------
    ! 9) restore C-grid
    !------------------
    if (ac) then
      if(verb) call stop_time ('(9) restore C-grid')
      ! ani: analysis increment: interpolate u,v-increments  to C grid
      ! fc : first guess       : restore original u,v-fields to C grid
      call uvtrafo(ani,0)  ! transform wind to rotated coordinates
      if (ani% grid% arakawa == 'C') then
        call A_to_C_grid (ani ,restore=.false.)
        call A_to_C_grid (fc  ,restore=.true. )
      endif
    endif
    !----------------------------
    ! 10) allocate analysis state
    !----------------------------
    if(verb) call stop_time ('(10) allocate analysis state')
    call construct (ana, template = ani)
    !------------------------------------------------------
    ! 11) temporarily derive vt, gh in forecast; save t, rh
    !------------------------------------------------------
    if(verb) call stop_time ('(11) derive vt, gh')
    allocate(to(size(fc%t,1),size(fc%t,2),size(fc%t,3),size(fc%t,4)))
    allocate(tv(size(fc%t,1),size(fc%t,2),size(fc%t,3),size(fc%t,4)))
    allocate(er(size(fc%t,1),size(fc%t,2),size(fc%t,3),size(fc%t,4)))
    allocate(rh(size(fc%t,1),size(fc%t,2),size(fc%t,3),size(fc%t,4)))
!   allocate(x (size(fc%t,1),size(fc%t,2),size(fc%t,3),size(fc%t,4)))
    if (.not. associated (fc% rh)) then
      if(.not.associated (fc% pf)) call finish ('finish_ana',&
                                                'fc% pf is missing')
      call allocate (fc, 'rh')
      fc% rh = rh_q (fc% q, fc% t, fc% pf)
    endif
    to = fc% t
    rh = fc% rh
!   x  = 0._wp
    !------------------------------------------------------------
    ! Beware: if water loading is to be taken into account below,
    ! it must be done consistently everywhere...
    !------------------------------------------------------------
!   if (fc% grid% ivctype == IVCTYPE_ICON) then
!     if (associated(fc% qcl)) x = x + fc% qcl
!     if (associated(fc% qci)) x = x + fc% qci
!     if (associated(fc% qr )) x = x + fc% qr
!     if (associated(fc% qs )) x = x + fc% qs
!     if (associated(fc% qg )) x = x + fc% qg
!   end if
!   fc% t = tv_t_q (fc% t, fc% q, x) !! did not inline on NEC SX
    fc% t = fc% t * (1.0_wp + RDDRM1 * fc% q)
!   deallocate (x)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! NEC SX: work around floating-point exception in gh_rh_tp
!NEC$ collapse
    fc% rh = gh_rh      (fc% rh, .true., 200._wp, fc% pf)
    !---------------------------------------------
    ! 12) analysis = analysis increment + forecast
    !---------------------------------------------
    if(verb) call stop_time ('(12) analysis increment + forecast')
    if(lvint) then
      !------------------------------------------------------
      ! spline interpolation of analysis to new vertical grid
      !------------------------------------------------------
      select case (ani% grid% vct)
      case (VCT_P_HYB)
        !-----------------------------------
        ! GME/HRM hybrid pressure coordinate
        !-----------------------------------
        ani = ani + fc            ! preliminary analysis
        ana% psr = ani % ps       ! analysis surface pressure
        call set_p (ana)
        call set_p (ani)
        call vert_intp (ani, ana) ! interpolate fields
        if (lanagp) then
          call set_tv_geo (ana)
          ana% t = ana% tv
          call deallocate (ani, 'geoh')
          call deallocate (ani, 'tv')
        endif
      case (VCT_Z_HYB, VCT_Z_GEN)
        !------------------------------
        ! ICON height hybrid coordinate
        !------------------------------
        if (licon) then
          !-----------------------------------------------------------------------
          ! full level geopotential analysis increments are provided by the EnVar.
          ! Temperature increments are derived from the hydrostatic equation.
          !-----------------------------------------------------------------------
          if (lanagp) then
            ani = ani + fc  ! preliminary analysis
            !------------------------------------------------
            ! Use geopotential matching given pressure levels
            !------------------------------------------------
            ani% pf               = fc % pf       ! restore fg pressure
            ani% ph               = fc % ph       !         fg pressure
            !------------------------------------------------------
            ! interpolate geopotential half levels from full levels
            !------------------------------------------------------
            if (.not.associated (ani% geof)) &
              call finish('finish_ana','geof analysis not present')
            if (.not.associated (ani% geoh)) call allocate (ani, 'geoh')
            ani% geoh(:,:,1,:) = 0._wp
            do k = 2, ke
              ani% geoh(:,:,k,:) =                           ani% geof (:,:,k-1,:)  &
                                 +     (ani% geof(:,:,k,:) - ani% geof (:,:,k-1,:)) &
                                 / log (ani% pf  (:,:,k,:) / ani% pf   (:,:,k-1,:)) &
                                 * log (ani% ph  (:,:,k,:) / ani% pf   (:,:,k-1,:))
            end do
            ani% geoh(:,:,ke+1,:) = Rd * ani% t (:,:,ke,:)                      &
                                  * log (ani% ps(:,:,1 ,:) / ani% ph(:,:,ke+1,:))
            ani% geoh(:,:,:,:) = ani% geoh(:,:,:,:) + ani% grid% hhl(:,:,:,:) * gacc
            !------------------------------
            ! keep lowest temperature layer
            !------------------------------
            ani% geoh(:,:,ke+1,:) = ani% geoh(:,:,ke,:)                              &
                                  - ani% t   (:,:,ke,:) / fc% t         (:,:,ke  ,:) &
                     * gacc * (ani% grid% hhl(:,:,ke,:) - ani% grid% hhl(:,:,ke+1,:))
            !-------------------------
            ! re-calculate full levels
            !-------------------------
            do k = 1, ke
              ani% geof  (:,:,k,:) = 0.5_wp * (ani% geoh(:,:,k,:) + ani% geoh(:,:,k+1,:))
            end do
            !-------------------------------------------
            ! Analysis increment <- change of thickness:
            ! (z(k)-z(k+1))/tv(k) = const. for p-levels
            !-------------------------------------------
            do k = 1, ke
              ani% t(:,:,k,:) = fc% t(:,:,k,:) *                           &
                (ani% geoh     (:,:,k,:) - ani% geoh     (:,:,k+1,:)) /     &
               ((ani% grid% hhl(:,:,k,:) - ani% grid% hhl(:,:,k+1,:)) * gacc)
            end do
          !---------------------------------------------------
          ! virtual temperature is used as primary information
          ! geopotential / pressure is derived
          !---------------------------------------------------
          else
            !----------------------------------------------------------------------
            ! Derive geopotential for preliminary analysis, needed by interpolation
            !----------------------------------------------------------------------
            ani                      = ani     + fc  ! preliminary analysis
            call set_geo (ani, geof=.true.)
            allocate (dz_(ie,je,nd))
            if (lfix_geosp) then
              !------------------------------------------------------------
              ! 1) Adjust geopotential of analysis:
              !    - lower boundary condition: surface geopotential
              !    - adjust layer thickness using temperature increments
              !    - relaxation of analysis geopotential heights
              !      towards model heights near model top
              ! 2) Derive pressure on half and full levels using
              !    analyzed surface pressure and adjusted layer thicknesses
              !------------------------------------------------------------
              k0 = minloc (abs (ani% grid% hhl(fc% grid% lb(1),fc% grid% lb(2),:,1) &
                                - h0), 1)
              k0 = p_min (k0)
              k0 = min (max (k0, 5), ke/2)         ! safeguard: 5 <= k0 <= ke/2
              allocate (fac(ie,je,nd))
              dz_ = 0._wp
              fac = 1._wp
              do k = ke, 1, -1
                if (k < k0) then
                  fac(:,:,:) = (ani% grid% hhl(:,:,k ,:) - ani% grid% hhl(:,:,1,:)) &
                             / (ani% grid% hhl(:,:,k0,:) - ani% grid% hhl(:,:,1,:))
                end if
                dz_(:,:,:) = dz_(:,:,:)                                            &
                           + (ani% grid% hhl(:,:,k,:) - ani% grid% hhl(:,:,k+1,:)) &
                           * gacc * (ani% t(:,:,k,:) / fc% t(:,:,k,:) - 1._wp)
                ani% geoh(:,:,k,:) = ani% geoh(:,:,k  ,:) + dz_(:,:,:) * fac(:,:,:)
                ani% geof(:,:,k,:) =(ani% geoh(:,:,k  ,:) +        &
                                     ani% geoh(:,:,k+1,:) ) * 0.5_wp
                ani% ph  (:,:,k,:) = ani% ph(:,:,k+1,:)             &
                                   * exp ((ani% geoh(:,:,k+1,:) -   &
                                           ani% geoh(:,:,k  ,:)   ) &
                                          / (Rd * ani% t(:,:,k,:))  )
                ani% pf  (:,:,k,:) = sqrt (ani% ph(:,:,k,:) * ani% ph(:,:,k+1,:))
              end do
              deallocate (fac)
            else
              !--------------------------------------------------------------
              ! Derive fictitious geopotential matching given pressure levels
              !--------------------------------------------------------------
              ani% pf               = fc % pf       ! restore fg pressure
              ani% ph               = fc % ph       !         fg pressure
              !----------------------------------------------------------------------
              ! Adjust geopotential of analysis (pressure coordinates!)
              ! Lower boundary condition: fictitious surface geopotential of analysis
              !----------------------------------------------------------------------
              dz_      (:,:,     :) = Rd * ani% t (:,:,ke,:)                      &
                                    * log (ani% ps(:,:,1 ,:) / ani% ph(:,:,ke+1,:))
              ani% geoh(:,:,ke+1,:) = ani% geoh(:,:,ke+1,:) + dz_(:,:,:)
#if 1
              !-------------------------------------------
              ! Analysis increment -> change of thickness:
              ! (z(k)-z(k+1))/tv(k) = const. for p-levels
              !-------------------------------------------
              do k = ke, 1, -1
                dz_(:,:,:) = dz_(:,:,:)                                            &
                           + (ani% grid% hhl(:,:,k,:) - ani% grid% hhl(:,:,k+1,:)) &
                           * gacc * (ani% t(:,:,k,:) / fc% t(:,:,k,:) - 1._wp)
                ani% geoh(:,:,k,:) = ani% geoh(:,:,k,:) + dz_(:,:,:)
              end do
              do k = 1, ke
                ani% geof(:,:,k,:) = 0.5_wp * ( ani% geoh(:,:,k  ,:) &
                                              + ani% geoh(:,:,k+1,:) )
              end do
#else
              !----------------------------------------------------
              ! Integration of hydrostatic equation for analysis
              ! (less stable, but self-consistent within analysis
              ! if water loading is taken into account everywhere.)
              !----------------------------------------------------
              do k = ke, 1, -1
                dz_      (:,:,  :) = Rd * ani% t (:,:,k  ,:)                   &
                                   * log (ani% ph(:,:,k+1,:) / ani% ph(:,:,k,:))
                ani% geof(:,:,k,:) = ani% geoh(:,:,k+1,:) + dz_(:,:,:) * 0.5_wp
                ani% geoh(:,:,k,:) = ani% geoh(:,:,k+1,:) + dz_(:,:,:)
              end do
#endif
            end if
            deallocate (dz_)
          endif
          !-----------------------------------------------------------------------
          ! Preserve near-surface temperature inversion features from the forecast
          !-----------------------------------------------------------------------
          if (lkeep_tinv) then
            call derive_tinv (fc, tinv)
            k = lbound (tinv,3)
            ani% t(:,:,k:,:) = ani% t(:,:,k:,:) - tinv(:,:,:,:)
          end if
          !------------
          ! Interpolate
          !------------
          call set_geo   (ana, geof=.true.)
          call vert_intp (ani, ana)                   ! interpolate fields
          !---------------------------------------------------------------------
          ! Restore near-surface temperature inversion features, adjust pressure
          !---------------------------------------------------------------------
          ana% ph(:,:,ke+1,:) = ana% ps(:,:,1,:)      ! restore surface pressure
          if (lkeep_tinv) then
             k = lbound (tinv,3)
             ana% t(:,:,k:,:) = ana% t(:,:,k:,:) + tinv(:,:,:,:)
!            !------------------------------------------------
!            ! Re-derive hydrostatic pressure, assuming
!            ! vertical grid is identical for input and output
!            ! -> tv(k)*ln(ph(k)/ph(k+1)) = const. for layer k
!            !------------------------------------------------
!            do k = ke, 1, -1
!              ana% ph  (:,:,k,:) = ana% ph(:,:,k+1,:)               &
!                                 * exp (fc%  t(:,:,k,:) /           &
!                                        ana% t(:,:,k,:)             &
!                                        * log (  fc% ph(:,:,k  ,:)  &
!                                               / fc% ph(:,:,k+1,:)) )
!              ana% pf  (:,:,k,:) = sqrt (ana% ph(:,:,k,:) * ana% ph(:,:,k+1,:))
!            end do
             !----------------------------------------------------------
             ! Re-derive hydrostatic pressure only where temperature was
             ! restored, keep pressure levels above
             !----------------------------------------------------------
             do k = lbound (tinv,3), ke
               ana% ph  (:,:,k+1,:)   = ana% ph(:,:,k,:)             &
                                  / exp (fc%  t(:,:,k,:) /           &
                                         ana% t(:,:,k,:)             &
                                         * log (  fc% ph(:,:,k  ,:)  &
                                                / fc% ph(:,:,k+1,:)) )
               ana% pf  (:,:,k,:) = sqrt (ana% ph(:,:,k,:) * ana% ph(:,:,k+1,:))
             end do
             ana% ph(:,:,1,:) = ana% pf(:,:,1,:)**2 & ! top half-level presssure
                              / ana% ph(:,:,2,:)
             deallocate (tinv)
          else
            ana% ph(:,:,1,:) = ana% pf(:,:,1,:)**2 & ! top half-level presssure
                             / ana% ph(:,:,2,:)
!           call set_ph    (ana)
          end if
        else
          !--------------------------
          ! COSMO hybrid z-coordinate
          !--------------------------
          if (lanagp) call finish ('finish_ana',                    &
            'lanagp=T not implemented for COSMO hybrid z-coordinate')
          ana% pf = ani% pf + fc % pf       ! analysis    pressure
          ani     = ani     + fc            ! preliminary analysis
          ani% pf = fc % pf                 ! old         pressure
          call vert_intp (ani, ana)         ! interpolate fields
          call allocate (ana, 'pp')
          ana% pp = ana% pf - ana% grid% p0 ! pressure disturbance
        endif
      end select
    else
      !---------------------------------------------
      ! no spline interpolation, just add increments
      !---------------------------------------------
      ana = ani + fc            ! simply add, do not interpolate
      select case (ani% grid% vct)
      case (VCT_P_HYB)
        ana% psr = ana% ps
      case default
        if (licon) then
!         call set_ph  (ana)
          call set_geo (ana, geof=.true.)   ! Needed for out_add=T
        endif
      end select
    endif
    !-------------------------------------
    ! 13) force negative RH values to zero
    !-------------------------------------
    if(verb) call stop_time ('(13) force negative RH to 0')
    if (rh0 <= 0._wp) then
      minrh = p_min (minval(ana% rh))
      if (dace% lpio) then
        write(6,*)
        write(6,*) 'minimum RH:', minrh
        write(6,*)
        write(6,*) 'values < 0 set to 0 .'
        write(6,*)
      endif
      ana% rh = max (ana% rh, 0._wp)
    endif
    !---------------------------------------------
    ! 14) derive spec.humid. and temperature
    !     from gen.humid.  and virtual temperature
    !---------------------------------------------
    if(verb) call stop_time ('(14) spec.humid. and temperature')
    call allocate (ana, 'q')
    select case (ana% grid% vct)
    case (VCT_P_HYB)
      call set_p  (ana)
    case default
      if (licon) then
        if (.not. associated (ana%pf)) &
        call set_pf (ana)
      endif
    end select
    tv = ana% t
    if(verb) call stop_time ('tq_tvgh')

#if 1
    ! Process one diamond at a time, since processing full fields
    ! at once requires lots of temporary storage
    do i = 1, size (tv,4)
       call tq_tvgh_vec (size (tv,1)*size (tv,2)*size (tv,3),   &
                         ana% t (:,:,:,i), ana% q (:,:,:,i),    &! output
                              tv(:,:,:,i), ana% rh(:,:,:,i),    &
                                           ana% pf(:,:,:,i),    &! input
                              er(:,:,:,i))                       ! error flag
    end do
#else
    ! Original version calling elemental subroutine, slow on NEC SX:
    call tq_tvgh        (ana% t,  ana% q,            &! output
                              tv, ana% rh, ana% pf,  &! input
                         er)                          ! error flag
#endif

!   call rh_gh (ana% rh, ana% rh)
!   call tq_tvrh_noxfail (ana% t,  ana% q,       &! output
!                         tv, ana% rh, ana% pf,  &! input
!                         i_fail= er)             ! error flag
    !---------------------------------
    ! handle error from q,t derivation
    !---------------------------------
    if (any(er < 0)) then
      write(0,*) 'finish_ana, tq_tvgh failed, times:',count(er < 0)
      if (lerror_tq >= 2) then
       do i = 1,size(er,3)
        if (any(er(:,:,i,:)<0)) then
         write(0,*) &
          'finish_ana, tq_tvgh failed, level, times:',i,count(er(:,:,i,:)<0)
         if (lerror_tq >= 3) then
          iel = minloc (er(:,:,i,:))
          ieb = ana% lb((/1,2,4/)) + iel - 1
          write(0,*) 'setup_psas, i,j,k,d,t,q,tv,rh,pf,ifail=',&
                     ieb(1),ieb(2),i,ieb(3)                   ,&
            ana% t  (ieb(1),ieb(2),i,ieb(3))                  ,&
            ana% q  (ieb(1),ieb(2),i,ieb(3))                  ,&
                 tv (iel(1),iel(2),i,iel(3))                  ,&
            ana% rh (ieb(1),ieb(2),i,ieb(3))                  ,&
            ana% pf (ieb(1),ieb(2),i,ieb(3))                  ,&
                 er (iel(1),iel(2),i,iel(3))
         endif
        endif
       end do
      endif
    endif
    call p_barrier
    if (any(er < 0)) then
      if (lerror_tq >= 4) then
        call finish ('finish_ana','tq_tvgh failed')
      else
        call message('finish_ana','tq_tvgh failed')
      endif
    endif
    !------------------------------------------------
    ! 15) special handling for stratospheric humidity
    !------------------------------------------------
    select case (fc% grid% model)
    case (MO_ICON, MO_GME, MO_IFS)
      if(verb) call stop_time ('(15) stratospheric humidity')
      call set_geo     (ana, geof=.true.)
      call set_strat_q (ana, ref=fc)
    end select
    !-------------------------------------------
    ! 16) restore temperature, relative humidity
    !-------------------------------------------
    if(verb) call stop_time ('(16) restore T, RH')
    ana% t  = t_tv_q (tv, ana% q)
    if(verb) call stop_time ('rh_q !! for NEC SX')
    ana% rh = rh_q   (ana% q, ana% t, ana% pf)
    !--------------------------------------------------------
    ! 17) copy some quantities from the forecast and ext. ana
    !--------------------------------------------------------
    if(verb) call stop_time ('(17) copy quantities from forecast and ext. ana')
    if ( l_extana ) then
      call prep_extana (extana, fc)
    else
      nullify (extana)
    endif
    if (associated (extana)) then
      call select_params (ana, fc, pass_fields, extana, read_extana, dealloc=.false.)
    else
      call select_params (ana, fc, pass_fields, dealloc=.false.)
    end if
    !-------------------------------------------------------------
    ! 18) provide zero increments for quantities not yet analysed.
    !     no cloud water for rel.humidity < 0.9
    !-------------------------------------------------------------
    if (associated (fc% qcl)) then
      if(verb) call stop_time ('(18) zero increments')
      call  allocate (ana, 'qcl')
      where (ana% rh < rh0_qcl0)
        ana% qcl = 0._wp
      elsewhere
        ana% qcl = fc% qcl
      endwhere
    end if
    !-----------------------------------------------------
    ! 19) restore temperature and rel.humidity in forecast
    !-----------------------------------------------------
    if(verb) call stop_time ('(19) restore T, RH in FC')
    fc% t  = to
    fc% rh = rh
    deallocate (to)
    deallocate (rh)
    deallocate (tv)
    deallocate (er)
    !-------------------------------------
    ! 20) re-derive ps from p (COSMO,ICON)
    !-------------------------------------
    if(verb) call stop_time ('(20) re-derive ps from p (COSMO,ICON)')
    if (ani% grid% vct /= VCT_P_HYB) call set_ps (ana)
    !-------------------------------------------------
    ! 21) derive t2m increment from lowest model level
    !     similarly for rh2m ("poor man's diagnostic")
    !     (note that model forecast may have td2m).
    !-------------------------------------------------
    if (associated (fc% t2m)) then
      call allocate (ani, 't2m')
      call allocate (ana, 't2m')
      if (version >= 3) then
        if(verb) call stop_time ('(21) derive t2m,rh2m increment')
        ani% t2m (:,:,1,:) = ana% t   (:,:,ke,:) -  fc% t   (:,:,ke,:)
        ana% t2m (:,:,1,:) =  fc% t2m (:,:, 1,:) + ani% t2m (:,:, 1,:)

        if (associated (fc% rh2m)) then
          call allocate (ana, 'rh2m')
          call allocate (ani, 'rh2m')
          ani% rh2m(:,:,1,:) = ana% rh  (:,:,ke,:) -  fc% rh  (:,:,ke,:)
          ana% rh2m(:,:,1,:) =  fc% rh2m(:,:, 1,:) + ani% rh2m(:,:, 1,:)
          ana% rh2m = min (max (ana% rh2m, 1.e-30_wp), 1._wp)

          !----------------------------------------------------
          ! Apply same increments to land-tile averages for now
          !----------------------------------------------------
          if (associated (fc% t2m_land) .and. associated (fc% rh2m_land)) then
            call allocate (ana, 't2m_land')
            call allocate (ana, 'rh2m_land')
            ana% t2m_land (:,:,1,:) = fc% t2m_land (:,:,1,:) + ani% t2m (:,:,1,:)
            ana% rh2m_land(:,:,1,:) = fc% rh2m_land(:,:,1,:) + ani% rh2m(:,:,1,:)
            ana% rh2m_land = min (max (ana% rh2m_land, 1.e-30_wp), 1._wp)
          end if

          call deallocate (ani, 'rh2m')
        end if
      else
        !-----------------------------------------------
        ! For the time being we do not analyse T2m,
        ! as long as we are forced to pass the FG as is.
        !-----------------------------------------------
        ani% t2m (:,:,1,:) = 0._wp
        ana% t2m (:,:,1,:) = fc% t2m (:,:, 1,:)
      end if
      call deallocate (ani, 't2m')
    end if
    !---------------------------
    ! 22) derive wind adjustment
    !---------------------------
    if (wind_adjust > 0) then
      if(verb) call stop_time ('(22) derive wind adjustment (COSMO,ICON)')
      call adjust_wind (fc, ana)
    end if

  end subroutine finish_ana

!==============================================================================

  subroutine derive_tinv (bg, tinv)
    !--------------------------------------------------------------
    ! Derive "near-surface temperature inversion" as departure from
    ! a constant lapse rate profile, whose structure shall be
    ! preserved when applying the analysis increments.
    ! The lapse rate is derived from background T at 500m and 1000m
    ! above the surface, consistent with ICON's temperature_intp.
    !--------------------------------------------------------------
    type(t_atm),           intent(in)  :: bg            ! Background state
    real(wp), allocatable, intent(out) :: tinv(:,:,:,:) ! Temperature inversion
    !----------------
    ! local variables
    !----------------
    integer               :: ke              ! number of levels
    integer               :: k, l            ! level index
    integer               :: k0              ! level index
    integer               :: k1, k2          ! bracketing level indices
    integer               :: i, j, d         ! loop indices
    integer,  allocatable :: l1 (:,:,  :)    ! auxiliary levels
    integer,  allocatable :: l2 (:,:,  :)    ! auxiliary levels
    real(wp)              :: t1, t2          ! Temperature @ zpbl1,zpbl2
    real(wp)              :: w1, w2          ! Interpolation weights
    real(wp), allocatable :: t0 (:,:,  :)    ! Reference temperature for extrap.
    real(wp), allocatable :: t_z(:,:,  :)    ! Temperature vertical gradient
    real(wp), allocatable :: z  (:,:,:,:)    ! Height above ground
    type(t_grid) ,pointer :: g               ! grid pointer
    integer               :: lb(4), ub(4)    ! array bounds

    g  => bg% grid
    lb =   g% lb(:)
    ub =   g% ub(:)
    ke =      ub(3)

    if (zpbl1 > zpbl2)            &
         call finish ("derive_tinv","zpbl1 > zpbl2")
    if (.not.associated (g% hhl)) &
         call finish ("derive_tinv","hhl not associated!")
    !--------------------------------------------------------------------
    ! Determine level indices bracketing zpbl2, starting near the surface
    !--------------------------------------------------------------------
    k2 = 0
    do k = ke, 2, -1
       if (minval (g% hhl(:,:,k,:) - g% hhl(:,:,ke+1,:)) >= zpbl2) then
          k2 = k-1
          exit
       end if
    end do
    if (k2 == 0) call finish ("derive_tinv","k2==0")

    allocate (z(lb(1):ub(1), lb(2):ub(2), k2:ke, lb(4):ub(4)))
    !-----------------------------------------------------------------------
    ! Derive layer mean heights above ground, level indices bracketing zpbl1
    !-----------------------------------------------------------------------
    do k = k2, ke
       z(:,:,k,:) = (g% hhl(:,:,k   ,:) + g% hhl(:,:,k+1,:)) * 0.5_wp &
                  -  g% hhl(:,:,ke+1,:)
    end do
    k1 = 0
    do k = k2, ke
       if (maxval (z(:,:,k,:)) < zpbl1) then
          k1 = k
          exit
       end if
    end do
    if (k1 == 0) call finish ("derive_tinv","k1==0")
    !--------------
    ! Sanity checks
    !--------------
    if (maxval (z(:,:,k1,:)) > zpbl1) call finish ("derive_tinv","z(k1) > zpbl1")
    if (minval (z(:,:,k2,:)) < zpbl2) call finish ("derive_tinv","z(k2) < zpbl2")

    allocate (l1(lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
    allocate (l2(lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
    l1 = 0
    l2 = 0
    !--------------------------------------------------
    ! Find layer indices immediately above zpbl1, zpbl2
    !--------------------------------------------------
    do d = lb(4), ub(4)
       do j = lb(2), ub(2)
          do i = lb(1), ub(1)
             do k = k1, k2, -1                  ! Scan from below
                if (z(i,j,k,d) > zpbl1) then
                   l1(i,j,d) = k                ! First one above
                   exit
                end if
             end do
             do k = k2, k1-1                    ! Scan from above
                if (z(i,j,k+1,d) < zpbl2) then
                   l2(i,j,d) = k                ! First one below
                   exit
                end if
             end do
          end do
       end do
    end do
    if (any (l1 == 0) .or. any (l2 == 0)) then
       call finish ("derive_tinv","determination of l1,l2 failed.")
    end if

    allocate (t0 (lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
    allocate (t_z(lb(1):ub(1), lb(2):ub(2), lb(4):ub(4)))
    !-----------------------------------------------------------------------
    ! Derive gradient between zpbl1 and zpbl2, reference temperature @ zpbl1
    !-----------------------------------------------------------------------
    do d = lb(4), ub(4)
       do j = lb(2), ub(2)
          do i = lb(1), ub(1)
             k  = l1(i,j,d)
             l  = l2(i,j,d)
             w1 = (z(i,j,k,d) - zpbl1) / (z(i,j,k,d) - z(i,j,k+1,d))
             w2 = (z(i,j,l,d) - zpbl2) / (z(i,j,l,d) - z(i,j,l+1,d))
             t1 = (1._wp-w1) * bg% t(i,j,k,d) + w1 * bg% t(i,j,k+1,d)
             t2 = (1._wp-w2) * bg% t(i,j,l,d) + w2 * bg% t(i,j,l+1,d)
             t0 (i,j,d) = t1
             t_z(i,j,d) = (t2 - t1) / (zpbl2 - zpbl1)
             t_z(i,j,d) = min (max (t_z(i,j,d), dtdz_min), dtdz_max)
          end do
       end do
    end do
    !------------------------------------------------------
    ! Derive near-surface inversion as difference between
    ! background state and extrapolated temperature profile
    !------------------------------------------------------
    k0 = minval (l1)
    allocate (tinv(lb(1):ub(1),lb(2):ub(2),k0:ke,lb(4):ub(4)))
    do d = lb(4), ub(4)
       do j = lb(2), ub(2)
          do i = lb(1), ub(1)
             do k = k0, ke
                if (z(i,j,k,d) < zpbl1) then
                   tinv(i,j,k,d) = bg% t(i,j,k,d)                              &
                                 - (t0(i,j,d) + t_z(i,j,d) * (z(i,j,k,d)-zpbl1))
                else
                   tinv(i,j,k,d) = 0._wp
                end if
             end do
          end do
       end do
    end do

    ! Debug: check where departure from gradient extrapolation is largest
!   block
!     integer :: j(4), pe
!     flush(6)
!     call p_barrier
!     do pe = 0, dace% npe-1
!      if (pe == dace% pe) then
!        print *, pe,"derive_tinv: min,max(tinv):", minval (tinv), maxval (tinv)
!        print *, pe,"derive_tinv: min,max(t_z) :", minval (t_z),  maxval (t_z)
!        print *, "lb:", lb
!        print *, "ub:", ub
!        print *
!        j    = minloc (tinv) + lb(:) - 1
!        j(3) = minloc (tinv (j(1),j(2),ke:k0:-1,j(4)), dim=1)
!        if (j(3) > 4) then             ! Possibly interesting inversion
!           print *, "lat,lon,hsurf=", real([ &
!                g%rlat(j(1),j(2),1,j(4))*180/pi,&
!                g%rlon(j(1),j(2),1,j(4))*180/pi,&
!                g%hhl (j(1),j(2),ke+1,j(4))])
!        end if
!        print *, "<< j:", j
!        print *, "  t0:", real (t0   (j(1),j(2),j(4)))
!        print *, " t_z:", real (t_z  (j(1),j(2),j(4)))
!        print *, "   z:", real (z    (j(1),j(2),ke:k2:-1,j(4)))
!        print *, "t_fc:", real (bg%t (j(1),j(2),ke:k2:-1,j(4)))
!        print *, "tinv:", real (tinv (j(1),j(2),ke:k0:-1,j(4)))
!        print *
!        j    = maxloc (tinv) + lb(:) - 1
!        j(3) = maxloc (tinv (j(1),j(2),ke:k0:-1,j(4)), dim=1)
!        if (j(3) > 4) then             ! Possibly interesting inversion
!           print *, "lat,lon,hsurf=", real([ &
!                g%rlat(j(1),j(2),1,j(4))*180/pi,&
!                g%rlon(j(1),j(2),1,j(4))*180/pi,&
!                g%hhl (j(1),j(2),ke+1,j(4))])
!        end if
!        print *, ">> j:", j
!        print *, "  t0:", real (t0   (j(1),j(2),j(4)))
!        print *, " t_z:", real (t_z  (j(1),j(2),j(4)))
!        print *, "   z:", real (z    (j(1),j(2),ke:k2:-1,j(4)))
!        print *, "t_fc:", real (bg%t (j(1),j(2),ke:k2:-1,j(4)))
!        print *, "tinv:", real (tinv (j(1),j(2),ke:k0:-1,j(4)))
!        print *
!      end if
!      flush(6)
!      call p_barrier
!     end do
!   end block

  end subroutine derive_tinv

!==============================================================================

  subroutine adjust_wind_ens (fc, ana)
    !---------------------------------------------------------------------
    ! Derive wind adjustment for the analysis ensemble in order to control
    ! the induced surface pressure tendency.  Adjustment modes:
    ! 0 : assume that increments should be divergence-free
    ! 1 : derive adjustments from analysed surface pressure tendencies
    !---------------------------------------------------------------------
    type(t_atm), intent(in)    :: fc (:)
    type(t_atm), intent(inout) :: ana(:)
    !----------------
    ! local variables
    !----------------
    integer     :: k
    type(t_atm) :: mean         ! Temporary for mean
    type(t_atm) :: fg_spread    ! Forecast ensemble spread
    type(t_atm) :: an_spread    ! Analysis ensemble spread
    type(t_atm) :: profile      ! Wind error profile
    real(wp)    :: w_fc, w_an   ! Weights

    if (wind_adjust <= 0) return

    if (adjust_mode > 0) then
       if (.not. associated (fc (1)% dpsdt)) &
            call finish ("adjust_wind_ens","fc% dpsdt not associated!")
       if (.not. associated (ana(1)% dpsdt)) &
            call finish ("adjust_wind_ens","ana% dpsdt not associated!")

       call construct  (mean,      template=fc(1), alloc='u v')
       call construct  (fg_spread, template=fc(1), alloc='u v')
       call construct  (an_spread, template=fc(1), alloc='u v')
       call mean_stdev (mean, fg_spread, fc)
       call mean_stdev (mean, an_spread, ana)
       call construct  (profile,   template=fc(1), alloc='u')
       w_an = min (max (w_spread_an, 0._wp), 1._wp)
       w_fc = 1._wp - w_an
       profile% u = sqrt ( w_fc * (fg_spread% u**2 + fg_spread% v**2) &
                         + w_an * (an_spread% u**2 + an_spread% v**2) )
       call destruct   (mean)
       call destruct   (fg_spread)
       call destruct   (an_spread)
    end if

    do k = 1, size (fc)
       if (fc(k)% members <= 0 .and. .not. adjust_det) cycle
       if (fc(k)% members >  0 .and. .not. adjust_ens) cycle
       if (dace% lpio) &
            write(*,'(a,i8)') "Adjusting wind, member =", fc(k)% member
       call adjust_wind (fc(k), ana(k), profile)
    end do

    if (adjust_mode > 0) then
       call destruct (profile)
    end if

  end subroutine adjust_wind_ens

!==============================================================================

  subroutine adjust_wind (fc, ana, profile)
    !--------------------------------------------------------------------
    ! Derive wind adjustment for the analysis state in order to reduce
    ! the induced surface pressure tendency.  The adjustment is derived
    ! from the requirement that the increments should be divergence-free.
    !--------------------------------------------------------------------
    type(t_atm), intent(in)             :: fc           ! Forecast state
    type(t_atm), intent(inout)          :: ana          ! Analysis state
    type(t_atm), intent(in),   optional :: profile      ! Wind error profile
    !----------------
    ! local variables
    !----------------
    integer               :: k, ke
    real(wp), allocatable :: phref_a(:) ! "Reference" half-level pressure
    real(wp), allocatable :: phref_b(:) ! "Reference" half-level pressure
    type(t_grid), pointer :: g
    type(t_atm)           :: xi         ! Control vector, temporary
    type(t_atm)           :: wk         ! Internal work space (u,v,div,vrt)
    type(t_atm)           :: b          ! RHS for adjustment calculation
    type(t_atm)           :: r, q       ! Temporaries (div, vrt)
    type(t_atm)           :: s, p       ! Temporaries (u, v)
    type(t_atm)           :: diag_fc    ! Diagnostic for background
    type(t_atm)           :: diag_ana   ! Diagnostic for background
    type(t_atm)           :: diag_adj   ! Diagnostic for adjustment
    real(wp)              :: max_fc     ! max of dpsdt, fc  (diag.)
    real(wp)              :: max_fc0    ! max of dpsdt, fc  (model)
    real(wp)              :: max_ana    ! max of dpsdt, ana (diag.)
    real(wp)              :: max_ana0   ! max of dpsdt, ana (LETKF)
    real(wp)              :: max_ani    ! max of dpsdt, ana increment
    real(wp)              :: max_ana_diff     ! max of dpsdt, ana diff between diagnostic and LETKF
    real(wp)              :: max_ana_diff_adj ! max of dpsdt, ana diff between diagnostic adjusted and LETKF
    real(wp)              :: max_adj    ! max of dpsdt, adjusted ana
    real(wp)              :: rms_fc     ! RMS of dpsdt, fc  (diag.)
    real(wp)              :: rms_fc0    ! RMS of dpsdt, fc  (model)
    real(wp)              :: rms_ana    ! RMS of dpsdt, ana (diag.)
    real(wp)              :: rms_ana0   ! RMS of dpsdt, ana (LETKF)
    real(wp)              :: rms_ani    ! RMS of dpsdt, ana increment
    real(wp)              :: rms_adj    ! RMS of dpsdt, adjusted ana
    real(wp)              :: gamma      ! Residual (normal equations)
    real(wp)              :: resne      ! Rel. residual (normal equations)
    real(wp)              :: alpha, beta
    real(wp)              :: gamma0
    character(3)          :: suf        ! Suffix of grib files
    real(wp)              :: dz0i       ! Inverse layer thickness
    real(wp)              :: sqz        ! Squeezing factor
    real(wp)              :: tmp, w     ! Temporary, weight
    integer               :: i, j, d    ! Indices

    if (wind_adjust <= 0) return
    if (fc% members <= 0 .and. .not. adjust_det) return
    if (fc% members >  0 .and. .not. adjust_ens) return

    g => fc% grid
    select case (g% gridtype)
    case (DWD6_ICON)
    case (WMO6_LATLON, WMO6_ROTLL)
    case default
       call finish ("adjust_wind","TODO (unsupported grid)")
    end select

!   suf = "grb"
    suf = "det"
    if (fc% members > 0) suf = char3(fc% member)

    call math_oper_init (g)
    call construct (xi, g, time=ana% time, ref_time=ana% ref_time)
    call construct (wk, g, time=ana% time, ref_time=ana% ref_time)
    call construct (b,  g, time=ana% time, ref_time=ana% ref_time)
    call allocate  (xi, "u v div vrt")
    call allocate  (wk, "u v div vrt dpsdt")
    call allocate  (b,  "div vrt")

    call construct (diag_fc,  g, time=fc%  time, ref_time=fc%  ref_time)
    call construct (diag_ana, g, time=ana% time, ref_time=ana% ref_time)
    call construct (diag_adj, g, time=ana% time, ref_time=ana% ref_time)
    call allocate  (diag_fc,  "dpsdt")
    call allocate  (diag_ana, "dpsdt")
    call allocate  (diag_adj, "dpsdt")
    if (diag_level > 1) then
       call allocate  (diag_adj, "u v div")
    end if
    !---------------------------------------
    ! Derive "reference" half-level pressure
    !---------------------------------------
    ke = g% nz
    allocate (phref_a(ke+1))
    allocate (phref_b(ke+1))
    do k = 1, ke+1
       phref_a(k) = sum (ana% ph(:,:,k,:))
       phref_b(k) = sum (fc % ph(:,:,k,:))
    end do
    phref_a = p_sum (phref_a) / g% nxny
    phref_b = p_sum (phref_b) / g% nxny
    !-----------------------------------------
    ! Derive diagnostics: preliminary analysis
    !-----------------------------------------
    xi = 0._wp
    call apply_a (ana, xi, refatm=ana, ph_ref=phref_a)
    diag_ana% dpsdt(:,:,1,:) = sum (xi% div, dim=3)
    max_ana = maxval (abs (diag_ana% dpsdt))
    max_ana = p_max (max_ana)
    rms_ana = sum (diag_ana% dpsdt ** 2)
    rms_ana = sqrt (p_sum (rms_ana) / g% nxny)
    if (associated (ana% dpsdt)) then
       max_ana0 = maxval (abs (ana% dpsdt))
       max_ana0 = p_max (max_ana0)
       rms_ana0 = sum (ana% dpsdt ** 2)
       rms_ana0 = sqrt (p_sum (rms_ana0) / g% nxny)
       max_ana_diff = maxval (abs (ana% dpsdt - diag_ana% dpsdt))
    end if
    !-------------------------------
    ! Derive diagnostics: background
    !-------------------------------
    b = 0._wp
    call apply_a (fc,  b,  refatm=fc,  ph_ref=phref_b)
    diag_fc % dpsdt(:,:,1,:) = sum (b % div, dim=3)
    max_fc  = maxval (abs (diag_fc% dpsdt))
    max_fc  = p_max (max_fc)
    rms_fc  = sum (diag_fc% dpsdt ** 2)
    rms_fc  = sqrt (p_sum (rms_fc) / g% nxny)
    if (associated (fc% dpsdt)) then
       max_fc0 = maxval (abs (fc% dpsdt))
       max_fc0 = p_max (max_fc0)
       rms_fc0 = sum (fc% dpsdt ** 2)
       rms_fc0 = sqrt (p_sum (rms_fc0) / g% nxny)
    end if
    !--------------------------------
    ! Diagnostics: analysis increment
    !--------------------------------
    max_ani = maxval (abs (diag_ana% dpsdt - diag_fc% dpsdt))
    max_ani = p_max (max_ani)
    rms_ani = sum ((diag_ana% dpsdt - diag_fc% dpsdt)** 2)
    rms_ani = sqrt (p_sum (rms_ani) / g% nxny)

    if (dace% lpio) then
       write(6,'()')
       if (associated (fc% dpsdt)) then
        write(6,'(A,f10.4,A)') "  max|dps/dt|  fc (model) =", max_fc0,  " Pa/s"
        write(6,'(A,f10.4,A)') "  rms(dps/dt)  fc (model) =", rms_fc0,  " Pa/s"
       end if
       write(6,'(A,f10.4,A)')  "  max|dps/dt|  fc (diag)  =", max_fc,   " Pa/s"
       write(6,'(A,f10.4,A)')  "  rms(dps/dt)  fc (diag)  =", rms_fc,   " Pa/s"
       if (associated (ana% dpsdt)) then
        write(6,'(A,f10.4,A)') "  max|dps/dt|  ana        =", max_ana0, " Pa/s"
        write(6,'(A,f10.4,A)') "  rms(dps/dt)  ana        =", rms_ana0, " Pa/s"
        write(6,'(A,f10.4,A)') "  max|dps/dt|  ana_diff   =", max_ana_diff, " Pa/s"
       end if
       write(6,'(A,f10.4,A)')  "  max|dps/dt|  ana (diag) =", max_ana,  " Pa/s"
       write(6,'(A,f10.4,A)')  "  rms(dps/dt)  ana (diag) =", rms_ana,  " Pa/s"
       write(6,'(A,f10.4,A)')  "  max|dps/dt|  inc        =", max_ani,  " Pa/s"
       write(6,'(A,f10.4,A)')  "  rms(dps/dt)  inc        =", rms_ani,  " Pa/s"
       write(6,'()')
    end if

    if (diag_level > 0) then
       ! TODO: make filenames depend on ensemble identification
       call write_grib (diag_fc,  ref=.false., mode='w',                    &
                        file=path_file (aux, trim (diag_base)//'_fg.' //suf))
       call write_grib (diag_ana, ref=.false., mode='w',                    &
                        file=path_file (aux, trim (diag_base)//'_ana.'//suf))
    end if

    if (wind_adjust > 1) then
       !-----------
       ! Derive RHS
       !-----------
       b% vrt = 0._wp
       b% div = b% div - xi% div

       if (adjust_mode > 0 .and. present    (profile)    &
                           .and. associated (fc % dpsdt) &
                           .and. associated (ana% dpsdt) ) then
          !------------------------------------------------------------
          ! Derive divergence adjustment corresponding to error profile
          !------------------------------------------------------------
          wk% dpsdt = ana% dpsdt - fc% dpsdt
          do k = 1, ke
             wk% u(:,:,k,:) = (ana% ph(:,:,k+1,:) - ana% ph(:,:,k,:)) &
                            * profile% u(:,:,k,:)
          end do
          wk% v(:,:,1,:) = sum (wk% u, dim=3)   ! Normalization factor
          do k = 1, ke
             wk% div(:,:,k,:) = wk% dpsdt(:,:,1,:) * wk% u(:,:,k,:) &
                                                   / wk% v(:,:,1,:)
          end do
          b% div = b% div + wk% div
       end if

       xi     = 0._wp
       select case (tolower (solver))
       case ("cgls")
          ! Allocate temporaries
          call alloc_dr (r)
          call alloc_dr (q)
          call alloc_uv (s)
          call alloc_uv (p)
          call initialize ()
          call cgls ()
          ! Cleanup
          call destruct (r)
          call destruct (q)
          call destruct (s)
          call destruct (p)
!      case ("lsqr")
          ! To be implemented...
       case default
          call finish ("adjust_wind","not implemented: " // trim (solver))
       end select

       if (wind_adjust == 4) then
          !----------------------------------------------------
          ! Reduce wind adjustment for strongly squeezed layers
          !----------------------------------------------------
          select case (g% vct)
          case (VCT_Z_HYB, VCT_Z_GEN)
             tmp = 1._wp / (sqz_trans(2) - sqz_trans(1))
             do       d = g% lb(4), g% ub(4)
                do    k = 1, ke
                   dz0i = 1._wp / (g% ak(k) - g% ak(k+1))
                   do j = g% lb(2), g% ub(2)
!DIR$ IVDEP
                   do i = g% lb(1), g% ub(1)
                      sqz = (g% hhl(i,j,k,d) - g% hhl(i,j,k+1,d)) * dz0i
                      sqz = min (max (sqz, sqz_trans(1)), sqz_trans(2))
                      w   = (sqz - sqz_trans(1)) * tmp
                      xi% u(i,j,k,d) = xi% u(i,j,k,d) * w
                      xi% v(i,j,k,d) = xi% v(i,j,k,d) * w
                   end do
                   end do
                end do
             end do
          case default
             call finish ("adjust_wind","wind_adjust=4 requires hybrid z")
          end select
       end if

       if (wind_adjust > 2) then
          ana% u = ana% u + xi% u
          ana% v = ana% v + xi% v
       end if
       !-------------------------------
       ! Derive diagnostics: adjustment
       !-------------------------------
       xi% div = 0._wp
       xi% vrt = 0._wp
       call apply_a (xi, xi, refatm=ana, ph_ref=phref_a)
       diag_adj% dpsdt(:,:,1,:) = sum (xi% div, dim=3)
       if (diag_level > 1) then
          diag_adj% u   = xi% u
          diag_adj% v   = xi% v
          diag_adj% div = xi% div
       end if
       !---------------------------------------------
       ! Diagnostics: adjusted analysis and increment
       !---------------------------------------------
       diag_ana% dpsdt = diag_ana% dpsdt + diag_adj% dpsdt
       max_adj = maxval (abs (diag_ana% dpsdt))
       if (associated (ana% dpsdt)) then
         max_ana_diff_adj = maxval (abs (ana% dpsdt - diag_ana% dpsdt))
       end if
       max_adj = p_max (max_adj)
       rms_adj = sum (diag_ana% dpsdt ** 2)
       rms_adj = sqrt (p_sum (rms_adj) / g% nxny)
       if (wind_adjust > 2) then
          max_ani = maxval (abs (diag_ana% dpsdt - diag_fc% dpsdt))
          max_ani = p_max (max_ani)
          rms_ani = sum ((diag_ana% dpsdt - diag_fc% dpsdt)** 2)
          rms_ani = sqrt (p_sum (rms_ani) / g% nxny)
       end if
       if (dace% lpio) then
          write(6,'()')
          if (wind_adjust > 2) then
            write(6,'(A,f10.4,A)')"  max|dps/dt|  inc_a =", max_ani, " Pa/s"
            write(6,'(A,f10.4,A)')"  rms(dps/dt)  inc_a =", rms_ani, " Pa/s"
          end if
          if (associated (ana% dpsdt)) then
            write(6,'(A,f10.4,A)') "  max|dps/dt|  ana_diff_adj   =", max_ana_diff_adj, " Pa/s"
          end if
          write(6,'(A,f10.4,A)')  "  max|dps/dt|  final =", max_adj, " Pa/s"
          write(6,'(A,f10.4,A)')  "  rms(dps/dt)  final =", rms_adj, " Pa/s"
          write(6,'()')
       end if

       if (diag_level > 0) then
          call write_grib (diag_adj, ref=.false., mode='w',                    &
                           file=path_file (aux, trim (diag_base)//'_adj.'//suf))
       end if
    end if

    call destruct (diag_fc)
    call destruct (diag_ana)
    call destruct (diag_adj)
    call destruct (xi)
    call destruct (wk)
    call destruct (b)
    call math_oper_cleanup ()
  contains
    subroutine cgls ()
      !--------------------------------------------
      ! Conjugate gradient for least squares (CGLS)
      ! Notation: C.Paige, M.Saunders, ACM 8 (p.57)
      !--------------------------------------------
      integer     :: ic                 ! Iteration counter
      real(wp)    :: gamm1, delta       ! Temporaries
      ic = 0
      if (dace% lpio) write(*,'(A10,4A16)')                  &
           " CGLS iter", "old s^2", "new s^2", "beta", "resne"
      do
         ic = ic + 1
         if (ic > max_iter) then
            if (dace% lpio) write(*,'(/,A)') " max_iter reached, exiting."
            exit
         end if
!        if (dace% lpio) print *, "CGLS iteration", ic
         q     = 0._wp                  ! q = A p
         call apply_a   (p, q, refatm=ana, ph_ref=phref_a)
         delta = dot_prod (q, q)
         delta = max (delta, EPSILON (delta))
         alpha = gamma / delta          ! alpha = gamma / |q|^2
         xi    = xi + alpha * p
         r     = r  - alpha * q         ! Update residual
         s     = 0._wp                  ! Residual of normal eqn.: s = A^T r
         call apply_a_t (r, s, refatm=ana, ph_ref=phref_a)
         gamm1 = gamma                  ! gamm1 = old norm^2 of residual s
         gamma = dot_prod (s, s)        ! gamma = new norm^2 of residual s
         resne = gamma / gamma0         ! Relative residual of normal equation
         beta  = gamma / gamm1
         p     = s + beta * p
         if (dace% lpio) write(*,'(i10,4f16.6)') ic, gamm1, gamma, beta, resne
         if (ic < min_iter) cycle
         if (resne < tol(1)) then
            if (dace% lpio) then
               write(*,'(/,1x,A,2f16.6)')                               &
                    "Convergence criterion 1 reached, resne < tol(1):", &
                    resne, tol(1)
            end if
            exit
         end if
         if (beta > tol(2) .and. tol(2) > 0) then
            if (dace% lpio) then
               write(*,'(/,1x,A,2f16.6)')                              &
                    "Convergence criterion 2 reached, beta > tol(2):", &
                    beta, tol(2)
            end if
            exit
         end if
      end do
    end subroutine cgls
    !-
    subroutine initialize ()
      ! Start vector
      xi        = 0._wp                                 ! x_0 = 0
      r         = b                                     ! r_0 = b
      s         = 0._wp                                 ! s_0 = A^T r_0
      call apply_a_t (r, s, refatm=ana, ph_ref=phref_a)
      p         = s                                     ! p_0 = s_0
      gamma0    = dot_prod (s, s)                       ! gamma0 = |s_0|^2
      gamma0    = max (gamma0, EPSILON (gamma0))
      gamma     = gamma0
      beta      = 0._wp
      resne     = 1._wp
    end subroutine initialize
    !-
    subroutine alloc_dr (x_)
      type(t_atm), intent(inout) :: x_
      call construct (x_, g)
      call allocate  (x_, "div vrt")
    end subroutine alloc_dr
    !-
    subroutine alloc_uv (x_)
      type(t_atm), intent(inout) :: x_
      call construct (x_, g)
      call allocate  (x_, "u v")
    end subroutine alloc_uv
    !-
    subroutine apply_a (x_, y_, refatm, ph_ref)
      type(t_atm),            intent(in)    :: x_
      type(t_atm),            intent(inout) :: y_
      type(t_atm),            intent(in)    :: refatm
      real(wp)   , optional,  intent(in)    :: ph_ref(:)

      real(wp) :: p0ref = 101325._wp
      integer  :: k, ke

      ke = refatm% grid% nz
      !print *, 'prima max e min di wk u e v: ', dace% pe, dace% pio, maxval(wk% u), minval(wk% v)
      if (associated (x_% u)) then
         do k = 1,ke
            wk% u(:,:,k,:) = x_% u(:,:,k,:)                              &
                           * (refatm% ph(:,:,k+1,:) - refatm% ph(:,:,k,:))
            wk% v(:,:,k,:) = x_% v(:,:,k,:)                              &
                           * (refatm% ph(:,:,k+1,:) - refatm% ph(:,:,k,:))
         end do
         !print *, 'dopo max e min di wk u e v: ', dace% pe, dace% pio, maxval(wk% u), minval(wk% v)
         call hor_div  (wk% u, wk% v, y_% div, g)
         call hor_curl (x_% u, x_% v, wk% vrt, g)
         if (present (ph_ref)) then
            do k = 1, ke
               y_% vrt(:,:,k,:) = wk% vrt(:,:,k,:) * (ph_ref(k+1) - ph_ref(k))
            end do
         else
            y_% vrt = wk% vrt * (p0ref/ke)
         end if
      else
         call finish ("apply_a","inconsistent arguments")
      end if
    end subroutine apply_a
    !-
    subroutine apply_a_t (x_, y_, refatm, ph_ref)
      type(t_atm),            intent(in)    :: x_
      type(t_atm),            intent(inout) :: y_
      type(t_atm),            intent(in)    :: refatm
      real(wp)   , optional,  intent(in)    :: ph_ref(:)

      real(wp) :: p0ref = 101325._wp
      integer  :: k, ke

      ke = refatm% grid% nz
      if (associated (y_% u)) then
         if (present (ph_ref)) then
            do k = 1, ke
               wk% vrt(:,:,k,:) = x_% vrt(:,:,k,:) * (ph_ref(k+1) - ph_ref(k))
            end do
         else
            wk% vrt = x_% vrt * (p0ref/ke)
         end if
         y_  = 0._wp
         call hor_div_t  (y_% u, y_% v, x_% div, g)
         do k = 1, ke
            y_% u(:,:,k,:) = y_% u(:,:,k,:)                              &
                           * (refatm% ph(:,:,k+1,:) - refatm% ph(:,:,k,:))
            y_% v(:,:,k,:) = y_% v(:,:,k,:)                              &
                           * (refatm% ph(:,:,k+1,:) - refatm% ph(:,:,k,:))
         end do
         call hor_curl_t (y_% u, y_% v, wk% vrt, g)
      else
         call finish ("apply_a_t","inconsistent arguments")
      end if
    end subroutine apply_a_t
    !-
    function dot_prod (a, b)
      type(t_atm), intent(in) :: a, b
      real(wp)                :: dot_prod

      dot_prod = 0._wp
      if (associated (a% u)) then
         if (.not.(                        associated (b% u) .and. &
                   associated (a% v) .and. associated (b% v))      ) then
            write(0,*) "dot_prod: associated (u,v)=",  &
                 associated (a% u), associated (b% u), &
                 associated (a% v), associated (b% v)
            call finish ("dot_prod","inconsistent arguments")
         end if
         dot_prod = dot_prod + sum (a% u * b% u) + sum (a% v * b% v)
      end if
      if (associated (a% div)) then
         if (.not.(                          associated (b% div) .and. &
                   associated (a% vrt) .and. associated (b% vrt))      ) then
            write(0,*) "dot_prod: associated (div,vrt)=",  &
                 associated (a% div), associated (b% div), &
                 associated (a% vrt), associated (b% vrt)
            call finish ("dot_prod","inconsistent arguments")
         end if
         dot_prod = dot_prod + sum (a% div * b% div) + sum (a% vrt * b% vrt)
      end if
      dot_prod = p_sum (dot_prod)
    end function dot_prod
    !-
  end subroutine adjust_wind

!==============================================================================

  subroutine read_nml_postmult ()
    !-------------------------
    ! Read namelist /POSTMULT/
    !-------------------------
    integer :: ierr
    logical :: init = .false.

    namelist /POSTMULT/ lkeep_tinv, lfix_geosp, zpbl1, zpbl2,      &
                        wind_adjust, diag_level, diag_base,        &
                        adjust_det, adjust_ens, adjust_mode,       &
                        w_spread_an, sqz_trans,                    &
                        solver, min_iter, max_iter, tol

    if (init) return
    init = .true.

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') ' Namelist /POSTMULT/'
      write(6,'()')
    endif
    !-------------
    ! set defaults
    !-------------
    lkeep_tinv = .false.                ! Preserve temp.inversions
    lfix_geosp = .false.                ! Debugging (do not use!)
    zpbl1      =   500._wp              ! AGL heights used in
    zpbl2      =  1000._wp              ! lapse rate computations
    wind_adjust = 0                     ! Wind adjustment mode (0=off)
    diag_level  = 1                     ! Level of diagnostics
    adjust_det  = .true.                ! Adjust deterministic analysis
    adjust_ens  = .false.               ! Adjust ensemble      members
    adjust_mode = 0                     ! Source of dpsdt for adjustment
    w_spread_an = 1.0_wp                ! Rel. weight ana/fg spread
    sqz_trans   = [ 0.5_wp, 0.8_wp ]    ! Mod. for squeezed levels
    solver      = "cgls"                ! Solver type
    min_iter    = 1                     ! Min. no. iterations
    max_iter    = 20                    ! Max. no. iterations
    tol(1)      = 0.1_wp                ! Convergence criteria: shrinkage
    tol(2)      = 0._wp                 ! Convergence criteria: other
    diag_base   = "diag"                ! Basename of diagnostics file
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('POSTMULT', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=POSTMULT, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_postmult','ERROR in namelist /POSTMULT/')
#else
        read (nnml ,nml=POSTMULT)
#endif
      end select
    end if

    if (dace% lpio) then
      write(6,'(A,L10)'    ) '  lkeep_tinv  =',lkeep_tinv
      write(6,'(A,L10)'    ) '  lfix_geosp  =',lfix_geosp
      write(6,'(A,F10.2,A)') '  zpbl1       =',zpbl1, ' m'
      write(6,'(A,F10.2,A)') '  zpbl2       =',zpbl2, ' m'
      write(6,'(A,I10)'    ) '  wind_adjust =',wind_adjust
      write(6,'(A,I10)'    ) '  diag_level  =',diag_level
      write(6,'(A,L10)'    ) '  adjust_det  =',adjust_det
      write(6,'(A,L10)'    ) '  adjust_ens  =',adjust_ens
      write(6,'(A,I10)'    ) '  adjust_mode =',adjust_mode
      write(6,'(A,F10.6)'  ) '  w_spread_an =',w_spread_an
      write(6,'(A,2F10.6)' ) '  sqz_trans   =',sqz_trans
      write(6,'(A,1X,A)'   ) '  solver      =',trim (solver)
      write(6,'(A,I10)'    ) '  min_iter    =',min_iter
      write(6,'(A,I10)'    ) '  max_iter    =',max_iter
      write(6,'(A,9F10.6)' ) '  tol         =',tol
      write(6,'(A,1X,A)'   ) '  diag_base   =',trim (diag_base)
      write(6,'()')
    end if
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (lkeep_tinv,  dace% pio)
    call p_bcast (lfix_geosp,  dace% pio)
    call p_bcast (zpbl1,       dace% pio)
    call p_bcast (zpbl2,       dace% pio)
    call p_bcast (wind_adjust, dace% pio)
    call p_bcast (diag_level,  dace% pio)
    call p_bcast (adjust_det,  dace% pio)
    call p_bcast (adjust_ens,  dace% pio)
    call p_bcast (adjust_mode, dace% pio)
    call p_bcast (w_spread_an, dace% pio)
    call p_bcast (sqz_trans,   dace% pio)
    call p_bcast (solver,      dace% pio)
    call p_bcast (min_iter,    dace% pio)
    call p_bcast (max_iter,    dace% pio)
    call p_bcast (tol,         dace% pio)
    call p_bcast (diag_base,   dace% pio)

    if (wind_adjust < 0 .or. wind_adjust > 4) &
         call finish ("read_nml_postmult","wind_adjust out of range 0..4")
    if (sqz_trans(1) >= sqz_trans(2)) &
         call finish ("read_nml_postmult","sqz_trans(1) >= sqz_trans(2)")
  end subroutine read_nml_postmult

!==============================================================================
end module mo_postmult
