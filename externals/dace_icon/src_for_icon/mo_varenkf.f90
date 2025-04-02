!
!+ Use the ensemble background in the variational scheme
!
MODULE mo_varenkf
!
! Description:
!   Use the ensemble background in the variational scheme.
!     Operators: Hi Hp X Cv Tr Ch Cht Trt Cvt Xt Hpt Hit
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
!  new module for VARenKF
! V1_37        2014-12-23 Andreas Rhodin
!  new namelist /VARENKF/
! V1_42        2015-06-08 Andreas Rhodin
!  optimisation for PSAS; vertical wavelet transform; allow w_ens_b+w_nmc_b/=1
! V1_45        2015-12-15 Harald Anlauf
!  use OpenMP workshare; use non-blocking communication; fix for TR15581
! V1_47        2016-06-06 Andreas Rhodin
!  some fixes for no.vertical gridpoints == 1 (no vertical localisation)
! V1_51        2017-02-24 Andreas Rhodin
!  namelist /VARENKF/: generalise flags l_fgerr,l_precond to integer
!
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2013
!==============================================================================
#include "tr15581.incf"
  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,       only: wp,           &! working precision kind parameter
                           i8             ! 8-byte integer
  use mo_exception,  only: finish         ! abort in case of error
  use mo_mpi_dace,   only: dace,         &! MPI group info
                           p_bcast        ! generic MPI broadcast routine
  use mo_wmo_tables, only: WMO3_ISOBARIC,&! grid type specification
                           DWD6_NONE      ! fake grid type for set of columns
  use mo_namelist,   only: position_nml, &! routine to position nml group
                           nnml,         &! namelist fortran unit number
                           POSITIONED     ! position_nml: OK return flag
  !---------------------
  ! model state and grid
  !---------------------
  use mo_atm_state,  only: t_atm,        &! atmospheric state derived type
                           construct,    &! derived type constructor routine
                           destruct,     &! derived type destructor  routine
                           allocate,     &! allocate   selected components
                           deallocate,   &! deallocate selected components
                           set_pointers, &! set pointers to array 'm'
                           operator(-),  &! atmosphere - atmosphere
                           operator(/),  &! atmosphere / number
!                          operator(+),  &! atmosphere + atmosphere
!                          operator(**), &! atmosphere ** number
!                          sqrt,         &! sqrt(atmosphere)
!                          print,        &! print atmospheric state
                           assignment(=),&! atmosphere = atmosphere
                           set_tv,       &! derive virtual temperature
                           set_rh,       &! derive relative humidity
                           set_p,        &! set pressure on full levels
                           set_geo,      &! set geopotential height
                           vert_intp,    &! vertical intp. to p.-levels
                           select_params  ! copy selected fields
  use mo_atm_grid,   only: t_grid,       &! atmospheric grid  derived type
                           construct,    &! derived type constructor routine
                           destruct       ! derived type destructor  routine
  use mo_physics,    only: gacc           ! gravity acceleration
  use mo_t_col,      only: t_cols,       &! data type to hold a bunch of columns
                           atm2col,      &! gather columns from atmosphere
!                          gather_cols,  &! gather columns from other processors
                           gather_cols_multi,&! gather columns from other processors
                           dealloc_cols, &! deallocate components of datatype t_cols
                           COL_P, COL_UV, COL_RH, COL_GEO, COL_TV2
  use mo_t_obs,      only: t_mcol,       &! model column index derived type
                           t_mcols        ! column index array container
  use mo_obs_sndrcv, only: p_mcol,       &! Pointers to model column descriptors
                           scatter_mcol, &! scatter model columns
                           alltoall_mcol  ! scatter model columns (alltoall)
  !-----------------------
  ! interpolation operator
  !-----------------------
  use mo_bg_err_io,  only: t_intop,      &! interpolation operator derived type
                           destruct       ! destructor  for t_intop
  !-----------------------
  ! localisation operators
  !-----------------------
  use mo_1dmra,      only: TR_ANA,       &! Analysis   (x -> w)
                           TR_SYN,       &! Synthesis  (x <- w)
                           TR_ADJ         ! Adjoint    (x -> w)
  use mo_lift_1d,    only: lift_1d        ! 1d wavelet transform
  use mo_varenkf_1d, only: t_C1,         &! vertical localisation operator
                           setup_c1,     &! set vert.loc.operator meta data
                           apply_vloc,   &! apply sqrt vert. lovalisation
                           apply_vloc_t, &! apply adj. sqrt vert. localisation
                           apply_vrint,  &! vertical (reverse) interpolation
                           apply_vint,   &! vertical interpolation
                           apply_vint_t, &! adjoint vertical interpolation
                           construct,    &! set vert.loc.operator meta data
                           destruct       ! deallocate t_C1 components
  use mo_varenkf_2d, only: t_C2,         &! horizontal localisation operator
                           construct,    &! set hor.loc.operator meta data
                           apply_c2,     &! apply sqrt 2D localisation operator
                           apply_c2t,    &! apply adjoint sqrt 2D loc. operator
                           destruct       ! deallocate t_C2 components
  use mo_random,     only: random_gauss   ! random number generator
  implicit none

  !================
  ! public entities
  !================
  private
  public :: t_Benkf          ! variational ensemble B matrix meta data type
  public :: Benkf            ! variational ensemble B matrix meta data
  public :: setup_varenkf    ! set up variational ensemble B matrix meta data
  public :: setup_C          ! re-initialise localisation parameters
  public :: varenkf_psas_opt ! re-distribute model grid columns of ensemble
  public :: v_lift_analysis  ! vertical wavelet analysis
  public :: v_lift_synth     ! vertical wavelet synthesis
  public :: v_lift_adj       ! adjoint vertical wavelet synthesis
  public :: apply_varenkf    ! apply B_enkf to atmospheric state vector
  public :: test_varenkf     ! simple test of VarEnKF B matrix operator
  public :: destruct         ! deallocate derived type components
  public :: read_nml_varenkf ! read vanelist /VARENKF/
  public :: w_ens_b          ! weight of ensemble B
  public :: w_nmc_b          ! weight of NMC      B (1-w_ens_b)
  public :: w_nmc_a          ! weight of NMC      B in anal.error calculation
  public :: i_fgerr          ! use ens_b for fg error
  public :: i_precond        ! use ens_b for preconditioner
  public :: l_psasopt        ! re-distribute model columns for PSAS HBH
  public :: l_anaerr         ! analysis error: use combined B_nmc/B_ens
  public :: n_vert           ! number of gridpoints for vertical wavelets
  public :: c_scale          ! vertical inter scale correlation parameter
  public :: c_vert           ! vertical correlation parameter
  public :: w_vscale         ! vertical wavelet localisation normalisation
  public :: set_vscale       ! set 'w_vscale' normalisation factors
  !===========
  ! Interfaces
  !===========
  interface destruct
    module procedure destruct_varenkf
  end interface destruct
!------------------------------------------------------------------------------
  !=========================
  ! Derived type definitions
  !=========================

  type t_Benkf
    !----------
    ! meta data
    !----------
    integer               :: n_ens            ! ensemble size
    !---------------------------------
    ! reference to background ensemble
    !---------------------------------
    type(t_grid) ,pointer :: grid             ! ensemble model grid
    type(t_grid) ,pointer :: grid0            ! grid for PSAS optimisation
    type (t_atm) ,pointer :: xb     (:)       ! background ensemble
    !-----------------------------------------------------------------
    ! ensemble deviations and mean (selected 3dvar control parameters)
    !-----------------------------------------------------------------
    integer               :: l_set            ! set to be used for PSAS
    integer               :: l_init           ! set actually initialised
    type (t_atm) ,pointer :: X      (:,:)     ! ensemble deviations
    type (t_atm)          :: mean   (0:1)     ! mean background
    !----------------------
    ! vertical localisation
    !----------------------
    type (t_C1)           :: vloc             ! vertical localisation data
    integer               :: k_vert           ! no.gridpoints in physical space
    integer               :: n_vert           ! no.gridpoints for wavelet
    integer               :: m_vert           ! no.gridpoints for frame
    integer               :: j_vert           ! total no. verticalgridpoints
    integer               :: l_vert           ! no. localisation gridpoints
    !------------------------
    ! horizontal localisation
    !------------------------
    type (t_C2)           :: hloc             ! horizontal localisation data
    !---------------------------------------
    ! interpolation to observation locations
    !---------------------------------------
    type (t_intop)        :: io               ! interpolation operator
  end type t_Benkf

!------------------------------------------------------------------------------
  !=================
  ! Module variables
  !=================

  type(t_Benkf) ,save :: Benkf          ! VarEnKF B matrix meta data

  real(wp) ,allocatable :: w_vscale (:) ! vertical scaling factors
  real(wp) ,allocatable :: a_vscale (:) ! vertical scaling factors

  !===================
  ! namelist /VARENKF/
  !===================
  !-----------------------------------
  ! partitioning of NMC and ensemble B
  !-----------------------------------
  real(wp)            :: w_ens_b   = 0._wp   ! weight of ensemble B
  real(wp)            :: w_nmc_b   = 1._wp   ! weight of NMC B (1-w_ens_b)
  real(wp)            :: w_nmc_a   = 1._wp   ! weight in anal.error calc.
  !--------------------------
  ! application of ensemble B
  !--------------------------
  real(wp)            :: lv_surf   = 0._wp   ! vert.localization scale, surface
  real(wp)            :: lv_top    = 0._wp   ! vert.localization scale, top
  logical             :: l_fgerr   = .false. ! depreciated
  logical             :: l_precond = .false. ! depreciated
  logical             :: l_psasopt = .true.  ! re-distribute model columns
  logical             :: l_anaerr  = .true.  ! analysis error: use w_nmc_b
  integer             :: i_fgerr   = 0       ! use ens_b for fg error
  integer             :: i_precond = 0       ! use ens_b for preconditioner
                                             !   values for i_fgerr, i_precond:
                                             !   0: do not use
  !--------------------------                !   1: localise in intp. space
  ! vertical wavelet approach                !   2: localise in obs.  space
  !--------------------------
  integer             :: n_vert    = 0       ! number of gridpoints
  real(wp)            :: c_scale   = 0.5_wp  ! inter scale correlation parameter
  real(wp)            :: c_vert    = 1.0_wp  ! vertical correlation parameter

  namelist /VARENKF/ w_ens_b, w_nmc_b,       &
                     lv_surf, lv_top,        &
                     l_fgerr, l_precond,     &
                     i_fgerr, i_precond,     &
                     l_psasopt, l_anaerr,    &
                     n_vert, c_scale, c_vert,&
                     w_nmc_a
!==============================================================================
contains
!==============================================================================

  subroutine read_nml_varenkf
  !------------------------
  ! read namelist /VARENKF/
  !------------------------
    integer :: ierr
    !-------------
    ! set defaults
    !-------------
    w_ens_b   =  0._wp  ! weight of ensemble B
    w_nmc_b   = -1._wp  ! weight of NMC      B (1-w_ens_b)
    w_nmc_a   = -1._wp  ! weight of NMC      B (analysis error calculation)
    lv_surf   =  0._wp  ! vert.localization scale, surface
    lv_top    =  0._wp  ! vert.localization scale, top
    l_fgerr   = .true.  ! depreciated
    l_precond = .true.  ! depreciated
    i_fgerr   = 1       ! use ens_b for fg error
    i_precond = 1       ! use ens_b for preconditioner
    l_psasopt = .true.  ! re-distribute model columns for PSAS HBH
    l_anaerr  = .true.  ! use w_nmc_b for analysis error
    n_vert    = 0       ! number of gridpoints
    c_scale   = 0.5_wp  ! inter scale correlation parameter
    c_vert    = 1.0_wp  ! vertical correlation parameter
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      call position_nml ('VARENKF', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=VARENKF, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_varenkf',          &
                                  'ERROR in namelist /VARENKF/')
#else
        read (nnml ,nml=VARENKF)
#endif
      end select
      !-------------------
      ! consistency checks
      !-------------------
      w_ens_b = max (0._wp, w_ens_b)
      if (w_nmc_b <  0._wp) then
        w_nmc_b = min (1._wp, max (0._wp, 1._wp - w_ens_b))
      endif
      if (w_nmc_a <  0._wp) then
        w_nmc_a = 1._wp; if (l_anaerr) w_nmc_a = w_nmc_b
      else
        w_nmc_a = min (w_nmc_a, 1._wp)
      end if
      if (w_ens_b == 0._wp) then
        l_fgerr   = .false.
        l_precond = .false.
        l_psasopt = .false.
      endif
      if (.not. l_fgerr  ) i_fgerr   = 0
      if (.not. l_precond) i_precond = 0
      !---------
      ! printout
      !---------
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') 'read namelist /VARENKF/'
      write(6,'()')
      write(6,*)  '  w_ens_b     =', w_ens_b
      write(6,*)  '  w_nmc_b     =', w_nmc_b
      write(6,*)  '  w_nmc_a     =', w_nmc_a
      write(6,'()')
      if (lv_surf > 0._wp) then
        write(6,*)'  lv_surf     =', lv_surf
      else
        write(6,*)'  lv_surf     =', " (taken from LETKF)"
      end if
      if (lv_top > 0._wp) then
        write(6,*)'  lv_top      =', lv_top
      else
        write(6,*)'  lv_top      =', " (taken from LETKF)"
      end if
      write(6,'()')
      write(6,*)  '  i_fgerr     =', i_fgerr
      write(6,*)  '  i_precond   =', i_precond
      write(6,*)  '  l_psasopt   =', l_psasopt
      write(6,*)  '  l_anaerr    =', l_anaerr
      write(6,'()')
      write(6,*)  '  n_vert      =', n_vert
      write(6,*)  '  c_scale     =', c_scale
      write(6,*)  '  c_vert      =', c_vert
      write(6,'()')
    endif
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (w_ens_b   ,dace% pio)
    call p_bcast (w_nmc_b   ,dace% pio)
    call p_bcast (w_nmc_a   ,dace% pio)
    call p_bcast (lv_surf   ,dace% pio)
    call p_bcast (lv_top    ,dace% pio)
    call p_bcast (i_fgerr   ,dace% pio)
    call p_bcast (i_precond ,dace% pio)
    call p_bcast (l_psasopt ,dace% pio)
    call p_bcast (l_anaerr  ,dace% pio)
    call p_bcast (n_vert    ,dace% pio)
    call p_bcast (c_scale   ,dace% pio)
    call p_bcast (c_vert    ,dace% pio)
  end subroutine read_nml_varenkf

!------------------------------------------------------------------------------

  subroutine destruct_varenkf (Benkf)
  type (t_Benkf) ,intent(inout) :: Benkf   ! VarEnKF B matrix meta data
  !----------------------------------------------
  ! deallocate components of derived type t_Benkf
  !----------------------------------------------
    call destruct (Benkf% grid)
    call destruct (Benkf% X)
    call destruct (Benkf% mean)
    call destruct (Benkf% hloc)
    call destruct (Benkf% vloc)
    call destruct (Benkf% io)
    deallocate    (Benkf% grid)
    deallocate    (Benkf% X)
    nullify       (Benkf% xb)
    if (associated  (Benkf% grid0)) then
      call destruct (Benkf% grid0)
      deallocate    (Benkf% grid0)
    endif
  end subroutine destruct_varenkf

!------------------------------------------------------------------------------

  subroutine setup_varenkf (Benkf, xb, mean, hloc, lv_surf_ref, lv_top_ref)
  !-------------------------------------------------------------
  ! set up the localised variational ensemble B matrix operators
  !-------------------------------------------------------------
  type(t_Benkf) ,intent(out)   :: Benkf       ! VarEnKF B matrix meta data
  type(t_atm)   ,pointer       :: xb(:)       ! ensemble deviations
  type(t_atm)   ,intent(inout) :: mean        ! ensemble mean
  real(wp)      ,intent(in)    :: hloc        ! horizontal localisation scale (m)
  real(wp)      ,intent(in)    :: lv_surf_ref ! vertical loc. scale at surface
  real(wp)      ,intent(in)    :: lv_top_ref  ! vertical loc. scale at model top
  optional                     :: lv_top_ref,&
                                  lv_surf_ref

    logical     :: ltv ! tv   is present
    logical     :: lrh ! rh   is present
    logical     :: lpf ! pf   is present
    logical     :: lph ! ph   is present
    logical     :: lgh ! geoh is present
    logical     :: lgf ! geof is present
    integer     :: i   ! ensemble member index
    type(t_atm) :: tmp ! temporary

    if (present (lv_surf_ref) .and. lv_surf <= 0._wp) lv_surf = lv_surf_ref
    if (present (lv_top_ref)  .and. lv_top  <= 0._wp) lv_top  = lv_top_ref
    if (lv_surf <= 0._wp .or. lv_top <= 0._wp .or. lv_surf > lv_top)        &
         call finish ("setup_varenkf","lv_surf, lv_top must be properly set")
    !----------------------------------------
    ! store references to ensemble deviations
    !----------------------------------------
    Benkf% n_ens  =  size (xb)
    Benkf% xb     => xb
    Benkf% l_set  =  1
    Benkf% l_init = -1
    !-------------------------------------
    ! check if required fields are present
    !-------------------------------------
!   if (.not.associated (mean% pf)) call finish ('setup_varenkf','pf missing')
    if (.not.associated (mean% t )) call finish ('setup_varenkf','t  missing')
    if (.not.associated (mean% q )) call finish ('setup_varenkf','q  missing')
    if (.not.associated (mean% u )) call finish ('setup_varenkf','u  missing')
    if (.not.associated (mean% v )) call finish ('setup_varenkf','v  missing')
    ltv = associated (mean% tv)
    lrh = associated (mean% rh)
    lpf = associated (mean% pf)
    lph = associated (mean% ph)
    lgh = associated (mean% geoh)
    lgf = associated (mean% geof)
    !----------------------------
    ! calculate control variables
    !----------------------------
    if (.not.lpf.or..not.lph) call set_p   (mean)
    if (.not.lpf.or..not.lph) call set_p   (xb)
    call set_tv (mean)
    call set_rh (mean)
    call set_tv (xb)
    call set_rh (xb)
    if (.not.lgf)             call set_geo (mean ,geof=.true.)
    if (.not.lgf)             call set_geo (xb   ,geof=.true.)
    !---------------------------------------------------------
    ! define grid for deviations as pressure coordinate system
    !---------------------------------------------------------
    nullify  (Benkf% grid0)
    allocate (Benkf% grid)
    call construct (Benkf% grid, mean% grid, levtyp=WMO3_ISOBARIC)
    !--------------------------------------------------------
    ! copy ensemble & mean (selected 3dvar control variables)
    !--------------------------------------------------------
    allocate           (Benkf% X      (Benkf% n_ens, 0:1)             )
    call construct     (Benkf% X(:,1), Benkf% grid, name='bg_dev'     )
    call select_params (Benkf% X(:,1), xb,         'u v tv rh pf geof')
    call construct     (Benkf% mean(1),Benkf% grid, name='mean'       )
    call select_params (Benkf% mean(1),mean,       'u v tv rh pf geof')
    !-------------------
    ! allocate temporary
    !-------------------
    call construct (tmp, Benkf% grid)
    call allocate (tmp, 'pf')
    !------------------
    ! loop over members
    !------------------
    do i = 1, Benkf% n_ens
      !--------------------------------------------------------
      ! interpolate to common (mean) pressure coordinate system
      !--------------------------------------------------------
      tmp% pf = mean% pf
!     Benkf% X(i)% m% i% ref = .false.
!     where (Benkf% X(i)% m% i% name == 'pf') Benkf% X(i)% m% i% ref = .true.
      call vert_intp (Benkf% X(i,1), tmp, inplace=.true.)
      Benkf% X(i,1) = tmp
      Benkf% X(i,1) = (Benkf% X(i,1) - Benkf% mean(1)) &
                    / sqrt (Benkf% n_ens - 1._wp)
      Benkf% X(i,1)% geof = Benkf% X(i,1)% geof / gacc   ! m^2/s^2 -> m
    end do

!call print (Benkf% X(1),comment='nachher')

!call print (Benkf% mean, comment='Benkf% mean')

!do i = 1, Benkf% n_ens
!  Benkf% X(i)% pf   = 0._wp
!  Benkf% X(i)% u    = 0._wp
!  Benkf% X(i)% v    = 0._wp
!  Benkf% X(i)% tv   = sqrt (1._wp / Benkf% n_ens) * (-1)**i
!  Benkf% X(i)% rh   = 0._wp
!  Benkf% X(i)% geof = 0._wp
!
!print *,dace% pe,'### Benkf%  X(i)% rh',i,(-1)**i, sqrt (1._wp / Benkf% n_ens) * (-1)**i
!
!end do
!
!   tmp = 0._wp
!   do i = 1, Benkf% n_ens
!     tmp = tmp + Benkf% X(i) ** 2
!   end do
!   tmp = sqrt(tmp)
!
!call print (tmp, comment='varenkf spread', verbose=.true.)
!call print (Benkf% X(1), comment='X(1)', verbose=.true.)
!call print (Benkf% X(2), comment='X(2)', verbose=.true.)

    !------------------------------
    ! set up localisation operators
    !------------------------------
    call setup_C (Benkf, 1, hloc, lv_surf, lv_top)
    !---------------------------------------
    ! set up vertical wavelet transformation
    !---------------------------------------
    Benkf% k_vert = Benkf% grid% nz
    Benkf% n_vert = n_vert
    Benkf% m_vert = n_vert * 2
    Benkf% j_vert = Benkf% k_vert + Benkf% m_vert
    if (Benkf% k_vert == 1) then
      Benkf% l_vert = 1
    else if (n_vert>0) then
      Benkf% l_vert = Benkf% m_vert
      call set_vscale
!     call v_lift_analysis
    else
      Benkf% l_vert = Benkf% vloc% nzr
    endif

    !---------
    ! clean up
    !---------
    call destruct (tmp)
    if (.not.ltv) call deallocate (mean, 'tv')
    if (.not.ltv) call deallocate (xb,   'tv')
    if (.not.lrh) call deallocate (mean, 'rh')
    if (.not.lrh) call deallocate (xb,   'rh')
    if (.not.lpf) call deallocate (mean, 'pf')
    if (.not.lpf) call deallocate (xb,   'pf')
    if (.not.lph) call deallocate (mean, 'ph')
    if (.not.lph) call deallocate (xb,   'ph')
    if (.not.lgh) call deallocate (mean, 'geoh')
    if (.not.lgh) call deallocate (xb,   'geoh')
    if (.not.lgf) call deallocate (mean, 'geof')
    if (.not.lgf) call deallocate (xb,   'geof')

  end subroutine setup_varenkf

!------------------------------------------------------------------------------

  subroutine setup_C (Benkf, l_init, hloc, lv_surf, lv_top)
  type(t_Benkf)      ,intent(inout) :: Benkf   ! VarEnKF B matrix meta data
  integer            ,intent(in)    :: l_init  ! ensemble mean
  real(wp) ,optional ,intent(in)    :: hloc    ! horizontal localisation scale
  real(wp) ,optional ,intent(in)    :: lv_top  ! vertical loc. scale at top
  real(wp) ,optional ,intent(in)    :: lv_surf ! vertical loc. scale at bottom
  !---------------------------------------------------------------------------
  ! set localisation parameters in 'Benkf'
  !
  ! If 'hloc', 'lv_top', 'lv_surf' are not given, the old parameters are kept,
  ! merely interpolation coefficients are adapted to re-distribution
  ! of grid-points given by 'Benkf% mean'
  !---------------------------------------------------------------------------

    target               :: Benkf
    real(wp)             :: r_loc
    type(t_atm), pointer :: mean
    type(t_grid),pointer :: refgrid

    !----------------------------------
    ! switch to required set of columns
    !----------------------------------
    if (l_init == Benkf% l_init) return
    Benkf% l_init =  l_init
    mean          => Benkf% mean (l_init)
    !--------------------------------------
    ! set up vertical localisation operator
    !--------------------------------------
    if (Benkf% grid% nz > 1) then
      call setup_c1 (Benkf% vloc, mean, lv_surf, lv_top, n_vert)
    endif
    !----------------------------------------
    ! set up horizontal localisation operator
    !----------------------------------------
    if (present (hloc)) then
      call construct (Benkf% hloc, mean% grid, hloc)
    else
      !----------------------------------------------------
      ! If "mean" has no real grid (collection of columns),
      ! pass the reference grid pointed to in Benkf% grid.
      !----------------------------------------------------
      nullify (refgrid)
      if (mean% grid% gridtype == DWD6_NONE) refgrid => Benkf% grid

      r_loc = Benkf% hloc% r_loc
      call destruct  (Benkf% hloc)
      call construct (Benkf% hloc, mean% grid, r_loc, refgrid=refgrid)
    endif

  end subroutine setup_C

!------------------------------------------------------------------------------

  subroutine varenkf_psas_opt (Benkf)
  type(t_benkf) ,intent(inout) ,target :: Benkf
  !----------------------------------------------------------------
  ! re-distribute model grid columns of VarEnKF ensemble (and mean)
  ! so that they match the distribution of observations.
  ! To be used for optimisation in the VarEnKF PSAS step.
  !
  ! 1) gather the required model columns on the PEs,
  !    as required for observation processing
  ! 2) set up grid meta data for the set of model columns
  ! 3) store columns in derived type t_atm
  ! 4) wavelet transform the original set
  ! 5) adjust the transposition and localisation meta data
  ! 6) wavelet transform the optimised set
  !----------------------------------------------------------------

!   type p_mcol
!     type (t_mcol) ,pointer :: p (:)
!   end type p_mcol
    !----------------
    ! local variables
    !----------------
    integer               :: i, j, k, d, l, m ! indices
    integer(i8)           :: latm         ! parameters required
    integer               :: nmcol        ! number of columns required
    integer               :: lb (4)       ! lbound    for selected columns
    integer               :: ub (4)       ! ubound    for selected columns
    real(wp) ,allocatable :: rlon (:,:,:) ! longitudes of selected columns
    real(wp) ,allocatable :: rlat (:,:,:) ! latitudes  of selected columns
    type(t_atm)  ,pointer :: x, y         ! pointer to source and dest.

    type (p_mcol)         :: mcols    (               dace% npe)  ! temporary
    type (t_cols)         :: cols_send(               dace% npe)
    type (t_cols)         :: cols_recv(0:dace% npe-1, dace% npe)
    type (t_mcols)        :: mc1      (               dace% npe)

    !------------------------------------
    ! only if optimisation is switched on
    !------------------------------------
    if (l_psasopt) then

    !------------------------------------------------------
    ! 0) set up grid meta data for the set of model columns
    !------------------------------------------------------
    nmcol = Benkf% io% mc% n
    lb    = 1
    ub    = (/nmcol, 1, Benkf% grid% nz, 1/)
    allocate (Benkf% grid0)
    allocate (rlon (nmcol,1,1))
    allocate (rlat (nmcol,1,1))
    do l = 1, nmcol
      i = Benkf% io% mc% c(l)% ijdtp (1)
      j = Benkf% io% mc% c(l)% ijdtp (2)
      d = Benkf% io% mc% c(l)% ijdtp (3)
      rlon (l,1,1) = Benkf% grid% rlon (i,j,1,d)
      rlat (l,1,1) = Benkf% grid% rlat (i,j,1,d)
    end do
    call construct (Benkf% grid0, gridtype=DWD6_NONE,               &
                    nx=ub(1), ny=ub(2), ke=ub(3), lon=rlon, lat=rlat)
!   Benkf% grid0% global = Benkf% grid% global
    !-----------------------------------------------------------------
    ! 1a) The meta data (column indices required for interpolation in
    !     the 'boxes' at the target PE) is gathered on the source PEs.
    !
    !     This is basically an alltoall communication pattern.
    !     Using alltoall communication directly may be more efficient.
    !-----------------------------------------------------------------
    mc1(dace% pe+1) = Benkf% io% mc       ! flat copy
    mc1% pe = (/(i,i=0,dace% npe-1)/)

#ifdef ALLTOALL_MCOL_ORIGINAL
    do i = 1, dace% npe
      nullify                   (mcols(i)%p)
      call scatter_mcol (mc1(i), mcols(i)%p, 1)
    end do
#else
    call alltoall_mcol  (mc1(:), mcols, tslot=1)
#endif
    !---------------------------------
    ! loop over mean, ensemble members
    !---------------------------------
    do m = 0, Benkf% n_ens
      if (m==0) then
        x => Benkf% mean(1)
        y => Benkf% mean(0)
      else
        x => Benkf% X(m,1)
        y => Benkf% X(m,0)
      endif
      !----------------------------------------------------------------
      ! 1b) On the source PE the required model columns are stored into
      !     an array of derived type t_cols for redistribution.
      !----------------------------------------------------------------
      latm = COL_P + COL_UV + COL_RH + COL_GEO + COL_TV2
      do i=1,dace% npe
        call atm2col (cols_send(i), x, mcols(i)%p, iatm=latm)
      end do
      !-----------------------------------------------------------------
      ! 1c) The model columns are transmitted to the target PEs
      !
      !     This is basically an alltoall communication pattern.
      !     Using alltoall communication directly may be more efficient.
      !-----------------------------------------------------------------
!!$      do i=1,dace% npe
!!$        call gather_cols (cols_send(i), cols_recv(:,i), root=mc1(i)% pe)
!!$      end do
      call gather_cols_multi (cols_send, cols_recv, root=mc1(:)% pe)
      !------------------------------------------------------
      ! 2) set up grid meta data for the set of model columns
      !------------------------------------------------------
      !    (moved out of loop over m)
      !---------------------------------------
      ! 3) store columns in derived type t_atm
      !---------------------------------------
      call construct (y, template=x, lb=lb, ub=ub, grid=Benkf% grid0)
      i = dace% pe + 1
      do j = 0, dace% npe-1
        l = 0
        do k=1,nmcol
          if (Benkf% io% mc% c(k)% ijdtp(5) == j) then
            l = l + 1
            y% pf   (k, 1, :, 1) = cols_recv(j,i)% col(l)% p
            y% geof (k, 1, :, 1) = cols_recv(j,i)% col(l)% geo
            y% tv   (k, 1, :, 1) = cols_recv(j,i)% col(l)% tv
            y% rh   (k, 1, :, 1) = cols_recv(j,i)% col(l)% rh
            y% u    (k, 1, :, 1) = cols_recv(j,i)% col(l)% u
            y% v    (k, 1, :, 1) = cols_recv(j,i)% col(l)% v
          endif
        end do
        if (j /= dace% pe) call dealloc_cols (cols_recv(j,i))
      end do
      call dealloc_cols (cols_send(i))
    end do

    !---------
    ! clean up
    !---------
    do i = 1, dace% npe
      deallocate (mcols(i)%p)
    end do

    endif

    !--------------------------------------
    ! 4) wavelet transform the original set
    !--------------------------------------
    call v_lift_analysis (Benkf, 1)

    if (.not. l_psasopt) return
    !-------------------------------------------------------
    ! 5) adjust the transposition and localisation meta data
    !-------------------------------------------------------
    Benkf% l_set = 0
    do i = 1, nmcol
      benkf% io% mc% c(i)% ijdtp = (/i,1,1,1,dace% pe/)
    end do
    call setup_c  (Benkf, 0)

    !---------------------------------------
    ! 6) wavelet transform the optimised set
    !---------------------------------------
    call v_lift_analysis (Benkf, 0)


  end subroutine varenkf_psas_opt

!------------------------------------------------------------------------------

  subroutine apply_varenkf (Benkf, a, z, li, lo)
  type(t_Benkf) ,intent(inout) :: Benkf  ! VarEnKF B matrix meta data
  type(t_atm)   ,intent(in)    :: z      ! gradient
  type(t_atm)   ,intent(inout) :: a      ! analysis increment
  integer       ,intent(in)    :: li     ! index for PSAS optimisation (in)
  integer       ,intent(in)    :: lo     ! index for PSAS optimisation (out)

    real(wp) ,allocatable :: t1 (:,:,:,:) ! temporary
    real(wp) ,allocatable :: t2 (:,:,:,:) ! temporary
    real(wp) ,allocatable :: t3 (:,:,:,:) ! temporary
    real(wp) ,allocatable :: t4 (:,:,:,:) ! temporary
    real(wp) ,allocatable :: t5 (:,:,:,:) ! temporary
    target                :: t1, t2, t4, t5
    real(wp) _POINTER     :: p4 (:,:,:,:)
    real(wp) _POINTER     :: p5 (:,:,:,:)
    integer               :: lb (4)       ! lower bound
    integer               :: ub (4)       ! upper bound
    integer               :: m            ! ensemble member index

    !---------------------
    ! allocate temporaries
    !---------------------
    !-----------------------------------------------
    ! model grid [optionally vertically transformed]
    !-----------------------------------------------
    lb     = Benkf% mean(li)% lb
    ub     = Benkf% mean(li)% ub
    ub (3) = Benkf% j_vert
    allocate (t1 (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
    !--------------------------
    ! vertical coefficient grid
    !--------------------------
    ub (3) = Benkf% l_vert
    allocate (t2 (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
    !---------------------------------------
    ! vertical + horizontal coefficient grid
    !---------------------------------------
    lb     = Benkf% hloc% gc% lb
    ub     = Benkf% hloc% gc% ub
    ub (3) = Benkf% l_vert
    allocate (t3 (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
    if (li == lo) then
      !---------------------
      ! adjoint on same grid
      !---------------------
      p4 => t2
      p5 => t1
    else
      !----------------------
      ! adjoint on model grid
      !-----------------------
      lb     = Benkf% mean(lo)% lb
      ub     = Benkf% mean(lo)% ub
      ub (3) = Benkf% j_vert
      allocate (t5 (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
      !--------------------------
      ! vertical coefficient grid
      !--------------------------
      ub (3) = Benkf% l_vert
      allocate (t4 (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))
      p4 => t4
      p5 => t5
    endif

    !---------------------------
    ! loop over ensemble members
    !---------------------------
    a    = 0._wp
    do m = 1, Benkf% n_ens
      !---------------------------------
      ! mutiply with ensemble deviations
      !---------------------------------
      call varenkf_X_t (t1, z, Benkf, m, li)
      !------------------------------------------
      ! check for correct localisation parameters
      !------------------------------------------
      call setup_C (Benkf, li)
      !----------------------
      ! vertical localisation
      !----------------------
      if (benkf% k_vert == 1) then
        t2 = t1
      else if (benkf% n_vert == 0) then  ! --->
        call apply_vloc_t (Benkf% vloc, t1, t2)
      else
        call apply_wloc_t (Benkf% vloc, t1, t2)
      endif
      !------------------------
      ! horizontal localisation
      !------------------------
      call apply_c2t (Benkf% hloc, t3, t2)
      !------------------------------------------
      ! check for correct localisation parameters
      !------------------------------------------
      call setup_C (Benkf, lo)
      !------------------------
      ! horizontal localisation
      !------------------------
      call apply_c2 (Benkf% hloc, p4, t3)
      !----------------------
      ! vertical localisation
      !----------------------
      if (benkf% k_vert == 1) then
        p5 = p4
      else if (benkf% n_vert == 0) then ! <---
        call apply_vloc (Benkf% vloc, p5, p4)
      else
        call apply_wloc (Benkf% vloc, p5, p4)
      endif
      !---------------------------------
      ! mutiply with ensemble deviations
      !---------------------------------
      call varenkf_X (p5, a, Benkf, m, lo)
    end do
    !---------
    ! clean up
    !---------

  end subroutine apply_varenkf

!------------------------------------------------------------------------------

  subroutine test_varenkf (Benkf)
  !+++++++++++++++++++++++++++++++++++++++++
  ! simple test of VarEnKF B matrix operator
  !+++++++++++++++++++++++++++++++++++++++++
  type (t_Benkf) ,intent(inout) :: Benkf   ! VarEnKF B matrix meta data

    type (t_atm) :: z, a ! gradient, analysis inctement
    integer      :: i    ! iteration count for run-time test
    integer      :: k
    !-----------------------------
    ! allocate argument and result
    !-----------------------------
    call construct (z, Benkf% grid, alloc = 'u v tv rh geof')
    call construct (a, Benkf% grid, alloc = 'u v tv rh geof')
    !-----------------
    ! pre-set argument
    !-----------------
    z = 0._wp
    k = (z% lb(3) + z% ub(3)) / 2
    if (1 >= z% lb(1) .and. 1 <= z% ub(1) .and. &
        1 >= z% lb(2) .and. 1 <= z% ub(2)       ) z% tv (1,1,k,1) = 1._wp
    !---------------------
    ! multiply with B-EnKF
    !---------------------
    do i = 1, 1
      call apply_varenkf (Benkf, a, z, 1, 1)
    end do
    !---------
    ! printout
    !---------
!   call print (z, comment = 'z (VarEnKF test input)')
!   call print (a, comment = 'a (VarEnKF test output)')
    !---------
    ! clean up
    !---------
    call destruct (a)
    call destruct (z)

  end subroutine test_varenkf

!------------------------------------------------------------------------------

  subroutine varenkf_X (w, d, Benkf, m, l)
  real(wp)       ,intent(in)    :: w (:,:,:,:) ! ensemble weights
  type(t_atm)    ,intent(inout) :: d           ! analysis increment
  type(t_Benkf)  ,intent(in)    :: Benkf       ! VarEnKF B matrix meta data
  integer        ,intent(in)    :: m           ! ensemble member index
  integer        ,intent(in)    :: l           ! set index
  !------------------------------------------------------------
  ! Part of the EnKF B matrix multiplication operator:
  !   Multiplication with the ensemble state
  !     output: analysis increment
  !       variables: wind component u
  !                  wind component v
  !                  temperature    t
  !                    or geopotential height (GME)
  !                    or pressure            (ICON)
  !                  relative humidity
  !     input:  weights of ensemble members
  !------------------------------------------------------------

#ifndef _CRAYFTN

!$omp parallel workshare
  d% u    = d% u    + Benkf% X(m,l)% u    * w
  d% v    = d% v    + Benkf% X(m,l)% v    * w
  d% tv   = d% tv   + Benkf% X(m,l)% tv   * w
  d% rh   = d% rh   + Benkf% X(m,l)% rh   * w
  d% geof = d% geof + Benkf% X(m,l)% geof * w
!$omp end parallel workshare

#else

!$omp parallel workshare
!DIR$ IVDEP
  d% u    = d% u    + Benkf% X(m,l)% u    * w
!DIR$ IVDEP
  d% v    = d% v    + Benkf% X(m,l)% v    * w
!DIR$ IVDEP
  d% tv   = d% tv   + Benkf% X(m,l)% tv   * w
!DIR$ IVDEP
  d% rh   = d% rh   + Benkf% X(m,l)% rh   * w
!DIR$ IVDEP
  d% geof = d% geof + Benkf% X(m,l)% geof * w
!$omp end parallel workshare

#endif

  end subroutine varenkf_X

!------------------------------------------------------------------------------

  subroutine varenkf_X_t (w, z, Benkf, m, l)
  real(wp)       ,intent(out) :: w (:,:,:,:) ! ensemble weights
  type(t_atm)    ,intent(in)  :: z           ! gradient
  type(t_Benkf)  ,intent(in)  :: Benkf       ! VarEnKF B matrix meta data
  integer        ,intent(in)  :: m           ! ensemble member index
  integer        ,intent(in)  :: l           ! set index

#ifndef __INTEL_COMPILER  /* ICE in Intel v15 */
!$omp parallel workshare
#endif

    w = z% u    * Benkf% X(m,l)% u  &
      + z% v    * Benkf% X(m,l)% v  &
      + z% tv   * Benkf% X(m,l)% tv &
      + z% rh   * Benkf% X(m,l)% rh &
      + z% geof * Benkf% X(m,l)% geof

#ifndef __INTEL_COMPILER
!$omp end parallel workshare
#endif

  end subroutine varenkf_X_t

!==========================================================================

  subroutine v_lift_analysis (Benkf, l)
  type(t_benkf) ,intent(inout) :: Benkf
  integer       ,intent(in)    :: l

    integer           :: i
    integer           :: j
    integer           :: d
    integer           :: m
    integer           :: k
    type (t_grid), pointer :: g
    real(wp) ,ALLOCATABLE  :: x (:,:,:,:)

    if (Benkf% n_vert == 0) return
    !---------------
    ! loop over sets
    !---------------
    select case (l)
    case (0)
      g => Benkf% grid0
    case (1)
      g => Benkf% grid
    end select
    !------------------
    ! loop over members
    !------------------
    do k = 1, size (Benkf% X,1)
      !--------------------
      ! loop over variables
      !--------------------
      do m = 1, size (Benkf% X(1,l)% m)
      if (.not. Benkf% X(k,l)% m(m)% i% alloc) cycle
        allocate (x (g% lb(1): g% ub(1), &
                     g% lb(2): g% ub(2), &
                     Benkf% j_vert     , &
                     g% lb(4): g% ub(4)) )
        !--------------------------------
        ! interpolation to regular p-grid
        !--------------------------------
        call apply_vrint (Benkf% vloc, Benkf% X(k,l)% m(m)% ptr, &
                          x(:,:,g% nz+1:g% nz+Benkf% n_vert,:))
        call apply_vint  (Benkf% vloc,                           &
                          x(:,:,1      :g% nz              ,:),  &
                          x(:,:,g% nz+1:g% nz+Benkf% n_vert,:)   )
        x                             (:,:,1:g% nz,:) = &
          Benkf% X(k,l)% m(m)% ptr - x(:,:,1:g% nz,:)
        !--------------------------------
        ! loop over horizontal gridpoints
        !--------------------------------
        do d = g% lb(4), g% ub(4)
        do j = g% lb(2), g% ub(2)
        do i = g% lb(1), g% ub(1)
          call lift_1d (x(i,j,g% nz+1:g% nz+Benkf% n_vert,d),       &
                        x(i,j,g% nz+1:g% nz+Benkf% m_vert,d), TR_ANA)
          !------------------
          ! wavelet transform
          !------------------
        end do
        end do
        end do
        !-------------------------------------
        ! replace field with wavelet transform
        !-------------------------------------
        deallocate (Benkf% X(k,l)% m(m)% ptr)
        MOVE_ALLOC (x, Benkf% X(k,l)% m(m)% ptr) ! Benkf% X(k,l)% m(m)% ptr => x
        Benkf% X(k,l)% m(m)% i% ub(3) = Benkf% j_vert
      end do
      Benkf% X(k,l)% ub(3) = Benkf% j_vert
      call set_pointers (Benkf% X(k,l))
    end do

  end subroutine v_lift_analysis

!==============================================================================

  subroutine v_lift_synth (w, x)
  type(t_atm) ,intent(inout) :: w ! IN
  type(t_atm) ,intent(inout) :: x ! OUT

    integer           :: i
    integer           :: j
    integer           :: d
    integer           :: m

    real(wp) ,allocatable :: y (:,:,:,:)

    allocate (y (x% lb(1) : x% ub(1), &
                 x% lb(2) : x% ub(2), &
                 Benkf% n_vert      , &
                 x% lb(4) : x% ub(4)) )

    !--------------------
    ! loop over variables
    !--------------------
    do m = 1, size (w% m)
      if (.not. w% m(m)% i% alloc) cycle
      if (.not. x% m(m)% i% alloc) call allocate (x, w% m(m)% i% name)

      !--------------------------------
      ! loop over horizontal gridpoints
      !--------------------------------
      do d = w% lb(4), w% ub(4)
      do j = w% lb(2), w% ub(2)
      do i = w% lb(1), w% ub(1)

        !------------------
        ! wavelet transform
        !------------------
        call lift_1d (y(i,j,:,d), w% m(m)% ptr (i,j,benkf% k_vert+1:,d), TR_SYN)

      end do
      end do
      end do

      !------------------------------
      ! iterpolation to physical grid
      !------------------------------
      call apply_vint (benkf% vloc, x% m(m)% ptr, y)

      !-----------------------------
      ! add grid scale perturbations
      !-----------------------------
      x% m(m)% ptr = x% m(m)% ptr + w% m(m)% ptr (:,:,:benkf% k_vert,:)

    end do

  end subroutine v_lift_synth
!==========================================================================
  subroutine v_lift_adj (w, x)
  type(t_atm) ,intent(inout) :: w
  type(t_atm) ,intent(in)    :: x

    integer           :: i
    integer           :: j
    integer           :: d
    integer           :: m

    real(wp) ,allocatable :: y (:,:,:,:)

!   !---------
!   ! printout
!   !---------
!   if(dace% lpio) then
!     write(6,*) '  v_lift_adj w% ub =',w% ub
!     write(6,*) '  v_lift_adj x% ub =',x% ub
!   endif

    allocate (y (x% lb(1) : x% ub(1), &
                 x% lb(2) : x% ub(2), &
                 Benkf% n_vert      , &
                 x% lb(4) : x% ub(4)) )

    !--------------------
    ! loop over variables
    !--------------------
    do m = 1, size (w% m)
      if (.not. x% m(m)% i% alloc) cycle
      if (.not. w% m(m)% i% alloc) call allocate (w, x% m(m)% i% name)

      !-----------------------------
      ! add grid scale perturbations
      !-----------------------------
      w% m(m)% ptr (:,:,:benkf% k_vert,:) = x% m(m)% ptr

      !------------------------------
      ! iterpolation to pressure grid
      !------------------------------
      call apply_vint_t (benkf% vloc, x% m(m)% ptr, y)

      !--------------------------------
      ! loop over horizontal gridpoints
      !--------------------------------
      do d = w% lb(4), w% ub(4)
      do j = w% lb(2), w% ub(2)
      do i = w% lb(1), w% ub(1)

        !------------------
        ! wavelet transform
        !------------------
        call lift_1d (y(i,j,:,d), w% m(m)% ptr (i,j,benkf% k_vert+1:,d), TR_ADJ)

      end do
      end do
      end do


    end do

  end subroutine v_lift_adj
!==========================================================================
  subroutine set_vscale
  !--------------------------------------------------------
  ! set vertical wavelet localisation normalisation factors
  !--------------------------------------------------------
    !----------
    ! variables
    !----------
    integer             :: i,j,k
    integer  ,parameter :: t = 100000  ! number of trials
    real(wp)            :: x  (benkf% m_vert)
    real(wp)            :: y  (benkf% m_vert)
    real(wp)            :: w  (benkf% m_vert)
    real(wp)            :: v  (benkf% m_vert)
    real(wp)            :: a
    integer             :: n

    n = benkf% n_vert

    !--------------
    ! preset arrays
    !--------------
    y = 0._wp
    a = 1._wp
    v = 1._wp
    j = 0
    k = n
    do
      v (j+1:j+k) = a
      a = a * c_scale
      j = j + k
      k = k / 2
      if (k==0) exit
    end do
    do i = 1, t
      call random_gauss (w)
      w = w * v
      j = 0
      k = n
      x = 0._wp
      do
        call lift_1d ( x(j+1:j+k),  w(j+1:), TR_SYN)
        j = j + k
        k = k / 2
        if (k==0) exit
      end do
      y = y + x**2
    end do
    y = sqrt (y / t)

    allocate (w_vscale (benkf% n_vert))
    allocate (a_vscale (benkf% n_vert))
    where (y==0) y = 1._wp
    w_vscale = 1._wp / y
    a_vscale = v

    if (dace% lpio) print *,'### set_vscale', y / v

  end subroutine set_vscale

!------------------------------------------------------------------------------

  subroutine apply_wloc (c1, a, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(out) :: a (:,:,:,:)
  real(wp)    ,intent(in)  :: v (:,:,:,:)

    integer             :: i1,i2,d
    integer             :: j,k

    real(wp) :: w (size(v,3))
    real(wp) :: x (size(v,3))

    !------------------
    ! scale interaction
    !------------------
    do d  = lbound(a,4), ubound(a,4)
    do i2 = lbound(a,2), ubound(a,2)
    do i1 = lbound(a,1), ubound(a,1)

      w = v (i1,i2,:,d)
      w = w * a_vscale
      j = 0
      k = benkf% n_vert
      x = 0._wp
      do
        call lift_1d ( x(j+1:j+k),  w(j+1:), TR_SYN)
        j = j + k
        k = k / 2
        if (k==0) exit
      end do
      x = x * w_vscale
      a (i1,i2,benkf% k_vert+1:,d) = x
    end do
    end do
    end do

    !------------------------
    ! grid scale disturbances
    !------------------------
    call apply_vloc (c1, a(:,:,:benkf% k_vert,:), &
                         a(:,:, benkf% k_vert+1:  &
                                benkf% k_vert+    &
                                benkf% n_vert,:)  )

  end subroutine apply_wloc

!------------------------------------------------------------------------------

  subroutine apply_wloc_t (c1, z, v)
  type (t_c1) ,intent(in)  :: c1
  real(wp)    ,intent(in)  :: z (c1%lb(1):,c1%lb(2):,:,:)
  real(wp)    ,intent(out) :: v (c1%lb(1):,c1%lb(2):,:,:)

    integer             :: i1,i2,d
    integer             :: j,k

    real(wp) :: w (size(v,3))
    real(wp) :: y (size(v,3))
    real(wp) :: x (size(v,3))

    !------------------------
    ! grid scale disturbances
    !------------------------
    v = 0._wp
    call apply_vloc_t (c1, z(:,:,:benkf% k_vert,:), &
                           v(:,:,:benkf% n_vert,:)  )
    v = v + z(:,:, benkf% k_vert+1:,:)
    !------------------
    ! scale interaction
    !------------------
    do d  = lbound(z,4), ubound(z,4)
    do i2 = lbound(z,2), ubound(z,2)
    do i1 = lbound(z,1), ubound(z,1)
      j = 0
      k = benkf% n_vert
      x = v (i1,i2,:,d)
      x = x * w_vscale
      y = 0._wp
      w = 0._wp
      do
        call lift_1d ( x(j+1:j+k),  y(j+1:), TR_ADJ)
        w = w + y
        j = j + k
        k = k / 2
        if (k==0) exit
      end do
      w = w * a_vscale
      v (i1,i2,:,d) = w
    end do
    end do
    end do

  end subroutine apply_wloc_t

!==========================================================================
end module mo_varenkf
