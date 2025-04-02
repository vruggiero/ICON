!
!+ IR emissivity model based on principal components
!
MODULE mo_ir_emis
!
! Description:
!   IR emissivity model based on principal components
!   from University of Wisconsin
!   adapted to DWD variational analysis system
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_28        2014/02/26 Andreas Rhodin
!  new module for usage of IR emissivity model based on principal components
! V1_29        2014/04/02 Andreas Rhodin
!  use PC fg coefficients for IR emissivity calculations
! V1_31        2014-08-21 Andreas Rhodin
!  specify fg errors, account for snow fraction, re-derive PC
! V1_48        2016-10-06 Andreas Rhodin
!  do not read variances derived from the atlas if not required
! V1_50        2017-01-09 Andreas Rhodin
!  base IR emissivity model on RTTOV12 sources; apply also for CRIS.
! V1_51        2017-02-24 Andreas Rhodin
!  change meaning of sink variable mode, make consistent with generalised humidity
!
! Code Description:
! Language: Fortran 2003
! Software Standards:
!
!==============================================================================
  !=============
  ! modules used
  !=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,        only: wp, sp         ! real precision kind parameters
  use mo_exception,   only: finish         ! abort in case of error
  use mo_mpi_dace,    only: dace,         &! MPI group info
                            p_gather,     &! generic MPI gather routine
                            p_scatterv,   &! generic MPI scatterv routine
                            p_sum,        &! MPI sum
                            p_bcast,      &! generic broadcast routine
                            p_bcast_ptr    ! broadcast pointers
  use mo_run_params,  only: data           ! constant data path
  use mo_time,        only: imm            ! get month from t-time
  use mo_namelist,    only: position_nml, &! routine to position nml group
                            nnml,         &! namelist fortran unit number
                            POSITIONED     ! position_nml: OK return flag
  !-------------------
  ! numerical routines
  !-------------------
  use mo_dec_matrix,  only: t_vector       ! decomposed vector
  use mo_algorithms,  only: index          ! sort, return index array
  use mo_matrix,      only: check_rs       ! eigenvectors of rs-matrix
  !-----------------------
  ! access to observations
  !-----------------------
  use mo_t_obs,       only: t_obs,        &! observation data type
                            t_spot,       &! report meta data type
                            set_sink       ! set sink variable parameters
  use mo_sink,        only: bc_xv          ! variable transformation
  use mo_rad,         only: t_rad_set,    &! radiance meta data type
                            rad_set,      &! radiance meta data
                            n_set,        &! number of valid rad_sets
                            set_indx,     &! get sat/instr/chan set index
                            chan_indx      ! get channel index in rad_set
  use mo_fdbk_tables, only: OT_RAD         ! radiance observation type
  use mo_rttov,       only: jplev,        &! no. of pressure levels
                            nsv            ! no. of single level inputs
  use mo_t_tovs,      only: t_tovs,       &! observation operator specific type
                            load,         &! load  t_tovs
                            store,        &! store t_tovs
                            destruct,     &! t_tovs  destructor routine
                            TTOVS_CEMI
  use mo_instrid,     only: rttov_instr    ! RTTOV number from instrument id
  !-----------------------
  ! access to gridded data
  !-----------------------
  use mo_wmo_tables,  only: WMO6_LATLON    ! latitude/longitude grid id.
  use mo_atm_state,   only: t_atm          ! atmospheric state derived type
  use mo_atm_grid,    only: t_grid,       &! grid meta data type
                            construct,    &! set up grid meta data
                            allocate,     &! allocate   grid type components
                            destruct       ! deallocate grid type components
  use mo_grid_intpol, only: grid_indices   ! get indices of neighbor gridpoints
  !----------------------------------------
  ! University Wisconsin IR emisivity atlas
  !----------------------------------------
  use mo_iratlas,     only: uwd,           &! IR emissivity atlas data
                            uwiremis_init, &! initialise the derived type
                            uwiremis_close_atlas, &! deallocate derived type
                            numpcs,        &! max. number of princ. components
                            numwave,       &! number of wavelenths
                            hsr_wavenum,   &! wave numbers
                            pcev,          &! principle component eigenvalues
                            pcm,           &! mean emissivity spectrum
                            sice_em,       &! sea ice emissivity spectrum
                            snow_em,       &! snow emissivity spectrum
                            bfemis_xgrid1, &! 1st longitude  in emiss. atlas
                            bfemis_ygrid1, &! 1st latitudes  in emiss. atlas
                            bfemis_gridres,&! resolution     of emiss. atlas
                            pca_stdv,      &! pc stdev (annual mean)
                            pcm_stdv        ! pc stdev (monthly mean)
#if (_RTTOV_VERSION >= 12) && !defined(__ICON__)
  use mod_uwiremis_atlas, only: rttov_uwiremis_angcorr, &
                                angcorrminzen   ! for the angular correction
#endif
   implicit none

  !================
  ! public entities
  !================
  private
  !------------
  ! subroutines
  !------------
  public :: read_nml_ir_emis ! read namelist /IREMIS/
  public :: init_ir_emis     ! initialise the module, read principal components
  public :: cleanup_ir_emis  ! deallocate module components
  public :: emis_pc          ! calculate emissivity from principal components
  public :: snw_frc_h        ! parameters to derive snow fraction from height
  public :: snw_frc_e        ! snow fraction (dummy parameter) error
  public :: n_pc             ! number of principle components to use
  !=================
  ! Module variables
  !=================
  integer  ,parameter   :: IASI  = 221  ! WMO instrument Id
  integer  ,parameter   :: CRIS  = 620  ! WMO instrument Id
  integer               :: npc   = 0    ! number of PCs used, 0 for not used
  real(wp) ,allocatable :: pcu (:,:)    ! principle component eigenvectors
  real(wp) ,allocatable :: pce   (:)    ! principle component eigenvalues
  real(wp) ,allocatable :: pces  (:)    ! location independent variance used
  real(wp) ,allocatable :: pcg   (:)    ! principle component global variances

  !-------------------
  ! namelist /IR_EMIS/               ! steering of IR PC emisivity model
  !-------------------
  integer  :: n_pc      = -1         ! number of principal components:
                                     !   -1: module not used, use FASTEM
                                     !    0: use mean emissivity profile
                                     !    1: mean emissivity + PC disturbances
  integer  :: fg_src    = 1          ! first guess source:
                                     !    0: use zero first guess
                                     !    1: use emissivity atlas
                                     !    2: cycled fg (not implemented)
  real(wp) :: fg_err(5) = 0._wp      ! first guess error contributions:
                                     !  (1): constant term
                                     !  (2): global   variance
                                     !  (3): spatial  variance
                                     !  (4): seasonal variance
                                     !  (5): PC       variance
                                     !  (6): cycled   (not implemented)
  real(wp) :: snw_frc_h(0:1) = 0._wp ! derive bg snow fraction from height
  real(wp) :: snw_frc_e      = 0._wp ! snow fraction error
  real(wp) :: wavenum_bnd(2) = (/0._wp,9999._wp/)
                                     ! bounds for PC decomposition
  logical  :: do_ang_corr = .false.  ! do angular correction
  integer  :: version     = 0        ! <12:old  12:from rttov 12 sources
  character(len=256) :: path = ''    ! path for atlas data

  namelist /IR_EMIS/ n_pc, fg_src, fg_err, snw_frc_h, snw_frc_e, wavenum_bnd, &
                     version, path, do_ang_corr

  !--------------------------------------------
  ! quantities derived from namelist parameters
  !--------------------------------------------
  logical  :: new_pc      = .false.  ! re-derive PC decomposition
  logical  :: read_pc     = .false.  ! derived flag:
                                     !    T: input of PC + emis.atlas required
                                     !    F: no further input required
  logical  :: stdv_pc     = .false.  ! input of emis.atlas std.deviations

  !-----------------------------------------------
  ! derived type for PC coeff. MPI scatter request
  !-----------------------------------------------
  type t_idx
   integer :: ib  ! destination box  index
   integer :: is  ! destination spot index
   integer :: pe  ! destination MPI rank
   integer :: i1  ! source longitude grid index
   integer :: i2  ! source latitude  grid index
  end type t_idx

!==============================================================================
contains
!==============================================================================

  subroutine init_ir_emis (bg, o, bi, ei, rad_set)
  type(t_atm)    ,intent(in)    :: bg        ! background atmospheric state
  type(t_obs)    ,intent(inout) :: o(:)      ! observation data
  type(t_vector) ,intent(inout) :: bi        ! interpolated background
  type(t_vector) ,intent(inout) :: ei        ! interpolated background error
  type(t_rad_set),intent(inout) :: rad_set(:)! radiance meta data
  !---------------------------------
  ! initialise the module
  ! read principal component vectors
  ! and first guess coefficients
  !---------------------------------

    logical               :: verbose = .false. ! verbosity flag
    integer               :: ierr              ! error return flag
    integer               :: month             ! month
    type(t_grid) ,pointer :: g                 ! emissivity atlas grid
    integer               :: npce              ! # of PC read

    !------------------------------------
    ! re-set flag to read PCs
    ! check if IR PC model is switched on
    !------------------------------------
    if (.not. read_pc) return
    read_pc = .false.
    !-------------------------------------------------------
    ! read the data into the mo_iratlas module on the I/O PE
    !-------------------------------------------------------
    npc  = min (n_pc, numwave)
    npce = npc
    if (dace% lpio) then
      if (new_pc) npce = numwave       ! read all PCs for reconstruction
      allocate (pcu  (npce, numwave))  ! PC eigenvectors
      allocate (pce  (npce         ))  ! PC eigenvalues (variances)
      month = imm (bg% time)
      call uwiremis_init(          &
                  trim(path)//'/', &! <- file path
                  month,           &! <- month for climatology
                  stdv_pc,         &! <- read stdev of PC atlas
                  verbose,         &! <- verbosity flag
                  do_ang_corr,     &! <- do angular correction
                  version,         &! <- atlas version
                  ierr,            &! -> error return variable
              npc=npce             )

      if (ierr/=0) call finish('init_ir_emis',         &
                               'error in uwiremis_init')
      pcu  = uwd% pcu (1:npce,:)    ! copy eigenvectors to this module
      pce  = pcev     (1:npce)      ! copy eigenvalues  to this module
      !----------------------------
      ! re-derive PCs within bounds
      !----------------------------
      if (new_pc) call pc_bnd
    else
      allocate (pcu  (npce, numwave))  ! PC eigenvectors
      allocate (pce  (npce         ))  ! PC eigenvalues (variances)
    endif
    allocate   (pces (npc          ))  ! PC variances used (static)
    !---------------------------
    ! broadcast scalar constants
    !---------------------------
    call p_bcast (uwd% single_inst ,dace% pio)
    call p_bcast (uwd% std_init    ,dace% pio)
    call p_bcast (uwd% do_ang_corr ,dace% pio)
    call p_bcast (uwd% nb_lats     ,dace% pio)
    call p_bcast (uwd% nb_lons     ,dace% pio)
    call p_bcast (uwd% nb_pack     ,dace% pio)
    call p_bcast (uwd% cv_lats     ,dace% pio)
    call p_bcast (uwd% cv_lons     ,dace% pio)
    call p_bcast (uwd% cv_pack     ,dace% pio)
    call p_bcast (uwd% igbp_lats   ,dace% pio)
    call p_bcast (uwd% igbp_lons   ,dace% pio)
    call p_bcast (uwd% nb_igbp     ,dace% pio)
    !--------------------------------------------------
    ! broadcast array components for angular correction
    !--------------------------------------------------
    if (uwd% do_ang_corr) then
      call p_bcast_ptr (uwd% igbp ,dace% pio)
      call p_bcast_ptr (uwd% p1d  ,dace% pio)
      call p_bcast_ptr (uwd% p2d  ,dace% pio)
      call p_bcast_ptr (uwd% p3d  ,dace% pio)
      call p_bcast_ptr (uwd% p1n  ,dace% pio)
      call p_bcast_ptr (uwd% p2n  ,dace% pio)
      call p_bcast_ptr (uwd% p3n  ,dace% pio)
    endif
    !---------------------------------------------------------------
    ! broadcast principle component eigenvectors + values to all PEs
    !---------------------------------------------------------------
    call p_bcast (pcu(:npc,:), dace% pio)
    call p_bcast (pce(:npc  ), dace% pio)
    !---------------------------------------------------
    ! calculate location independent part of PC variance
    !---------------------------------------------------
    pces = fg_err(1)**2              &! constant term
         + fg_err(5)**2 * pce(:n_pc)  ! PC eigenvalues
    if (dace% lpio) then
      write (6,'()')
      write (6,'(a)') repeat ('-',79)
      write (6,'()')
      write (6,'(a)')      '  initialise mo_ir_emis'
      write (6,'()')
      write (6,'(a,10f8.2,(/15x,10f8.2))') '    pces      =',sqrt (pces)
      write (6,'()')
    endif
    !-----------------------------------------
    ! set up grid meta data for the atlas data
    !-----------------------------------------
    allocate (g)
    call construct (g, gridtype= WMO6_LATLON,               &
                             nx= uwd% nb_lons,              &
                             ny= uwd% nb_lats,              &
                            lo1= bfemis_xgrid1  / 1000._wp, &
                            la1= bfemis_ygrid1  / 1000._wp, &
                             di= bfemis_gridres / 1000._wp, &
                             dj=-bfemis_gridres / 1000._wp  )

    !-----------------------------------------------
    ! distribute emissivity flag (misuse 'lsm' flag)
    !-----------------------------------------------
    call allocate (g,'lsm')
    if (dace% lpio) g% lsm (:,:,1,1) = transpose (uwd% bfemis_flag)
    call p_bcast (g% lsm (:,:,1,1), dace% pio)

    !---------------------------------
    ! set up sensor specific meta data
    !---------------------------------
    call pc_chan (rad_set)

    !------------------------------------------
    ! distribute the atlas data to observations
    !------------------------------------------
    call pc2obs (g, o, bi, ei)

    !------------------------------
    ! set angular correction factor
    !------------------------------
    if (uwd% do_ang_corr) call  set_cemi (o, rad_set)

    !---------
    ! clean up
    !---------
    call destruct (g)
    deallocate    (g)
    if (uwd% do_ang_corr) then
      deallocate (uwd% igbp)
      deallocate (uwd% p1d )
      deallocate (uwd% p2d )
      deallocate (uwd% p3d )
      deallocate (uwd% p1n )
      deallocate (uwd% p2n )
      deallocate (uwd% p3n )
    endif

    !-------------------------------
    ! clean up the mo_iratlas module
    !-------------------------------
    if (dace% lpio) call uwiremis_close_atlas

  end subroutine init_ir_emis

!--------------------------------------------------------------------------

  subroutine pc2obs (g, o, bi, ei)
  !------------------------------------------------------
  ! scatter the PC coefficients to observation data types
  !------------------------------------------------------
  type(t_grid)   ,intent(in)    :: g    ! emissivity atlas grid
  type(t_obs)    ,intent(inout) :: o(:) ! observation data
  type(t_vector) ,intent(inout) :: bi   ! interpolated background
  type(t_vector) ,intent(inout) :: ei   ! interpolated background error

    real(wp)                 :: w   (4)           ! weight
    integer                  :: idx (4,4)         ! [Point, index]
    integer                  :: jk  (4)
    integer                  :: np                ! no of neigbours
    integer                  :: ib
    integer                  :: is
    integer                  :: n, n0             ! number of IASI FOVs
    integer                  :: n1 (1:dace% npe)
    integer                  :: i, j, k
    integer                  :: i1, i2, i3
    real(wp)                 :: dlon, dlat
    type(t_idx) ,allocatable :: tix  (:)
    type(t_idx) ,allocatable :: tix0 (:)
    logical                  :: land
    real(wp)    ,allocatable :: coef0 (:)
    real(wp)    ,allocatable :: coef  (:)
    real(wp)    ,allocatable :: var0  (:)
    real(wp)    ,allocatable :: var   (:)
    integer                  :: npv               ! no of profile inputs
    integer                  :: ncv               ! total no of inputs
    integer                  :: npcl              ! no of PC inputs (this spot)
    integer                  :: npcf              ! no of PC bg coefficients
    real(wp)                 :: snw_frc           ! snow fraction

    if (any(rad_set(1:n_set)%gopts%lev_mode > 0)) call finish('pc2obs', &
         'Radiance datasets with lev_mode > 0 not implemented so far')

    npv  = 2*jplev

    !--------------------------------------
    ! set index array size for observations
    !--------------------------------------
    n = 0
    do ib = 1, size(o)
      if (o(ib)% pe /= dace% pe) cycle
      do is = 1, o(ib)% n_spot
        if (o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        do i = o(ib)% spot(is)% o% i + 1, &
               o(ib)% spot(is)% o% i +    &
               o(ib)% spot(is)% o% n
          select case (o(ib)% body(i)% lev_sig)
          case (IASI, CRIS)
            n = n + 1
            exit
          end select
        end do
      end do
    end do

    !-------------------------------------------
    ! set up column index array for observations
    !-------------------------------------------
    allocate (tix (n))
    n = 0
    do ib = 1, size(o)
      if (o(ib)% pe /= dace% pe) cycle
      do is = 1, o(ib)% n_spot
        if (o(ib)% spot(is)% hd% obstype /= OT_RAD) cycle
        do i = o(ib)% spot(is)% o% i + 1, &
               o(ib)% spot(is)% o% i +    &
               o(ib)% spot(is)% o% n
          select case (o(ib)% body(i)% lev_sig)
          case (IASI, CRIS)
            dlat = o(ib)% spot(is)% col% c% dlat
            dlon = o(ib)% spot(is)% col% c% dlon
            call grid_indices (dlon,     & ! <-- geodetic longitude
                               dlat,     & ! <-- geodetic latitude
                               g,        & ! <-- grid data type
                               idx,      & ! --> Grid point indices
                               w,        & ! --> Weights
                               np        ) ! --> number of points returned
            jk = index (w)
            land = .false.
            do k = np,1,-1
              i1 = idx(jk(k),1)
              i2 = idx(jk(k),2)
              if (g% lsm (i1,i2,1,1) > 0) then
                land = .true.
                n = n + 1
                tix(n)% ib = ib
                tix(n)% is = is
                tix(n)% pe = dace% pe
                tix(n)% i1 = i1
                tix(n)% i2 = i2
                exit
              endif
            end do
            if (.not.land) then
              !-------------------------------
              ! no coefficients for this point
              !-------------------------------
              ncv  = o(ib)% spot(is)% i% n
              npcl = ncv - nsv - npv
              i2   = o(ib)% spot(is)% i% i + ncv
              i1   = i2 - npcl + 1
              if (npcl /= npc) call finish('pc2obs','npcl /= npc')
              bi% s(ib)% x(i1:i2) = 0._wp
              ei% s(ib)% x(i1:i2) = 0._wp
            endif
            exit
          end select
        end do
      end do
    end do

    !-------------
    ! send request
    !-------------
    call p_gather (n, n1, dace% pio)
    n0 = p_sum (n)
    if (dace% lpio) then
      allocate (tix0(n0))
    else
      allocate (tix0(0))
    endif
    call p_gather_idx (tix(1:n), tix0, dace% pio)

    !-----------------------
    ! prepare data to return
    !-----------------------
    npcf = min (npc, numpcs)
    allocate   (pcg   (npcf)     )
    if (dace% lpio) then
      allocate (coef0 (npcf * n0))
      allocate (var0  (npcf * n0))
      pcg = sum (pca_stdv (:,:npcf), dim=1) / size (pca_stdv, 1)
      do i = 1, n0
        !------------------------------------------------
        ! Find the emissivity PC coefs and their variance
        !------------------------------------------------
        j = uwd% bfemis_lut(tix0(i)%i2, (tix0(i)%i1))
        k = (i-1) * npcf
        if (j > 0) then
          coef0 (k+1:k+npcf) =  uwd% pca_coef (j,:npcf) * uwd% pca_sfac(:npcf) &
                                                        + uwd% pca_offs(:npcf)
          var0  (k+1:k+npcf) =      (pcg      (  :npcf) * fg_err(2)) ** 2 &
                             +      (pca_stdv (j,:npcf) * fg_err(3)) ** 2 &
                             +      (pcm_stdv (j,:npcf) * fg_err(4)) ** 2
        else
          coef0 (k+1:k+npcf) = 0._wp
          var0  (k+1:k+npcf) = 0._wp
        end if
      end do
    else
      allocate (coef0 (0))
      allocate (var0  (0))
    endif

    !------------
    ! return data
    !------------
    n1 = n1 * npcf
    allocate (coef (n*npcf))
    allocate (var  (n*npcf))
    call p_scatterv (coef0, n1, coef, dace% pio)
    call p_scatterv (var0,  n1, var,  dace% pio)
    call p_bcast    (pcg,             dace% pio)
    !--------------------------
    ! fill into data structures
    !--------------------------
    do j = 1, n
      k    = (j-1) * npcf
      ib   = tix(j)% ib
      is   = tix(j)% is
      ncv  = o(ib)% spot(is)% i% n
      npcl = ncv - nsv - npv
      i2   = o(ib)% spot(is)% i% i + ncv
      i1   = i2 - npcl + 1
      i3   = i1 + npcf - 1
      if (npcl == 0) cycle
      if (npcl /= npc) call finish ('pc2obs','npcl /= npc')
      !------------------------------------------------
      ! set background principal component coefficients
      !------------------------------------------------
      bi% s(ib)% x(i1:i2) = 0._wp
      select case (fg_src)
      case (1)
        bi% s(ib)% x(i1:i3) = coef (k+1:k+npcf)
      end select
      !-----------------------------
      ! set background snow fraction
      !-----------------------------
      snw_frc = 0._wp
      if (o(ib)% spot(is)% hs_bg <= snw_frc_h(0)) then
        snw_frc = 0._wp
      else if (o(ib)% spot(is)% hs_bg >= snw_frc_h(1)) then
        snw_frc = 1._wp
      else
        snw_frc = o(ib)% spot(is)% hs_bg - snw_frc_h(0) &
                        / ( snw_frc_h(1) - snw_frc_h(0) )
      endif
      call set_sink (o(ib)% spot(is), o(ib), 3, snw_frc, snw_frc_e, npv+5, &
                     0._wp, 1._wp, 0.5_wp, 0.5_wp, 2                       )
      !----------------------
      ! set coefficient error
      !----------------------
      ei% s(ib)% x(i1:i2) = pces (:npcl) ! spacially invariant term
      ei% s(ib)% x(i1:i3) =       ei% s(ib)% x(i1:i3) + var (k+1:k+npcf)
      ei% s(ib)% x(i1:i2) = sqrt (ei% s(ib)% x(i1:i2))
    end do

  end subroutine pc2obs

!--------------------------------------------------------------------------

  subroutine pc_bnd
  !----------------------------------------
  ! derive new PC eigenvector + eigenvalues
  ! for wavenumbers within bounds
  !----------------------------------------

    integer               :: i, j, k
    integer               :: i1, i2
    real(wp) ,allocatable :: cov (:,:)
    real(wp) ,allocatable :: evc (:,:)
    real(wp) ,allocatable :: evl (:)
    integer               :: npce
    integer               :: npcn
    integer               :: npc1, npc2, npc3

    !--------------------
    ! derive array bounds
    !--------------------
    i1 = 1; i2 = numwave
    do i = 1, numwave
      if (hsr_wavenum(i) < wavenum_bnd(1)) i1 = i
      if (hsr_wavenum(i) > wavenum_bnd(2)) then
        i2 = i
        exit
      endif
    end do
    npcn = i2 - i1 + 1

    !-------------------------
    ! re-construct covariances
    !-------------------------
    npce = size (pce)
    allocate (cov (i1:i2, i1:i2))
    allocate (evc (i1:i2,  npcn))
    allocate (evl         (npcn))
    do i = i1, i2
      do j = i1, i2
        cov (i,j) = 0._wp
        do k = 1, npce
          cov (i,j) = cov (i,j) + pcu (k,i) * pce (k) * pcu (k,j)
        end do
      end do
    end do

    !-------------------------
    ! eigenvalue decomposition
    !-------------------------
    call check_rs (cov, evc=evc, evl=evl)

    !----------------------
    ! copy back leading PCs
    !----------------------
    npc1 = 1
    npc2 = n_pc
    if (fg_src == 1) npc1 = min (n_pc, numpcs) + 1
    npc3 = min (npc2 - npc1 + 1, npcn)

    deallocate (pcu, pce)
    allocate   (pcu (npc2,  numwave))
    allocate   (pce (npc2))
    pcu (:npc1-1,:) = uwd% pcu (1:npc1-1,:)
    pcu ( npc1: ,:) = 0._wp
    pce             = 0._wp
    do k = 1, npc3
      i = k + npc1 - 1
      j = npcn + 1 - k
      pcu (i, i1:i2) = evc (:,j)
      pce (i)        = evl   (j)
    end do

  end subroutine pc_bnd

!==========================================================================

  subroutine cleanup_ir_emis
  !-----------------------------
  ! deallocate module components
  !-----------------------------

    if (allocated  (pcu ))  deallocate (pcu )
    if (allocated  (pce ))  deallocate (pce )
    if (allocated  (pces))  deallocate (pces)

  end subroutine cleanup_ir_emis

!==========================================================================

  subroutine pc_chan (rad_set)
  !-------------------------------------------------------
  ! interpolate principle component eigenvectors
  ! and mean emissivities to wave numbers of used channels
  !-------------------------------------------------------
  type (t_rad_set) ,intent(inout) :: rad_set(:)  ! radiance meta data
  target :: rad_set

    !---------------------------------------------------
    ! Smallest wavenumber        of IASI
    ! Largest wavenumber         of IASI
    ! Wavenumber difference between IASI channels (1/cm)
    !---------------------------------------------------
    real(wp) ,parameter :: iasi_start_wvn =  645._wp
    real(wp) ,parameter :: iasi_end_wvn   = 2760._wp
    real(wp) ,parameter :: iasi_dwvn      =  0.25_wp

    real(wp) ,parameter :: cris_start_wvn_long  = 650._wp
    real(wp) ,parameter :: cris_end_wvn_long    = 1095._wp
    real(wp) ,parameter :: cris_dwvn_long       = 0.625_wp

    real(wp) ,parameter :: cris_start_wvn_mid   = 1210._wp
    real(wp) ,parameter :: cris_end_wvn_mid     = 1750._wp
    real(wp) ,parameter :: cris_dwvn_mid        = 1.25_wp

    real(wp) ,parameter :: cris_start_wvn_short = 2155._wp
    real(wp) ,parameter :: cris_end_wvn_short   = 2550._wp
    real(wp) ,parameter :: cris_dwvn_short      = 2.5_wp

    integer                  :: i, ii, ic, i1, i2, j
    integer                  :: n_chan
    integer                  :: n_instr
    type(t_rad_set) ,pointer :: rs
    real(wp)                 :: w1, w2

    !----------------------------------
    ! cycle over satellites/instruments
    !----------------------------------
    do i = 1, size (rad_set)
      rs => rad_set(i)

      !-------------------------------------------
      ! check if PC emissivity model is applicable
      !-------------------------------------------
      n_chan  = rs% n_chan
      n_instr = rs% n_instr
      if (n_chan  == 0)                  cycle
      if (n_instr == 0)                  cycle
      if (npc     == 0)                  cycle
      if (all ( rs% instr_wmo /= IASI &
          .and. rs% instr_wmo /= CRIS )) cycle
      !--------------------
      ! allocate components
      !--------------------
      if (.not.associated (rs% emis_pc)) then
        allocate (rs% wavenum        (n_chan)); rs% wavenum   =  0._wp
        allocate (rs% emis_land      (n_chan)); rs% emis_land =  0._wp
        allocate (rs% emis_snow      (n_chan)); rs% emis_snow =  0._wp
        allocate (rs% emis_sice      (n_chan)); rs% emis_sice =  0._wp
        allocate (rs% i1             (n_chan)); rs% i1        = -1
        allocate (rs% w1             (n_chan)); rs% w1        =  0._wp
        allocate (rs% emis_pc   (npc, n_chan)); rs% emis_pc   =  0._wp
      endif

      !----------------------
      ! loop over instruments
      !----------------------
      do ii = 1, n_instr

        !-------------------
        ! loop over channels
        !-------------------
        select case (rs% instr_wmo (ii))
        case (IASI, CRIS)
          do ic = rs% o_ch_i(ii) + 1, rs% o_ch_i(ii) + rs% n_ch_i(ii)

            !--------------------------------------
            ! derive wavenumber from channel number
            !--------------------------------------
            select case (rs% instr_wmo (ii))
            case (IASI)
              !-----
              ! IASI
              !-----
              rs% wavenum (ic) = (rs% chan(ic)-1) * iasi_dwvn + iasi_start_wvn
            case (CRIS)
              !-----
              ! CRIS
              !-----
              select case (rs%chan(ic))
              case ( 27:713)    ! longwave band
                rs% wavenum (ic) = (rs% chan(ic)-1) * cris_dwvn_long  + cris_start_wvn_long
              case ( 716:1142)  ! midwave band
                rs% wavenum (ic) = (rs% chan(ic)-1) * cris_dwvn_mid   + cris_start_wvn_mid
              case (1147:1301)  ! shortwave band
                rs% wavenum (ic) = (rs% chan(ic)-1) * cris_dwvn_short + cris_start_wvn_short
              end select
            end select

            !-------------------------------------
            ! interpolate PC and mean emissivities
            !-------------------------------------
            i1 = floor ( (rs% wavenum (ic) - hsr_wavenum(1)) / 5._wp) + 1
            i1 = max (0, min (numwave, i1))
            select case (i1)
            case default
              w2 = (rs% wavenum (ic) - hsr_wavenum(i1)) / 5._wp
              i2 = i1 + 1
            case (0)
              i1 = 1
              i2 = i1 + 1
              w2 = 0._wp
            case (numwave)
              i1 = numwave - 1
              i2 = numwave
              w2 = 1._wp
            end select
            w1 = 1._wp - w2
            rs% emis_land   (ic) = w1 * pcm     (i1) + w2 * pcm     (i2)
            rs% emis_snow   (ic) = w1 * snow_em (i1) + w2 * snow_em (i2)
            rs% emis_sice   (ic) = w1 * sice_em (i1) + w2 * sice_em (i2)
            rs% i1          (ic) = i1
            rs% w1          (ic) = w1
            do j = 1, npc
              rs% emis_pc (j,ic) = w1 * pcu   (j,i1) + w2 * pcu   (j,i2)
            end do
          end do
        end select
      end do
    end do
  end subroutine pc_chan

!----------------------------------------------------------------------------

  subroutine set_cemi (o, rad_set)
  !--------------------------------
  ! store angular correction factor
  !--------------------------------
  type(t_obs)      ,intent(inout) :: o(:)        ! observation meta data
  type (t_rad_set) ,intent(in)    :: rad_set(:)  ! radiance meta data

    integer               :: ib, is, i, j, k, i1
    integer               :: satid, grid, iset, instr, ichan, rtins, igbp_type
    real(wp)              :: satzen, solzen, dlat, dlon, angcorr(2), w1
    integer               :: gridx, gridy
    integer               :: ilat, ilon
    type(t_spot) ,pointer :: si
    type(t_tovs)          :: tovs
    logical               :: first

    integer, parameter :: igbp_gridres =      50     !    0.05 deg
    integer, parameter :: igbp_ygrid1  =   89950     !   89.95 deg
    integer, parameter :: igbp_xgrid1  = -179950     ! -179.95 deg

    do ib = 1, size(o)
      if (o(ib)% pe /= dace% pe) cycle
      do is = 1, o(ib)% n_spot
        si => o(ib)% spot(is)
        if (si% hd% obstype /= OT_RAD) cycle
        satid  = si% hd% satid                  ! satellite id
        grid   = si% hd% grid_id                ! grid id
        iset   = set_indx (rad_set (1:n_set), satid=satid, grid=grid)
        first  = .true.

        satzen = si% stzen
        solzen = si% sozen
        dlat   = si% col% c% dlat
        dlon   = si% col% c% dlon
        if (dlon < -180._wp) dlon = dlon + 360.0_wp
        if (dlon >= 180._wp) dlon = dlon - 360.0_wp
        ilat = NINT (dlat * 1000._wp)
        ilon = NINT (dlon * 1000._wp)

        if (iset <= 0) call finish('set_cemi','iset <= 0')
        do k = 1, si% o% n
          i = k + si% o% i
          instr = int (o(ib)% body(i)% lev_sig)
          select case (instr)
          case (IASI, CRIS)
            rtins = rttov_instr (instr, satid)
            ichan = nint(o(ib)% olev (i))
            j     = chan_indx (rtins, ichan, rad_set(iset))
            i1    = rad_set(iset)% i1(j)
            w1    = rad_set(iset)% w1(j)

            gridy = nint (abs (igbp_ygrid1 - ilat) * 1._wp / igbp_gridres) + 1
            gridx = nint (abs (igbp_xgrid1 - ilon) * 1._wp / igbp_gridres) + 1
            igbp_type = uwd% igbp (gridy, gridx)

#if (_RTTOV_VERSION >= 12) && !defined(__ICON__)
            if (satzen > angcorrminzen) then
              call rttov_uwiremis_angcorr (     &
                  uwd% p1d(i1:i1+1,:), uwd% p2d(i1:i1+1,:), uwd% p3d(i1:i1+1,:), &
                  uwd% p1n(i1:i1+1,:), uwd% p2n(i1:i1+1,:), uwd% p3n(i1:i1+1,:), &
                  solzen,                       &! in
                  satzen,                       &! in
                  igbp_type,                    &! in
                  angcorr)                       ! out
              if (first) then
                first = .false.
                call load (o(ib), si, tovs, tovs_io=0)
!                if (associated(tovs%cemi)) deallocate(tovs%cemi)
                allocate(tovs%cemi(tovs%nchan))
                tovs%init  = tovs%init + TTOVS_CEMI
                tovs%ncemi = tovs%nchan
              endif
              tovs% cemi (k) = w1 * angcorr(1) + (1._wp-w1) * angcorr(2)
!print *,'###',k,i1,w1,angcorr,tovs% cemi(k)
            endif
#endif
          end select
        end do
        if (.not.first) then
          call store (o(ib), si, tovs, tovs_io=TTOVS_CEMI)
          call destruct         (tovs)
        endif
      end do
    end do

  end subroutine set_cemi

!============================================================================

  subroutine emis_pc (spot, cemi, obs, snowfrac, stype, coef, emiss, de_dc)
  !--------------------------------------------------------------
  ! calculate emisivities from principal component coefficients
  !   observational data provided in derived types t_spot / t_obs
  !--------------------------------------------------------------
  type(t_spot)       ,intent(in)  :: spot         ! FOV meta data
  real(sp)           ,intent(in)  :: cemi(:)      ! angular correction factor
  type(t_obs)        ,intent(in)  :: obs          ! observation meta data
  real(wp)           ,intent(in)  :: snowfrac     ! snow fraction
  integer            ,intent(in)  :: stype        ! surface type
  real(wp)           ,intent(in)  :: coef  (:)    ! PC coefficients
  real(wp)           ,intent(out) :: emiss   (:)  ! emissivities
  real(wp) ,optional ,intent(out) :: de_dc (0:,:) ! Jakobi matrix d emis / d pc

    integer            :: satid, grid, iset, ichan, instr, rtins
    integer            :: i               ! observation index in report
    integer            :: j               ! observation index in 'rad_set'
    integer            :: k               ! observation index in 'body'
    real(wp)           :: em              ! preliminary emissivities (snow free)
    integer ,parameter :: bc    = 1       ! use log transform
    real(wp),parameter :: xmaxr = 0.01_wp ! nonlinear range at upper bound
    real(wp),parameter :: xminr = 0.50_wp ! nonlinear range at lower bound
    real(wp)           :: d               ! Jakobian of nonlinear transform
    real(wp)           :: c               ! angular correction factor
    !-------------
    ! set defaults
    !-------------
    emiss = 0._wp
    if (present (de_dc)) de_dc = 0._wp

    !----------------------------------
    ! return if PC approach is not used
    !----------------------------------
    if (npc <= 0) return

    !-----------------------------
    ! get sensor meta data indices
    !-----------------------------
    satid  = spot% hd% satid                  ! satellite id
    grid   = spot% hd% grid_id                ! grid id
    iset   = set_indx (rad_set (1:n_set), satid=satid, grid=grid)
    if (iset <= 0) call finish('emis_pc','iset <= 0')

    !------------------
    ! copy back results
    !------------------
    do i = 1, spot% o% n
      k     = spot% o% i + i
      instr = int (obs% body (k)% lev_sig)
      select case (instr)
      case (IASI, CRIS)
        rtins = rttov_instr (instr, satid)
        ichan =              nint(obs% olev (k))
        j     = chan_indx (rtins, ichan, rad_set(iset))

        select case (stype)
        !--------------------------------
        ! sea ice: use presribed spectrum
        !--------------------------------
        case (2)
          emiss(i) = rad_set(iset)% emis_sice (j)

        !-----------------------------
        ! sea: keep 0 for FASTEM model
        !-----------------------------
        case default

        !---------------------------
        ! land surface: use PC model
        !---------------------------
        case (0)
          !-----------------------------------------------------
          ! preliminary emissivities from pc (linear regression)
          !-----------------------------------------------------
          c  = cemi (i)
          d  = 1._wp - snowfrac
          em = c * sum (coef * rad_set(iset)% emis_pc (:,j)) &
                             + rad_set(iset)% emis_land (j)
          if (present(de_dc)) then
            !-----------------------------------------------
            ! derivative with respect to the PC coefficients
            ! index 0:   with respect to snow fraction
            !-----------------------------------------------
            de_dc (1:,i) = rad_set(iset)% emis_pc (:,j) * d * c
            de_dc (0 ,i) = rad_set(iset)% emis_snow (j) &
                         - em
          endif
          em = d * em + snowfrac * rad_set(iset)% emis_snow(j)
          !------------------------------------------------------------
          ! final emissivities (nonlinear transform for emissivity < 1)
          !------------------------------------------------------------
          call bc_xv (emiss(i), em, 0._wp, xminr, xmaxr, 1._wp, bc, d)
          if (present(de_dc)) de_dc (:,i) = de_dc (:,i) * d
        end select
      end select
    end do

  end subroutine emis_pc

!==========================================================================

  subroutine read_nml_ir_emis
  !------------------------
  ! read namelist /IR_EMIS/
  !------------------------

    integer :: ierr
    !------------
    ! set default
    !------------
    n_pc        = 6       ! 6 PC by default
    version     = 0
    do_ang_corr = .false. ! do angular correction
    path        = data
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write (6,'()')
      write (6,'(a)') repeat ('-',79)
      write (6,'()')
      call position_nml ('IR_EMIS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=IR_EMIS, iostat=ierr)
        if (ierr/=0) call finish ('read_nml_ir_emis',          &
                                  'ERROR in namelist /IR_EMIS/')
#else
        read (nnml ,nml=IR_EMIS)
#endif
        write (6,'(a)') '  namelist /IR_EMIS/ read'
      case default
        write (6,'(a)') '  namelist /IR_EMIS/ not present'
        n_pc = -1    ! switch off IR emissivity model
      end select
      write (6,'()')
      !-------------------
      ! consistency checks
      !-------------------
      fg_err = max (0._wp, fg_err)
      n_pc   = min (n_pc,  numwave)
      if (n_pc <= 0)                  fg_src  = 0
      if (n_pc <= 0)                  fg_err  = 0._wp
      if (n_pc >  0)                  read_pc = .true.
      if (any (fg_err(2:4) /= 0._wp)) stdv_pc = .true.
      new_pc = wavenum_bnd(1) > hsr_wavenum(1) .or. &
               wavenum_bnd(2) < hsr_wavenum(numwave)
      !---------
      ! printout
      !---------
      write(6,'(a,i5)')    '    n_pc        = ',n_pc
      write(6,'(a,i5)')    '    fg_src      = ',fg_src
      write(6,'(a,6f8.2)') '    fg_err      = ',fg_err
      write(6,'(a,2f8.2)') '    snw_frc_h   = ',snw_frc_h
      write(6,'(a, f8.2)') '    snw_frc_e   = ',snw_frc_e
      write(6,'(a,2f8.2)') '    wavenum_bnd = ',wavenum_bnd
      write(6,'(a,4x,l1)') '    read_pc     = ',read_pc
      write(6,'(a,4x,l1)') '    stdv_pc     = ',stdv_pc
      write(6,'(a,4x,l1)') '    new_pc      = ',new_pc
      write(6,'(a,4x,l1)') '    do_ang_corr = ',do_ang_corr
      write(6,'(a,i5)')    '    version     = ',version
      write(6,'(a,a)')     '    path        = ',trim(path)
      write(6,'()')
      !-----------------
      ! invalid settings
      !-----------------
      select case (fg_src)
      case (0,1)
      case default
        call finish ('read_nml_ir_emis','fg_src/=1,2  not implemented')
      end select
!     if (fg_err(6) /= 0._wp ) call finish ('read_nml_ir_emis',        &
!                                           'fg_err(6) not implemented')
    endif
    !----------
    ! broadcast
    !----------
    call p_bcast (n_pc        ,dace% pio)
    call p_bcast (fg_src      ,dace% pio)
    call p_bcast (fg_err      ,dace% pio)
    call p_bcast (snw_frc_h   ,dace% pio)
    call p_bcast (snw_frc_e   ,dace% pio)
    call p_bcast (wavenum_bnd ,dace% pio)
    call p_bcast (read_pc     ,dace% pio)
    call p_bcast (stdv_pc     ,dace% pio)
    call p_bcast (new_pc      ,dace% pio)
    call p_bcast (version     ,dace% pio)
    call p_bcast (do_ang_corr ,dace% pio)
    call p_bcast (path        ,dace% pio)
  end  subroutine read_nml_ir_emis

!==========================================================================
#define DERIVED type(t_idx)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_idx
#include "p_gather_derived.incf"
#undef  DERIVED
!==========================================================================
end module mo_ir_emis
