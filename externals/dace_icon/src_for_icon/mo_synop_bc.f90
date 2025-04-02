!
!+ Bias correction for SYNOP data
!
MODULE mo_synop_bc
!
! Description:
!   Bias correction for SYNOP data.
!   Also used for aggregation of precipitation and wind gust data.
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress, Elisabeth Bauernschubert
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_48        2016-10-06 Andreas Rhodin
!  SYNOP bias correction, used for aggregation of precip.
! V1_49        2016-10-25 Harald Anlauf
!  t_spot: generalize bc_airep -> bc_index
! V1_50        2017-01-09 Harald Anlauf
!  Fix 'bug' in commented code
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin   DWD  2016  original source
! Harald Anlauf    DWD  2016
! E.Bauernschubert DWD  2020
!-----------------------------------------------------------

!=============
! Modules used
!=============

  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_kind,       only: wp,              &! working precision kind parameter
                           sp,              &! single precision kind parameter
                           i8                ! 8 byte integer kind parameter
  use mo_physics,    only: pi                ! 3.1415...
  use mo_exception,  only: finish            ! abort in case of error
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_sum,           &! sum over PEs
                           p_max,           &! max over PEs
                           p_bcast,         &! generic MPI broadcast routine
                           p_gather,        &! generic MPI gather routine
                           p_barrier         ! MPI barrier
  use mo_run_params, only: ana_time,        &! analysis time
                           flag_biasc_synop  ! SYNOP bias correction
  use mo_time,       only: days,            &! derive days (real)
                           hours,           &! derive hours
                           ihh,             &! derive daily hours
                           operator(-),     &! calculate time difference
                           cyyyymmddhhmm,   &! string from time
                           time_cyyyymmddhhmm! time from string
  use slatec_module, only: sort              ! sort routine
  use mo_fdbk_tables,only: OT_SYNOP,        &! SYNOP report type ID
                           OT_DRIBU,        &! DRIBU report type ID
                           VN_Z,            &! varno for geopotential height
                           VN_PS,           &! varno for surface pressure
                           VN_T2M,          &! varno for 2 meter temperature
                           VN_RH2M,         &! varno for 2 metre relative humidity
                           VN_N,            &! total cloud amount
                           VN_RR,           &! varno for rain rate
                           VN_RAD_GL,       &! global  radiation
                           VN_RAD_DF,       &! diffuse radiation
                           VN_GUST,         &! varno for gust
                           VN_TMIN, VN_TMAX  ! varno for tmin, tmax
  use mo_dec_matrix, only: t_vector          ! vector data type
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_t_datum,    only: rvind             ! invalid value
  use mo_t_use,      only: decr_use,        &! decrease the state of a datum
                           STAT_DISMISS,    &! observation status flags
                           STAT_PASSIVE,    &! ...
                           STAT_ACTIVE_0I,  &! ...
                           STAT_ACCEPTED,   &! observation status flags
                           CHK_BIASCOR,     &! flag for no bias correction
                           CHK_NO_OBS        ! no observation in report
  use mo_obs_set,    only: t_obs_set         ! observation data derived type
  use mo_obs_tables, only: rept_use          ! use table entry
  use mo_biasc_io,   only: t_bcor_head,     &! file header derived type
                           new_bc_head,     &! construct new file header
                           open_bc_read,    &! open  file for reading
                           open_bc_write,   &! open  file for writing
                           close_bc,        &! close file
                           bc_paths,        &! full pathnames
                           bc_obstyp,       &! observation type in files
                           nbcf,            &! number of files in list
                           verbose           ! Verbosity level of bias corr.
  use mo_wigos,      only: t_wsi,           &! WIGOS id data type
                           wsi_mode,        &! WIGOS station id mode
                           wsi_prefix,      &! Prefix of "fake station id"
                           wsi_decode,      &! Derive WSI from given string
                           wsi_encode,      &! Derive statid,stathash from t_wsi
                           wsi_to_text       ! Convert t_wsi to text format
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,        only:                  &! NetCDF f90 interface
                           nf90_def_dim,    &! define dimensions
                           nf90_def_var,    &! define variables
                           nf90_put_att,    &! define attribute
                           nf90_enddef,     &! end of definition mode
                           nf90_put_var,    &! write variable
                           nf90_strerror,   &! derive error character string
                           NF90_FLOAT,      &! float     type id
                           NF90_INT,        &! integer   type id
                           NF90_SHORT,      &! short int type id
                           NF90_CHAR,       &! character type id
                           NF90_GLOBAL,     &! global attribute flag
                           NF90_NOERR,      &! NetCDF return code for no error
                           NF90_FILL_FLOAT   ! NetCDF fillvalue
  use mo_t_netcdf,   only: ncid,            &! NetCDF file id to use
                           stanc,           &! NetCDF start  parameter to use
                           counc,           &! NetCDF count  parameter to use
                           strnc,           &! NetCDF stride parameter to use
                           get_attr,        &! get attribute
                           get_var,         &! get variable
                           get_dim           ! get dimension
  implicit none

!================
! Public entities
!================

  private
  !--------------
  ! derived types
  !--------------
  public :: t_synop_bc       ! synop bias correction data
! public :: t_bc             ! component of t_synop_bc
! public :: t_plat_bc        ! component of t_synop_bc
  !------------
  ! subroutines
  !------------
  public :: synop_bc_init      ! initialize module: read namelist, biascor.coeff.
  public :: synop_bc_bfg       ! apply bias correction, called before fg-check
! public :: synop_bc_afg       ! update bias correction coefs., after first guess
  public :: synop_bc_aan       ! update bias correction coefs., after analysis
  public :: read_synop_bc_file ! read synop bias correction file
  !----------------------------------------------------
  ! namelist parameters to be set in namelist SYNOP_OBS
  !----------------------------------------------------
  public :: biascor_mode     ! mode used for updating
  public :: reg_param        ! regularization parameter for non-linear bc
  public :: t_decay          ! accumulation decay time (days)
  public :: n_required       ! number of entries required for correction
  public :: bc_fallback      ! fallback if biasc-file not present
  public :: aggregate        ! aggregate rain rate observations
  public :: bc_t2m           ! apply b/c for t2m (nonlinear bc)?
  public :: bc_rh2m          ! apply b/c for rh2m (nonlinear bc)?
  public :: bc_synop         ! apply b/c for SYNOP Land?
  public :: bc_ship          ! apply b/c for SYNOP Ship?
  public :: bc_buoy          ! apply b/c for DRIBU?
  public :: bc_metar         ! apply b/c for METAR?
  public :: chk_zstation     ! check consistency of reported station height
  public :: z0_bc_synop      ! threshold parameter for SYNOP Land
  public :: z0_bc_ship       ! threshold parameter for SYNOP Ship
  public :: z0_bc_buoy       ! threshold parameter for DRIBU
  public :: z0_bc_metar      ! threshold parameter for METAR
  public :: compress_ag      ! internally "compress" aggregated data?
  public :: BC_NOBC,BC_UP,BC_FG,BC_AN,BC_VARBC,BC_INVAL ! values for biascor_mode

!=========================
! Derived type definitions
!=========================

  !----------------------
  ! station specific data
  !----------------------
  type t_bc
    integer  :: n         = 0      ! number of entries
    real(wp) :: ob        = 0._wp  ! obs-bg  mean deviation
    real(wp) :: oa        = 0._wp  ! obs-ana mean deviation
    real(wp) :: ob_ob     = 0._wp  ! obs-bg  variance
    real(wp) :: oa_oa     = 0._wp  ! obs-ana variance
    real(wp) :: n_ac      = 0._wp  ! number of accumulated entries
    real(wp) :: ob_ac     = 0._wp  ! obs-bg  mean     accumulated
    real(wp) :: oa_ac     = 0._wp  ! ana-bg  mean     accumulated
    real(wp) :: ob_ob_ac  = 0._wp  ! obs-bg  variance accumulated
    real(wp) :: oa_oa_ac  = 0._wp  ! ana-bg  variance accumulated
    real(wp) :: o_err     = 0._wp  ! nominal obs error
    real(wp) :: b_err     = 0._wp  ! nominal bg  error
    real(wp) :: o_err_ac  = 0._wp  ! nominal obs error
    real(wp) :: b_err_ac  = 0._wp  ! nominal bg  error
    real(wp) :: bc_b      = 0._wp  ! bias correction background value
    real(wp) :: bc_a      = 0._wp  ! bias correction analysis   value
!   real(wp) :: bc_b_err  = 0._wp  ! bias correction background error
!   real(wp) :: bc_a_err  = 0._wp  ! bias correction analysis   error
  end type t_bc

  !----------------------------
  ! aggregation of observations
  !----------------------------
  type t_ag
    integer  :: ibc         =  -1     ! element index in array
    integer  :: hours       =   0     ! time of last entry
    real(sp) :: rr_1 (0:23) = -99._sp !  1h rain rates
    real(sp) :: rr_3 (0: 7) = -99._sp !  3h rain rates
    real(sp) :: rr_6 (0: 3) = -99._sp !  6h rain rates
    real(sp) :: rr12 (0: 1) = -99._sp ! 12h rain rates
    real(sp) :: gu_1 (0:23) = -99._sp !  1h gusts
    real(sp) :: gu_3 (0: 7) = -99._sp !  3h gusts
    real(sp) :: gu_6 (0: 3) = -99._sp !  6h gusts
    real(sp) :: ma_1 (0:23) = -99._sp !  1h tmax
    real(sp) :: ma_6 (0: 3) = -99._sp !  6h tmax
    real(sp) :: mi_1 (0:23) = -99._sp !  1h tmin
    real(sp) :: mi_6 (0: 3) = -99._sp !  6h tmin
    real(sp) :: gr_1 (0:23) = -99._sp !  1h global  radiation
    real(sp) :: gr_3 (0: 7) = -99._sp !  3h global  radiation
    real(sp) :: gr_6 (0: 3) = -99._sp !  6h global  radiation
    real(sp) :: dr_1 (0:23) = -99._sp !  1h diffuse radiation
    real(sp) :: dr_3 (0: 7) = -99._sp !  3h diffuse radiation
    real(sp) :: dr_6 (0: 3) = -99._sp !  6h diffuse radiation
  end type t_ag

  !----------------------------------------------------------
  ! number of coefficients (for T2M and RH2M bias correction)
  !----------------------------------------------------------
  integer, parameter :: ncoeff    = 5 * 2  ! total number of coefficients; 5 trigonometric basis functions
                                           ! for diurnal cycle, 2 polynomials for cloudcover

  !----------------------------------------------------------
  ! Coefficients (for T2M and RH2M bias correction)
  !----------------------------------------------------------
  type t_coeff
    real(wp) :: t2m_bc_b (ncoeff) = 0._wp  ! coefficients for nonlinear bias correction for t2m  (background)
    real(wp) :: t2m_bc_a (ncoeff) = 0._wp  ! coefficients for nonlinear bias correction for t2m  (analysis)
    real(wp) :: rh2m_bc_b(ncoeff) = 0._wp  ! coefficients for nonlinear bias correction for rh2m (background)
    real(wp) :: rh2m_bc_a(ncoeff) = 0._wp  ! coefficients for nonlinear bias correction for rh2m (analysis)
  end type t_coeff

  !------------------------------------------
  ! observed quantity indices and varno array
  !------------------------------------------
  integer, parameter :: nob     = 3     ! number of observed quantities
  integer, parameter :: OB_PS   = 1     ! surface pressure
  integer, parameter :: OB_T2M  = 2     ! 2 meter temperature
  integer, parameter :: OB_RH2M = 3     ! 2 metre relative humidity (varno 58)
  integer, parameter :: vvarno(3) = (/VN_PS, VN_T2M, VN_RH2M/)

  !-----------------------
  ! platform specific data
  !-----------------------
  type t_plat_bc
    character(len=8) :: statid  = ''       ! SYNOP station id
    integer  (i8)    :: hash    = 0        ! integer representation of name
    integer          :: code    = -1       ! CMA codetype
    real(sp)         :: z       = -999._sp ! station height [m]
    real(sp)         :: lat     = -999._sp ! (last known) latitude
    real(sp)         :: lon     = -999._sp ! (last known) longitude
    type(t_wsi)      :: wsi                ! WIGOS station id
  end type t_plat_bc

  !---------------------
  ! bias correction data
  !---------------------
  type t_synop_bc
    type(t_bcor_head)        :: h                 ! file header data
    integer                  :: biascor_mode = 0  ! mode used for updating
    integer                  :: n            = 0  ! number of entries
    type (t_plat_bc),pointer :: plat (:)          ! platform specific metadata
    type (t_ag)     ,pointer :: ag   (:)          ! observations to aggregate
    type (t_coeff)  ,pointer :: coeff(:)          ! coefficients for nonlinear bias correction
    type (t_bc)     ,pointer :: data (:,:)        ! data (nob,:)
    logical                  :: compressed   = .false.
  end type t_synop_bc

  !-------------------------------------
  ! derived type used for cross-checking
  !-------------------------------------
  type t_tmp
    character(len=8) :: statid ! station name
    integer  (i8)    :: hash   ! integer representation of name
    integer          :: code   ! CMA codetype
    integer          :: pe     ! processor of report
    integer          :: i      ! old index
    integer          :: j      ! new index
    real(sp)         :: lat    ! latitude
    real(sp)         :: lon    ! longitude
    type(t_wsi)      :: wsi    ! WIGOS id
  end type t_tmp

  !-----------------------------------
  ! values for biascor_mode definition
  !-----------------------------------
  integer ,parameter :: BC_INVAL = -9  ! not set so far
  integer ,parameter :: BC_NOBC  =  0  ! no bias correction
  integer ,parameter :: BC_UP    =  1  ! only update bias correction file
  integer ,parameter :: BC_FG    =  2  ! apply  bias corr. from first guess
  integer ,parameter :: BC_AN    =  3  ! apply  bias corr. from analysis
  integer ,parameter :: BC_VARBC =  4  ! variational bias correction


!=================
! Module variables
!=================

  type(t_synop_bc),save :: synop_bc                 ! bias correction data
  type(t_synop_bc),save :: empty_bc                 ! empty bias correction
  !-----------------
  ! namelist entries
  !-----------------
  integer               :: biascor_mode = BC_INVAL  ! mode used for updating
  real(wp)              :: reg_param    = 50        ! regularization parameter for non-linear bc
  real(wp)              :: t_decay      =   -30._wp ! accumulation decaytime
  integer               :: n_required   =    50     ! # of entries required
  logical               :: bc_fallback  = .false.   ! biasc-file not present
  logical               :: aggregate    = .true.    ! aggregate observations
  logical               :: bc_t2m       = .false.   ! apply b/c for t2m (nonlinear bc)?
  logical               :: bc_rh2m      = .false.   ! apply b/c for rh2m (nonlinear bc)?
  logical               :: bc_synop     = .false.   ! apply b/c for SYNOP Land?
  logical               :: bc_ship      = .false.   ! apply b/c for SYNOP Ship?
  logical               :: bc_buoy      = .false.   ! apply b/c for DRIBU?
  logical               :: bc_metar     = .false.   ! apply b/c for METAR?
  logical               :: chk_zstation = .false.   ! check station height
  real(wp)              :: z0_bc_synop  =   999._wp ! threshold for SYNOP Land
  real(wp)              :: z0_bc_ship   =   999._wp ! threshold for SYNOP Ship
  real(wp)              :: z0_bc_buoy   =   999._wp ! threshold for DRIBU
  real(wp)              :: z0_bc_metar  =   999._wp ! threshold for METAR
  logical               :: compress_ag  = .true.    ! "compress" ag data?

contains
!==============================================================================

  subroutine synop_bc_init
  !---------------------------------------
  ! Initialize this module:
  ! read bias-correction coefficient files
  !---------------------------------------

    integer :: i    ! loop index
    integer :: ierr ! error return value


    !--------------------------------------------------------------------
    ! set defaults for namelist parameters depending on 'ga3_biasc_synop'
    ! if not set so far
    !--------------------------------------------------------------------
    if (biascor_mode == BC_INVAL) then
      select case (flag_biasc_synop)
      case (-1)
        biascor_mode = BC_NOBC
        bc_fallback  = .false.
      case ( 0)
        biascor_mode = BC_FG
        bc_fallback  = .true.
      case ( 1)
        biascor_mode = BC_FG
        bc_fallback  = .false.
      end select
    endif

    if (flag_biasc_synop < 0) then
      bc_synop = .false.
      bc_ship  = .false.
      bc_buoy  = .false.
      bc_metar = .false.
    end if

    if (biascor_mode /= 0) then
      !---------------------------------------
      ! read bias-correction coefficient files
      !---------------------------------------
      ierr = -1
      if (flag_biasc_synop /= 0 .or. .not. bc_fallback) then
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_SYNOP) then
            call read_synop_bc_file (synop_bc, bc_paths(i))
            ierr = 0
            exit
          endif
        end do
      end if
      !-----------------------------------------------
      ! create empty bias-correction coefficient files
      !-----------------------------------------------
      if (ierr /= 0) then
        if (.not.bc_fallback) &
          call finish ('synop_bc_init','coefficient file not present !')
        call new_bc_head (synop_bc% h, OT_SYNOP, t_decay=(/abs(t_decay)/))
        synop_bc% biascor_mode = biascor_mode
        synop_bc% n            = 0
        allocate (synop_bc% plat (0))
        allocate (synop_bc% ag   (0))
        allocate (synop_bc% data (nob,0))
        if (bc_t2m .or. bc_rh2m) allocate (synop_bc% coeff (0))
        if (dace% lpio) write(6,'(a,a)')'    created empty file ', &
                                             trim (synop_bc%h% path)
      endif
      !-------------------------------
      ! prepare bias correction to use
      !-------------------------------
      select case (biascor_mode)
      case (BC_UP)
        synop_bc% data% bc_b = 0
        synop_bc% data% bc_a = 0
      case (BC_FG)
        where (synop_bc% data% n_ac >= n_required)
          synop_bc% data% bc_b = synop_bc% data% ob_ac
          synop_bc% data% bc_a = synop_bc% data% ob_ac
        elsewhere
          synop_bc% data% bc_b = NF90_FILL_FLOAT
          synop_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_AN)
        where (synop_bc% data% n_ac >= n_required)
          synop_bc% data% bc_b = synop_bc% data% oa_ac
          synop_bc% data% bc_a = synop_bc% data% oa_ac
        elsewhere
          synop_bc% data% bc_b = NF90_FILL_FLOAT
          synop_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_VARBC)
        synop_bc% data% bc_b = synop_bc% data% bc_a
      end select
      ! In case of t2m and rh2m bias correction we use the coefficients in synop_bc% coeff.
      if (bc_t2m) then
         do i = 1, synop_bc% n
            synop_bc% coeff(i)% t2m_bc_b  = synop_bc% coeff(i)% t2m_bc_a
         end do
      end if
      if (bc_rh2m) then
         do i = 1, synop_bc% n
            synop_bc% coeff(i)% rh2m_bc_b = synop_bc% coeff(i)% rh2m_bc_a
         end do
      end if
    endif

  end subroutine synop_bc_init

!------------------------------------------------------------------------------

  subroutine read_synop_bc_file (bc, file, inst)
  !--------------------------------
  ! read synop bias correction file
  !--------------------------------
  type(t_synop_bc)    ,intent(out) :: bc    ! derived type variable to fill
  character(len=*)    ,intent(in)  :: file  ! name of file to read
  logical   ,optional ,intent(in)  :: inst  ! read instantaneous data as well

    integer                    :: i, ierr
    logical                    :: linst
    real(sp)      ,allocatable :: rr(:,:)
    type(t_bc)        ,pointer :: p (:,:)   ! pointer to statistics file entries
    integer                    :: wsilen    ! WSI length (dimension)
    integer                    :: wsimode   ! WSI mode found in b/c file
    integer(i8)                :: wsihash
    character(10)              :: statid    ! "Fake" station ID
    character(32), allocatable :: wsichr(:) ! Station ID (expanded WSI)
    type(t_wsi)                :: wsi

    linst = .false.; if (present(inst)) linst = inst

    if (dace% lpio) then
      !-------------------------------
      ! open NetCDF file, read header
      !-------------------------------
      bc% h% path = file
      call open_bc_read (bc% h)

      !----------------
      ! read dimensions
      !----------------
      ncid  = bc% h% ncid
      call get_dim (bc% n ,'platform')

      !--------------------
      ! allocate components
      !--------------------
      allocate (bc% plat (    bc% n))
      allocate (bc% ag   (    bc% n))
      allocate (bc% data (nob,bc% n))
      if (bc_t2m .or. bc_rh2m) allocate (bc% coeff (bc% n))

      !---------------
      ! read variables
      !---------------
      if (bc% n > 0) then
        stanc = 1
        counc = 0
        strnc = 1

        !--------
        ! 1h data
        !--------
        allocate (rr (24, bc% n))

        call get_var (rr ,'rr_1'  )
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% rr_1 = rr(:,i); enddo

        call get_var (rr ,'gust_1')
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gu_1 = rr(:,i); enddo

        call get_var (rr ,'gsr_1', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gr_1 = rr(:,i); enddo

        call get_var (rr ,'dsr_1', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% dr_1 = rr(:,i); enddo

        call get_var (rr ,'tmax_1', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% ma_1 = rr(:,i); enddo

        call get_var (rr ,'tmin_1', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% mi_1 = rr(:,i); enddo

        deallocate (rr)
        !--------
        ! 3h data
        !--------
        allocate (rr ( 8, bc% n))

        call get_var (rr ,'rr_3'  )
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% rr_3 = rr(:,i); enddo

        call get_var (rr ,'gust_3')
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gu_3 = rr(:,i); enddo

        call get_var (rr ,'gsr_3', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gr_3 = rr(:,i); enddo

        call get_var (rr ,'dsr_3', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% dr_3 = rr(:,i); enddo

        deallocate (rr)
        !--------
        ! 6h data
        !--------
        allocate (rr ( 4, bc% n))

        call get_var (rr ,'rr_6'  )
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% rr_6 = rr(:,i); enddo

        call get_var (rr ,'gust_6')
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gu_6 = rr(:,i); enddo

        call get_var (rr ,'gsr_6', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% gr_6 = rr(:,i); enddo

        call get_var (rr ,'dsr_6', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% dr_6 = rr(:,i); enddo

        call get_var (rr ,'tmax_6', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% ma_6 = rr(:,i); enddo

        call get_var (rr ,'tmin_6', ierr=ierr)
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% mi_6 = rr(:,i); enddo

        deallocate (rr)
        !---------
        ! 12h data
        !---------
        allocate (rr ( 2, bc% n))

        call get_var (rr ,'rr_12' )
        where (rr == rvind) rr = -99._sp
        do i = 1, bc% n; bc% ag(i)% rr12 = rr(:,i); enddo

        deallocate (rr)

        call get_var (bc% ag  % hours      ,'hours'   )
        call get_var (bc% plat% statid     ,'statid'  )
!       bc% plat% code = -1
!       bc% plat% z    = -999._sp
!       bc% plat% lat  = -999._sp
!       bc% plat% lon  = -999._sp
        if (bc% h% version /= "00.00") then
          call get_var (bc% plat% code     ,'codetype')
          if (chk_zstation) then
           call get_var(bc% plat% z        ,'z_station')
           where (bc% plat% z == NF90_FILL_FLOAT)
             bc% plat% z   = -999._sp
           end where
          end if
          call get_var (bc% plat% lat      ,'lat'      )
          call get_var (bc% plat% lon      ,'lon'      )
          where (bc% plat% lat == NF90_FILL_FLOAT)
             bc% plat% lat = -999._sp
             bc% plat% lon = -999._sp
          end where
          !-----------------------------------------------------------------
          ! SHIP: do not keep previous position, but remember the actual one
          !-----------------------------------------------------------------
!         where (bc% plat% code == 21)
!            bc% plat% lat = -999._sp
!            bc% plat% lon = -999._sp
!         end where
          select case (bc% h% version)
          case ("00.01")
             !----------------------------------------------------------------
             ! Version "00.01": Bias correction is/was only implemented for PS
             !----------------------------------------------------------------
             p => bc% data(:1,:)
             call get_var (p% n_ac     ,'n_ac'    )
             call get_var (p% ob_ac    ,'ob_ac'   )
             call get_var (p% oa_ac    ,'oa_ac'   )
             call get_var (p% ob_ob_ac ,'ob_ob_ac')
             call get_var (p% oa_oa_ac ,'oa_oa_ac')
             call get_var (p% o_err_ac ,'o_err_ac')
             call get_var (p% b_err_ac ,'b_err_ac')
             call get_var (p% bc_b     ,'bc_b'    )
             call get_var (p% bc_a     ,'bc_a'    )
             if (linst) then
                call get_var (p% n     ,'n'    )
                call get_var (p% ob    ,'ob'   )
                call get_var (p% oa    ,'oa'   )
                call get_var (p% ob_ob ,'ob_ob')
                call get_var (p% oa_oa ,'oa_oa')
                call get_var (p% o_err ,'o_err')
                call get_var (p% b_err ,'b_err')
             endif
          case ("01.00")
             !---------------------------------------------------------------------------
             ! Version "01.00": Additionally, a nonlinear bias correction is implemented
             ! for t2m and rh2m. The statistics for the classical approach are written to
             ! the derived type bc% data% ...
             !---------------------------------------------------------------------------
             call get_var (bc% data% n_ac     ,'n_ac'    )
             call get_var (bc% data% ob_ac    ,'ob_ac'   )
             call get_var (bc% data% oa_ac    ,'oa_ac'   )
             call get_var (bc% data% ob_ob_ac ,'ob_ob_ac')
             call get_var (bc% data% oa_oa_ac ,'oa_oa_ac')
             call get_var (bc% data% o_err_ac ,'o_err_ac')
             call get_var (bc% data% b_err_ac ,'b_err_ac')
             call get_var (bc% data% bc_b     ,'bc_b'    )
             call get_var (bc% data% bc_a     ,'bc_a'    )
             if (linst) then
                call get_var (bc% data% n     ,'n'    )
                call get_var (bc% data% ob    ,'ob'   )
                call get_var (bc% data% oa    ,'oa'   )
                call get_var (bc% data% ob_ob ,'ob_ob')
                call get_var (bc% data% oa_oa ,'oa_oa')
                call get_var (bc% data% o_err ,'o_err')
                call get_var (bc% data% b_err ,'b_err')
             endif
             !------------------------------------------------------------------------
             ! If version "01.00" of bias_SYNOP file and nonlinear bias correction for
             ! t2m or rh2m is switched on, then get the correct coefficients from file
             !------------------------------------------------------------------------
!          if (bc% h% version == "01.00") then
             if (bc_t2m) then
                allocate (rr ( ncoeff, bc% n))
                call get_var (rr ,'bc_t2m_b' )
                do i = 1, bc% n; bc% coeff(i)% t2m_bc_b = rr(:,i); enddo
                call get_var (rr ,'bc_t2m_a' )
                do i = 1, bc% n; bc% coeff(i)% t2m_bc_a = rr(:,i); enddo
                deallocate (rr)
             end if
             if (bc_rh2m) then
                allocate (rr ( ncoeff, bc% n))
                call get_var (rr ,'bc_rh2m_b' )
                do i = 1, bc% n; bc% coeff(i)% rh2m_bc_b = rr(:,i); enddo
                call get_var (rr ,'bc_rh2m_a' )
                do i = 1, bc% n; bc% coeff(i)% rh2m_bc_a = rr(:,i); enddo
                deallocate (rr)
             end if
             end select
!          end if
        end if
        do i = 1, bc% n
          bc% plat(i)% hash = transfer (bc% plat(i)% statid, bc% plat(i)% hash)
          bc% ag  (i)% ibc  = i
        end do

        if (wsi_mode > 0) then
          call get_dim (wsilen, 'wsichars', ierr)
          if (ierr /= NF90_NOERR) then
            write(*,*) "Note: wsi is missing in ", trim (file)
          else
            if (wsilen /= 32) call finish ("read_synop_bc","wsichars/=32")
            call get_attr (wsimode, 'wsi', 'wsi_mode')
            if (wsimode /= wsi_mode)                                   &
                 write(*,'(a,i0,a,i0,a)') " Note: wsi_mode differs: ", &
                 wsimode, " (found) /= ", wsi_mode, " (current)"
            if (wsimode <= 0 .or. wsimode > 3)                      &
                 call finish ("read_synop_bc","invalid wsi:wsi_mode")

            allocate (wsichr(bc% n))
            call get_var (wsichr, 'wsi')
            do i = 1, bc% n
               if (bc% plat(i)% statid(1:1) /= wsi_prefix) cycle
               if (wsichr(i)(1:2) == "0-") then
                  call wsi_decode (wsichr(i), wsi, ierr)
                  if (ierr == 0 .and. wsi% wigii < 20000) then
                     call wsi_encode (wsi, statid, wsihash)
                     !-------------------------------------------------------
                     ! Retain stored/used fake station ID for ECMWF-like mode
                     ! unless wsi_mode is changed from non-ECMWF-like:
                     !-------------------------------------------------------
                     if (wsi_mode /= 1 .or. wsimode /= 1) &
                          bc% plat(i)% statid = statid
                     bc%      plat(i)% hash   = wsihash
                     bc%      plat(i)% wsi    = wsi
                  else
                     write(0,*) "read_synop_bc: bad WSI: ", trim (wsichr(i))
                  end if
               end if
            end do
          end if
        end if
      endif
      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif

    !----------
    ! broadcast
    !----------
    call p_bcast (bc% n, dace% pio)
    if (.not.dace% lpio) then
      allocate (bc% plat (    bc% n))
      allocate (bc% ag   (    bc% n))
      allocate (bc% data (nob,bc% n))
      if (bc_t2m .or. bc_rh2m) then
        allocate (bc% coeff(  bc% n))
      endif
    endif
    call p_bcast_ag   (bc% ag  , dace% pio)
    call p_bcast_plat (bc% plat, dace% pio)
    call p_bcast_data (bc% data, dace% pio)
    if (bc_t2m .or. bc_rh2m) then
      call p_bcast_coeff(bc% coeff,dace% pio)
    endif

  end subroutine read_synop_bc_file

!------------------------------------------------------------------------------

  subroutine destruct_bc_file (bc)
  type(t_synop_bc) ,intent(inout) :: bc    ! derived type variable to fill

    !----------------------
    ! deallocate components
    !----------------------
    deallocate (bc% plat)
    deallocate (bc% ag  )
    if (bc_t2m .or. bc_rh2m) deallocate (bc% coeff)
    deallocate (bc% data)

    bc = empty_bc
  end subroutine destruct_bc_file

!------------------------------------------------------------------------------

  subroutine write_synop_bc_file (bc)
  !---------------------------------
  ! write synop bias correction file
  !---------------------------------
  type(t_synop_bc) ,intent(inout) :: bc

    integer :: status                            ! NetCDF return value
    integer :: dimid (2)                         ! NetCDF dimension id
    integer :: dimcoef                           ! NetCDF dimension id for bc coefficients
    integer :: dimch, dimwch                     ! NetCDF dimension id
    integer :: dimh24, dimh8, dimh4, dimh2       ! NetCDF dimension id
    integer :: statid, hours, wsi                ! NetCDF variable ids
    integer :: codetype, z_station, lat, lon     ! .
    integer :: varno                             ! NetCDF variable ids
    integer :: rr_1, rr_3, rr_6, rr_12           ! NetCDF variable ids
    integer :: gu_1, gu_3, gu_6                  ! NetCDF variable ids
    integer :: gr_1, gr_3, gr_6                  ! NetCDF variable ids
    integer :: dr_1, dr_3, dr_6                  ! NetCDF variable ids
    integer :: ma_1, mi_1, ma_6, mi_6            ! NetCDF variable ids
    integer :: n, n_ac                           ! .
    integer :: ob, oa                            ! .
    integer :: ob_ob, oa_oa                      ! .
    integer :: ob_ac, oa_ac                      ! .
    integer :: ob_ob_ac, oa_oa_ac                ! .
    integer :: bc_b, b_err , b_err_ac ! bc_b_err ! .
    integer :: bc_a, o_err , o_err_ac ! bc_a_err ! .
    integer :: bc_t2m_b, bc_rh2m_b, bc_t2m_a, bc_rh2m_a
    !++++++++++++++++++++++++++++++++++++++++++
    ! workaround for NEC SX9 nf90_put_var(text)
    !++++++++++++++++++++++++++++++++++++++++++
    character(len=8),allocatable :: statids(:)
    character(32)   ,allocatable :: wsichar(:) ! Station ID (expanded WSI)
    real(wp)        ,allocatable :: rr   (:,:)
    integer                      :: i

    if (dace% lpio .and. bc%n > 0) then
      !-------------------------------
      ! open NetCDF file, write header
      !-------------------------------
!     if (bc_t2m .or. bc_rh2m) then
        bc% h% version     = "01.00"         ! File format version
!     else
!       bc% h% version     = "00.01"         ! File format version
!     endif
      bc% h% t_decay(2:) = bc% h% t_decay(1) ! Write only the first t_decay value
      call open_bc_write (bc% h, OT_SYNOP)

      !---------------
      ! set attributes
      !---------------
      status = nf90_put_att (bc%h% ncid, NF90_GLOBAL, 'biascor_mode', &
                                                       biascor_mode   )
      !------------------
      ! define dimensions
      !------------------
      call def_dim ('observation', nob,            dimid(1))
      call def_dim ('platform',    size(bc% plat), dimid(2))
      call def_dim ('coeff',       ncoeff,         dimcoef )
      call def_dim ('chars',       8,              dimch   )
      if (wsi_mode > 0)                                    &
      call def_dim ('wsichars',    len (wsichar),  dimwch  )
      call def_dim ('h_24',        24,             dimh24  )
      call def_dim ('h_8',         8,              dimh8   )
      call def_dim ('h_4',         4,              dimh4   )
      call def_dim ('h_2',         2,              dimh2   )

      !-----------------
      ! define variables
      !-----------------
      status = nf90_def_var (bc%h% ncid ,'statid' ,NF90_CHAR, &
                             (/dimch,dimid(2)/), statid)
      status = nf90_put_att (bc%h% ncid , statid, 'longname', &
                             'station id as character string' )
      if (wsi_mode > 0) then
       status = nf90_def_var(bc%h% ncid ,'wsi'    ,NF90_CHAR, &
                             [dimwch,dimid(2)],   wsi         )
       status = nf90_put_att(bc%h% ncid , wsi,    'longname', &
                             'WIGOS station id as string'     )
       status = nf90_put_att(bc%h% ncid , wsi,    'wsi_mode', &
                             wsi_mode                         )
      end if

      status = nf90_def_var (bc%h% ncid ,'rr_12', NF90_FLOAT, &
                             (/dimh2, dimid(2)/), rr_12       )
      status = nf90_def_var (bc%h% ncid ,'rr_6' , NF90_FLOAT, &
                             (/dimh4, dimid(2)/), rr_6        )
      status = nf90_def_var (bc%h% ncid ,'rr_3' , NF90_FLOAT, &
                             (/dimh8, dimid(2)/), rr_3        )
      status = nf90_def_var (bc%h% ncid ,'rr_1' , NF90_FLOAT, &
                             (/dimh24,dimid(2)/), rr_1        )

      status = nf90_def_var (bc%h% ncid ,'gust_6', NF90_FLOAT, &
                             (/dimh4, dimid(2)/), gu_6        )
      status = nf90_def_var (bc%h% ncid ,'gust_3', NF90_FLOAT, &
                             (/dimh8, dimid(2)/), gu_3        )
      status = nf90_def_var (bc%h% ncid ,'gust_1', NF90_FLOAT, &
                             (/dimh24,dimid(2)/), gu_1        )

      status = nf90_def_var (bc%h% ncid ,'gsr_6', NF90_FLOAT, &
                             (/dimh4, dimid(2)/), gr_6        )
      status = nf90_def_var (bc%h% ncid ,'gsr_3', NF90_FLOAT, &
                             (/dimh8, dimid(2)/), gr_3        )
      status = nf90_def_var (bc%h% ncid ,'gsr_1', NF90_FLOAT, &
                             (/dimh24,dimid(2)/), gr_1        )

      status = nf90_def_var (bc%h% ncid ,'dsr_6', NF90_FLOAT, &
                             (/dimh4, dimid(2)/), dr_6        )
      status = nf90_def_var (bc%h% ncid ,'dsr_3', NF90_FLOAT, &
                             (/dimh8, dimid(2)/), dr_3        )
      status = nf90_def_var (bc%h% ncid ,'dsr_1', NF90_FLOAT, &
                             (/dimh24,dimid(2)/), dr_1        )

      status = nf90_def_var (bc%h% ncid ,'tmax_6', NF90_FLOAT, &
                             (/dimh4, dimid(2)/), ma_6        )
      status = nf90_def_var (bc%h% ncid ,'tmin_6', NF90_FLOAT, &
                             (/dimh4, dimid(2)/), mi_6        )
      status = nf90_def_var (bc%h% ncid ,'tmax_1', NF90_FLOAT, &
                             (/dimh24,dimid(2)/), ma_1        )
      status = nf90_def_var (bc%h% ncid ,'tmin_1', NF90_FLOAT, &
                             (/dimh24,dimid(2)/), mi_1        )

      call def_var1 ('hours'     ,hours     ,'hours')
      call def_var1s('codetype'  ,codetype  ,'CMA codetype')
      call def_var1f('z_station' ,z_station ,'station height used')
      call def_var1f('lat'       ,lat       ,'last known latitude')
      call def_var1f('lon'       ,lon       ,'last known longitude')
      call def_vars ('varno'     ,varno     ,'type of the observed quantity')
      call def_var  ('n'         ,n         ,'entries')
      call def_var  ('ob'        ,ob        ,'obs-bg  mean')
      call def_var  ('oa'        ,oa        ,'obs-ana mean')
      call def_var  ('ob_ob'     ,ob_ob     ,'obs-bg  variance')
      call def_var  ('oa_oa'     ,oa_oa     ,'obs-ana variance')
      call def_var  ('n_ac'      ,n_ac      ,'entries accumulated')
      call def_var  ('ob_ac'     ,ob_ac     ,'obs-bg  mean accumulated')
      call def_var  ('oa_ac'     ,oa_ac     ,'obs-ana mean accumulated')
      call def_var  ('ob_ob_ac'  ,ob_ob_ac  ,'obs-bg  variance accumulated')
      call def_var  ('oa_oa_ac'  ,oa_oa_ac  ,'obs-ana variance accumulated')
      call def_var  ('o_err'     ,o_err     ,'nominal obs error squared')
      call def_var  ('b_err'     ,b_err     ,'nominal bg  error squared')
      call def_var  ('o_err_ac'  ,o_err_ac  ,               &
                     'nominal obs error squared accumulated')
      call def_var  ('b_err_ac'  ,b_err_ac  ,               &
                     'nominal bg  error squared accumulated')
      call def_var  ('bc_b'      ,bc_b      ,'bias correction background value')
      call def_var  ('bc_a'      ,bc_a      ,'bias correction analysis   value')
      ! If synop bias correction coefficients are used for t2m or rh2m
      ! write them to bc_t2m_b and bc_rh2m_b in bias-file
      if (bc_t2m) then
         call def_varbc('bc_t2m_b'    ,bc_t2m_b  ,'bias correction background coefficients for t2m')
         call def_varbc('bc_t2m_a'    ,bc_t2m_a  ,'bias correction analysis   coefficients for t2m')
      end if
      if (bc_rh2m) then
         call def_varbc('bc_rh2m_b'   ,bc_rh2m_b ,'bias correction background coefficients for rh2m')
         call def_varbc('bc_rh2m_a'   ,bc_rh2m_a ,'bias correction analysis   coefficients for rh2m')
      end if
!!    call def_var  ('bc_b_err ' ,bc_b_err  ,                  &
!!                   'bias correction background error squared')
!!    call def_var  ('bc_a_err ' ,bc_a_err  ,                 &
!!                   'bias correction analysis error squared')
      status = nf90_enddef  (bc%h% ncid)
      !----------------
      ! write variables
      !----------------
      allocate (statids (bc% n))
      statids = bc% plat% statid
!+++  status  = nf90_put_var (bc%h% ncid, statid, bc% plat% statid)
      status  = nf90_put_var (bc%h% ncid, statid, statids)
      deallocate (statids)

      if (wsi_mode > 0) then
         allocate (wsichar(bc% n))
         do i = 1, bc% n
            if (bc% plat(i)% wsi% valid) then
               call wsi_to_text (bc% plat(i)% wsi, wsichar(i))
            else
               wsichar(i) = ""
            end if
         end do
         status = nf90_put_var (bc%h% ncid, wsi, wsichar)
         deallocate (wsichar)
      end if

      !--------
      ! 6h data
      !--------
      allocate   (rr ( 2, bc% n))

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% rr12; enddo
      where (rr < -0.11_wp) rr = rvind
      call put_var  ('rr_12' ,rr_12 ,rr); deallocate (rr)

      allocate   (rr ( 4, bc% n))

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% rr_6; enddo
      where (rr < -0.11_wp) rr = rvind
      call put_var  ('rr_6'  ,rr_6  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gu_6; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('gust_6'  ,gu_6  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gr_6; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('gsr_6'  ,gr_6  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% dr_6; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('dsr_6'  ,dr_6  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% ma_6; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('tmax_6' ,ma_6  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% mi_6; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('tmin_6' ,mi_6  ,rr)

      deallocate (rr)
      !--------
      ! 3h data
      !--------
      allocate   (rr ( 8, bc% n))

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% rr_3; enddo
      where (rr < -0.11_wp) rr = rvind
      call put_var  ('rr_3'  ,rr_3  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gu_3; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('gust_3'  ,gu_3  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gr_3; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('gsr_3'  ,gr_3  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% dr_3; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('dsr_3'  ,dr_3  ,rr)

      deallocate (rr)
      !--------
      ! 1h data
      !--------
      allocate   (rr (24, bc% n))

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% rr_1; enddo
      where (rr < -0.11_wp) rr = rvind
      call put_var  ('rr_1'  ,rr_1  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gu_1; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('gust_1'  ,gu_1  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% gr_1; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('gsr_1'  ,gr_1  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% dr_1; enddo
      where (rr < 0.) rr = rvind
      call put_var  ('dsr_1'  ,dr_1  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% ma_1; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('tmax_1'  ,ma_1  ,rr)

      do i = 1, bc% n; rr(:,i) = bc% ag(i)% mi_1; enddo
      where (rr < 0._wp) rr = rvind
      call put_var  ('tmin_1'  ,mi_1  ,rr)

      deallocate (rr)

      call put_var1 ('hours'    ,hours    ,bc% ag  % hours)
      call put_var1 ('codetype' ,codetype ,bc% plat% code)
      where (bc% plat% z == -999._sp)
         bc% plat% z   = NF90_FILL_FLOAT
      end where
      where (bc% plat% lat == -999._sp)
         bc% plat% lat = NF90_FILL_FLOAT
         bc% plat% lon = NF90_FILL_FLOAT
      end where
      call put_var1f('z_station',z_station,bc% plat% z)
      call put_var1f('lat'      ,lat      ,bc% plat% lat)
      call put_var1f('lon'      ,lon      ,bc% plat% lon)
      call put_var1 ('varno'    ,varno    ,vvarno)
      call put_vari ('n'        ,n        ,bc% data% n)
      call put_var  ('ob'       ,ob       ,bc% data% ob)
      call put_var  ('oa'       ,oa       ,bc% data% oa)
      call put_var  ('ob_ob'    ,ob_ob    ,bc% data% ob_ob)
      call put_var  ('oa_oa'    ,oa_oa    ,bc% data% oa_oa)
      call put_var  ('n_ac'     ,n_ac     ,bc% data% n_ac)
      call put_var  ('ob_ac'    ,ob_ac    ,bc% data% ob_ac)
      call put_var  ('oa_ac'    ,oa_ac    ,bc% data% oa_ac)
      call put_var  ('ob_ob_ac' ,ob_ob_ac ,bc% data% ob_ob_ac)
      call put_var  ('oa_oa_ac' ,oa_oa_ac ,bc% data% oa_oa_ac)
      call put_var  ('o_err'    ,o_err    ,bc% data% o_err )
      call put_var  ('b_err'    ,b_err    ,bc% data% b_err )
      call put_var  ('o_err_ac' ,o_err_ac ,bc% data% o_err_ac)
      call put_var  ('b_err_ac' ,b_err_ac ,bc% data% b_err_ac)
      call put_var  ('bc_b'     ,bc_b     ,bc% data% bc_b)
      call put_var  ('bc_a'     ,bc_a     ,bc% data% bc_a)

      !----------------------------------------------------
      ! Write coefficients for t2m and rh2m bias correction
      !----------------------------------------------------
      if (bc_t2m) then
         allocate   (rr ( ncoeff, bc% n))
         do i = 1, bc% n; rr(:,i) = bc% coeff(i)% t2m_bc_b; enddo
         call put_var  ('bc_t2m_b'  ,bc_t2m_b  ,rr)
         do i = 1, bc% n; rr(:,i) = bc% coeff(i)% t2m_bc_a; enddo
         call put_var  ('bc_t2m_a'  ,bc_t2m_a  ,rr)
         deallocate (rr)
      end if
      if (bc_rh2m) then
         allocate   (rr ( ncoeff, bc% n))
         do i = 1, bc% n; rr(:,i) = bc% coeff(i)% rh2m_bc_b; enddo
         call put_var  ('bc_rh2m_b'  ,bc_rh2m_b  ,rr)
         do i = 1, bc% n; rr(:,i) = bc% coeff(i)% rh2m_bc_a; enddo
         call put_var  ('bc_rh2m_a'  ,bc_rh2m_a  ,rr)
         deallocate (rr)
      end if

      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif

  contains
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var1 (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

    end subroutine put_var1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var1f (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    real(sp)         ,intent(in) :: values (:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

    end subroutine put_var1f
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_vari (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

    end subroutine put_vari
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    real(wp)         ,intent(in) :: values (:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

    end subroutine put_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_dim (name, size, dimid)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(in)  :: size
    integer          ,intent(out) :: dimid

      status = nf90_def_dim (bc%h% ncid ,name ,size ,dimid)

      if (status /= NF90_NOERR) then
        write(0,*)   'write_synop_bc_file: def_dim '//trim(name)//' : size =',size
        call finish ('write_synop_bc_file: def_dim '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

    end subroutine def_dim
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_FLOAT, dimid, varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if
    end subroutine def_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var1 (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_INT, dimid(2), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1 '//trim(name),&
                      trim(nf90_strerror(status))                 )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1 '//trim(name),&
                      trim(nf90_strerror(status))                 )
      end if
    end subroutine def_var1
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var1f (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_FLOAT, dimid(2), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1f '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1f '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if
    end subroutine def_var1f
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var1s (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_SHORT, dimid(2), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1s '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_var1s '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if
    end subroutine def_var1s
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_vars (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_SHORT, dimid(1), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_vars '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_vars '//trim(name),&
                      trim(nf90_strerror(status))                  )
      end if
    end subroutine def_vars
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_varbc (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_FLOAT, (/dimcoef,dimid(2)/), varid)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_varbc '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_synop_bc_file: def_varbc '//trim(name),&
                      trim(nf90_strerror(status))                )
      end if
    end subroutine def_varbc
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine write_synop_bc_file

!==============================================================================

  subroutine synop_bc_bfg (obs)
  !---------------------------------------------------------------------
  ! SYNOP bias correction routine to be called before first guess check.
  ! Check for missing entries in bc file, extend file
  ! Apply bias correction to observations.
  !---------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    !----------------
    ! local variables
    !----------------
    integer ,parameter       :: mp = 30000   ! max. number of new synop stations
    type(t_tmp)              :: mis_pe (mp)  ! missing entries in bc file
    type(t_tmp)              :: mis_all(mp)  ! missing entries in bc file
    integer(i8)              :: hash         ! integer representation of name
    integer                  :: code         ! CMA codetype
    integer                  :: n            ! counter / PE
    integer                  :: m            ! counter all PEs
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j, k      ! loop indices
    integer                  :: ibc          ! bias correction index
    integer                  :: clcover      ! cloudcover (for 2tm bias correction)
    integer                  :: state        ! state to set if no BC available
    real(wp)                 :: rhour        ! time of observation
    real(wp)                 :: rawobs       ! raw observation
    real(wp)                 :: bcor         ! bias correction
    real(wp)                 :: scalf        ! scaling factor (dz to dp)
    real(wp)                 :: z0_bc        ! bias correction threshold param.
    real(wp)                 :: tmp          ! temporary
    type(t_plat_bc) ,pointer :: plat (:)     ! temporary
    type(t_ag)      ,pointer :: ag   (:)     ! temporary
    type(t_coeff)   ,pointer :: coeff(:)     ! temporary
    type(t_bc)      ,pointer :: data (:,:)   ! temporary
!   logical                  :: first        ! first obs. (z or ps) of this spot
    real(wp)                 :: M_trig(5)    ! array with trigonometric basis functions
                                             ! for T2M and RH2M bias correction (diurnal variation)
    real(wp)                 :: M_bfct(5*2)  ! array with trigonometric and polynomial basis functions
                                             ! for T2M and RH2M bias correction (diurnal variation + cloudcover)
    real(wp)                 :: v_coeff(5*2) ! array with coefficients for t2m bias correction
    real(wp) ,parameter      :: scal = 2 * pi / 24 ! 2*pi/24; scaling factor

    !---------------------------------------------------------
    ! relate bias correction file entries to reports in 3dvar.
    ! check for missing entries in bc file.
    ! extend file if necessary.
    !---------------------------------------------------------
    if (biascor_mode /= 0) then
      !-----------------------------------------------
      ! relate bias correction file entries to reports
      ! check for missing entries on this PE
      !-----------------------------------------------
      n = 0
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
spot:   do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP .and. &
              obs% o(ib)% spot(is)% hd% obstype /= OT_DRIBU       ) cycle
          obs% o(ib)% spot(is)% bc_index = 0
          !------------------------------------------
          ! Skip stations with non-unique identifiers
          !------------------------------------------
          if (obs% o(ib)% spot(is)% hd% codetype == 21 .or. &
              obs% o(ib)% spot(is)% hd% codetype == 24      ) then
             if (obs% o(ib)% spot(is)% statid    == "SHIP") cycle spot
          end if
          !-----------------------------------
          ! Check for reasonable report status
          !-----------------------------------
          if (obs% o(ib)% spot(is)% use% state <= STAT_DISMISS) cycle spot

          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
!           !+++++++++++++++++++++++++++++++++++++++++++++
!           ! TODO: check for reasonable observed quantity
!           !+++++++++++++++++++++++++++++++++++++++++++++
!           select case (obs% o(ib)% varno(i))
!           case (VN_Z, VN_PS)
!           case default
!             cycle
!           end select

!           hash = transfer (obs% o(ib)% spot(is)% statid, hash)
            hash =           obs% o(ib)% spot(is)% stat_hash
            code =           obs% o(ib)% spot(is)% hd% codetype
            !----------------------------------------------------------
            ! Some stations alternatively report as manual or automatic
            ! Choose representative codetype.
            !----------------------------------------------------------
            select case (code)
            case (11:14)
               code = 11        ! SYNOP LAND
            case (21:24)
               code = 21        ! SYNOP SHIP
!           case (140,17)       ! METAR, METAR-USA (recoded)
!           case (20)           ! coastal (CMAN)
!           case (165)          ! BUOY
            end select
            do j = 1, size(synop_bc% plat)
              if   (hash == synop_bc% plat(j)% hash) then
                if (code /= synop_bc% plat(j)% code) then
                   !---------------------------------------------------------
                   ! Assume existing entry of unknown codetype was SYNOP Land
                   ! (to keep previously existing aggregations).
                   !---------------------------------------------------------
                   if (synop_bc% plat(j)% code /= -1) cycle
                   if (                   code /= 11) cycle
                   synop_bc% plat(j)% code = code
                end if
                obs% o(ib)% spot(is)% bc_index = j
                cycle spot
              endif
            end do
            do j = 1, n
              if (mis_pe (j)% hash == hash .and. &
                  mis_pe (j)% code == code       ) then
                obs% o(ib)% spot(is)% bc_index = -j
                cycle spot
              endif
            end do
            !--------------------------------
            ! add to list of missing stations
            !--------------------------------
            n = n + 1
            if (n > mp) call finish('synop_bc_bfg','n > mp')
            obs% o(ib)% spot(is)% bc_index = - n
            mis_pe (n)% statid  = obs% o(ib)% spot(is)% statid
            mis_pe (n)% hash    = obs% o(ib)% spot(is)% stat_hash
            mis_pe (n)% wsi     = obs% o(ib)% spot(is)% wsi
            mis_pe (n)% code    = code
            mis_pe (n)% i       = n
            mis_pe (n)% j       = 0
            mis_pe (n)% lat     = obs% o(ib)% spot(is)% col% c% dlat
            mis_pe (n)% lon     = obs% o(ib)% spot(is)% col% c% dlon
            cycle spot
          end do
        end do spot
      end do
      !----------------------------
      ! Update codetype where known
      !----------------------------
      synop_bc% plat% code = p_max (synop_bc% plat% code)
      !---------------------------------------
      ! cross-check missing entries on all PEs
      !---------------------------------------
      m = p_sum (n)
      if (m > mp) call finish('synop_bc_bfg','m > mp')
      if (m > 0) then
        call p_gather_tmp (mis_pe (1:n), mis_all(1:m), dace% pio)
        if (dace% lpio) then
          k             = 1
          i             = 1
          mis_all(i)% j = k
          mis_pe (k)    = mis_all(i)
outer:    do i = 2, m
            do j = 1, k
              if (mis_all(i)% hash == mis_pe(j)% hash .and. &
                  mis_all(i)% code == mis_pe(j)% code       ) then
                mis_all(i)% j = j
                cycle outer
              endif
            end do
            k = k + 1
            mis_all(i)% j = k
            mis_pe (k)    = mis_all(i)
          end do outer
          !------------
          ! extend file
          !------------
          i = size (synop_bc% plat)
          plat => synop_bc% plat
          ag   => synop_bc% ag
          if (bc_t2m .or. bc_rh2m) coeff => synop_bc% coeff
          data => synop_bc% data
          allocate (synop_bc% plat (    i+k))
          allocate (synop_bc% ag   (    i+k))
          if (bc_t2m .or. bc_rh2m) allocate (synop_bc% coeff(    i+k))
          allocate (synop_bc% data (nob,i+k))
          synop_bc% n = i+k
          synop_bc% plat(  :i) = plat
          synop_bc% ag  (  :i) = ag
          if (bc_t2m .or. bc_rh2m) synop_bc% coeff(  :i) = coeff
          synop_bc% data(:,:i) = data
          do j = 1, k
            synop_bc% plat (i+j)% hash   = mis_pe (j)% hash
            synop_bc% plat (i+j)% statid = mis_pe (j)% statid
            synop_bc% plat (i+j)% wsi    = mis_pe (j)% wsi
            synop_bc% plat (i+j)% code   = mis_pe (j)% code
            synop_bc% ag   (i+j)% ibc    = i+j
            synop_bc% plat (i+j)% lat    = mis_pe (j)% lat
            synop_bc% plat (i+j)% lon    = mis_pe (j)% lon
          end do
          deallocate (plat)
          deallocate (ag  )
          if (bc_t2m .or. bc_rh2m) deallocate (coeff)
          deallocate (data)
          mis_all (1:m)% j = mis_all (1:m)% j + i
        endif
        !--------------------------------------
        ! redistribute information to other PEs
        !--------------------------------------
        call p_scatter_tmp (mis_all(1:m), mis_pe(1:n), dace% pio)
        if (dace% lpio) m = i+k
        call p_bcast (m, dace% pio)
        if (.not.dace% lpio) then
          deallocate (synop_bc% plat)
          deallocate (synop_bc% ag  )
          if (bc_t2m .or. bc_rh2m) deallocate (synop_bc% coeff)
          deallocate (synop_bc% data)
          allocate   (synop_bc% plat (    m))
          allocate   (synop_bc% ag   (    m))
          if (bc_t2m .or. bc_rh2m) allocate   (synop_bc% coeff(    m))
          allocate   (synop_bc% data (nob,m))
        endif
        call p_bcast_plat (synop_bc% plat, dace% pio)
        call p_bcast_ag   (synop_bc% ag,   dace% pio)
        if (bc_t2m .or. bc_rh2m) call p_bcast_coeff(synop_bc% coeff,dace% pio)
        call p_bcast_data (synop_bc% data, dace% pio)
        synop_bc% n = m
        !-----------------------------------
        ! relate reports to new file entries
        !-----------------------------------
        do ib = 1, size (obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP .and. &
                obs% o(ib)% spot(is)% hd% obstype /= OT_DRIBU       ) cycle
            if (obs% o(ib)% spot(is)% bc_index < 0)      &
                obs% o(ib)% spot(is)% bc_index = mis_pe( &
               -obs% o(ib)% spot(is)% bc_index)% j
          end do
        end do
      endif
      !---------------------
      ! aggregate rain rates
      !---------------------
      call aggregate_rr (obs)
    endif
    !---------------------------------------
    ! Temporarily "compress" aggregated data
    !---------------------------------------
    call compress (synop_bc)

    !-----------------------------------
    ! apply bias-correction within 3dvar
    !-----------------------------------
    if (biascor_mode >= BC_FG) then
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP .and. &
              obs% o(ib)% spot(is)% hd% obstype /= OT_DRIBU       ) cycle

          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                      cycle
          !--------------------------------------
          ! Try to remember "last known position"
          !--------------------------------------
          synop_bc% plat(ibc)% lat = obs% o(ib)% spot(is)% col% c% dlat
          synop_bc% plat(ibc)% lon = obs% o(ib)% spot(is)% col% c% dlon

          select case (obs% o(ib)% spot(is)% hd% codetype)
          case (11,14,20)       ! SYNOP Land (incl. coast)
             if (.not. bc_synop)                             cycle
             z0_bc = z0_bc_synop
          case (21,24)          ! SHIP
             if (.not. bc_ship)                              cycle
             z0_bc = z0_bc_ship
          case (140,17)         ! METAR, METAR-USA
             if (.not. bc_metar)                             cycle
             z0_bc = z0_bc_metar
          case (165)            ! DRIBU
             if (.not. bc_buoy)                              cycle
             z0_bc = z0_bc_buoy
          case default          ! not implemented
             cycle
          end select

          !-----------------------------------------------------------------
          ! Retrieve cloudcover information for T2M and RH2M bias correction
          !-----------------------------------------------------------------
          ! Check for cloudcover
          ! If no cloudcover available use cloudcover 8
          clcover = 8
          do k = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
             select case (obs% o(ib)% varno(k))
             case (VN_N)
                clcover = obs% o(ib)% body(k)% o
                ! Write cloudcover information in spot% phase for verification
                obs% o(ib)% spot(is)% phase = clcover
                exit
             case default
                cycle
             end select
          end do

          !--------------------------------------------------------------------------
          ! Construct basis functions for T2M and RH2M bias correction
          ! (5 trigonometric for diurnal variations and 2 polynomial for cloud cover)
          !--------------------------------------------------------------------------
          ! Array M_trig with trigonometric basis functions.
          rhour = ihh (obs% o(ib)% spot(is)% hd% time)
          M_trig(1) = 1
          do k = 1, 2
             M_trig(2*k)   = sin(k*scal*rhour)      ! scaling with scal = 2*pi/24
             M_trig(2*k+1) = cos(k*scal*rhour)
          end do
          ! extended by polynomial basis functions for cloud cover
          M_bfct = (/ M_trig, (9-clcover) * M_trig /)

!         first = .true.
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            if (obs% o(ib)% body(i)% use% state <= STAT_DISMISS) cycle
            select case (obs% o(ib)% varno(i))
            case (VN_Z)
               scalf = 1._wp
               j = OB_PS
            case (VN_PS)
               scalf = obs% o(ib)% spot(is)% pz_bg
               j = OB_PS
            case (VN_T2M)  ! case for T2M bias correction
               j = OB_T2M
            case (VN_RH2M) ! case for RH2M bias correction
               j = OB_RH2M
            case default
               cycle
            end select

            if      (synop_bc% plat(ibc)% z == -999._sp               ) then
               synop_bc% plat(ibc)% z   = obs% o(ib)% spot(is)% z
!!             write(0,*) dace% pe, "synop_bc_bfg: station, codetype, height: ", &
!!                  obs% o(ib)% spot(is)% statid, " : ", &
!!                  obs% o(ib)% spot(is)% hd% codetype, &
!!                  obs% o(ib)% spot(is)% z
            else if (synop_bc% plat(ibc)% z /= obs% o(ib)% spot(is)% z) then
               write(0,*) "synop_bc_bfg: WARNING: non-constant station height: ", &
                    obs% o(ib)% spot(is)% statid, " : ", real (obs% o(ib)% spot(is)% z), &
                    " /= ", real (synop_bc% plat(ibc)% z), "  codetypes:", &
                    obs% o(ib)% spot(is)% hd% codetype, synop_bc% plat(ibc)% code
!!write(0,*) "synop_bc_bfg: ", obs% o(ib)% spot(is)% statid, obs% o(ib)% spot(is)% z, &
!!                  obs% o(ib)% varno(i),   obs% o(ib)% body(i)% lev_typ, &
!!                  obs% o(ib)% body(i)% o, obs% o(ib)% body(i)% use% state
               !--------------------------
               !+++ TODO: reset statistics
               !--------------------------
            end if

            rawobs = obs% o(ib)% body(i)% o  - &
                     obs% o(ib)% body(i)% bc
            !------------------------------------------------------------------
            ! We currently use the same n_required for every observed quantity:
            !------------------------------------------------------------------
            if (synop_bc% data(j,ibc)% n_ac < n_required) then
               bcor = 0._wp
               if (.not. (obs% o(ib)% varno(i) == VN_Z .or.   &
                          obs% o(ib)% varno(i) == VN_PS       ) .or.  z0_bc /= 999.) then
                  state = rept_use (obs% o(ib)% spot(is)% hd% obstype)% use (CHK_BIASCOR)
                  call decr_use (obs% o(ib)% body(i)% use, &
                                 check = CHK_BIASCOR,      &
                                 state = state,            &
                                 lflag = .true.            )
               endif
            else
               select case (obs% o(ib)% varno(i))
               case (VN_Z,VN_PS)
                  !--------------------------------------------------------
                  ! If modulus of bias is below threshold, treat as 0, else
                  ! derive effective bias correction (subtract threshold)
                  !--------------------------------------------------------
                  bcor = synop_bc% data(j,ibc)% bc_b
                  tmp  = max (abs (bcor) - z0_bc, 0._wp)
                  bcor = sign (tmp, bcor)
                  bcor = bcor * scalf
               case (VN_T2M)
                  if (bc_t2m) then
                     v_coeff = synop_bc% coeff(ibc)% t2m_bc_b
                  else
                     v_coeff = 0
                  end if
                  ! calculate bias correction with basis functions for trig. interpolation
                  bcor = dot_product ( M_bfct, v_coeff )
                  ! Store calculated bias correction in bc_b
                  synop_bc% data(j,ibc)% bc_b = bcor
               case (VN_RH2M)
                  if (bc_rh2m) then
                     v_coeff = synop_bc% coeff(ibc)% rh2m_bc_b
                  else
                     v_coeff = 0
                  end if
                  ! calculate bias correction with basis functions for trig. interpolation
                  bcor = dot_product ( M_bfct, v_coeff )
                  ! Store calculated bias correction in bc_b
                  synop_bc% data(j,ibc)% bc_b = bcor
                  ! In case of rh2m close to 1: reduce bcor if necessary
                  bcor = max(min(bcor, 1 - rawobs), rawobs - 1)
               end select
            end if
            obs% o(ib)% body(i)% bc =        - bcor
            obs% o(ib)% body(i)% o  = rawobs - bcor
          end do
        end do
      end do
    endif

  end subroutine synop_bc_bfg

!------------------------------------------------------------------------------

  subroutine aggregate_rr (obs)
  !---------------------
  ! aggregate rain rates
  !---------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    integer                 :: i, j, k, n, m
    integer                 :: ib        ! box index
    integer                 :: is        ! spot index
    integer                 :: ibc       ! plane index
    real(wp)                :: rhour     ! time of observation
    integer                 :: ihour     ! time of observation
    integer                 :: dhour     ! offset to last storage
    integer                 :: jhour     ! hour of the day
    integer                 :: j_1, j_3, j_6, j_12
    integer                 :: g_1, g_3, g_6
    integer                 :: a_1, a_6  ! indices: tmax
    integer                 :: i_1, i_6  ! indices: tmin
    integer                 :: s_1, s_3, s_6
    integer                 :: d_1, d_3, d_6
    logical                 :: lrr, lgu  ! flag present entries
    logical                 :: lma, lmi  ! flag present entries
    logical                 :: lrd
    logical    ,allocatable :: mask(:)
    type(t_ag) ,allocatable :: agl(:), agg(:)
    type(t_ag) ,pointer     :: ag        !
    real(sp)                :: rrr       ! temporary
    real(sp)                :: gsr, dsr
    integer                 :: nrr, jrr
    integer                 :: ngr, ndr
    real(sp)                :: rgu       ! temporary
    real(sp)                :: rmi, rma  ! temporaries
    integer                 :: ngu
    integer                 :: state     ! observation state to set
    !-----------------------------------
    ! cycle over SYNOP land observations
    !-----------------------------------
    allocate (mask (size (synop_bc% ag))); mask = .false.
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP) cycle
        !-----------------------------------------------------
        ! Skip ships for the time being (for multiple reasons)
        !-----------------------------------------------------
        if (obs% o(ib)% spot(is)% hd% codetype == 21 .or. &
            obs% o(ib)% spot(is)% hd% codetype == 24      ) cycle

        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                      cycle
        ag => synop_bc% ag(ibc)

        !-------------------
        ! derive time offset
        !-------------------
!       rhour = hours (obs% o(ib)% spot(is)% actual_time)
        rhour = hours (obs% o(ib)% spot(is)% hd% time)
        ihour = nint  (rhour)
        dhour = ihour - ag% hours
        jhour = modulo (ihour, 24)
!       if ((rhour - ihour) /= 0._wp) cycle
!       if (abs (rhour - ihour) > 0.16_wp) cycle  ! allow up to ~9 minutes (US BUFR)
        if (abs (rhour - ihour) > 0.34_wp) cycle  ! allow up to ~20 minutes (Hungary)

        !-----------------------------
        ! check for rain-rates present
        !-----------------------------
        j_1 = 0; j_3 = 0; j_6 = 0; j_12 = 0; lrr = .false.
        g_1 = 0; g_3 = 0; g_6 = 0;           lgu = .false.
        a_1 = 0;          a_6 = 0;           lma = .false.
        i_1 = 0;          i_6 = 0;           lmi = .false.
        s_1 = 0; s_3 = 0; s_6 = 0;           lrd = .false.
        d_1 = 0; d_3 = 0; d_6 = 0

        do i = 1, obs% o(ib)% spot(is)% o% n
          j = i + obs% o(ib)% spot(is)% o% i

          select case (obs% o(ib)% varno(j))
          case (VN_RR)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                j_1  = j
                lrr  = .true.
              case ( 3)
                if (mod (jhour, 3) /= 0) cycle
                j_3  = j
                lrr  = .true.
              case ( 6)
                if (mod (jhour, 6) /= 0) cycle
                j_6  = j
                lrr  = .true.
              case (12)
                if (mod (jhour, 12) /= 0) cycle
                j_12 = j
                lrr  = .true.
              end select
            endif

          case (VN_RAD_GL)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                s_1  = j
                lrd  = .true.
              case ( 3)
                if (mod (jhour, 3) /= 0) cycle
                s_3  = j
                lrd  = .true.
              case ( 6)
                if (mod (jhour, 6) /= 0) cycle
                s_6  = j
                lrd  = .true.
              end select
            endif

          case (VN_RAD_DF)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                d_1  = j
                lrd  = .true.
              case ( 3)
                if (mod (jhour, 3) /= 0) cycle
                d_3  = j
                lrd  = .true.
              case ( 6)
                if (mod (jhour, 6) /= 0) cycle
                d_6  = j
                lrd  = .true.
              end select
            endif

          case (VN_GUST)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                g_1  = j
                lgu  = .true.
              case ( 3)
                if (mod (jhour, 3) /= 0) cycle
                g_3  = j
                lgu  = .true.
              case ( 6)
                if (mod (jhour, 6) /= 0) cycle
                g_6  = j
                lgu  = .true.
              end select
            endif

          case (VN_TMAX)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                a_1  = j
                lma  = .true.
              case (12)                         ! Aggregated 12h value every 6h
                if (mod (jhour, 6) /= 0) cycle
                lma  = .true.
              end select
            endif

          case (VN_TMIN)
            if (obs% o(ib)% body (j)% o /= rvind) then
              select case (nint(obs% o(ib)% olev(j)))
              case ( 1)
                i_1  = j
                lmi  = .true.
              case (12)                         ! Aggregated 12h value every 6h
                if (mod (jhour, 6) /= 0) cycle
                lmi  = .true.
              end select
            endif
          end select
        end do

        !------------------------
        ! remove outdated entries
        !------------------------
        if (.not. (lrr .or. lrd .or. lgu .or. lma .or. lmi)) cycle
        mask (ibc) = .true.
        do i = 0, min (dhour, 24)-1
          j = modulo (jhour -i, 24)
                               ag% rr_1 (j    ) = -99._sp
          if (mod (j, 3) == 0) ag% rr_3 (j / 3) = -99._sp
          if (mod (j, 6) == 0) ag% rr_6 (j / 6) = -99._sp
          if (mod (j,12) == 0) ag% rr12 (j /12) = -99._sp

                               ag% gr_1 (j    ) = -99._sp
          if (mod (j, 3) == 0) ag% gr_3 (j / 3) = -99._sp
          if (mod (j, 6) == 0) ag% gr_6 (j / 6) = -99._sp

                               ag% dr_1 (j    ) = -99._sp
          if (mod (j, 3) == 0) ag% dr_3 (j / 3) = -99._sp
          if (mod (j, 6) == 0) ag% dr_6 (j / 6) = -99._sp

                               ag% gu_1 (j    ) = -99._sp
          if (mod (j, 3) == 0) ag% gu_3 (j / 3) = -99._sp
          if (mod (j, 6) == 0) ag% gu_6 (j / 6) = -99._sp

                               ag% ma_1 (j    ) = -99._sp
                               ag% mi_1 (j    ) = -99._sp
          if (mod (j, 6) == 0) ag% ma_6 (j / 6) = -99._sp
          if (mod (j, 6) == 0) ag% mi_6 (j / 6) = -99._sp
        end do

        !-------------------
        ! insert new entries
        !-------------------
        j = jhour
                                         ag% hours = ihour
        if (j_1 /= 0)                    ag% rr_1 (j   ) = max (obs% o(ib)% body (j_1 )% o, 0._sp)
        if (j_3 /= 0 .and. mod(j, 3)==0) ag% rr_3 (j/ 3) = max (obs% o(ib)% body (j_3 )% o, 0._sp)
        if (j_6 /= 0 .and. mod(j, 6)==0) ag% rr_6 (j/ 6) = max (obs% o(ib)% body (j_6 )% o, 0._sp)
        if (j_12/= 0 .and. mod(j,12)==0) ag% rr12 (j/12) = max (obs% o(ib)% body (j_12)% o, 0._sp)

        if (s_1 /= 0)                    ag% gr_1 (j   ) =      obs% o(ib)% body (s_1 )% o
        if (s_3 /= 0 .and. mod(j, 3)==0) ag% gr_3 (j/ 3) =      obs% o(ib)% body (s_3 )% o
        if (s_6 /= 0 .and. mod(j, 6)==0) ag% gr_6 (j/ 6) =      obs% o(ib)% body (s_6 )% o

        if (d_1 /= 0)                    ag% dr_1 (j   ) =      obs% o(ib)% body (d_1 )% o
        if (d_3 /= 0 .and. mod(j, 3)==0) ag% dr_3 (j/ 3) =      obs% o(ib)% body (d_3 )% o
        if (d_6 /= 0 .and. mod(j, 6)==0) ag% dr_6 (j/ 6) =      obs% o(ib)% body (d_6 )% o

        if (g_1 /= 0)                    ag% gu_1 (j   ) =      obs% o(ib)% body (g_1 )% o
        if (g_3 /= 0 .and. mod(j, 3)==0) ag% gu_3 (j/ 3) =      obs% o(ib)% body (g_3 )% o
        if (g_6 /= 0 .and. mod(j, 6)==0) ag% gu_6 (j/ 6) =      obs% o(ib)% body (g_6 )% o
        if (a_1 /= 0)                    ag% ma_1 (j   ) =      obs% o(ib)% body (a_1 )% o
        if (i_1 /= 0)                    ag% mi_1 (j   ) =      obs% o(ib)% body (i_1 )% o
      end do
    end do

    !---------------------
    ! send results to PE 0
    !---------------------
    n = count (mask)
    allocate (agl (n))
    if (n > 0) agl(:) = pack (synop_bc% ag, mask)
    m = p_sum (n)
    if (.not.dace% lpio) m = 0
    allocate (agg (m))
    call p_gather_ag (agl, agg, dace% pio)
    deallocate (mask)
    !---------------------------------------
    ! on PE 0 combine changes from other PEs
    !---------------------------------------
    if (dace% lpio) then
      allocate (mask (size (synop_bc% ag)))
      mask = .false.
      do k = 1, m
        ibc = agg (k)% ibc
        ag => synop_bc% ag(ibc)
        if (mask(ibc)) then
          !---------------------------------------------
          ! multiple changes: remove outdated time slots
          !---------------------------------------------
          ihour = agg (k)% hours
          dhour = ihour - ag% hours
          jhour = modulo (ihour, 24)
          do i = 0, min (dhour, 24)-1
            j = modulo (jhour -i, 24)
                                 ag% rr_1 (j    ) = -99._sp
            if (mod (j, 3) == 0) ag% rr_3 (j / 3) = -99._sp
            if (mod (j, 6) == 0) ag% rr_6 (j / 6) = -99._sp
            if (mod (j,12) == 0) ag% rr12 (j /12) = -99._sp
                                 ag% gu_1 (j    ) = -99._sp
            if (mod (j, 3) == 0) ag% gu_3 (j / 3) = -99._sp
            if (mod (j, 6) == 0) ag% gu_6 (j / 6) = -99._sp
                                 ag% ma_1 (j    ) = -99._sp
                                 ag% mi_1 (j    ) = -99._sp
            if (mod (j, 6) == 0) ag% ma_6 (j / 6) = -99._sp
            if (mod (j, 6) == 0) ag% mi_6 (j / 6) = -99._sp
          end do
          !----------------
          ! combine changes
          !----------------
          j = jhour
                            ag% rr_1 (j   ) = max (ag% rr_1 (j   ), agg(k)% rr_1 (j   ))
          if (mod(j, 3)==0) ag% rr_3 (j/ 3) = max (ag% rr_3 (j/ 3), agg(k)% rr_3 (j/ 3))
          if (mod(j, 6)==0) ag% rr_6 (j/ 6) = max (ag% rr_6 (j/ 6), agg(k)% rr_6 (j/ 6))
          if (mod(j,12)==0) ag% rr12 (j/12) = max (ag% rr12 (j/12), agg(k)% rr12 (j/12))
                            ag% gu_1 (j   ) = max (ag% gu_1 (j   ), agg(k)% gu_1 (j   ))
          if (mod(j, 3)==0) ag% gu_3 (j/ 3) = max (ag% gu_3 (j/ 3), agg(k)% gu_3 (j/ 3))
          if (mod(j, 6)==0) ag% gu_6 (j/ 6) = max (ag% gu_6 (j/ 6), agg(k)% gu_6 (j/ 6))
                            ag% ma_1 (j   ) = max (ag% ma_1 (j   ), agg(k)% ma_1 (j   ))
                            ag% mi_1 (j   ) = max (ag% mi_1 (j   ), agg(k)% mi_1 (j   ))
!         if (mod(j, 6)==0) ag% ma_6 (j/ 6) = max (ag% ma_6 (j/ 6), agg(k)% ma_6 (j/ 6))
!         if (mod(j, 6)==0) ag% mi_6 (j/ 6) = max (ag% mi_6 (j/ 6), agg(k)% mi_6 (j/ 6))
                            ag% hours = ihour
        else
          !-------------------------------
          ! first change so far: just copy
          !-------------------------------
          ag = agg (k)
          mask (ibc) = .true.
        endif
      end do
      deallocate (mask)
      !-----------------------
      ! derive missing entries
      !-----------------------
      rhour = hours (ana_time)
      ihour = nint  (rhour)
      if ((rhour - ihour) == 0._wp) then
      jhour = modulo (ihour, 24)

        do i = 1, size (synop_bc% ag)
          ag => synop_bc% ag(i)
          if (ag% hours >= ihour) then
            if (modulo (ihour, 3) == 0) then

              !----------------------
              ! derive rr_3 from rr_1 (rain)
              ! derive gr_3 from gr_1 (pyranometer)
              ! derive dr_3 from dr_1 (pyranometer)
              ! derive gu_3 from gu_1 (gust)
              !----------------------
              nrr = 0
              rrr = 0._sp
              ngr = 0
              gsr = 0._sp
              ndr = 0
              dsr = 0._sp
              ngu = 0
              rgu = 0._sp
              do k = 0, 2
                j = modulo (ihour-k, 24)

                if (ag% rr_1 (j) >= 0._sp) then
                  rrr = rrr + ag% rr_1 (j)
                  nrr = nrr + 1
                else
                  jrr = j
                endif

                if (ag% gr_1 (j) >= 0._sp) then
                  gsr = gsr + ag% gr_1 (j)
                  ngr = ngr + 1
                endif

                if (ag% dr_1 (j) >= 0._sp) then
                  dsr = dsr + ag% dr_1 (j)
                  ndr = ndr + 1
                endif

                if (ag% gu_1 (j) >= 0._sp) then
                  rgu = max (rgu, ag% gu_1 (j))
                  ngu = ngu + 1
                endif
              end do

              j = jhour
              if (nrr == 3 .and. ag% rr_3 (j/3) < 0._sp) ag% rr_3 (j/3) = rrr
              if (ngr == 3 .and. ag% gr_3 (j/3) < 0._sp) ag% gr_3 (j/3) = gsr
              if (ndr == 3 .and. ag% dr_3 (j/3) < 0._sp) ag% dr_3 (j/3) = dsr
              if (ngu == 3 .and. ag% gu_3 (j/3) < 0._sp) ag% gu_3 (j/3) = rgu

              if (modulo (ihour, 6) == 0) then

                !----------------------
                ! derive rr_6 from rr_3
                ! derive gr_6 from gr_3
                ! derive dr_6 from dr_3
                ! derive gu_6 from gu_3
                !----------------------
                k = modulo (j-3,24)
                if (ag% rr_6 (j/6) <  0._sp .and. &
                    ag% rr_3 (j/3) >= 0._sp .and. &
                    ag% rr_3 (k/3) >= 0._sp       ) then
                    ag% rr_6 (j/6) = ag% rr_3 (j/3) + ag% rr_3 (k/ 3)
                endif
                if (ag% gr_6 (j/6) <  0._sp .and. &
                    ag% gr_3 (j/3) >= 0._sp .and. &
                    ag% gr_3 (k/3) >= 0._sp       ) then
                    ag% gr_6 (j/6) = ag% gr_3 (j/3) + ag% gr_3 (k/ 3)
                endif
                if (ag% dr_6 (j/6) <  0._sp .and. &
                    ag% dr_3 (j/3) >= 0._sp .and. &
                    ag% dr_3 (k/3) >= 0._sp       ) then
                    ag% dr_6 (j/6) = ag% dr_3 (j/3) + ag% dr_3 (k/ 3)
                endif
                if (ag% gu_6 (j/6) <  0._sp .and. &
                    ag% gu_3 (j/3) >= 0._sp .and. &
                    ag% gu_3 (k/3) >= 0._sp       ) then
                    ag% gu_6 (j/6) = max (ag% gu_3 (j/3), ag% gu_3 (k/ 3))
                endif

                !----------------------
                ! derive mi_6 from mi_1
                ! derive ma_6 from ma_1
                !----------------------
                k = modulo (j-6,24) / 6
                if (ag% mi_6 (j/6) <  0._sp .and. &
                    ag% mi_1 (j  ) >  0._sp       ) then
                   rmi = minval (ag% mi_1(k+1:k+5))
                   if (rmi > 0._sp) ag% mi_6 (j/6) = min (rmi, ag% mi_1 (j))
                end if
                if (ag% ma_6 (j/6) <  0._sp .and. &
                    ag% ma_1 (j  ) >  0._sp       ) then
                   rmi = minval (ag% ma_1(k+1:k+5))
                   rma = maxval (ag% ma_1(k+1:k+5))
                   if (rmi > 0._sp) ag% ma_6 (j/6) = max (rma, ag% ma_1 (j))
                end if

                if (modulo (ihour, 12) == 0) then

                  k = modulo (j-6,24)
                  !----------------------
                  ! derive rr12 from rr_6
                  !----------------------
                  if (ag% rr12 (j/12) <  0._sp .and. &
                      ag% rr_6 (j/ 6) >= 0._sp .and. &
                      ag% rr_6 (k/ 6) >= 0._sp       ) then
                    ag% rr12 (j/12) = ag% rr_6 (j/ 6) + ag% rr_6 (k/ 6)
                  endif
                  !----------------------
                  ! derive rr_6 from rr12
                  !----------------------
                  if (ag% rr12 (j/12) >= 0._sp .and. &
                      ag% rr_6 (j/ 6) <  0._sp .and. &
                      ag% rr_6 (k/ 6) >= 0._sp       ) then
                      ag% rr_6 (j/ 6) = ag% rr12 (j/12) - ag% rr_6 (k/ 6)
                  endif

                endif                     ! modulo (ihour, 12) == 0
                !----------------------
                ! derive rr_3 from rr_6
                !----------------------
                k = modulo (j-3,24)
                if (ag% rr_6 (j/6) >= 0._sp .and. &
                    ag% rr_3 (j/3) <  0._sp .and. &
                    ag% rr_3 (k/3) >= 0._sp       ) then
                    ag% rr_3 (j/3) = ag% rr_6 (j/6) - ag% rr_3 (k/ 3)
                endif

              endif                       ! modulo (ihour, 6) == 0

              !----------------------
              ! derive rr_1 from rr_3
              !----------------------
              if (nrr == 2 .and. ag% rr_3 (j/3) >= 0._sp) ag% rr_1 (jrr) = ag% rr_3 (j/3) - rrr

            endif                         ! modulo (ihour, 3) == 0
          endif
        end do
      endif
    endif
    !-----------------------------
    ! broadcast changes to all PEs
    !-----------------------------
    call p_bcast_ag (synop_bc% ag, dace% pio)
    deallocate (agl, agg)

    !------------------------------------------
    ! store derived entries in observation data
    !------------------------------------------
    if (aggregate) then
      !-----------------------------------
      ! cycle over SYNOP land observations
      !-----------------------------------
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP) cycle
          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                      cycle
          ag => synop_bc% ag(ibc)
          !-------------------
          ! derive time offset
          !-------------------
!         rhour = hours (obs% o(ib)% spot(is)% actual_time)
          rhour = hours (obs% o(ib)% spot(is)% hd% time)
          ihour = nint  (rhour)
          jhour = modulo (ihour, 24)
!         if ((rhour - ihour) /= 0._wp) cycle
!         if (abs (rhour - ihour) > 0.16_wp) cycle  ! allow up to ~9 minutes (US BUFR)
          if (abs (rhour - ihour) > 0.34_wp) cycle  ! allow up to ~20 minutes (Hungary)
          !-----------------------------
          ! check for rain-rates missing
          !-----------------------------
          do i = 1, obs% o(ib)% spot(is)% o% n
            j = i + obs% o(ib)% spot(is)% o% i
            if (obs% o(ib)% body (j)% o == rvind) then
              !-----------------------
              ! derive missing entries
              !-----------------------
              state = STAT_DISMISS
              select case (obs% o(ib)% varno(j))
              case default
                cycle

              case (VN_RR)
                select case (nint(obs% o(ib)% olev(j)))
                case ( 1)
                  if (ag% rr_1(jhour) >= 0._sp) then
                    obs% o(ib)% body (j)% o = ag% rr_1 (jhour)
                    state = STAT_PASSIVE
                  endif
                case ( 3)
                  if (mod (jhour, 3) == 0) then
                    if (ag% rr_3(jhour/3) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% rr_3 (jhour/3)
                      state = STAT_PASSIVE
                    endif
                  endif
                case ( 6)
                  if (mod (jhour, 6) == 0) then
                    if (ag% rr_6(jhour/6) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% rr_6 (jhour/6)
                      state = STAT_PASSIVE
                    endif
                  endif
                case (12)
                  if (mod (jhour, 12) == 0) then
                    if (ag% rr12(jhour/12) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% rr12 (jhour/12)
                      state = STAT_PASSIVE
                    endif
                  endif
                end select

              case (VN_RAD_GL)
                select case (nint(obs% o(ib)% olev(j)))
                case ( 1)
                  if (ag% gr_1(jhour) >= 0._sp) then
                    obs% o(ib)% body (j)% o = ag% gr_1 (jhour)
                    state = STAT_PASSIVE
                  endif
                case ( 3)
                  if (mod (jhour, 3) == 0) then
                    if (ag% gr_3(jhour/3) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% gr_3 (jhour/3)
                      state = STAT_PASSIVE
                    endif
                  endif
                case ( 6)
                  if (mod (jhour, 6) == 0) then
                    if (ag% gr_6(jhour/6) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% gr_6 (jhour/6)
                      state = STAT_PASSIVE
                    endif
                  endif
                end select

              case (VN_RAD_DF)
                select case (nint(obs% o(ib)% olev(j)))
                case ( 1)
                  if (ag% dr_1(jhour) >= 0._sp) then
                    obs% o(ib)% body (j)% o = ag% dr_1 (jhour)
                    state = STAT_PASSIVE
                  endif
                case ( 3)
                  if (mod (jhour, 3) == 0) then
                    if (ag% dr_3(jhour/3) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% dr_3 (jhour/3)
                      state = STAT_PASSIVE
                    endif
                  endif
                case ( 6)
                  if (mod (jhour, 6) == 0) then
                    if (ag% dr_6(jhour/6) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% dr_6 (jhour/6)
                      state = STAT_PASSIVE
                    endif
                  endif
                end select

              case (VN_GUST)
                select case (nint(obs% o(ib)% olev(j)))
                case ( 3)
                  if (mod (jhour, 3) == 0) then
                    if (ag% gu_3(jhour/3) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% gu_3 (jhour/3)
                      state = STAT_PASSIVE
                    endif
                  endif
                case ( 6)
                  if (mod (jhour, 6) == 0) then
                    if (ag% gu_6(jhour/6) >= 0._sp) then
                      obs% o(ib)% body (j)% o = ag% gu_6 (jhour/6)
                      state = STAT_PASSIVE
                    endif
                  endif
                end select

              case (VN_TMAX)
                select case (nint(obs% o(ib)% olev(j)))
                case (12)
                  if (mod (jhour, 6) == 0) then
                     k = modulo (jhour-6, 24)
                     if (obs% o(ib)% body (j)% o == rvind .and. &
                         ag% ma_6(jhour/6)       >  0._sp .and. &
                         ag% ma_6(    k/6)       >  0._sp       ) then
                        obs% o(ib)% body (j)% o = max (ag% ma_6(jhour/6), &
                                                       ag% ma_6(    k/6)  )
                        state = STAT_PASSIVE
                     end if
                  endif
                end select

              case (VN_TMIN)
                select case (nint(obs% o(ib)% olev(j)))
                case (12)
                  if (mod (jhour, 6) == 0) then
                     k = modulo (jhour-6, 24)
                     if (obs% o(ib)% body (j)% o == rvind .and. &
                         ag% mi_6(jhour/6)       >  0._sp .and. &
                         ag% mi_6(    k/6)       >  0._sp       ) then
                        obs% o(ib)% body (j)% o = min (ag% mi_6(jhour/6), &
                                                       ag% mi_6(    k/6)  )
                        state = STAT_PASSIVE
                     end if
                  endif
                end select
              end select

              call decr_use (obs% o(ib)% body(j)% use, &
                             check = CHK_NO_OBS,       &
                             state = state,            &
                             lflag = .true.            )
              obs% o(ib)% body(j)% use% check = CHK_NO_OBS
            endif
          end do
        end do
      end do
    endif

  end subroutine aggregate_rr

!==============================================================================

  subroutine synop_bc_aan (obs, y)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  type (t_vector)  ,intent(in)    :: y    ! background, observation space

    integer             :: i, j, k   ! index
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ibc       ! plane index
    real(wp)            :: bg        ! background
    real(wp)            :: an        ! analysis
    real(wp)            :: o         ! raw observation
    real(wp)            :: dob       ! obs - bg
    real(wp)            :: doa       ! obs - ana
    real(wp)            :: eo        ! observation error
    real(wp)            :: eb        ! background  error
    real(wp)            :: dt        ! time since last update of statistics
    real(wp)            :: f         ! weight for accumulated statistics
    real(wp)            :: scali     ! inverse scaling factor (dp to dz)
    integer             :: ier       ! error return flag
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp  (:,:) ! pointer to statistics file entries
    type(t_bc), pointer :: p         ! pointer to statistics file entries
    logical             :: first     ! first obs. (z or ps) of this spot?
    logical             :: count     ! count obs in statistics
    logical             :: cor       ! bias corrected?
    integer             :: np_synop  ! number of synop land processed
    integer             :: nc_synop  ! number of synop land corrected
    integer             :: np_ship   ! number of synop ship processed
    integer             :: nc_ship   ! number of synop ship corrected
    integer             :: np_buoy   ! number of dribu      processed
    integer             :: nc_buoy   ! number of dribu      corrected
    integer             :: np_metar  ! number of metar      processed
    integer             :: nc_metar  ! number of metar      corrected
    integer             :: ns        ! number of stations in statistics
    real, allocatable   :: tmpz  (:) ! latest station height
    real, allocatable   :: tmplon(:) ! latest longitude
    real, allocatable   :: tmplat(:) ! latest latitude
    real(wp)            :: M_trig(5) ! array with trigonometric basis
                                     ! functions for T2M bias correction
    real(wp)            :: M_bfct(5*2)  ! array with trigonometric and polynomial basis functions
                                        ! for T2M and RH2M bias correction (diurnal variation + cloudcover)
    real(wp) ,parameter :: scal = 2 * pi / 24 ! 2*pi/24; scaling factor
    integer             :: clcover      ! cloudcover (for t2m and rh2m bias correction)
    real(wp)            :: rhour        ! time of observation
    real(wp)            :: dbias        ! delta bias for update of coeff. for t2m bias correction
    real(wp)            :: dcoeff(5*2)  ! delta of coefficients for t2m bias correction
    type(t_coeff), pointer :: p_coeff   ! pointer to coefficients for t2m and rh2m bias correction

    if (biascor_mode == 0) return

    !-----------------
    ! set sums to zero
    !-----------------
    pp => synop_bc% data(:,:)
    pp% n     = 0
    pp% ob    = 0._wp  ! obs-bg  mean deviation
    pp% ob_ob = 0._wp  ! obs-bg  variance
    pp% oa    = 0._wp  ! obs-ana mean deviation
    pp% oa_oa = 0._wp  ! obs-ana variance
    pp% o_err = 0._wp  ! nominal obs error
    pp% b_err = 0._wp  ! nominal bg  error

    np_synop  = 0      ! number of synop land processed
    nc_synop  = 0      ! number of synop land corrected
    np_ship   = 0      ! number of synop ship processed
    nc_ship   = 0      ! number of synop ship corrected
    np_buoy   = 0      ! number of dribu      processed
    nc_buoy   = 0      ! number of dribu      processed
    np_metar  = 0      ! number of metar      processed
    nc_metar  = 0      ! number of metar      processed
    !-----------------------
    ! Update bias statistics
    !-----------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype /= OT_SYNOP .and. &
            obs% o(ib)% spot(is)% hd% obstype /= OT_DRIBU       ) cycle
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                      cycle
        first = .true.
        count = .true.

        !-----------------------------------------------------------------
        ! Cloudcover information for T2M and RH2M bias correction
        !-----------------------------------------------------------------
        clcover = obs% o(ib)% spot(is)% phase
        if (clcover < 0 .or. clcover > 8) clcover = 8

        !--------------------------------------------------------------------------
        ! Construct basis functions for T2M and RH2M bias correction
        ! (5 trigonometric for diurnal variations and 2 polynomial for cloud cover)
        !--------------------------------------------------------------------------
        ! Array M_trig with trigonometric basis functions.
        rhour = ihh (obs% o(ib)% spot(is)% hd% time)
        M_trig(1) = 1
        do k = 1, 2
           M_trig(2*k)   = sin(k*scal*rhour)      ! scaling with scal = 2*pi/24
           M_trig(2*k+1) = cos(k*scal*rhour)
        end do
        ! extended by polynomial basis functions for cloud cover
        M_bfct = (/ M_trig, (9-clcover) * M_trig /)

        do i = obs% o(ib)% spot(is)% o% i + 1,                       &
               obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
          if (obs% o(ib)% body(i)% use% state <= STAT_DISMISS) cycle
          select case (obs% o(ib)% varno(i))
          case (VN_Z)
             scali = 1._wp
             j = OB_PS
          case (VN_PS)
             scali = 1._wp / obs% o(ib)% spot(is)% pz_bg
             j = OB_PS
          case (VN_T2M)
             scali = 1._wp
             j = OB_T2M
          case (VN_RH2M)
             scali = 1._wp
             j = OB_RH2M
          case default
             cycle
          end select

          select case (obs% o(ib)% body(i)% use% state)
          case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
            if (first .OR. j == OB_T2M .OR. j == OB_RH2M) then  ! only once for Z/PS!
               bg  = obs% o(ib)% body(i)% bg
               an  = y% s(ib)% x(i)
               o   = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
               dob = (o - bg)                * scali    ! convert to height
               doa = (o - an)                * scali
               eo  = obs% o(ib)% body(i)% eo * scali
               eb  = obs% o(ib)% body(i)% eb * scali
               p   => synop_bc% data(j,ibc)
               p% n      = p% n      + 1
               p% ob     = p% ob     + dob
               p% ob_ob  = p% ob_ob  + dob * dob
               p% oa     = p% oa     + doa
               p% oa_oa  = p% oa_oa  + doa * doa
               p% o_err  = p% o_err  + eo
               p% b_err  = p% b_err  + eb
               if (j == OB_PS) first = .false.
            end if
            ! In case of t2m or rh2m update coefficients for bias statistic via 3dvar
            ! dcoeff = H_t (alpha + H_t^t H_t)^-1  dbias
            ! where H_t = M_bfct (trig. basis functions for time, polynomial for cloudcover)
            ! dbias = raw observation minus first guess minus bias
            if (j == OB_T2M .and. bc_t2m) then
               dbias = o - bg + obs% o(ib)% body(i)% bc  ! dob = (o - bg)
               dcoeff  = M_bfct * 1/(reg_param + dot_product( M_bfct, M_bfct )) * dbias
               p_coeff => synop_bc% coeff(ibc)
               p_coeff% t2m_bc_a  = p_coeff% t2m_bc_b  + dcoeff      ! update coefficients
            end if
            if (j == OB_RH2M .and. bc_rh2m) then
               dbias   = o - bg + obs% o(ib)% body(i)% bc  ! dob = (o - bg)
               dcoeff  = M_bfct * 1/(reg_param + dot_product( M_bfct, M_bfct )) * dbias
               p_coeff => synop_bc% coeff(ibc)
               p_coeff% rh2m_bc_a  = p_coeff% rh2m_bc_b + dcoeff     ! update coefficients
            end if
          end select
          if (count) then
            cor = (obs% o(ib)% body(i)% bc /= 0._sp)
            select case (obs% o(ib)% spot(is)% hd% codetype)
            case (11:14,20)
               np_synop          = np_synop + 1
               if (cor) nc_synop = nc_synop + 1
            case (21:24)
               np_ship           = np_ship  + 1
               if (cor) nc_ship  = nc_ship  + 1
            case (140,17)
               np_metar          = np_metar + 1
               if (cor) nc_metar = nc_metar + 1
            case (165)
               np_buoy           = np_buoy  + 1
               if (cor) nc_buoy  = nc_buoy  + 1
            end select
            count = .false.      ! Count each report only once
          end if
        end do ! i
      end do   ! is
    end do     ! ib

    np_synop = p_sum (np_synop)
    nc_synop = p_sum (nc_synop)
    np_ship  = p_sum (np_ship)
    nc_ship  = p_sum (nc_ship)
    np_buoy  = p_sum (np_buoy)
    nc_buoy  = p_sum (nc_buoy)
    np_metar = p_sum (np_metar)
    nc_metar = p_sum (nc_metar)
    !----------------
    ! sum up over PEs
    !----------------
    synop_bc% data% n      = p_sum (synop_bc% data% n)
    synop_bc% data% ob     = p_sum (synop_bc% data% ob)
    synop_bc% data% ob_ob  = p_sum (synop_bc% data% ob_ob)
    synop_bc% data% oa     = p_sum (synop_bc% data% oa)
    synop_bc% data% oa_oa  = p_sum (synop_bc% data% oa_oa)
    synop_bc% data% o_err  = p_sum (synop_bc% data% o_err)
    synop_bc% data% b_err  = p_sum (synop_bc% data% b_err)

    !-----------------------------------------
    ! update valid station height, coordinates
    ! to "last used" or "last known" position
    !-----------------------------------------
    ns        = synop_bc% n
    allocate (tmpz(ns), tmplat(ns), tmplon(ns))
    tmpz  (:) = synop_bc% plat% z
    tmplat(:) = synop_bc% plat% lat
    tmplon(:) = synop_bc% plat% lon

    if (dace% lpio) then
       tmpz   = -999._sp
       tmplat = -999._sp
       tmplon = -999._sp
    end if
    tmpz      = p_max (tmpz)
    tmplat    = p_max (tmplat)
    tmplon    = p_max (tmplon)
    if (dace% lpio) then
       where (tmpz   /= -999._sp) synop_bc% plat% z   = tmpz
       where (tmplat /= -999._sp) synop_bc% plat% lat = tmplat
       where (tmplon /= -999._sp) synop_bc% plat% lon = tmplon
    end if
    deallocate (tmpz, tmplat, tmplon)

    !-----------------------------------------------
    ! Remove statistics which is not to be collected
    !-----------------------------------------------
!   if (.not. (bc_t2m .or. bc_rh2m)) then
!     pp(2:3,:) = t_bc()
!   end if

    !-----------------------------------------
    ! weight factor for accumulated statistics
    !-----------------------------------------
    dt = days (ana_time - time_cyyyymmddhhmm (synop_bc% h% last_date))
    if (synop_bc% h% t_decay(1) > 0._wp .and. dt > 0._wp) then
      f = exp ( - dt / synop_bc% h% t_decay(1))
    else
      f = 1._wp
    endif

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => synop_bc% data(:,:)
    !-------------------------------
    ! rescale accumulated statistics
    !-------------------------------
    pp% n_ac       = pp% n_ac     * f
    pp% ob_ac      = pp% ob_ac    * pp% n_ac
    pp% ob_ob_ac   = pp% ob_ob_ac * pp% n_ac
    pp% oa_ac      = pp% oa_ac    * pp% n_ac
    pp% oa_oa_ac   = pp% oa_oa_ac * pp% n_ac
    pp% o_err_ac   = pp% o_err_ac * pp% n_ac
    pp% b_err_ac   = pp% b_err_ac * pp% n_ac
    !-----------
    ! accumulate
    !-----------
    pp% n_ac       = pp% n_ac     + pp% n
    pp% ob_ac      = pp% ob_ac    + pp% ob
    pp% ob_ob_ac   = pp% ob_ob_ac + pp% ob_ob
    pp% oa_ac      = pp% oa_ac    + pp% oa
    pp% oa_oa_ac   = pp% oa_oa_ac + pp% oa_oa
    pp% o_err_ac   = pp% o_err_ac + pp% o_err
    pp% b_err_ac   = pp% b_err_ac + pp% b_err
    !--------
    ! rescale
    !--------
    where (pp% n > 0)
      pp% ob       = pp% ob       / pp% n
      pp% ob_ob    = pp% ob_ob    / pp% n
      pp% oa       = pp% oa       / pp% n
      pp% oa_oa    = pp% oa_oa    / pp% n
      pp% o_err    = pp% o_err    / pp% n
      pp% b_err    = pp% b_err    / pp% n
    endwhere

    where (pp% n_ac > 0)
      pp% ob_ac    = pp% ob_ac    / pp% n_ac
      pp% ob_ob_ac = pp% ob_ob_ac / pp% n_ac
      pp% oa_ac    = pp% oa_ac    / pp% n_ac
      pp% oa_oa_ac = pp% oa_oa_ac / pp% n_ac
      pp% o_err_ac = pp% o_err_ac / pp% n_ac
      pp% b_err_ac = pp% b_err_ac / pp% n_ac
    endwhere

    !---------------------------
    ! Uncompress aggregated data
    !---------------------------
    call uncompress (synop_bc)

    !-----
    ! sort
    !-----
    allocate (iperm (synop_bc% n))
    call sort (synop_bc% plat% statid, iperm, 1, ier)
    synop_bc% plat = synop_bc% plat (  iperm)
    synop_bc% ag   = synop_bc% ag   (  iperm)
    synop_bc% data = synop_bc% data (:,iperm)
    if (bc_t2m .or. bc_rh2m) synop_bc% coeff = synop_bc% coeff (iperm)
    deallocate (iperm)

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,'(a)') '  SYNOP height/bias correction'
      write(6,'()')
      write(6,'(a)') '    SYNOP bias correction statistics'
      write(6,'()')
      write(6,42)    'SYNOP LAND', nc_synop, np_synop, nc_synop*100./max(np_synop,1)
      write(6,42)    'SYNOP SHIP', nc_ship , np_ship , nc_ship *100./max(np_ship ,1)
      write(6,42)    'DRIBU     ', nc_buoy , np_buoy , nc_buoy *100./max(np_buoy ,1)
      write(6,42)    'METAR     ', nc_metar, np_metar, nc_metar*100./max(np_metar,1)
42    format (4x,a,": ",i8,2x,"out of",i8,2x,"(",f5.1," %)")
    endif

    call remove_inactive (synop_bc)
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'()')
      write(6,'(a,a)')    '    last_date = ',synop_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    synop_bc% h%    t_decay(1)
      write(6,'(a,f13.4)')'    dt        =',                    dt
      write(6,'(a,f13.4)')'    f         =',                    f
      write(6,'(a,i8)')   '    n         =',sum(synop_bc% data% n)
      write(6,'(a,f11.2)')'    n_ac      =',sum(synop_bc% data% n_ac)
      write(6,'()')
      if (biascor_mode >= BC_FG) call write_bias_diag (synop_bc)

      if (verbose > 0) then
        !---------------------------------
        ! Verbosity of bias correction:
        ! 1 = print only non-empty entries
        ! 2 = print all
        !---------------------------------
        ! PS statistics
        write(6,'(a)') '    Statistics for PS:'
        write(6,'(a)') '    statid   codetype  n_ac     bias   stddev'
        write(6,'()')
        do i = 1, synop_bc% n
          if (synop_bc% data(1,i)% n_ac > 0._wp)                &
            write(6,'(4x,a,1x,i5,3(2x,f7.2))')                  &
                                synop_bc% plat(i)% statid,      &
                                synop_bc% plat(i)% code,        &
                                synop_bc% data(1,i)% n_ac,      &
                                synop_bc% data(1,i)% bc_a,      &
              sqrt (max (0._wp, synop_bc% data(1,i)% ob_ob_ac   &
                              - synop_bc% data(1,i)% ob_ac ** 2))
        end do
        write(6,'()')
        ! T2M statistics
        write(6,'(a)') '    Statistics for T2M:'
        write(6,'(a)') '    statid   codetype  n_ac     bias   stddev'
        write(6,'()')
        do i = 1, synop_bc% n
          if (synop_bc% data(2,i)% n_ac > 0._wp)                &
            write(6,'(4x,a,1x,i5,3(2x,f7.2))')                  &
                                synop_bc% plat(i)% statid,      &
                                synop_bc% plat(i)% code,        &
                                synop_bc% data(2,i)% n_ac,      &
                                synop_bc% data(2,i)% bc_a,      &
              sqrt (max (0._wp, synop_bc% data(2,i)% ob_ob_ac   &
                              - synop_bc% data(2,i)% ob_ac ** 2))
        end do
        write(6,'()')
        ! T2M statistics
        write(6,'(a)') '    Statistics for RH2M:'
        write(6,'(a)') '    statid   codetype  n_ac     bias   stddev'
        write(6,'()')
        do i = 1, synop_bc% n
          if (synop_bc% data(3,i)% n_ac > 0._wp)                &
            write(6,'(4x,a,1x,i5,3(2x,f7.2))')                  &
                                synop_bc% plat(i)% statid,      &
                                synop_bc% plat(i)% code,        &
                                synop_bc% data(3,i)% n_ac,      &
                                synop_bc% data(3,i)% bc_a,      &
              sqrt (max (0._wp, synop_bc% data(3,i)% ob_ob_ac   &
                              - synop_bc% data(3,i)% ob_ac ** 2))
        end do
        write(6,'()')
      end if
    endif

    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_synop_bc_file (synop_bc)
    call destruct_bc_file (synop_bc)

  contains

    subroutine write_bias_diag (bc)
      type(t_synop_bc) ,intent(in) :: bc
      !------------------------------------------------------------------
      ! Write SYNOP bias diagnostics aggregated to a 60°x60° regular grid
      !------------------------------------------------------------------
      integer  :: i, j, k
      integer  :: n_s(6,3), n_b(6,3)
      real(wp) :: b_s(6,3), b_b(6,3)
      n_s = 0     ; n_b = 0
      b_s = 0._wp ; b_b = 0._wp
      do k = 1, bc% n
         if (bc% data(1,k)% n == 0) cycle
         i = int ((bc% plat(k)% lon + 180._sp) / 60._sp) + 1
         j = int ((bc% plat(k)% lat +  90._sp) / 60._sp) + 1
         if (i < 1 .or. i > 6) cycle
         if (j < 1 .or. j > 3) cycle
         select case (bc% plat(k)% code)
         case (11)
            n_s(i,j) = n_s(i,j) + 1
            b_s(i,j) = b_s(i,j) + bc% data(1,k)% ob
         case (165)
            n_b(i,j) = n_b(i,j) + 1
            b_b(i,j) = b_b(i,j) + bc% data(1,k)% ob
         end select
      end do
      !----------
      ! Normalize
      !----------
      where (n_s > 0) b_s = b_s / n_s
      where (n_b > 0) b_b = b_b / n_b
      if (sum (n_s) + sum (n_b) > 0) then
         write(6,'(A)')           '    Summary of regionally averaged bias &
                                  & [gpm]  (60°x60° grid)        obs.count:'
         write(6,'()')
         write(6,'(11x,7(i4,"°",:,"..."))')         (-180+60*i,i=0,6)
         write(6,'()')
         write(6,'(A,6f8.2,1x,6i5)') '  SYNOP NH:', b_s(:,3), n_s(:,3)
         write(6,'(A,6f8.2,1x,6i5)') '  SYNOP TR:', b_s(:,2), n_s(:,2)
         write(6,'(A,6f8.2,1x,6i5)') '  SYNOP SH:', b_s(:,1), n_s(:,1)
         write(6,'()')
         write(6,'(A,6f8.2,1x,6i5)') '  DRIBU NH:', b_b(:,3), n_b(:,3)
         write(6,'(A,6f8.2,1x,6i5)') '  DRIBU TR:', b_b(:,2), n_b(:,2)
         write(6,'(A,6f8.2,1x,6i5)') '  DRIBU SH:', b_b(:,1), n_b(:,1)
         write(6,'()')
      end if
    end subroutine write_bias_diag

  end subroutine synop_bc_aan

!------------------------------------------------------------------------------

  subroutine remove_inactive (bc)
    type(t_synop_bc) ,intent(inout) :: bc
    !----------------------------------------------------
    ! Remove inactive stations when accumulated effective
    ! number of observations falls below threshold
    ! *and* last aggregation is older than 7 days.
    !----------------------------------------------------
    integer             :: idx(bc% n)       ! Index list
    integer             :: i                ! Loop index
    integer             :: n                ! Entry count
    integer             :: ana_hour         ! Time of analysis
    integer             :: min_hour         ! "Age" limit for aggregated values
    type(t_synop_bc)    :: tmp
    real(wp), parameter :: thresh = 0.5_wp

    if (bc% n <= 0) return

    ana_hour = int (hours (ana_time))
    min_hour = ana_hour - 7 * 24

    idx = 0
    n   = 0
    do i = 1, bc% n
       if (        bc% ag  (  i)% hours >= min_hour .or. &
           maxval( bc% data(:,i)% n_ac) >  thresh        ) then
          n      = n + 1
          idx(n) = i
       end if
    end do
    if (n == bc% n) return

    if (dace% lpio) then
       write(6,'(/,A,i0)') '    inactive entries removed: ', bc% n - n
       if (verbose > 1) then
          write(6,'(/,a)') '    statid   codetype  n_ac(PS)  n_ac(T2M)  n_ac(RH2M)  age(agg)'
          do i = 1, bc% n
             if (        bc% ag    (i)% hours <  min_hour .and. &
                 maxval( bc% data(:,i)% n_ac) <= thresh         ) then
                write(6,'(4x,a,1x,i5,1x,f8.2,2x,f8.2,3x,f8.2,8x,i5)')  &
                                bc% plat(  i)% statid, &
                                bc% plat(  i)% code,   &
                                bc% data(1,i)% n_ac,   &
                                bc% data(2,i)% n_ac,   &
                                bc% data(3,i)% n_ac,   &
                     ana_hour - bc% ag  (  i)% hours
             end if
          end do
       end if
    end if

    allocate (tmp% plat(    n))
    allocate (tmp% ag  (    n))
    allocate (tmp% data(nob,n))
    tmp% plat = bc% plat(  idx(1:n))
    tmp% ag   = bc% ag  (  idx(1:n))
    tmp% data = bc% data(:,idx(1:n))
    deallocate (bc% plat, bc% data, bc% ag)
    bc% n     =  n
    bc% plat  => tmp% plat
    bc% ag    => tmp% ag
    bc% data  => tmp% data

    if (bc_rh2m .or. bc_t2m) then
       allocate    (tmp% coeff(n))
       tmp% coeff =  bc% coeff(idx(1:n))
       deallocate   (bc% coeff)
       bc% coeff => tmp% coeff
    endif

  end subroutine remove_inactive

!------------------------------------------------------------------------------

  subroutine compress (bc)
    type(t_synop_bc), intent(inout) :: bc
    !--------------------------------------
    ! Internally "compress" aggregated data
    !--------------------------------------

    integer          :: pm                ! no. processors for "compressed" data
    integer          :: lb(0:dace% npe-1) ! Lower bounds
    integer          :: ub(0:dace% npe-1) ! Upper bounds
    type(t_synop_bc) :: tmp               ! Temporary

    if (.not. compress_ag) return
    if (bc% compressed) return
    bc% compressed = .true.

    call get_compress_params (bc% n, pm, lb, ub)

    if (associated (bc% ag)) then
       allocate (tmp% ag(lb(dace% pe):ub(dace% pe)))
       tmp% ag =  bc% ag(lb(dace% pe):ub(dace% pe))
       deallocate (bc% ag)
       bc% ag => tmp% ag
    end if

  end subroutine compress
!------------------------------------------------------------------------------
  subroutine uncompress (bc)
    type(t_synop_bc), intent(inout) :: bc
    !--------------------------
    ! "Recover" aggregated data
    !--------------------------
    integer          :: pm                ! no. processors for "compressed" data
    integer          :: lb(0:dace% npe-1) ! Lower bounds
    integer          :: ub(0:dace% npe-1) ! Upper bounds
    type(t_synop_bc) :: tmp               ! Temporary

    if (.not. bc% compressed) return
    bc% compressed = .false.

    call get_compress_params (bc% n, pm, lb, ub)

    tmp = bc

    allocate         (bc% ag  (bc% n))
    call p_gather_ag (tmp% ag, bc% ag, dace% pio)
    call p_bcast_ag  (         bc% ag, dace% pio)
    deallocate       (tmp% ag)

  end subroutine uncompress
!------------------------------------------------------------------------------
  subroutine get_compress_params (n, pm, lb, ub)
    integer, intent(in)  :: n         ! No. elements to "compress"
    integer, intent(out) :: pm        ! no. processors keeping "compressed" data
    integer, intent(out) :: lb(0:)    ! Lower bounds
    integer, intent(out) :: ub(0:)    ! Upper bounds

    integer  :: pe            ! Processor index
    integer  :: off           ! Offset
    integer  :: nt            ! Chunk size

    pm  = min (dace% npe, n)
    off = 0
    do pe = 0, pm-1
       nt     = (n - (pe+1)) / dace% npe + 1
       lb(pe) = off + 1
       ub(pe) = off + nt
       off    = ub(pe)
    end do
    ub(pm:) = off
    lb(pm:) = off + 1
    if (off /= n) then
       if (dace% lpio) then
          write(0,*) "n :", n
          write(0,*) "lb:", lb
          write(0,*) "ub:", lb
       end if
       call p_barrier ()
       call finish ("get_compress_params","internal bounds checking error")
    end if
  end subroutine get_compress_params

!==============================================================================

! subroutine synop_bc_afg (obs, lwrite)
! !-----------------------------------------------------------------------
! ! Bias correction routine to be called after analysis
! ! Update bias correction statistics.
! ! Write updated correction coefficient file.
! !-----------------------------------------------------------------------
! type (t_obs_set) ,intent(inout) :: obs    ! observation
! logical          ,intent(in)    :: lwrite ! write bc file
!
! end subroutine synop_bc_afg

!==============================================================================
!----------------------------------------------------------------------------
! subroutine p_gather_tmp (sendbuffer,receivebuffer,root,[comm],[recvcounts])
!----------------------------------------------------------------------------
#define DERIVED type(t_tmp)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_tmp
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  p_gather_DERIVED
!==============================================================================
!--------------------------------------------------------------------
! subroutine p_scatter_tmp (sendbuf,recvbuf,root,[comm],[sendcounts])
!--------------------------------------------------------------------
#define DERIVED type(t_tmp)
#undef  MPI_TYPE
#define p_scatter_DERIVED p_scatter_tmp
#include "p_scatter_derived.incf"
#undef  DERIVED
!==============================================================================
!----------------------------------------------------
! subroutine p_bcast_plat (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_plat_bc),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_plat
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  p_bcast_DERIVED
!==============================================================================
!----------------------------------------------------
! subroutine p_bcast_ag (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_ag),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_ag
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  p_bcast_DERIVED
!==============================================================================
!---------------------------------------------------------------------------
! subroutine p_gather_ag (sendbuffer,receivebuffer,root,[comm],[recvcounts])
!---------------------------------------------------------------------------
#define DERIVED type(t_ag)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_ag
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  p_gather_DERIVED
!==============================================================================
!----------------------------------------------------
! subroutine p_bcast_data (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_bc),dimension(:,:)
#define VECTOR
#define RANK 2
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_data
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  RANK
#undef  p_bcast_DERIVED
!==============================================================================
!==============================================================================
!----------------------------------------------------
! subroutine p_bcast_coeff (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_coeff),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_coeff
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  p_bcast_DERIVED
!==============================================================================
end module mo_synop_bc
