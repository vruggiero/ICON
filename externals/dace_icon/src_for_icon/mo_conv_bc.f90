!
!+ Bias correction for individual conventional platforms
!
MODULE mo_conv_bc
!
! Description:
!   Bias correction for individual conventional platforms.
!   Currently: Aircraft (temperature) bias correction
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_13        2011/11/01 Andreas Rhodin
!  Bias correction for individual conventional platforms
! V1_19        2012-04-16 Andreas Rhodin
!  options to read instantaneous statistics as well
!              use data even if biascorrection is not done
! V1_20        2012-06-18 Harald Anlauf
!  verbosity level for bias correction
! V1_22        2013-02-13 Harald Anlauf
!  stubs for aircraft flight-track reconstruction
! V1_23        2013-03-26 Harald Anlauf
!  Increase fine-graining of diagnostic output for aircraft bias correction
! V1_26        2013/06/27 Andreas Rhodin
!  replace STAT_ACTIVE_0 by STAT_NOTACTIVE
! V1_27        2013-11-08 Andreas Rhodin
!  t_obs: provide components for bias correction predictors
! V1_28        2014/02/26 Harald Anlauf
!  aircraft_bc_init: ignore existing bias corr. when flag_biasc_airep==0
! V1_29        2014/04/02 Robin Faulwetter
!  Do not write bias_AIREP*, if no AIREP data and statistics are available
! V1_30        2014-04-09 Andreas Rhodin
!  make recalculation of aircraft phase independent from PE configuration
! V1_48        2016-10-06 Robin Faulwetter
!  Implemented level-based thinning
! V1_49        2016-10-25 Harald Anlauf
!  t_spot: generalize bc_airep -> bc_index
! V1_50        2017-01-09 Harald Anlauf
!  aircraft_bc_aan: fix typo in accumulation
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-----------------------------------------------------------

!=============
! Modules used
!=============

  use mo_kind,       only: wp,              &! working precision kind parameter
                           i8                ! 8 byte integer kind parameter
  use mo_exception,  only: finish            ! abort in case of error
  use mo_time,       only: days,            &! derive days (real)
                           seconds,         &! derive seconds (real)
                           operator(-),     &! calculate time difference
                           cyyyymmddhhmm,   &! string from time
                           time_cyyyymmddhhmm! time from string
! use mo_t_obs,      only: t_obs             ! observation data type
  use mo_t_obs,      only: t_spot            ! observation data type
  use mo_obs_set,    only: t_obs_set         ! observation data derived type
  use mo_obs_tables, only: rept_use          ! use table entry
  use mo_dec_matrix, only: t_vector          ! vector data type
  use mo_mpi_dace,   only: dace,            &! MPI group info
                           p_sum,           &! sum over PEs
                           p_bcast,         &! generic MPI broadcast routine
                           p_gather          ! generic MPI gather routine
  use mo_run_params, only: ana_time,        &! analysis time
!                          run_type,        &! haupt=0, vor=1, ass=2
                           flag_biasc_airep  ! aircraft bias correction
  use mo_biasc_io,   only: t_bcor_head,     &! file header derived type
                           new_bc_head,     &! construct new file header
                           open_bc_read,    &! open  file for reading
                           open_bc_write,   &! open  file for writing
                           close_bc,        &! close file
                           bc_paths,        &! full pathnames
                           bc_obstyp,       &! observation type in files
                           nbcf,            &! number of files in list
                           verbose           ! Verbosity level of bias corr.
  use mo_fdbk_tables,only: VN_T,            &! temperature observation flag
                           VN_RH,           &! relative humidity       flag
                           OT_AIREP          ! Aircraft report type ID
  use mo_t_use,      only: CHK_BIASCOR,     &! flag for no bias correction
                           decr_use,        &! decrease the state of a datum
                           STAT_PASSIVE,    &! observation status flags
                           STAT_ACTIVE_0I,  &!
                           STAT_ACCEPTED     !
  use netcdf,        only:                  &! NetCDF f90 interface
                           nf90_def_dim,    &! define dimensions
                           nf90_def_var,    &! define variables
                           nf90_put_att,    &! define attribute
                           nf90_enddef,     &! end of definition mode
                           nf90_put_var,    &! write variable
                           nf90_strerror,   &! derive error character string
                           NF90_FLOAT,      &! float     type id
                           NF90_CHAR,       &! character type id
                           NF90_GLOBAL,     &! global attribute flag
                           NF90_NOERR,      &! NetCDF return code for no error
                           NF90_FILL_FLOAT   ! NetCDF fillvalue
  use mo_t_netcdf,   only: ncid,            &! NetCDF file id to use
                           stanc,           &! NetCDF start  parameter to use
                           counc,           &! NetCDF count  parameter to use
                           strnc,           &! NetCDF stride parameter to use
                           get_var,         &! get variable
                           get_dim           ! get dimension
  use slatec_module, only: sort              ! sort routine
  implicit none

!================
! Public entities
!================

  private
  !--------------
  ! derived types
  !--------------
  public :: t_aircraft_bc    ! aircraft bias correction data
  public :: t_plat_bc, t_bc  ! components of t_aircraft_bc
  !------------
  ! subroutines
  !------------
  public :: aircraft_bc_init ! initialize module: read namelist, biascor.coeff.
  public :: aircraft_bc_bfg  ! apply bias correction, called before fg-check
  public :: aircraft_bc_aan  ! update bias correction coefs., after analysis
  public :: read_bcor_file   ! read aircraft bias correction file
  !----------------------------------------------------
  ! namelist parameters to be set in namelist AIREP_OBS
  !----------------------------------------------------
  public :: biascor_mode     ! mode used for updating
  public :: t_decay          ! accumulation decay time (days)
  public :: force_level      ! force phase to level above this pressure
  public :: fr_land_bc       ! land fraction demanded for statistics
  public :: n_required       ! number of entries required for correction
  public :: bc_fallback      ! fallback if biasc-file not present
  public :: scan_tracks      ! scan tracks, set (ascend/descend/flight) phase
  public :: time_sep         ! time required for separated tracks
  public :: rate_asc         ! ascend  rate required
  public :: rate_desc        ! descend rate required
  public :: rate_lev         ! level   rate required
  public :: BC_NOBC,BC_UP,BC_FG,BC_AN,BC_VARBC ! values for biascor_mode

!=========================
! Derived type definitions
!=========================

  !---------------------------------------------------------
  ! BUFR MPHAI (8004) : PHASE OF AIRCRAFT FLIGHT code figure
  !---------------------------------------------------------
  integer,parameter :: BUFR_UNSTEADY = 2
  integer,parameter :: BUFR_LEVEL    = 3
  integer,parameter :: BUFR_LEVEL_HW = 4
  integer,parameter :: BUFR_ASCEND   = 5
  integer,parameter :: BUFR_DESCEND  = 6
  integer,parameter :: BUFR_MISSING  = 7

  !-----------------------------
  ! aircraft phase specific data
  !-----------------------------
  type t_bc
    integer  :: n         = 0      ! number of entries
    real(wp) :: ob        = 0._wp  ! obs-bg  mean deviation
    real(wp) :: oa        = 0._wp  ! obs-ana mean deviation
    real(wp) :: ab        = 0._wp  ! ana-obs mean deviation
    real(wp) :: ob_ob     = 0._wp  ! obs-bg  variance
    real(wp) :: oa_oa     = 0._wp  ! obs-ana variance
    real(wp) :: ab_ab     = 0._wp  ! ana-obs variance
    real(wp) :: ob_oa     = 0._wp  !
    real(wp) :: ob_ab     = 0._wp  !
    real(wp) :: oa_ab     = 0._wp  !
    real(wp) :: n_ac      = 0._wp  ! number of accumulated entries
    real(wp) :: ob_ac     = 0._wp  ! obs-bg  mean     accumulated
    real(wp) :: oa_ac     = 0._wp  ! ana-bg  mean     accumulated
    real(wp) :: ab_ac     = 0._wp  ! ana-obs mean     accumulated
    real(wp) :: ob_ob_ac  = 0._wp  ! obs-bg  variance accumulated
    real(wp) :: oa_oa_ac  = 0._wp  ! ana-bg  variance accumulated
    real(wp) :: ab_ab_ac  = 0._wp  ! ana-obs variance accumulated
    real(wp) :: ob_oa_ac  = 0._wp  !
    real(wp) :: ob_ab_ac  = 0._wp  !
    real(wp) :: oa_ab_ac  = 0._wp  !
    real(wp) :: o_err     = 0._wp  ! nominal obs error
    real(wp) :: b_err     = 0._wp  ! nominal bg  error
    real(wp) :: o_err_ac  = 0._wp  ! nominal obs error
    real(wp) :: b_err_ac  = 0._wp  ! nominal bg  error
    real(wp) :: bc_b      = 0._wp  ! bias correction background value
    real(wp) :: bc_a      = 0._wp  ! bias correction analysis   value
!   real(wp) :: bc_b_err           ! bias correction background error
!   real(wp) :: bc_a_err           ! bias correction analysis   error
  end type t_bc

  !------------------------------
  ! aircraft flight phase indices
  !------------------------------
  integer, parameter :: nph      = 4  ! number of aircraft phases
  integer ,parameter :: PH_ASC   = 1  ! ascending
  integer ,parameter :: PH_DES   = 2  ! descending
  integer ,parameter :: PH_LEV   = 3  ! flight level
  integer ,parameter :: PH_UNSP  = 4  ! un-specified

  !--------------------------
  ! observed quantity indices
  !--------------------------
  integer, parameter :: nob   = 2     ! number of observed quantities
  integer ,parameter :: OB_T  = 1     ! temperature
  integer ,parameter :: OB_RH = 2     ! relative humidity

  !-----------------------
  ! platform specific data
  !-----------------------
  type t_plat_bc
    character(len=8) :: statid = ''   ! aircraft registration number
    integer  (i8)    :: hash   = 0    ! integer representation of name
  end type t_plat_bc

  !------------------------------
  ! Aircraft bias correction data
  !------------------------------
  type t_aircraft_bc
    type(t_bcor_head)        :: h                 ! file header data
    integer                  :: biascor_mode = 0  ! mode used for updating
    integer                  :: n            = 0  ! number of entries
    type (t_plat_bc),pointer :: plat (:)          ! platform specific metadata
    type (t_bc)     ,pointer :: data (:,:,:)      ! (nph,nob,:) data
  end type t_aircraft_bc

  !-------------------------------------
  ! derived type used for cross-checking
  !-------------------------------------
  type t_tmp
    character(len=8) :: statid ! aircraft registration number
    integer  (i8)    :: hash   ! integer representation of name
    integer          :: pe     ! processor       of report
    integer          :: i      ! old index
    integer          :: j      ! new index
  end type t_tmp

  !--------------------------------------------
  ! derived type used for flight-track assembly
  !--------------------------------------------
  type t_track
!   type(t_plat_bc)  :: plat          ! platform specific metadata
    character(len=8) :: statid  = ''  ! aircraft registration number
    integer  (i8)    :: hash    = 0   ! integer representation of name
    integer          :: time_oa = 0   ! obs.time - ana.time    [seconds]
    integer          :: pe            ! processor       of report
    integer          :: ib            ! observation box  index
    integer          :: is            ! observation spot index
    integer          :: phase         ! flight phase
    real             :: level         ! height                      [Pa]
    integer          :: id      = 0   ! observation ID from source/subset
  end type t_track

  !------------------------
  ! values for biascor_mode
  !------------------------
  integer ,parameter :: BC_NOBC  = 0  ! no bias correction
  integer ,parameter :: BC_UP    = 1  ! only update bias corrrection file
  integer ,parameter :: BC_FG    = 2  ! apply  bias corr. from first guess
  integer ,parameter :: BC_AN    = 3  ! apply  bias corr. from analysis
  integer ,parameter :: BC_VARBC = 4  ! variational bias correction

!=================
! Module variables
!=================

  type(t_aircraft_bc),save :: aircraft_bc              ! bias correction data
  !-----------------
  ! namelist entries
  !-----------------
  integer                  :: biascor_mode = BC_NOBC   ! mode used for updating
  real(wp)                 :: t_decay      =   -30._wp ! accumulation decaytime
  real(wp)                 :: force_level  = 31000._wp ! force phase to level
  real(wp)                 :: fr_land_bc   =    0.6_wp ! land fraction
  real(wp)                 :: rate_asc     =     0._wp ! ascend  rate required
  real(wp)                 :: rate_desc    =     0._wp ! descend rate required
  real(wp)                 :: rate_lev     =     0._wp ! level   rate required
  integer                  :: time_sep     =   900     ! time separation (s)
  integer                  :: n_required   =    50     ! # of entries required
  logical                  :: bc_fallback  = .false.   ! biasc-file not present
  logical                  :: scan_tracks  = .true.    ! scan tracks, set phase
  !----------------------------------
  ! Aircraft flight track information
  !----------------------------------
  type(t_track), pointer :: tracks(:) => NULL ()  ! Aircraft flight tracks
  integer                :: ntracks   =  0        ! No. of spots

contains
!==============================================================================

  subroutine aircraft_bc_init
  !---------------------------------------
  ! Initialize this module:
  ! read bias-correction coefficient files
  !---------------------------------------

    integer :: i    ! loop index
    integer :: ierr ! error return value

    if (biascor_mode /= 0) then
      !---------------------------------------
      ! read bias-correction coefficient files
      !---------------------------------------
      ierr = -1
      if (flag_biasc_airep /= 0 .or. .not. bc_fallback) then
        do i = 1, nbcf
          if (bc_obstyp(i) == OT_AIREP) then
            call read_bcor_file (aircraft_bc, bc_paths(i))
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
          call finish ('aircraft_bc_init','coefficient file not present !')
        call new_bc_head (aircraft_bc% h, OT_AIREP, t_decay=(/abs(t_decay)/))
        aircraft_bc% biascor_mode = biascor_mode
        aircraft_bc% n            = 0
        allocate (aircraft_bc% plat (0))
        allocate (aircraft_bc% data (nph,nob,0))
        if (dace% lpio) write(6,'(a,a)')'    created empty file ', &
                                          trim (aircraft_bc%h% path)
      endif
      !-------------------------------
      ! prepare bias correction to use
      !-------------------------------
      select case (biascor_mode)
      case (BC_UP)
        aircraft_bc% data% bc_b = 0
        aircraft_bc% data% bc_a = 0
      case (BC_FG)
        where (aircraft_bc% data% n_ac >= n_required)
          aircraft_bc% data% bc_b = aircraft_bc% data% ob_ac
          aircraft_bc% data% bc_a = aircraft_bc% data% ob_ac
        elsewhere
          aircraft_bc% data% bc_b = NF90_FILL_FLOAT
          aircraft_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_AN)
        where (aircraft_bc% data% n_ac >= n_required)
          aircraft_bc% data% bc_b = aircraft_bc% data% oa_ac
          aircraft_bc% data% bc_a = aircraft_bc% data% oa_ac
        elsewhere
          aircraft_bc% data% bc_b = NF90_FILL_FLOAT
          aircraft_bc% data% bc_a = NF90_FILL_FLOAT
        endwhere
      case (BC_VARBC)
        aircraft_bc% data% bc_b = aircraft_bc% data% bc_a
      end select
    endif

  end subroutine aircraft_bc_init

!------------------------------------------------------------------------------

  subroutine read_bcor_file (bc, file, inst)
  !-----------------------------------
  ! read aircraft bias correction file
  !-----------------------------------
  type(t_aircraft_bc) ,intent(out) :: bc    ! derived type variable to fill
  character(len=*)    ,intent(in)  :: file  ! name of file to read
  logical   ,optional ,intent(in)  :: inst  ! read instantanious data as well

    integer :: i
    logical :: linst

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
      allocate (bc% plat (        bc% n))
      allocate (bc% data (nph,nob,bc% n))

      !---------------
      ! read variables
      !---------------
      if (bc% n > 0) then
        stanc = 1
        counc = 0
        strnc = 1
        call get_var (bc% plat% statid   ,'statid'  )
        call get_var (bc% data% n_ac     ,'n_ac'    )
        call get_var (bc% data% ob_ac    ,'ob_ac'   )
        call get_var (bc% data% oa_ac    ,'oa_ac'   )
        call get_var (bc% data% ab_ac    ,'ab_ac'   )
        call get_var (bc% data% ob_ob_ac ,'ob_ob_ac')
        call get_var (bc% data% oa_oa_ac ,'oa_oa_ac')
        call get_var (bc% data% ab_ab_ac ,'ab_ab_ac')
        call get_var (bc% data% ob_oa_ac ,'ob_oa_ac')
        call get_var (bc% data% ob_ab_ac ,'ob_ab_ac')
        call get_var (bc% data% oa_ab_ac ,'oa_ab_ac')
        call get_var (bc% data% o_err_ac ,'o_err_ac')
        call get_var (bc% data% b_err_ac ,'b_err_ac')
        call get_var (bc% data% bc_a     ,'bc_a'    )
        if (linst) then
          call get_var (bc% data% n     ,'n'    )
          call get_var (bc% data% ob    ,'ob'   )
          call get_var (bc% data% oa    ,'oa'   )
          call get_var (bc% data% ab    ,'ab'   )
          call get_var (bc% data% ob_ob ,'ob_ob')
          call get_var (bc% data% oa_oa ,'oa_oa')
          call get_var (bc% data% ab_ab ,'ab_ab')
          call get_var (bc% data% ob_oa ,'ob_oa')
          call get_var (bc% data% ob_ab ,'ob_ab')
          call get_var (bc% data% oa_ab ,'oa_ab')
          call get_var (bc% data% o_err ,'o_err')
          call get_var (bc% data% b_err ,'b_err')
        endif
!       bc% plat% hash = transfer (bc% plat% statid, bc% plat% hash)
        do i = 1, bc% n
          bc% plat(i)% hash = transfer (bc% plat(i)% statid, bc% plat(i)% hash)
        end do
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
      allocate (bc% plat (        bc% n))
      allocate (bc% data (nph,nob,bc% n))
    endif
    call p_bcast_plat (aircraft_bc% plat, dace% pio)
    call p_bcast_data (aircraft_bc% data, dace% pio)

  end subroutine read_bcor_file

!------------------------------------------------------------------------------

  subroutine remove_inactive (bc)
    type(t_aircraft_bc) ,intent(inout) :: bc
    !----------------------------------------------------
    ! Remove inactive stations when accumulated effective
    ! number of observations falls below threshold
    !----------------------------------------------------
    integer             :: idx(bc% n)       ! Index list
    integer             :: i, j, k
    integer             :: n
    type(t_aircraft_bc) :: tmp
    real(wp), parameter :: thresh = 0.5_wp

    if (bc% n <= 0) return

    idx = 0
    n   = 0
    do i = 1, bc% n
       if (sum (bc% data(:,:,i)% n_ac) > thresh) then
          n      = n + 1
          idx(n) = i
       end if
    end do
    if (n == bc% n) return

    if (dace% lpio) then
       write(6,'(/,A,i0)') '    inactive entries removed: ', bc% n - n
       if (verbose > 1) then
          write(6,'(/,a)') '    statid      n(T)   n(rh)'
          do i = 1, bc% n
             if (sum (bc% data(:,:,i)% n_ac) <= thresh) then
                write(6,'(4x,a,2f8.3)') aircraft_bc% plat    (i)% statid, &
                                   sum (aircraft_bc% data(:,1,i)% n_ac),  &
                                   sum (aircraft_bc% data(:,2,i)% n_ac)
             end if
          end do
       end if
    end if

    j = size (bc% data, dim=1)
    k = size (bc% data, dim=2)
    allocate (tmp% plat(    n))
    allocate (tmp% data(j,k,n))
    tmp% plat = bc% plat(    idx(1:n))
    tmp% data = bc% data(:,:,idx(1:n))
    deallocate (bc% plat, bc% data)
    bc% n     =  n
    bc% plat  => tmp% plat
    bc% data  => tmp% data

  end subroutine remove_inactive

!------------------------------------------------------------------------------

  subroutine write_bcor_file (bc)
  !------------------------------------
  ! write aircraft bias correction file
  !------------------------------------
  type(t_aircraft_bc) ,intent(inout) :: bc

    integer :: status                            ! NetCDF return value
    integer :: dimid (3)                         ! NetCDF dimension id
    integer :: dimch                             ! NetCDF dimension id
    integer :: statid                            ! NetCDF variable ids
    integer :: n                                 ! .
    integer :: ob, oa, ab                        ! .
    integer :: ob_ob, oa_oa, ab_ab               ! .
    integer :: ob_oa, ob_ab, oa_ab               ! .
    integer :: n_ac                              ! .
    integer :: ob_ac, oa_ac, ab_ac               ! .
    integer :: ob_ob_ac, oa_oa_ac, ab_ab_ac      ! .
    integer :: ob_oa_ac, ob_ab_ac, oa_ab_ac      ! .
    integer :: bc_b, b_err , b_err_ac ! bc_b_err ! .
    integer :: bc_a, o_err , o_err_ac ! bc_a_err ! .
    !++++++++++++++++++++++++++++++++++++++++++
    ! workaround for NEC SX9 nf90_put_var(text)
    !++++++++++++++++++++++++++++++++++++++++++
    character(len=8),allocatable :: statids(:)

    if (dace% lpio .and. bc%n > 0) then
      !-------------------------------
      ! open NetCDF file, write header
      !-------------------------------
      bc% h% t_decay(2:) = bc% h% t_decay(1) ! Write only one the first t_decay value
      call open_bc_write (bc% h, OT_AIREP)

      !---------------
      ! set attributes
      !---------------
      status = nf90_put_att (bc%h% ncid, NF90_GLOBAL, 'biascor_mode', &
                                                       biascor_mode   )
      !------------------
      ! define dimensions
      !------------------
      call def_dim ('phase',       nph,            dimid(1))
      call def_dim ('observation', nob,            dimid(2))
      call def_dim ('platform',    size(bc% plat), dimid(3))
      call def_dim ('chars',       8,              dimch   )
      !-----------------
      ! define variables
      !-----------------
      status = nf90_def_var (bc%h% ncid ,'statid' ,NF90_CHAR, &
                             (/dimch,dimid(3)/), statid)
      status = nf90_put_att (bc%h% ncid , statid, 'longname',&
                             'station id as character string')
      call def_var ('n'         ,n         ,'entries')
      call def_var ('ob'        ,ob        ,'obs-bg  mean')
      call def_var ('oa'        ,oa        ,'obs-ana mean')
      call def_var ('ab'        ,ab        ,'ana-bg  mean')
      call def_var ('ob_ob'     ,ob_ob     ,'obs-bg  variance')
      call def_var ('oa_oa'     ,oa_oa     ,'obs-ana variance')
      call def_var ('ab_ab'     ,ab_ab     ,'ana-bg  variance')
      call def_var ('ob_oa'     ,ob_oa     ,'obs-bg  obs-ana covariance')
      call def_var ('ob_ab'     ,ob_ab     ,'obs-bg  ana-bg  covariance')
      call def_var ('oa_ab'     ,oa_ab     ,'obs-ana ana-bg  covariance')
      call def_var ('n_ac'      ,n_ac      ,'entries accumulated')
      call def_var ('ob_ac'     ,ob_ac     ,'obs-bg  mean accumulated')
      call def_var ('oa_ac'     ,oa_ac     ,'obs-ana mean accumulated')
      call def_var ('ab_ac'     ,ab_ac     ,'ana-bg  mean accumulated')
      call def_var ('ob_ob_ac'  ,ob_ob_ac  ,'obs-bg  variance accumulated')
      call def_var ('oa_oa_ac'  ,oa_oa_ac  ,'obs-ana variance accumulated')
      call def_var ('ab_ab_ac'  ,ab_ab_ac  ,'ana-bg  variance accumulated')
      call def_var ('ob_oa_ac'  ,ob_oa_ac  ,'obs-bg  obs-ana covariance acc.')
      call def_var ('ob_ab_ac'  ,ob_ab_ac  ,'obs-bg  ana-bg  covariance acc.')
      call def_var ('oa_ab_ac'  ,oa_ab_ac  ,'obs-ana ana-bg  covariance acc.')
      call def_var ('bc_b'      ,bc_b      ,'first guess bias correction')
      call def_var ('bc_a'      ,bc_a      ,'analysis bias correction')
      call def_var ('o_err'     ,o_err     ,'nominal obs error squared')
      call def_var ('b_err'     ,b_err     ,'nominal bg  error squared')
      call def_var ('o_err_ac'  ,o_err_ac  ,               &
                    'nominal obs error squared accumulated')
      call def_var ('b_err_ac'  ,b_err_ac  ,               &
                    'nominal bg  error squared accumulated')
!     call def_var ('bc_b'      ,bc_b      ,'bias correction background value')
!     call def_var ('bc_a'      ,bc_a      ,'bias correction analysis   value')
!     call def_var ('bc_b_err ' ,bc_b_err  ,                  &
!                   'bias correction background error squared')
!     call def_var ('bc_a_err ' ,bc_a_err  ,                 &
!                   'bias corrrection analysis error squared')
      status = nf90_enddef  (bc%h% ncid)
      !----------------
      ! write variables
      !----------------
      allocate (statids (bc% n))
      statids = bc% plat% statid
!+++  status  = nf90_put_var (bc%h% ncid, statid, bc% plat% statid)
      status  = nf90_put_var (bc%h% ncid, statid, statids)
      deallocate (statids)
      call put_vari ('n'        ,n        ,bc% data% n)
      call put_var  ('ob'       ,ob       ,bc% data% ob)
      call put_var  ('oa'       ,oa       ,bc% data% oa)
      call put_var  ('ab'       ,ab       ,bc% data% ab)
      call put_var  ('ob_ob'    ,ob_ob    ,bc% data% ob_ob)
      call put_var  ('oa_oa'    ,oa_oa    ,bc% data% oa_oa)
      call put_var  ('ab_ab'    ,ab_ab    ,bc% data% ab_ab)
      call put_var  ('ob_oa'    ,ob_oa    ,bc% data% ob_oa)
      call put_var  ('ob_ab'    ,ob_ab    ,bc% data% ob_ab)
      call put_var  ('oa_ab'    ,oa_ab    ,bc% data% oa_ab)
      call put_var  ('n_ac'     ,n_ac     ,bc% data% n_ac)
      call put_var  ('ob_ac'    ,ob_ac    ,bc% data% ob_ac)
      call put_var  ('oa_ac'    ,oa_ac    ,bc% data% oa_ac)
      call put_var  ('ab_ac'    ,ab_ac    ,bc% data% ab_ac)
      call put_var  ('ob_ob_ac' ,ob_ob_ac ,bc% data% ob_ob_ac)
      call put_var  ('oa_oa_ac' ,oa_oa_ac ,bc% data% oa_oa_ac)
      call put_var  ('ab_ab_ac' ,ab_ab_ac ,bc% data% ab_ab_ac)
      call put_var  ('ob_oa_ac' ,ob_oa_ac ,bc% data% ob_oa_ac)
      call put_var  ('ob_ab_ac' ,ob_ab_ac ,bc% data% ob_ab_ac)
      call put_var  ('oa_ab_ac' ,oa_ab_ac ,bc% data% oa_ab_ac)
      call put_var  ('bc_b'     ,bc_b     ,bc% data% bc_b)
      call put_var  ('bc_a'     ,bc_a     ,bc% data% bc_a)
      call put_var  ('o_err '   ,o_err    ,bc% data% o_err )
      call put_var  ('o_err_ac' ,o_err_ac ,bc% data% o_err_ac)
      call put_var  ('b_err '   ,b_err    ,bc% data% b_err )
      call put_var  ('b_err_ac' ,b_err_ac ,bc% data% b_err_ac)

      !-----------
      ! close file
      !-----------
      call close_bc (bc% h)
    endif

  contains
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_vari (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    integer          ,intent(in) :: values (:,:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_vari
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine put_var (name, varid, values)
    character(len=*) ,intent(in) :: name
    integer          ,intent(in) :: varid
    real(wp)         ,intent(in) :: values (:,:,:)

      status = nf90_put_var (bc%h% ncid ,varid, values)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: put_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine put_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_dim (name, size, dimid)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(in)  :: size
    integer          ,intent(out) :: dimid

      status = nf90_def_dim (bc%h% ncid ,name ,size ,dimid)

      if (status /= NF90_NOERR) then
        write(0,*)   'write_bcor_file: def_dim '//trim(name)//' : size =',size
        call finish ('write_bcor_file: def_dim '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

    end subroutine def_dim
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    subroutine def_var (name, varid, longname)
    character(len=*) ,intent(in)  :: name
    integer          ,intent(out) :: varid
    character(len=*) ,intent(in)  :: longname

      status = nf90_def_var (bc%h% ncid ,name ,NF90_FLOAT, dimid, varid)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if

      status = nf90_put_att (bc%h% ncid ,varid, 'longname', longname)

      if (status /= NF90_NOERR) then
        call finish ('write_bcor_file: def_var '//trim(name),&
                      trim(nf90_strerror(status))            )
      end if
    end subroutine def_var
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine write_bcor_file

!==============================================================================

  subroutine aircraft_bc_bfg (obs)
  !------------------------------------------------------------------------
  ! Radiance bias correction routine to be called before first guess check.
  ! Check for missing entries in bc file, extend file
  ! Apply bias correction to observations.
  !------------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation

    !----------------
    ! local variables
    !----------------
    integer ,parameter       :: mp = 10000   ! max. number of new planes
    type(t_tmp)              :: mis_pe (mp)  ! missing entries in bc file
    type(t_tmp)              :: mis_all(mp)  ! missing entries in bc file
    integer(i8)              :: hash         ! integer representation of name
    integer                  :: n            ! counter / PE
    integer                  :: m            ! counter all PEs
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    integer                  :: i, j, k      ! loop indices
    type(t_plat_bc) ,pointer :: plat (:)     ! temporary
    type (t_bc)     ,pointer :: data (:,:,:) ! temporary
    real(wp)                 :: rawobs       ! raw observation
    real(wp)                 :: bcor         ! bias correction
    integer                  :: ibc          ! bias correction index
    integer                  :: ip           ! phase index
    integer                  :: state        ! state to set if no BC available

    !---------------------------------------------------------
    ! relate bias correction file entries to reports in 3dvar.
    ! check for missing entries in bc file.
    ! extend file if necessary.
    !---------------------------------------------------------
    if (biascor_mode /= 0) then
      !----------------------------------------------------------------------
      ! gather aircraft flight-track information and (re-)assign flight phase
      !----------------------------------------------------------------------
      if (scan_tracks) then
        call aircraft_tracks_gather (obs)
        call aircraft_tracks_scan   (obs)
        call aircraft_tracks_delete ()
      endif
      !-----------------------------------------------
      ! relate bias correction file entries to reports
      ! check for missing entries on this PE
      !-----------------------------------------------
      n = 0
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
spot:   do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype  /= OT_AIREP) cycle
          if (obs% o(ib)% spot(is)% hd% codetype == 146     ) cycle ! MODES
          obs% o(ib)% spot(is)% bc_index = 0
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            select case (obs% o(ib)% varno(i))
            case (VN_T, VN_RH)
            case default
              cycle
            end select
            hash = transfer (obs% o(ib)% spot(is)% statid, hash)
            do j = 1, size(aircraft_bc% plat)
              if  (hash == aircraft_bc% plat(j)% hash) then
                obs% o(ib)% spot(is)% bc_index = j
                cycle spot
              endif
            end do
            do j = 1, n
              if (mis_pe (j)% hash == hash) then
                obs% o(ib)% spot(is)% bc_index = -j
                cycle spot
              endif
            end do
            n = n + 1
            if (n > mp) call finish('aircraft_bc_bfg','n > mp')
            obs% o(ib)% spot(is)% bc_index = - n
            mis_pe (n)% statid = obs% o(ib)% spot(is)% statid
            mis_pe (n)% hash   = hash
            mis_pe (n)% i      = n
            mis_pe (n)% j      = 0
            cycle spot
          end do
        end do spot
      end do
      !---------------------------------------
      ! cross-check missing entries on all PEs
      !---------------------------------------
      m = p_sum (n)
      if (m > mp) call finish('aircraft_bc_bfg','m > mp')
      if (m > 0) then
        call p_gather_tmp (mis_pe (1:n), mis_all(1:m), dace% pio)
        if (dace% lpio) then
          k             = 1
          i             = 1
          mis_all(i)% j = k
          mis_pe (k)    = mis_all(i)
outer:    do i = 2, m
            do j = 1, k
              if (mis_all(i)% hash == mis_pe(j)% hash) then
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
          i = size (aircraft_bc% plat)
          plat => aircraft_bc% plat
          data => aircraft_bc% data
          allocate (aircraft_bc% plat (i+k))
          allocate (aircraft_bc% data (nph,nob,i+k))
          aircraft_bc% n = i+k
          aircraft_bc% plat(    :i) = plat
          aircraft_bc% data(:,:,:i) = data
          do j = 1, k
            aircraft_bc% plat (i+j)% hash   = mis_pe (j)% hash
            aircraft_bc% plat (i+j)% statid = mis_pe (j)% statid
          end do
          deallocate (plat, data)
          mis_all (1:m)% j = mis_all (1:m)% j + i
        endif
        !--------------------------------------
        ! redistribute information to other PEs
        !--------------------------------------
        call p_scatter_tmp (mis_all(1:m), mis_pe(1:n), dace% pio)
        if (dace% lpio) m = i+k
        call p_bcast (m, dace% pio)
        if (.not.dace% lpio) then
          deallocate (aircraft_bc% plat)
          deallocate (aircraft_bc% data)
          allocate   (aircraft_bc% plat (m))
          allocate   (aircraft_bc% data (nph,nob,m))
        endif
        call p_bcast_plat (aircraft_bc% plat, dace% pio)
        call p_bcast_data (aircraft_bc% data, dace% pio)
        aircraft_bc% n = m
        !-----------------------------------
        ! relate reports to new file entries
        !-----------------------------------
        do ib = 1, size (obs% o)
          if (obs% o(ib)% pe /= dace% pe) cycle
          do is = 1, obs% o(ib)% n_spot
            if (obs% o(ib)% spot(is)% hd% obstype  /= OT_AIREP) cycle
            if (obs% o(ib)% spot(is)% hd% codetype == 146     ) cycle ! MODES
            if (obs% o(ib)% spot(is)% bc_index < 0)      &
                obs% o(ib)% spot(is)% bc_index = mis_pe( &
               -obs% o(ib)% spot(is)% bc_index)% j
          end do
        end do
      endif
    endif

    !-----------------------------------
    ! apply bias-correction within 3dvar
    !-----------------------------------
    if (biascor_mode >= BC_FG) then
      state = rept_use (OT_AIREP)% use (CHK_BIASCOR)
      do ib = 1, size (obs% o)
        if (obs% o(ib)% pe /= dace% pe) cycle
        do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype  /= OT_AIREP) cycle
          if (obs% o(ib)% spot(is)% hd% codetype == 146     ) cycle ! MODES
          ibc = obs% o(ib)% spot(is)% bc_index
          if (ibc <= 0)                                       cycle
          select case (obs% o(ib)% spot(is)% phase)
          case (BUFR_LEVEL,BUFR_LEVEL_HW)   ! Level flight
            ip = PH_LEV
          case (BUFR_ASCEND)     ! Ascending (ASC)
            ip = PH_ASC
          case (BUFR_DESCEND)     ! Descending (DES)
            ip = PH_DES
          case (BUFR_MISSING)     ! Missing value
            ip = PH_UNSP
          case default ! Unsteady, reserved
            ip = 0
          end select
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            select case (obs% o(ib)% varno(i))
            case (VN_T)
              j = OB_T
            case (VN_RH)
!             j = OB_RH
              cycle
            case default
              cycle
            end select
            if (ip > 0) then
              k = ip
              if (k==PH_UNSP .and. obs%o(ib)% olev(i) < force_level) k = PH_LEV
              rawobs = obs% o(ib)% body(i)% o  - &
                       obs% o(ib)% body(i)% bc
              bcor = aircraft_bc% data(k,j,ibc)% bc_b
              if (bcor == NF90_FILL_FLOAT) then
                bcor = 0._wp
                call decr_use (obs% o(ib)% body(i)% use, &
                               check = CHK_BIASCOR,      &
                               state = state,            &
                               lflag = .true.            )
              endif
              obs% o(ib)% body(i)% bc =        - bcor
              obs% o(ib)% body(i)% o  = rawobs - bcor
            else
              call decr_use (obs% o(ib)% body(i)% use, &
                             check = CHK_BIASCOR,      &
                             state = state,            &
                             lflag = .true.            )
            endif
            exit
          end do
        end do
      end do
    endif
  end subroutine aircraft_bc_bfg

!==============================================================================

  subroutine aircraft_bc_aan (obs, y)
  !-----------------------------------------------------------------------
  ! Bias correction routine to be called after analysis
  ! Update bias correction statistics.
  ! Write updated correction coefficient file.
  !-----------------------------------------------------------------------
  type (t_obs_set) ,intent(inout) :: obs  ! observation
  type (t_vector)  ,intent(in)    :: y    ! background, observation space

    integer             :: i, j, k   ! index
!   integer             :: n         ! number of aircrafts in statistics
    integer             :: ib        ! box index
    integer             :: is        ! spot index
    integer             :: ip        ! phase index
    integer             :: ibc       ! plane index
    real(wp)            :: bg        ! background
    real(wp)            :: an        ! analysis
    real(wp)            :: o         ! raw observation
    real(wp)            :: eo        ! observation error
    real(wp)            :: eb        ! background  error
    real(wp)            :: dob       ! obs - bg
    real(wp)            :: doa       ! obs - ana
    real(wp)            :: dab       ! ana - bg
    real(wp)            :: dt        ! time since last update of statistics
    real(wp)            :: f         ! weight for accumulated statistics
    integer             :: ier       ! error return flag
    integer,    pointer :: iperm (:) ! permutation index array for sorting
    type(t_bc), pointer :: pp(:,:,:) ! pointer to statistics file entries
    type(t_bc), pointer :: p     (:) ! pointer to statistics file entries

    if (biascor_mode == 0) return

!   if (run_type < 2)      return    ! update file only in analysis cycle

!   n = size(aircraft_bc% plat)

    !-----------------
    ! set sums to zero
    !-----------------
    pp => aircraft_bc% data(:,:,:)
    pp% n     = 0
    pp% ob    = 0._wp  ! obs-bg  mean deviation
    pp% oa    = 0._wp  ! obs-ana mean deviation
    pp% ab    = 0._wp  ! ana-obs mean deviation
    pp% ob_ob = 0._wp  ! obs-bg  variance
    pp% oa_oa = 0._wp  ! obs-ana variance
    pp% ab_ab = 0._wp  ! ana-obs variance
    pp% ob_oa = 0._wp  !
    pp% ob_ab = 0._wp  !
    pp% oa_ab = 0._wp  !
    pp% o_err = 0._wp  !
    pp% b_err = 0._wp  !

    !-----------------------
    ! Update bias statistics
    !-----------------------
    do ib = 1, size (obs% o)
      if (obs% o(ib)% pe /= dace% pe) cycle
      do is = 1, obs% o(ib)% n_spot
        if (obs% o(ib)% spot(is)% hd% obstype  /= OT_AIREP) cycle
        if (obs% o(ib)% spot(is)% hd% codetype == 146     ) cycle ! MODES
        ibc = obs% o(ib)% spot(is)% bc_index
        if (ibc <= 0)                                       cycle
        if (obs% o(ib)% spot(is)% sl_bg < fr_land_bc)       cycle
        select case (obs% o(ib)% spot(is)% phase)
        case (BUFR_LEVEL,BUFR_LEVEL_HW)   ! Level flight
          ip = PH_LEV
        case (BUFR_ASCEND)     ! Ascending (ASC)
          ip = PH_ASC
        case (BUFR_DESCEND)     ! Descending (DES)
          ip = PH_DES
        case (BUFR_MISSING)     ! Missing value
          ip = PH_UNSP
        case default ! Unsteady, reserved
          ip = 0
        end select
        if (ip > 0) then
          do i = obs% o(ib)% spot(is)% o% i + 1,                       &
                 obs% o(ib)% spot(is)% o% i + obs% o(ib)% spot(is)% o% n
            select case (obs% o(ib)% varno(i))
            case (VN_T)
              j = OB_T
            case (VN_RH)
              j = OB_RH
            case default
              cycle
            end select
            k = ip
            if (k==PH_UNSP .and. obs%o(ib)% olev(i) < force_level) k = PH_LEV
            p => aircraft_bc% data(k,:,ibc)
            select case (obs% o(ib)% body(i)% use% state)
            case (STAT_PASSIVE, STAT_ACTIVE_0I:STAT_ACCEPTED)
              bg = obs% o(ib)% body(i)% bg
              an = y% s(ib)% x(i)
              o  = obs% o(ib)% body(i)% o - obs% o(ib)% body(i)% bc
              eo = obs% o(ib)% body(i)% eo
              eb = obs% o(ib)% body(i)% eb
              dob  = o  - bg
              doa  = o  - an
              dab  = an - bg
              p(j)% n      = p(j)% n      + 1
              p(j)% ob     = p(j)% ob     + dob
              p(j)% oa     = p(j)% oa     + doa
              p(j)% ab     = p(j)% ab     + dab
              p(j)% ob_ob  = p(j)% ob_ob  + dob * dob
              p(j)% oa_oa  = p(j)% oa_oa  + doa * doa
              p(j)% ab_ab  = p(j)% ab_ab  + dab * dab
              p(j)% ob_oa  = p(j)% ob_oa  + dob * doa
              p(j)% ob_ab  = p(j)% ob_ab  + dob * dab
              p(j)% oa_ab  = p(j)% oa_ab  + doa * dab
              p(j)% o_err  = p(j)% o_err  + eo
              p(j)% b_err  = p(j)% b_err  + eb
            end select
          end do
        endif
      end do
    end do

    !----------------
    ! sum up over PEs
    !----------------
    aircraft_bc% data% n      = p_sum (aircraft_bc% data% n)
    aircraft_bc% data% ob     = p_sum (aircraft_bc% data% ob)
    aircraft_bc% data% oa     = p_sum (aircraft_bc% data% oa)
    aircraft_bc% data% ab     = p_sum (aircraft_bc% data% ab)
    aircraft_bc% data% ob_ob  = p_sum (aircraft_bc% data% ob_ob)
    aircraft_bc% data% oa_oa  = p_sum (aircraft_bc% data% oa_oa)
    aircraft_bc% data% ab_ab  = p_sum (aircraft_bc% data% ab_ab)
    aircraft_bc% data% ob_oa  = p_sum (aircraft_bc% data% ob_oa)
    aircraft_bc% data% ob_ab  = p_sum (aircraft_bc% data% ob_ab)
    aircraft_bc% data% oa_ab  = p_sum (aircraft_bc% data% oa_ab)
    aircraft_bc% data% o_err  = p_sum (aircraft_bc% data% o_err )
    aircraft_bc% data% b_err  = p_sum (aircraft_bc% data% b_err )

    !-----------------------------------------
    ! weight factor for accumulated statistics
    !-----------------------------------------
    dt = days (ana_time - time_cyyyymmddhhmm (aircraft_bc% h% last_date))
    if (aircraft_bc% h% t_decay(1) > 0._wp .and. dt > 0._wp) then
      f = exp ( - dt / aircraft_bc% h% t_decay(1))
    else
      f = 1._wp
    endif

    !-----------------------------------
    ! Update accumulated bias statistics
    !-----------------------------------
    pp => aircraft_bc% data(:,:,:)
    !-------------------------------
    ! rescale accumulated statistics
    !-------------------------------
    pp% n_ac       = pp% n_ac     * f
    pp% ob_ac      = pp% ob_ac    * pp% n_ac
    pp% oa_ac      = pp% oa_ac    * pp% n_ac
    pp% ab_ac      = pp% ab_ac    * pp% n_ac
    pp% ob_ob_ac   = pp% ob_ob_ac * pp% n_ac
    pp% oa_oa_ac   = pp% oa_oa_ac * pp% n_ac
    pp% ab_ab_ac   = pp% ab_ab_ac * pp% n_ac
    pp% ob_oa_ac   = pp% ob_oa_ac * pp% n_ac
    pp% ob_ab_ac   = pp% ob_ab_ac * pp% n_ac
    pp% oa_ab_ac   = pp% oa_ab_ac * pp% n_ac
    pp% o_err_ac   = pp% o_err_ac * pp% n_ac
    pp% b_err_ac   = pp% b_err_ac * pp% n_ac
    !-----------
    ! accumulate
    !-----------
    pp% n_ac       = pp% n_ac     + pp% n
    pp% ob_ac      = pp% ob_ac    + pp% ob
    pp% oa_ac      = pp% oa_ac    + pp% oa
    pp% ab_ac      = pp% ab_ac    + pp% ab
    pp% ob_ob_ac   = pp% ob_ob_ac + pp% ob_ob
    pp% oa_oa_ac   = pp% oa_oa_ac + pp% oa_oa
    pp% ab_ab_ac   = pp% ab_ab_ac + pp% ab_ab
    pp% ob_oa_ac   = pp% ob_oa_ac + pp% ob_oa
    pp% ob_ab_ac   = pp% ob_ab_ac + pp% ob_ab
    pp% oa_ab_ac   = pp% oa_ab_ac + pp% oa_ab
    pp% o_err_ac   = pp% o_err_ac + pp% o_err
    pp% b_err_ac   = pp% b_err_ac + pp% b_err
    !--------
    ! rescale
    !--------
    where (pp% n > 0)
      pp% ob       = pp% ob       / pp% n
      pp% oa       = pp% oa       / pp% n
      pp% ab       = pp% ab       / pp% n
      pp% ob_ob    = pp% ob_ob    / pp% n
      pp% oa_oa    = pp% oa_oa    / pp% n
      pp% ab_ab    = pp% ab_ab    / pp% n
      pp% ob_oa    = pp% ob_oa    / pp% n
      pp% ob_ab    = pp% ob_ab    / pp% n
      pp% oa_ab    = pp% oa_ab    / pp% n
      pp% o_err    = pp% o_err    / pp% n
      pp% b_err    = pp% b_err    / pp% n
    endwhere

    where (pp% n_ac > 0)
      pp% ob_ac    = pp% ob_ac    / pp% n_ac
      pp% oa_ac    = pp% oa_ac    / pp% n_ac
      pp% ab_ac    = pp% ab_ac    / pp% n_ac
      pp% ob_ob_ac = pp% ob_ob_ac / pp% n_ac
      pp% oa_oa_ac = pp% oa_oa_ac / pp% n_ac
      pp% ab_ab_ac = pp% ab_ab_ac / pp% n_ac
      pp% ob_oa_ac = pp% ob_oa_ac / pp% n_ac
      pp% ob_ab_ac = pp% ob_ab_ac / pp% n_ac
      pp% oa_ab_ac = pp% oa_ab_ac / pp% n_ac
      pp% o_err_ac = pp% o_err_ac / pp% n_ac
      pp% b_err_ac = pp% b_err_ac / pp% n_ac
    endwhere

    !-----
    ! sort
    !-----
    allocate (iperm (aircraft_bc% n))
    call sort (aircraft_bc% plat% statid, iperm, 1, ier)
    aircraft_bc% plat = aircraft_bc% plat (    iperm)
    aircraft_bc% data = aircraft_bc% data (:,:,iperm)
    deallocate (iperm)

    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'( )')
      write(6,'(a)')      '  Aircraft temperature bias correction'
    end if

    call remove_inactive (aircraft_bc)
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write(6,'( )')
      write(6,'(a,a)')    '    last_date = ',aircraft_bc% h% last_date
      write(6,'(a,a)')    '    ana_date  = ',cyyyymmddhhmm(ana_time)
      write(6,'(a,f11.2)')'    t_decay   =',    aircraft_bc% h%    t_decay(1)
      write(6,'(a,f13.4)')'    dt        =',                       dt
      write(6,'(a,f13.4)')'    f         =',                       f
      write(6,'(a,i8)')   '    n         =',sum(aircraft_bc% data% n)
      write(6,'(a,f11.2)')'    n_ac      =',sum(aircraft_bc% data% n_ac)
      write(6,'( )')
      if (verbose > 0) then
        !---------------------------------
        ! Verbosity of bias correction:
        ! 1 = print only non-empty entries
        ! 2 = print all
        !---------------------------------
        write(6,'(a)') &
             '           ------ entries ------ --- bias correction --- --- innovation rmse ---'
        write(6,'(a)') &
             'statid     ASC   DES   LEV  UNSP   ASC   DES   LEV  UNSP   ASC   DES   LEV  UNSP'
        write(6,'( )')
        do i = 1, aircraft_bc% n
          if (verbose > 1 .or. sum (aircraft_bc% data(:,1,i)% n_ac) > 0._wp) &
          write(6,'(a,4f6.0,8f6.2)')&
                               aircraft_bc% plat    (i)% statid,    &
                               aircraft_bc% data(:,1,i)% n_ac,      &
                               aircraft_bc% data(:,1,i)% bc_a,      &
             sqrt (max (0._wp, aircraft_bc% data(:,1,i)% ob_ob_ac   &
                             - aircraft_bc% data(:,1,i)% ob_ac ** 2))
        end do
        write(6,'( )')
      end if
    endif
    !------------------------------------------
    ! Write updated correction coefficient file
    !------------------------------------------
    call write_bcor_file (aircraft_bc)
  end subroutine aircraft_bc_aan

!------------------------------------------------------------------------------
  subroutine aircraft_tracks_gather (obs)
    !-----------------------------------------
    ! Gather aircraft flight-track information
    ! and distribute to all processors
    !-----------------------------------------
    type (t_obs_set) ,intent(inout) :: obs  ! observations
    !----------------
    ! local variables
    !----------------
    integer                  :: ib           ! box index
    integer                  :: is           ! spot index
    real                     :: level        ! height
    integer                  :: i            ! loop index
    integer                  :: ier          ! error status
    integer                  :: m            ! auxiliary
    integer                  :: n            ! auxiliary
    type(t_spot)   ,pointer  :: s            ! auxiliary spot pointer
    type(t_track)  ,pointer  :: t            ! auxiliary track pointer
    type(t_track)  ,pointer  :: tmp(:)       ! auxiliary array for sorting
    integer, allocatable     :: idx(:)       ! index arrays for sorting
    integer                  :: mtracks      ! Current size of tracks(:)
    integer, parameter       :: m0 = 100     ! Initial no. of spots

    if (associated (tracks)) call aircraft_tracks_delete ()
    mtracks = m0
    allocate (tracks(mtracks))

    do ib = 1, size (obs% o)
       if (obs% o(ib)% pe /= dace% pe) cycle
spot:  do is = 1, obs% o(ib)% n_spot
          if (obs% o(ib)% spot(is)% hd% obstype  /= OT_AIREP) cycle
          if (obs% o(ib)% spot(is)% hd% codetype == 146     ) cycle ! MODES
          s => obs% o(ib)% spot(is)
          do i = s% o% i + 1, s% o% i + s% o% n
            !--------------------------------------------------
            ! Use level of first valid observation at this spot
            ! (may need checking against plevel if available).
            !--------------------------------------------------
            level = obs%o(ib)% olev(i)
            !------------------------------
            ! Increase array size if needed
            !------------------------------
            if (ntracks >= mtracks) then
               tmp => tracks
               m = nint (mtracks*1.5) + 1
               allocate (tracks(m))
               tracks(1:mtracks) = tmp
               mtracks = m
               deallocate (tmp)
            end if
            !------------------------------
            ! Store essential aircraft data
            !------------------------------
            ntracks    =  ntracks + 1
            t          => tracks(ntracks)
            t% statid  =  s% statid
            t% hash    =  transfer (t% statid, t% hash) ! Should consider endianness..
            t% time_oa =  nint (seconds (s% actual_time - ana_time))
            t% pe      =  dace% pe
            t% ib      =  ib
            t% is      =  is
            t% phase   =  s% phase
            t% level   =  level
            t% id      =  s% id
            cycle spot
          end do
       end do spot
    end do
    !-----------------------------
    ! Gather data on I/O processor
    !-----------------------------
    n = p_sum (ntracks)

!print*,dace% pe,'aircraft_tracks_gather: sending to I/O pe, ntracks=', ntracks
!if (dace% lpio) then
!   print *, "total track entries: n=", n
!end if

    allocate (tmp(n))
    call p_gather_track (tracks(1:ntracks), tmp(1:n), dace% pio)
    deallocate (tracks)
    allocate (tracks(n))
    !-------------------------------------------------------------------
    ! Sort aircraft flight tracks according to registration id and time.
    !-------------------------------------------------------------------
    if (dace% lpio) then
       allocate (idx(n))
       !--------------------------
       ! First sort w.r.t. 3rd key
       !--------------------------
       call sort (tmp% id,      idx, 1, ier)
       tmp = tmp(idx)
       !--------------------------
       ! Then  sort w.r.t. 2nd key
       !--------------------------
       call sort (tmp% time_oa, idx, 1, ier)
       tmp = tmp(idx)
       !-------------------------
       ! Then sort w.r.t. 1st key
       !-------------------------
       call sort (tmp% statid,  idx, 1, ier)
!      call sort (tmp% hash  ,  idx, 1, ier)  !    Faster, but diff.order
       tracks = tmp(idx)
       deallocate (idx)
    end if
    deallocate (tmp)
    call p_bcast_track (tracks(1:n), dace% pio)
    ntracks = n

  end subroutine aircraft_tracks_gather
!------------------------------------------------------------------------------
  subroutine aircraft_tracks_scan (obs)
    type (t_obs_set) ,intent(inout) :: obs  ! observations
    !----------------
    ! local variables
    !----------------
    integer                  :: i, j    ! loop index
    integer                  :: i0      ! first index of current plane
    integer                  :: n       ! auxiliary variable (counter)
    integer                  :: xphase(ntracks)  ! new flight phase
    real                     :: rate  (ntracks)  ! ascend/descend rate
!   character(len=8)         :: statid  ! station id
    type(t_track)  ,pointer  :: t       ! auxiliary track pointer
    type(t_track)  ,pointer  :: t1      ! auxiliary track pointer
    type(t_track)  ,pointer  :: t2      ! auxiliary track pointer
    type(t_spot)   ,pointer  :: s       ! auxiliary spot pointer
    logical                  :: change  ! true for phase change
    integer(i8)              :: hash    ! hash of station id

    if (ntracks == 0) return
    if (.not. associated (tracks)) &
         call finish ("aircraft_tracks_scan","tracks not allocated!")

    !---------------------------
    ! loop over all observations
    !---------------------------
    xphase = tracks(:ntracks)% phase
    rate   = 0.
    do i = 1, ntracks
      t  => tracks(i)
      if (t% phase == BUFR_MISSING) then
        if (t% level <= force_level) then
            xphase(i)   = BUFR_LEVEL
        else
          if (i < ntracks) then
            t1 => tracks(i+1)
            if (t1% statid == t% statid .and.      &
                t1% time_oa - t% time_oa < time_sep) then
                rate(i) = (t1% level   - t% level  ) &
                  / max(1, t1% time_oa - t% time_oa)
              if (rate(i) < -rate_asc      ) then
                xphase(i)   = BUFR_ASCEND
              else if (rate(i) > rate_desc ) then
                xphase(i)   = BUFR_DESCEND
              else if (abs(rate(i)) < rate_lev .and. t% level <= 31000.) then
                xphase(i)   = BUFR_LEVEL
              else
                xphase(i)   = t1% phase
              endif
            endif
          endif
        endif
        if (xphase(i) == BUFR_MISSING) then
          if (i > 1) then
            t2 => tracks(i-1)
            if (t% statid == t2% statid      .and. &
                t% time_oa - t2% time_oa < time_sep) then
                rate(i) = (t% level   - t2% level  ) &
                  / max(1, t% time_oa - t2% time_oa)
              if (rate(i) < -rate_asc      ) then
                xphase(i)   = BUFR_ASCEND
              else if (rate(i) > rate_desc ) then
                xphase(i)   = BUFR_DESCEND
              else if (abs(rate(i)) < rate_lev .and. t% level <= 31000.) then
                xphase(i)   = BUFR_LEVEL
              else
                xphase(i)   = t2% phase
              endif
            endif
          endif
        endif

        if (t% pe == dace% pe) then
          !----------------------------------------------------
          ! overwrite phase if observation is on this processor
          !----------------------------------------------------
          s => obs% o(t% ib)% spot(t% is)
          s% phase = xphase(i)
        end if
      endif
    end do

    !---------
    ! printout
    !---------
    if (dace% lpio .and. verbose > 2) then
      hash   = 0
      change = .false.
      n  = 0
      i0 = 1
      do i = 1, ntracks
        t  => tracks(i)
        if (t% hash /= hash) then
          n      = n + 1
          if (change) then
            write (6,*)
            do j = i0, i-1
              t1 => tracks(j)
              write(6,*) ' AIRCRAFT-PHASE',n,t1%phase /= xphase(j), t1%statid,&
              xphase(j), t1%phase, t1%time_oa - tracks(i0)%time_oa, t1%level, &
              rate(j)
            end do
          endif
          change = .false.
          hash   = t% hash
          i0     = i
        endif
        if (xphase(i) /= t% phase) change = .true.
      end do
      write(6,*)
      write(6,*) "aircraft_tracks_scan: number of different aircraft =", n
    endif

  end subroutine aircraft_tracks_scan
!------------------------------------------------------------------------------
  subroutine aircraft_tracks_delete ()
    if (associated (tracks)) deallocate (tracks)
    nullify (tracks)
    ntracks = 0
  end subroutine aircraft_tracks_delete
!------------------------------------------------------------------------------

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
!----------------------------------------------------
! subroutine p_bcast_data (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_bc),dimension(:,:,:)
#define VECTOR
#define RANK 3
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_data
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  RANK
#undef  p_bcast_DERIVED
!==============================================================================
!----------------------------------------------------------------------------
! subroutine p_gather_track (sendbuffer,receivebuffer,root,[comm],[recvcounts])
!----------------------------------------------------------------------------
#define DERIVED type(t_track)
#undef  MPI_TYPE
#define p_gather_DERIVED p_gather_track
#include "p_gather_derived.incf"
#undef  DERIVED
#undef  p_gather_DERIVED
!----------------------------------------------------
! subroutine p_bcast_track (buffer, p_source, [comm])
!----------------------------------------------------
#define DERIVED type(t_track),dimension(:)
#define VECTOR
#undef  MPI_TYPE
#define p_bcast_DERIVED p_bcast_track
#include "p_bcast.incf"
#undef  DERIVED
#undef  VECTOR
#undef  p_bcast_DERIVED
!==============================================================================
end module mo_conv_bc
