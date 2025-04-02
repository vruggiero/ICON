!
!+ Data type holding IASI SATEM observations and operations thereon
!
MODULE mo_satem
!
! Description:
! Definition of data type holding IASI SATEM observations
! and operations thereon
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  scan observation input file by generic routine in module mo_obs_netcdf
! V1_3         2008/12/08 Detlef Pingel
!  allow different integration schemes for geopotential height
! V1_4         2009/03/26 Detlef Pingel
!  work around for NEC-SX bug
! V1_9         2010/04/20 Andreas Rhodin
!  TSK_SHRINK in subroutines process: pass parameter 'state' to shrink_report
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  technical changes
! V1_20        2012-06-18 Andreas Rhodin
!  changed comments and names of public routines in module mo_p_output
! V1_28        2014/02/26 Andreas Rhodin
!  changed interface to new_int
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Detlef Pingel  DWD 08/2008 original code
!==============================================================================

!-------------------------------------------------
! uncomment to exit in case of unknown BUFR codes:
!-------------------------------------------------

  !=============
  ! modules used
  !=============
  !------------------------
  ! general purpose modules
  !------------------------
  use mo_kind,          only: wp,                &! working precision
                              dp                  ! double precision
  use mo_exception,     only: finish              ! abort routine
  use mo_mpi_dace,      only: dace,              &! MPI group info
                              p_bcast             ! generic broadcast routine
  use mo_namelist,      only: position_nml,      &! position namelist
                              POSITIONED,        &! return value
                              nnml                ! namelist file unit number
  use mo_physics,       only: gacc,              &! gravity acceleration
                              Rgas=>R,           &! gas constant
                              tv_t_q              ! get tv from t, q
  use mo_run_params,    only: path_file           ! concatenate path/filename
  use mo_p_output,      only: oline,             &! output line buffer
                              iol,               &! indexof next line to write
                        nl => nextline,          &! increment line number
                              flush_buf_pio       ! write buffer on I/O PE
  use mo_cntrlvar,      only: tq_tvgh             ! generalized humidity transform
  use mo_algorithms,    only: init_splinex,      &! routines for
                              splint              ! spline interpolation
  !---------------------------------------
  ! atmospheric state data type definition
  !---------------------------------------
  use mo_atm_state,     only: t_atm               ! atmospheric state data type
  use mo_time,          only: t_time,            &! date&time data type
                              init_time,         &! initialise time data type
                              cyyyymmdd,chhmmss   ! conversion routines
  !-----------------
  ! matrix data type
  !-----------------
  use mo_dec_matrix,    only: t_vector_segm       ! segment data type
  !----------------------
  ! observation data type
  !----------------------

  use mo_obs_tables,    only: decr_rpt_use,      &! degrade status of report
                              check_report_0,    &! basic checks on reports
                              check_report_1,    &! basic checks on reports
                              rept_use            ! status flag table
  use mo_t_use,         only: STAT_DISMISS,      &! status flag: dismiss report
                              STAT_ACTIVE,       &!            : active
                              CHK_INSDAT,        &!            : insufficient data
                              CHK_NONE,          &!            : no check fired
                              CHK_DOMAIN,        &!            : inconsis. profiles
                              t_use               ! data type to hold state
  use mo_t_datum,       only: t_datum             ! date+time derived type
  use mo_obs_set,       only: t_obs_block         ! obs data type
  use mo_t_obs,         only: t_obs,             &! observation data type
                              t_spot,            &! observation meta data type
                              t_head,            &! component of t_spot
                              new_spot,          &! reserve memory for meta data
                              set_xuv,           &! set unit vectors, solar zenith
                              new_obs,           &! reserve memory for observations
                              new_int,           &! reserve memory for interp.space
                              new_par,           &! reserve memory for parameters
                              shrink_report,     &! remove passive observations
                              OBS_H,             &! interp. observation type: geop.
                              ITY_ICOL,          &! interpolation type: column
                              source,            &! table of observation files
                              n_source,          &! number of     source files
             SATEM_obstype => SATEM,             &! SATEM obs. operator ident.
                              FT_NETCDF,         &! input file type
                              set_vqc_insitu,    &! bounds for vqc routine
                              set_int_insitu,    &! set interpolation space
  !---------------------
  ! predefined constants
  !---------------------
                              TSK_INIT,          &! task value: initialise module
                              TSK_READ,          &! read observations
                              TSK_SET_CHR,       &! set observation characteristics
                              TSK_SETUP_COLS,    &! setup model columns required
                              TSK_SETUP_FULL,    &! final setup of PSAS-space
                              TSK_SETUP_FUL0,    &! setup descriptionof PSAS-space
                              TSK_K,             &! derive tl model
                              TSK_R,             &! setup observational error
                              TSK_SHRINK          ! release unused obs. in report
  use mo_fdbk_tables,   only: VN_Q,              &!         q              code
                              VN_T,              &!         T              code
                              VN_Z,              &!         geopotential   code
                              OT_SATEM            ! report type value for SATEM
  use mo_t_col,         only: t_cols,            &! model columns data type
                              COL_TV,            &! required obervations virt. t
                              COL_GEOH,          &! required obervations geop. h.l.
                              COL_PH,            &! required obervations geop. h.l.
                              COL_RH              ! required obervations rel. hum.
  use mo_satid,         only: mnem                ! derive: mnemonic <-> satid
  use mo_grid_intpol,   only: idx_init            ! determine model indices for location
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,      only: nf90_open,             &!
                         nf90_close,            &!
                         nf90_inq_dimid,        &!
                         nf90_inq_varid,        &!
                         nf90_Inquire_Dimension,&!
                         nf90_get_var,          &!
                         nf90_strerror,         &!
                         NF90_NOWRITE,          &! mode flag to open a dataset
                         NF90_NOERR              ! status return value: no error

  implicit none
!==============================================================================
  !----------------
  ! public entities
  !----------------
  private
  public :: t_satem                   ! SATEM data type
  public :: process_satem             ! routine to process IASI SATEM data
  public :: read_satem_netcdf         ! read IASI SATEM data (fov,metadata,profiles)
  public :: read_IASI_l2_flags_netcdf ! read IASI SATEM data (fov,metadata)
!==============================================================================
  real(wp), parameter     :: invalid = -999._wp    ! invalid value
  integer,  parameter     :: nmxl    = 10          ! maximal number of thickness layers
  !-----------------------------------
  ! private IASI SATEM data type
  !-----------------------------------
  type t_satem_level                               ! obs for one IASI SATEM level
    real(wp)              :: pn      = invalid     ! PRESSURE (VERT.LOCATION)
    real(wp)              :: tdbt    = invalid     ! TEMPERATURE/DRY BULB TEMPERATURE
    real(wp)              :: tdbt2   = invalid     ! TEMPERATURE/DRY BULB TEMPERATURE
    real(wp)              :: var_t   = invalid     ! OBS. ERROR TEMPERATURE
    real(wp)              :: var_q   = invalid     ! OBS. ERROR SPEC. HUMIDITY
    real(wp)              :: pn0     = invalid     ! PRESSURE (VERT.LOCATION)
    real(wp)              :: mixr    = invalid     ! MIXING RATIO (humidity)
    real(wp)              :: mixr2   = invalid     ! MIXING RATIO (humidity)
  end type t_satem_level

  type t_satem                                     ! SATEM data type
    character(len=64)     :: file        = ''      ! Input file name
    type (t_time)         :: time                  ! satem date
    integer               :: sat_id      = 0       ! SATELLITE IDENTIFIER
    integer               :: ogc         = 0       ! ORIGINATING/GENERATING CENTRE
    integer               :: sain        = 0       ! SATELLITE INSTRUMENTS
    integer               :: sacl        = 0       ! SATELLITE CLASSIFIKATION
    integer               :: orbnu       = 0       ! ORBIT NUMBER
    integer               :: scaln       = 0       ! SCAN LINE NUMBER
    integer               :: fievn       = 0       ! FIELD VIEW NUMBER
    integer               :: sind        = 0       ! SUPERADIABATIC INDICATOR
    integer               :: lase        = 0       ! LAND/SEA QUALIFIER
    integer               :: dani        = 0       ! DAY/NIGHT QUALIFIER
    integer               :: i4          = 0       ! SATELLITE DATA PROCESSING TECHNIQUE USED
    integer               :: sgli        = 0       ! SUN-GLINT INDICATOR
    integer               :: FLG_CLDFRM  = 0       ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
    integer               :: FLG_CLDSUM  = 0       ! INSTRUMENT DETECTING CLOUDS
    integer               :: FLG_IASIBAD = 0       ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
    integer               :: FLG_040195  = 0       ! QUALITY AND COMPLETENESS OF RETRIEVAL
    integer               :: FLG_040196  = 0       ! RETRIEVAL CHOICE INDICATOR
    integer               :: FLG_040197  = 0       ! SATELLITE MANOEUVRE INDICATOR
    integer               :: FLG_040198  = 0       ! SELECTION OF BACKGROUND STATE
    integer               :: vtsa        = 0       ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
    integer               :: noob        = 0       ! NUMBER OF OBSERVATIONS
    integer               :: noob0       = 0       ! NUMBER OF OBSERVATIONS
    integer               :: vtsa0       = 0       ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
    real(wp)              :: lah         = invalid ! LATITUDE (HIGH ACCURACY)
    real(wp)              :: loh         = invalid ! LONGITUDE (HIGH ACCURACY)
    real(wp)              :: szan        = invalid ! SATELLITE ZENITH ANGLE
    real(wp)              :: da          = invalid ! BEARING OR AZIMUTH
    real(wp)              :: soza        = invalid ! SOLAR ZENITH ANGLE
    real(wp)              :: das         = invalid ! SOLAR AZIMUTH
    real(wp)              :: ppp         = invalid ! PRESSURE
    real(wp)              :: tsts        = invalid ! SKIN TEMPERATURE
    real(wp)              :: hp          = invalid ! HEIGHT OF STATION
    integer               :: nlevel      =  0      ! number of samples
  end type t_satem

  !---------------
  ! NetCDF indices
  !---------------
  integer                 :: ncid        ! NetCDF file id
  character(len=128)      :: ncfile      ! NetCDF file name
  integer                 :: status      ! NetCDF return variable from functions
  integer                 :: stanc  (3)  ! NetCDF start  parameter
  integer                 :: counc  (3)  ! NetCDF count  parameter
  integer                 :: strnc  (3)  ! NetCDF stride parameter

  !=========
  ! namelist
  !=========
  !--------
  ! general
  !--------
  logical,save            :: layer_thickness = .true.       ! obs: layer thickness (T) <-> t, q (F)
  logical                 :: smooth_satem    = .false.      ! smooth retrieved t, q profiles
  integer                 :: verbose         = 1            ! 0: silent; >=1 verbose
  !----------------------------------
  ! input, SATEM layer specification
  !----------------------------------
  real(wp)                :: satem_init_levs(nmxl) = -1._wp ! boundaries [Pa] of thickness layers
  integer,save            :: n_clear_satem_b = -100000000   ! first clear spot to be processed
  integer,save            :: n_clear_satem_e =  100000000   ! last clear spot to be processed
  integer                 :: obs_layer_method               ! method to derive geop. from t
  !--------------------
  ! observational error
  !--------------------
  real(wp)                :: obs_err_abs_t = 1._wp          ! obs. error temperature
  real(wp)                :: obs_err_rel_q = 1._wp          ! obs. error relative himidity

  namelist /IASI_SATEM/  verbose,                                        &
                         satem_init_levs, smooth_satem, obs_err_abs_t,   &
                         obs_err_rel_q, layer_thickness,n_clear_satem_b, &
                         n_clear_satem_e, obs_layer_method
  !================================
  ! private constants and variables
  !================================
  integer         ,save   :: n_layer                    ! levels to define increments
  integer         ,save   :: satem_int_size         = 0 ! length of data type t_satem
  integer         ,save   :: satem_levels_int_size  = 0 ! length of data type t_satem_level
  integer         ,save   :: satem_byte_size        = 0 ! length of data type t_satem
  integer         ,save   :: satem_levels_byte_size = 0 ! length of data type t_satem_level
  real(wp),pointer,save   :: satem_levs(:)              ! levels to define increments
  real(wp),pointer,save   :: satem_levs_mid(:)          ! levels to define increments

!==============================================================================
contains
!==============================================================================
  subroutine read_IASI_l2_flags_netcdf (satem,ierr)
    type (t_satem),pointer               :: satem(:)            ! observations data type
    integer               ,intent(out)   :: ierr
    !------------------------------------------------------------------------
    ! IASI SATEM meta data and flags are read from a NetCDF file
    !------------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    integer    :: j              ! index variable
    integer    :: i_file         ! index variable
    integer    :: entry          ! reort ID derived from source/record

    !-------------------------------
    ! extra variables for necdf file
    !-------------------------------
    integer    :: nsatems     ! no. of SATEM profiles in input file
    integer    :: i_level_t   ! no. of levels temperature
    integer    :: i_level_m   ! no. of levels mixing ratio (rel. humidity)
    integer    :: sat_id      ! SATELLITE IDENTIFIER
    integer    :: yyyy        ! YEAR
    integer    :: mo          ! MONTH
    integer    :: dd          ! DAY
    integer    :: hh          ! HOUR
    integer    :: mi          ! MINUTE
    integer    :: ss          ! SECOND
    integer    :: ogc         ! ORIGINATING/GENERATING CENTRE
    integer    :: sain        ! SATELLITE INSTRUMENTS
    integer    :: sacl        ! SATELLITE CLASSIFIKATION
    integer    :: orbnu       ! ORBIT NUMBER
    integer    :: scaln       ! SCAN LINE NUMBER
    integer    :: fievn       ! FIELD VIEW NUMBER
    integer    :: sind        ! SUPERADIABATIC INDICATOR
    integer    :: lase        ! LAND/SEA QUALIFIER
    integer    :: dani        ! DAY/NIGHT QUALIFIE
    integer    :: i4          ! SATELLITE DATA PROCESSING TECHNIQUE USED
    integer    :: sgli        ! SUN-GLINTEGER INDICATOR
    integer    :: FLG_CLDFRM  ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
    integer    :: FLG_CLDSUM  ! INSTRUMENT DETECTING CLOUDS
    integer    :: FLG_IASIBAD ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
    integer    :: FLG_040195  ! QUALITY AND COMPLETENESS OF RETRIEVAL
    integer    :: FLG_040196  ! RETRIEVAL CHOICE INDICATOR
    integer    :: FLG_040197  ! SATELLITE MANOEUVRE INDICATOR
    integer    :: FLG_040198  ! SELECTION OF BACKGROUND STATE
    integer    :: vtsa        ! VERTICAL SIGNIFICANCE SATELL.OBSERV.)
    integer    :: noob        ! NUMBER OF OBSERVATIONS
    integer    :: noob0       ! NUMBER OF OBSERVATIONS
    integer    :: vtsa0       ! VERTICAL SIGNIFICANCE SATELL.OBSERV.)
    real(wp)   :: lah         ! LATITUDE (HIGH ACCURACY)
    real(wp)   :: loh         ! LONGITUDE (HIGH ACCURACY)
    real(wp)   :: szan        ! SATELLITE ZENITH ANGLE
    real(wp)   :: da          ! BEARING OR AZIMUTH
    real(wp)   :: soza        ! SOLAR ZENITH ANGLE
    real(wp)   :: das         ! SOLAR AZIMUTH
    real(wp)   :: ppp         ! PRESSURE
    real(wp)   :: tsts        ! SKIN TEMPERATURE
    real(wp)   :: hp          ! HEIGHT OF STATION

    !-----------------------------------
    ! read from files (from one PE only)
    !-----------------------------------
    ierr=-1
    do i_file = 1, n_source

       !-------------------------
       ! skip inappropriate files
       !-------------------------
       if (source(i_file)% filetype /= FT_NETCDF) cycle
       if (source(i_file)% obstype  /= OT_SATEM)  cycle
       if (source(i_file)% pe       /= dace% pe)  cycle

       write(6,'(a,a)') ' reading ',source(i_file)% file
       !-----------------
       ! open NetCDF file
       !-----------------
       ncfile = path_file (source(i_file)% path, source(i_file)% file)
       status = nf90_open (ncfile, NF90_NOWRITE, ncid)
       if (status /= NF90_NOERR) call error ('nf90_open')
       ierr=status

       !-------------------------------
       ! preset total number of reports
       !-------------------------------
       entry = sum (source(1:i_file-1)% entries)

       !------------------------------------------
       ! read number of satem and number of levels
       !------------------------------------------
       call get_dim (nsatems  , 'BUFR_records')
       call get_dim (i_level_t, 'Loop_000_maxlen')
       call get_dim (i_level_m, 'Loop_001_maxlen')

       print * ,' satems read: ' , nsatems

       ! allocate SATEM data structure
       allocate(satem(nsatems))

       !-------------------------
       ! read netcdf file, header
       !-------------------------
       ! loop over satem profiles
       do j=1,nsatems
          entry = entry + 1
          ! read header info of satem profile nsatems
          stanc  = j          ! read from satem profile j
          counc  = 1
          strnc  = 1

          call get_var_int    ( sat_id      ,'MI1I2')  ! SATELLITE IDENTIFIER
          call get_var_int    ( ogc         ,'MOGC')   ! ORIGINATING/GENERATING CENTRE
          call get_var_int    ( sain        ,'MSAIN')  ! SATELLITE INSTRUMENTS
          call get_var_int    ( sacl        ,'MSACL')  ! SATELLITE CLASSIFIKATION
          call get_var_int    ( yyyy        ,'MJJJ')   ! YEAR
          call get_var_int    ( mo          ,'MMM')    ! MONTH
          call get_var_int    ( dd          ,'MYY')    ! DAY
          call get_var_int    ( hh          ,'MGG')    ! HOUR
          call get_var_int    ( mi          ,'NGG')    ! MINUTE
          call get_var_int    ( ss          ,'MSEC')   ! SECOND
          call get_var_int    ( orbnu       ,'NORBNU') ! ORBIT NUMBER
          call get_var_int    ( scaln       ,'NSCALN') ! SCAN LINE NUMBER
          call get_var_int    ( fievn       ,'NFIEVN') ! FIELD VIEW NUMBER
          call get_var_int    ( sind        ,'MSIND')  ! SUPERADIABATIC INDICATOR
          call get_var_int    ( lase        ,'MLASE')  ! LAND/SEA QUALIFIER
          call get_var_int    ( dani        ,'MDANI')  ! DAY/NIGHT QUALIFIE
          call get_var_int    ( i4          ,'MI4')    ! SATELLITE DATA PROCESSING TECHNIQUE USED
          call get_var_int    ( sgli        ,'MSGLI')  ! SUN-GLINT INDICATOR
          call get_var_int    ( FLG_CLDFRM  ,'040192') ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
          call get_var_int    ( FLG_CLDSUM  ,'040193') ! INSTRUMENT DETECTING CLOUDS
          call get_var_int    ( FLG_IASIBAD ,'040194') ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
          call get_var_int    ( FLG_040195  ,'040195') ! QUALITY AND COMPLETENESS OF RETRIEVAL
          call get_var_int    ( FLG_040196  ,'040196') ! RETRIEVAL CHOICE INDICATOR
          call get_var_int    ( FLG_040197  ,'040197') ! SATELLITE MANOEUVRE INDICATOR
          call get_var_int    ( FLG_040198  ,'040198') ! SELECTION OF BACKGROUND STATE
          call get_var_int    ( vtsa        ,'MVTSA')  ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          call get_var_int    ( noob        ,'NNOOB')  ! NUMBER OF OBSERVATIONS
          call get_var_int    ( noob0       ,'NNOOB0') ! NUMBER OF OBSERVATIONS
          call get_var_int    ( vtsa0       ,'MVTSA0') ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          call get_var_real   ( lah         ,'MLAH')   ! LATITUDE (HIGH ACCURACY)
          call get_var_real   ( loh         ,'MLOH')   ! LONGITUDE (HIGH ACCURACY)
          call get_var_real   ( szan        ,'MSZAN')  ! SATELLITE ZENITH ANGLE
          call get_var_real   ( da          ,'MDA')    ! BEARING OR AZIMUTH
          call get_var_real   ( soza        ,'MSOZA')  ! SOLAR ZENITH ANGLE
          call get_var_real   ( das         ,'MDAS')   ! SOLAR AZIMUTH
          call get_var_real   ( ppp         ,'MPPP')   ! PRESSURE
          call get_var_real   ( tsts        ,'MTSTS')  ! SKIN TEMPERATURE
          call get_var_real   ( hp          ,'MHP')    ! HEIGHT OF STATION

          !-------------------------
          ! store in satem data type
          !-------------------------
          satem(j)% file        = source(i_file)% file    ! INPUT FILE NAME
          satem(j)% sat_id      = sat_id       ! SATELLITE IDENTIFIER
          satem(j)% ogc         = ogc          ! ORIGINATING/GENERATING CENTRE
          satem(j)% sain        = sain         ! SATELLITE INSTRUMENTS
          satem(j)% sacl        = sacl         ! SATELLITE CLASSIFIKATION
          satem(j)% orbnu       = orbnu        ! ORBIT NUMBER
          satem(j)% scaln       = scaln        ! SCAN LINE NUMBER
          satem(j)% fievn       = fievn        ! FIELD VIEW NUMBER
          satem(j)% sind        = sind         ! SUPERADIABATIC INDICATOR
          satem(j)% lase        = lase         ! LAND/SEA QUALIFIER
          satem(j)% dani        = dani         ! DAY/NIGHT QUALIFIER
          satem(j)% i4          = i4           ! SATELLITE DATA PROCESSING TECHNIQUE USED
          satem(j)% sgli        = sgli         ! SUN-GLINT INDICATOR
          satem(j)% vtsa        = vtsa         ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          satem(j)% noob        = noob         ! NUMBER OF OBSERVATIONS
          satem(j)% noob0       = noob0        ! NUMBER OF OBSERVATIONS
          satem(j)% vtsa0       = vtsa0        ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          satem(j)% lah         = lah          ! LATITUDE (HIGH ACCURACY)
          satem(j)% loh         = loh          ! LONGITUDE (HIGH ACCURACY)
          satem(j)% szan        = szan         ! SATELLITE ZENITH ANGLE
          satem(j)% da          = da           ! BEARING OR AZIMUTH
          satem(j)% soza        = soza         ! SOLAR ZENITH ANGLE
          satem(j)% das         = das          ! SOLAR AZIMUTH
          satem(j)% ppp         = ppp          ! PRESSURE
          satem(j)% tsts        = tsts         ! SKIN TEMPERATURE
          satem(j)% hp          = hp           ! HEIGHT OF STATION
          satem(j)% nlevel      = i_level_t    ! NO. OF LEVELS TEMPERATURE
          satem(j)% FLG_CLDFRM  = FLG_CLDFRM   ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
          !-----------------------
          ! cloud flag FLG_CLDFRM:
          !-----------------------
          !  bit  0    Cloudy but no height assignment possible
          !  bit  1    One cloud formation with ATOVS height assignment
          !  bit  2    Two cloud formations with ATOVS height assignment
          !  bit  3    Three cloud formations with ATOVS height assignment
          !  bit  4    One cloud formation with NWP forecast height assignment
          !  bit  5    Two cloud formations with NWP forecast height assignment
          !  bit  6    Three cloud formations with NWP forecast assignment
          !  bit  7    One cloud formation with climatological height assignment
          !  bit  8    Two cloud formations with climatological height assignment
          !  bit  9    Three cloud formations with climatological height assignment
          !  bit 10    The height assignment of ?rst cloud formation is ambiguous
          !  bit 11    The height assignment of second cloud formation is ambiguous
          !  bit 12    The height assignment of third cloud formation is ambiguous
          !  bit 13-15 not used
          satem(j)% FLG_CLDSUM  = FLG_CLDSUM   ! INSTRUMENT DETECTING CLOUDS
          satem(j)% FLG_IASIBAD = FLG_IASIBAD  ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
          satem(j)% FLG_040195  = FLG_040195   ! QUALITY AND COMPLETENESS OF RETRIEVAL
          satem(j)% FLG_040196  = FLG_040196   ! RETRIEVAL CHOICE INDICATOR
          satem(j)% FLG_040197  = FLG_040197   ! SATELLITE MANOEUVRE INDICATOR
          satem(j)% FLG_040198  = FLG_040198   ! SELECTION OF BACKGROUND STATE

          !-------------
          ! convert time
          !-------------
          call init_time (satem(j)% time, yyyy, mo, dd, hh, mi, ss)
       end do ! read satem profiles

       !------------------
       ! close NetCDF file
       !------------------
       status = nf90_close (ncid)
       if (status /= NF90_NOERR) call error ('nf90_close')
    enddo

  end subroutine read_IASI_l2_flags_netcdf
!========================================================================



  subroutine read_satem_netcdf (obs)
    type (t_obs)   ,intent(inout) :: obs            ! observations data type
    !------------------------------------------------------------------------
    ! IASI SATEM data are read from a NetCDF file and stored in the observation
    ! data type OBS.
    !------------------------------------------------------------------------
    integer       , parameter     :: mxvob = 100    ! max. no. of IASI SATEM layers
    !----------------
    ! local variables
    !----------------
    type(t_spot)                  :: spt            ! report meta data variable
    type(t_spot)                  :: empty          ! default initialised
    type(t_head)                  :: head           ! data usually stored in BUFR header
    type(t_use)                   :: use            ! state of the report
    type(t_satem)                 :: satem          ! IASI SATEM meta data
    type(t_satem_level), pointer  :: satem_level(:) =>NULL() ! IASI SATEM level data to store
    integer                       :: n              ! number of levels read
    integer                       :: i              ! index variable
    integer                       :: j              ! index variable
    integer                       :: i_file         ! index variable
    integer                       :: entry          ! reort ID derived from source/record

    !-------------------------------
    ! extra variables for necdf file
    !-------------------------------
    integer    :: nsatems     ! no. of SATEM profiles in input file
    integer    :: i_level_t   ! no. of levels temperature
    integer    :: i_level_m   ! no. of levels mixing ratio (rel. humidity)
    integer    :: sat_id      ! SATELLITE IDENTIFIER
    integer    :: yyyy        ! YEAR
    integer    :: mo          ! MONTH
    integer    :: dd          ! DAY
    integer    :: hh          ! HOUR
    integer    :: mi          ! MINUTE
    integer    :: ss          ! SECOND
    integer    :: ogc         ! ORIGINATING/GENERATING CENTRE
    integer    :: sain        ! SATELLITE INSTRUMENTS
    integer    :: sacl        ! SATELLITE CLASSIFIKATION
    integer    :: orbnu       ! ORBIT NUMBER
    integer    :: scaln       ! SCAN LINE NUMBER
    integer    :: fievn       ! FIELD VIEW NUMBER
    integer    :: sind        ! SUPERADIABATIC INDICATOR
    integer    :: lase        ! LAND/SEA QUALIFIER
    integer    :: dani        ! DAY/NIGHT QUALIFIE
    integer    :: i4          ! SATELLITE DATA PROCESSING TECHNIQUE USED
    integer    :: sgli        ! SUN-GLINTEGER INDICATOR
    integer    :: FLG_CLDFRM  ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
    integer    :: FLG_CLDSUM  ! INSTRUMENT DETECTING CLOUDS
    integer    :: FLG_IASIBAD ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
    integer    :: FLG_040195  ! QUALITY AND COMPLETENESS OF RETRIEVAL
    integer    :: FLG_040196  ! RETRIEVAL CHOICE INDICATOR
    integer    :: FLG_040197  ! SATELLITE MANOEUVRE INDICATOR
    integer    :: FLG_040198  ! SELECTION OF BACKGROUND STATE
    integer    :: vtsa        ! VERTICAL SIGNIFICANCE SATELL.OBSERV.)
    integer    :: noob        ! NUMBER OF OBSERVATIONS
    integer    :: noob0       ! NUMBER OF OBSERVATIONS
    integer    :: vtsa0       ! VERTICAL SIGNIFICANCE SATELL.OBSERV.)
    real(wp)   :: lah         ! LATITUDE (HIGH ACCURACY)
    real(wp)   :: loh         ! LONGITUDE (HIGH ACCURACY)
    real(wp)   :: szan        ! SATELLITE ZENITH ANGLE
    real(wp)   :: da          ! BEARING OR AZIMUTH
    real(wp)   :: soza        ! SOLAR ZENITH ANGLE
    real(wp)   :: das         ! SOLAR AZIMUTH
    real(wp)   :: ppp         ! PRESSURE
    real(wp)   :: tsts        ! SKIN TEMPERATURE
    real(wp)   :: hp          ! HEIGHT OF STATION
    real(wp)   :: pn (mxvob)  ! PRESSURE (VERT.LOCATION)
    real(wp)   :: tdbt(mxvob) ! TEMPERATURE/DRY BULB TEMPERATURE
    real(wp)   :: pn0(mxvob)  ! PRESSURE (VERT.LOCATION)
    real(wp)   :: mixr (mxvob)! MIXING RATIO (humidity)

    !--------------------------------------
    ! check for usage of SATEM
    !--------------------------------------
    if (rept_use(OT_SATEM)% use(CHK_NONE) <= STAT_DISMISS) return

    !-----------------------------------
    ! read from files (from one PE only)
    !-----------------------------------
    do i_file = 1, n_source

       !-------------------------
       ! skip inappropriate files
       !-------------------------
       if (source(i_file)% filetype /= FT_NETCDF) cycle
       if (source(i_file)% obstype  /= OT_SATEM)  cycle
       if (source(i_file)% pe       /= dace% pe)  cycle

       write(6,'(a,a)') ' reading ',source(i_file)% file

       !-----------------
       ! open NetCDF file
       !-----------------
       ncfile = path_file (source(i_file)% path, source(i_file)% file)
       status = nf90_open (ncfile, NF90_NOWRITE, ncid)
       if (status /= NF90_NOERR) call error ('nf90_open')

       !-------------------------------
       ! preset total number of reports
       !-------------------------------
       entry = sum (source(1:i_file-1)% entries)

       !---------------------------------
       ! read number of satem
       !  and number of levels
       !---------------------------------
       call get_dim     ( nsatems , 'BUFR_records')
       print * ,' satems read: ' , nsatems
       nsatems = min (nsatems, rept_use(OT_SATEM)% max_proc)
       print * ,' satems used: ' , nsatems
       call get_dim     ( i_level_t   , 'Loop_000_maxlen')
       call get_dim     ( i_level_m   , 'Loop_001_maxlen')
       print * ,' level2 temperature  levels   : ' , i_level_t
       print * ,' level2 mixing ratio levels   : ' , i_level_m

       if ((i_level_t > mxvob).or.(i_level_m > mxvob)) &
            call finish('read_satem_netcdf','iob > mxvob')

       ! initialize counter for clear spot profils:
       i=0

       !-------------------------
       ! read netcdf file, header
       !-------------------------
       ! loop over satem profiles
       do j=1,nsatems
          entry = entry + 1
          ! read header info of satem profile nsatems
          stanc  = j          ! read from satem profile j
          counc  = 1
          strnc  = 1

          call get_var_int    ( sat_id      ,'MI1I2')  ! SATELLITE IDENTIFIER
          call get_var_int    ( ogc         ,'MOGC')   ! ORIGINATING/GENERATING CENTRE
          call get_var_int    ( sain        ,'MSAIN')  ! SATELLITE INSTRUMENTS
          call get_var_int    ( sacl        ,'MSACL')  ! SATELLITE CLASSIFIKATION
          call get_var_int    ( yyyy        ,'MJJJ')   ! YEAR
          call get_var_int    ( mo          ,'MMM')    ! MONTH
          call get_var_int    ( dd          ,'MYY')    ! DAY
          call get_var_int    ( hh          ,'MGG')    ! HOUR
          call get_var_int    ( mi          ,'NGG')    ! MINUTE
          call get_var_int    ( ss          ,'MSEC')   ! SECOND
          call get_var_int    ( orbnu       ,'NORBNU') ! ORBIT NUMBER
          call get_var_int    ( scaln       ,'NSCALN') ! SCAN LINE NUMBER
          call get_var_int    ( fievn       ,'NFIEVN') ! FIELD VIEW NUMBER
          call get_var_int    ( sind        ,'MSIND')  ! SUPERADIABATIC INDICATOR
          call get_var_int    ( lase        ,'MLASE')  ! LAND/SEA QUALIFIER
          call get_var_int    ( dani        ,'MDANI')  ! DAY/NIGHT QUALIFIE
          call get_var_int    ( i4          ,'MI4')    ! SATELLITE DATA PROCESSING TECHNIQUE USED
          call get_var_int    ( sgli        ,'MSGLI')  ! SUN-GLINT INDICATOR
          call get_var_int    ( FLG_CLDFRM  ,'040192') ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
          call get_var_int    ( FLG_CLDSUM  ,'040193') ! INSTRUMENT DETECTING CLOUDS
          call get_var_int    ( FLG_IASIBAD ,'040194') ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
          call get_var_int    ( FLG_040195  ,'040195') ! QUALITY AND COMPLETENESS OF RETRIEVAL
          call get_var_int    ( FLG_040196  ,'040196') ! RETRIEVAL CHOICE INDICATOR
          call get_var_int    ( FLG_040197  ,'040197') ! SATELLITE MANOEUVRE INDICATOR
          call get_var_int    ( FLG_040198  ,'040198') ! SELECTION OF BACKGROUND STATE
          call get_var_int    ( vtsa        ,'MVTSA')  ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          call get_var_int    ( noob        ,'NNOOB')  ! NUMBER OF OBSERVATIONS
          call get_var_int    ( noob0       ,'NNOOB0') ! NUMBER OF OBSERVATIONS
          call get_var_int    ( vtsa0       ,'MVTSA0') ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          call get_var_real   ( lah         ,'MLAH')   ! LATITUDE (HIGH ACCURACY)
          call get_var_real   ( loh         ,'MLOH')   ! LONGITUDE (HIGH ACCURACY)
          call get_var_real   ( szan        ,'MSZAN')  ! SATELLITE ZENITH ANGLE
          call get_var_real   ( da          ,'MDA')    ! BEARING OR AZIMUTH
          call get_var_real   ( soza        ,'MSOZA')  ! SOLAR ZENITH ANGLE
          call get_var_real   ( das         ,'MDAS')   ! SOLAR AZIMUTH
          call get_var_real   ( ppp         ,'MPPP')   ! PRESSURE
          call get_var_real   ( tsts        ,'MTSTS')  ! SKIN TEMPERATURE
          call get_var_real   ( hp          ,'MHP')    ! HEIGHT OF STATION

          !-------------------------
          ! store in satem data type
          !-------------------------
          satem% file        = source(i_file)% file    ! INPUT FILE NAME
          satem% sat_id      = sat_id       ! SATELLITE IDENTIFIER
          satem% ogc         = ogc          ! ORIGINATING/GENERATING CENTRE
          satem% sain        = sain         ! SATELLITE INSTRUMENTS
          satem% sacl        = sacl         ! SATELLITE CLASSIFIKATION
          satem% orbnu       = orbnu        ! ORBIT NUMBER
          satem% scaln       = scaln        ! SCAN LINE NUMBER
          satem% fievn       = fievn        ! FIELD VIEW NUMBER
          satem% sind        = sind         ! SUPERADIABATIC INDICATOR
          satem% lase        = lase         ! LAND/SEA QUALIFIER
          satem% dani        = dani         ! DAY/NIGHT QUALIFIER
          satem% i4          = i4           ! SATELLITE DATA PROCESSING TECHNIQUE USED
          satem% sgli        = sgli         ! SUN-GLINT INDICATOR
          satem% vtsa        = vtsa         ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          satem% noob        = noob         ! NUMBER OF OBSERVATIONS
          satem% noob0       = noob0        ! NUMBER OF OBSERVATIONS
          satem% vtsa0       = vtsa0        ! VERTICAL SIGNIFICANCE (SATELL.OBSERV.)
          satem% lah         = lah          ! LATITUDE (HIGH ACCURACY)
          satem% loh         = loh          ! LONGITUDE (HIGH ACCURACY)
          satem% szan        = szan         ! SATELLITE ZENITH ANGLE
          satem% da          = da           ! BEARING OR AZIMUTH
          satem% soza        = soza         ! SOLAR ZENITH ANGLE
          satem% das         = das          ! SOLAR AZIMUTH
          satem% ppp         = ppp          ! PRESSURE
          satem% tsts        = tsts         ! SKIN TEMPERATURE
          satem% hp          = hp           ! HEIGHT OF STATION
          satem% nlevel      = i_level_t    ! NO. OF LEVELS TEMPERATURE
          satem% FLG_CLDFRM  = FLG_CLDFRM   ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
          satem% FLG_CLDSUM  = FLG_CLDSUM   ! INSTRUMENT DETECTING CLOUDS
          satem% FLG_IASIBAD = FLG_IASIBAD  ! VALIDATION FLAG FOR IASI LEVEL 1 PRODUCT
          satem% FLG_040195  = FLG_040195   ! QUALITY AND COMPLETENESS OF RETRIEVAL
          satem% FLG_040196  = FLG_040196   ! RETRIEVAL CHOICE INDICATOR
          satem% FLG_040197  = FLG_040197   ! SATELLITE MANOEUVRE INDICATOR
          satem% FLG_040198  = FLG_040198   ! SELECTION OF BACKGROUND STATE

          !----------------------------------
          !!!!! process clear spots only !!!!
          !----------------------------------
          if (satem% FLG_CLDFRM > 0) cycle   ! CLOUD FORMATION AND HEIGHT ASSIGNMENT

          !-------------
          ! convert time
          !-------------
          call init_time (satem% time, yyyy, mo, dd, hh, mi, ss)

          !---------------------------------------------------
          ! store header data, insert DWD dbkz, perform checks
          !---------------------------------------------------
          head% modtype     = SATEM_obstype ! observation operator type
          head% obstype     = OT_SATEM      ! observation operator type
          head% buf_type    = 3             ! BUFR type
          head% buf_subtype = 223           ! BUFR subtype
          head% source      = i_file        ! source file number
          head% record      = j             ! record in file
          head% time        = satem% time   ! date/time stamp
          head% dbkz        = 1794          ! DBKZ/internal ID
          head% id          = entry         ! ID derived from source/record
!         !------------
!         ! not yet set
!         !------------
!         head% codetype    = obsid% codetype
!         head% db_time     = db_time(is)
!         head% idbk        = report_subti
!         head% center      = s1cent(is)
!         head% subcenter   = s1cents(is)

          call check_report_0 (use, head,1)
          if (use% state <= STAT_DISMISS) cycle

          !------------------------------------------
          ! for diagnistics:
          ! select clear sky profile
          ! no. [n_clear_satem_b ... n_clear_satem_e]
          !------------------------------------------
          i=i+1
          if (i < n_clear_satem_b) cycle
          if (i > n_clear_satem_e) exit

          !--------------------------------------
          ! read level-2 data (profiles)
          ! store to t_satem_level data structure
          !--------------------------------------
          n = i_level_t
          if(associated(satem_level)) deallocate(satem_level)
          allocate(satem_level(n))
          stanc  = (/1,j,0/)                      ! read from satem j
          counc  = (/n,1,0/)                      ! n data
          strnc  = (/1,1,0/)                      !
          call get_var_real_1   ( pn  , 'MPN')    ! PRESSURE (VERT.LOCATION)
          call get_var_real_1   ( tdbt, 'MTDBT')  ! TEMPERATURE/DRY BULB TEMPERATURE
          call get_var_real_1   ( pn0 , 'MPN0')   ! PRESSURE (VERT.LOCATION)
          call get_var_real_1   ( mixr, 'MMIXR')  ! MIXING RATIO (REL. HUMIDITY)
          satem_level% pn    = pn   (1:n)
          satem_level% tdbt  = tdbt (1:n)
          satem_level% pn0   = pn0  (1:n)
          satem_level% mixr  = mixr (1:n)

          !-------------------------------------------
          ! Final preparation of observation type data
          ! Storage into components of 'obs' and 'spot'
          !-------------------------------------------
          spt            = empty
          spt% use       = use
          spt% hd        = head
          spt% hd% satid = sat_id
          spt% ident     = sat_id
          spt% statid    = mnem(sat_id)
          call check_store_satem (satem, satem_level, spt, obs)
          call flush_buf_pio
          deallocate (satem_level)
       end do ! read satem profiles

       !------------------
       ! close NetCDF file
       !------------------
       status = nf90_close (ncid)
       if (status /= NF90_NOERR) call error ('nf90_close')
    enddo

  end subroutine read_satem_netcdf
!============================================================================
! Auxiliary routines to read NetCDF files follow
!----------------------------------------------------------------------------
  subroutine get_dim (dim ,name, stat)
  !-------------------------------
  ! get dimension from Netcdf file
  !-------------------------------
  integer       ,intent(out) :: dim
  character(len=*)  ,intent(in)  :: name
  integer, optional ,intent(out) :: stat
  integer :: dimid
  status = nf90_inq_dimid   (ncid, name, dimid)
  if (present (stat)) stat = status
  if (status /= NF90_NOERR) then
    if (present (stat)) return
    call error ('nf90_inq_dimid ('//name//')')
  endif
  status = nf90_Inquire_Dimension (ncid, dimid, len=dim)
  if (status /= NF90_NOERR) call error ('nf90_inq_dimlen ('//name//')')
  end subroutine get_dim
!----------------------------------------------------------------------------
  subroutine get_var_int (int ,name)
  !---------------------
  ! get 0-D int variable
  !---------------------
  integer      ,intent(out) :: int
  character(len=*) ,intent(in)  :: name
  integer :: varid
  status = nf90_inq_varid (ncid, name, varid)
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, int, stanc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::int ('//name//')')
  end subroutine get_var_int
!----------------------------------------------------------------------------
  subroutine get_var_real (x ,name)
  !---------------------
  ! get 0-D real variable
  !---------------------
  real(dp)         ,intent(out) :: x
  character(len=*) ,intent(in)  :: name
  integer :: varid
  status = nf90_inq_varid (ncid, name, varid)
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, x, stanc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::double ('//name//')')
  end subroutine get_var_real
!----------------------------------------------------------------------------
  subroutine get_var_real_1 (x ,name)
  !----------------------
  ! get 1-D real variable
  !----------------------
  real(dp)         ,intent(out) :: x (:)
  character(len=*) ,intent(in)  :: name
  integer  :: varid
  status = nf90_inq_varid (ncid, name, varid)    ! get varid
  if (status /= NF90_NOERR) call error ('nf90_inq_varid ('//name//')')
  status = nf90_get_var (ncid, varid, x, stanc, counc, strnc)
  if (status /= NF90_NOERR) call error ('nf90_get_var::double ('//name//')')
  end subroutine get_var_real_1
!----------------------------------------------------------------------------
  subroutine error (string)
  !-------------------------
  ! abort on error condition
  !-------------------------
  character(len=*) string
  character(len=80) :: str
  str = nf90_strerror (status)
  write (0,*) 'read_satem_netcdf: '//trim(string)//', file='//trim(ncfile)
  write (0,*) trim(str)
  write (0,*) 'start  =',stanc
  write (0,*) 'count  =',counc
  write (0,*) 'stride =',strnc
  call finish ('read_satem_netcdf',trim(string)//', file='//trim(ncfile))
  end subroutine error
!==============================================================================
  subroutine process_satem (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                            state)
    integer            ,intent(in)             :: task    ! what to do
    type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
    type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
    type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
    type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
    type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
    type(t_vector_segm),intent(inout),optional :: y       ! observed quantity
    real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
    type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
    integer            ,intent(in)   ,optional :: state   ! status flag

    !----------------------------------------------------------------------------
    ! Main observation processing routine
    ! similar for all observation operators
    !----------------------------------------------------------------------------
    !----------------
    ! local variables
    !----------------
    type(t_satem_level) ,pointer :: satem_levels(:) ! satem level data type
    type(t_satem)                :: satem           ! satem data type
    integer                      :: ii, io          ! indices
    integer                      :: i,j,k,l,n,i2    ! indices
    integer                      :: tsk             ! task, local copy
    logical                      :: change          ! true for observations changed
    integer                      :: jtv, jrh, jh    ! indices
    real(wp)                     :: p               ! pressure
    real(wp)                     :: q               ! specific humidity
    real(wp)                     :: t               ! temperature
    real(wp)                     :: gh              ! generalized humidity
    real(wp)                     :: tv              ! virtual temperature
    real(wp)                     :: h               ! geopotential
    real(wp)                     :: dt_tv, dt_gh    ! gradients
    real(wp)                     :: dq_tv, dq_gh    ! gradients
    integer                      :: ifail           ! error code
    integer                      :: ind_levs_i      ! indices to locate layer
    integer                      :: ind_levs_j      ! indices to locate layer
    real(wp),allocatable         :: Hnew(:)
    real(wp),allocatable         :: yn  (:)
    real(wp),allocatable         :: R(:,:)          ! layer thickness: R explicit
    real(wp),allocatable         :: Rnew(:)         !                  R packed
    real(wp),allocatable         :: gp_cov(:)       !                  covariances

    nullify (satem_levels)

    !==============================
    ! observation non_specific part
    !==============================
    !----------------
    ! tsk == TSK_INIT
    ! read namelist
    !----------------
    tsk = task
    if (iand (TSK_INIT,tsk) /= 0) then
       call read_satem_nml
       tsk=tsk-TSK_INIT
    endif
    if (tsk==0) return

    !-----------------
    ! tsk == TSK_READ:
    ! read data
    !-----------------
    if (iand (TSK_READ,tsk) /= 0) then
       call read_satem_netcdf (obs% o)
       tsk=tsk-TSK_READ
    endif
    if (tsk==0) return

    !---------------------------------------------
    ! TSK_SET_CHR: set observation characteristics
    !---------------------------------------------
    if (iand (TSK_SET_CHR,tsk) /= 0) then
       tsk=tsk-TSK_SET_CHR
    endif
    if (tsk==0) return

    !==========================
    ! observation specific part
    !==========================
    if (.not.present(spot)) call finish('process_satem','spot not present')
    if (spot% hd% modtype /= SATEM_obstype) return

    !==================================================================
    ! tsk == TSK_SHRINK:
    ! release unused observations in the report
    !==================================================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
       call shrink_report (spot, obs%o, state, change)
       tsk = tsk - TSK_SHRINK
       if (tsk == 0) return
    endif

    !========================================================
    ! setup model columns
    !========================================================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then
       if (layer_thickness) then
          call idx_init (                          &
               spot% col% c,                       &! <-  column descriptor
               spot% col% h,                       &!  -> interpolation coefs.
               obs% o% mc,                         &! <-> model column descriptors
               COL_GEOH+COL_PH,                    &! <-  fields required
               0,                                  &! <-  tracers required
               atm% grid,                          &! <-  model grid
               spot% i_time,                       &! <-  time slot
               spot% w_time                        )! <-  time interpolation weight
       else
          call idx_init (                          &
               spot% col% c,                       &! <-  column descriptor
               spot% col% h,                       &!  -> interpolation coefs.
               obs% o% mc,                         &! <-> model column descriptors
               COL_TV+COL_RH,                      &! <-  fields required
               0,                                  &! <-  tracers required
               atm% grid,                          &! <-  model grid
               spot% i_time,                       &! <-  time slot
               spot% w_time                        )! <-  time interpolation weight
       endif
       if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)
       tsk = tsk - TSK_SETUP_COLS

    endif
    if(tsk==0) goto 888

    !========================================================
    ! final setup of PSAS interpolation space
    !========================================================
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
       !----------------------------
       ! request interpolation space
       !----------------------------
       if (layer_thickness) then
!         call new_int (obs% o, spot, n_layer+1)         ! for layer thickness:
          obs% o% lev  (spot%i%i+1:spot%i%i+spot%i%n)= & ! different sizes of
               log(satem_levs(1:n_layer+1))              ! observation and
          obs% o% t_int(spot%i%i+1:spot%i%i+spot%i%n)= & ! interpolation space
               OBS_H
       endif
       tsk = tsk - TSK_SETUP_FUL0

    endif

    !========================================================
    ! set interpolation space: observed quantities and levels
    !========================================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
       tsk=tsk-TSK_SETUP_FULL
    endif
    if(tsk==0) goto 888

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
       if (dace% pe==obs% o% pe) then

          !-----------------------------------------------------------------------
          ! set observational error
          ! (for layer thickness only. Observational errors
          ! for t and q have been specified via /SATEM/ namelist).
          ! The layer thickness errors are set according to ECMWF
          ! research manual 1 "ECMWF data assimilation - Scientific documentation"
          !-----------------------------------------------------------------------
          if (any(obs% o% varno (spot% o% i+1:spot% o% i+ spot% o% n) == VN_Z)) then
             if (layer_thickness) then
                allocate (R(n_layer,n_layer))
                allocate (gp_cov(n_layer))
                allocate (Rnew(n_layer*n_layer))
                R =  0._wp
                do i=1,n_layer
                   R(i,i)=1._wp
                enddo
                if (spot% stclf == 0) then

                   !-------------
                   ! clear spot
                   !-------------
                   gp_cov(1)= 27._wp; gp_cov(2)= 17._wp; gp_cov(3)= 25._wp ! covariances
                   gp_cov(4)= 53._wp; gp_cov(5)= 33._wp; gp_cov(6)= 36._wp ! [gpm]
                   gp_cov(7)= 74._wp
                   R(2,1) =  0.22_wp; R(3,1) = -0.18_wp; R(3,2) =  0.27_wp ! correlations
                   R(4,3) = -0.16_wp; R(5,3) =  0.16_wp; R(6,4) =  0.15_wp !
                   R(7,6) =  0.27_wp                                       !
                else

                   !-------------
                   ! cloudy spot
                   !-------------
                   gp_cov(1)= 33._wp; gp_cov(2)= 20._wp; gp_cov(3)= 31._wp ! covariances
                   gp_cov(4)= 56._wp; gp_cov(5)= 33._wp; gp_cov(6)= 36._wp ! [gpm]
                   gp_cov(7)= 74._wp
                   R(2,1) =  0.22_wp; R(3,1)= -0.34_wp; R(3,2)=  0.20_wp   ! correlations
                   R(4,2) = -0.25_wp; R(4,3)= -0.29_wp; R(6,4)=  0.20_wp   !
                endif
                do i = 1, n_layer
                   R(i,:) = R(i,:) * gp_cov(i)
                   R(:,i) = R(:,i) * gp_cov(i)
                enddo
                do i = 1, n_layer
                   do j = i + 1, n_layer
                      R(i,j) = R(j,i)
                   enddo
                enddo
                !----------------------------------------------------
                ! store R matrix in 1d vector (sparse representation)
                !----------------------------------------------------
                do n = 1, n_layer
                   Rnew((n-1)*n_layer+1:n*n_layer) = R(:,n)
                enddo
                deallocate (R, gp_cov)
             endif
          endif

          !---------------------------------
          ! store R in sparse representation
          !---------------------------------
          if (any(obs%o%varno(spot%o%i+1:spot%o%i+spot%o%n)==VN_T) .or.  &    ! load additional
              any(obs%o%varno(spot%o%i+1:spot%o%i+spot%o%n)==VN_Q))      &    ! SATEM info if
              call load_satem (obs% o, spot, satem, satem_levels)              ! needed (t,q)
          k = obs% R% ia (spot% o% i+1)
          do j=1,spot% o% n
             obs% R% ia (spot% o% i + j) = k
             if (obs% o% varno (spot% o% i+j) == VN_T) then                   ! temperature:
                obs% R% ja     (k) = spot% o% i + j                          ! row index
                obs% R% packed (k) = satem_levels(nint(j/2._wp))% var_t        !
                k = k + 1                                                      !
             endif
             if (obs% o% varno (spot% o% i+j) == VN_Q) then                   ! spec. humidity:
                obs% R% ja     (k) = spot% o% i + j                                         ! row index
                obs% R% packed (k) = satem_levels(nint(j/2._wp))% var_q        !
                k = k + 1                                                      !
             endif
             if (obs% o% varno (spot% o% i+j) == VN_Z) then                   ! thickness layers:
                ind_levs_j =                                                  &! index to identify
                     sum(minloc(abs(satem_levs_mid-obs%o%olev(spot%o%i+j))))   ! layer of observation
                do i=1,spot% o% n                                              !
                   !++ sx8 workaround: ++++
                   i2=spot%o%i+i
                   !+++++++++++++++++++++++
                   ind_levs_i =                                               &! index to identify
                        sum(minloc(abs(satem_levs_mid-obs%o%olev(i2))))        ! layer of observation
                   l = (ind_levs_j - 1) * n_layer + ind_levs_i                 !
                   if (Rnew(l) /= 0._wp) then                                  !
                      obs% R% packed (k) = Rnew(l)                             !
                       obs% R% ja (k) =  i2                                     ! row index
                      k = k + 1                                                !
                   endif                                                       !
                end do                                                         !
             endif                                                             !
          end do
          obs% R% ia (spot% o% i + spot% o% n + 1) = k                         ! column index
          if (associated(satem_levels)) deallocate (satem_levels)

          !-----------------
          ! setup vqc bounds
          !-----------------
          call set_vqc_insitu(spot, obs% o)
       endif
       tsk = tsk - TSK_R
    endif
    if(tsk==0) goto 888

    !==================
    ! set up H operator
    !==================
    if (iand (TSK_K,tsk) /= 0) then
       if (spot% pe_eval == dace% pe) then
          n = spot% o% n
          allocate (Hnew(spot% o% n * spot% i% n))
          allocate (yn  (spot% o% n))
          Hnew = 0._wp
          j  = 1
          ii = spot%i% i + 1
          p  = -huge(p)
          do i = 1, spot% o% n
             io = spot%o% i + i
             if (p /= obs% o% olev(io)) then
                p  = obs% o% olev(io)
                tv = -1._wp
                h  = -1._wp
             endif
             select case (obs% o% varno(io))
             case (VN_T, VN_Q)

                !------------------------------
                ! determine fg values t, q from
                ! 3D-Var parameters tv, gh,
                ! + gradients
                !------------------------------
                if (tv < 0._wp) then
                   jtv = j
                   jrh = j+1
                   tv  = xi%x (ii)
                   gh  = xi%x (ii+1)
                   call tq_tvgh (                     &!
                        t,                            &! ->  t
                        q,                            &! ->  q
                        tv,                           &! <-  tv
                        gh,                           &! <-  gh
                        p,                            &! <-  p
                        ifail,                        &! -> error code
                        dt_tv, dt_gh, dq_tv, dq_gh)    ! -> gradient
                   if (ifail < 0) &
                        call finish ('process_satem (TSK_K)','tq_tvgh failed')
                   ii = ii + 2
                   j  = j  + 2
                endif
                if (obs% o% varno(io) == VN_T) then
                   yn(i) = t
                   if (present(y)) y%x (io) = t
                   Hnew ((jtv - 1) * n + i) = dt_tv
                   Hnew ((jrh - 1) * n + i) = dt_gh
                else if (obs% o% varno(io) == VN_Q) then
                   yn(i) = q
                   if (present(y)) y%x (io) = q
                   Hnew ((jtv - 1) * n + i) = dq_tv
                   Hnew ((jrh - 1) * n + i) = dq_gh
                endif
             case (VN_Z)

                !----------------------------------------------------
                ! determine layer thickness values as corresponding
                ! differences of model geopotential background.
                ! The H matrix is then
                !     / -1  1  0  0  0 ...          \   \
                !     |  0 -1  1  0  0 ...          |   |
                !     |  0  0 -1  1  0 ...          |   obs. space:
                ! H =                  ...              size = n_obs
                !     |                ... -1  1  0 |   |
                !     \                ...  0 -1  1 /   /
                !
                !     \---- interpol. space:  ------/
                !           size = n_layer+1
                !----------------------------------------------------
                if (h < 0._wp) then
                   ii = spot%i% i +                                        & ! indices to identify
                        sum(minloc(abs(obs%o%olev(io)-satem_levs_mid)))      ! layer of observation
                   jh = j +                                                & !
                        sum(minloc(abs(obs%o%olev(io)-satem_levs_mid)))-1    !
                   h  = xi%x (ii+1) - xi%x (ii)
                endif
                if (present(y)) y%x (io) = h
                Hnew ((jh - 1) * n + i ) = -1._wp      ! entries of H
                Hnew ((jh    ) * n + i ) =  1._wp      !
                yn(i) = xi%x (ii+1) - xi%x (ii)        ! fg layer thickness
             case default
                write (0,*) 'process_satem, TSK_K: invalid obstype',&
                     obs% o% varno(io)
                call finish('process_satem, TSK_K','invalid observation type')
             end select
          end do

          !----------------------------
          ! store H matrix, first guess
          !----------------------------
          k = obs% H% ia (spot% i% i + 1)
          do j=1,spot% i% n                              ! columns
             obs% H% ia (spot% i% i +j) = k              ! column index
             do i=1,spot% o% n                           ! rows
                l = (j - 1) * n + i
                if (Hnew(l) /= 0._wp) then
                   obs% H% packed (k) = Hnew(l)          ! H entries
                   obs% H% ja (k) = spot% o% i + i       ! row index
                   k = k + 1
                endif
             end do
          end do
          obs% H% ia (spot% i% i + spot% i% n + 1) = k       ! column index
!         obs% xi% x (spot% i% i+1:spot% i% i+spot% i% n) = &! fg
!              xi% x (spot% i% i+1:spot% i% i+spot% i% n)
          obs% yi% x (spot% o% i+1:spot% o% i+spot% o% n) = yn (1:spot%o% n)
          deallocate (Hnew)
          deallocate (yn)
       endif
       tsk = tsk - TSK_K
       if (tsk == 0) goto 888
    endif

!   !=============
!   ! unknown task
!   !=============
!   if (tsk /= 0) then
!      write (0,*)  'process_satem:  unknown task',tsk
!      call finish ('process_satem','unknown task')
!   endif

    !========================
    ! release local variables
    !========================
888 continue

  end subroutine process_satem
!==============================================================================
!==============================================================================
  !---------------------------
  ! Private auxiliary routines
  !---------------------------
!------------------------------------------------------------------------------
  subroutine set_size
  !------------------------------------------------------------------
  ! store sizes of derived data types T_SATEM and T_SATEM_LEVELS
  ! (in termes of size of component OBS% PAR) into
  ! private module variables SATEM_INT_SIZE and SATEM_LEVELS_INT_SIZE
  !------------------------------------------------------------------
    type (t_satem)         :: satem
    type (t_satem_level)   :: satem_levels
    type (t_obs)  :: obs
    if (satem_int_size == 0) then
      satem_int_size         = size (transfer (satem, obs% par))
      satem_levels_int_size  = size (transfer (satem_levels, obs% par))
      satem_byte_size        = size (transfer (satem, (/' '/) ))
      satem_levels_byte_size = size (transfer (satem_levels, (/' '/) ))
    endif
  end subroutine set_size
!------------------------------------------------------------------------------
  subroutine store_satem (obs, spot, satem, satem_levels)
  type (t_obs)          ,intent(inout) :: obs              ! data of all observations
  type (t_spot)         ,intent(inout) :: spot             ! meta data of this observation
  type (t_satem)        ,intent(in)    :: satem            ! satem meta data
  type (t_satem_level)  ,intent(in)    :: satem_levels (:) ! satem data
  !-----------------------------------------------------------------------------
  ! Store the data from variables SATEM and SATEM_LEVELS in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------------
    integer ,pointer :: par (:)
    integer          :: n, m
    if (satem_int_size == 0) call set_size
    n = spot%col%nlev
    m = satem_int_size + n * satem_levels_int_size
    if (spot% p% i < 0) call new_par (obs, m, spot=spot)
    if (m < spot% p% n) spot% p% n = m
    if (m > spot% p% n) call finish('store_satem','m > spot% p% n')
    par => obs % par (spot% p% i+1 : spot% p% i + spot% p% n)
    par (1 :satem_int_size)  = transfer(satem           ,par)
    par (satem_int_size+1 :) = transfer(satem_levels(:n),par)
  end subroutine store_satem
!------------------------------------------------------------------------------
  subroutine load_satem (obs, spot, satem, satem_levels)
  type (t_obs)          ,intent(in)  :: obs              ! data of all observations
  type (t_spot)         ,intent(in)  :: spot             ! meta data of this observation
  type (t_satem)        ,intent(out) :: satem            ! satem meta data
  type (t_satem_level)  ,pointer     :: satem_levels (:) ! satem data
  !------------------------------------------------------------------
  ! Load the data from component PAR of OBS from position provided by
  ! SPOT. Store into SATEM and SATEM_LEVELS.
  ! allocate SATEM_LEVELS with size required.
  !------------------------------------------------------------------
    allocate (satem_levels (spot% col% nlev))
    if (satem_int_size == 0) call set_size
    satem  = transfer (obs% par (spot% p% i+1 : spot% p% i + satem_int_size), satem)
    satem_levels = transfer (obs% par (spot% p% i+1 + satem_int_size : &
                               spot% p% i   + spot% p% n), satem_levels)
  end subroutine load_satem
!==============================================================================
  subroutine read_satem_nml
  !---------------------------
  ! read namelist /IASI_SATEM/
  !---------------------------
    integer          :: ierr                 ! error indicator
    integer          :: i                    ! index
    integer          :: ind_levs(nmxl)       ! indices to pick up
    integer          :: n_levs  = 0          ! no. of levels to confine satem layers
    real(wp),pointer :: old_levs(:)          ! levels to confine satem layers (temp.)
    real(wp),pointer :: satem_levs_tmp(:)    ! levels to confine satem layers (temp.)

    !-------------
    ! set defaults
    !-------------
    verbose              = 0                 ! verbosity: be chatty for verbose>0
    smooth_satem         = .false.           ! smooth initial t, q retrieval
    satem_init_levs( : ) = (/100000._wp,  &  ! SATEM layer height pressure levels [Pa]
                              70000._wp,  &  ! according to ECMWF research manual 1
                              50000._wp,  &  ! "ECMWF data assimilation -
                              30000._wp,  &  !  Scientific documentation"
                              10000._wp,  &  !
                               5000._wp,  &  !
                               3000._wp,  &  !
                               1000._wp,  &  !
                                  0._wp,  &  !
                                  0._wp   /) !
    obs_err_abs_t        = 1._wp             ! nominal absolute obs. error of t
    obs_err_rel_q        = 1._wp             ! nominal relative obs. error of q
    layer_thickness      = .true.            ! true: consider layer thickness [gpm] as
                                             ! observation data instead of t, q
    n_clear_satem_b      = -100000000        ! first clear spot to be processed
    n_clear_satem_e      =  100000000        ! last  "

    obs_layer_method     = 2                 ! method to derive geopotential from retrieved
                                             ! temperature: 1 (spline integration)
                                             !              2 (linear integration)
                                             !              3 (staircase t averaging)
    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
       call position_nml ('IASI_SATEM', status=ierr)
       select case (ierr)
       case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=IASI_SATEM, iostat=ierr)
          if (ierr/=0) call finish ('read_satem_nml','ERROR in namelist /IASI_SATEM/')
#else
          read (nnml ,nml=IASI_SATEM)
#endif
       end select
    endif

    !-----------------------------------------
    ! set satem layer parameters (number,
    ! layer boundary levels, layer midpoints),
    ! broadcast layer parameters
    !-----------------------------------------
    if (dace% lpio) then
       allocate(satem_levs_tmp(size(satem_init_levs)))
       satem_levs_tmp = satem_init_levs
       n_levs = 0
       do i = 1, size(satem_levs_tmp)
          if (satem_levs_tmp(i) <= 0._wp)  cycle
          n_levs = n_levs + 1
          ind_levs(n_levs) = i
       enddo
       !------------------------------
       ! finish if less than 2 layer
       ! boundary levels are specified
       !------------------------------
       if (n_levs<2)                      &
            call finish ('read_satem_nml',&
            'ERROR: specify at least one lower and upper bound for thickness layer')
       old_levs  => satem_levs_tmp
       n_layer = n_levs - 1
    endif
    call p_bcast (n_layer,dace% pio)
    allocate (satem_levs (n_layer+1))
    if (dace% lpio) then
       satem_levs = old_levs(ind_levs(1:n_layer+1))
    endif
    call p_bcast (satem_levs, dace% pio)
    allocate (satem_levs_mid(n_layer))
    do i = 1, n_layer
       satem_levs_mid(i) = &
            (satem_levs(i) + satem_levs(i+1))/2._wp
    enddo

    !-------------------------------------
    ! broadcast other namelist parameters
    !-------------------------------------
    call p_bcast (smooth_satem     ,dace% pio)
    call p_bcast (verbose          ,dace% pio)
    call p_bcast (layer_thickness  ,dace% pio)
    call p_bcast (obs_err_abs_t    ,dace% pio)
    call p_bcast (obs_err_rel_q    ,dace% pio)
    call p_bcast (n_clear_satem_b  ,dace% pio)
    call p_bcast (n_clear_satem_e  ,dace% pio)
    call p_bcast (obs_layer_method ,dace% pio)

  end subroutine read_satem_nml
!==============================================================================
  !-------------------------------------------------
  ! Called from routine reading BUFR or NetCDF files
  !   Final preparation of observation type data
  !   Storage into components of 'obs'
  !-------------------------------------------------
  subroutine check_store_satem (satem, satem_levels, spot, obs)
    integer             ,parameter     :: mxvob = 100       ! max. no. levels per satemprofile
    type(t_satem)       ,intent(inout) :: satem             ! obs. data of satem profile
    type(t_satem_level) ,pointer       :: satem_levels (:)  ! obs. data on satem level
    type(t_spot)        ,intent(inout) :: spot              ! spot data type
    type(t_obs)         ,intent(inout) :: obs               ! observation data type
    integer                      :: ind(mxvob)              ! indices to pick up
    integer                      :: k, i, j                 ! indices
    integer                      :: nlevel                  ! number of levels used
    type(t_satem_level) ,pointer :: old (:)                 ! temporary
    type(t_spot)        ,pointer :: spt                     ! temporary
    real(wp)                     :: w, ww                   ! temporary weight
    integer                      :: k1,k2,k3                ! indices
    integer                      :: id                      ! observation id
    type(t_datum)                :: bod                     ! datum
    real(wp)                     :: pn    = invalid         ! temp. pressure
    real(wp)                     :: tdbt  = invalid         ! temp. temperature
    real(wp)                     :: tdbt2 = invalid         ! temp. temperature var.
    real(wp)                     :: mixr  = invalid         ! temp. rel. humidity
    real(wp)                     :: mixr2 = invalid         ! temp.  rel. humidity var.
    real(wp)                     :: q                       ! specific humidity
    real(wp)                     :: t                       ! temperature
    integer                      :: ind_layers(nmxl)        ! indices to pick up
    integer                      :: nobs_layers (nmxl) = 0  ! no. of obs in levels
    integer                      :: i_int, n_int            ! indices of spline integration
    real(wp)                     :: geop_i                  ! tmp. geop. of layer fraction
    real(wp) ,pointer            :: int_geop(:)             !      "          "     "
    real(wp)                     :: l_th_j                  !      "          "     "
    integer                      :: it                      ! no. of steps in spline integr.
    real(wp)                     :: del, sum, tnm, x        ! temporaries in spline integr.
    real(wp) ,pointer            :: rho_inv(:), rho_inv_2(:)! tmp. rho^-1 vector
    real(wp)                     :: rho_inv_i, rho_inv_i1   !   "    "    values
    real(wp)                     :: rho_inv_x               !   "    "    values
    real(wp) ,pointer            :: l_th(:)                 ! tmp. layer thickness
    real(wp)                     :: w_top, w_bot            ! tmp. wheights in geop. integration
    integer                      :: i_top, i_bot            ! tmp. indices    "         "
    integer                      :: l_top, l_bot            ! top, bot. levels of t retrieval
                                                            ! profile to be used for SATEM
    nlevel = 0

    !----------------------------------------
    ! sanity check for data in SATEM profile,
    ! count valid level, condense array
    !----------------------------------------
    do k = 1, satem% nlevel
       if (( satem_levels(k)% tdbt <= 0._wp     ) .or. &
            (satem_levels(k)% tdbt >= 1000._wp  ) .or. &
            (satem_levels(k)% pn   == invalid)) cycle
       nlevel = nlevel + 1
       ind(nlevel) = k
    enddo
    old  => satem_levels
    allocate (satem_levels (nlevel))
    satem% nlevel = nlevel
    satem_levels = old (ind(1:nlevel))

    !-----------------------------------------------
    ! specify observation errors for t, q (namelist)
    ! (layer thickness obs. errors are specified in
    ! subroutine process_satem)
    !-----------------------------------------------
    do k=1, satem% nlevel
       satem_levels(k)% var_t =                                   & ! temperature
            max (0._wp, obs_err_abs_t**2)                           !
       satem_levels(k)% var_q =                                   & ! rel. humidity
            max (0._wp, (obs_err_rel_q*satem_levels(k)% mixr)**2)   !
    enddo

    !-----------------------------------------------------------------------------
    ! keep selected levels and nonempty spots.
    ! First assignement of satem% nlevel as no. of t,q-levels in
    ! retrieved IASI profile. When layer thickness is used as observation,
    ! satem% nlevel is assigned the number of thickness layers lateron (see below)
    !-----------------------------------------------------------------------------
    i = 0
    if (smooth_satem) then
       nlevel = nlevel - 2
       i = 1
    endif
    if (nlevel <= 0) then
       call decr_rpt_use (spot, CHK_INSDAT, comment='nlevel <= 0')
       return
    endif
    deallocate (old)
    old => satem_levels
    allocate (satem_levels (nlevel))
    satem% nlevel = nlevel

!   satem_levels = old (ind(1+i:satem% nlevel-i)) !? ind refers to previous'old'
    satem_levels = old (1+i:nlevel+i)             !!
    ind = (/(i,i=1,mxvob)/)                       !! fix ind for subsequent use

    !-----------------------------------------------
    ! smooth observational SATEM temperatures (tdbt)
    !-----------------------------------------------
    if (smooth_satem) then
       satem_levels(1)% tdbt2 = 0._wp
       do i = 2, satem% nlevel + 1
          ww   = 0._wp
          tdbt  = 0._wp
          tdbt2 = 0._wp
          mixr  = 0._wp
          mixr2 = 0._wp
          pn    = 0._wp
          k1 = ind (i-1)
          k2 = ind (i)
          k3 = ind (i+1)
          do k = k1, k2
             w     = (old(k)%pn - old(k1)%pn) / (old(k2)%pn - old(k1)%pn)
             ww    = ww    + w
             pn    = pn    + w * old(k)% pn
             tdbt  = tdbt  + w * log(old(k)% tdbt)
             tdbt2 = tdbt2 + w * log(old(k)% tdbt**2)
             mixr  = mixr  + w * log(old(k)% mixr)
             tdbt2 = mixr2 + w * log(old(k)% mixr**2)
          end do
          do k = k2 + 1, k3
             w    = (old(k)%pn - old(k3)%pn) / (old(k2)%pn - old(k3)%pn)
             ww   = ww   + w
             pn    = pn    + w * old(k)% pn
             tdbt  = tdbt  + w * log(old(k)% tdbt)
             tdbt2 = tdbt2 + w * log(old(k)% tdbt**2)
             mixr  = mixr  + w * log(old(k)% mixr)
             mixr2 = mixr2 + w * log(old(k)% mixr**2)
          end do
          satem_levels(i-1)% pn    =  pn    / ww
          satem_levels(i-1)% tdbt  =  exp(tdbt  / ww)
          satem_levels(i-1)% tdbt2 =                                       &
               max (0._wp, exp(tdbt2 / ww) - satem_levels(i-1)% tdbt **2)
          satem_levels(i-1)% mixr  =  exp(mixr  / ww)
          satem_levels(i-1)% mixr2 =                                       &
               max (0._wp, exp(mixr2 / ww) - satem_levels(i-1)% mixr **2)
       end do
    else
       satem_levels% tdbt2 = 0._wp
       satem_levels% mixr2 = 0._wp
    endif
    deallocate (old)

    !---------
    ! printout
    !---------
    if (verbose > 0) then
       call nl
       call nl; write (oline(iol),'(a,a)')     ' statid  : ',spot% statid
       call nl; write (oline(iol),'(a,a)')     ' date    : ',cyyyymmdd(satem% time)
       call nl; write (oline(iol),'(a,a)')     ' time    : ',chhmmss  (satem% time)
       call nl; write (oline(iol),'(a,2f10.3)')' lat/lon : ',satem% lah, satem% loh
       call nl; write (oline(iol),'(a)')
       call nl; write (oline(iol),'(a)') &
            ' selected SATEM levels and temperatures:'
       call nl; write (oline(iol),'(a)') &
            '   k         p       lon       lat       tdbt       dtdbt2       mixr       dmixr2'
       do k = 1, satem% nlevel
          call nl; write (oline(iol),'(i4,f12.3,2f10.3,4f12.4)') k,satem_levels(k)% pn,&
               satem% loh, satem% lah, &
               satem_levels(k)% tdbt, sqrt(satem_levels(k)% tdbt2), &
               satem_levels(k)% mixr, sqrt(satem_levels(k)% mixr2)
       end do
       call nl; write (oline(iol),'(a)')
       call nl; write (oline(iol),'(a)') repeat('-',79)
    endif

    if (layer_thickness) then
       !-------------------------
       ! for layer thickness:
       ! t --> tv
       !-------------------------
       do i = 1, satem% nlevel
          t  = satem_levels(i)% tdbt
          q  = satem_levels(i)% mixr
          satem_levels(i)% tdbt = tv_t_q(t, q) ! <-- store virt. temp. in tdbt !!!!!
       enddo

       !-------------------------------
       ! limits of t-profile used for
       ! spline-interpolation of rho^-1
       !-------------------------------
       l_top = 15                 ! highest t level
       l_bot = satem% nlevel      ! lowest  t level

       !--------------------------------------------
       ! count and indicate non-empty layers.
       ! satem% nlevel now is the number of
       ! layers (not no. of levels in t,q-retrieval)
       !--------------------------------------------
       nobs_layers = 0                         ! index to count non-empty layers
       do i = 1, satem% nlevel
          do j = 1, n_layer
             if ((satem_levels(i)%pn <= satem_levs(j  ))  .and.   &
                  (satem_levels(i)%pn > satem_levs(j+1))) then
                nobs_layers(j) = nobs_layers(j) + 1
             endif
          enddo
       enddo

       i = 0
       ind_layers  = 0
       do j = 1, n_layer
          if (nobs_layers(j) == 0) cycle
          i = i + 1
          ind_layers(i) = j
       enddo
       satem% nlevel = i                       ! number of non-empty layers

       !---------------------------------------------
       ! derive layer thickness [m] from t retrieval:
       !---------------------------------------------
       allocate(l_th(satem% nlevel))              ! tmp. vector to hold thickness
       select case (obs_layer_method)
       case(1)
          !---------------------------------
          ! spline interpolation of rho^-1,
          ! nmerical integration:
          !---------------------------------
          n_int = 5                               ! use 2**n_i_int intervals
                                                  ! for integration
          ! determine profile of rho^-1
          allocate(rho_inv(l_top:l_bot))
          allocate(rho_inv_2(l_top:l_bot))

          rho_inv = satem_levels(l_top:l_bot)% tdbt * Rgas / &
               satem_levels(l_top:l_bot)% pn

          ! initialize spline routine with 2nd derivatives
          call init_splinex(                     &
               satem_levels(l_top:l_bot)% pn,    &
               rho_inv,                          &
               rho_inv_2                         &
               )
          ! loop over SATEM layer thickness levels
          do i = 1, satem% nlevel
             call splint (                       &
                  satem_levels(l_top:l_bot)% pn, &
                  rho_inv,                       &
                  rho_inv_2,                     &
                  satem_levs(ind_layers(i)),     &
                  rho_inv_i)
             call splint (                       &
                  satem_levels(l_top:l_bot)% pn, &
                  rho_inv,                       &
                  rho_inv_2,                     &
                  satem_levs(ind_layers(i)+1),   &
                  rho_inv_i1)
             ! loop over refined integration intervals
             do i_int = 1, n_int
                if (i_int == 1) then
                   geop_i = 0.5_wp                                                   &
                        * (satem_levs(ind_layers(i)+1) - satem_levs(ind_layers(i)))  &
                        * (rho_inv_i + rho_inv_i1 )
                else
                   it = 2**(i_int-2)
                   tnm = it
                   del = (satem_levs(ind_layers(i)+1)-satem_levs(ind_layers(i)))/tnm
                   x = satem_levs(ind_layers(i)) + 0.5_wp * del

                   sum = 0._wp
                   ! midpoint rule in each integration interval
                   do j = 1, it
                      call splint (satem_levels(l_top:l_bot)% pn, &
                           rho_inv,                               &
                           rho_inv_2,                             &
                           x,                                     &
                           rho_inv_x)
                      sum = sum + rho_inv_x
                      x = x + del
                   enddo
                   geop_i = 0.5_wp                                                          &
                        * (geop_i + (satem_levs(ind_layers(i)+1)-satem_levs(ind_layers(i))) &
                        * sum/tnm)
                endif
             enddo
             ! tmp. layer thickness result [m]
             l_th(i) = - geop_i / gacc
          enddo
          deallocate(rho_inv,rho_inv_2)
       case (2)
          !----------------------------------
          ! linear interpolation of rho^-1,
          ! analytic integration:
          !----------------------------------
          allocate(int_geop(l_top:l_bot))
          int_geop = 0._wp
          ! geop. [m] between two neighbouring t values in retrieved profile
          do i = l_top, l_bot-1
             int_geop(i) = 0.5_wp *                                                &
                  (  satem_levels(i  )% tdbt/satem_levels(i  )% pn                 &
                  +  satem_levels(i+1)% tdbt/satem_levels(i+1)% pn)                &
                  * (satem_levels(i+1)%pn-satem_levels(i)%pn)*(Rgas/gacc)
          enddo
          ! loop over SATEM levels
          do j = 1, satem% nlevel
             ! find indices of upper, lower t value in each layer
             do i = l_top, l_bot-1
                if ((satem_levels(i)%  pn <= satem_levs(ind_layers(j)+1))    .and. &
                    (satem_levels(i+1)%pn >  satem_levs(ind_layers(j)+1)))         &
                     i_top = i+1
                if ((satem_levels(l_top)% pn >= satem_levs(ind_layers(j)  )) .and. &
                    (satem_levels(l_top)% pn <  satem_levs(ind_layers(j)+1)))      &
                     i_top = l_top
                if ((satem_levels(i  )%pn <= satem_levs(ind_layers(j)))      .and. &
                    (satem_levels(i+1)%pn >  satem_levs(ind_layers(j))))           &
                     i_bot = i
                if ((satem_levels(l_bot)% pn <  satem_levs(ind_layers(j)  )) .and. &
                    (satem_levels(l_bot)% pn >= satem_levs(ind_layers(j)+1)))      &
                     i_bot = l_bot
             enddo
             ! sum up geopoential in core section of each layer
             l_th_j = 0._wp
             do k = i_top, i_bot-1
                l_th_j = l_th_j + int_geop(k)
             enddo
             ! do some corrections for upper and lower edge of SATEM layer
             if (i_top > l_top) then
                w_top = (satem_levs(ind_layers(j)+1)-satem_levels(i_top)%pn) /           &
                     (satem_levels(i_top-1)%pn-satem_levels(i_top)%pn)
                l_th_j = l_th_j + w_top * int_geop(i_top-1)
             endif
             if (i_top == l_top)                                                         &
                  l_th_j = l_th_j + (satem_levels(l_top )% tdbt/satem_levels(l_top)% pn) &
                  * (satem_levels(l_top)%pn-satem_levs(ind_layers(j)+1))*(Rgas/gacc)
             if (i_bot < l_bot) then
                w_bot = (satem_levs(ind_layers(j))-satem_levels(i_bot)%pn) /             &
                     (satem_levels(i_bot+1)%pn-satem_levels(i_bot)%pn)
                l_th_j = l_th_j + w_bot * int_geop(i_bot)
             endif
             if (i_bot == l_bot) then
                l_th_j = l_th_j + satem_levels(l_bot )% tdbt/satem_levels(l_bot)% pn     &
                     * (-satem_levels(l_bot)%pn+satem_levs(ind_layers(j)))*(Rgas/gacc)
             endif
             ! tmp. layer thickness result [m]
             l_th(j) = l_th_j
          enddo
          deallocate(int_geop)
       case(3)
          !-------------------
          ! staircase rho^-1:
          ! average t in layer
          !-------------------
          l_th=0._wp
          ! determine t average in layer:
          do i = l_top, l_bot
             do j = 1, satem% nlevel
                if ((satem_levels(i)%pn <= satem_levs(ind(j)  ))  .and.   &
                    (satem_levels(i)%pn >  satem_levs(ind(j)+1)))  then
                   l_th(j) = l_th(j) + satem_levels(i)% tdbt
                endif
             enddo
          enddo
          ! loop over SATEM levels
          do i = 1, satem% nlevel
             ! tmp. layer thickness result [m]
             l_th(i) =                                                       &
                  l_th(i) /  nobs_layers(ind_layers(i)) * Rgas / gacc *      &
                  log(satem_levs(ind_layers(i))/satem_levs(ind_layers(i)+1))
          enddo
       case default
          call finish('check_store_satem','unknown integration method')
       end select

    endif

    !====================================
    ! store data in observation data type
    !====================================
    if (satem_int_size == 0) call set_size

    !----------
    ! meta data
    !----------
    spot% col% c% dlat        = satem% lah          ! LATITUDE (HIGH ACCURACY)
    spot% col% c% dlon        = satem% loh          ! LONGITUDE (HIGH ACCURACY)
    spot% col%    nlev        = satem% nlevel       ! NO. OF LEVELS
    spot%         int_type    = ITY_ICOL            ! INTERPOLATION TYPE: COLUMN
    spot%         ident       = satem% sat_id       ! SATELLITE IDENTIFIER
    spot%         actual_time = satem% time         ! DATE/TIME STAMP
    spot%         ps          = satem% ppp          ! PRESSURE
    spot%         z           = satem% hp           ! HEIGHT OF STATION
    spot%         stzen       = satem% szan         ! SATELLITE ZENITH ANGLE
    spot%         stclf       = satem% FLG_CLDFRM   ! CLOUD FORMATION AND HEIGHT ASSIGNMENT
    spot%         stlsf       = satem% lase         ! LAND/SEA QUALIFIER
    spot%         stret       = satem% FLG_040196   ! RETRIEVAL CHOICE INDICATOR
    spot%         sttyp       = satem% sain         ! SATELLITE INSTRUMENTS
    spot%         statid      = mnem(satem% sat_id) ! SATELLITE IDENTIFIER MNEMONIC

    call set_xuv  (spot)
    !--------------------------
    ! report selection (filter)
    !--------------------------
    call check_report_1 (spot)

    if (spot% use% state > STAT_DISMISS) then
       call new_spot (obs,1, set_id=.true.)
       spt => obs% spot (obs% n_spot)
       id      = spt% id
       spt     = spot
       spt% id = id
       !-------------
       ! observations
       !-------------
       bod % mn          = 'IASI SATEM'
       bod % use % state = STAT_ACTIVE
       bod % use % check = CHK_NONE

       if (layer_thickness) then
          !-------------------------------
          ! observations = layer thickness
          !-------------------------------
          call new_obs (obs, satem% nlevel , spot=spt)
          spt% nr             =  21                          ! no. of non-zero entries in R
                                                             ! according to ECMWF research manual 1
                                                             ! "ECMWF data assimilation -
                                                             !  Scientific documentation"
          do i=1,spt% o% n
             obs % varno (spt% o% i+i ) = VN_Z              ! varno h
             ! layer thickness from average t:
             obs %  body (spt% o% i+i)    = bod
             obs %  body (spt% o% i+i)% o = l_th(i)
             ! store midpoint pressure as obs pressure level:
             obs %  olev (spt% o% i+i) = &
                  (satem_levs(ind_layers(i))+satem_levs(ind_layers(i)+1))/2._wp
          enddo
          call new_int (obs, spt, n_layer+1)
       else
          !--------------------
          ! observations = t, q
          !--------------------
          call new_obs (obs, 2 * satem% nlevel, spot=spt)
          spt% nr = 2 * spt% col% nlev                                       ! t,q: R diagonal
          obs    % varno (spt% o% i+1 : spt% o% i + spt% o% n: 2)   = VN_T  ! varno t
          obs    % varno (spt% o% i+2 : spt% o% i + spt% o% n: 2)   = VN_Q  ! varno q
          do i = 1, 2
             obs % body  (spt% o% i+i : spt% o% i + spt% o% n: 2)   = bod    ! datum
             obs % olev  (spt% o% i+i : spt% o% i + spt% o% n: 2)   = &      ! obs pressure levels
                  satem_levels(1:spt% col% nlev)% pn
          enddo
          obs    % body  (spt% o% i+1 : spt% o% i + spt% o% n: 2)%o = &      ! obs t
               satem_levels(1:spt% col% nlev)% tdbt                          !
          obs    % body  (spt% o% i+2 : spt% o% i + spt% o% n: 2)%o = &      ! obs q
               satem_levels(1:spt% col% nlev)% mixr                          !
          call set_int_insitu (spt, obs)
          call store_satem (obs, spt, satem, satem_levels)
       endif

       !------------------------------------------------
       ! store information read in observation data type
       !------------------------------------------------
    endif

  end subroutine check_store_satem
!==============================================================================
end module mo_satem
