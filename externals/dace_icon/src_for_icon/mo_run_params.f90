!
!+ read 3D-Var runtime parameters from namelist /RUN/
!
MODULE mo_run_params
!
! Description:
!   General runtime parameters read from namelist /RUN/:
!     - Method to use
!     - MPI decomposition
!     - analysis time, reference time of forecast
!     - input/output paths
!     - input/output file names
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
! V1_2         2008/12/04 Andreas Rhodin
!  new namelist parameters: opt_fields, urun_fields
! V1_4         2009/03/26 Andreas Rhodin
!  changed namelist parameters (use sea ice mask: fr_ice)
! V1_5         2009/05/25 Andreas Rhodin
!  new namelist parameters: fdbk_basename, model
! V1_7         2009/08/24 Andreas Rhodin
!  new namelist parameter : range_ana (set time range value for COSMO/LETKF)
! V1_8         2009/12/09 Andreas Rhodin
!  define GRIB fields for COSMO
!  new namelist parameters: interp_strato,
!                           file_strato (interpolate IFS stratosphere)
!  function path_file: replace string "_EXP_" by experiment number
!                      replace  _ANA_TIMEMM_ or _ANA_TIMEMMSS_ by analysis time
!                      same for _FCR_TIMEMM_,   _FCR_TIMEMMSS_
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Andreas Rhodin
!  adjust run_type,nex if nex/2**14/=0 (misspecified in namelist)
! V1_13        2011/11/01 Andreas Rhodin
!  changes for GRIB2 API, LETKF, printout of host id
! V1_15        2011/12/06 Andreas Rhodin
!  change logical namelist parameter interp_strato to integer value
! V1_19        2012-04-16 Andreas Rhodin
!  changes for GME LETKF
! V1_20        2012-06-18 Andreas Rhodin
!  subroutine set_nprocs: make nprocs optional
! V1_22        2013-02-13 Harald Anlauf
!  changes for ICON / GRIB-API / CDI
! V1_23        2013-03-26 Andreas Rhodin
!  replace forecast lead time string in path_file
! V1_26        2013/06/27 Harald Anlauf
!  Changes for GRIB2/GRIB_API/ICON
! V1_27        2013-11-08 Harald Anlauf
!  for ICON: always require qr qs fr_ice, remove t2m,td2m from default settings
! V1_28        2014/02/26 Harald Anlauf
!  Fixes for TEMP verification of IFS fields;
!  Namelist /run/: default ga3_biasc_airep to 2099010100 (off)
! V1_29        2014/04/02 Harald Anlauf
!  Revise printing of host to MPI rank assignment
! V1_31        2014-08-21 Andreas Rhodin
!  preparations for MEC:
!    change default in namelist /run/ : aux = output ,
!    new optional parameters anatime,fcreftime,leadtime to function path_file
! V1_35        2014-11-07 Harald Anlauf
!  implement blocking of gather operation to limit memory usage;
!  to use, see namelist run::io_max_gather
! V1_40        2015-02-27 Harald Anlauf
!  Diagnose compiler vendor and version
! V1_43        2015-08-19 Harald Anlauf
!  set defaults for MEC.
!  path_file: _ANDATE_ (yyyymmdd).
!  Namelist /RUN/, dealloc_fields: remove fields from analysis file.
! V1_44        2015-09-30 Hendrik Reich
!  add 'cosmo_refatm' (reference atmosphere to be used for COSMO) to namelist /RUN/
! V1_45        2015-12-15 Harald Anlauf
!  Namelist /RUN/, ready_det: ready file for deterministic analysis
! V1_46        2016-02-05 Andreas Rhodin
!  fallback to default COSMO ivctype and irefatm if not specified in GRIB1
! V1_47        2016-06-06 Harald Anlauf
!  Enable GRIB encoding of selected analysis fields with 24 bits
!  write COSMO vcp to GRIB, changes for refatm=2 (COMET)
! V1_48        2016-10-06 Andreas Rhodin
!  support COSMO ivctype=3,4 (sleeve coordinates)
!  namelist /run/: new parameter 'otime_shift'
!                  no default initialisation for 'adjust_sst_snow'
! V1_49        2016-10-25 Alexander Cress
!  Online bias correction for SCATT
! V1_50        2017-01-09 Andreas Rhodin
!  new namelist parameter otime_include, new defaults grib_edition,grib_library
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002       original code
!                      2006-2008  changes
!---------------------------------------------------------------

  !=============
  ! Modules used
  !=============
  use mo_kind,          only: i8, wp            ! integer*8 kind parameter
  use mo_namelist,      only: nnml,            &! file unit of namelist file
                              position_nml,    &! position namelist group
                              POSITIONED        ! return flag
  use mo_mpi_dace,      only: dace,            &! DACE parallelisation info
                              p_parallel,      &! flag for parallel mode
                              p_bcast           ! overloaded broadcast routine
  use mo_time,          only: t_time,          &! time data type
                              init_time,       &! set time data type
                              now,             &! returns actual time
                              zero_time,       &! zero time interval
                              operator(+),     &! add      times
                              operator(-),     &! subtract times
                              operator(/),     &! divide   times
                              operator(/=),    &! compare  times
                              chhhmm,          &! return string
                              chhhmmss,        &! return string
                              cddhhmmss,       &! return string
                              cdddhhmmss,      &! return string
                              cyyyymmdd,       &! return string
                              cyyyymmddhh,     &! return string
                              cyyyymmddhhmm,   &! return string
                              cyyyymmddhhmmss, &! return string
                              time_yyyymmddhh, &! derive time from integer
                              sign,            &! compare times
!                             ihh,             &! derive hours
                              hours,           &! time to hours (real)
                              p_bcast           ! broadcast derived type
  use mo_dace_string,   only: char5,           &! integer -> char(len=5)
                              char3,           &! integer -> char(len=3)
                              eval_string       ! evaluate strings
  use mo_gribtables,    only: dwd_to_v3d        ! convert to 3dvar convention
  use mo_grib_handling, only: set_grib,        &! set GRIB edition and library
                              grib_api_version,&! get grib_api_version
                              gribtabversion,  &! GRIB2 minimum tablesVersion
                              grib_packing      ! GRIB packing or compression
  use mo_exception,     only: finish            ! return on error condition
  use mo_fortran_units, only: get_unit_number, &! request an unused unit number
                              return_unit_number! return the unit number
  use mo_system,        only: getlogin,        &! get user name
                              gethostname       ! get host name
  use mo_version,       only: var3d_version,   &! return version number
                              var3d_date,      &!        version date
                              var3d_time,      &!        version time
                              compiler_info     ! compiler vendor and version
  use mo_wmo_tables,    only: DWD6_ICOSAHEDRON,&! WMO table 6: grid-type
                              DWD6_ICON,       &!
                              WMO6_GAUSSIAN,   &!
                              WMO6_ROTLL
  implicit none

  !================
  ! Public entities
  !================
  private
  !----------------------------------
  ! subroutine to read namelist /RUN/
  !----------------------------------
  public :: nml_run_flags
  !---------------------------------
  ! subroutine to set nproc1, nproc2
  !---------------------------------
  public :: set_nprocs
  !------------------------------------------
  ! namelist variables and derived parameters
  !------------------------------------------
  public :: method, nproc1, nproc2, nproc_repro, barrier, check_sync
  public :: ana_time, fc_ref_time, fc_time, fc_hours, obs_timeshift
  public :: otime_exclude
  public :: data, iopath, input, output, obsinput, aux
  public :: fg_file, blacklists, oldanerr_file, invar_det
  public :: ana_file, ana_err_file, grid_file !, out_file
  public :: p_readgrib, p_readbufr, npe_read_obs
  public :: nex, run_type, runtype, pass_fields
  public :: read_fields, opt_fields, fdbk_addvar, urun_fields, dealloc_fields, gribout_24bit
  public :: expseed, sst_fields, ana_date, model
  public :: fdbk_basename
  public :: interp_strato  ! flag to interpolate IFS stratosphere
  public :: file_strato    ! stratosphere file name
  public :: grid_file_strato! stratosphere grid file
  public :: abort_strato   ! flag to abort if file is not present
! public :: grib_edition   ! GRIB edition to write (-1=input,1,2)
! public :: grib_library   ! GRIB-API used: 1=GRIBEX  2=GRIB2-API
  public :: read_obs_first ! flag to read obs. before background
  public :: lcost_obs      ! estimate cost (for reading) from number of obs
  public :: io_max_gather  ! gather buffer size for I/O
  public :: io_max_scatter ! scatter buffer size for I/O
  public :: cosmo_refatm   ! reference atmosphere     to be used for COSMO:1or2
  public :: cosmo_ivctype  ! vertical coordinate type to be used for COSMO:1to4
  public :: cosmo_svc      ! large-scale, small-scale decay rate for ivctype 3,4
  public :: cosmo_nfltvc   ! number of filter applications  used for ivctype 3,4
  public :: l_ke_in_gds    ! explicit GDS entry for number of model levels
  public :: ready_det      ! ready file, determ. analysis
  public :: nature_run     ! run to generate 'truth', no: bc, fg-check, ana.incr.
  public :: chk_icon_grid  ! check for correct order of neighbour references
  public :: debug_rss      ! Memory usage debugging level
  public :: node_id        ! Node index for rank
  public :: geoid_format   ! file format of geoid (ASCII or netCDF)
  public :: geoid_file     ! geoid file name
  !----------------------------------------------------
  ! switch behaviour of analysis system at a given time
  !----------------------------------------------------
  public :: flag_biasc_rad     ! radiance online bias correction
  public :: flag_biasc_airep   ! aircraft bias correction
  public :: flag_biasc_synop   ! SYNOP    bias correction
  public :: flag_biasc_gpsgb   ! GPSGB    bias correction
  public :: flag_biasc_scatt   ! SCATT    bias correction
  public :: flag_biasc_wlidar  ! WLIDAR  bias correction
  !--------------------
  ! general information
  !--------------------
  public :: host               ! host name
  public :: user               ! user name
  public :: run_time           ! time of run of this program
  public :: center             ! generating center id (default DWD)
  public :: subcenter          ! subcenter id         (default none)
  public :: proc_ana_err       ! generating process for analysis error
  public :: proc_ana           ! generating process for analysis
  public :: proc_ana_ens       ! generating process for analysis ensemble
  public :: range_ana          ! time range parameter for analysis
  public :: ensemble_id        ! ensemble id for GRIB output
  !---------
  ! routines
  !---------
  public :: path_file          ! concatenate: path / file . suffix
  public :: basename           ! strip leading directories from path
  public :: process_grid       ! derive generating process from grid
  public :: cleanup_run_params ! clean up this module

  !------------------------
  ! steering of the program
  !------------------------
  character(len=12)  :: method       = 'PSAS' ! 'PSAS'  (default 3dvar mode)
                                              ! 'PSAS+LETKF'
                                              ! 'LETKF'
                                              ! 'ENVAR', 'ENVAR+LETKF'
                                              ! 'GMESTAT', 'MEC'
  integer            :: expseed       = 0     ! nondefault random number seed
  integer            :: interp_strato = 0     ! interpolate stratosphere
  logical            :: abort_strato  =.false.! abort if file is not present
  logical            :: read_obs_first=.false.! read obs. before background
  logical            :: lcost_obs     =.false.! get cost (for reading)from #obs
  logical            :: nature_run    =.false.! generate 'truth'
  integer            :: chk_icon_grid = -1    ! check neighbour references 0/1=off/on
  !-----------------------
  ! specification of model
  !-----------------------
  character(len=8)   :: model        = 'ICON' ! or 'COSMO', 'IFS','GME'
  integer            :: cosmo_refatm = -1     ! reference atmosphere for COSMO
  integer            :: cosmo_ivctype= -1     ! reference atmosphere for COSMO
  real(wp)           :: cosmo_svc(2) = 0._wp  ! large,small-scale decay rate
  integer            :: cosmo_nfltvc = 0      ! number of filter applications
  logical            :: l_ke_in_gds  = .false.! explicit GDS entry for number of model levels
  !-------------------------------------------
  ! fields to read and write from/to GRIB file
  !-------------------------------------------
  character(len=512) :: read_fields   = ' '   ! fields to read from forecast
  character(len=512) :: sst_fields    = ' '   ! read from SST analysis (00UTC)
  character(len=512) :: opt_fields    = ' '   ! read optionally
  character(len=512) :: fdbk_addvar   = ' '   ! fdbk output optionally
  character(len=512) :: urun_fields   = ' '   ! unspecified runtype
  character(len=512) :: pass_fields   = ' '   ! pass to analysis (relabel)
  character(len=128) :: dealloc_fields= ' '   ! do not write to analysis
  character(len=128) :: gribout_24bit = ' '   ! fields to write with 24 bits
  integer            :: run_type      =  3    ! haupt=0, (vor=1), ass=2, test=3
  character(len=8)   :: runtype       = 'forecast' ! set to '' to read analyses
  !--------------------
  ! MPI parallelisation
  !--------------------
  integer            :: nproc1        = -1     ! number of PEs in direction 1
  integer            :: nproc2        = -1     ! number of PEs in direction 2
  integer            :: nproc_repro   = -1     ! fict. PEs for repeatable runs
  logical            :: barrier       = .true. ! call mpibarrier for diagnostics
  logical            :: check_sync    = .false.! check syncronisation(stop_time)
  integer            :: p_readgrib    = -1     ! PE used to read GRIB
  integer            :: p_readbufr    = -1     ! PE used to read BUFR
  integer            :: npe_read_obs  = 999999 ! number of PEs to read obsv
  integer(i8)        :: io_max_gather = -1     ! I/O gather buffer size [bytes]
  integer(i8)        :: io_max_scatter= -1     ! I/O scatter buffer size [bytes]
  integer            :: debug_rss     =  0     ! Memory usage debugging level
  !------
  ! dates
  !------
  type(t_time) ,save :: ana_time              ! analysis time
  type(t_time) ,save :: ana_date              ! analysis date (0 UTC)
  type(t_time) ,save :: fc_ref_time           ! reference time (forecast start)
  type(t_time) ,save :: fc_time               ! forecast interval
  type(t_time) ,save :: obs_timeshift         ! shift of observation time window
  real(wp)           :: fc_hours       = -1   ! forecast interval (hours)
  integer            :: yyyymmddhh_ana =  0   ! analysis time
  integer            :: yyyymmddhh_ref =  0   ! time of forecast start
  integer(i8)        :: time_ana       =  0   ! analysis time  (ccccmmddhhmmss)
  integer(i8)        :: time_ref       =  0   ! forecast start (ccccmmddhhmmss)
  integer(i8)        :: runtime        =  0   ! program run time YYYYMMDDHHMM[SS]
  integer            :: otime_shift    = -9999! observation window shift (hhmm)
  integer            :: otime_include  =  0   ! -1:start,  +1: end of interval
  logical            :: otime_exclude  =.true.! T for ICON, F for COSMO
  !----------------
  ! directory paths
  !----------------
  character(len=256) :: data           = ''   ! constant input data
  character(len=256) :: iopath         = ''   ! path for input/output dirs
  character(len=256) :: input          = ''   ! input             directory
  character(len=256) :: output         = ''   ! output            directory
  character(len=256) :: obsinput       = ''   ! observation input directory
  character(len=256) :: aux            = ''   ! additional output directory
  !-----------------
  ! input file names
  !-----------------
  character(len=128) :: fg_file        = ''   ! forecast file
  character(len=128) :: blacklists(5)  = ''   ! blacklists
  character(len=128) :: oldanerr_file  = ''   ! analysis error
  character(len=128) :: invar_det      = ''   ! invariant fields
  character(len=128) :: grid_file      = ''   ! grid metadata (ICON)
  character(len=128) :: file_strato    = ''   ! stratosphere file name
  character(len=128) :: grid_file_strato = '' ! grid metadata for file_strato
  character(len=128) :: geoid_format   = 'ascii'       ! type/version of geoid
  character(len=128) :: geoid_file     = 'ww15mgh.grd' ! geoid file name
  !----------------------------
  ! output file characteristics
  !----------------------------
  character(len=128) :: ana_file       = ''   ! analysis
  character(len=128) :: ana_err_file   = ''   ! analysis error
! character(len=128) :: out_file       = ''   ! protocol
  character(len=128) :: fdbk_basename  = ''   ! feedback file
  character(len=128) :: ready_det      = ''   ! ready file, determ. analysis
  integer            :: grib_edition   =  2   ! GRIB edition to write (-1,1,2)
  integer            :: grib_library   =  2   ! API: 1=(C)GRIBEX  2=GRIB2-API
  !--------------------------------------------------------
  ! identification of center, process, experiment, ensemble
  !--------------------------------------------------------
  integer            :: center       =  78 ! DWD
  integer            :: subcenter    = 255 ! none
  integer            :: proc_ana_err =   0 ! process: analysis error
  integer            :: proc_ana     =   0 ! process: analysis
  integer            :: proc_ana_ens =   0 ! process: analysis ensemble
  integer            :: range_ana    =   0 ! time range parameter for analysis
  integer            :: nex          =  0     ! experiment number
  integer            :: ensemble_id  =   1 ! ensemble id for GRIB output
  !----------------------------------------------------
  ! switch behaviour of analysis system at a given time
  !----------------------------------------------------
  integer            :: ga3_biasc_rad   = 2099010100 ! radiance online biascor.
  integer            :: ga3_biasc_airep = 2099010100 ! aircraft bias correction
  integer            :: ga3_biasc_synop = 2099010100 ! SYNOP    bias correction
  integer            :: ga3_biasc_gpsgb = 2099010100 ! GPSGB    bias correction
  integer            :: ga3_biasc_scatt = 2099010100 ! SCATT    bias correction
  integer            :: ga3_biasc_wlidar= 2099010100 ! WLIDAR   bias correction
  type(t_time) ,save :: date_biasc_rad
  type(t_time) ,save :: date_biasc_airep
  type(t_time) ,save :: date_biasc_synop
  type(t_time) ,save :: date_biasc_gpsgb
  type(t_time) ,save :: date_biasc_scatt
  type(t_time) ,save :: date_biasc_wlidar
  integer            :: flag_biasc_rad   = -1        ! radiance online biascor.
  integer            :: flag_biasc_airep = -1        ! aircraft bias correction
  integer            :: flag_biasc_synop = -1        ! SYNOP    bias correction
  integer            :: flag_biasc_gpsgb = -1        ! GPSGB    bias correction
  integer            :: flag_biasc_scatt = -1        ! SCATT    bias correction
  integer            :: flag_biasc_wlidar= -1        ! WLIDAR   bias correction

  !===============
  ! namelist /RUN/
  !===============
  namelist /RUN/ method, nproc1, nproc2, nproc_repro, nature_run,            &
                 model, barrier, runtime, check_sync, chk_icon_grid,         &
                 yyyymmddhh_ana, yyyymmddhh_ref, fc_hours,                   &
                 iopath, data, input, output, aux, obsinput,                 &
                 fg_file, blacklists, oldanerr_file, grid_file,              &
                 ana_file, ana_err_file, invar_det, ready_det,               &
                 p_readgrib, nex, run_type, npe_read_obs, lcost_obs,         &
                 proc_ana, proc_ana_err, proc_ana_ens, center, subcenter,    &
                 runtype, pass_fields, read_fields,                          &
                 opt_fields, fdbk_addvar, urun_fields, sst_fields,           &
                 expseed, fdbk_basename, range_ana,                          &
                 ensemble_id, time_ana, time_ref, otime_shift,               &
                 interp_strato, file_strato, grid_file_strato, abort_strato, &
                 grib_edition, grib_library, read_obs_first, otime_include,  &
                 ga3_biasc_rad, ga3_biasc_airep, ga3_biasc_synop,            &
                 ga3_biasc_gpsgb, ga3_biasc_scatt, ga3_biasc_wlidar,         &
                 cosmo_refatm, cosmo_ivctype, cosmo_svc, cosmo_nfltvc,       &
                 l_ke_in_gds, dealloc_fields, gribout_24bit, gribtabversion, &
                 grib_packing, io_max_gather, io_max_scatter, debug_rss,     &
                 geoid_format, geoid_file

  ! no longer used: out_file
  !---------------
  ! run parameters
  !---------------
  type(t_time) ,save              :: run_time           ! time of this model run
  character (len=16)              :: user = ""
  character (len=16)              :: host = ""
  character (len=16) ,allocatable :: hosts  (:)
  integer ,protected ,allocatable :: node_id(:)         ! Node index for rank

contains
!------------------------------------------------------------------------------
  subroutine nml_run_flags (check)
  logical ,optional, intent(in) :: check ! flag: check for presence of files
  !--------------------
  ! read namelist /RUN/
  !--------------------
!   integer             :: pe
    integer             :: ierr
    integer             :: yyyy,mm,dd,hh,mi !,ss
    integer             :: iunit, ios
    integer             :: major, minor, rev
    character(len=10)   :: c_time
    logical             :: chk
    character(len=96)   :: vendor, ftn_version  ! Fortran compiler information

    if (dace% lpio) then
      !--------------
      ! read namelist
      !--------------
      call position_nml ('RUN', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=RUN, iostat=ierr)
        if (ierr/=0) call finish ('nml_run_flags','ERROR in namelist /RUN/')
#else
        read (nnml ,nml=RUN)
#endif
      case default
        call finish('nml_run_flags','namelist /RUN/ must be present.')
      end select
      !-------------------
      ! consistency checks
      !-------------------
      if (nproc1<=0 .and. nproc2<=0) call set_nprocs (nproc1, nproc2, dace% npe)
      if (nproc1<=0 .and. nproc2> 0) nproc1 = dace% npe / nproc2
      if (nproc1> 0 .and. nproc2<=0) nproc2 = dace% npe / nproc1
      nproc1 = max (nproc1,1)
      nproc2 = max (nproc2,1)
      if (nproc1 * nproc2 /= dace% npe) then
        write (0,*)
        write (0,*) 'nproc1    =',nproc1
        write (0,*) 'nproc2    =',nproc2
        write (0,*) 'dace% npe =',dace% npe
        write (0,*)
        call finish ('nml_run_flags','nproc1 * nproc2 /= dace% npe')
      endif
      if (p_readgrib >= dace% npe) p_readgrib = dace% pio
      if (p_readgrib <  0        ) p_readgrib = dace% pio
      if (p_readbufr >= dace% npe) p_readbufr = dace% pio
      if (p_readbufr <  0        ) p_readbufr = dace% pio
      if (nex / 2**14 /= 0) then
        write (6,'()')
        write (6,'(a,i6,i12)') '  namelist /run/: run_type, nex      read :'&
                                                 ,run_type, nex
        if (run_type < 0 .or. run_type > 2) run_type = nex / 2**14
        nex = mod (nex, 2**14)
        write (6,'(a,i6,i12)') '                  run_type, nex forced to :'&
                                                 ,run_type, nex
        write (6,'()')
      endif
      if (method == 'VARENKF') method = 'ENVAR+LETKF' ! VARENKF is alias

      if (nature_run) then
        !--------------------------------------------------
        ! nature run: degrade to 'PSAS', no bias correction
        !--------------------------------------------------
        method = 'PSAS'
        ga3_biasc_rad   = 2099010100
        ga3_biasc_airep = 2099010100
        ga3_biasc_synop = 2099010100
        ga3_biasc_gpsgb = 2099010100
        ga3_biasc_scatt = 2099010100
        ga3_biasc_wlidar= 2099010100
      endif

      !-------------------------------
      ! set defaults for ICON or COSMO
      !-------------------------------
      call dwd_to_v3d (pass_fields)   ! convert from DWD to 3dvar convention

      select case (model)
      case ('GME','IFS')

        select case (method)
        case ('PSAS')
          call eval_string (pass_fields,'qci t2m td2m qr qs '//pass_fields)
        case default
          call eval_string (pass_fields,'qci qr qs '//pass_fields)
        end select
        call eval_string   (opt_fields, 'qr qs fr_ice '//opt_fields)
        call eval_string   (read_fields,'ps psr tsurf t u v q qcl fr_ice '//&
                            read_fields                                     )
        call eval_string   (urun_fields, urun_fields)
        call eval_string   (sst_fields,  sst_fields)
        call eval_string   (read_fields, pass_fields, mode=1)

      case ('ICON')                                          ! ++++ ADJUST ++++

        select case (method)
!       case ('PSAS')
!         call eval_string (pass_fields,'qci qr qs t2m td2m '//pass_fields)
        case ('MEC')
          pass_fields = ''
        case default
          call eval_string (pass_fields,'qci qr qs '//pass_fields)
          call eval_string   (opt_fields, 'ps '       //opt_fields)
          call eval_string   (read_fields,'pf tsurf t u v q qcl fr_ice '//&
                              read_fields                                 )
          call eval_string   (urun_fields, urun_fields)
          call eval_string   (sst_fields,  sst_fields)
          call eval_string   (read_fields, pass_fields, mode=1)
        end select

      case ('COSMO')                                         ! ++++ ADJUST ++++

        select case (method)
        case ('ADJUST')
        case ('GMESTAT')
        case ('MEC')
          pass_fields = ''
        case default
          call eval_string (pass_fields,'qci tsurf w qr qs '//pass_fields)
          call eval_string (opt_fields, 'qg fr_ice hsurf &
               &lsm z0 soiltyp plcov lai rootdp vio3 hmo3 t_snow w_i &
               &qv_s w_snow t_s t_so w_so freshsnw for_e for_d rho_snow '//opt_fields)
          call eval_string (urun_fields, 'pp tsurf t u v w q qcl qci qr qs qg fr_ice hsurf &
               &lsm z0 soiltyp plcov lai rootdp vio3 hmo3 t_snow w_i &
               &qv_s w_snow t_s t_so w_so freshsnw for_e for_d rho_snow '//urun_fields)
          call eval_string (sst_fields,  sst_fields)
          call eval_string (read_fields, 'pp t u v w q qcl qg fr_ice hsurf &
              &lsm z0 soiltyp plcov lai rootdp vio3 hmo3 t_snow w_i &
              &qv_s w_snow t_s t_so w_so freshsnw for_e for_d rho_snow '//read_fields)
          call eval_string (read_fields, pass_fields, mode=1)
        end select

      end select
      call eval_string   (fdbk_addvar, fdbk_addvar)

    endif

    !-------------------
    ! broadcast namelist
    !-------------------
    if (p_parallel) then
      call p_bcast (method          ,dace% pio)
      call p_bcast (model           ,dace% pio)
      call p_bcast (cosmo_refatm    ,dace% pio)
      call p_bcast (cosmo_ivctype   ,dace% pio)
      call p_bcast (cosmo_svc       ,dace% pio)
      call p_bcast (cosmo_nfltvc    ,dace% pio)
      call p_bcast (l_ke_in_gds     ,dace% pio)
      call p_bcast (nproc1          ,dace% pio)
      call p_bcast (nproc2          ,dace% pio)
      call p_bcast (nproc_repro     ,dace% pio)
      call p_bcast (barrier         ,dace% pio)
      call p_bcast (check_sync      ,dace% pio)
      call p_bcast (yyyymmddhh_ana  ,dace% pio)
      call p_bcast (yyyymmddhh_ref  ,dace% pio)
      call p_bcast (time_ana        ,dace% pio)
      call p_bcast (time_ref        ,dace% pio)
      call p_bcast (runtime         ,dace% pio)
      call p_bcast (fc_hours        ,dace% pio)
      call p_bcast (otime_shift     ,dace% pio)
      call p_bcast (otime_include   ,dace% pio)
      call p_bcast (otime_exclude   ,dace% pio)
      call p_bcast (data            ,dace% pio)
      call p_bcast (iopath          ,dace% pio)
      call p_bcast (input           ,dace% pio)
      call p_bcast (obsinput        ,dace% pio)
      call p_bcast (output          ,dace% pio)
      call p_bcast (aux             ,dace% pio)
      call p_bcast (fg_file         ,dace% pio)
      call p_bcast (ana_file        ,dace% pio)
      call p_bcast (ana_err_file    ,dace% pio)
!     call p_bcast (out_file        ,dace% pio)
      call p_bcast (invar_det       ,dace% pio)
      call p_bcast (ready_det       ,dace% pio)
      call p_bcast (grid_file       ,dace% pio)
      call p_bcast (fdbk_basename   ,dace% pio)
      call p_bcast (blacklists      ,dace% pio)
      call p_bcast (oldanerr_file   ,dace% pio)
      call p_bcast (p_readgrib      ,dace% pio)
      call p_bcast (p_readbufr      ,dace% pio)
      call p_bcast (npe_read_obs    ,dace% pio)
      call p_bcast (lcost_obs       ,dace% pio)
      call p_bcast (io_max_gather   ,dace% pio)
      call p_bcast (io_max_scatter  ,dace% pio)
      call p_bcast (nex             ,dace% pio)
      call p_bcast (run_type        ,dace% pio)
      call p_bcast (runtype         ,dace% pio)
      call p_bcast (pass_fields     ,dace% pio)
      call p_bcast (read_fields     ,dace% pio)
      call p_bcast (opt_fields      ,dace% pio)
      call p_bcast (fdbk_addvar     ,dace% pio)
      call p_bcast (urun_fields     ,dace% pio)
      call p_bcast (sst_fields      ,dace% pio)
      call p_bcast (dealloc_fields  ,dace% pio)
      call p_bcast (gribout_24bit   ,dace% pio)
      call p_bcast (proc_ana        ,dace% pio)
      call p_bcast (proc_ana_err    ,dace% pio)
      call p_bcast (range_ana       ,dace% pio)
      call p_bcast (proc_ana_ens    ,dace% pio)
      call p_bcast (center          ,dace% pio)
      call p_bcast (subcenter       ,dace% pio)
      call p_bcast (ensemble_id     ,dace% pio)
      call p_bcast (expseed         ,dace% pio)
      call p_bcast (interp_strato   ,dace% pio)
      call p_bcast (file_strato     ,dace% pio)
      call p_bcast (grid_file_strato,dace% pio)
      call p_bcast (abort_strato    ,dace% pio)
      call p_bcast (grib_edition    ,dace% pio)
      call p_bcast (grib_library    ,dace% pio)
      call p_bcast (gribtabversion  ,dace% pio)
      call p_bcast (grib_packing    ,dace% pio)
      call p_bcast (read_obs_first  ,dace% pio)
      call p_bcast (ga3_biasc_rad   ,dace% pio)
      call p_bcast (ga3_biasc_airep ,dace% pio)
      call p_bcast (ga3_biasc_synop ,dace% pio)
      call p_bcast (ga3_biasc_gpsgb ,dace% pio)
      call p_bcast (ga3_biasc_scatt ,dace% pio)
      call p_bcast (ga3_biasc_wlidar,dace% pio)
      call p_bcast (nature_run      ,dace% pio)
      call p_bcast (chk_icon_grid   ,dace% pio)
      call p_bcast (debug_rss       ,dace% pio)
      call p_bcast (geoid_format    ,dace% pio)
      call p_bcast (geoid_file      ,dace% pio)
    endif

    !--------------
    ! get hostnames
    !--------------
    call get_hostnames ()
    if (debug_rss > 0)  &
        call get_nodes ()

    !----------------------------------
    ! derive actual time, machine, etc.
    !----------------------------------
    call getlogin    (user)
    if (runtime == 0) then
      run_time = now()
    else          ! yyyymmddhhmm[ss]
      if (runtime > 999999999999_i8) then
!       ss   = mod (runtime , 100_i8)
        runtime =   runtime / 100_i8
      endif
      yyyy =      runtime / 100000000_i8
      mm   = mod (runtime / 1000000_i8   ,100_i8)
      dd   = mod (runtime / 10000_i8     ,100_i8)
      hh   = mod (runtime / 100_i8       ,100_i8)
      mi   = mod (runtime                ,100_i8)
      call init_time (run_time, yyyy=yyyy, mo=mm, dd=dd, hh=hh, mi=mi)
      if (yyyy > 9999 .or. yyyy < 1000) then
        if (dace% lpio) &
          write (0,'(a,i0)') '  runtime possibly out of range?  yyyy = ', yyyy
      end if
    endif
    call p_bcast (run_time, dace% pio)
    call p_bcast (user,     dace% pio)
    !---------------------
    ! derive analysis time
    !---------------------
    if (yyyymmddhh_ana /= 0 .and. time_ana == 0) &
      time_ana = yyyymmddhh_ana * 10000_i8
    if (time_ana /= 0) then
      call init_time (ana_time, date=time_ana)
      ana_date = ana_time
      ana_date% secs = 0
    else
      call finish('nml_run_flags','yyyymmddhh_ana or time_ana must be present in /RUN/')
    endif
    !-------------------------------------
    ! set analysis date dependent switches
    !-------------------------------------
    date_biasc_rad   = time_yyyymmddhh (ga3_biasc_rad)
    date_biasc_airep = time_yyyymmddhh (ga3_biasc_airep)
    date_biasc_synop = time_yyyymmddhh (ga3_biasc_synop)
    date_biasc_gpsgb = time_yyyymmddhh (ga3_biasc_gpsgb)
    date_biasc_scatt = time_yyyymmddhh (ga3_biasc_scatt)
    date_biasc_wlidar= time_yyyymmddhh (ga3_biasc_wlidar)
    flag_biasc_rad   = sign (1, ana_time - date_biasc_rad  )
    flag_biasc_airep = sign (1, ana_time - date_biasc_airep)
    flag_biasc_synop = sign (1, ana_time - date_biasc_synop)
    flag_biasc_gpsgb = sign (1, ana_time - date_biasc_gpsgb)
    flag_biasc_scatt = sign (1, ana_time - date_biasc_scatt)
    flag_biasc_wlidar= sign (1, ana_time - date_biasc_wlidar)
    !------------------------------
    ! derive time of forecast start
    !------------------------------
    if (yyyymmddhh_ref /= 0 .and. time_ref == 0) &
      time_ref = yyyymmddhh_ref * 10000_i8
    if (fc_hours <  0._wp .and. time_ref == 0) fc_hours = 3._wp   ! default
    if (fc_hours >= 0._wp) call init_time (fc_time, ss=nint (fc_hours*3600))
    if (time_ref /= 0) then
      call init_time (fc_ref_time, date=time_ref)
      if (fc_hours >= 0._wp) then
        if (fc_ref_time + fc_time /= ana_time)                &
          call finish('nml_run_flags',                        &
                      'time_ref and fc_hours are inconsistent')
      else
        fc_time  = ana_time - fc_ref_time                   ! derive interval
        fc_hours = hours (fc_time)
      endif
    else if (fc_hours >= 0._wp) then                        ! derive start time
      fc_ref_time = ana_time - fc_time
    else
      call finish('nml_run_flags',&
                  'yyyymmddhh_ref or fc_hours must be present in /RUN/')
    endif

    !---------------------------------
    ! shift ob observation time window
    !---------------------------------
    if (otime_shift == -9999) then
      select case (method)
      case ('MEC')
        obs_timeshift = zero_time
        otime_exclude  = .false.
      case default
        select case (model)
        case ('GME','IFS','ICON')
          obs_timeshift = fc_time / 2._wp
          otime_exclude = .true.
        case default
          obs_timeshift = zero_time
          otime_exclude  = .false.
        end select
      end select
    else
      call init_time (obs_timeshift, hhmmss = 100 * otime_shift)
    endif
    select case (otime_include)
    case (1)
      otime_exclude  = .false.
    case (-1)
      otime_exclude = .true.
    end select
    otime_include = 1; if (otime_exclude) otime_include = -1

    !----------
    ! set paths
    !----------
    if (data     =='') data     = '../data' ;data   = path_file ('./'   ,data)
    if (iopath   =='') iopath   = '../io'   ;iopath = path_file ('./'   ,iopath)
    if (input    =='') input    = 'input'   ;input  = path_file (iopath ,input)
    if (output   =='') output   = 'output'  ;output = path_file (iopath ,output)
    if (obsinput =='') obsinput =  input
    if (aux      =='') aux      =  output   ;aux    = path_file (output ,aux)
    !---------------------
    ! set input file names
    !---------------------
    c_time = cyyyymmddhh(fc_ref_time)
    if (oldanerr_file =='') oldanerr_file = 'i_ga_err_'//c_time
    if (method =='GMESTAT' .or. &
        method =='MEC'    ) oldanerr_file = '/none/'
    !----------------------
    ! set output file names
    !----------------------
    c_time = cyyyymmddhh(ana_time)
    if (ana_err_file == '') ana_err_file = 'i_ga_err_'//c_time
    if (ana_file     == '') ana_file     = 'ixxxa_ga_'//c_time
!   if (out_file     == '') out_file     = 'out_'  //c_time
    !-----------------------
    ! add paths to filenames
    !-----------------------
    fg_file       = path_file (input,  fg_file)
    oldanerr_file = path_file (input,  oldanerr_file)
    ana_err_file  = path_file (output, ana_err_file)
    ana_file      = path_file (output, ana_file)
!   out_file      = path_file (output, out_file)
    ready_det     = path_file (output, ready_det)
    invar_det     = path_file (input,  invar_det)
    grid_file     = path_file (input,  grid_file)
    file_strato   = path_file (input,  file_strato)
    grid_file_strato= path_file (input, grid_file_strato)
    !---------------------------------------------------------------------
    ! for 3dvar-PSAS: check if files are readable, remove old output files
    !---------------------------------------------------------------------
    chk = (method(1:4) == 'PSAS'  .or. method(1:5) == 'ENVAR')
    if (present(check)) chk = check
    if (dace% lpio .and. chk) then
      iunit = get_unit_number()

      open (iunit, file=fg_file, status='old', action='read', iostat=ios)
      if (ios/=0) call finish ('nml_run_flags','cannot open '//trim(fg_file))
      close (iunit)

      if (oldanerr_file/='/none/') then
        open (iunit, file=oldanerr_file, status='old',action='read',iostat=ios)
        if (ios/=0) call finish ('nml_run_flags',&
                                 'cannot open '//trim(oldanerr_file))
        close (iunit)
      endif

      open (iunit, file=ana_err_file, action='write', iostat=ios)
      if (ios/=0) call finish ('nml_run_flags',&
                               'cannot write '//trim(ana_err_file))
      close (iunit, status='delete')

      open (iunit, file=ana_file, action='write', iostat=ios)
      if (ios/=0) call finish ('nml_run_flags',&
                               'cannot write '//trim(ana_file))
      close (iunit, status='delete')

      if (ready_det /= "") then
        open (iunit, file=ready_det, action='write', iostat=ios)
        if (ios/=0) call finish ('nml_run_flags',&
                                 'cannot write '//trim(ready_det))
        close (iunit, status='delete')
      end if

      call return_unit_number (iunit)
    endif
    !-------------------------------
    ! check if linked with GRIB2 API
    !-------------------------------
#ifndef GRIB_API
    grib_edition   =  1   ! GRIB edition to write
    grib_library   =  1   ! GRIB API to use: 1=(C)GRIBEX  2=GRIB2-API
#else
    call set_grib (edition=grib_edition, library=grib_library)
#endif
    call grib_api_version (major, minor, rev)
    !------------------------------
    ! retrieve compiler information
    !------------------------------
    call compiler_info (vendor=vendor, version=ftn_version)
    !---------
    ! printout
    !---------
    if (dace% lpio) then
      write (6,'(a       )') repeat('-',79)
      write (6,'(        )')
      write (6,'(a       )') '  DWD 3D-Var/LETKF'
      write (6,'(        )')
      write (6,'(a       )') '  Version: '//var3d_version()               &
                                          //var3d_date()//' '//var3d_time()
      write (6,'(        )')
      write (6,'(a,a     )') '  host               = ',trim (host)
!     do pe = 0, dace% npe-1
!       write (6,'(a,a,i5)') '  host (pe)      = ',hosts (pe), pe
!     end do
    end if
    call print_nodes ()
    if (dace% lpio) then
      write (6,'(        )')
      write (6,'(a,a     )') '  run time (gmt)     = ',cyyyymmddhhmm(run_time)
      write (6,'(a,a     )') '  user               = ',user
      write (6,'(a,i10   )') '  center             = ',center
      write (6,'(a,i10   )') '  subcenter          = ',subcenter
      write (6,'(a,i10   )') '  process (analysis) = ',proc_ana
      write (6,'(a,i10   )') '    (analysis error) = ',proc_ana_err
      write (6,'(a,i10   )') ' (analysis ensemble) = ',proc_ana_ens
      write (6,'(a,i10   )') '  range_ana          = ',range_ana
      write (6,'(a,i10   )') '  ensemble_id        = ',ensemble_id
      write (6,'(        )')
      write (6,'(a,i10   )') '  nex                = ',nex
      write (6,'(a,i10   )') '  run_type           = ',run_type
      write (6,'(a,a     )') '  runtype            = ',runtype
      write (6,'(a,a     )') '  method             = ',method
      write (6,'(a,a     )') '  model              = ',model
      write (6,'(a,i10   )') '  cosmo_refatm       = ',cosmo_refatm
      write (6,'(a,i10   )') '  cosmo_ivctype      = ',cosmo_ivctype
      write (6,'(a,i10   )') '  cosmo_nfltvc       = ',cosmo_nfltvc
      write (6,'(a,2f11.0)') '  cosmo_svc          = ',cosmo_svc
      write (6,'(a,l1    )') '  l_ke_in_gds        = ',l_ke_in_gds
      write (6,'(        )')
      write (6,'(a,i10   )') '  nprocs             = ',dace% npe
      write (6,'(a,i10   )') '  nproc1             = ',nproc1
      write (6,'(a,i10   )') '  nproc2             = ',nproc2
      write (6,'(a,i10   )') '  nproc_repro        = ',nproc_repro
      write (6,'(a,l1    )') '  barrier            = ',barrier
      write (6,'(a,l1    )') '  check_sync         = ',check_sync
      write (6,'(        )')
      if (mod (fc_ref_time% secs, 3600) == 0 .and. &
          mod (ana_time%    secs, 3600) == 0       ) then
       write(6,'(a,a     )') '  reference time     = ',cyyyymmddhh  (fc_ref_time)
       write(6,'(a,a     )') '  analysis  time     = ',cyyyymmddhh  (ana_time)
      else
       write(6,'(a,a     )') '  reference time     = ',cyyyymmddhhmm(fc_ref_time)
       write(6,'(a,a     )') '  analysis  time     = ',cyyyymmddhhmm(ana_time)
      end if
      write (6,'(a,a     )') '  analysis  date     = ',cyyyymmddhh(ana_date)
      write (6,'(a,a     )') '  observ. time shift =      ',chhhmm(obs_timeshift)
      write (6,'(a,i2    )') '  otime_include      = ',otime_include
      write (6,'(a,f5.2  )') '  forecast hours     =      ',fc_hours
      write (6,'(        )')
      write (6,'(a,a     )') '  data     path      = ',trim(data)
      write (6,'(a,a     )') '  input    path      = ',trim(input)
      write (6,'(a,a     )') '  output   path      = ',trim(output)
      write (6,'(a,a     )') '  aux      path      = ',trim(aux)
      write (6,'(a,a     )') '  obsinput path      = ',trim(obsinput)
      write (6,'(        )')
      write (6,'(a,a     )') '  fg_file            = ',trim(fg_file)
      write (6,'(a,a     )') '  invar_det          = ',trim(invar_det)
      write (6,'(a,a     )') '  grid_file          = ',trim(grid_file)
      write (6,'(a,a     )') '  oldanerr_file      = ',trim(oldanerr_file)
      write (6,'(a,a     )') '  ana_err_file       = ',trim(ana_err_file)
      write (6,'(a,a     )') '  ana_file           = ',trim(ana_file)
!     write (6,'(a,a     )') '  out_file           = ',trim(out_file)
      write (6,'(a,a     )') '  ready_det          = ',trim(ready_det)
      write (6,'(a,a     )') '  fdbk_basename      = ',trim(fdbk_basename)
      write (6,'(a,a     )') '  pass_fields        = ',trim(pass_fields)
      write (6,'(a,a     )') '  read_fields        = ',trim(read_fields)
      write (6,'(a,a     )') '  opt_fields         = ',trim( opt_fields)
      write (6,'(a,a     )') '  sst_fields         = ',trim( sst_fields)
      write (6,'(a,a     )') '  urun_fields        = ',trim(urun_fields)
      write (6,'(a,a     )') '  dealloc_fields     = ',trim(dealloc_fields)
      write (6,'(a,a     )') '  gribout_24bit      = ',trim(gribout_24bit)
      write (6,'(a,a     )') '  file_strato        = ',trim(file_strato)
      write (6,'(a,a     )') '  grid_file_strato   = ',trim(grid_file_strato)
      write (6,'(a,i0    )') '  interp_strato      = ',interp_strato
      write (6,'(a,l1    )') '  abort_strato       = ',abort_strato
      write (6,'(a,l1    )') '  read_obs_first     = ',read_obs_first
      write (6,'(a,l1    )') '  nature_run         = ',nature_run
      write (6,'(a,i0    )') '  chk_icon_grid      = ',chk_icon_grid
      write (6,'(a,i0    )') '  debug_rss          = ',debug_rss
      write (6,'(a,i0    )') '  expseed            = ',expseed
      if (io_max_gather  >= 0) &
       write(6,'(a,i0    )') '  io_max_gather      = ',io_max_gather
      if (io_max_scatter >= 0) &
       write(6,'(a,i0    )') '  io_max_scatter     = ',io_max_scatter
      write (6,'(        )')
      write (6,'(a,a,1x,a)') '  compiler & version = ',trim(vendor), &
                                                       trim(ftn_version)
      write (6,'(a,i2    )') '  grib_edition       = ',grib_edition
      write (6,'(a,i2    )') '  grib_library       = ',grib_library
      write (6,'(a,i2    )') '  gribtabversion     = ',gribtabversion
      write (6,'(a,i2,2(".",i0))') &
                             '  grib_api version   = ',major, minor, rev
      write (6,'(a,a     )') '  grib_packing       = ',trim(grib_packing)
      write (6,'(a,a     )') '  geoid_format       = ',trim(geoid_format)
      write (6,'(a,a     )') '  geoid_file         = ',trim(geoid_file)
      write (6,'(a,a5,i2 )') '  rttov packages     : '           &
#if defined (_RTTOV_VERSION)
                                                      ,"RTTOV", _RTTOV_VERSION &
#endif
                                                      ;
      write (6,'()')
      write (6,'(a,i10,i4)') '  ga3_biasc_rad      = ',ga3_biasc_rad,   &
                                                      flag_biasc_rad
      write (6,'(a,i10,i4)') '  ga3_biasc_airep    = ',ga3_biasc_airep, &
                                                      flag_biasc_airep
      write (6,'(a,i10,i4)') '  ga3_biasc_synop    = ',ga3_biasc_synop, &
                                                      flag_biasc_synop
      write (6,'(a,i10,i4)') '  ga3_biasc_gpsgb    = ',ga3_biasc_gpsgb, &
                                                      flag_biasc_gpsgb
      write (6,'(a,i10,i4)') '  ga3_biasc_scatt    = ',ga3_biasc_scatt, &
                                                      flag_biasc_scatt
      write (6,'(a,i10,i4)') '  ga3_biasc_wlidar   = ',ga3_biasc_wlidar, &
                                                      flag_biasc_wlidar
      write (6,'(        )')
      write (6,'(a       )') repeat('-',79)
    endif
    if (allocated(hosts)) deallocate (hosts)

#ifdef GRIB_API
    select case (grib_packing)
    case ("","grid_simple","grid_second_order")
    case ("grid_ccsds")
      if (major * 100 + minor < 204)                                         &
           call finish ("nml_run_flags",                                     &
                        "ECCODES >= 2.4 required for " // trim (grib_packing))
    case default
      call finish ("nml_run_flags","unsupported packing: "//trim (grib_packing))
    end select
#endif

  end subroutine nml_run_flags
!------------------------------------------------------------------------------
  elemental function path_file (pathname, filename, suffix, iens, &
                                anatime, fcreftime, leadtime      )
  character(len=*) ,intent(in)           :: pathname  ! path
  character(len=*) ,intent(in)           :: filename  ! file
  character(len=*) ,intent(in) ,optional :: suffix    ! file suffix
  integer          ,intent(in) ,optional :: iens      ! ensemble member
  type(t_time)     ,intent(in) ,optional :: anatime   ! analyis time
  type(t_time)     ,intent(in) ,optional :: fcreftime ! forecast reference time
  type(t_time)     ,intent(in) ,optional :: leadtime  ! forecast lead time
  character(len=256)                     :: path_file ! concatenated pathname
  !--------------------------------------------------------------------------
  ! This routine concatenates path, basename and optionally suffix to derive
  !  a complete file pathname. Special character strings within the pathname
  !  are replaced.
  ! If filename does not begin with '/' or '.'
  !   pathname and filename are concateneted
  !   a '/' is added in between if required.
  ! If filename contains special character strings these strings are
  !  replaced by the analysis, forecast reference time, etc.
  ! If suffix is present the string is appended to the filename.
  ! If filename contains _ENS_ this string is replaced by the (3digit)
  !  ensemble member id if present. Otherwise iens is appended to the
  !  file name as a suffix.
  ! The following special character strings are replaced:
  !
  !  _ANDATE_        analysis date             yyyymmdd
  !  _ANA_TIME_      analysis time             yyyymmddhh
  !  _ANA_TIMEMM_    analysis time             yyyymmddhhmm
  !  _ANA_TIMEMMSS_  analysis time             yyyymmddhhmmss
  !  _FCR_TIME_      forecast reference time   yyyymmddhh
  !  _FCR_TIMEMM_    forecast reference time   yyyymmddhhmm
  !  _FCR_TIMEMMSS_  forecast reference time   yyyymmddhhmmss
  !  VVVMM           forecast lead time               hhhmm
  !  VVVMMSS         forecast lead time               hhhmmss
  !  DDVVMMSS        forecast lead time              ddhhmmss
  !  DDDVVMMSS       forecast lead time             dddhhmmss
  !  _EXP_           experiment id             iii
  !
  ! Values for analysis time, forecast reference time, forecast lead time are
  !  taken from namelist /run/ or if present from the respective optional
  !  parameters.
  !--------------------------------------------------------------------------

    character(len=256) :: fname   ! file name
    integer            :: i       ! index in character string
    integer            :: j       ! index to iterate replacements
    type(t_time)       :: t_ana   ! analysis time to use
    type(t_time)       :: t_fcref ! forecast reference time to use
    type(t_time)       :: t_lead  ! forecast lead time to use
    !---------------------------------
    ! for empty argument return blanks
    !---------------------------------
    if (filename == '') then
      path_file = ''
      return
    endif
    !------------------------------
    ! concatenate path and basename
    !------------------------------
    if (pathname      /= ''  .and. &
        filename(1:1) /= '/' .and. &
        filename(1:1) /= '.'       ) then
      i = len_trim (pathname)
      if (pathname(i:i) == '/') then
        fname = trim(pathname)//filename
      else
        fname = trim(pathname)//'/'//filename
      endif
    else
      fname = filename
    endif
    !----------------------------------------------------
    ! set times to use from optional argument or namelist
    !----------------------------------------------------
    t_ana   = ana_time;    if (present (anatime))   t_ana   = anatime
    t_fcref = fc_ref_time; if (present (fcreftime)) t_fcref = fcreftime
    t_lead  = fc_time;     if (present (leadtime))  t_lead  = leadtime
    !-------------------------------------------
    ! replace date-string or exp.id. in filename
    ! try 2 times
    !-------------------------------------------
    do j = 1, 2
      i     = index (fname, '_ANDATE_')
      if (i /= 0) fname (i:i+7)  = cyyyymmdd       (t_ana)
      i     = index (fname, '_ANA_TIME_')
      if (i /= 0) fname (i:i+9)  = cyyyymmddhh     (t_ana)
      i     = index (fname, '_ANA_TIMEMM_')
      if (i /= 0) fname (i:i+11) = cyyyymmddhhmm   (t_ana)
      i     = index (fname, '_ANA_TIMEMMSS_')
      if (i /= 0) fname (i:i+13) = cyyyymmddhhmmss (t_ana)
      i     = index (fname, '_FCR_TIME_')
      if (i /= 0) fname (i:i+9)  = cyyyymmddhh     (t_fcref)
      i     = index (fname, '_FCR_TIMEMM_')
      if (i /= 0) fname (i:i+11) = cyyyymmddhhmm   (t_fcref)
      i     = index (fname, '_FCR_TIMEMMSS_')
      if (i /= 0) fname (i:i+13) = cyyyymmddhhmmss (t_fcref)
      i     = index (fname, '_EXP_')
      if (i /= 0) fname (i:i+4)  = char5           (nex)
      i     = index (fname, 'DDDVVMMSS')
      if (i /= 0) fname (i:i+8)  = cdddhhmmss      (t_lead)
      i     = index (fname, 'DDVVMMSS')
      if (i /= 0) fname (i:i+7)  = cddhhmmss       (t_lead)
      i     = index (fname, 'VVVMMSS')
      if (i /= 0) fname (i:i+6)  = chhhmmss        (t_lead)
      i     = index (fname, 'VVVMM')
      if (i /= 0) fname (i:i+4)  = chhhmm          (t_lead)
    end do
    !--------------------------
    ! append suffix to basename
    !--------------------------
    if (present(suffix)) then
      i = len_trim (fname)
      fname (i+1:) = suffix
    endif
    !---------------------------------
    ! add / replace ensemble member id
    !---------------------------------
    if (present (iens)) then
      i     = index (fname, '_ENS_')
      if (i /= 0) then
        fname (i:i+3) = char3 (iens)
        fname (i+4:)  = fname (i+6:)
      else
        fname = trim(fname)//'.'//char3 (iens)
      endif
    endif
    path_file = fname
  end function path_file
!------------------------------------------------------------------------------
  elemental function basename (pathname, suffix)
    character(len=*)             ,intent(in) :: pathname  ! full path
    character(len=*) ,optional   ,intent(in) :: suffix    ! suffix pattern
    character(len=len(pathname))             :: basename
    !------------------------------------------------------
    ! Derive filename with any leading directories removed,
    ! optionally removing trailing suffix
    !------------------------------------------------------
    integer :: i, j
    i = index (pathname, '/', back=.true.)      ! Find last '/'
    j = 0
    if (present (suffix)) then
       if (suffix /= "") j = index (pathname, suffix, back=.true.) - 1
    end if
    if (j <= i) j = len_trim (pathname)
    basename = pathname(i+1:j)
  end function basename
!------------------------------------------------------------------------------
  subroutine process_grid (grid, ni, ensemble)
  integer ,intent(in) :: grid      ! WMO grid id
  integer ,intent(in) :: ni        ! resolution parameter
  logical ,intent(in) :: ensemble  ! flag for ensemble analysis
  !----------------------------------------------------
  ! derive generating process from grid specification.
  ! set parameters proc_ana, proc_ana_err, proc_ana_ens
  ! if not given by the namelist.
  !----------------------------------------------------
    integer :: p_ana, p_err
    !-----------------------------------------------
    ! return if parameters are given by the namelist
    !-----------------------------------------------
    if (ensemble) then
      if (proc_ana_ens /= 0)                     return
    else
      if (proc_ana /= 0 .and. proc_ana_err /= 0) return
    endif
    !----------------------------------
    ! derive process id from resolution
    !----------------------------------
    select case (grid)
    case (DWD6_ICON)
       p_ana = 1        ! ICON global icosahedral grid (resolution-independent)
       p_err = 26       ! ICON analysis error
    case (DWD6_ICOSAHEDRON)
       select case (ni)
       case (16)
          p_ana = 171
          p_err = 189
       case (32)
          p_ana = 141
          p_err = 190
       case (48)
          p_ana = 143
          p_err = 191
       case (64)
          p_ana = 145
          p_err = 192
       case (96)
          p_ana = 147
          p_err = 193
       case (128)
          p_ana = 149
          p_err = 194
       case (192)
          p_ana = 173
          p_err = 195
       case (256)
          p_ana = 175
          p_err = 196
       case (384)
          p_ana = 206
          p_err = 210
       case default
          write(0,*)  'process_grid:  invalid ni =',ni
          call finish('process_grid','invalid ni')
       end select
    case (WMO6_ROTLL)
       ! Temporary logic, may need revision
       p_ana = 137                                      ! COSMO-DE analysis
       p_err =   0
       if (proc_ana     == 0) proc_ana     = p_ana
       if (proc_ana_ens == 0) proc_ana_ens = proc_ana
    case (WMO6_GAUSSIAN)
       select case (ni)
       case (16)
          p_ana = 177
          p_err = 189
       case (32)
          p_ana = 161
          p_err = 190
       case (48)
          p_ana = 163
          p_err = 191
       case (64)
          p_ana = 165
          p_err = 192
       case (96)
          p_ana = 167
          p_err = 193
       case (128)
          p_ana = 169
          p_err = 194
       case (192)
          p_ana = 179
          p_err = 195
       case (256)
          p_ana = 181
          p_err = 196
       case default
          write(0,*)  'process_grid:  invalid ng =',ni
          call finish('process_grid','invalid ng for Gaussian grid')
       end select
    case default
      write(0,*)  'process_grid:  no process defined for grid =',grid
      call finish('process_grid','no process defined for grid')
    end select
    !------------------------------------------------------
    ! set parameters if not given by the namelist, printout
    !------------------------------------------------------
    if (ensemble) then
      if (proc_ana_ens == 0) proc_ana_ens = p_ana
      if(dace% lpio) then
        write(6,'(a)')
        write(6,'(a,i5)') '  proc_ana_ens =', proc_ana_ens
        write(6,'(a)')
      endif
    else
      if (proc_ana     == 0) proc_ana     = p_ana
      if (proc_ana_err == 0) proc_ana_err = p_err
      if(dace% lpio) then
        write(6,'(a)')
        write(6,'(a,i5)') '  process      =', proc_ana
        write(6,'(a,i5)') '  proc_ana_err =', proc_ana_err
        write(6,'(a)')
      endif
    endif
  end subroutine process_grid
!------------------------------------------------------------------------------
  subroutine set_nprocs (nproc1, nproc2, nprocs)
  integer ,intent(out)          :: nproc1 ! no.of processors in direction 1
  integer ,intent(out)          :: nproc2 ! no.of processors in direction 2
  integer ,intent(in) ,optional :: nprocs ! total number of processors
  !------------------------------------------------------------------------
  ! set nproc1,nproc2 (number of processors for model domain decomposition)
  ! for the case that nproc1,nproc2 are not specified in the namelist
  ! so that nproc1 * nproc2 = nprocs
  !------------------------------------------------------------------------
    integer :: ie, idiff, isqrt, np1, np2, n
    ie     = huge(ie)
    n      = dace% npe; if (present(nprocs)) n = nprocs
    isqrt  = nint (sqrt (real(n)))
    do np1 = isqrt,1,-1
      np2 = n / np1
      idiff  = n - np1 * np2
      if (idiff < 0) cycle
      if (idiff < ie) then
        ie     = idiff
        nproc1 = np1
        nproc2 = np2
        if (ie == 0) exit
      endif
    end do
  end subroutine set_nprocs
!------------------------------------------------------------------------------
  subroutine get_hostnames ()
    !--------------
    ! get hostnames
    !--------------
    integer :: pe
    call gethostname (host)
    allocate (hosts (0 : dace% npe-1))
    if (p_parallel) then
      do pe = 0, dace% npe-1
        hosts(pe) = host
        call p_bcast (hosts(pe), pe)
      end do
    else
      hosts = host
    endif
  end subroutine get_hostnames
!------------------------------------------------------------------------------
  subroutine get_nodes ()
    !-----------------------------------------
    ! Derive association of MPI ranks to nodes
    ! (identified by their hostnames)
    !-----------------------------------------
    integer            :: i                     ! loop index
    integer            :: n                     ! # tasks on node
    integer            :: pe                    ! first rank on node
    integer            :: node                  ! Logical node number
    character(len=16)  :: nodename

    if (.not. allocated (hosts)) call get_hostnames ()

    if (allocated (node_id)) return

    allocate (node_id(0:dace% npe-1))
    node_id = 0

    !-------------------------------------------------
    ! Initialize: associate rank 0 with logical node 1
    !-------------------------------------------------
    node        = 1
    pe          = 0
    node_loop: do
       node_id(pe) = node
       nodename    = hosts(pe)
       n           = 1
       do i = pe+1, dace% npe-1
          if (node_id(i) > 0) cycle
          if (hosts(i) == nodename) then
             n          = n + 1
             node_id(i) = node
          end if
       end do
       !-----------------------------
       ! Search for next logical node
       !-----------------------------
       do i = pe+1, dace% npe-1
          if (node_id(i) == 0) then
             node = node + 1
             pe   = i
             cycle node_loop
          end if
       end do
       exit node_loop
    end do node_loop
  end subroutine get_nodes
!------------------------------------------------------------------------------
  subroutine print_nodes ()
    !----------------------------------------
    ! Print association of MPI ranks to nodes
    ! (identified by their hostnames)
    !----------------------------------------
    integer            :: i                     ! loop index
    integer            :: n                     ! # tasks on node
    integer            :: pe                    ! first rank on node
    integer            :: node                  ! Logical node number
    integer            :: node_id(0:dace% npe-1)
    integer, parameter :: nline = 10            ! max. ranks/line
    integer            :: rank   (dace% npe)
    character(len=16)  :: nodename

    if (.not. allocated (hosts)) call get_hostnames ()
    node_id = 0

    if (dace% lpio) then
       write(6,'()')
       write(6,'(2x,a,8x,a4,3x,a)') "hostname", "node","MPI ranks:"
       !-------------------------------------------------
       ! Initialize: associate rank 0 with logical node 1
       !-------------------------------------------------
       node        = 1
       pe          = 0
       node_loop: do
          node_id(pe) = node
          nodename    = hosts(pe)
          n           = 1
          rank(n)     = pe
          do i = pe+1, dace% npe-1
             if (node_id(i) == 0) then
                if (hosts(i) == nodename) then
                   n          = n + 1
                   rank(n)    = i
                   node_id(i) = node
                end if
             end if
          end do
          do i = 1, n, nline
             if (i == 1) then
                write(6,'(2x,a16,i4," : ",99i5)') &
                     nodename, node,              &
                     rank(i:min(i+nline-1,n))
             else
                write(6,'(25x,99i5)') &
                     rank(i:min(i+nline-1,n))
             end if
          end do
          !-----------------------------
          ! Search for next logical node
          !-----------------------------
          do i = pe+1, dace% npe-1
             if (node_id(i) == 0) then
                node = node + 1
                pe   = i
                cycle node_loop
             end if
          end do
          exit node_loop
       end do node_loop
    end if
  end subroutine print_nodes
!------------------------------------------------------------------------------
  subroutine cleanup_run_params
  !---------------------
  ! clean up this module
  !---------------------
    if (allocated(hosts  )) deallocate (hosts)
    if (allocated(node_id)) deallocate (node_id)
  end subroutine cleanup_run_params
!------------------------------------------------------------------------------
end module mo_run_params
