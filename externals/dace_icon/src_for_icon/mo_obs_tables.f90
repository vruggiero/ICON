!
!+ Tables with characteristics of the report types
!
MODULE mo_obs_tables
!
! Description:
!   This module holds tables characterising the report types and statistics
!   on their usage in the 3D-Var and LETKF:
!
!   type (t_rept_char) t_rept_char:
!     Characteristics of the report type.
!     Associated module type and Datenbankkennzahlen.
!
!   type (t_rept_stat) rept_stat,
!   dism_stat, pass_stat, reje_stat:
!     Report type usage statistics
!
!   type (t_rept_stat) obsv_stat:
!     Observation variable usage statistics
!
!   type (t_rept_use) rept_use:
!     Observation type usage and selection flags (namelist /REPORT/)
!
! Current Maintainer: DWD, Harald Anlauf, Alexander Cress
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Harald Anlauf
!  read_nml_checks: fix error message; Bump npr to 100000
! V1_4         2009/03/26 Andreas Rhodin
!  Option to read NetCDF files in parallel
! V1_5         2009/05/25 Harald Anlauf
!  write_pending, write_head: optimize for NEC SX-9
!  check_report_1: set CHK_DOMAIN for bad horizontal coordinates
!  introduce consistency check for status(CHK_DOMAIN) <= STAT_DISMISS
! V1_7         2009/08/24 Andreas Rhodin
!  add DBKZ=10385 (BUOY, new BUFR format)
! V1_8         2009/12/09 Andreas Rhodin
!  make processing of use_cntr, moni_cntr more intuitive
! V1_9         2010/04/20 Andreas Rhodin
!  account for status flag STAT_DEFAULT=0 .
!  option to exclude are near boundaries in out of domain check
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  Module variable rept_char: add known dbkz numbers for GPSRO
! V1_12        2010/09/03 Alexander Cress
!  Activate CODAR reports (dbkz 528, for Antarctic balloon experiment)
! V1_13        2011/11/01 Andreas Rhodin
!  changes for oceansat-2, windprofilers
!  extend statistics for used variables
!  implement checks on cutoff time
! V1_15        2011/12/06 Andreas Rhodin
!  implement data structures to trace invalid observation operators
! V1_19        2012-04-16 Andreas Rhodin
!  changes for RADAR observation operator (COSMO LETKF)
! V1_20        2012-06-18 Andreas Rhodin
!  for LETKF: adjust time_b,time_e defaults; dont check for "synop" time
! V1_22        2013-02-13 Andreas Rhodin
!  changes for wind profilers, Slant Total Delay operater
!  revised check for data base time < analysis window (20 min hardcoded)
!  write comments for REJECTED reports, don't write double entries
! V1_23        2013-03-26 Andreas Rhodin
!  define module type RADAR
! V1_26        2013/06/27 Andreas Rhodin
!  replace STAT_ACTIVE_0 by STAT_NOTACTIVE
! V1_28        2014/02/26 Andreas Rhodin
!  make update_statistics public
! V1_29        2014/04/02 Alexander Cress
!  Handle MODE-S data
! V1_31        2014-08-21 Andreas Rhodin
!  fix bug for GMESTAT feedback file input (introduced in V1_29)
! V1_35        2014-11-07 Alexander Cress
!  changes for TEMP SHIP BUFR reports
! V1_36        2014-11-13 Gerhard Paul
!  changes for new SYNOP BUFR reports
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_43        2015-08-19 Harald Anlauf
!  print_rept_stat: widen format for observation statistics
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_45        2015-12-15 Andreas Rhodin
!  MEC: generalise 'check_report'; COSMO LETKF: set CHK_OPERATOR to PASSIVE
! V1_47        2016-06-06 Andreas Rhodin
!  set ZTD/STD passive by default, less restrictive quality check for COSMO-MEC
! V1_48        2016-10-06 Andreas Rhodin
!  handle shift of assimilation interval in ICON/COSMO
! V1_51        2017-02-24 Andreas Rhodin
!  chanhes for COMET (OT_SOIL); fixes for GPSRO and SYNOP (version=1) operators
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  2004-2008  original source
! Gerhard Paul    2008       add SKY Datenbankkennzahlen
! Harald Anlauf   2008       various fixes
!-------------------------------------------------------------------
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!============================================================================

  !-------------
  ! modules used
  !-------------
  use mo_kind,      only: wp                 ! kind parameters
  use mo_mpi_dace,  only: dace,             &! MPI group info
                          p_bcast,          &! broadcast routine
                          p_send,           &! generic MPI send routine
                          p_recv,           &! generic MPI receive routine
                          p_sum,            &! generic MPI sum routine
                          p_gather           ! generic MPI gather routine
!                         p_gatherv          ! generic MPI gatherv routine
  use mo_exception, only: finish             ! abort routine
  use mo_namelist,  only: position_nml,     &! routine to position nml group
                          nnml,             &! namelist fortran unit number
                          POSITIONED         ! position_nml: OK    return flag
  use mo_obstypes,  only: obstype_dbkz,     &! BUFR/CMA types from DWD-dbkz
                          t_obsid            ! derived type for table entry
  use mo_t_obs,     only: t_obs,            &! observation derived type
                          t_spot,           &! report meta data type
                          t_head,           &! component of t_spot
                          invalid,          &! invalid observation value
                          SYNOP, TEMP, TOVS,&! observation module numbers
                          GPSRO, AMV, AIREP,&!
                          GPSGB, RADAR,     &!
                          COSMO, WLIDAR,    &!
                          SCATT, SOIL, MWR, &!
                          dace_op,          &! use dace obs. operators (MEC)
                          n_dace_op          ! length of non-empty entries in dace_op
  use mo_wmo_tables,only: WMOA_SURF_LAND,   &! Surface data - land
                          WMOA_SURF_SEA,    &! Surface data - sea
                          WMOA_VERT_SOUND,  &! Vertical soundings
                          WMOA_SINGL_SAT,   &! Single level upper-air data
                          WMOA_SINGL_LEV     ! Single level upper-air data
  use mo_dwd_tables,only: dbkz_mnem          ! get description from dbkz
  use mo_t_use,     only: t_use,            &! data type to hold report state
                          decr_use,         &! decrease the state of a report
                          n_chk,            &!
                          chk_key,          &! id for 'checks'
                          stat_key,         &! id for 'states'
                          stat_mnem,        &! mnemonics for 'states':
                          n_stat,           &! number of status values
                          STAT_ABORT,       &!
                          STAT_FORGET,      &!
                          STAT_DISMISS,     &!
                          STAT_MERGED,      &!
                          STAT_OBS_ONLY,    &!
                          STAT_NOTACTIVE,   &!
                          STAT_ACTIVE,      &!
                          STAT_ACTIVE_0I,   &!
                          STAT_ACTIVE_1,    &!
                          STAT_ACCEPTED,    &!
                          STAT_PASSIVE,     &!
                          STAT_PAS_REJ,     &!
                          STAT_REJECTED,    &!
                          STAT_DEFAULT,     &!
                          chk,              &! mnem. + states for 'checks':
                          CHK_NONE, CHK_TIME, CHK_AREA, CHK_HEIGHT, &
                          CHK_THIN, CHK_RULE, CHK_SUBTYP,           &
                          CHK_DOMAIN, CHK_NOTUSED, CHK_INSDAT,      &
                          CHK_OPERATOR, CHK_SURF
  use mo_fdbk_tables,                       &!
                    only: OT_SYNOP,         &! observation type values:
                          OT_AIREP,         &!
                          OT_SATOB,         &!
                          OT_DRIBU,         &!
                          OT_TEMP,          &!
                          OT_PILOT,         &!
                          OT_SATEM,         &!
                          OT_PAOB,          &!
                          OT_SCATT,         &!
                          OT_RAD,           &!
                          OT_GPSRO,         &!
                          OT_GPSGB,         &!
                          OT_RADAR,         &!
                          OT_POWER,         &!
                          OT_SOIL,          &!
                          OT_OBJECT,        &!
                          OT_LIGHTN,        &!
                          OT_WLIDAR,        &!
                          OT_MWR,           &!
                          n_ot,             &! number of obs.types
                          n_vn,             &! number of varnos
                          varno,            &! table of variable numbers
                          init_fdbk_tables   ! initialise the tables
  use mo_obs_rules, only: get_rule           ! get observation processing rules
  use mo_time,      only: t_time,           &! time data type
                          operator(+),      &! add times
                          operator(-),      &! subtract times
                          operator(/),      &! divide time
                          operator(<),      &! compare times
                          operator(>),      &! compare times
                          operator(<=),     &! compare times
                          operator(>=),     &! compare times
                          operator(/=),     &! compare times
                          cyyyymmddhh,      &! derive 'yyyymmddhh'
                          cyyyymmddhhmm,    &! derive 'yyyymmddhhmm'
                          chhmmss,          &! derive 'hhmmss' from time
                          ihh,              &! derive hh (hours) from time
                          init_time,        &! routine to set time variable
                          zero_time,        &! default initialisation value
                          invalid_time       ! invalid time value
  use mo_run_params,only: ana_time,         &! analysis time
                          run_time,         &! time of the model run
                          fc_ref_time,      &! reference time (forecast start)
                          obs_timeshift,    &! shift of observation time window
                          otime_exclude,    &! T for compatibility < r13181
                          aux,              &! path for auxiliary output
                          path_file,        &! concatenate path/filename
                          host,             &! hostname
                          user,             &! user name
                          nex,              &! experiment number
                          run_type,         &! assimilation cycle, forecasts,..
                          fc_ref_time,      &! forecast start time
                          fc_hours,         &! forecast duration
                          method             ! '3dvar', 'LETKF',...
  use mo_usstd,     only: h_p_usstd          ! h(gpm) from p US std.atm.
  use mo_fortran_units,&
                    only: get_unit_number    ! get a free unit number
!                         return_unit_number ! release unit number after use
  implicit none

  !================
  ! public entities
  !================
  private
  !----------------------------
  ! report type characteristics
  !----------------------------
  public :: rept_char          ! report type characteristics table
  public :: t_rept_char        ! table entry type definition
  public :: idb_dbk            ! function: table index (Datenbankkennzahl)
  !-----------------------
  ! report type statistics
  !-----------------------
  public :: t_rept_stat        ! table entry type definition
  public :: rept_stat          ! observation type statistics table
  public :: rept_stat0         ! empty obs. type statistics table entry
  public :: dism_stat          ! table for statistics of dismissed reports
  public :: pass_stat          ! table for statistics of passive   reports
  public :: reje_stat          ! table for statistics of rejected  reports
  public :: oonly_stat         ! table for statistics of obsv.only reports
  public :: print_rept_stat    ! routine to print observation type statistics
  public :: derive_rept_stat   ! re-derive observation type statistics
  public :: gather_rept_stat   ! gather report type statistics on I/O PE
  !----------------------------------------------------
  ! report type usage and selection (namelist /REPORT/)
  !----------------------------------------------------
  public :: rept_use           ! report type usage table
  public :: t_rept_use         ! table entry type definition
  public :: read_nml_report    ! read namelist /REPORT/
  public :: print_rept_use     ! print result of /REPORT/
  public :: set_rept_use       ! set default values
  !------------------------------------
  ! perform checks on report types,
  ! update report status and statistics
  !------------------------------------
  public :: check_report_0     ! preset state variable, simple checks
  public :: check_report_1     ! standard checks
  public :: dism_report_0      ! dismiss report
  public :: decr_rpt_use       ! change use-flags of report
  public :: update_statistics  ! update report statistics
  !-----------------------------------------------------------
  ! write file on report type usage (ASCII: 1 line per report)
  !-----------------------------------------------------------
  public :: write_report       ! write report to status file
  public :: write_pending      ! write pending reports to status file
  !-------------------------------------------------------------
  ! observation types (to be revised for ground based GPS/Radar)
  !-------------------------------------------------------------
  public :: obstyp, t_obstyp   ! table, derived type for observation types

  !==========
  ! Constants
  !==========
  !---------------------------------------
  ! table of CMA observation operator type
  !---------------------------------------
  type t_obstyp
    character(len=8)  :: name ! mnemonic
    integer           :: key  ! key
    character(len=32) :: desc ! description
  end type t_obstyp

  type (t_obstyp) ,parameter :: obstyp(n_ot)= &
    (/ t_obstyp ('SYNOP'  ,OT_SYNOP ,''), &
       t_obstyp ('AIREP'  ,OT_AIREP ,''), &
       t_obstyp ('SATOB'  ,OT_SATOB ,''), &
       t_obstyp ('DRIBU'  ,OT_DRIBU ,''), &
       t_obstyp ('TEMP'   ,OT_TEMP  ,''), &
       t_obstyp ('PILOT'  ,OT_PILOT ,''), &
       t_obstyp ('SATEM'  ,OT_SATEM ,''), &
       t_obstyp ('PAOB '  ,OT_PAOB  ,''), &
       t_obstyp ('SCATT'  ,OT_SCATT ,''), &
       t_obstyp ('RAD'    ,OT_RAD   ,''), &
       t_obstyp ('GPSRO'  ,OT_GPSRO ,''), &
       t_obstyp ('GPSGB'  ,OT_GPSGB ,''), &
       t_obstyp ('RADAR'  ,OT_RADAR ,''), &
       t_obstyp ('POWER'  ,OT_POWER ,''), &
       t_obstyp ('SOIL'   ,OT_SOIL  ,''), &
       t_obstyp ('OBJECT' ,OT_OBJECT,''), &
       t_obstyp ('LIGHTN' ,OT_LIGHTN,''), &
       t_obstyp ('WLIDAR' ,OT_WLIDAR,''), &
       t_obstyp ('MWR'    ,OT_MWR   ,'') /)

  integer ,parameter :: ndb    = 40 ! number of Datenbankkennzahlen

  !==========================================================
  ! observation type characteristics  table entry  definition
  !==========================================================
  type t_rept_char
    !------------
    ! description
    !------------
    character(len=8) :: mnem       ! mnemonic
    character(len=64):: desc       ! verbal description
    !---------------
    ! identification
    !---------------
    integer          :: obstype    ! observation type
    integer          :: mod        ! f90 module number
    integer          :: bufrtyp(2) ! WMO BUFR types
    integer          :: dbk(ndb)   ! DWD Datenbankkennzahlen
    !----------------
    ! characteristics
    !----------------
  end type t_rept_char

  !=====================================================
  ! observation type statistics  table entry  definition
  !=====================================================
  type t_rept_stat
    integer :: processed = 0 ! no. reports processed
    integer :: merged    = 0 ! no. reports merged with others (TEMP A,B,C,D)
    integer :: obs_only  = 0 ! observation operator not applied
    integer :: dismissed = 0 ! no. reports dismissed (neither assi. nor monit.)
    integer :: passive   = 0 ! no. passive  reports (monitored)
    integer :: rejected  = 0 ! no. rejected reports (currently monitored)
    integer :: active    = 0 ! no. active   reports (assimilated)
    integer :: accepted  = 0 ! no. accepted reports (w_vqc > 0.5)
  end type t_rept_stat

  !============================================
  ! observation type use table entry definition
  !============================================
  integer, parameter :: nc        = 20 ! size of generating center array
  integer ,parameter ::  UKN_CNTR = -1 ! unknown (not specified in input)
  integer ,parameter ::  ANY_CNTR = -2 !
  integer ,parameter :: NONE_CNTR = -3 !

  integer ,parameter :: OU_NONE       = 0     !   not used
  integer ,parameter :: OU_LETKF      = 1     !   only store obs/fg values
  integer ,parameter :: OU_INP        = 2     !   need input values (for BC)
  integer ,parameter :: OU_FWD        = 3     !   only use forward operator
  integer ,parameter :: OU_VAR        = 4     !   use forward and Jakobian

  type t_rept_use
   integer      :: init           ! initialize (use) the respective module
   integer      :: use (1: n_chk) ! state values for specific checks
   integer      :: max_proc       ! max. number of processed observations
   integer      :: max_act        ! max. number of active observations
   integer      :: ni             ! resolution flag for simple thinning
   real(wp)     :: height_t       ! top    bound [hPa]
   real(wp)     :: height_b       ! bottom bound [hPa]
   real(wp)     :: height_top     ! top    bound [m]
   real(wp)     :: height_bot     ! bottom bound [m]
   real(wp)     :: lat_nb         ! northern bound
   real(wp)     :: lat_sb         ! southern bound
   real(wp)     :: lon_eb         ! eastern bound
   real(wp)     :: lon_wb         ! western bound
   real(wp)     :: excl_bnd       ! exclude observations at lateral boundaries
   real(wp)     :: min_dist       ! minimum distance between observations [m]
   real(wp)     :: minvdist       ! minimum vertical distance           [hPa]
   type(t_time) :: time_b         ! begin of valid time range
   type(t_time) :: time_e         ! end   of valid time range
   type(t_time) :: time_cutoff    ! cut-off time
   logical      :: lhour          ! assimilate at this hour of the day
   logical      :: time_s1        ! apply check to section 1 (synop) time
   real(wp)     :: fr_land        ! land fraction required
   real(wp)     :: fr_sea         ! sea  fraction required
   real(wp)     :: fr_noice       !      fraction required
   real(wp)     :: min_tsurf      ! minimum surface temperature allowed   [C]
   integer      :: deriv_p        ! derive p from h if not present
   integer      :: moni_cntr (nc) ! generating centers (monitoring)
   integer      :: use_cntr  (nc) ! generating centers (assimilation)
   logical      :: read1pe        ! read NetCDF file from 1 PE
   integer      :: inv_op (6)     ! observation operator validity check
  end type t_rept_use

  !==================================================
  ! observation type characteristics table definition
  !==================================================
  integer            :: i            ! used in implied do loop

  type (t_rept_char) :: rept_char (n_ot) = (/ &
    !  1:
    t_rept_char (  'SYNOP',&! mnemonic
                   'SYNOP',&! verbal description
                 OT_SYNOP ,&! observation type
                    SYNOP ,&! f90 module number
                 (/WMOA_SURF_LAND,WMOA_SURF_SEA/), &! WMO BUFR types
                 (/0,1,5,9,128,131,156,166,256,    &!
                   265,384,10000,10005,10128,      &! synop land manual new bufr, synop land autom. new bufr
                   10015,10143,10150,              &! synop land manual, autom. (WIGOS)
                   10158,                          &! cman coastal stations (USA)
                   170,10170,                      &! SWIS (Road weather stations)
                   10256,10384,                    &! synop ship manual new bufr, synop ship autom. new bufr
                   182,                            &! synop snow observations
                   (-1,i=24,ndb)/))                &! DWD Datenbankkennzahlen
                 ,         &!
    !  2:
    t_rept_char (  'AIREP',&! mnemonic
                   'AIREP',&! verbal description
                 OT_AIREP ,&! observation type
                    AIREP ,&! f90 module number
                 (/WMOA_SINGL_LEV,-1/)            ,&! WMO BUFR types
                 (/528,529,530,                    &! dbkz: Codar,Airep,Amdar
                   532,533,534,                    &! dbkz: ACARS: Germany LH; USA; Acar:UKmetoffice
                   535,538,                        &! dbkz: ACARS: China ; others(Korea,Japan,..)
                   10532,                          &! dbkz: ACARS: single-level (unified format)
                   10533,                          &! dbkz: TAMDAR and AFIRS
                   542,10534,                      &! dbkz: MODES
                  (-1,i=13,ndb)/))                 &! Datenbankkennz.
                 ,         &!
    !  3:
    t_rept_char (  'SATOB',&! mnemonic
                   'SATOB (AMV)',&! verbal description
                 OT_SATOB ,&! observation type
                    AMV   ,&! f90 module number
                 (/WMOA_SINGL_SAT,-1/)            ,&! WMO BUFR types
                 (/1672,1673,1674,1675,1677,       &! dbkz: wind(2),cloud,temp.upper air humidity
                   1704,                           &! dbkz: EUMETSAT                 geostationary sat.
                   1705,                           &! dbkz: NOAA/NESIDS              geostationary sat.
                   1706,                           &! dbkz: NOAA/NESIDS MODIS(ASCII) polar sat.
                   1707,                           &! dbkz: CHINA                    geostationary sat.
                   1708,                           &! dbkz: JAPAN                    geostationary sat.
                   1709,                           &! dbkz: NOAA/NESIDS NOAAxx       polar sat.
                   1710,                           &! dbkz: NOAA/NESIDS MODIS(BUFR)  polar sat.
                   1720,                           &! dbkz: EUMETSAT Sentinel-3 A/B  polar sat.
                   1711,1712,1713,1714,            &! dbkz: cloud,temp.,upper air humidity,ozone
                  (-1,i=18,ndb)/))&
                 ,         &!
    !  4:
    t_rept_char (  'DRIBU',&! mnemonic
                   'DRIBU (buoys and scatterometer)' ,&! verbal description
                 OT_DRIBU ,&! observation type
                    SYNOP ,&! f90 module number
                 (/WMOA_SURF_SEA,-1/), &! WMO BUFR types
                 (/385,10385,(-1,i=3,ndb)/)) &
                 ,         &!
    !  5:
    t_rept_char (  'TEMP' ,&! mnemonic
                   'TEMP' ,&! verbal description
                 OT_TEMP  ,&! observation type
                    TEMP  ,&! f90 module number
                 (/WMOA_VERT_SOUND,-1/),              &! WMO BUFR types
                 (/516, 517, 518, 519,                &! dbkz:TEMP      A,B,C,D mobile
                   520, 521, 522, 523,                &! dbkz:TEMP      A,B,C,D
                   524, 525, 526, 527,                &! dbkz:TEMP      A,B,C,D mobile
                   536,                               &! dbkz:TEMP      A,B,C,D merged <1983
                   10520, 10521,                      &! dbkz:TEMP BUFR
                   10526, 10527, 10574,               &! dbkz:TEMP BUFR high resolution
                   10776, 10777,                      &! dbkz:TEMP SHIP BUFR
                   10782, 10783, 10785,               &! dbkz:TEMP SHIP BUFR high resolution
                   10516, 10517, 10570,               &!      TEMP mobile BUFR
                   776, 777, 778, 779,                &!      TEMP SHIP A,B,C,D
                   792,                               &!      TEMP SHIP A,B,C,D merged <1983
                   780, 781, 782, 783, 10780,         &!      TEMP DROP A,B,C,D,BUFR
                   (-1,i=37,ndb)/))                   &!
                 ,         &!
    !  6:
    t_rept_char (  'PILOT',&! mnemonic
                   'PILOT',&! verbal description
                 OT_PILOT ,&! observation type
                    TEMP  ,&! f90 module number
                 (/WMOA_VERT_SOUND,-1/),              &! WMO BUFR types
                 (/508, 509, 510, 511,                &! dbkz:PILOT      A,B,C,D geopt.height
                   512, 513, 514, 515,                &! dbkz:PILOT      A,B,C,D
                   548, 549, 550, 551, 552,           &! dbkz:WINDPROF   USA
                   553, 554,                          &! dbkz:WINDPROF   u,v;t,w Lindenberg,Europe
                   555,                               &! dbkz:WINDPROF   Japan
                   556,                               &! dbkz:WINDPROF   RASS Deutschl
                   10553,                             &! dbkz:WINDPROF   RASS Deutschl
                   584,                               &! dbkz:SCADA, wind power turbines
                   10600,                             &! dbkz:VAD WIND PROFILE
                   764, 765, 766, 767,                &!      PILOT SHIP A,B,C,D geopt.height
                   768, 769, 770, 771,(-1,i=29,ndb)/))&!      PILOT SHIP A,B,C,D
                 ,         &!
    !  7:
    t_rept_char (  'SATEM',&! mnemonic
                   'SATEM (not implemented)',&! verbal description
                 OT_SATEM ,&! observation type
                    -1    ,&! f90 module number
                    -1    ,&! WMO BUFR types
                    -1    )&! DWD Datenbankkennzahlen
                 ,         &!
    !  8:
    t_rept_char (  'PAOB' ,&! mnemonic
                   'PAOB pressure observations derived from sat.',&!
                 OT_PAOB  ,&! observation type
                    SYNOP ,&! f90 module number
                  (/253,-1/),          &! WMO BUFR types
                  (/386,(-1,i=2,ndb)/))&! DWD Datenbankkennzahlen
                 ,         &!
    !  9:
    t_rept_char (  'SCATT',&! mnemonic
                   'Scatterometer', &! verbal description
                 OT_SCATT ,&! observation type
                    SCATT ,&! f90 module number
                  [WMOA_SURF_SEA,-1]                     ,&! WMO BUFR types
                  [1697,1698,1699,1700,1770              ,&! Scatterometer
                   1701,1702,1780                        ,&! Altimeter
                   (-1,i=9,ndb)])                         &! DWD Datenbankkennzahlen
                 ,         &!
    ! 10:
    t_rept_char (  'RAD'   ,&! mnemonic
                   'Radiances',&! verbal description
                 OT_RAD    ,&! observation type
                    TOVS   ,&! f90 module number
                    -1     ,&! WMO BUFR types
                    -1    ) &! DWD Datenbankkennzahlen
                 ,          &!
    !-------------------------------------
    ! Up to No. 10: (ECMWF) CMA convention.
    ! 3D-Var local definitions follow:
    !-------------------------------------
    ! 11:
    t_rept_char (  'GPSRO',&! mnemonic
                   'GPS Radio occultations (bending angles)',&! description
                 OT_GPSRO ,&! observation type
                    GPSRO ,&! f90 module number
!                   -1    ,&! WMO BUFR types
                 (/WMOA_VERT_SOUND,-1/),    &! WMO BUFR types
                 (/1694,1695,(-1,i=3,ndb)/))&! DWD Datenbankkennzahlen
                 ,         &!
    ! 12:
    t_rept_char (  'GPSGB',&! mnemonic
                   'Ground Based GPS observations',& ! verbal description
                 OT_GPSGB ,&! observation type
                    GPSGB ,&! f90 module number
                    (/WMOA_SURF_LAND, -1/), &! WMO BUFR types
                    (/94,95,(-1,i=3,ndb)/) )&! DWD Datenbankkennzahlen
                 ,         &!
    ! 13:
    t_rept_char (  'RADAR',&! mnemonic
                   'Radar',&! verbal description
                 OT_RADAR ,&! observation type
                    RADAR ,&! f90 module number
                    -1    ,&! WMO BUFR types
                    -1    )&! DWD Datenbankkennzahlen
                 ,         &!
    ! 14:
    t_rept_char (  'POWER',&! mnemonic
                   'Power',&! verbal description
                 OT_POWER ,&! observation type
                    00000 ,&! f90 module number !++ not yet defined +++!
                    -1    ,&! WMO BUFR types
                    -1    )&! DWD Datenbankkennzahlen
                 ,         &!
    ! 15:
    t_rept_char (  'SOIL' ,&! mnemonic
                   'Soil' ,&! verbal description
                 OT_SOIL  ,&! observation type
                    SOIL  ,&! f90 module number
                     -1   ,&! WMO BUFR types
                  [ 1699,          &! ASCAT (level 2)
                    (-1,i=2,ndb)]) &! DWD Datenbankkennzahlen
                 ,         &!
    ! 16:
    t_rept_char (  'OBJECT',&! mnemonic
                   'Object',&! verbal description
                 OT_OBJECT ,&! observation type
                    000000 ,&! f90 module number
                    -1     ,&! WMO BUFR types
                    -1     )&! DWD Datenbankkennzahlen
                 ,          &!
    ! 17:
    t_rept_char (  'LIGHTN'   ,&! mnemonic
                   'Lightning',&! verbal description
                 OT_LIGHTN    ,&! observation type
                    000000    ,&! f90 module number
                    -1        ,&! WMO BUFR types
                    -1        )&! DWD Datenbankkennzahlen
                 ,          &!
    ! 18:
    t_rept_char (  'WLIDAR'    ,&! mnemonic
                   'Wind Lidar',&! verbal description
                 OT_WLIDAR     ,&! observation type
                    WLIDAR     ,&! f90 module number
                 (/WMOA_VERT_SOUND,-1/),&! WMO BUFR types
                 (/1815,(-1,i=2,ndb)/)) &! DWD Datenbankkennzahlen
                 ,          &!
    ! 19:
    t_rept_char (  'MWR',   &! mnemonic
                   'Microwave Radiometer',&! verbal description
                 OT_MWR    ,&! observation type
                    MWR    ,&! f90 module number
                    -1     ,&! WMO BUFR types
                    -1     )&! DWD Datenbankkennzahlen
    !-------------------------------------
                                /)

  ! This must be after the actual declaration (Intel compiler bug)
  protected          :: rept_char

  type (t_rept_use) ,target  ,save :: rept_use (n_ot)

  type (t_rept_stat)         ,save :: rept_stat0 ! default status variable
  type (t_rept_stat) ,target ,save :: rept_stat (0:ndb,  n_ot)
  type (t_rept_stat) ,target ,save :: obsv_stat (0:n_vn, n_ot)

  integer, save ::  dism_stat (n_chk, 0:n_ot) = 0
  integer, save ::  pass_stat (n_chk, 0:n_ot) = 0
  integer, save ::  reje_stat (n_chk, 0:n_ot) = 0
  integer, save :: oonly_stat (n_chk, 0:n_ot) = 0

  logical, save       :: seen_abort = .false.
  character(len=32)   :: abort_msg  = ""
  !===================================================
  ! module variables (write report file 'monREPRT.nc')
  !===================================================
  integer  ,parameter :: npr = 100000   ! Maximum no. of pending reports
  integer  ,save      :: ipr = 0        ! Current no. of pending reports
  integer  ,save      :: iu  = 0        ! unit number for write_report

  !------------------------
  ! Report title and format
  !------------------------
  character(len=*) ,parameter :: title = &
       '# type    station    date    time     lat.     lon.  &
       &code  dbkz fil.    rec.        id check       status      flags'
  character(len=*) ,parameter :: fm = &
       '(3(1x,a),2f9.3,2i6,i4,1x,i8,1x,i9,1x,2a,b32.32,1x,a,1x,i1,1x,a)'

  character(len=200),    save :: pending_report (npr)
  integer                     :: lenr(npr) = 0    ! Lengths of pend. reports

  !-----------
  ! interfaces
  !-----------
  interface write_report
    module procedure write_report
    module procedure write_head
  end interface

!==============================================================================
contains
!==============================================================================
  subroutine print_rept_stat (comment, unit, reports, vars)
  !----------------------------------------
  ! print observation type statistics table
  !----------------------------------------
  character(len=*) ,intent(in) ,optional :: comment  ! print this comment
  integer          ,intent(in) ,optional :: unit     ! write to Fortran unit
  logical          ,intent(in) ,optional :: reports  ! write report statistics
  logical          ,intent(in) ,optional :: vars     ! " single observation "

    !----------------
    ! local variables
    !----------------
    integer            :: irv         ! report or varno
    integer            :: iun         ! unit number
    integer            :: io          ! report type index
    integer            :: id          ! subtype index
    integer            :: check       ! check
    integer            :: n_processed ! no. processed
    integer            :: idbk        ! Datenbankkennzahl
    character (len=9)  :: mnem, dummy ! mnemonic
    character (len=80) :: desc        ! description
    type(t_obsid)      :: obsid
    logical            :: lr, lv      ! copy of optional arguments

    character(*) ,parameter :: fh(2) = ['(2x,a8,7a8,2a8,1x,a)   ',&! rep./obs.
                                        '(2x,a8,7a10,a7,a6,1x,a)'] ! header format
    character(*) ,parameter :: fh2   =  '(2x,"dbkz code")'         ! 2nd line
    character(*) ,parameter :: fb(2) = ['(1x,a9,7i8,2i8,1x,a)   ',&! rep./obs.
                                        '(1x,a9,7i10,i7,i6,1x,a)'] ! body   format
    type (t_rept_stat),pointer  :: stat (:,:)

    if (.not.dace% lpio) return
    call init_fdbk_tables
    !----------------------------
    ! process optional parameters
    !----------------------------
    iun = 6      ;if (present (unit   )) iun = unit
    lr  = .true. ;if (present (reports)) lr  = reports
    lv  = .true. ;if (present (vars   )) lv  = vars
    do irv = 1, 2
      if (irv == 1) then
        if (.not. lr) cycle
        stat => rept_stat
      else
        if (.not. lv) cycle
        stat => obsv_stat
      endif

      !---------------------------
      ! update report type summary
      !---------------------------
      stat(0,:)% processed  = sum (stat(1:,:)% processed  , dim=1)
      stat(0,:)% merged     = sum (stat(1:,:)% merged     , dim=1)
      stat(0,:)% obs_only   = sum (stat(1:,:)% obs_only   , dim=1)
      stat(0,:)% dismissed  = sum (stat(1:,:)% dismissed  , dim=1)
      stat(0,:)% passive    = sum (stat(1:,:)% passive    , dim=1)
      stat(0,:)% rejected   = sum (stat(1:,:)% rejected   , dim=1)
      stat(0,:)% active     = sum (stat(1:,:)% active     , dim=1)
      stat(0,:)% accepted   = sum (stat(1:,:)% accepted   , dim=1)

      !-------------
      ! write header
      !-------------
!     write(iun,'(a)') repeat('-',79)
      write(iun,'()')
      if (irv == 1) then
        if (present(comment)) then
          write(iun,'(4x,a)') 'Observation Report Statistics '//trim(comment)
        else
          write(iun,'(4x,a)') 'Observation Report Statistics'
        endif
      else
        if (present(comment)) then
          write(iun,'(4x,a)') 'Observation Statistics '//trim(comment)
        else
          write(iun,'(4x,a)') 'Observation Statistics'
        endif
      endif
      write(iun,'()')
      if (all (stat(0,:)% processed == 0)) then
        write(iun,'(a)') '  No statistics available'
        write(iun,'( )')
        cycle
      endif
      write(iun,fh(irv)) &
                    'obstype','proces.','merged','obsonly','dism.','passive',&
                    'reject.','active' ,'acc.'  ,'lost' ,'description'
      write(iun,fh2)

      !----------------------------------------
      ! loop over observation types, write body
      !----------------------------------------
      do io = 1, n_ot
        write(iun,'("  '//repeat('-',128)//'")')
        write(iun,fb(irv)) rept_char(  io)% mnem ,&
                           stat(0,io)% processed ,&
                           stat(0,io)% merged    ,&
                           stat(0,io)% obs_only  ,&
                           stat(0,io)% dismissed ,&
                           stat(0,io)% passive   ,&
                           stat(0,io)% rejected  ,&
                           stat(0,io)% active    ,&
                           stat(0,io)% accepted  ,&
                      lost(stat(0,io))           ,&
                      trim(rept_char(  io)% desc)
        if (irv == 1) then
          !----------------------------
          ! loop over Datenbankkennzahl
          !----------------------------
          do id = 1, ndb
            n_processed = stat(id,io)% processed
            idbk        = rept_char(io)% dbk(id)
            if (n_processed == 0) cycle             ! skip empty lines
            if (idbk   >= 0) then              ! Datenbankkennzahl present
              obsid = obstype_dbkz (idbk)
              write (mnem,'(i5,1x,i3)') idbk, obsid% codetype
              call dbkz_mnem (idbk,dummy,desc)
            else                               ! no Datenbankkennzahl present
              if (id == 1) exit
              mnem   = '   other'
              desc   = 'no valid Datenbankkennzahl'
            endif
            write(iun,fb(irv))             mnem        ,&
                                           n_processed ,&
                              stat(id,io)% merged      ,&
                              stat(id,io)% obs_only    ,&
                              stat(id,io)% dismissed   ,&
                              stat(id,io)% passive     ,&
                              stat(id,io)% rejected    ,&
                              stat(id,io)% active      ,&
                              stat(id,io)% accepted    ,&
                         lost(stat(id,io))             ,&
                         trim(                  desc)
            if (idbk < 0) exit
          end do
        else
          !--------------------
          ! loop over variables
          !--------------------
          do id = 1, n_vn
            n_processed = stat(id,io)% processed
            idbk        = varno% e (id)% value
            mnem        = varno% e (id)% name
            desc        = varno% e (id)% description
            if (n_processed == 0) cycle             ! skip empty lines
            write(iun,fb(irv))                  mnem        ,&
                                                n_processed ,&
                              stat(id,io)% merged      ,&
                              stat(id,io)% obs_only    ,&
                              stat(id,io)% dismissed   ,&
                              stat(id,io)% passive     ,&
                              stat(id,io)% rejected    ,&
                              stat(id,io)% active      ,&
                              stat(id,io)% accepted    ,&
                         lost(stat(id,io))             ,&
                         trim(                  desc)
            if (idbk < 0) exit
          end do
        endif
      end do
      !--------------------------
      ! write summary and trailer
      !--------------------------
      write(iun,'("  '//repeat('-',128)//'")')
      write(iun,fb(irv))               'total   '  ,&
                   sum      (stat(0,:)% processed) ,&
                   sum      (stat(0,:)% merged)    ,&
                   sum      (stat(0,:)% obs_only)  ,&
                   sum      (stat(0,:)% dismissed) ,&
                   sum      (stat(0,:)% passive)   ,&
                   sum      (stat(0,:)% rejected)  ,&
                   sum      (stat(0,:)% active)    ,&
                   sum      (stat(0,:)% accepted)  ,&
                   sum (lost(stat(0,:)))
      if (irv == 1) then
        !--------------------------------
        ! statistics on DISMISSED reports
        !--------------------------------
        do io = 0, n_ot
          if (sum(dism_stat (:,io)) == 0) cycle
          mnem = 'UNKNOWN'; if (io>0) mnem = rept_char(io)% mnem
          write(iun,'()')
          write(iun,'(a,a,a)') '  ',trim(mnem),' reports DISMISSED due to:'
          write(iun,'()')
          do check = 1, n_chk
            if (dism_stat (check,io) /= 0) write (iun,'(4x,a12,i8)') &
              chk(check)% mnem, dism_stat (check,io)
          end do
          write (iun,'(4x,a12,i8)') 'total',sum(dism_stat (:,io))
        end do
        !-------------------------------
        ! statistics on OBS_ONLY reports
        !-------------------------------
        do io = 0, n_ot
          if (sum(oonly_stat (:,io)) == 0) cycle
          mnem = 'UNKNOWN'; if (io>0) mnem = rept_char(io)% mnem
          write(iun,'()')
          write(iun,'(a,a,a)') '  ',trim(mnem),' reports OBS_ONLY due to:'
          write(iun,'()')
          do check = 1, n_chk
            if (oonly_stat (check,io) /= 0) write (iun,'(4x,a12,i8)') &
              chk(check)% mnem, oonly_stat (check,io)
          end do
          write (iun,'(4x,a12,i8)') 'total',sum(oonly_stat (:,io))
        end do
        !------------------------------
        ! statistics on PASSIVE reports
        !------------------------------
        do io = 0, n_ot
          if (sum(pass_stat (:,io)) == 0) cycle
          mnem = 'UNKNOWN'; if (io>0) mnem = rept_char(io)% mnem
          write(iun,'()')
          write(iun,'(a,a,a)') '  ',trim(mnem),' reports PASSIVE due to:'
          write(iun,'()')
          do check = 1, n_chk
            if (pass_stat (check,io) /= 0) write (iun,'(4x,a12,i8)') &
              chk(check)% mnem, pass_stat (check,io)
          end do
          write (iun,'(4x,a12,i8)') 'total',sum(pass_stat (:,io))
        end do
        !-------------------------------
        ! statistics on REJECTED reports
        !-------------------------------
        do io = 0, n_ot
          if (sum(reje_stat (:,io)) == 0) cycle
          mnem = 'UNKNOWN'; if (io>0) mnem = rept_char(io)% mnem
          write(iun,'()')
          write(iun,'(a,a,a)') '  ',trim(mnem),' reports REJECTED due to:'
          write(iun,'()')
          do check = 1, n_chk
            if (reje_stat (check,io) /= 0) write (iun,'(4x,a12,i8)') &
              chk(check)% mnem, reje_stat (check,io)
          end do
          write (iun,'(4x,a12,i8)') 'total',sum(reje_stat (:,io))
        end do
      endif
      !--------------
      ! write trailer
      !--------------
      write(iun,'()')
      write(iun,'(a)') repeat('-',79)
    end do
  end subroutine print_rept_stat
!------------------------------------------------------------------------------
  subroutine derive_rept_stat (obs, reports, vars, step)
  !--------------------------------------------------------------------
  ! re-derive observation type statistics from observation type content
  !--------------------------------------------------------------------
  type (t_obs) ,intent(in)           :: obs (:) ! observation data
  logical      ,intent(in) ,optional :: reports ! statistics for reports
  logical      ,intent(in) ,optional :: vars    ! stats.for single observations
  integer      ,intent(in) ,optional :: step    ! for specific LETKF step only

    integer                     :: ib, is, i, j
    integer                     :: state, check, obstype, idbk
    logical                     :: lr, lv       ! copy of optional arguments
    integer                     :: ls           ! copy of optional argument
    type (t_rept_stat) ,pointer :: stat
    type (t_rept_stat) ,pointer :: stats

    !----------------------------
    ! process optional parameters
    !----------------------------
    lr = .true.; if (present (reports)) lr = reports
    lv = .true.; if (present (vars   )) lv = vars
    ls = 0;      if (present (step   )) ls = step

    !----------------------------------------
    ! set observation type statistics to zero
    !----------------------------------------
    if (lv) then
      obsv_stat = rept_stat0
    endif
    if (lr) then
       rept_stat = rept_stat0
       dism_stat = 0
       pass_stat = 0
       reje_stat = 0
      oonly_stat = 0
    endif

    !---------------------
    ! re-derive statistics
    !---------------------
    call init_fdbk_tables
    do ib = 1, size (obs)
      if (obs(ib)% pe /= dace% pe) cycle
      do is = 1, obs(ib)% n_spot
        obstype =  obs(ib)% spot(is)% hd% obstype
        !--------
        ! reports
        !--------
        if (lr) then
          idbk  =  obs(ib)% spot(is)% hd% idbk
          state =  obs(ib)% spot(is)% use% state
          check =  obs(ib)% spot(is)% use% check
          stat  => rept_stat (   0, obstype)
          stats => rept_stat (idbk, obstype)
          select case (state)
          case (STAT_ABORT)
          case (STAT_FORGET)
          case (STAT_DISMISS)
            stat % dismissed = stat % dismissed + 1
            stats% dismissed = stats% dismissed + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
            dism_stat   (check, obstype) = &
              dism_stat (check, obstype) + 1
          case (STAT_MERGED)
            stat % merged    = stat % merged    + 1
            stats% merged    = stats% merged    + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
          case (STAT_PASSIVE, STAT_PAS_REJ, STAT_NOTACTIVE)
            stat % passive   = stat % passive   + 1
            stats% passive   = stats% passive   + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
            pass_stat   (check, obstype) = &
              pass_stat (check, obstype) + 1
          case (STAT_REJECTED)
            stat % rejected  = stat % rejected  + 1
            stats% rejected  = stats% rejected  + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
            reje_stat   (check, obstype) = &
              reje_stat (check, obstype) + 1
          case (STAT_OBS_ONLY)
            stat % obs_only  = stat % obs_only  + 1
            stats% obs_only  = stats% obs_only  + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
            oonly_stat   (check, obstype) = &
              oonly_stat (check, obstype) + 1
          case (STAT_ACTIVE_0I:STAT_ACTIVE_1)
            stat % active    = stat % active    + 1
            stats% active    = stats% active    + 1
            stat % processed = stat % processed + 1
            stats% processed = stats% processed + 1
          case default
            write (0,*) 'derive_rept_stat:  invalid report state',state
            call finish('derive_rept_stat','invalid report state')
          end select
        endif
        !----------
        ! variables
        !----------
        if (lv) then
          do i = obs(ib)% spot(is)% o%i + 1 ,                  &
                 obs(ib)% spot(is)% o%i + obs(ib)% spot(is)% o%n
            if (present(step) .and. ls /= obs(ib)% body(i)%set% ekf_pass) cycle
            do j = 1, n_vn
              if (varno% e (j)% value == obs(ib)% varno(i)) exit
            end do
            state =  obs(ib)% body(i)% use% state
            stat  => obsv_stat (0, obstype)
            stats => obsv_stat (j, obstype)
            select case (state)
            case (STAT_ABORT)
            case (STAT_FORGET)
            case (STAT_OBS_ONLY)
              stat % obs_only  = stat % obs_only  + 1
              stats% obs_only  = stats% obs_only  + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_DISMISS)
              stat % dismissed = stat % dismissed + 1
              stats% dismissed = stats% dismissed + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_MERGED)
              stat % merged    = stat % merged    + 1
              stats% merged    = stats% merged    + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_PASSIVE, STAT_PAS_REJ, STAT_NOTACTIVE)
              stat % passive   = stat % passive   + 1
              stats% passive   = stats% passive   + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_REJECTED)
              stat % rejected  = stat % rejected  + 1
              stats% rejected  = stats% rejected  + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_ACTIVE_0I:STAT_ACTIVE_1,STAT_ACCEPTED+1:)
              stat % active    = stat % active    + 1
              stats% active    = stats% active    + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case (STAT_ACCEPTED)
              stat % accepted  = stat % accepted  + 1
              stats% accepted  = stats% accepted  + 1
              stat % processed = stat % processed + 1
              stats% processed = stats% processed + 1
            case default
              write (0,*) 'derive_rept_stat:  invalid variable state',state
              call finish('derive_rept_stat','invalid variable state')
            end select
          end do
        endif
      end do
    end do

  end subroutine derive_rept_stat
!------------------------------------------------------------------------------
  elemental function lost (entry)
  integer                        :: lost
  type (t_rept_stat) ,intent(in) :: entry
    lost = entry% processed &
         - entry% merged    &
         - entry% obs_only  &
         - entry% dismissed &
         - entry% passive   &
         - entry% rejected  &
         - entry% active    &
         - entry% accepted
  end function lost
!------------------------------------------------------------------------------
  subroutine read_nml_report
  !----------------------------
  ! read namelist 'REPORT'
  !----------------------------
    !-------------------
    ! namelist variables
    !-------------------
    integer           :: init
    character(len=12) :: check
    character(len=12) :: type
    character(len=12) :: use
    integer           :: max_proc
    integer           :: max_act
    integer           :: ni
    real(wp)          :: height_t
    real(wp)          :: height_b
    real(wp)          :: height_top
    real(wp)          :: height_bot
    real(wp)          :: lat_nb
    real(wp)          :: lat_sb
    real(wp)          :: lon_eb
    real(wp)          :: lon_wb
    real(wp)          :: excl_bnd    ! distance excluded at boundaries [degree]
    integer           :: time_b      ! begin of valid time range       [+-hhmm]
    integer           :: time_e      ! end   of valid time range       [+-hhmm]
    integer           :: time_cutoff ! cut-off time                      [hhmm]
    integer           :: time_s1     ! apply check to section 1 (synop) time
    real(wp)          :: min_dist    ! min. distance between observations  [km]
    real(wp)          :: minvdist    ! min. vertical distance             [hPa]
    integer           :: hours (8)
    real(wp)          :: fr_land     ! land fraction reqired
    real(wp)          :: fr_sea      ! sea  fraction reqired
    real(wp)          :: fr_noice    !      fraction reqired
    real(wp)          :: min_tsurf   ! minimum surface temperature allowed
    integer           :: deriv_p     ! derive p from h if not present
    integer           :: read1pe     ! read NetCDF files from 1 PE:   [0=F,1=T]
    integer           :: inv_op (6)  ! observation operator validity check

    integer           :: moni_cntr(nc)! generating centers (monitoring)
    integer           :: use_cntr (nc)! generating centers (assimilation)

    namelist /REPORT/ check, type, use, max_act, max_proc, ni, min_dist,  &
                      height_t, height_b, height_top, height_bot,         &
                      lat_nb, lat_sb, lon_eb, lon_wb,                     &
                      time_b, time_e, hours, minvdist, fr_land, fr_sea,   &
                      fr_noice, min_tsurf, init, deriv_p, moni_cntr,      &
                      use_cntr, read1pe, excl_bnd, time_cutoff, inv_op,   &
                      time_s1
    !---------------------------------------------------------
    ! Auxiliary interface for specification of argument intent
    !---------------------------------------------------------
    INTERFACE
      SUBROUTINE p_bcast_derivedtype (buffer, count, source, comm)
      import :: t_rept_use
      TYPE(t_rept_use) ,INTENT(inout) :: buffer(*) ! variable to bcast
      INTEGER          ,INTENT(in)    :: count     ! len(byte) of variable
      INTEGER          ,INTENT(in)    :: source    ! source processor index
      INTEGER          ,INTENT(in)    :: comm      ! communicator
      END SUBROUTINE p_bcast_derivedtype
    END interface
    !----------------
    ! local variables
    !----------------
    integer          :: i, ierr
    logical          :: first = .true.
    integer          :: it1, it2        ! report type indices
    integer          :: iuse, ichk
    integer          :: l               ! size of data type in byte
#if defined(__ibm__)
    integer          :: ios             ! iostat return parameter
#endif

    !----------------------------
    ! call this routine only once
    !----------------------------
    if (.not. first) then
      if (dace% lpio) then
        write(6,'(a)') repeat('-',79)
        write(6,'(a)')
        write(6,'(a)') &
               ' namelist /REPORT/ already read: skipping read_nml_report'
        write(6,'(a)')
      endif
      return
    endif

    !-------------------------------------
    ! set default values in table rept_use
    !-------------------------------------
    call set_rept_use
    call print_rept_use ('Default')

    !-------------------------------
    ! read namelist group repeatedly
    !-------------------------------
    do
      !----------------------------------------
      ! preset namelist variables with defaults
      !----------------------------------------
      init        = OU_VAR
      check       = ''
      type        = ''
      use         = ''
      max_act     = huge (0)
      max_proc    = huge (0)
      ni          = huge (0)
      height_t    = huge (0._wp)
      height_b    = huge (0._wp)
      height_top  = huge (0._wp)
      height_bot  = huge (0._wp)
      lat_nb      = huge (0._wp)
      lat_sb      = huge (0._wp)
      lon_eb      = huge (0._wp)
      lon_wb      = huge (0._wp)
      excl_bnd    = huge (0._wp)
      min_dist    = huge (0._wp)
      minvdist    = huge (0._wp)
      time_b      = huge (0)
      time_e      = huge (0)
      time_cutoff = huge (0)
      time_s1     = -1
      hours       = -1
      fr_land     = huge (0._wp)
      fr_sea      = huge (0._wp)
      fr_noice    = huge (0._wp)
      min_tsurf   = huge (0._wp)
      deriv_p     = huge (0)
      read1pe     = -1
      inv_op      = -1
      moni_cntr   = -9
      use_cntr    = -9
      !--------------
      ! read namelist
      !--------------
      if (dace% lpio) then
        call position_nml ('REPORT' ,lrewind=first ,status=ierr)
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
          read (nnml ,nml=REPORT, iostat=ios)
          if (ios/=0) call finish ('read_nml_report',          &
                                   'ERROR in namelist /REPORT/')
#else
          read (nnml ,nml=REPORT)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (dace% lpio) then
        if (first) then
          write(6,'(a)') repeat('-',79)
          write(6,'(a)')
          if (ierr == POSITIONED) then
             write(6,'(a)') ' read namelist /REPORT/:'
          else
             write(6,'(a,i0)') ' namelist /REPORT/ not found: ierr = ', ierr
          end if
          write(6,'(a)')
        endif
      endif
      first = .false.
      if (ierr /= POSITIONED) exit
      !----------------------
      ! determine report type
      !----------------------
      if (type == '') then
        it1 = 1
        it2 = n_ot
      else
        it1 = 0
        do i = 1, n_ot
          if (rept_char(i)% mnem == type) then
            it1 = i
            it2 = i
            exit
          endif
        end do
        if (it1 == 0) call finish ('read_nml_rept_use',&
          'invalid value for variable TYPE in namelist /REPORT/: '//type)
      endif
      !-------------------
      ! determine use flag
      !-------------------
      iuse = stat_key (use)
      if (use==' ') iuse = 0
      if (iuse < 0) call finish ('read_nml_rept_use',&
        'invalid value for variable USE in namelist /REPORT/: '//use)
      !-------------
      ! set use flag
      !-------------
      if (use /= '') then
        ichk = chk_key (check)
        if (check==' ') ichk = CHK_NONE
        if (ichk < 0) call finish ('read_nml_rept_use',&
          'invalid value for variable CHECK in namelist /REPORT/: '//check)
        rept_use (it1:it2)% use (ichk) = iuse
      endif
      if (dace% lpio) then
        if (type/=' ') write(6,'("  type = ",a10," ,use = ",a10)') type, use
      endif
      !--------------
      ! set init flag
      !--------------
      if (it1==it2) rept_use (it1)% init = init
      rept_use (OT_RADAR)% init = min (rept_use (OT_RADAR)% init, OU_FWD)
      !----------------
      ! set other flags
      !----------------
      if (read1pe==0.or.read1pe==1) rept_use (it1:it2)% read1pe   = read1pe==1
      if (max_act  /= huge (0)    ) rept_use (it1:it2)% max_act   = max_act
      if (max_proc /= huge (0)    ) rept_use (it1:it2)% max_proc  = max_proc
      if (ni       /= huge (0)    ) rept_use (it1:it2)% ni        = ni
      if (deriv_p  /= huge (0)    ) rept_use (it1:it2)% deriv_p   = deriv_p
      if (height_t /= huge (0._wp)) rept_use (it1:it2)% height_t  = height_t
      if (height_b /= huge (0._wp)) rept_use (it1:it2)% height_b  = height_b
      if (height_top/=huge (0._wp)) rept_use (it1:it2)% height_top= height_top
      if (height_bot/=huge (0._wp)) rept_use (it1:it2)% height_bot= height_bot
      if (lat_nb   /= huge (0._wp)) rept_use (it1:it2)% lat_nb    = lat_nb
      if (lat_sb   /= huge (0._wp)) rept_use (it1:it2)% lat_sb    = lat_sb
      if (lon_eb   /= huge (0._wp)) rept_use (it1:it2)% lon_eb    = lon_eb
      if (lon_wb   /= huge (0._wp)) rept_use (it1:it2)% lon_wb    = lon_wb
      if (excl_bnd /= huge (0._wp)) rept_use (it1:it2)% excl_bnd  = excl_bnd
      if (min_dist /= huge (0._wp)) rept_use (it1:it2)% min_dist = min_dist &
                                                                 * 1000._wp
      if (minvdist /= huge (0._wp)) rept_use (it1:it2)% minvdist = minvdist &
                                                                 *  100._wp
      if (fr_land  /= huge (0._wp)) rept_use (it1:it2)% fr_land  = fr_land
      if (fr_sea   /= huge (0._wp)) rept_use (it1:it2)% fr_sea   = fr_sea
      if (fr_noice /= huge (0._wp)) rept_use (it1:it2)% fr_noice = fr_noice
      if (min_tsurf/= huge (0._wp)) rept_use (it1:it2)% min_tsurf= min_tsurf
      if (time_s1  ==       1     ) rept_use (it1:it2)% time_s1  = .true.
      if (time_s1  ==       0     ) rept_use (it1:it2)% time_s1  = .false.
      if (time_b   /= huge (0)    ) then
        call init_time (rept_use (it1:it2)% time_b, hh=      time_b/100,&
                                                    mi= mod (time_b,100))
        rept_use (it1:it2)% time_b = rept_use (it1:it2)% time_b + ana_time
      endif
      if (time_e   /= huge (0)    ) then
        call init_time (rept_use (it1:it2)% time_e, hh=      time_e/100,&
                                                    mi= mod (time_e,100))
        rept_use (it1:it2)% time_e = rept_use (it1:it2)% time_e + ana_time
      endif
      if (time_cutoff /= huge (0)    ) then
        call init_time (rept_use (it1:it2)% time_cutoff,     &
                                   hh=      time_cutoff/100, &
                                   mi= mod (time_cutoff,100))
        rept_use (it1:it2)% time_cutoff = rept_use (it1:it2)% time_cutoff &
                                        + ana_time
!     else
!       rept_use (it1:it2)% time_cutoff = invalid_time
      endif
      if (any (hours >= 0        )) rept_use (it1:it2)% lhour    = &
                                    any (hours == ihh (ana_time))
      if (any (moni_cntr /= -9)) then
        do i = it1, it2
          rept_use (i)% moni_cntr = moni_cntr
        end do
      end if
      if (any ( use_cntr /= -9)) then
        do i = it1, it2
          rept_use (i)%  use_cntr =  use_cntr
        end do
      end if
      if (any (inv_op /= -1)) then
      do i = it1, it2
        rept_use (i)% inv_op = inv_op
      end do
      end if
      !--------------------------------------------------------------------
      ! Take over default settings for height_top, height_bot when not set,
      ! while height_t, height_b pressure level bounds are set.
      !--------------------------------------------------------------------
      do i = it1, it2
         if (rept_use(i)% height_top == 99999._wp .and. &
             rept_use(i)% height_t   >      0._wp       ) then
            rept_use (i)% height_top = h_p_usstd (rept_use(i)% height_t * 100)
         end if
         if (rept_use(i)% height_bot ==  -999._wp .and. &
             rept_use(i)% height_b   <   1100._wp       ) then
            rept_use (i)% height_bot = h_p_usstd (rept_use(i)% height_b * 100)
         end if
      end do
      !-------------------
      ! consistency checks
      !-------------------
      rept_use(:)% use(CHK_DOMAIN) = min (STAT_DISMISS, rept_use(:)% use(CHK_DOMAIN))
      !----------
      ! broadcast
      !----------
      l = size (transfer (rept_use ,(/' '/)))
      call p_bcast_derivedtype (rept_use, l, dace% pio, dace% comm)
    end do
    !---------
    ! printout
    !---------
    call print_rept_use ('Namelist')

  end subroutine read_nml_report
!------------------------------------------------------------------------------
  subroutine set_rept_use

    integer          :: i
    integer          :: init
    type(t_time)     :: halfhour

    !---------------------------------------
    ! set default values in table 'rept_use'
    !---------------------------------------
    ! effect of failure of certain checks
    !------------------------------------
    do i=1,size(rept_use)
      rept_use(i) % use = chk% default_state
    end do

    !-------------------------------
    ! activate selected report types
    !-------------------------------
    select case (method)
    case default
      init = OU_NONE                           ! operators not used
    case ('LETKF')
      init = OU_LETKF                          ! only store obs/fg values
    case ('MEC', 'GMESTAT', 'VERI_ENS', 'VERI_ENSFC')
      init = OU_FWD                            ! use forward operator
    case ('PSAS', 'ENVAR', 'PSAS+LETKF', 'ENVAR+LETKF')
      init = OU_VAR                            ! use forward operator and Jakobian
    end select

    rept_use (:)       % use (CHK_NONE) = STAT_FORGET  !
    rept_use (OT_TEMP) % use (CHK_NONE) = STAT_ACTIVE  ! activate TEMPs
    rept_use (OT_PILOT)% use (CHK_NONE) = STAT_ACTIVE  ! activate PILOTs
    rept_use (OT_SYNOP)% use (CHK_NONE) = STAT_ACTIVE  ! activate SYNOPs
    rept_use (OT_AIREP)% use (CHK_NONE) = STAT_ACTIVE  ! activate AIREPs
    rept_use (OT_SCATT)% use (CHK_NONE) = STAT_ACTIVE  ! activate SCATTerometer
    rept_use (OT_DRIBU)% use (CHK_NONE) = STAT_ACTIVE  ! activate DRIBUs
    rept_use (OT_SATOB)% use (CHK_NONE) = STAT_ACTIVE  ! activate AMVs
    rept_use (OT_GPSGB)% use (CHK_NONE) = STAT_PASSIVE ! passive  ZTD/STD
    rept_use (:)       % init           = init         ! module initialisation ?
!    rept_use (OT_RAD)  % init           = OU_NONE      ! RADIANCES not used by default
    rept_use (OT_RADAR)% init           = OU_NONE      ! RADAR     not used by default

    select case (method)
    !--------------------------------------------
    ! additional observation types in COSMO LETKF
    !--------------------------------------------
    case ('LETKF')
      rept_use (OT_GPSGB)% use (CHK_NONE)     = STAT_ACTIVE
      rept_use           % use (CHK_OPERATOR) = STAT_PASSIVE
    !----------------------------------------
    ! additional observation types in GMESTAT
    !----------------------------------------
    case ('GMESTAT','VERI_ENS')
      rept_use (OT_GPSRO)% use (CHK_NONE) = STAT_ACTIVE
    !------------------------------------------
    ! additional observation types in COSMO MEC
    !------------------------------------------
    case ('MEC')
      rept_use  (OT_GPSGB)% use (CHK_NONE) = STAT_ACTIVE
      rept_use  (OT_SOIL )% use (CHK_NONE) = STAT_ACTIVE
      rept_char (OT_SYNOP)% mod = COSMO
      rept_char (OT_DRIBU)% mod = COSMO
      rept_char (OT_AIREP)% mod = COSMO
      rept_char (OT_PILOT)% mod = COSMO
      rept_char (OT_TEMP )% mod = COSMO
      rept_char (OT_GPSGB)% mod = GPSGB         ! do not use old COSMO code

      if (n_dace_op > 0) then
         do i=1,n_dace_op
           select case(trim(dace_op(i)))
           case('SYNOP')
              rept_char (OT_SYNOP)% mod = SYNOP
           case('DRIBU')
              !rept_char (OT_DRIBU)% mod = SYNOP
              print*, 'in subroutine set_rept_use: DRIBU not in modtype'
           case('AIREP')
              rept_char (OT_AIREP)% mod = AIREP
           case('PILOT')
              !rept_char (OT_PILOT)% mod = TEMP
              print*, 'in subroutine set_rept_use: PILOT not in modtype'
           case('TEMP')
              rept_char (OT_TEMP)% mod  = TEMP
           case default
           end select
        end do
      end if

    !----------------------------------------
    ! less restrictive settings for COSMO MEC
    !----------------------------------------
      rept_use (:)% use (CHK_HEIGHT ) = STAT_PASSIVE
      rept_use (:)% use (CHK_SURF   ) = STAT_REJECTED
      rept_use (:)% use (CHK_NOTUSED) = STAT_PASSIVE
      rept_use (:)% use (CHK_THIN   ) = STAT_PASSIVE
      rept_use (:)% use (CHK_AREA   ) = STAT_PASSIVE
    end select

    !-----------------------------
    ! set default bounds on checks
    !-----------------------------
    rept_use (:)% max_act     = huge (0)
    rept_use (:)% max_proc    = huge (0)
    rept_use (:)% height_t    =    0._wp
    rept_use (:)% min_dist    =    0._wp
    rept_use (:)% minvdist    =    0._wp
    rept_use (:)% height_b    = 1100._wp
    rept_use (:)% height_top  =99999._wp
    rept_use (:)% height_bot  = -999._wp
    rept_use (:)% lat_nb      =   90._wp
    rept_use (:)% lat_sb      =  -90._wp
    rept_use (:)% lon_eb      =  180._wp
    rept_use (:)% lon_wb      = -180._wp
    rept_use (:)% excl_bnd    =    0._wp
    rept_use (:)% time_cutoff = invalid_time
    rept_use (:)% lhour       = .true.
    rept_use (:)% fr_land     = 0._wp
    rept_use (:)% fr_sea      = 0._wp
    rept_use (:)% fr_noice    = 0._wp
    rept_use (:)% min_tsurf   = -999._wp
    rept_use (:)% deriv_p     = 0
    rept_use (:)% read1pe     = .true.
    do i=1,size(rept_use)
      rept_use (i)% inv_op    = -1
    end do
    !---------------------------------------
    ! report type specific defaults: surface
    !---------------------------------------
    rept_use (OT_RAD)  % fr_sea   = 1._wp
    rept_use (OT_PILOT)% deriv_p  = 1            ! derive pressure
    rept_use (OT_AIREP)% deriv_p  = 1            !   from height
    !------------------------------------------
    ! report type specific defaults: time frame
    !------------------------------------------
    call init_time (halfhour,  mi=30)
    rept_use (:)         % time_b   = fc_ref_time + obs_timeshift
    rept_use (:)         % time_e   = ana_time    + obs_timeshift
    select case (method)
    case ('LETKF','MEC','GMESTAT')
      rept_use (:)       % time_s1  = .false.
    case default
      rept_use (OT_SYNOP)% time_b   = ana_time - halfhour
      rept_use (OT_SYNOP)% time_e   = ana_time + halfhour
      rept_use (OT_SATOB)% time_b   = ana_time - halfhour
      rept_use (OT_SATOB)% time_e   = ana_time + halfhour
      rept_use (OT_DRIBU)% time_b   = ana_time - halfhour
      rept_use (OT_DRIBU)% time_e   = ana_time + halfhour
      rept_use (:)       % time_s1  = .true.
    end select
    !----------------------------------------
    ! : minimum distance between observations
    !----------------------------------------
    rept_use (OT_SATOB)% min_dist = 100000._wp
    rept_use (OT_RAD)  % min_dist = 100000._wp
    rept_use (OT_AIREP)% min_dist =  50000._wp
    rept_use (OT_SATOB)% minvdist =   4000._wp
    rept_use (OT_AIREP)% minvdist =   4000._wp
    !--------------------------------
    ! generated centers accepted
    !
    ! -1: not specified in input file
    ! -2: any
    ! -3: none
    !--------------------------------
    do i=1,size(rept_use)
      rept_use(i)% moni_cntr(:) = -9
      rept_use(i)%  use_cntr(:) = -9
      rept_use(i)% moni_cntr(1) = ANY_CNTR
      rept_use(i)%  use_cntr(1) = ANY_CNTR
    end do

  end subroutine set_rept_use
!------------------------------------------------------------------------------
  subroutine print_rept_use (text)
  character(len=*) ,intent(in) :: text

    integer          :: i, j
    character(len=6) :: chhmmss_cutoff (n_ot)
    character(len=9) :: c
    character(len=8) :: c_cntr (n_ot)

    if (dace% lpio) then
      where (invalid_time /= rept_use % time_cutoff)
        chhmmss_cutoff = chhmmss(rept_use %time_cutoff)
      elsewhere
        chhmmss_cutoff = '******'
      endwhere
      write (6,'()')
      write (6,'(1x,a,a)') text," settings in table 'rept_use' :"
      write (6,'(" type        : ",99(a8  ,1x))')         rept_char  % mnem
      write (6,'(" max_proc    : ",99(i6  ,3x))')         rept_use   % max_proc
      write (6,'(" max_act     : ",99(i6  ,3x))')         rept_use   % max_act
      write (6,'(" ni          : ",99(i6  ,3x))')         rept_use   % ni
      write (6,'(" deriv_p     : ",99(i6  ,3x))')         rept_use   % deriv_p
      write (6,'(" height_t    : ",99(f8.1,1x))')         rept_use   % height_t
      write (6,'(" height_b    : ",99(f8.1,1x))')         rept_use   % height_b
      write (6,'(" height_top  : ",99(f8.1,1x))')         rept_use   % height_top
      write (6,'(" height_bot  : ",99(f8.1,1x))')         rept_use   % height_bot
      write (6,'(" lat_nb      : ",99(f8.1,1x))')         rept_use   % lat_nb
      write (6,'(" lat_sb      : ",99(f8.1,1x))')         rept_use   % lat_sb
      write (6,'(" lon_eb      : ",99(f8.1,1x))')         rept_use   % lon_eb
      write (6,'(" lon_wb      : ",99(f8.1,1x))')         rept_use   % lon_wb
      write (6,'(" excl_bnd    : ",99(f8.1,1x))')         rept_use   % excl_bnd
      write (6,'(" time_b      : ",99(a8  ,1x))') chhmmss(rept_use   % time_b)
      write (6,'(" time_e      : ",99(a8  ,1x))') chhmmss(rept_use   % time_e)
      write (6,'(" time_cutoff : ",99(a8  ,1x))') chhmmss_cutoff
      write (6,'(" time_s1     : ",99(l6,3x))')           rept_use   % time_s1
      write (6,'(" fr_land     : ",99(f8.1,1x))')         rept_use   % fr_land
      write (6,'(" fr_sea      : ",99(f8.1,1x))')         rept_use   % fr_sea
      write (6,'(" fr_noice    : ",99(f8.1,1x))')         rept_use   % fr_noice
      write (6,'(" min_tsurf   : ",99(f8.1,1x))')         rept_use   % min_tsurf
      write (6,'(" lhour       : ",99(l6,3x))')           rept_use   % lhour
      write (6,'(" read1pe     : ",99(l6,3x))')           rept_use   % read1pe
      write (6,'(" inv_op      : ",99(6i1,3x))')      (/ (rept_use(i)% inv_op,  &
                                                         i=1,size(rept_use))  /)
      write (6,'(" init        : ",99(i6  ,3x))')         rept_use   % init
      write (6,'(a)') ' checks      :'
      do i=1,n_chk
        write (6,'(1x,a10,  "  : ",99(a8  ,1x))') chk(i)% mnem,               &
                                                  stat_mnem(rept_use(:)%use(i))
      end do
      c = 'used cntr'
      do i=1,nc
       if (all (rept_use(:)% use_cntr(i) < NONE_CNTR)) cycle
        do j=1,n_ot
          select case (rept_use(j)% use_cntr(i))
          case (:-4)
            c_cntr (j) = ''
          case (NONE_CNTR)
            c_cntr (j) = '   none '
          case (ANY_CNTR)
            c_cntr (j) = '    any '
          case (UKN_CNTR)
            c_cntr (j) = 'unknown '
          case default
            write(c_cntr(j),'(i7,1x)') rept_use(j)% use_cntr(i)
          end select
        end do
        write (6,'(1x,a10,  "  : ",99(a8  ,1x))') c,c_cntr
        c = ''
      end do
      c = 'moni cntr'
      do i=1,nc
        if (all (rept_use(:)% moni_cntr(i) < NONE_CNTR)) cycle
        do j=1,n_ot
          select case (rept_use(j)% moni_cntr(i))
          case (:-4)
            c_cntr (j) = ''
          case (NONE_CNTR)
            c_cntr (j) = '   none '
          case (ANY_CNTR)
            c_cntr (j) = '    any '
          case (UKN_CNTR)
            c_cntr (j) = 'unknown '
          case default
            write(c_cntr(j),'(i7,1x)') rept_use(j)% moni_cntr(i)
          end select
        end do
        write (6,'(1x,a10,  "  : ",99(a8  ,1x))') c,c_cntr
        c = ''
      end do
      write (6,'()')
    endif

  end subroutine print_rept_use
!------------------------------------------------------------------------------

  subroutine check_report_0 (state, head, nsubset, keep)
  !--------------------------------------------------------------------------
  ! preset report state variable
  !
  ! and perform a number of simple checks:
  !   CHK_SUBTYP  valid report subtype : Datenbankkennzahl
  !   CHK_TIME    time in valid range  : section 1 (synop)  time (if enabled)
  !                                      data base (cutoff) time
  !--------------------------------------------------------------------------
  type(t_use)  ,intent(inout)        :: state   ! state variable to set
  type(t_head) ,intent(in)           :: head    ! report header
  integer      ,intent(in)           :: nsubset ! number of subsets
  logical      ,intent(in) ,optional :: keep    ! keep state

    integer                     :: report_type   ! report type
    integer                     :: subtype_index ! subtype index
    type (t_rept_use)  ,pointer :: use
    type (t_rept_stat) ,pointer :: stat
    type (t_rept_stat) ,pointer :: stats
!   integer                     :: report_subtype
    integer                     :: my_use
    integer                     :: rule_num      ! rule no. setting 'use'
    integer                     :: statold       ! status before rejection
    character(len=12)           :: csubset =' subsets:   '
    character(len=20)           :: comment
    integer                     :: ns            ! number of subsets
    logical                     :: lkeep
    !---------------
    ! default status
    !---------------
    type(t_time) ,parameter :: m20 = t_time(0, 20*60) ! 20 minutes
    comment          = ''
    report_type      = head% obstype
    subtype_index    = head% idbk
    use              => rept_use                     (report_type)
    stat             => rept_stat                  (0,report_type)
    stats            => rept_stat      (subtype_index,report_type)
!   report_subtype   =  rept_char(report_type)% dbk (subtype_index)
    lkeep = .false.; if (present(keep)) lkeep = keep
    if (lkeep) then
      call decr_use (state, use% use(CHK_NONE), CHK_NONE)
    else
      state% check     =          CHK_NONE
      state% state     = use% use(CHK_NONE)
      state% flags     = 0
    endif
    !--------------
    ! simple checks
    !--------------
    if (stat% passive + stat% rejected + stat% active + 1 > use% max_proc) &
                          call decr_use (state, STAT_DISMISS, CHK_THIN)
    if (otime_exclude) then
      if   (                 use% time_s1) then
        if (head% time    <  use% time_b      ) call set_state (CHK_TIME)
        if (head% time    >= use% time_e      ) call set_state (CHK_TIME)
      endif
    else
      if   (                 use% time_s1) then
        if (head% time    <= use% time_b      ) call set_state (CHK_TIME)
        if (head% time    >  use% time_e      ) call set_state (CHK_TIME)
      endif
    endif
    if   (           .not. use% lhour       ) call set_state (CHK_TIME)
    if   (head% db_time /= zero_time        ) then
      if (invalid_time  /= use% time_cutoff .and. &
          head% db_time >  use% time_cutoff ) call set_state (CHK_TIME)
      if (head% db_time <  use% time_b - m20) then
        comment = 'db_time='//cyyyymmddhhmm(head% db_time)
        call decr_use  (state, STAT_DISMISS, CHK_TIME)
      endif
    endif
    !----------------
    ! kennzahl-filter
    !----------------
    call get_rule (head% modtype, head% buf_type, head% buf_subtype, &
                   head% dbkz, '',                                   &
                   obstype  = head% obstype,                         &
                   codetype = head% codetype,                        &
                   use      = my_use,                                &
                   rule_num = rule_num                               )
    select case (my_use)
    case (STAT_ACTIVE:)
    case (STAT_ABORT:STAT_ACTIVE-1)
      statold = state% state
      call decr_use (state, my_use, CHK_SUBTYP)
      if (statold /= state% state .and. rule_num > 0) then
        write (comment,'("rule: ",i3)') rule_num
      end if
    case default
      call finish('check_report_0','invalid value of my_use')
    end select

    !------------------
    ! generating center
    !------------------
    if (all (use%  use_cntr /= ANY_CNTR .and.   &
             use%  use_cntr /= head% center  )) &
      call decr_use (state, STAT_PASSIVE, CHK_NOTUSED)

    if (all (use% moni_cntr /= ANY_CNTR     .and. &
             use%  use_cntr /= head% center .and. &
             use% moni_cntr /= head% center  ))   &
      call decr_use (state, STAT_DISMISS, CHK_NOTUSED)

    !------------------
    ! update statistics
    !------------------
    ns = nsubset
    if (state% state > STAT_DISMISS) ns = 1
!call check_statistics('check_report_0: before')
    stat % processed = stat % processed + ns
    stats% processed = stats% processed + ns
    call update_statistics (stat, stats, int(state% state), ns)
!call check_statistics('check_report_0: after')
    select case (int(state% state))
    case (STAT_ABORT)
      write (csubset(10:12),'(i3)') ns
      call write_report (head, state, comment=csubset)
    case (STAT_DISMISS)
      if(comment=='') then
         write (csubset(10:12),'(i3)') ns
         comment = csubset
      end if
      call write_report (head, state, comment=comment)
      dism_stat   (state%check, report_type) = &
        dism_stat (state%check, report_type) + ns
    end select

  contains

    subroutine set_state (check)
    integer ,intent(in) :: check
      call decr_use (state, use% use (check), check)
    end subroutine set_state

  end subroutine check_report_0
!------------------------------------------------------------------------------
  subroutine dism_report_0 (state, head, check)
  type(t_use)  ,intent(out) :: state ! state variable to set
  type(t_head) ,intent(in)  :: head  ! report header
  integer      ,intent(in)  :: check ! reason for dismissal
    type (t_rept_stat) ,pointer :: stat
    type (t_rept_stat) ,pointer :: stats

    stat  => rept_stat          (0,head% obstype)
    stats => rept_stat (head% idbk,head% obstype)
    stat % dismissed = stat % dismissed + 1
    stats% dismissed = stats% dismissed + 1
    call decr_use (state, STAT_DISMISS, check)
    call write_report (head, state)
    dism_stat   (check, head% obstype) = &
      dism_stat (check, head% obstype) + 1
  end subroutine dism_report_0
!------------------------------------------------------------------------------

  subroutine check_report_1 (report)
  type(t_spot) ,intent(inout)        :: report        ! report
  !----------------------------------------------------------------------------
  ! perform a number of simple checks:
  !   CHK_TIME   time in valid range            (actual time)
  !   CHK_AREA   location in valid area         (specified area)
  !   CHK_DOMAIN coordinates in range           (-90..90; -180..360)
  !   CHK_HEIGHT location in valid height range (specified height range)
  !   CHK_RULE   complex rule                   (type/dbkz/statid/satid/latlon)
  !   CHK_THIN   max number of reports          (for testing)
  !----------------------------------------------------------------------------

    integer                     :: report_type
    integer                     :: subtype_index
    type (t_rept_use)  ,pointer :: use
    type (t_rept_stat) ,pointer :: stat
    type (t_rept_stat) ,pointer :: stats
    integer                     :: new_use
    integer                     :: old_use
    integer                     :: rule_use
    integer                     :: rule_num      ! rule no. setting 'use'
    real(wp)                    :: dlon, dlat, height
    character(len=20)           :: comment
    real(wp), parameter         :: TOL = 1.e-12_wp   ! rounding in BUFR decoding

    !--------------------
    ! set local variables
    !--------------------
    report_type   =  report% hd% obstype
    subtype_index =  report% hd% idbk
    use           => rept_use                (report_type)
    stat          => rept_stat             (0,report_type)
    stats         => rept_stat (subtype_index,report_type)
    comment       =  ''

    !---------------
    ! perform checks
    !---------------
    old_use = report% use% state

    !----------------------------------
    ! max. no of active+passive reports
    !----------------------------------
    if (stat% passive + stat% rejected + stat% active > use% max_proc) &
      call set_dismiss (CHK_THIN, '')
    !-----------------
    ! valid station id
    !-----------------
    if (report% statid == '') call set_state (CHK_INSDAT, 'no station id')
    !-----------------------------------
    ! validity of horizontal coordinates
    ! adjust longitudes above 180°E
    !-----------------------------------
    dlat    = report% col% c% dlat
    dlon    = report% col% c% dlon
    if (     dlon < -180._wp  &
        .or. dlon >  360._wp) then
                              call set_state (CHK_DOMAIN,'longitude')
    else if (dlon >  180._wp) then
      dlon = dlon - 360._wp
      report% col% c% dlon = dlon
    end if
    !---------------------------------
    ! catch rounding issues near poles
    !---------------------------------
    if (abs (dlat) >  90._wp) then
      if (abs (dlat) > 90._wp + TOL) then
                              call set_state (CHK_DOMAIN,'latitude' )
      else
        if (dlat > 90._wp) then
          dlat =  90._wp           ! Assign to North pole
        else
          dlat = -90._wp           ! Assign to South pole
        end if
        report% col% c% dlat = dlat
      end if
    end if
    !-----
    ! area
    !-----
    if (use% lat_nb >= use% lat_sb) then
      !-------------------------------------
      ! exclude observations OUT of the area
      !-------------------------------------
      if (dlat > use% lat_nb .or.  &
          dlat < use% lat_sb) call set_state (CHK_AREA,'lat OUT of area')
    else
      !-------------------------------------
      ! exclude observations WITHIN the area
      !-------------------------------------
      if (dlat > use% lat_nb .and. &
          dlat < use% lat_sb) call set_state (CHK_AREA,'lat WITHIN area')
    endif
    if (dlon > use% lon_eb) dlon = dlon - 360._wp
    if (dlon < use% lon_wb) dlon = dlon + 360._wp
    if (use% lon_eb >= use% lon_wb) then
      !-------------------------------------
      ! exclude observations OUT of the area
      !-------------------------------------
      if (dlon < use% lon_wb .or.  &
          dlon > use% lon_eb) call set_state (CHK_AREA,'lon OUT of area')
    else
      !-------------------------------------
      ! exclude observations WITHIN the area
      !-------------------------------------
      if (dlon < use% lon_wb .and. &
          dlon > use% lon_eb) call set_state (CHK_AREA,'lon WITHIN area')
    endif
    !-------
    ! height
    !-------
    if (report% ps /= invalid) then
       height = report% ps / 100.0_wp
       if (use% height_b >= use% height_t) then
          !------------------------------------------
          ! exclude observations OUTSIDE height range
          !------------------------------------------
          if (height > use% height_b .or.  &
              height < use% height_t     ) &
            call set_state (CHK_HEIGHT,'OUTSIDE height range')
       else
          !-----------------------------------------
          ! exclude observations WITHIN height range
          !-----------------------------------------
          if (height > use% height_b .and. &
              height < use% height_t     ) &
            call set_state (CHK_HEIGHT,'WITHIN height range')
       end if
    end if
    !------------
    ! actual time
    !------------
    if (otime_exclude) then
      if (report% actual_time <  use% time_b .or. &
          report% actual_time >= use% time_e      ) then
        call set_state (CHK_TIME,'head: '//cyyyymmddhhmm(report% hd% time))
      endif
    else
      if (report% actual_time <= use% time_b .or. &
          report% actual_time >  use% time_e      ) then
        call set_state (CHK_TIME,'head: '//cyyyymmddhhmm(report% hd% time))
      endif
    endif

!   if (report% actual_time /= report% hd% time) then
!     comment = 'head: '//cyyyymmddhhmm(report% hd% time)
!     call set_state (CHK_TIME)
!   endif

    !------
    ! rules
    !------
    call get_rule ( report% hd%     modtype,     &
                    report% hd%     buf_type,    &
                    report% hd%     buf_subtype, &
                    report% hd%     dbkz,        &
                    report%         statid,      &
     obstype =      report% hd%     obstype,     &
    codetype =      report% hd%     codetype,    &
       satid = int (report% hd%     satid),      &
      center =      report%         center_id,   &
         lon =      report% col% c% dlon,        &
         lat =      report% col% c% dlat,        &
         use =      rule_use,                    &
    rule_num =      rule_num                     )

    select case (rule_use)
    case (STAT_ACTIVE:)
    case (STAT_ABORT:STAT_ACTIVE-1)
      if (rule_use < report% use% state) then
         if (rule_num > 1) then
            write (comment,'("rule: ",i3)') rule_num
         else
            comment = 'type,dbkz,statid,satid,lat/lon'
         end if
      end if
      call decr_use (report% use, rule_use, CHK_RULE)
    case default
      call finish('check_report_1','invalid value of rule_use')
    end select

    !--------------
    ! update tables
    !--------------
    new_use = report% use% state
    if (new_use /= old_use) then
!call check_statistics('check_report_1: before')
      call update_statistics (stat, stats, new_use,  1)
      call update_statistics (stat, stats, old_use, -1)
!call check_statistics('check_report_1: after')
      select case (new_use)
      case (STAT_ABORT)
        call write_report (report, comment=comment)
      case (STAT_DISMISS)
        call write_report (report, comment=comment)
        dism_stat   (report% use% check, report_type) = &
          dism_stat (report% use% check, report_type) + 1
      end select
    endif

  contains

    subroutine set_state (check, com)
    integer          ,intent(in) :: check
    character(len=*) ,intent(in) :: com
      if (use% use (check) < report% use% state) comment = com
      call decr_use (report% use, use% use (check), check)
    end subroutine set_state

    subroutine set_dismiss (check, com)
    integer          ,intent(in) :: check
    character(len=*) ,intent(in) :: com
      if (STAT_DISMISS < report% use% state) comment = com
      call decr_use (report% use, STAT_DISMISS, check)
    end subroutine set_dismiss

  end subroutine check_report_1
!------------------------------------------------------------------------------
  subroutine decr_rpt_use (report, check, use, comment)
  type(t_spot)    ,intent(inout)        :: report        ! report
  integer         ,intent(in)           :: check         ! check
  integer         ,intent(in) ,optional :: use           ! new use
  character(len=*),intent(in) ,optional :: comment

    integer                     :: report_type
    integer                     :: subtype_index
    integer                     :: new_use
    integer                     :: old_use
    type (t_rept_stat) ,pointer :: stat
    type (t_rept_stat) ,pointer :: stats

    report_type   =  report% hd% obstype
    subtype_index =  report% hd% idbk
    stat          => rept_stat             (0,report_type)
    stats         => rept_stat (subtype_index,report_type)

    new_use = STAT_DEFAULT
    if (present (use))           new_use = use
    if (new_use == STAT_DEFAULT) new_use = rept_use (report_type)% use (check)
    old_use = report% use% state

    call decr_use (report% use, new_use, check)
    new_use = report% use% state

    if (new_use /= old_use) then
      report% comment = ''; if (present(comment)) report% comment = comment
!call check_statistics('decr_rpt_use: before')
      call update_statistics (stat, stats, new_use,  1)
      call update_statistics (stat, stats, old_use, -1)
!call check_statistics('decr_rpt_use: after')
      select case (new_use)
      case (STAT_ABORT)
        call write_report (report, comment)
      case (STAT_MERGED)
        call write_report (report, comment)
      case (STAT_DISMISS)
        call write_report (report, comment)
        dism_stat   (check, report_type) = &
          dism_stat (check, report_type) + 1
      end select
    endif

  end subroutine decr_rpt_use
!------------------------------------------------------------------------------
  subroutine update_statistics (stat, stats, use, incr)
  type(t_rept_stat) ,intent(inout) :: stat
  type(t_rept_stat) ,intent(inout) :: stats
  integer           ,intent(in)    :: use
  integer           ,intent(in)    :: incr
    select case (use)
    case (STAT_ABORT)
    case (STAT_OBS_ONLY)
      stat % obs_only  = stat % obs_only  + incr
      stats% obs_only  = stats% obs_only  + incr
    case (STAT_FORGET)
      stat % processed = stat % processed - incr
      stats% processed = stats% processed - incr
    case (STAT_DISMISS)
      stat % dismissed = stat % dismissed + incr
      stats% dismissed = stats% dismissed + incr
    case (STAT_MERGED)
      stat % merged    = stat % merged    + incr
      stats% merged    = stats% merged    + incr
    case (STAT_PASSIVE, STAT_PAS_REJ, STAT_NOTACTIVE)
      stat % passive   = stat % passive   + incr
      stats% passive   = stats% passive   + incr
    case (STAT_REJECTED)
      stat % rejected  = stat % rejected  + incr
      stats% rejected  = stats% rejected  + incr
    case (STAT_ACTIVE_0I:STAT_ACTIVE_1)
      stat % active    = stat % active    + incr
      stats% active    = stats% active    + incr
    case default
      if (use < 1 .or. use > n_stat) then
        write(0,*)  'update_statistics','invalid status:',use
        call finish('update_statistics','invalid status')
      else
        call finish('update_statistics','invalid status: '//stat_mnem(use))
      endif
    end select
  end subroutine update_statistics
!------------------------------------------------------------------------------
#if !defined (__SUNPRO_F95)
! Work around Sun f95 compiler bug; check_statistics is currently not needed

  subroutine check_statistics (comment)
  character(len=*) ,intent(in) :: comment
    integer :: i, j, l
    do i=ubound(rept_stat,1),lbound(rept_stat,1),-1
    do j=lbound(rept_stat,2),ubound(rept_stat,2)
      l = lost (rept_stat(i,j))
      if (l/=0) then
        write (0,*) dace% pe,'check_statistics: n,rpt,lost=',i,j,l,comment
        write (0,*) dace% pe,'check_statistics: n         =',rept_stat(i,j)
        call finish('check_statistics',comment)
      endif
    end do
    end do
  end subroutine check_statistics

#endif
!==============================================================================
  subroutine write_report (report, comment)
  !---------------------------------------------------
  ! Write an entry to report file 'monREPRT.nc'.
  ! Entries are collected in a bufr on each non-IO PE.
  ! They are communicated and written out later.
  !---------------------------------------------------
  type(t_spot)     ,intent(in)           :: report  ! report to print
  character(len=*) ,intent(in) ,optional :: comment ! comment
    if(present(comment)) then
      call write_head (report% hd,           &
                       report% use,          &
                       report% statid,       &
                       report% col% c% dlat, &
                       report% col% c% dlon, &
                       report% actual_time,  &
                       report% corme,        &
                               comment       )
    else
      call write_head (report% hd,           &
                       report% use,          &
                       report% statid,       &
                       report% col% c% dlat, &
                       report% col% c% dlon, &
                       report% actual_time,  &
                       report% corme,        &
                       report% comment       )
    endif
  end subroutine write_report
!------------------------------------------------------------------------------
  subroutine write_pending
  !------------------------------------------------------
  ! Write pending entries to report file 'monREPRT.nc'.
  ! Entries were collected in a buffer on each non-IO PE.
  ! They are now communicated and written out.
  !------------------------------------------------------
    integer :: pe, i, jpr, kpr (0:dace% npe-1)
    logical ::             lpr (0:dace% npe-1)
!   integer        ,allocatable :: ibuf(:)
!   integer        ,allocatable :: lens(:,:)
    character(len=*), parameter :: fmt_lost = &
         "('##### ',i10,&
         &' reports lost (buffer full, increase npr in mo_obs_tables) #####')"

!#if defined(_ECMWF_)
!if(.not.dace% lpio) write(6,*) 'write_pending: CURRENTLY DISABLED'
!ipr = 0
!return     !++++++ workaround for bug in  xlf 9.10.6  ?  +++++++++++
!#endif

    if (dace% lpio) then
      if (iu == 0) call open_report ()
      !----------------------------------------------------------
      ! Write pending reports from I/O processor and clear buffer
      !----------------------------------------------------------
      jpr = min(ipr,npr)
      do i = 1, jpr
        write (iu,'(A)') pending_report(i)(1:lenr(i))
      end do
      if (ipr>npr) write (iu,fmt_lost) ipr-npr
      ipr = 0
      if (seen_abort) then
        call finish ('write_pending','state=ABORT, comment='//trim(abort_msg))
      end if
    end if

    if (dace% npe == 1) return

FTRACE_BEGIN("write_pending:gather")
    call p_gather (ipr,        kpr, dace% pio)
    call p_gather (seen_abort, lpr, dace% pio)

!   if (dace% lpio) then
!      allocate (ibuf(npr * dace% npe))
!   else
!      allocate (ibuf(0))
!   end if
!   call p_gatherv (lenr, ibuf, root=dace% pio)
!   if (dace% lpio) then
!      allocate (lens(npr,0:dace% npe-1))
!      lens = reshape (ibuf, shape (lens))
!   end if
!   deallocate (ibuf)

    if (.not.dace% lpio) then
      !---------------------------------------------
      ! Send pending reports from non-I/O processors
      !---------------------------------------------
      jpr = min(ipr,npr)
      if (jpr > 0) then
          call p_send (pending_report(1:jpr), dace% pio, p_tag=42)
          call p_send (          lenr(1:jpr), dace% pio, p_tag=43)
      end if
      ipr = 0
      if (seen_abort) call p_send (abort_msg, dace% pio, 3)
    else
      !-----------------------------------------
      ! Collect pending reports on I/O processor
      !-----------------------------------------
      do pe = 0, dace% npe-1
        if (pe == dace% pio) cycle
        ipr = kpr(pe)
        jpr = min (ipr, npr)
        if (jpr > 0) then
          call p_recv (pending_report(1:jpr), pe,   p_tag=42)
          call p_recv (          lenr(1:jpr), pe,   p_tag=43)
!         lenr(1:jpr) = lens(1:jpr, pe)
        end if
        do i = 1, jpr
          write (iu,'(A)') pending_report(i)(1:lenr(i))
        end do
        if (ipr>npr) write (iu,fmt_lost) ipr-npr
        if (lpr(pe)) then
          call p_recv (abort_msg, pe, 3)
          call finish ('write_pending', &
                       'state=ABORT, comment='//trim(abort_msg))
        end if
      end do
      ipr = 0
    endif
FTRACE_END  ("write_pending:gather")
  end subroutine write_pending
!------------------------------------------------------------------------------
  subroutine write_head(head, status, statid, dlat, dlon, actual_time, corme, &
                        comment)
  !---------------------------------------------------
  ! Write an entry to report file 'monREPRT.nc'.
  ! Entries are collected in a bufr on each non-IO PE.
  ! They are communicated and written out later.
  !---------------------------------------------------
  type(t_head)     ,intent(in)           :: head         ! Report header
  type(t_use)      ,intent(in)           :: status       ! Report status
  character(len=*) ,intent(in) ,optional :: statid       ! station name
  real(wp)         ,intent(in) ,optional :: dlat         ! latitude
  real(wp)         ,intent(in) ,optional :: dlon         ! longitude
  type(t_time)     ,intent(in) ,optional :: actual_time  ! actual time
  integer          ,intent(in) ,optional :: corme        ! correction message
  character(len=*) ,intent(in) ,optional :: comment      ! comment

    character(len=12)              :: time, atime
    character(len=10)              :: name
    character(len=40)              :: com
    real(wp)                       :: lat, lon
    integer                        :: corm
    character(len(pending_report)) :: pr
    !--------------------
    ! optional parameters
    !--------------------
FTRACE_BEGIN("write_head:body")
    time =                                             cyyyymmddhhmm(head%  time)
    atime= ''        ;if (present(actual_time)) atime= cyyyymmddhhmm(actual_time)
    name = 'unknown' ;if (present(statid))      name = statid
    com  = ''        ;if (present(comment))     com  = comment
    lat  = -999._wp  ;if (present(dlat))        lat  = dlat
    lon  = -999._wp  ;if (present(dlon))        lon  = dlon
    corm = 0         ;if (present(corme))       corm = corme

    if (abs (lat) > 999._wp) lat = -999._wp
    if (abs (lon) > 999._wp) lon = -999._wp

    if (name == '')          name = 'unknown'

    write (pr,fm)   rept_char(  head% obstype)% mnem,   &!  8+1 +
                    name,                               &! 10+1 +
                    time,                               &! 12+1 +
                    lat,                                &!  9   +
                    lon,                                &!  9   +
                    head% codetype,                     &!  6   +
                    head% dbkz,                         &!  6   +
                    head% source,                       &!  4   +
                    head% record,                       &!  8+1 +
                    head% id,                           &!  9+1 +
                    chk(status% check)% mnem,           &! 12+1 +
                    stat_mnem (status% state),          &! 12   +
                               status% flags,           &! 32   +
                    atime,                              &! 12+1 +
                    corm,                               &!  1+1 +
                    trim(com)                            ! 40+1 = 199

    ipr = ipr + 1
    if (ipr <= npr) then
       pending_report(ipr) = pr
       lenr(ipr) = len_trim (pr)
    endif
    if (status% state == STAT_ABORT) then
       seen_abort = .true.
       abort_msg  = com
       if (ipr > npr) then
          pending_report(npr) = pr
          lenr(npr) = len_trim (pr)
       end if
    end if
FTRACE_END  ("write_head:body")
  end subroutine write_head
!------------------------------------------------------------------------------
  subroutine open_report ()
  !--------------------------------
  ! open report file 'monREPRT.nc'.
  !--------------------------------

    iu = get_unit_number()
    open (iu, file=path_file(aux,'monREPRT.nc'), status='replace')
    write (iu,'(a,a  )') '#'
    write (iu,'(a,a  )') '# Observation Flags'
    write (iu,'(a,a  )') '#'
    write (iu,'(a,i10)') '# experiment     = ',nex
    write (iu,'(a,i10)') '# run_type       = ',run_type
    write (iu,'(a,f5.2)')'# forecast hours =      ',fc_hours
    if (mod (fc_ref_time% secs, 3600) == 0 .and. &
        mod (ana_time%    secs, 3600) == 0       ) then
     write(iu,'(a,a  )') '# analysis  time = ',cyyyymmddhh  (ana_time)
     write(iu,'(a,a  )') '# reference time = ',cyyyymmddhh  (fc_ref_time)
    else
     write(iu,'(a,a  )') '# analysis  time = ',cyyyymmddhhmm(ana_time)
     write(iu,'(a,a  )') '# reference time = ',cyyyymmddhhmm(fc_ref_time)
    end if
    write (iu,'(a,a  )') '# run time (gmt) = ',cyyyymmddhhmm(run_time)
    write (iu,'(a,a  )') '# host           = ',host
    write (iu,'(a,a  )') '# user           = ',user
    write (iu,'(a)')     '#'
    write (iu,'(a)')        title
    write (iu,'(a)')     '#'
  end subroutine open_report
!==============================================================================
  function idb_dbk (dbkz, obt) result (idbk)
  !--------------------------------------------------------------------
  ! Determine the index idbk in table rept_char and rept_stat(idbk,obt)
  ! for given observation type (obt) and Datenbankkennzahl.
  !--------------------------------------------------------------------
  integer ,intent(in) :: dbkz ! Datenbankkennzahl
  integer ,intent(in) :: obt  ! observation type
  integer             :: idbk ! index in table rept_stat
    integer :: i
    !---------------------
    ! search dbkz in table
    !---------------------
    do i=1,ndb
      if (rept_char(obt)%dbk(i) == dbkz .or. rept_char(obt)%dbk(i) == -1) then
        idbk = i
        return
      endif
    end do
    !-------------------
    ! abort if not found
    !-------------------
#if defined (__SUNPRO_F95)
    ! Work around Sun f95 compiler bug: cannot have any write() here -> ICE!
#else
    write (0,*)'idb_dbk: dbkz for obt not found, ndb too small?'
    write (0,*)'idb_dbk: dbkz =', dbkz
    write (0,*)'idb_dbk: obt  =', obt, rept_char(obt)% mnem
    write (0,*)'idb_dbk: dbk()=', rept_char(obt)% dbk
    write (0,*)'idb_dbk: ndb  =', ndb
#endif
    call finish ('idb_dbk','dbkz for obt not found, ndb too small?')
  end function idb_dbk
!==============================================================================
  subroutine gather_rept_stat
  !--------------------------------
  ! gather flag statistics  on PE 0
  !--------------------------------

    type(t_rept_stat),allocatable :: rept (:,:) ! copy of remote rept_stat
    type(t_rept_stat),allocatable :: obsv (:,:) ! copy of remote rept_stat
    integer                       :: pe         ! processor element index

    if (.not.dace% lpio) then
      call p_send_rept_stat (rept_stat, dace% pio)
      call p_send_rept_stat (obsv_stat, dace% pio)
      rept_stat = rept_stat0
      obsv_stat = rept_stat0
    endif

    if (dace% lpio) then
      allocate (rept (size(rept_stat,1),size(rept_stat,2)))
      allocate (obsv (size(obsv_stat,1),size(obsv_stat,2)))
      do pe = 0, dace% npe-1
        if (pe /= dace% pio) then
          call p_recv_rept_stat (rept, pe)
          rept_stat% processed = rept_stat% processed + rept% processed
          rept_stat% merged    = rept_stat% merged    + rept% merged
          rept_stat% obs_only  = rept_stat% obs_only  + rept% obs_only
          rept_stat% dismissed = rept_stat% dismissed + rept% dismissed
          rept_stat% passive   = rept_stat% passive   + rept% passive
          rept_stat% rejected  = rept_stat% rejected  + rept% rejected
          rept_stat% active    = rept_stat% active    + rept% active
          rept_stat% accepted  = rept_stat% accepted  + rept% accepted
          call p_recv_rept_stat (obsv, pe)
          obsv_stat% processed = obsv_stat% processed + obsv% processed
          obsv_stat% merged    = obsv_stat% merged    + obsv% merged
          obsv_stat% obs_only  = obsv_stat% obs_only  + obsv% obs_only
          obsv_stat% dismissed = obsv_stat% dismissed + obsv% dismissed
          obsv_stat% passive   = obsv_stat% passive   + obsv% passive
          obsv_stat% rejected  = obsv_stat% rejected  + obsv% rejected
          obsv_stat% active    = obsv_stat% active    + obsv% active
          obsv_stat% accepted  = obsv_stat% accepted  + obsv% accepted
        endif
      end do
      deallocate (rept, obsv)
    endif

    dism_stat = p_sum  (dism_stat); if (.not.dace% lpio)  dism_stat = 0
    pass_stat = p_sum  (pass_stat); if (.not.dace% lpio)  pass_stat = 0
    reje_stat = p_sum  (reje_stat); if (.not.dace% lpio)  reje_stat = 0
   oonly_stat = p_sum (oonly_stat); if (.not.dace% lpio) oonly_stat = 0

  end subroutine gather_rept_stat
!==============================================================================
! send rept_stat
!---------------
#define VECTOR
#define RANK 2
#define DERIVED type(t_rept_stat),dimension(:,:)
#define p_send_DERIVED p_send_rept_stat
#undef  MPI_TYPE
#include "p_send.incf"
!==============================================================================
! recv rept_stat
!---------------
#define VECTOR
#define RANK 2
#define DERIVED type(t_rept_stat),dimension(:,:)
#define p_recv_DERIVED p_recv_rept_stat
#undef  MPI_TYPE
#include "p_recv.incf"
!==============================================================================
end module mo_obs_tables
