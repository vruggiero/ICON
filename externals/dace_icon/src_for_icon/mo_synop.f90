!
!+ surface observation operator (SYNOP,SHIP,BUOY,PAOB) specific routines
!
MODULE mo_synop
!
! Description:
!   Surface observation operator (SYNOP,SHIP,BUOY,PAOB) specific routines.
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
! V1_4         2009/03/26 Alexander Cress
!  Read scatterometer data from NetCDF
! V1_5         2009/05/25 Alexander Cress
!  fix scatterometer bias correction
! V1_6         2009/06/10 Alexander Cress
!  subroutine mleqc_qscat: check for invalid index to proposed wind
! V1_7         2009/08/24 Andreas Rhodin
!  include dbkz=10385 in array 'kzahl' (BUOY, new BUFR format)
! V1_8         2009/12/09 Andreas Rhodin
!  optimizations, cleanup, workaround bug in xlf V12.1
! V1_9         2010/04/20 Harald Anlauf
!  qscat_bcor: improve numerical consistency between platforms/compilers
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  read_synop_netcdf: reduce diagnostic output for netcdf_verb==0
! V1_13        2011/11/01 Harald Anlauf
!  changes for BUOYs, METARs, OCEANSAT-2, TEMP merge; optimisations for NEC SX
! V1_14        2011/11/08 Harald Anlauf
!  Cleanup quality control for OSCAT
! V1_15        2011/12/06 Andreas Rhodin
!  remove unused component t_datum% u
! V1_16        2011/12/09 Alexander Cress
!  SCATTEROMETER: store fov in spti% phase
! V1_20        2012-06-18 Harald Anlauf
!  read_synop_netcdf: if pmsl (mpppp) is reported as 0, set to invalid value
!  Oscat: reject data if windspeed is <3m/s or >30m/s
! V1_22        2013-02-13 Harald Anlauf
!  Disambiguation of synops and ships via bufr_type
!  Enable ships reporting as "SHIP"
!  Bugfix for "new buoys" (KZ=10385) in the Antarctic
!  Separate treatment of ship and buoys
!  Additional diagnostics for double entries
! V1_24        2013/04/26 Andreas Rhodin
!  SCATTerometer: store satid in spot%hd%satid (required for correct thinning)
! V1_27        2013-11-08 Harald Anlauf
!  Dismiss SYNOPs with lat=0,lon=0
! V1_28        2014/02/26 Andreas Rhodin
!  new interface to new_int
! V1_29        2014/04/02 Harald Anlauf
!  Reduce default verbosity for scatterometer data
! V1_31        2014-08-21 Andreas Rhodin
!  changes for ECMWF SYNOP (BUFR2)NetCDF input
! V1_37        2014-12-23 Alexander Cress
!  preference BUFR SYNOP reports
! V1_40        2015-02-27 Harald Anlauf
!  changes for new SHIP BUFR reports (A.Cress)
!  bugfix for reading scatterometer reports in NetCDF (R.Faulwetter)
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for COSMO MEC
! V1_45        2015-12-15 Harald Anlauf
!  Implement windspeed observations for Jason-2; Preparations for SARAL/Altika
! V1_46        2016-02-05 Harald Anlauf
!  Improve handling of SYNOP wind obs; enhance check of multiple reports;
!  check against proper missing value for station height;
!  Fix DRIBU station id; extend SYNOP/DRIBU code for monitoring of T2m
! V1_47        2016-06-06 Harald Anlauf
!  scatterometer bias correction;
!  extend altimetry data QC (Jason-2) to ECMWF recommendations;
!  handle Jason-3, Sentinel 3A
! V1_48        2016-10-06 Andreas Rhodin
!  passive monitoring for RR,GUST,T2M,FF,DD
!  RH2M assimilation: option to use model first guess
! V1_49        2016-10-25 Andreas Rhodin
!  extend sanity check for reported SYNOP pressure; bug fix for MEC
! V1_50        2017-01-09 Harald Anlauf
!  add checks for altimeter; changes to run 3dvar with COSMO/COMET data.
! V1_51        2017-02-24 Andreas Rhodin
!  changes to pass observations from COSMO fof-files through the 3dvar
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD 2000-2008  original source
! Oliver Schmid   DWD       2005       new obs data type
! Gerhard Paul    DWD       2008       input in NetCDF format
! Harald Anlauf   DWD       2008       modifications for SX8
!==============================================================================

!-------------------------------------------------
! uncomment to exit in case of unknown BUFR codes:
!
!#define CHECKCODES
!-------------------------------------------------

  !=============
  ! modules used
  !=============
  use mo_kind,       only: wp, sp, i2, i8  ! real, integer kind parameters
  use mo_namelist,   only: position_nml,  &! position namelist
                           nnml,          &! namelist Fortran unit number
                           POSITIONED      ! ok    code from position_nml
  use mo_mpi_dace,   only: dace,          &! MPI group info
                           p_bcast         ! broadcast routine
  use mo_dace_string,only: char3           ! conversion: int -> char(len=3)
  use mo_obs_set,    only: t_obs_block     ! obs. data type
!                          t_obs_set       ! observation data derived type
  use mo_exception,  only: finish          ! abort on error condition
  use mo_atm_state,  only: t_atm           ! atmospheric state data type
  use mo_t_datum,    only: t_datum,       &! data type for one observed datum
                           inv_datum,     &! constant for invalid datum
                           set_datum,     &! set t_datum% o  (observed value)
                           set_qbits,     &! set t_datum% qc (quality bits)
                           cmp_datum,     &! compare data
                           merge_datum,   &! merge data
                           count_datum,   &! count relative data content
                           print,         &! generic print routine
                           rvind,         &! invalid value
                           SRC_DER,       &! derived  value flag
                           QC_OK,         &! quality flag: OK
                           QC_MISS,       &!               missing
                           QC_INCON,      &!               inconsistent
                           QC_NOUSE,      &!               do not use
                           QC_CLIM         !               climatological range
  use mo_t_use,      only: t_use,         &! status variable data type
                           use_0,         &! default values of type use
                           decr_use,      &! decrease the state of a datum
                           STAT_ACTIVE,   &!
                           STAT_DISMISS,  &!
                           STAT_PASSIVE,  &!
                           STAT_OBS_ONLY, &!
                           CHK_NOTUSED,   &!
                           CHK_CORR,      &!
                           CHK_CORRERR,   &!
                           CHK_REDUNDANT, &!
                           CHK_DBLERR,    &!
                           CHK_INSDAT,    &!
                           CHK_DOMAIN,    &!
                           CHK_TIME,      &!
                           CHK_HEIGHT,    &! height range
                           CHK_THIN,      &! thinning flag
                           CHK_QI,        &! quality index
!                          CHK_BIASCOR,   &! flag for no bias correction
                           CHK_BLACKLIST, &! blacklisting flag
                           CHK_NOIMPL      !
  use mo_wmo_tables, only: WMO0_ECMWF,    &! generating center
                           WMO0_DWD        !
  use mo_obs_tables, only: decr_rpt_use,  &!
                           idb_dbk,       &! index in table rept_stat
!                          rept_char,     &! observation type characteristics
                           check_report_0,&! init. flags, standard checks
                           check_report_1  !
  use mo_t_obs,      only: t_obs,         &! observation data type
                           t_spot,        &! spot data type
                           t_head,        &! observation data type
                           derive_dbkz,   &! derive DBKZ if not present
                           new_par,       &! request space in parameter array
                           new_spot,      &! request space in metadata  array
                           new_obs,       &! request space in observation array
                           new_int,       &! request space in int.obs.    array
                           set_xuv,       &! set unit vectors, solar zenith
                           invalid,       &! invalid observation value (real)
                           empty_spot,    &! empty spot
                           set_vqc_insitu,&! subroutine to set VQC bounds
                           shrink_report, &! remove passive observations
                           monitor_ff,    &! flag to monitor wind speed
                           monitor_dd,    &! flag to monitor wind direction
                           source,        &! list   of Report source files
                           SYNOP,         &! module flag value
!                          CHR_ID,        &! H is the identity operator
                           CHR_NONL,      &! H is a non-linear operator
                           TSK_INIT,      &! task flag: initialize modules
                           TSK_READ,      &!   read observations
                           TSK_SET_CHR,   &!   set observation characteristics
                           TSK_SETUP_COLS,&!   determine model columns required
                           TSK_SETUP_FULL,&!   setup description of PSAS-space
                           TSK_SETUP_FUL0,&!   not used for SYNOP
                           TSK_SHRINK,    &!  release unused obs. in report
                           TSK_R,         &!   setup observational errors
                           TSK_K,         &!   evaluate linear operator
                           TSK_Y,         &!   evaluate nonlinear operator
                           ITY_ICOL,      &! interpolation type: column
                           netcdf_verb,   &! verbosity of NetCDF decoding
                           OBS_H,OBS_U,OBS_V,OBS_T,OBS_TV,OBS_RH,OBS_DUM
  use mo_wigos,      only: operator (==)   ! Same WIGOS station ID?
  use mo_obs_trange, only: apply_trange    ! apply time range operator
  use mo_synop_bc   ,only: biascor_mode,  &! mode used for bias correction
                           reg_param,     &! regularization parameter for non-linear bc
                           t_decay,       &! biasc. accum. decay time (days)
                           n_required,    &! number of entries required for bc
                           bc_fallback,   &! action if biasc-file not present
                           aggregate,     &! aggregate rain rate observations
                           bc_t2m,        &! apply b/c for t2m (nonlinear bc)?
                           bc_rh2m,       &! apply b/c for rh2m (nonlinear bc)?
                           bc_synop,      &! apply b/c for SYNOP Land?
                           bc_ship,       &! apply b/c for SYNOP Ship?
                           bc_buoy,       &! apply b/c for DRIBU?
                           bc_metar,      &! apply b/c for METAR?
                           chk_zstation,  &! check station height consistency
                           z0_bc_synop,   &! threshold parameter for SYNOP Land
                           z0_bc_ship,    &! threshold parameter for SYNOP Ship
                           z0_bc_buoy,    &! threshold parameter for DRIBU
                           z0_bc_metar,   &! threshold parameter for METAR
                           compress_ag     ! "compress" aggregated data?
  use mo_fdbk_tables,only: VN_U, VN_U10M, &!          wind component code
                           VN_V, VN_V10M, &!                         code
                           VN_RH,VN_RH2M, &!          rh             code
                           VN_Z,          &!          geopotential   code
                           VN_FF,         &!          wind speed     code
                           VN_DD,         &!          wind direction code
                           VN_GUST,       &!          gust speed     code
                           VN_T,          &!          temperature    code
                           VN_TD,         &!          dewpoint temp. code
                           VN_T2M,        &!          2m temperature code
                           VN_TD2M,       &!          2m dewpoint t. code
                           VN_TSEA,       &!          sea/water temp.code
                           VN_TRTR,       &!          time range (h) code
                           VN_RR,         &!          rain rate      code
                           VN_P,          &!                pressure code
                           VN_PS,         &!        surface pressure code
                           VN_HEIGHT,     &!          height         code
                           VN_HOSAG,      &!     height above ground code
                           VN_DEPTH,      &!     depth below surface code
                           VN_PTEND,      &!       pressure tendency code
                           VN_N, VN_N_L,  &! cloud amount
                           VN_N_M, VN_N_H,&! cloud amount
                           VN_NH,         &! cloud height
                           VN_CEIL,       &! cloud ceiling
                           VN_VV,         &! visibility
                           VN_RAD_GL,     &! global    radiation
                           VN_RAD_DF,     &! diffusive radiation
                           VN_RAD_LW,     &! long wave radiation
                           VN_SDEPTH,     &! snow depth
                           VN_WW,         &! present weather code    code
                           VN_GCLG,       &! general cloud group     code
                           VN_ICLG,       &! individual cloud layer group
                           VN_NUM,        &! ordinal (channel) number   (  )
                           OT_SYNOP,      &! data type specification
                           OT_DRIBU,      &! data type specification
                           VN_TMIN,VN_TMAX !
  use mo_t_col,      only: t_cols,        &!   model columns data type
                           t_col,         &!
                           COL_TV,COL_Q,  &!
                           COL_RH,        &!
                           COL_UV,        &!
                           COL_CLC,       &!
                           COL_GEO,       &!
                           COL_RANGE,     &!
                           nt              ! number of time range slots
  use mo_cosmo_conv, only: octa_percent,  &! convert cloud cover (%) to octa
                           present_clc     ! flag indicating if CLC is present
  use mo_time,       only: init_time,     &! initialise time data type
                           operator(<),   &! compare  times
                           operator(>=),  &! compare  times
                           operator(==),  &! compare  times
                           operator (-),  &! subtract times
                           hours,         &! derive hours   from time variable
                           minutes,       &! derive minutes from time variable
                           chhmm,         &!
                           print,         &! generic print routine
                           cyyyymmddhhmmss ! derive string from time
  use mo_dec_matrix, only: t_vector_segm, &! vector segment data type
                           mp              ! kind parameter for matrix elements
  use mo_physics,    only: esw_t,         &! sat.vapour pressure over water
                           t_esw,         &! inverse esw_t
                           desw_dt,       &! derivative: d esw / d t
                           fd_uv,         &! derive ff, dd from u,v
                           gacc,          &! gravity acceleration
                           R,             &! gas constant
!                          c_p,           &! specific heat [J /kg /K]
                           lapse_cl,      &! Climatological lapse rate [K/m]
                           d2r             ! pi/180, factor degree  -> radians
  use mo_cntrlvar,   only: trh_tvgh        ! generalized humidity transform
  use mo_usstd,      only: p_h_usstd,     &! p from h(gpm) US std.atm.
                           p_r             ! reference pressure
  use mo_grid_intpol,only: idx_init
  use mo_bufr_dwd,   only: t_bufr,        &! BUFR record data type
                           inv_bufr,      &! indicator for invalid value
                           nvind,         &! integer missing value indicator
#ifdef  CHECKCODES
                           bufr_get_entry_texts, &!
                           bufr_get_entry_units, &!
#endif
                           bufr_print_sections,  &! print sections 0-4
                           bufr_print_subset,    &! print body of BUFR record
                           bufr_get_character     ! decode character entry
! use mo_obs_rules,  only: get_rule,      &! routine to get a rule
!                          iud,           &! undefined integer value
!                          rud,           &! undefined real    value
!                          t_set           ! result data type
  use mo_obstypes,only:    t_obsid,       &! observation id table entry
                           obstype_dbkz,  &! derive obsids from dbkz
                           obstype_bufr    ! derive obsids from bufr type
  use mo_p_output,   only: oline,         &! output line buffer
                           iol,           &! index of next line to write
                           nextline        ! increment line number
  use mo_obs_err,    only: synop_obs_err   ! set SYNOP,DRIBU,PAOB obs.err.
  use mo_run_params, only: ana_time,      &! analysis time
                           method          ! check for 'MEC','GMESTAT'
  !--------------------------------------
  ! obs data header info read from NetCDF
  !--------------------------------------
  use mo_head_netcdf,only: ncid,          &! NetCDF file id
                           dimids_max,    &! max number of NetCDF dimension ids
                           imissing,      &! NetCDF _FillValue for integer
                           rmissing,      &! NetCDF _FillValue for reals
                           s2ikz,         &! DWD-internal classifier
                           s1cat,         &! data category
                           s1catls,       &! local data sub category
                           s1cent,        &! data centre
                           stime,         &! header observation time (section1)
                           db_time,       &! data bank time
                           s1cents,       &! data sub centre
                           s1updat,       &! update sequence no.
                           mlah,          &! latitude
                           mloh,          &! longitude
                           obs_time,      &! body observation time
                           ystidn,        &! any type of station identifier as variable
                           istidn,        &! WMO numeric station number combined
                           lwsi,          &! WIGOS station identifier valid
                           wsi,           &! WIGOS station id stored representation
                           wsihash         ! Station id hash
  !---------------------
  ! netCDF f90 interface
  !---------------------
  use netcdf,        only: nf90_Inquire,          &!
                           nf90_Inquire_Dimension,&!
                           nf90_Inquire_Variable, &!
                           nf90_inq_varid,        &!
                           nf90_get_var,          &!
                           NF90_MAX_NAME,         &!
                           NF90_NOERR,            &!
                           NF90_INT,              &!
                           NF90_BYTE,             &!
                           NF90_SHORT,            &!
                           NF90_FLOAT,            &!
                           NF90_DOUBLE             !
  implicit none
!------------------------------------------------------------------------------
  !================
  ! public entities
  !================
  private
  public :: process_synop
  public :: read_synop_bufr
  public :: read_synop_netcdf ! read observations from netCDF file
  public :: t_synop           ! SYNOP observation          data type
  public :: inva_syn          ! invalid SYNOP obs          t_synop constant
  public :: check_store_synop ! subroutine passed to BUFR read routine
  public :: dtdzp             ! temperature gradient for p extrapolation
  public :: zbs               ! temperature gradient extrapolation
  public :: use_ps_model      ! invalid value: use model surface pressure
  public :: version           ! operator version
  public :: read_synop_nml    ! read namelist /SYNOP_OBS/
!------------------------------------------------------------------------------
  !======================
  ! data type definitions
  !======================
  type t_synop
    !-------------------------------------------------------------
    integer        :: p_used ! p_s, p_gp, p_msl
    type (t_datum) :: total  ! quality control flag
    type (t_datum) :: z_ls   ! height of land surface     [m]
    !-------------------------------------------------------------
    type (t_datum) :: zr     ! height of station reported [m]
    type (t_datum) :: z      ! height of station used     [m]
    type (t_datum) :: ps     ! pressure at station height [Pa]
    !-------------------------------------------------------------
    type (t_datum) :: p_msl  ! pressure mean sea level    [Pa]
    type (t_datum) :: z_msl  ! height   mean sea level    [m]       ! always 0
    !-------------------------------------------------------------
    type (t_datum) :: p      ! pressure level             [Pa]
    type (t_datum) :: gp     ! geopotential               [gpm]
    !-------------------------------------------------------------
    type (t_datum) :: p_as   ! pressure     to assimilate [Pa]
    type (t_datum) :: z_as   ! geopotential to assimilate [gpm]
    !-------------------------------------------------------------
    type (t_datum) :: p_ref  ! reference pressure         [Pa]
    type (t_datum) :: h_ts   ! height of t/td/rh sensor   [m]
    type (t_datum) :: t      ! temperature                [K]
    type (t_datum) :: td     ! dewpoint temperature       [K]
    type (t_datum) :: rh     ! relative humidity          [ ]
    !-------------------------------------------------------------
    type (t_datum) :: h_ws   ! height of wind sensor      [m]
    type (t_datum) :: ff     ! wind speed                 [m/s]
    type (t_datum) :: dd     ! wind direction             [deg]
    type (t_datum) :: uu     ! wind component             [m/s]
    type (t_datum) :: vv     ! wind component             [m/s]
    !-------------------------------------------------------------
    type (t_datum) :: tsea   ! sea/water temperature      [K]
    type (t_datum) :: h_tsea ! depth of sensor below surf.[m]
    !-------------------------------------------------------------
    real           :: p_tend= rvind ! 3 hour pressure change     [Pa]
    real           :: gust1 = -99.  ! maximum wind speed(gusts)  [m/s] 1 hour
    real           :: gust3 = -99.  ! maximum wind speed(gusts)  [m/s] 3 hour
    real           :: gust6 = -99.  ! maximum wind speed(gusts)  [m/s] 6 hour
    real           :: tmax1 = rvind ! maximum temperature        [K]  -1 hour
    real           :: tmin1 = rvind ! minimum temperature        [K]  -1 hour
    real           :: tmax  = rvind ! maximum temperature        [K] -12 hour
    real           :: tmin  = rvind ! minimum temperature        [K] -12 hour
    integer        :: n     = -1    ! cloud cover (total)        [%]
!   integer        :: vs    = -1    ! vertical significance (surface observ.)
!   real           :: nn    = -99.  ! (lowest) cloud amount      [octas]
    real           :: nh    = rvind ! height of base of cloud    [m]
    real           :: ceil  = rvind ! cloud ceiling above ground [m]
    real           :: clcl  = -99.  ! low    cloud amount        [octas]
    real           :: clcm  = -99.  ! middle cloud amount        [octas]
    real           :: clch  = -99.  ! high   cloud amount        [octas]
    integer        :: gwg   = nvind ! general cloud group
    integer        :: icl1  = nvind ! individual cloud layer 1
    integer        :: icl2  = nvind ! individual cloud layer 2
    integer        :: icl3  = nvind ! individual cloud layer 3
    integer        :: icl4  = nvind ! individual cloud layer 4
    integer        :: ww    = nvind ! present weather
    integer        :: w     = nvind ! past weather (2) i.e.  6 hour
    real           :: vis   = rvind ! horizontal visibility      [m]
    integer        :: e     = nvind ! state of ground (with or without snow)
    real           :: sdepth= rvind ! total snow depth           [m]
    real           :: rr1  = -99.   ! total precipitation        [kg/m**2]  -1 hour
    real           :: rr3  = -99.   ! total precipitation        [kg/m**2]  -3 hour
    real           :: rr6  = -99.   ! total precipitation        [kg/m**2]  -6 hour
    real           :: rr12 = -99.   ! total precipitation        [kg/m**2] -12 hour
    real           :: rr24 = -99.   ! total precipitation        [kg/m**2] -24 hour
    real           :: lwr1 = -99.   ! long wave radiation        [J/m**2]   -1 hour
    real           :: gsr1 = -99.   ! global solar radiation     [J/m**2]   -1 hour
    real           :: gsr3 = -99.   ! global solar radiation     [J/m**2]   -3 hour
    real           :: gsr6 = -99.   ! global solar radiation     [J/m**2]   -6 hour
    real           :: dsr1 = -99.   ! diffuse solar radiation    [J/m**2]   -1 hour
    real           :: dsr3 = -99.   ! diffuse solar radiation    [J/m**2]   -3 hour
    real           :: dsr6 = -99.   ! diffuse solar radiation    [J/m**2]   -6 hour
    integer        :: qobl = nvind  ! Quality of buoy location
    integer        :: nvs  =  -1    ! Speed of motion (SHIP)     [m/s]
    !-------------------------------------------------------------
  end type t_synop

  integer, parameter :: P_NONE = 0
  integer, parameter :: P_S    = 1
  integer, parameter :: P_GP   = 2
  integer, parameter :: P_MSL  = 3

  type (t_synop) ,parameter :: inva_syn =                          &
        t_synop(P_NONE   ,                                         &
                inv_datum,inv_datum,inv_datum,inv_datum,inv_datum, &
                inv_datum,inv_datum,inv_datum,inv_datum,inv_datum, &
                inv_datum,inv_datum,inv_datum,inv_datum,inv_datum, &
                inv_datum,inv_datum,inv_datum,inv_datum,inv_datum, &
                inv_datum,inv_datum,inv_datum,                     &
!               inv_datum,inv_datum,inv_datum,inv_datum,           &! ff2,..,vv2
                rvind    ,-99.     ,-99.     ,-99.               , &! p_tend,..,gust6
                rvind    ,rvind    ,rvind    ,rvind    ,-1       , &! tmax/tmin, n
                rvind    ,rvind    ,-99.     ,-99.     ,-99.     , &! nh,..,clch
                nvind    ,nvind    ,nvind    ,nvind    ,nvind    , &! gwg,...,icl4
                nvind    ,nvind    ,rvind    ,nvind    ,rvind    , &! ww,..,sdepth
                -99.     ,-99.     ,-99.     ,-99.     ,-99.     , &! rr1,..,rr24
                -99.     ,-99.     ,-99.     ,-99.               , &! lwr1,..,gsr6
                -99.     ,-99.     ,-99.                         , &! dsr1,..,dsr6
                nvind    ,-1                                       )

! integer, parameter :: synop_int_size = size (transfer (inva_syn,(/0/)))
  integer            :: synop_int_size = 0
  real(sp),parameter :: use_ps_model = 989898._sp ! use model surface pressure
  integer, parameter :: CBH_CLR = 16000 ! cloud base height assigned at clear sky
                                        !   (either in obs or in model equival.
                                        !    should be < 20920)
!------------------------------------------------------------------------------
  !===============================
  ! List of "Datenbank-Kennzahlen"
  !===============================
  integer, parameter :: kz_synop(7) = [   0, 128, 156, 10000, 10128, 10015, 10143 ]
  integer, parameter :: kz_nat  (2) = [   5,      10005 ] ! German national BUFR template
  integer, parameter :: kz_swis (2) = [ 170,      10170 ] ! Road weather stations
  integer, parameter :: kz_ship (4) = [ 256, 384, 10256, 10384 ]
  integer, parameter :: kz_metar(1) = [   1 ]
  integer, parameter :: kz_buoy (2) = [ 385, 10385 ]
  integer, parameter :: kz_cman (1) = [ 10158 ]
  integer, parameter :: kz_metar_usa(1) = [ 131 ]
  integer, parameter :: kz_seemhem (2) = [ 166, 10150 ]
! integer, parameter :: kz_scat (5) = [ 1697, 1699, 1700, 1701, 1702 ]
!------------------------------------------------------------------------------
  !-------------------
  ! Namelist SYNOP_OBS
  !-------------------
  logical  :: use_t2     = .false.   ! 2 m temperature
  logical  :: use_rh     = .false.   ! 2 m humidity
  logical  :: mon_td     = .true.    ! 2 m dewpoint temperature (monitoring)
  logical  :: mon_tr     = .true.    ! time range variables (tmin,tmax,gust,rain)
  logical  :: mon_cl     = .true.    ! cloud amount, etc.
  logical  :: mon_cwoo   = .false.   ! obs only cloud + weather var. (no oper.)
  logical  :: mon_tsea   = .false.   ! sea temperature monitoring (SHIP,BUOY)
  logical  :: mon_snow   = .false.   ! monitor snow depth (SYNOP)
  logical  :: mon_ps     = .false.   ! monitor off-synoptic ps
  logical  :: mon_t2     = .false.   ! monitor off-synoptic t2m/td2m
  logical  :: use_vs     = .true.    ! 10 m wind over sea
  logical  :: use_vl     = .false.   ! 10 m wind over land (tropics, z < 150 m)
  logical  :: use_vl_gl  = .false.   ! 10 m wind over land (everywhere)
  logical  :: use_p      = .false.   ! p(h) is observable
  logical  :: use_h      = .true.    ! h(p) is observable
  logical  :: ag_rr_1h   = .false.   ! aggregate 1h rain rate
  logical  :: ag_rr_3h   = .true.    ! aggregate 3h rain rate
  logical  :: ag_rr_6h   = .true.    ! aggregate 6h rain rate
  logical  :: ag_rd_1h   = .false.   ! aggregate 1h radiation
  logical  :: ag_rd_3h   = .true.    ! aggregate 3h radiation
  logical  :: ag_rd_6h   = .true.    ! aggregate 6h radiation
  logical  :: ag_gust_3h = .true.    ! aggregate 3h gust
  logical  :: ag_gust_6h = .true.    ! aggregate 6h gust
  logical  :: ag_tmaxmin = .true.    ! aggregate tmax/tmin
  logical  :: ag_combine = .true.    ! combine rr and gust from different BUFRs
  integer  :: vprc_cloud = 1         ! version of cloud obs pre-processing
  real(wp) :: chk_ps(2)  = 0._wp     ! consistency check p_station vs. p_msl
  integer  :: chk_strict = 0         ! consistency check strictness level (0-3)
  logical  :: prt_data   = .false.   ! print data
  integer  :: zls_to_z   = 3         ! set z from land surface   height (see
  integer  :: zmsl_to_z  = 4         ! set z from mean sea level height  below)
  real(wp) :: p_d_msl    = 100       ! max difference between ps and pmsl [hPa]
  real(wp) :: dtdzp      = lapse_cl  ! temperature gradient for p extrapolation
  real(wp) :: ssdmax_v10m= -1._wp    ! max. SSO standard dev. for 10m wind [m]
  real(wp) :: dz_max_v10m= -1._wp    ! max. height diff. to model surface [m]
  real(wp) :: ssdmax_rh2m= -1._wp    ! max. SSO standard dev. for 2m humi.[m]
  real(wp) :: dz_max_rh2m= -1._wp    ! max. height diff. to model surface [m]
  real(wp) :: ssdmax_t2m = -1._wp    ! max. SSO standard dev. for 2m temp.[m]
  real(wp) :: dz_max_t2m = -1._wp    ! max. height diff. to model surface [m]
  integer  :: p_land (3) = (/1,2,3/) ! pressure to use: ...
  integer  :: p_ship (3) = (/3,0,0/) ! 0=none 1=p_s, 2=p_gp, 3=p_msl
  integer  :: p_buoy (3) = (/3,0,0/) !
  integer  :: p_metar(3) = (/1,3,0/) !
  integer  :: chk_dbl    = 0         ! extended check for double entries
  real(wp) :: chk_dbl_dt =  0._wp    ! tolerance for time difference [min]
  integer  :: chk_cldinfo= 0         ! mask for cloud information checks:
                                     ! bit 0: treat missing vertical significance
                                     ! bit 1: handle "no cloud" in cloud type
  integer  :: chk_ship   = 0         ! extended check for ship reports
                                     ! 0=none, 1=check if buoy or rig,2=dup.check
  integer  :: chk_zb_mode= 0         ! consistency check for barometer height
                                     ! 0=none, 1=reject suspicious, 2=correct
  real(sp) :: chk_zb_max =  2._sp    ! barometer height suspicious below 2m
  real(sp) :: dzs_zb_max = 20._sp    ! max. diff. |z_station-z_barometer| [m]
  integer  :: version    = 0         ! operator version:
                                     ! 0: old rh2m from 10m model layer
                                     ! 1: from model 2m for synop-land and buoy
                                     ! 2: always rh2m from 2m model diagnostics
                                     ! 3: revised adjoint for t2m, rh2m
                                     ! 4: 3 + interpolation between 2m/10m
                                     ! 5: 3 + use 2m land-tile average values
                                     ! 6: 4 + use 2m land-tile average values
  integer  :: verbose    = 1         ! Verbosity level of consistency checks
  integer  :: v_buoy     = 0         ! operator version for wind, buoys:
                                     ! 0: use 10m wind
                                     ! 1: use logarithmic wind profile
  real(wp) :: buoy_h0_v    = 4.0_wp  ! buoy: default height of wind sensor  [m]
  real(sp) :: buoy_hr_v(2) = [2, 10] ! buoy: valid height range of w.sensor [m]
  real(sp) :: buoy_hr_t(2) = [2, 10] ! buoy: valid height range of t.sensor [m]
  real(sp) :: hr_t2m   (2) = [0, 99] ! land: valid height range of t.sensor [m]
  integer  :: chk_tmaxmin  = 0       ! mask for Tmax/Tmin checks (0-7)
  logical  :: bugfix_t2m   = .true.  ! bugfix for t2m, obs.operator setup
  logical  :: bugfix_rh2m  = .true.  ! bugfix for rh2m, obs.operator setup
  !
  ! values for zls_to_z or zmsl_to_z :
  !   0: use station height
  !   1: use surface height if station height is not present
  !   2: use surface height if surface height < station height
  !   3: use surface height if surface height > station height
  !   4: use surface height if surface height is present
  !   5: use surface height
  !
  namelist  /SYNOP_OBS/ use_t2, mon_td, use_rh, use_vs, use_vl, use_p, use_h, &
                        mon_tr, mon_cl, mon_cwoo, mon_tsea, p_d_msl, p_land,  &
                        p_ship, p_buoy, p_metar, zls_to_z, zmsl_to_z, dtdzp,  &
                        chk_ps, chk_strict, chk_dbl, prt_data, verbose,       &
                        use_vl_gl, version, ag_rr_1h, ag_rr_3h, ag_rr_6h,     &
                        ag_rd_1h, ag_rd_3h, ag_rd_6h, ag_gust_3h, ag_gust_6h, &
                        ag_tmaxmin, ag_combine, vprc_cloud,                   &
                        biascor_mode, t_decay, n_required, bc_fallback,       &
                        v_buoy, buoy_h0_v, buoy_hr_v, buoy_hr_t, hr_t2m,      &
                        ssdmax_v10m, dz_max_v10m, ssdmax_rh2m, dz_max_rh2m,   &
                        chk_zb_mode, chk_zb_max, dzs_zb_max, chk_dbl_dt,      &
                        bc_synop, bc_ship, bc_buoy, bc_metar, chk_zstation,   &
                        z0_bc_synop, z0_bc_ship, z0_bc_buoy, z0_bc_metar,     &
                        compress_ag, chk_ship, mon_ps, mon_t2, chk_cldinfo,   &
                        reg_param, bc_t2m, bc_rh2m, ssdmax_t2m, dz_max_t2m,   &
                        chk_tmaxmin, bugfix_t2m, mon_snow, bugfix_rh2m
  !-------------------
  ! derived quantities
  !-------------------
  real(wp) :: p_stat = 10000._wp ! p_d_msl [Pa]
  real(wp) :: z_stat =   800._wp ! p_d_msl [gpm]
  !-------------------------------------------
  ! Module parameters used in NetCDF interface
  !-------------------------------------------
  integer, parameter :: TY_I  = 1       ! type of field: integer
  integer, parameter :: TY_F  = 2       ! type of field: real (float)
  integer, parameter :: NDIM1 = 1       ! dimension of field: 1
  integer, parameter :: NDIM2 = 2       ! dimension of field: 2
!==============================================================================
contains ! subroutines
!=====================
  subroutine process_synop (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
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

    !----------------
    ! local variables
    !----------------
    integer                :: ni         ! actual number of coef. for interpol.
    integer                :: i, k       ! loop index
    integer                :: tsk        ! task
    integer(i8)            :: iatm       ! model column flag
    integer       ,pointer :: ty (:)     ! unpacked types
    logical                :: change     ! argument from shrink_report
    integer                :: iii, ioi   !+++ work around NEC SX bug
    integer                :: ifail      ! error code
    real(wp)               :: ff, dd     ! wind speed, direction
    real(wp)               :: t,  tv     ! temporaries (temperature)
    real(wp)               :: rh, gh     ! temporaries (humidity)
    real(wp)               :: ff_x, ff_y ! derivatives
    real(wp)               :: dd_x, dd_y ! derivatives
    integer                :: i_t,i_r,i_z! interpolation space indices
    integer                :: i_u,i_v,i_d! interpolation space indices
    integer                :: itv        ! interpolation space indices
    integer                :: o_t,o_r,o_z! observation   space indices
    integer                :: o_u,o_v,o_f! observation   space indices
    integer                :: o_td,o_ps  ! observation   space indices
    integer                :: o_d        ! observation   space indices
    real(wp)               :: e, td      ! temporaries for d td / d rh
    real(wp)               :: de_dt, dtd_de, dtd_dt, dtd_dr
    real(wp)               :: dt_tv      ! derivatives
    real(wp)               :: dt_gh      ! derivatives
    real(wp)               :: drh_tv     ! derivatives
    real(wp)               :: drh_gh     ! derivatives
    real(wp)    ,parameter :: rhm = 1.e-10_wp
    real(wp)               :: sf         ! scaling factor, wind speed profile
    real(wp)               :: lev        ! level (temporary)
    real(wp)               :: p          ! pressure (temporary)
    real(wp)               :: wmax       ! max. weight of all columns
    logical                :: lw         ! logarithmic/scaled wind profile?
    !---------------------------------------
    ! variables valid for a given time range
    !---------------------------------------
    type(t_col) ,pointer :: c                 ! temporary
    integer              :: ic                ! column index
    real(wp)             :: vmax_10m    (nt)  ! max. 10m wind speed
    real(wp)             :: tmin_2m     (nt)  ! min.  2m temperature
    real(wp)             :: tmax_2m     (nt)  ! max.  2m temperature
    real(wp)             :: tot_prec    (nt)  ! total precipitation
    real(wp)             :: aswdir_s    (nt)  ! downw.direct SW rad
    real(wp)             :: aswdifd_s   (nt)  ! down. diffusive r.
    !---------------------------------
    ! variables to overtake from model
    !---------------------------------
    real(wp)             :: clct              ! total  cloud cover
    real(wp)             :: clcl              ! low    cloud cover
    real(wp)             :: clcm              ! medium cloud cover
    real(wp)             :: clch              ! high   cloud cover
    real(wp)             :: ceil              ! cloud ceiling (at 1 g.pt.)
    real(wp)             :: vis               ! hor. visibility (at 1 g.pt.)
    real(wp)             :: hsurf_1gp         ! model orography (at 1 g.pt.)
    real(wp), allocatable:: clc_1gp  (:)      ! cloud cover on model levels
    real(wp), allocatable:: hml_1gp  (:)      ! height of model half levels
    !--------------------------------------
    ! some constants
    !   hard coded, values are not critical
    !--------------------------------------
!   real(wp) ,parameter    :: R     = 287.04_wp  ! [J /kg /K] gas constant
!   real(wp) ,parameter    :: gamma = gacc / c_p ! [K /m] adiabatic T gradient
!   real(wp) ,parameter    :: cp_R  = c_p  / R   ! c_p / R
    !---------------------------
    ! process optional arguments
    !---------------------------
    !  TSK_INIT       ! initialize modules
    !  TSK_READ       ! read observations
    !  TSK_SET_CHR    ! set observation characteristics
    !  TSK_SETUP_COLS ! setup columns
    !  TSK_SETUP_FUL0 ! setup interpolation space
    !  TSK_SETUP_FULL ! setup description of PSAS-space
    !  TSK_R          ! setup observational error
    !  TSK_Y          ! run forward operator
    !  TSK_K          ! evaluate linear operator
    !==================================================
    ! skip tasks not required for this observation type
    !==================================================
    tsk = task
    tsk = iand (task, not (&
          TSK_READ))        ! Input is read by separate BUFR reading routine
    if (tsk == 0) return
    !===============
    ! initialisation
    !===============
    if (iand (TSK_INIT,tsk) /= 0) then
      call read_synop_nml
      tsk = tsk - TSK_INIT
      if (tsk == 0) return
    endif

    !=============================================
    ! TSK_SET_CHR: set observation characteristics
    !=============================================
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      tsk = tsk - TSK_SET_CHR
      spot% int_type  = ITY_ICOL
      spot% cost      = 1._wp
      spot% nr        = spot% o%n
      spot% char      = CHR_NONL
      if (tsk == 0) return
    endif

    !==========================================
    ! tsk == TSK_SHRINK:
    ! release unused observations in the report
    !==========================================
    if (iand (TSK_SHRINK,tsk) /= 0) then
      call shrink_report (spot, obs%o, state, change)
      if (.not. change .or. spot% i%n==0) then
        tsk = tsk - TSK_SHRINK
        if (tsk == 0) return
      endif
    endif

    !===========
    ! PSAS-space
    !===========
    if (iand (TSK_SETUP_FUL0+TSK_SHRINK,tsk) /= 0) then
!     if (dace% pe==obs% o% pe) then
        ioi = spot% o% i
        i_t = 0; i_r = 0; i_u = 0; i_v = 0; i_z = 0; i_d = 0; itv = 0
        o_z = 0; o_ps = 0; o_t = 0
        do i = 1, spot% o% n
          select case (obs% o% varno (ioi+i))
          case (VN_U, VN_U10M, VN_V, VN_V10M, VN_FF, VN_DD)
            if (i_u == 0) then
              i_u = 1
            endif
            if (i_v == 0) then
              i_v = 1
            endif
            o_v   = i
          case (VN_RH, VN_RH2M)
            if (i_r == 0) then
              i_r = 1
              o_r = i
            endif
            ! Emulate old bug for stations reporting rh2m but not t2m
            if (i_t == 0 .and. version >= 3 .and. bugfix_rh2m) then
              i_t = 1
              o_t = i
            end if
          case (VN_TD, VN_TD2M)
            if (i_r == 0) then
              i_r = 1
              o_r = i
            endif
            if (i_t == 0) then
              i_t = 1
              o_t = i
            endif
          case (VN_T2M, VN_T)
            if (i_t == 0) then
              i_t = 1
              o_t = i
            endif
            ! Emulate old bug for stations reporting t2m but not td2m/rh2m
            if (i_r == 0 .and. version >= 3 .and. bugfix_t2m) then
              i_r = 1
              o_r = i
            end if
          case (VN_Z)
            if (i_z == 0) then
              i_z = 1
            endif
            o_z = i
          case (VN_PS, VN_P)
            if (i_z == 0) then
              i_z = 1
            endif
            o_ps = i
          case default
            cycle
          end select
        end do

        ni = 0
        if (i_z > 0) then
          ni = ni + 1; i_z = ni
        endif
        if (i_t > 0) then
          ni = ni + 1; i_t = ni
        endif
        if (i_r > 0) then
          ni = ni + 1; i_r = ni
        endif
        if (i_u > 0) then
          ni = ni + 1; i_u = ni
        endif
        if (i_v > 0) then
          ni = ni + 1; i_v = ni
        endif
        if (i_d > 0) then
          ni = ni + 1; i_d = ni
        endif

!       if (ni == 0) then
!         ni  = ni + 1
!         i_d = ni
!       endif

        if (iand (TSK_SHRINK,tsk) /= 0) then
           k = 0
           if (ni == 0) then
             spot% i% n = ni
           else
             do i = 1, spot% i% n
               select case (obs% o% t_int (spot%i%i+i))
               case (OBS_DUM)
                 k = i_d
               case (OBS_U)
                 k = i_u
               case (OBS_V)
                 k = i_v
               case (OBS_T)
                 k = i_t
               case (OBS_TV)
                 k = i_t
               case (OBS_RH)
                 k = i_r
               case (OBS_H)
                 k = i_z
               end select
               if (k > 0) then
                 obs% o% t_int (spot%i%i+k) = obs% o% t_int (spot%i%i+i)
                 obs% o% lev   (spot%i%i+k) = obs% o% lev   (spot%i%i+i)
               endif
             end do
             spot% i% n = ni
           endif
        else
          call new_int (obs% o, spot, ni)
          iii = spot% i% i
          ioi = spot% o% i

          if (i_d > 0) then
            obs% o% t_int (iii+i_d) = OBS_DUM
            obs% o% lev   (iii+i_d) = log (100000._wp)
          endif

          if (i_u > 0 .and. i_v > 0) then
            if (any (obs% o% body (ioi+o_v)% lev_typ == [VN_HEIGHT,VN_HOSAG])) then
              if (obs% o% body (ioi+o_v)% plev <= 0._sp) then
                lev = cols% col(spot% col% h% imc(1,1))% s% ps
              else
                lev = real (obs% o% body (ioi+o_v)% plev, wp)
              endif
            else
              lev = obs% o% olev (ioi+o_v)
            end if
            obs% o% t_int (iii+i_u) = OBS_U
            obs% o% t_int (iii+i_v) = OBS_V
            obs% o% lev   (iii+i_u) = log (lev)
            obs% o% lev   (iii+i_v) = log (lev)
          else if (i_u > 0 .or. i_v > 0) then
            call finish ('process_synop(TSK_SETUP_FUL0)','i_u>0 /= i_v>0')
          endif

          if (i_t > 0) then
            if (any (obs% o% body (ioi+o_t)% lev_typ == [VN_HEIGHT,VN_HOSAG])) then
              if (obs% o% body (ioi+o_t)% plev <= 0._sp) then
                lev = cols% col(spot% col% h% imc(1,1))% s% ps
              else
                lev = real (obs% o% body (ioi+o_t)% plev, wp)
              endif
            else
              lev = obs% o% olev (ioi+o_t)
            end if
            if (version < 3) then
             obs%o% t_int (iii+i_t) = OBS_T     ! old formulation
            else
             obs%o% t_int (iii+i_t) = OBS_TV
            end if
            obs% o% lev   (iii+i_t) = log (lev)
          endif
          if (i_r > 0) then
            if (any (obs% o% body (ioi+o_r)% lev_typ == [VN_HEIGHT,VN_HOSAG])) then
              if (obs% o% body (ioi+o_r)% plev <= 0._sp) then
                lev = cols% col(spot% col% h% imc(1,1))% s% ps
              else
                lev = real (obs% o% body (ioi+o_r)% plev, wp)
              endif
            else
              lev = obs% o% olev (ioi+o_r)
            end if
            obs% o% t_int (iii+i_r) = OBS_RH
            obs% o% lev   (iii+i_r) = log (lev)
          endif
          if (i_z > 0) then
            if (o_z > 0) then
              lev = obs% o% olev (ioi+o_z)
            else
              if (obs% o% body (ioi+o_ps)% o > 0._sp) then
                lev = real (obs% o% body (ioi+o_ps)% o, wp)
              else
                obs% o% body (ioi+o_ps)% o = rvind
                if (obs% o% body (ioi+o_ps)% plev <= 0._sp) then
                  lev = cols% col(spot% col% h% imc(1,1))% s% ps
                else
                  lev = real (obs% o% body (ioi+o_ps)% plev, wp)
                endif
              endif
            endif
            obs% o% t_int (iii+i_z) = OBS_H
            obs% o% lev   (iii+i_z) = log (lev)
          endif
        endif
!     endif
      if (iand (TSK_SETUP_FUL0,tsk) /= 0) tsk = tsk - TSK_SETUP_FUL0
      if (iand (TSK_SHRINK    ,tsk) /= 0) tsk = tsk - TSK_SHRINK
      if (tsk == 0) return
    endif

    !========================
    ! determine model columns
    !========================
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then

      select case (version)
      case default
         iatm = COL_TV + COL_Q  + COL_UV
      case (3:)
         iatm = COL_TV + COL_RH + COL_UV
      end select
      if (present_clc) iatm = iatm + COL_CLC + COL_GEO
      select case (method)
      case ('MEC', 'GMESTAT', 'VERI_ENS')
        iatm = iatm + COL_RANGE
      end select

      call idx_init (       &
             spot% col% c,  &! <-  column descriptor
             spot% col% h,  &!  -> interpolation coefficients
             obs% o% mc,    &! <-> model column descriptors
             iatm,          &! <-  fields required
             0,             &! <-  tracers required
             atm% grid,     &! <-  model grid
             spot% i_time,  &! <-  time slot
             spot% w_time   )! <-  time interpolation weight

      if (spot% col% h% imc(1,1) == 0) call decr_rpt_use (spot, CHK_DOMAIN)

      tsk = tsk - TSK_SETUP_COLS
      if (tsk == 0) return
    endif

    !=================================================================
    ! load temp data from observation data type for further processing
    !=================================================================
!   n = spot% o% n
!!! call load_synop (obs% o, spot, o)

    !===========================================
    ! setup description of PSAS-space
    ! observed values were set up while reading.
    !===========================================
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk = tsk - TSK_SETUP_FULL
      if (tsk == 0) return
    endif

    !===============================================
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    !===============================================
    if (iand (TSK_R,tsk) /= 0) then
!     ob => obs% o% obs   (spot% o% i + 1 : spot% o% i + spot% o% n)
      ty => obs% o% varno (spot% o% i + 1 : spot% o% i + spot% o% n)
!     lv => obs% o% lev   (spot% i% i + 1 : spot% o% i + spot% i% n)
      if (obs% o% pe == dace% pe) then
        !=================================
        ! setup SYNOP observational errors
        !=================================
        call synop_obs_err (obs, spot, ty)
        !=================
        ! setup vqc bounds
        !=================
        call set_vqc_insitu (spot, obs% o)
        !====================================================
        ! Check difference of station height to model surface
        !====================================================
        if (      dz_max_v10m       >  0._wp         .and.  &
                  spot% hd% obstype == OT_SYNOP      .and.  &
            (     spot% z           == -999._wp .or.        &
             abs (spot% z - spot% gp_bg/gacc) > dz_max_v10m)) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_U,VN_U10M,VN_V,VN_V10M,VN_FF,VN_DD)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
        if (      dz_max_rh2m       >  0._wp         .and.  &
                  spot% hd% obstype == OT_SYNOP      .and.  &
            (     spot% z           == -999._wp .or.        &
             abs (spot% z - spot% gp_bg/gacc) > dz_max_rh2m)) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_RH,VN_RH2M,VN_TD,VN_TD2M)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
        if (      dz_max_t2m        >  0._wp         .and. &
                  spot% hd% obstype == OT_SYNOP      .and. &
            (     spot% z           == -999._wp .or.       &
             abs (spot% z - spot% gp_bg/gacc) > dz_max_t2m)) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_T,VN_T2M)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
        !==============================================================
        ! Apply check on SSO standard deviation (ssd_bg must be known).
        ! Optionally adjust observation error (not yet implemented).
        !==============================================================
        if (ssdmax_v10m > 0._wp .and. spot% hd% obstype == OT_SYNOP   &
                                .and. spot% ssd_bg      >  ssdmax_v10m) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_U,VN_U10M,VN_V,VN_V10M,VN_FF,VN_DD)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
        if (ssdmax_rh2m > 0._wp .and. spot% hd% obstype == OT_SYNOP   &
                                .and. spot% ssd_bg      >  ssdmax_rh2m) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_RH,VN_RH2M,VN_TD,VN_TD2M)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
        if (ssdmax_t2m  > 0._wp .and. spot% hd% obstype == OT_SYNOP  &
                                .and. spot% ssd_bg      >  ssdmax_t2m) then
          ioi = spot% o% i
          do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_T,VN_T2M)
              call decr_use (obs% o% body(ioi+i)% use, STAT_PASSIVE, check=CHK_HEIGHT)
            end select
          end do
        end if
      endif
      tsk = tsk - TSK_R
      if (tsk == 0) return
    endif

    !==============================
    ! tsk == TSK_Y
    ! evaluate observation operator
    !==============================
    if (iand (TSK_Y,tsk) /= 0) then
      iii = spot% i% i
      ioi = spot% o% i
      !------------
      ! set indices
      !------------
      i_t = 0; i_r = 0; i_u = 0; i_v = 0; i_z = 0; itv = 0
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii+i))
        case (OBS_H)
          i_z = i
        case (OBS_RH)
          i_r = i
        case (OBS_U)
          i_u = i
        case (OBS_V)
          i_v = i
        case (OBS_T)
          i_t = i
        case (OBS_TV)
          itv = i
        case default
          call finish ('process_synop(TSK_Y)','invalid t_int')
        end select
      end do

      !-----------------------------------------
      ! interpolate time range variables for MEC
      !-----------------------------------------
      select case (method)
      case ('MEC', 'GMESTAT', 'VERI_ENS')

        allocate ( clc_1gp (cols% ke) )
        allocate ( hml_1gp (cols% ke) )
        vmax_10m  = 0._wp
        tmin_2m   = 0._wp
        tmax_2m   = 0._wp
        tot_prec  = 0._wp
        aswdir_s  = 0._wp
        aswdifd_s = 0._wp
        clct      = 0._wp
        clcl      = 0._wp
        clcm      = 0._wp
        clch      = 0._wp
        ceil      = 0._wp
        vis       = 0._wp
        hsurf_1gp = 0._wp
        clc_1gp   = 0._wp
        hml_1gp   = 0._wp
        wmax      = 0._wp

        do i = 1, size(spot% col% h% imc,1)

          ic = spot% col% h% imc(i,1)
          if (ic==0) exit
          c => cols% col(ic)

          vmax_10m  = vmax_10m  + spot% col% h% w(i) * c% s% vmax_10m
          tmin_2m   = tmin_2m   + spot% col% h% w(i) * c% s% tmin_2m
          tmax_2m   = tmax_2m   + spot% col% h% w(i) * c% s% tmax_2m
          tot_prec  = tot_prec  + spot% col% h% w(i) * c% s% tot_prec
          aswdir_s  = aswdir_s  + spot% col% h% w(i) * c% s% aswdir_s
          aswdifd_s = aswdifd_s + spot% col% h% w(i) * c% s% aswdifd_s
          clct      = clct      + spot% col% h% w(i) * c% s% clct
          clcl      = clcl      + spot% col% h% w(i) * c% s% clcl
          clcm      = clcm      + spot% col% h% w(i) * c% s% clcm
          clch      = clch      + spot% col% h% w(i) * c% s% clch
          !   note: horizontal averaging of ceiling is highly questionable!
!         ceil      = ceil      + spot% col% h% w(i) * c% s% ceiling
!         hsurf     = hsurf     + spot% col% h% w(i) * c% s% geosp / gacc
          !   --> use ceiling of grid pt with largest interpolation weight
          !      (alternative: min of ceiling AGL over all g.p. with weight > 0)
          if (spot% col% h% w(i) > wmax) then
            wmax      = spot% col% h% w(i)
            ceil      = c% s% ceiling
            vis       = c% s% vis
            hsurf_1gp = c% s% geosp / gacc
            if ((associated(c% clc)) .and. (associated(c% geo))) then
              clc_1gp = c%    clc
              hml_1gp = c%    geo   / gacc
            endif
          endif
        end do
        if (.not. ((associated(c% clc)) .and. (associated(c% geo))))           &
          clc_1gp = -99._wp

      end select

      !---------------------------------------------------------
      ! scaling factor for near-surface wind observation (DRIBU)
      !---------------------------------------------------------
      sf =  1._wp
      lw = (spot% hd% obstype == OT_DRIBU) .and.  &
           (v_buoy > 0) .and. (spot% z0_bg > 0._wp)
      if (lw) then
         do i = 1, spot% o% n
            select case (obs% o% varno (ioi+i))
            case (VN_FF,VN_U,VN_V)
               sf = get_scalf (spot, obs% o% body(ioi+i), obs% o% olev(ioi+i))
               exit
            end select
         end do
      end if
      !------------------------------
      ! evaluate observation operator
      !------------------------------
      obs% o% body (ioi+1:ioi+spot% o% n)% op_na = 0
      ff = -1._wp
      rh = -1._wp; if (i_r > 0) rh = xi% x(iii+i_r)
      t  = -1._wp; if (i_t > 0) t  = xi% x(iii+i_t)
      tv = -1._wp; if (itv > 0) tv = xi% x(iii+itv)
      !-------------------------------------------------------------------
      ! Allow verification mode to proceed (analyses need not provide t2m)
      !-------------------------------------------------------------------
      if (version >= 3 .and. tv < 0._wp) itv = 0

      if (itv > 0) then
         p  = exp (obs% o% lev(iii+itv))
         gh = rh
         call trh_tvgh (t,                            &!  -> t
                        rh,                           &!  -> rh
                        tv,                           &! <-  tv
                        gh,                           &! <-  gh
                        p,                            &! <-  p
                        ifail                       )  !  -> error code
         if (ifail < 0) then
            write(0,*) spot% statid, ": process_synop(TSK_Y): t,rh,tv,gh,p =", &
                 real([t,rh,tv,gh,p])
            call finish ('process_synop(TSK_Y)','trh_tvgh failed')
         end if
      end if
      do i = 1, spot% o% n
        select case (obs% o% varno (ioi+i))
        case (VN_FF)
          if (ff < 0._wp) call fd_uv (ff, dd, xi%x (iii+i_u), xi%x (iii+i_v))
          y%x (ioi+i) = ff * sf
        case (VN_DD)
          if (ff < 0._wp) call fd_uv (ff, dd, xi%x (iii+i_u), xi%x (iii+i_v))
          y%x (ioi+i) = dd
        case (VN_T2M, VN_T)
          y%x (ioi+i) = t
          if (0._wp >= t) then
            obs% o% body(ioi+i)% op_na = 1
            y% x (ioi+i) = rvind
          endif
        case (VN_RH, VN_RH2M)
          y%x (ioi+i) = rh
!         if (0._wp > rh) y% x (ioi+i) = rvind
        case (VN_TD, VN_TD2M)
          td = -1._wp
          if (t > 0._wp) td = t_esw (min (1._wp, max (rhm, rh)) * esw_t (t))
          y%x (ioi+i) = td
          if (0._wp >= td) y% x (ioi+i) = rvind
        case (VN_U,  VN_U10M)
          y%x (ioi+i) =  xi%x (iii+i_u) * sf
        case (VN_V,  VN_V10M)
          y%x (ioi+i) =  xi%x (iii+i_v) * sf
        case (VN_Z)
          y%x (ioi+i) =  xi%x (iii+i_z)
        case (VN_PS, VN_P)
          y%x (ioi+i) = (xi%x (iii+i_z) - obs%o%olev(ioi+i)) * spot% pz_bg &
                                        + obs%o%body(ioi+i)% o             &
                                        - obs%o%body(ioi+i)% bc
        case (VN_TSEA)
          y%x (ioi+i) = spot% tw_bg
          if (spot% tw_bg <= 0._wp) then
            obs% o% body(ioi+i)% op_na = 1
            y%x (ioi+i) = rvind
          end if
        case (VN_SDEPTH)
          !--------------------------------------------
          ! passive monitoring, skipping glacier points
          !--------------------------------------------
          y%x (ioi+i) = spot% hs_bg
          if (spot% hs_bg < 0._wp .or. spot% hs_bg > 39.999_wp) then
            obs% o% body(ioi+i)% op_na = 1
            y%x (ioi+i) = rvind
          end if
        case (VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW)
          !------------------------------------------------------
          ! time range variables, handled by routine apply_trange
          !------------------------------------------------------
          y%x (ioi+i) = rvind
        case (VN_N, VN_NH, VN_CEIL, VN_N_L, VN_N_M, VN_N_H, VN_VV )
          !------------------------------------------------------
          ! cloud related variables, handled by routine apply_clc
          !------------------------------------------------------
          y%x (ioi+i) = rvind
        case (VN_WW, VN_GCLG, VN_ICLG )
          !---------
          ! obs only
          !---------
          y%x (ioi+i) = rvind
        case default
          call finish ('process_synop(TSK_Y)',                        &
                       'invalid varno: '//char3(obs% o% varno (ioi+i)))
        end select
      end do
      !---------------------------------
      ! set time range variables for MEC
      !---------------------------------
      select case (method)
      case ('MEC', 'GMESTAT', 'VERI_ENS')
        call apply_trange (spot, obs, y% x(spot% o% i + 1         :  &
                                           spot% o% i + spot% o% n), &
                           vmax_10m, tmin_2m, tmax_2m,               &
                           tot_prec, aswdir_s, aswdifd_s             )

        call apply_clc (spot, obs, y, clct, clcl, clcm, clch,        &
                        ceil, vis, hsurf_1gp, clc_1gp, hml_1gp )

        deallocate ( clc_1gp )
        deallocate ( hml_1gp )
      end select

      tsk = tsk - TSK_Y
      if (tsk == 0) return
    endif

    !================
    ! tsk == TSK_K
    ! derive Jakobian
    !================
    if (iand (TSK_K,tsk) /= 0) then
      ioi  = spot% o% i
      iii  = spot% i% i
      !------------
      ! set indices
      !------------
      i_t = 0; i_r = 0; i_u = 0; i_v = 0; i_z = 0; itv = 0
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii +i))
        case (OBS_H)
          i_z = i
        case (OBS_RH)
          i_r = i
        case (OBS_U)
          i_u = i
        case (OBS_V)
          i_v = i
        case (OBS_T)
          i_t = i
        case (OBS_TV)
          itv = i
        case default
          write (0,*)  'process_synop(TSK_K): t_int=', obs% o% t_int(iii+i)
          call finish ('process_synop(TSK_K)','invalid t_int')
        end select
      end do

      o_t=0; o_r=0; o_u=0; o_v=0; o_z=0; o_f=0; o_d=0; o_td=0; o_ps=0
      do i = 1, spot% o%n
        select case (obs% o% varno (ioi+i))
        case (VN_U,  VN_U10M)
          o_u  = i
        case (VN_V,  VN_V10M)
          o_v  = i
        case (VN_RH, VN_RH2M)
          o_r  = i
        case (VN_T2M, VN_T)
          o_t  = i
        case (VN_TD, VN_TD2M)
          o_td = i
        case (VN_Z)
          o_z  = i
        case (VN_PS, VN_P)
          o_ps = i
        case (VN_FF)
          o_f  = i
        case (VN_DD)
          o_d  = i
        end select
      end do

      ff = -1._wp
      rh = -1._wp; if (i_r > 0) rh = xi% x(iii+i_r)
      t  = -1._wp; if (i_t > 0) t  = xi% x(iii+i_t)
      tv = -1._wp; if (itv > 0) tv = xi% x(iii+itv)
      if (itv > 0) then
         !----------------------------------------
         ! Adjoint code requires t2m and td2m/rh2m
         !----------------------------------------
         if (version >= 3 .and. tv < 0._wp) then
            call finish ('process_synop(TSK_K)','tv not available')
         end if

         p  = exp (obs% o% lev(iii+itv))
         gh = rh
         call trh_tvgh (t,                            &!  -> t
                        rh,                           &!  -> rh
                        tv,                           &! <-  tv
                        gh,                           &! <-  gh
                        p,                            &! <-  p
                        ifail,                        &!  -> error code
                        dt_tv, dt_gh, drh_tv, drh_gh)  !  -> gradient
         if (ifail < 0) then
            write(0,*) spot% statid, ": process_synop(TSK_K): t,rh,tv,gh,p =", &
                 real([t,rh,tv,gh,p])
            call finish ('process_synop(TSK_K)','trh_tvgh failed')
         end if
      end if
      !------------------------------------
      ! derivatives d ff / d u , d ff / d v
      !------------------------------------
      if (o_f > 0 .or. o_d > 0) then
        call fd_uv (ff, dd, obs% xi% x(iii+i_u), obs% xi% x(iii+i_v), &
                    ff_x, ff_y, dd_x, dd_y                            )
      endif
      !---------------------------------------------------------
      ! scaling factor for near-surface wind observation (DRIBU)
      !---------------------------------------------------------
      sf = 1._wp
      lw = (spot% hd% obstype == OT_DRIBU) .and.  &
           (v_buoy > 0) .and. (spot% z0_bg > 0._wp)
      i  = max (o_f, o_u)
      if (lw .and. i > 0)                                                &
           sf = get_scalf (spot, obs% o% body(ioi+i), obs% o% olev(ioi+i))
      !-------------------------------------
      ! derivatives d td / t t , d td / d rh
      !-------------------------------------
      if (o_td > 0) then
        if (rh >= 1._wp) then
          dtd_dt = 1._wp
          dtd_dr = 0._wp
          td     = t
        else if (rh < rhm) then
          dtd_dt = 0._wp
          dtd_dr = 0._wp
          td     = t_esw (rhm * esw_t (t))
        else
          e      =          esw_t   (t)
          de_dt  =          desw_dt (t)
          td     =          t_esw   (rh * e)
          dtd_de = 1._wp  / desw_dt (td)
          dtd_dt = dtd_de * de_dt *  rh
          dtd_dr = dtd_de * e
        endif
        if (0._wp >= td) then
          td     = rvind
          dtd_dt = 0._wp
          dtd_dr = 0._wp
        endif
      endif
      !-------------
      ! set Jakobian
      !-------------
      k = obs% H% ia (iii + 1)
!NEC$ nomove
      do i = 1, spot% i%n
        select case (obs% o% t_int (iii+i))
        case (OBS_H)
          if (o_z > 0) then
            obs% H% ja     (k)     = ioi+o_z        ! row    index
            obs% H% packed (k)     = 1._mp          ! coefficient
            k = k + 1
          endif
          if (o_ps> 0) then
            obs% H% ja     (k)     = ioi+o_ps       ! row    index
            obs% H% packed (k)     = spot% pz_bg    ! coefficient
            k = k + 1
          endif
        case (OBS_RH)
         if (itv > 0) then
          !-------------------------
          ! Control variables: tv,gh
          !-------------------------
          if (o_r > 0) then
            obs% H% ja     (k)     = ioi+o_r        ! row    index
            obs% H% packed (k)     = drh_gh         ! coefficient
            k = k + 1
          endif
          if (o_t > 0) then
            obs% H% ja     (k)     = ioi+o_t        ! row    index
            obs% H% packed (k)     = dt_gh          ! coefficient
            k = k + 1
          endif
          if (o_td > 0) then
            obs% H% ja     (k)     = ioi+o_td       ! row    index
            obs% H% packed (k)     = dtd_dr * drh_gh + dtd_dt * dt_gh
            k = k + 1
          endif
         else
          if (o_r > 0) then
            obs% H% ja     (k)     = ioi+o_r        ! row    index
            obs% H% packed (k)     = 1._mp          ! coefficient
            k = k + 1
          endif
          if (o_td > 0) then
            obs% H% ja     (k)     = ioi+o_td       ! row    index
            obs% H% packed (k)     = dtd_dr         ! coefficient
            k = k + 1
          endif
         end if
        case (OBS_U)
          if (o_u > 0) then
            obs% H% ja     (k)     = ioi+o_u        ! row    index
            obs% H% packed (k)     = real (sf,mp)   ! coefficient
            k = k + 1
          endif
          if (o_f > 0) then
            obs% H% ja     (k)     = ioi+o_f        ! row    index
            obs% H% packed (k)     = ff_x * sf      ! coefficient
            k = k + 1
          endif
          if (o_d > 0) then
            obs% H% ja     (k)     = ioi+o_d        ! row    index
            obs% H% packed (k)     = dd_x           ! coefficient
            k = k + 1
          endif
        case (OBS_V)
          if (o_v > 0) then
            obs% H% ja     (k)     = ioi+o_v        ! row    index
            obs% H% packed (k)     = real (sf,mp)   ! coefficient
            k = k + 1
          endif
          if (o_f > 0) then
            obs% H% ja     (k)     = ioi+o_f        ! row    index
            obs% H% packed (k)     = ff_y * sf      ! coefficient
            k = k + 1
          endif
          if (o_d > 0) then
            obs% H% ja     (k)     = ioi+o_d        ! row    index
            obs% H% packed (k)     = dd_y           ! coefficient
            k = k + 1
          endif
        case (OBS_T)
          if (o_t > 0) then
            obs% H% ja     (k)     = ioi+o_t        ! row    index
            obs% H% packed (k)     = 1._mp          ! coefficient
            k = k + 1
          endif
          if (o_td > 0) then
            obs% H% ja     (k)     = ioi+o_td       ! row    index
            obs% H% packed (k)     = dtd_dt         ! coefficient
            k = k + 1
          endif
        case (OBS_TV)
          if (o_t > 0) then
            obs% H% ja     (k)     = ioi+o_t        ! row    index
            obs% H% packed (k)     = dt_tv          ! coefficient
            k = k + 1
          endif
          if (o_r > 0) then
            obs% H% ja     (k)     = ioi+o_r        ! row    index
            obs% H% packed (k)     = drh_tv         ! coefficient
            k = k + 1
          endif
          if (o_td > 0) then
            obs% H% ja     (k)     = ioi+o_td       ! row    index
            obs% H% packed (k)     = dtd_dt * dt_tv + dtd_dr * drh_tv
            k = k + 1
          endif
        case default
          write (0,*)  'process_synop(TSK_K): t_int=', obs% o% t_int(iii+i)
          call finish ('process_synop(TSK_K)','invalid t_int')
        end select
        obs% H% ia (iii + i + 1) = k
      end do

      !------------------------------
      ! evaluate observation operator
      !------------------------------
!NEC$ nomove
      do i = 1, spot% o% n
        select case (obs% o% varno (ioi+i))
        case (VN_FF)
          obs%yi%x (ioi+i) = ff * sf
        case (VN_DD)
          obs%yi%x (ioi+i) = dd
        case (VN_TD,  VN_TD2M)
          obs%yi%x (ioi+i) = td
        case (VN_T2M, VN_T)
          obs%yi%x (ioi+i) = t
          if (0._wp >= t) then
            obs% o% body(ioi+i)% op_na = 1
            obs%yi%x (ioi+i) = rvind
          endif
        case (VN_RH, VN_RH2M)
          obs%yi%x (ioi+i) = rh
        case (VN_U,  VN_U10M)
          obs%yi%x (ioi+i) = obs%xi%x (iii+i_u) * sf
        case (VN_V,  VN_V10M)
          obs%yi%x (ioi+i) = obs%xi%x (iii+i_v) * sf
        case (VN_Z)
          obs%yi%x (ioi+i) = obs%xi%x (iii+i_z)
        case (VN_PS, VN_P)
          obs%yi%x (ioi+i) = (obs%xi%x (iii+i_z) - obs%o%olev(ioi+i)) * spot% pz_bg &
                                                 + obs%o%body(ioi+i)% o             &
                                                 - obs%o%body(ioi+i)% bc
        case (VN_TSEA)
          obs%yi%x (ioi+i) = spot% tw_bg
          if (spot% tw_bg <= 0._wp) then
            obs% o% body(ioi+i)% op_na = 1
            obs%yi%x (ioi+i) = rvind
          end if
        case (VN_SDEPTH)
          obs%yi%x (ioi+i) = spot% hs_bg
          if (spot% hs_bg < 0._wp .or. spot% hs_bg > 39.999_wp) then
            obs% o% body(ioi+i)% op_na = 1
            obs%yi%x (ioi+i) = rvind
          end if
        case (VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW)
          !------------------------------------------------------
          ! time range variables, handled by module mo_obs_trange
          !------------------------------------------------------
          obs%yi%x (ioi+i) = rvind
        case (VN_N, VN_NH, VN_CEIL, VN_N_L, VN_N_M, VN_N_H, VN_VV)
          !---------------------------------------------
          ! passive monitoring, implementation postponed
          !---------------------------------------------
          obs%yi%x (ioi+i) = rvind
        case (VN_WW, VN_GCLG, VN_ICLG )
          !---------
          ! obs only
          !---------
          obs%yi%x (ioi+i) = rvind
        case default
          call finish ('process_synop(TSK_K)',                        &
                       'invalid varno: '//char3(obs% o% varno (ioi+i)))
        end select
      end do

      tsk = tsk - TSK_K
      if (tsk == 0) return
    endif

    !===========
    ! left tasks
    !===========
    if (tsk /= 0) then
      if (dace% lpio) write (6,*) 'process_synop: unknown task',tsk
      call finish ('process_synop','unknown task')
    endif
    !--------------------------------------------------------------------------
  contains
    !--------------------------------------------------------------------------
    pure function get_scalf (spot, body, lev) result (scalf)
      !--------------------------------------------------
      ! Derive scaling factor of logarithmic wind profile
      !--------------------------------------------------
      type(t_spot),  intent(in) :: spot   ! Observation metadata
      type(t_datum), intent(in) :: body   ! Observation data
      real(wp),      intent(in) :: lev    ! Level of observation
      real(wp)                  :: scalf  ! Scaling factor

      real(wp)                  :: h_ws   ! Height of wind sensor

      scalf = 1._wp

      select case (spot% hd% obstype)
      case (OT_DRIBU)
         if (v_buoy == 0) return
         if (body% lev_typ == VN_HEIGHT) then
            h_ws = lev
         else
            h_ws = buoy_h0_v
         end if
      case default
         !h_ws = 10._wp
         return
      end select

      if (spot% z0_bg > 0._wp .and. spot% z0_bg < h_ws) then
         scalf = log (h_ws / spot% z0_bg) / log (10._wp / spot% z0_bg)
      end if

    end function get_scalf
    !--------------------------------------------------------------------------
  end subroutine process_synop

!------------------------------------------------------------------------------

  subroutine apply_clc (spot, obs, y, clct, clcl, clcm, clch,                 &
                        ceil, vis, hsurf_1gp, clc_1gp, hml_1gp )
    !--------------------------------------------------------------------
    ! apply obs operator to compute model equivalents for cloud variables
    !--------------------------------------------------------------------

    type(t_spot)        ,intent(in)    :: spot     ! report header meta data
    type(t_obs_block)   ,intent(in)    :: obs      ! observation data
    type(t_vector_segm) ,intent(inout) :: y        ! model equivalent
    real(wp)            ,intent(in)    :: clct     ! total  cloud cover
    real(wp)            ,intent(in)    :: clcl     ! low    cloud cover
    real(wp)            ,intent(in)    :: clcm     ! medium cloud cover
    real(wp)            ,intent(in)    :: clch     ! high   cloud cover
    real(wp)            ,intent(in)    :: ceil     ! cloud ceiling ASL at 1 g.pt
    real(wp)            ,intent(in)    :: vis      ! horiz. visibility at 1 g.pt
    real(wp)            ,intent(in)    :: hsurf_1gp! model orography   at 1 g.pt
    real(wp)            ,intent(in)    :: clc_1gp (:) ! cloud cover on model lev
    real(wp)            ,intent(in)    :: hml_1gp (:) ! height of model levels

    real(wp)            ,allocatable   :: hhl (:) ! height of model half levels

      integer  :: i, j    ! observation indices
      integer  :: k       ! model level index
      integer  :: ke      ! number of model (main) levels
      real(wp) :: cil     ! cloud ceiling     (derived from model cloud column)
      real(wp) :: cbh     ! cloud base height (derived from model cloud column)
      logical  :: lcol(2) ! need to compute values from model cloud column

    ! note: missing value = -99. for clc?, ceil (-> subr. atm2col, mo_t_col.f90)
    !                            and clc (-> subr. process_synop)

      lcol = .false.

      !   assign values to 'y', based on input from 2-D model fields
      do i = 1, spot% o% n
        j = spot% o% i + i
        select case (obs% o% varno (j))
        case (VN_N)
          if (clct  >= 0._wp) y% x (j) = octa_percent (clct)
        case (VN_N_L)
          if (clcl  >= 0._wp) y% x (j) = octa_percent (clcl)
        case (VN_N_M)
          if (clcm  >= 0._wp) y% x (j) = octa_percent (clcm)
        case (VN_N_H)
          if (clch  >= 0._wp) y% x (j) = octa_percent (clch)
        case (VN_CEIL)
          ! 'ceil'  is defined ASL (above  sea   level) in ICON, thresh: > 50 %
          ! the obs is defined AGL (above ground level) in Synop, thr.: > 4 oct.
          if (ceil  >= 0._wp) y% x (j) = ceil - hsurf_1gp
          if (ceil  <  0._wp) lcol(1) = .true.
        case (VN_VV)
          if (vis   >= 0._wp) y% x (j) = vis
        case (VN_NH)
          lcol(2) = .true.
        end select
      end do

      !   if possible + required, compute CBH + ceiling from model cloud column
      if ((minval(clc_1gp) >= 0._wp) .and. ((lcol(1)) .or. (lcol(2)))) then
        !   compute height of model half levels
        ke  = size( clc_1gp )
        allocate ( hhl (ke) )
        hhl (ke+1) = hsurf_1gp
        do k = ke , 1 , -1
          hhl (k)  = hml_1gp(k) + (hml_1gp(k) - hhl(k+1))
        enddo
        !   compute cloud base height and ceiling from model cloud column
        cil = real( CBH_CLR , wp )
        cbh = real( CBH_CLR , wp )
        do k = 1, ke
          !   note: in ICON, a threshold of 0.5 is used for field 'ceil'
          if (clc_1gp(k) > 0.55_wp)  cil = min( cil , hhl(k+1) - hsurf_1gp )
          if (clc_1gp(k) > 0.01_wp)  cbh = min( cbh , hhl(k+1) - hsurf_1gp )
        enddo
        deallocate ( hhl )
        !   assign values of ceiling (if needed) and CBH to 'y'
        do i = 1, spot% o% n
          j = spot% o% i + i
          select case (obs% o% varno (j))
          case (VN_CEIL)
            if (lcol(1))  y% x (j) = cil
          case (VN_NH)
                          y% x (j) = cbh
          end select
        end do
      endif

  end subroutine apply_clc

!==============================================================================

  subroutine read_synop_bufr (bufr, spt, obs, lkeep, cc)

  type (t_bufr) ,intent(inout)        :: bufr  ! BUFR record to decode
  type (t_spot) ,intent(inout)        :: spt   ! meta information to set
  type (t_obs)  ,intent(inout)        :: obs   ! observations data type to set
  logical       ,intent(out)          :: lkeep ! accept observation ?
  integer       ,intent(in) ,optional :: cc    ! part of year ccyy

    !----------------
    ! index variables
    !----------------
    integer        :: is               ! sub-set    index
    integer        :: ie               ! entry in sub-set
    integer        :: id               ! descriptor index
    integer        :: i

    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
    integer            :: ival  ! value decoded from BUFR (integer)
    real(sp)           :: rval  ! value decoded from BUFR (real)
    character(len=8)   :: ymnem ! value decoded from BUFR (string)
!   integer            :: itype !
    logical            :: diff, m12, m21
    type (t_synop)     :: s1

    type(t_synop)      :: s                ! SYNOP observation read
    integer            :: yyyy,mo,dd,hh,mi ! actual time read
    character(len=32)  :: com
    logical            :: unused = .false.
!   integer            :: ip
    real(wp)           :: dt1, dt2         ! difference to analysis time

    !---------------------
    ! loop over data, copy
    !---------------------
    lkeep = .false.
    if (bufr% sec3% num_subsets /= 1) then
      call decr_rpt_use (spt, CHK_NOIMPL, &
                         comment='read_synop_bufr: cannot handle subsets')
      return
    endif
    if (bufr% sec0% edition == 4) then
      spt% corme = bufr% sec1% update      ! Correction message?
    end if
    is = 1
    call construct_synop (s)
    do ie=1,bufr% nbufdat(is)
      !-------------
      ! decode datum
      !-------------
      ival  = bufr% ibufdat (is,ie)
      !----------------------------------------
      ! if valid, copy datum from BUFR record
      !   to structures of type t_spot, t_datum
      !----------------------------------------
      if(ival /= inv_bufr) then
        !---------------------------------------
        ! decode real/character datum, mnemonics
        !---------------------------------------
        id    = bufr% idescidx(is,ie)
!       itype = bufr% itype(id)
        ymnem = bufr% ymnem(id)
        rval  = rvind
        if (bufr% is_char(id) == 0) then
          IF (ival /= inv_bufr) &
            rval = ival * bufr% scale(id)
        endif
        !----------------------
        ! copy individual datum
        !----------------------
        select case (ymnem)     ! BUFR data descriptors

        case ('MCORME') ! COR MESSAGE MARK
          spt% corme = ival
        !------------------------
        ! class 01 identification
        !------------------------
        case ('YDDDD') ! SHIP OR MOBILE LAND STATION IDENTIFIER
          call bufr_get_character (bufr, ival, spt% statid)
        case ('YCCCC') ! ICAO LOCATION INDICATOR
          call bufr_get_character (bufr, ival, spt% statid)
        case ('MABNN') ! BUOY/PLATFORM IDENTIFIER
          if(spt% hd% obstype == OT_DRIBU) then
            write(spt%  statid,'(i5.5)')       ival
          else
            write(spt%  statid,'("BP_",i5.5)') ival
          endif
        case ('YSSOSN') ! SHORT STATION OR SITE NAME
          if (spt% statid=='') call bufr_get_character (bufr,ival,spt% statid)
        case ('MDS') ! direction of motion
        case ('NVS') ! speed of motion [m/s]
        case ('MII') ! WMO   block number
          spt%  ident  = 1000      * int(rval)
        case ('NIII') ! WMO station number
          spt%  ident  = spt% ident + int(rval)
        !-------------------------
        ! class 02 instrumentation
        !-------------------------
        case ('NIX') ! type of station
        !-------------------------
        ! class 04 location (time)
        !-------------------------
        case ('MJJJ') ! year
          yyyy = ival
          if(yyyy<100 .and. present(cc)) yyyy = yyyy + cc * 100
        case ('MMM') ! month
          mo   = ival
        case ('MYYQ')   ! Q-bits for following value              [code table]
        case ('MYY') ! day
          dd   = ival
        case ('MGG') ! hour
          hh   = ival
        case ('NGG') ! minute
          mi   = ival
        !---------------------------------
        ! class 05 location (horizontal-1)
        !---------------------------------
        case ('MLAH')  ! latitude (high accuracy)                [deg]
          spt% col% c% dlat = rval
        case ('MLALAQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MLALA') ! latitude (coarse accuracy)              [deg]
          if (spt% col% c% dlat == rvind) spt% col% c% dlat = rval
        !---------------------------------
        ! class 06 location (horizontal-2)
        !---------------------------------
        case ('MLOH')  ! longitude (high accuracy)                [deg]
          spt% col% c% dlon = rval
        case ('MLOLOQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MLOLO') ! longitude (coarse accuracy)              [deg]
          if (spt% col% c% dlon == rvind) spt% col% c% dlon =rval
        !-----------------------------
        ! class 07 location (vertical)
        !-----------------------------
        case ('MHP') ! height of station                    [m]
          call set_datum (s% zr ,rval ,spt% corme)
        case ('MPN') ! pressure (vert.location)             [Pa]
          call set_datum (s% p ,rval ,spt% corme)
        case ('MPNQ')! Q-bits for following value           [code table]
!!!       call set_qbits (s% p ,ival)
        !---------------------------------
        ! class 08 significance qualifiers
        !---------------------------------
        !----------------------------------------
        ! class 10 vertical elements and pressure
        !----------------------------------------
        case ('MHHA')  ! height of land surface                  [m]
          call set_datum (s% z_ls ,rval ,spt% corme)
        case ('MPPP')  ! pressure                                [Pa]
          call set_datum (s% ps ,rval ,spt% corme)
        case ('MPPPQ') ! Q-bits for following value              [code table]
          call set_qbits (s% ps ,ival)
        case ('MPPPP') ! pressure reduced to mean sea level      [Pa]
          call set_datum (s% p_msl ,rval ,spt% corme)
        case ('MPPPPQ')! Q-bits for following value              [code table]
           call set_qbits (s% p_msl ,ival)
        case ('MPHPH') ! altimeter setting (QNH) (METAR)         [Pa]
           call set_datum (s% p_msl ,rval ,spt% corme)
        case ('NHNHN') ! geopotential                            [m**2/s**2]
          call set_datum (s% gp    ,rval ,spt% corme)
          s% gp% o = s% gp% o / gacc                        ! -> [gpm]
        case ('NHNHNQ')
          call set_qbits (s% gp    ,ival)
        !-----------------------------
        ! class 11 wind and turbulence
        !-----------------------------
        case ('NDD') ! wind direction        at 10 m   [deg]!->[rad]
          call set_datum (s% dd ,rval ,spt% corme) ! * d2r
        case ('NDDQ')
          call set_qbits (s% dd ,ival)
        case ('NFF') ! wind speed            at 10 m        [m/s]
          call set_datum (s% ff ,rval ,spt% corme)
        case ('NFFQ')
          call set_qbits (s% ff ,ival)
        !---------------------
        ! class 12 temperature
        !---------------------
        case ('MTTT')  ! dry bulb  temperature at  2 m        [K]
          call set_datum (s% t ,rval ,spt% corme)
        case ('MTTTQ')
          call set_qbits (s% t ,ival)
        case ('MTFTF') ! wet bulb  temperature at  2 m        [K]
        case ('MTDTD') ! dew point temperature at  2 m        [K]
          call set_datum (s% td ,rval ,spt% corme)
        case ('MTDTDQ')
          call set_qbits (s% td ,ival)
        !--------------------------------------------
        ! class 13 hydrographic/hydrological elements
        !--------------------------------------------
        case ('MUUU')  ! relative humidity     at  2 m   [%]->[ ]
          call set_datum (s% rh ,rval * 0.01_sp ,spt% corme)
        case ('MUUUQ')
          call set_qbits (s% rh ,ival)
        !----------------
        ! BUFR4 templates
        !----------------
        case ('MHOSNN') ! height of station ground above mean sea  [M]
          call set_datum (s% z_ls ,rval ,spt% corme)
        case ('MHOBNN') ! height of barometer above mean sea level [M]
          call set_datum (s% zr ,rval ,spt% corme)
        case ('MTDBTQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% t ,ival)
        case ('MTDBT')  ! temperature/dry bulb temperature         [K]
          call set_datum (s% t ,rval ,spt% corme)
        case ('NDNDNQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% dd ,ival)
        case ('NDNDN')  ! wind direction                           [DEGREE_TRUE]
          call set_datum (s% dd ,rval ,spt% corme)
        case ('NFNFNQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% ff ,ival)
        case ('NFNFN')  ! wind speed                               [M/S]
          call set_datum (s% ff ,rval ,spt% corme)
        case ('MTDNHQ') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
!!!       call set_qbits (s% td ,ival)
        case ('MTDNH')  ! dew-point temperature                    [K]
          call set_datum (s% td ,rval ,spt% corme)
        case ('NHHHN')  ! geopotential height                      [GPM]
          call set_datum (s% gp    ,rval ,spt% corme)
#ifdef  CHECKCODES
        !-----------------------
        ! check for unused codes
        !-----------------------
        case ('MOBITQ') ! overall quality bits                    [code table]
        case ('MADDF')  ! associated field significance           [code table]
        case ('MA1')    ! WMO region number                       [code table]
        case ('NIW')    ! type of instrument.for wind measurement [flag table]
        case ('MADDF0') ! associated field significance           [code table]
        case ('NPPPQ')  ! Q-bits for following value              [code table]
        case ('NPPP')   ! 3 hour pressure change                  [Pa]
        case ('NA')     ! characteristic of pressure tendency     [code table]
        case ('MVVQ')   ! Q-bits for following value              [code table]
        case ('MVV')    ! horizontal visibility                   [m]
        case ('NWWQ')   ! Q-bits for following value              [code table]
        case ('NWW')    ! present weather                         [code table]
        case ('MW1Q')   ! Q-bits for following value              [code table]
        case ('MW1')    ! past weather (1)                        [code table]
        case ('MW2Q')   ! Q-bits for following value              [code table]
        case ('MW2')    ! past weather (2)                        [code table]
        case ('MNQ')    ! Q-bits for following value              [code table]
        case ('MN')     ! cloud cover (total)                     [%]
        case ('NHQ')    ! Q-bits for following value              [code table]
        case ('NH')     ! height of base of cloud                 [m]
        case ('MCCQ')   ! Q-bits for following value              [code table]
        case ('MCC')    ! cloud type                              [code table]
        case ('MCC0Q')  ! Q-bits for following value              [code table]
        case ('MCC0')   ! cloud type                              [code table]
        case ('MCC1Q')  ! Q-bits for following value              [code table]
        case ('MCC1')   ! cloud type                              [code table]
        case ('NFXGU')  ! maximum wind speed(gusts)               [m/s]
        case ('MGGTR1') ! durat.of time relat.to following value  [hour]
        case ('MRRRQ')  ! Q-bits for following value              [code table]
        case ('MRRR')   ! total precipitation/total water equiv.  [kg/m**2]
        case ('MNHQ')   ! Q-bits for following value              [code table]
        case ('MNH')    ! cloud amount                            [code table]
        case ('MVTSU')  ! vertical significance (surface observ.) [code table]
        case ('MCC2')   ! cloud type                              [code table]
        case ('NH0Q')   ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('NH0')    ! height of base of cloud                 [m]
        case ('MNH0Q')  ! Q-bits for following value              [code table]
        case ('MNH0')   ! cloud amount                            [code table]
        case ('NSSS')   ! total snow depth                        [m]
        case ('MGGACT') ! actual hour of observation              [hour]
        case ('NGGACT') ! actual minute of observation            [minute]
        case ('MPW')    ! period of waves                         [s]
        case ('MHW')    ! height of waves                         [m]
        case ('MCOTR')  ! duration of time relat.to follow.value  [code table]
        case ('MGGTR3') ! durat.of time relat.to following value  [hour]
        case ('MSSSS')  ! global radiation integr.o.period specif.[J/m**2]
        case ('YSUPL')  ! 008 characters                          [CCITT IA5]
        case ('NFXME')  ! maximum wind speed(10 min mean wind)    [m/s]
        case ('MMOSTM') ! method of sea-surface temperature measu.[code table]
        case ('MTSQ')   ! Q-bits for following value              [code table]
        case ('MTS')    ! sea/water temperature                   [K]
        case ('MGGTR2') ! durat.of time relat.to following value  [hour]
        case ('MSSSM0') ! total sunshine                          [minute]
        case ('ME')     ! state of ground (with or without snow)  [code table]
        case ('MGGTR4') ! durat.of time relat.to following value  [hour]
        case ('MDFSR')  ! diff.solar radiation integr.o.per.spec. [J/m**2]
        case ('MPWPW')  ! period of wind waves                    [s]
        case ('MHWHW')  ! height of wind waves                    [m]
        case ('MGGTR5') ! durat.of time relat.to following value  [hour]
        case ('MLWR')   ! long wave radiat.,integr./period specif.[J/m**2]
        case ('ME0')    ! state of ground (with or without snow)  [code table]
        case ('MGGTR0Q')! Q-bits for following value              [code table]
        case ('MGGTR0') ! durat.of time relat.to following value  [hour]
        case ('MTNTNQ') ! Q-bits for following value              [code table]
        case ('MTNTN')  ! minimum temp., height/period specified  [K]
        case ('MSSSMQ') ! Q-bits for following value              [code table]
        case ('MSSSM')  ! total sunshine                          [minute]
        case ('MTXTXQ') ! Q-bits for following value              [code table]
        case ('MTXTX')  ! maximum temp., height/period specified  [K]
        case ('MGGTRQ') ! Q-bits for following value              [code table]
        case ('MGGTR')  ! durat.of time relat.to following value  [hour]
        case ('MTGTG')  ! ground minimum temp., past 12 hours     [K]
        case ('MVTSU0') ! vertical significance (surface observ.) [code table]
        case ('MNH1')   ! cloud amount                            [code table]
        case ('MCC3')   ! cloud type                              [code table]
        case ('MCT')    ! cloud top description                   [code table]
        case ('NDW12')  ! direction of swell waves                [degree true]
        case ('NHT')    ! height of top of cloud                  [m]
        case ('NP24Q')  ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('NP24')   ! 24 hour pressure change                 [Pa]
        case ('MEEE')   ! evaporation/evapotranspiration          [kg/m**2]
        case ('NIE')    ! instrument/crop for evapo(transpi)ration[code table]
        case ('MPW12')  ! period of swell waves                   [s]
        case ('MHW12')  ! height of swell waves                   [m]
        case ('MMOWTM') ! method of wet-bulb temperature measurem.[code table]
        case ('Loop000':'Loop999')  ! Start of Loop
        case ('Lcnt000':'Lcnt999')  ! Loop Counter
        case ('MDREP')  ! delayed descriptor replication factor      numeric
        !------
        ! METAR
        !------
        case ('NWWIP')  ! intensity or proximity of weather       [CODE_TABLE]
        case ('YWWP')   ! weather phenomena (METAR)               [CCITT_IA5]
        case ('NWWIP0') ! intensity or proximity of weather       [CODE_TABLE]
        case ('YWWP0')  ! weather phenomena (METAR)               [CCITT_IA5]
        case ('YWWD')   ! descriptor of weather                   [CCITT_IA5]
        case ('YWWD0')  ! descriptor of weather                   [CCITT_IA5]
        !-----------------------
        ! BUFR4 template (DRIBU)
        !-----------------------
        case ('MWRSA')  ! WMO REGION SUB-AREA                      [NUMERIC]
        case ('MTODB')  ! TYPE OF DATA BUOY                        [CODE_TABLE]
        case ('MTISI')  ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('MJJJ0')  ! YEAR                                     [YEAR]
        case ('MMM0')   ! MONTH                                    [MONTH]
        case ('MYY0')   ! DAY                                      [DAY]
        case ('MGG0')   ! HOUR                                     [HOUR]
        case ('NGG0')   ! MINUTE                                   [MINUTE]
        case ('MQBST')  ! QUALITY OF BUOY SATELLITE TRANSMISSION   [CODE_TABLE]
        case ('MQOBL')  ! QUALITY OF BUOY LOCATION                 [CODE_TABLE]
        case ('MLQC')   ! LOCATION QUALITY CLASS (..66 CONFIDENCE) [CODE_TABLE]
        case ('NDRTY')  ! DROGUE TYPE                              [CODE_TABLE]
        case ('MDCI')   ! DEPTH CORRECTION INDICATOR               [CODE_TABLE]
        case ('NZNZN')  ! DEPTH BELOW SEA/WATER SURFACE            [M]
        case ('MTN00Q') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTN00')  ! SEA/WATER TEMPERATURE                    [K]
        case ('MDREP0') ! DELAYED DESCRIPTOR REPLICATION FACTOR    [NUMERIC]
        case ('MTISI1') ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('NGGTP')  ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MTISI3') ! TIME SIGNIFICANCE                        [CODE_TABLE]
        case ('NOMDP')  ! OPERATOR OR MANUFACTURER DEF. PARAMETER  [NUMERIC]
        case ('MDREP1') ! DELAYED DESCRIPTOR REPLICATION FACTOR    [NUMERIC]
        case ('NDRDE')  ! DROGUE DEPTH                             [M]
        case ('NK2')    ! METHOD OF SALINITY/DEPTH MEASUREMENT     [CODE_TABLE]
        case ('MSNSN')  ! SALINITY                                 [CODE_TABLE]
        case ('MHAWAS1')! HEIGHT OF SENSOR ABOVE WATER SURFACE     [M]
        case ('MANTYP') ! ANEMOMETER TYPE                          [CODE_TABLE]
        !------------------------
        ! BUFR4 templates (SYNOP)
        !------------------------
        case ('MGGTP')  ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP1') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP2') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP3') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP4') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('NCI')    ! SEA ICE CONCENTRATION                    [CODE_TABLE]
        case ('MAMTI')  ! AMOUNT AND TYPE OF ICE                   [CODE_TABLE]
        case ('MICSI')  ! ICE SITUATION                            [CODE_TABLE]
        case ('MICDE')  ! ICE DEVELOPMENT                          [CODE_TABLE]
        case ('MSI')    ! ICE DEPOSIT (THICKNESS)                  [M]
        case ('MRS')    ! RATE OF ICE ACCREATION                   [CODE_TABLE]
        case ('NGGTP0') ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MTNTNHQ')! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTNTNH') ! MINIMUM TEMP., HEIGHT/PERIOD SPECIFIED   [K]
        case ('YSOSN')  ! STATION OR SITE NAME                     [CCITT_IA5]
        case ('MHOSEN') ! HEIGHT OF SENSOR ABOVE LOCAL GROUND      [M]
        case ('MVTSU2') ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)  [CODE_TABLE]
        case ('MGGTP5') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('NGGTP1') ! TIME PERIOD OR DISPLACEMENT              [MINUTE]
        case ('MGGTP6') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP0') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGGTP7') ! TIME PERIOD OR DISPLACEMENT              [HOUR]
        case ('MGLSR')  ! GLOBAL SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]
        case ('MTXTXHQ')! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MTXTXH') ! MAXIMUM TEMP., HEIGHT/PERIOD SPECIFIED   [K]
        case ('MRR24Q') ! Q-BITS FOR FOLLOWING VALUE               [CODE TABLE]
        case ('MRR24')  ! TOTAL PRECIPITATION, PAST 24 HOURS       [KG/M**2]
        case ('MDSRH')  ! DIFFUSE SOLAR RADIAT. (HIGH ACC.) INTEGR [J/M**2]
        case ('MEEEV')  ! EVAPORATION/EVAPOTRANSPIRATION           [KG/M**2]
        case ('MTGTGH') ! GROUND MINIMUM TEMP., PAST 12 HOURS      [K]
        !--------------------
        ! print unknown codes
        !--------------------
        case default
          call bufr_get_entry_texts (bufr)
          call bufr_get_entry_units (bufr)
          write(0,*) bufr% itype(id),bufr% ifxy(id)
          write(0,'(a,a,a)') ymnem, bufr% ytext(id),      &
                          '['//trim(bufr% yunit(id))//']'
          unused = .true.
#endif
        end select
      end if
    end do
    !--------------------------------------
    ! in case of unknown codes:
    ! print suspicious BUFR record and exit
    !--------------------------------------
    if (unused) then
      call bufr_print_sections (bufr)
      call bufr_print_subset   (bufr, is)
      call finish ('read_synop_bufr','code(s) not implemented')
    endif
    !-----------------------------------------------
    ! if station name is not set, use station number
    !-----------------------------------------------
    if (spt% statid   =='')   write(spt% statid,'(i5.5)') spt% ident
    !----------------------------------------------------
    ! 'QSCAT' or 'ASCAT' for Scatterometer coded as DRIBU
    !----------------------------------------------------
    s% z = s% zr                  ! Preset used height from reported height
    select case (spt% hd% dbkz)
    case (1697)
      spt% statid = 'QSCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1698,1699)
      spt% statid = 'ASCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1700)
      spt% statid = 'OSCAT'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1701)
      spt% statid = 'JASON-2'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1702)
      spt% statid = 'SARAL'
      s% z% o  = 0._sp
      s% z% qc = 0_i2
    case (1)
    !-----------------------------------------------------------
    ! metar: reported pressure is truncated, add fractional part
    !+++ disabled as height/bias correction is the way to go +++
    !-----------------------------------------------------------
!     if (s% p_msl% qc == 0) then
!       ip = int (s% p_msl% o)
!       if (mod(ip,100)==0) then
!         s% p_msl% o = s% p_msl% o + 50._sp    ! German METARs.+0.5hPa
!       else
!         s% p_msl% o = s% p_msl% o + 16.925_sp ! Internat.METARs+.16..
!       endif
!     endif
    end select
    !----------------
    ! standard checks
    !----------------
    lkeep = .true.
    call init_time (spt% actual_time, yyyy, mo, dd, hh, mi)
    call check_report_1 (spt)
    if (spt% use% state <= STAT_DISMISS) lkeep = .false.

    if (lkeep) then
      !-------------------------
      ! check for double entries
      !-------------------------
      select case (spt% hd% obstype)
      case (OT_SYNOP, OT_DRIBU)
      do i=1,obs% n_spot
       if  ( obs% spot(i)% hd% obstype == spt% hd% obstype ) then
        if ((obs% spot(i)% statid      == spt% statid) .or.&
            (obs% spot(i)% ident       /= 0 .and.    &
             obs% spot(i)% ident       == spt% ident )     ) then
         if( obs% spot(i)% actual_time == spt% actual_time ) then
          !-------------------------
          ! same station id and time
          !-------------------------
          call load_synop (obs, obs% spot(i), s1)
          call cmp_synop (s1, s, diff, m12, m21)
          if (obs%spot(i)% corme>0 .or. spt% corme>0) then
            !-------------------
            ! correction reports
            !-------------------
            if      (spt% corme > obs%spot(i)% corme) then       ! s corrects s1
              write (com,'("corme = ",i0)') spt% corme
              call decr_rpt_use (obs% spot(i), CHK_CORR, comment=com)
              call merge_synop (s1, s)
              call check_store_synop (s1, spt, obs, lkeep, repl=i)
              lkeep = .false.
              s     = s1        ! Retain saved data for diagnostic printout
            else if (spt% corme < obs%spot(i)% corme) then       ! s1 corrects s
              spt% corme = obs%spot(i)% corme
              write (com,'("corme = ",i0)') spt% corme
              call decr_rpt_use (obs% spot(i), CHK_CORR, comment=com)
              call merge_synop (s, s1)
              call check_store_synop (s , spt, obs, lkeep, repl=i)
              lkeep = .false.
            else
              write (com,'("corme = ",i0)') obs%spot(i)% corme
              call decr_rpt_use (obs% spot(i) ,CHK_CORRERR ,STAT_DISMISS, com)
              write (com,'("corme = ",i0)') spt% corme
              call decr_rpt_use (spt          ,CHK_CORRERR ,STAT_DISMISS, com)
              lkeep = .false.
            endif
          else
            !----------------------
            ! no correction reports
            !----------------------
            if (diff) then
              !-----------------------------
              ! conflicting double occurence
              !-----------------------------
              call decr_rpt_use (obs% spot(i) ,CHK_DBLERR ,STAT_DISMISS)
              call decr_rpt_use (spt          ,CHK_DBLERR ,STAT_DISMISS)
              lkeep = .false.
            elseif (m21) then
              !------------------------------
              ! non conflicting messages,
              ! the latter extends the former
              !------------------------------
              call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS)
              call check_store_synop (s, spt, obs, lkeep, repl=i)
              lkeep = .false.
            elseif (m12 .or. spt% hd% obstype /= OT_DRIBU) then
              !---------------------------------
              ! non conflicting messages,
              ! the former may extend the latter
              !---------------------------------
              call decr_rpt_use (spt  ,CHK_REDUNDANT ,STAT_DISMISS)
              lkeep = .false.
            else  ! .not.(m12.or.m21) .and. obstype == OT_DRIBU
              !---------------------------------
              ! non conflicting messages, buoys:
              ! compare quality of buoy location
              !---------------------------------
!------------------------------------------------------------------
! Quality of buoy location:
!  0 Reliable     (location was made over two satellite passes)
!  1 Latest known (no location over the corresponding
!  2 Dubious      (location made over one pass only; a second
!                  solution is possible in 5 per cent of the cases)
!  3 Missing value
!------------------------------------------------------------------
              if (s% qobl == s1% qobl) then
                 call decr_rpt_use (spt          ,CHK_REDUNDANT ,STAT_DISMISS)
                 lkeep = .false.
              else if (abs (s% qobl) > abs (s1% qobl)) then
                 call decr_rpt_use (spt          ,CHK_QI        ,STAT_DISMISS)
                 lkeep = .false.
              else
                 call decr_rpt_use (obs% spot(i) ,CHK_QI        ,STAT_DISMISS)
                 call check_store_synop (s, spt, obs, lkeep, repl=i)
                 lkeep = .false.
              end if
            endif
          endif
         !-----------------------------------
         ! same station id but different time
         !-----------------------------------
         else if (chk_dbl > 0) then
          dt1 = abs(minutes(obs% spot(i)% actual_time - ana_time))
          dt2 = abs(minutes(     spt    % actual_time - ana_time))
          if (dt1 <= dt2) then
            call decr_rpt_use (spt  ,CHK_THIN ,STAT_DISMISS)
          else
            call decr_rpt_use (obs% spot(i) ,CHK_THIN ,STAT_DISMISS)
            call check_store_synop (s, spt, obs, lkeep, repl=i)
          endif
          lkeep = .false.
         endif    ! time
        endif     ! statid
        if (.not.lkeep) exit
        !-----------------------------------
        ! different station at same location
        !-----------------------------------
        if (chk_dbl > 0                                     .and. &
            obs% spot(i)% col% c% dlat == spt% col% c% dlat .and. &
            obs% spot(i)% col% c% dlon == spt% col% c% dlon       ) then
          dt1 = abs(minutes(obs% spot(i)% actual_time - ana_time))
          dt2 = abs(minutes(     spt    % actual_time - ana_time))
          if (dt1 <= dt2) then
            call decr_rpt_use (spt  ,CHK_THIN ,STAT_DISMISS)
          else
            call decr_rpt_use (obs% spot(i) ,CHK_THIN ,STAT_DISMISS)
            call check_store_synop (s, spt, obs, lkeep, repl=i)
          endif
          lkeep = .false.
        endif
       endif      ! obstype
       if (.not.lkeep) exit
      end do
      end select
      !--------------------
      ! no double occurence
      !--------------------
      if (lkeep) call check_store_synop (s, spt, obs, lkeep)
    endif

    if (prt_data) call print_synop (spt, s)

  end subroutine read_synop_bufr
!------------------------------------------------------------------------------
  subroutine check_store_synop (synop, spot, obs, lkeep, repl)
  type(t_synop),intent(inout)        :: synop ! SYNOP level information
  type(t_spot) ,intent(inout)        :: spot  ! meta data of this observation
  type(t_obs)  ,intent(inout)        :: obs   ! data of all observations
  logical      ,intent(out)          :: lkeep ! flag: keep or reject
  integer      ,intent(in) ,optional :: repl  ! observation to replace

    !----------------
    ! local variables
    !----------------
    target                  :: synop
    integer, parameter      :: mv = 40    ! max. no. variables used
    integer                 :: typ (mv)   ! observation type
    real(wp)                :: lev (mv)   ! level of observation
    type(t_datum)           :: bod (mv)   ! body entry
    integer                 :: id
    type (t_spot)  ,pointer :: s
    integer                 :: n          ! number of valid observations
    integer                 :: i
    type (t_datum) ,pointer :: zs         ! link to synop% z_ls or synop% z_msl
    integer                 :: zs_to_z    ! copy of zls_to_z    or zmsl_to_z
    integer                 :: p_use(3)   ! copy of p_land, p_ship or p_buoy
    logical                 :: land, ship
    logical                 :: buoy, scat
    logical                 :: metar
    logical                 :: cman
    real(wp)                :: p_chk
    real(wp)                :: z_chk
    real(wp)                :: z_ref
    real(wp)                :: e_chk
    integer                 :: l_chk
    logical                 :: first = .true.
    logical                 :: use_v
    real(wp)                :: clct       ! total cloud cover (here: octas)
    logical                 :: ok
    real(sp)                :: tmin, tmax

    n     = 0
    lkeep = .false.
    land  = .false.
    ship  = .false.
    buoy  = .false.
    scat  = .false.
    metar = .false.
    cman  = .false.

    !----------------
    ! set preferences
    !----------------
    select case (spot% hd% buf_type)
    case (0)
      !------------------
      ! for land stations
      !------------------
      land    =  .true.
      p_use   =  p_land       ! pressure levels to use
      zs_to_z =  zls_to_z     ! station height  to use
      zs      => synop% z_ls
      select case (spot% hd% codetype)
      case (20)
        cman  = .true.        ! coastal station
      end select
    case (1, 12)
      !----------------------
      ! for ships (and buoys)
      !----------------------
      ship    =  .true.
      p_use   =  p_ship       ! pressure levels to use
      zs_to_z =  zmsl_to_z    ! station height  to use
      zs      => synop% z_msl
    case default
      call finish('check_store_synop','BUFR-type not in 0,1,12')
    end select
    select case (spot% hd% dbkz)
    case (1)
      !------
      ! METAR
      !------
      land    =  .true.
      metar   =  .true.
      p_use   =  p_metar
      zs_to_z =  zls_to_z
      zs      => synop% z_ls
    case (kz_buoy(1), kz_buoy(2))       ! (385, 10385)
      !-----------------------------------------------------------
      ! Hack for BUOYs, since BUFR subtype appears incorrectly set
      !-----------------------------------------------------------
      buoy    = .true.
      p_use   =  p_buoy       ! pressure levels to use
      zs_to_z =  zmsl_to_z    ! station height  to use
      zs      => synop% z_msl
    case default
!     scat    =  any (spot% hd% dbkz == kz_scat)
    end select

    !-------------------------------------------------------
    ! Consistency check: barometer height vs. station height
    ! Some stations erroneously report the height difference
    !-------------------------------------------------------
    if (land .and. chk_zb_mode > 0) then
      if (     synop% zr% qc == 0     .and. synop% z_ls% qc == 0          .and. &
               synop% zr% o  >= 0._sp .and. synop% zr  % o  <= chk_zb_max .and. &
          abs (synop% zr% o  -              synop% z_ls% o) >  dzs_zb_max) then
        if (chk_zb_mode == 1) then
           synop% zr% qc = QC_INCON
           synop% z % qc = QC_INCON
        else
           synop% z% o   = synop% z_ls% o + synop% zr% o
           synop% z% src = SRC_DER
           !print*, "setting z =", synop% z% o
        end if
        if (verbose > 0) then
          call nextline
          write(oline(iol),'(a,1x,a,1x,a,2(f10.2,1x),4x,3(a,f10.2,6x),a,i6,a)') &
               ' SYNOP-zb-check:', spot% statid, chhmm(spot% actual_time),    &
               spot% col% c% dlon, spot% col% c% dlat, 'zs=', synop% z_ls% o, &
               'zb=', synop% zr% o, 'z=', synop% z% o, 'z%qc=', synop% z% qc, &
               ' #CHK_ZB'
        end if
      end if
    end if

    !--------------------------------------------------
    ! set station height z from land surface height z_ls
    !--------------------------------------------------
    select case (zs_to_z)
    case (0)
      !--------------------------
      ! always use station height
      !--------------------------
    case (1)
      !----------------------------------------------------
      ! use surface height if station height is not present
      !----------------------------------------------------
      if (    synop% z  % qc /= 0 ) synop% z = zs
    case (2)
      !------------------------------------------------------
      ! use surface height if surface height < station height
      !------------------------------------------------------
      if (    synop% z  % qc /= 0 ) synop% z = zs
      if (           zs % qc == 0 &
        .and.        zs % o  <    &
              synop% z  % o       ) synop% z = zs
    case (3)
      !------------------------------------------------------
      ! use surface height if surface height > station height
      !------------------------------------------------------
      if (    synop% z  % qc /= 0 ) synop% z = zs
      if (           zs % qc == 0 &
        .and.        zs % o  >    &
              synop% z  % o       ) synop% z = zs
    case (4)
      !------------------------------------------------
      ! use surface height if surface height is present
      !------------------------------------------------
      if (    zs% qc == 0 ) synop% z = zs
    case (5)
      !--------------------------
      ! always use surface height
      !--------------------------
      synop% z = zs
    case default
      call finish ('check_store_synop','zls_to_z not in 0:5')
    end select
    !--------------------------------
    ! store station or surface height
    !--------------------------------
    if (synop% z% qc == 0) spot% z = synop% z% o

    !-----------------------------------
    ! Sanity check for reported pressure
    !-----------------------------------
    if (synop% ps    % o <= 0._sp ) synop% ps    % qc = QC_INCON
    if (synop% p     % o <= 0._sp ) synop% p     % qc = QC_INCON
    if (synop% p_msl % o <= 0._sp ) synop% p_msl % qc = QC_INCON

    if (buoy .or. ship) then
       if ( synop% ps   % qc == 0         .and. &
           (synop% ps   % o  <   85000._sp .or. &
            synop% ps   % o  >  110000._sp     )) synop% ps   % qc = QC_INCON
       if ( synop% p_msl% qc == 0         .and. &
           (synop% p_msl% o  <   85000._sp .or. &
            synop% p_msl% o  >  110000._sp     )) synop% p_msl% qc = QC_INCON
    else
       if ( synop% ps   % qc == 0         .and. &
           (synop% ps   % o  <   30000._sp .or. &
            synop% ps   % o  >  120000._sp     )) synop% ps   % qc = QC_INCON
       ! SYNOP land or METAR (QNH)
       if ( synop% p_msl% qc == 0         .and. &
           (synop% p_msl% o  <   80000._sp .or. &
            synop% p_msl% o  >  110000._sp     )) synop% p_msl% qc = QC_INCON
    end if

    !-------------------------------------------------------------
    ! Sanity check for reported temperature, dew-point temperature
    !-------------------------------------------------------------
    if (synop% t%  qc == QC_OK .and. (synop% t%  o < 100._sp .or. &
                                      synop% t%  o > 400._sp     )) then
        synop% t%  qc =  QC_CLIM  ! Temperature out of expected range
    end if
    if (synop% td% qc == QC_OK .and. (synop% td% o < 100._sp .or. &
                                      synop% td% o > 400._sp     )) then
        synop% td% qc =  QC_CLIM  ! Dew-point temperature out of range
    end if
    if (synop% rh% qc == QC_OK .and.  synop% rh% o < 0._sp) then
        synop% rh% qc =  QC_MISS  ! Unphysical rel.humidity
    end if

    if (buoy .or. ship) then
       if (synop% tsea% qc == QC_OK .and. (synop% tsea% o < 200._sp .or. &
                                           synop% tsea% o > 350._sp     )) then
           synop% tsea% qc =  QC_CLIM  ! Sea/water temperature out of range
       end if
    end if

    !-----------------------------------------
    ! check consistency of p_station and p_msl
    !-----------------------------------------
    l_chk = 0
    z_ref = -999._wp
    p_chk = -999._wp
    z_chk = -999._wp
    e_chk = 0._wp
    if (chk_ps(1) /= 0._wp) then
      if   (synop% ps  % qc == 0 .and. &     ! pressure at station height [Pa]
            synop% p   % qc == 0 .and. &     ! pressure level             [Pa]
            synop% gp  % qc == 0 .and. &     ! geopotential               [gpm]
            synop% z   % qc == 0 .and. &     ! height of station          [m]
            synop% t   % qc == 0       ) then
        l_chk = P_GP
        p_chk = synop% p%  o
        z_ref = synop% gp% o
        z_chk = zbs (real ( synop% p%  o, wp), &
                     real ( synop% ps% o, wp), &
                     real ( synop% t%  o, wp), &
                     gacc * synop% z%  o       ) / gacc

      endif
      if   (synop% ps   % qc == 0 .and. &
            synop% p_msl% qc == 0 .and. &
            synop% z    % qc == 0 .and. &
            synop% t    % qc == 0       ) then
        l_chk = P_MSL
        p_chk = synop% p_msl% o
        z_ref = 0._wp
        z_chk = zbs (real ( synop% p_msl% o, wp), &
                     real ( synop% ps%    o, wp), &
                     real ( synop% t%     o, wp), &
                     gacc * synop% z%     o       ) / gacc
      endif
      if (l_chk /= 0) then
        e_chk = sqrt ( (chk_ps(1)                      )**2 + &
                       (chk_ps(2)*(synop% z% o - z_ref))**2   )
        if ( abs (z_ref - z_chk) > e_chk ) then
          synop% z% qc = QC_INCON
        endif
      else
        if (chk_ps(1) < 0._wp) then
          synop% z% qc = QC_INCON
        endif
      endif
      if (synop% z% qc == QC_INCON) then
        if (btest (chk_strict,0)) synop% p_msl% qc = QC_INCON  ! Distrust p_msl
        if (btest (chk_strict,1)) synop% ps%    qc = QC_INCON  ! Distrust p_s
      end if
      !---------------------------------------------------------
      ! Plausibility checks: (diurnal) air temperature variation
      !---------------------------------------------------------
      if (chk_tmaxmin > 0) then
        if (btest (chk_tmaxmin, 2) .and. ag_tmaxmin) then
           !---------------------------------------
           ! Tmax/Tmin: 12-hourly vs. 1-hourly data
           ! (also from different/merged reports).
           ! Allow some small tolerance (~ 0.1 K).
           ! If inconsistent, try to recover later
           ! from aggregation of 1-hourly data.
           !---------------------------------------
           if (synop% tmax1 /= rvind) then
              tmax = synop% tmax
              if (tmax /= rvind .and. tmax > 0._sp                 &
                                .and. tmax < synop% tmax1 - 0.11_sp) then
                 synop% tmax = -99._sp  ! Provide placeholder for aggregation
              end if
           end if
           if (synop% tmin1 /= rvind) then
              tmin = synop% tmin
              if (tmin /= rvind .and. tmin > 0._sp                 &
                                .and. tmin > synop% tmin1 + 0.11_sp) then
                 synop% tmin = -99._sp
              end if
           end if
        end if

        if (btest (chk_tmaxmin, 0) .and. synop% t% qc == QC_OK) then
           !----------------------------------
           ! Tmax/Tmin vs. T2M: 12-hourly data
           ! allow tolerance for North-America
           !----------------------------------
           tmax = synop% tmax
           tmin = synop% tmin
           if (tmax <  0._sp) tmax = rvind
           if (tmin <  0._sp) tmin = rvind
           if (tmax /= rvind .and. (tmax < synop% t% o - 0.21_sp .or. &
                                    tmax > synop% t% o + 60.0_sp)     ) then
              synop% tmax = rvind
           end if
           if (tmin /= rvind .and. (tmin > synop% t% o + 0.21_sp .or. &
                                    tmin < synop% t% o - 60.0_sp)     ) then
              synop% tmin = rvind
           end if
           !---------------------------------
           ! Tmax/Tmin vs. T2M: 1-hourly data (slightly tighter limits)
           !---------------------------------
           tmax = synop% tmax1
           tmin = synop% tmin1
           if (tmax /= rvind .and. (tmax < synop% t% o - 0.11_sp .or. &
                                    tmax > synop% t% o + 30.0_sp)     ) then
              synop% tmax1 = rvind
           end if
           if (tmin /= rvind .and. (tmin > synop% t% o + 0.11_sp .or. &
                                    tmin < synop% t% o - 30.0_sp)     ) then
              synop% tmin1 = rvind
           end if
        end if

        if (btest (chk_tmaxmin, 1)) then
           !------------------------------
           ! Tmax vs. Tmin: 12-hourly data
           !------------------------------
           tmax = synop% tmax
           tmin = synop% tmin
           if (tmax <  0._sp) tmax = rvind
           if (tmin <  0._sp) tmin = rvind
           if (tmax /= rvind .and. tmin /= rvind .and. &
                (tmax < tmin .or. tmax > tmin + 60._sp)) then
              synop% tmax = rvind
              synop% tmin = rvind
           end if
           !-----------------------------
           ! Tmax vs. Tmin: 1-hourly data
           !-----------------------------
           tmax = synop% tmax1
           tmin = synop% tmin1
           if (tmax /= rvind .and. tmin /= rvind .and. &
                (tmax < tmin .or. tmax > tmin + 30._sp)) then
              synop% tmax1 = rvind
              synop% tmin1 = rvind
           end if
        end if
      end if
      !--------------------
      ! diagnostic printout
      !--------------------
      if ((synop% p_as% qc == 0 .or. &
           synop% ps  % qc == 0 .or. &
           synop% p   % qc == 0      ) .and. verbose > 0) then
        if (first) then
           call nextline
           write(oline(iol), '()')
           call nextline
           write(oline(iol), '(A)') " SYNOP-ps-check: &
                &station    time       lon        lat        z%o       ps%o &
                &chk z% p%qc   p_chk     z_ref     z_chk     e_chk"
           first = .false.
        end if
        if (l_chk == 0 .or. synop% z% qc /= 0 .or. verbose > 1) then
          call nextline
          write(oline(iol), '(a,1x,a,1x,a,4(f10.2,1x),3i3,4f10.2,a)')         &
           ' SYNOP-ps-check:', spot% statid, chhmm(spot% actual_time),        &
           spot% col% c% dlon, spot% col% c% dlat, synop% z% o, synop% ps% o, &
           l_chk, synop% z% qc, synop% ps% qc, p_chk, z_ref, z_chk, e_chk,    &
           ' #CHK_PS'
        end if
      endif
    endif
    if (metar) then
       !-------------------------
       ! Convert QNH to p_station
       !-------------------------
       if (synop% p_msl% qc == QC_OK .and. &
           synop% z    % qc == QC_OK .and. &
           synop% ps   % qc /= QC_OK       ) then
          synop% ps% o   = synop% p_msl% o * p_h_usstd (real (synop% z% o,wp)) &
                                           / p_r
          synop% ps% qc  = QC_OK
          synop% ps% src = SRC_DER
       end if
    end if
    if(prt_data) print *, repeat('=',64)
    !------------------------------
    ! select pressure to assimilate
    !------------------------------
    do i=1,3
      select case (p_use(i))
      !---------------------------
      ! pressure at station height
      !---------------------------
      case (P_S)
        if (synop% p_as% qc /= 0 .and. &
            synop% ps  % qc == 0 .and. &
            synop% z   % qc == 0       ) then
          synop% p_as   = synop% ps
          synop% z_as   = synop% z
          synop% p_used = p_use(i)
          if(prt_data) print *,'p_used: P_S'
        endif
      !------------------------------------
      ! pressure at standard pressure level
      !------------------------------------
      case (P_GP)
        if   (synop% p_as% qc             /= 0 .and. &
              synop% p   % qc             == 0 .and. &
              synop% gp  % qc             == 0 .and. &
              synop% z   % qc             == 0 .and. &
          abs(synop% z% o - synop% gp% o) <= z_stat  ) then
          synop% p_as   = synop% p
          synop% z_as   = synop% gp
          synop% p_used = p_use(i)
          if(prt_data) print *,'p_used: P_GP (z-gp ok)'
        endif
        if     (synop% p_as% qc             /= 0 .and. &
                synop% p   % qc             == 0 .and. &
                synop% gp  % qc             == 0 .and. &
                synop% ps  % qc             == 0 .and. &
            abs(synop% p% o - synop% ps% o) <= p_stat  ) then
          synop% p_as   = synop% p
          synop% z_as   = synop% gp
          synop% p_used = p_use(i)
          if(prt_data) print *,'p_used: P_GP (p-ps ok)'
        endif
      !---------------------------
      ! pressure at mean sea level
      !---------------------------
      case (P_MSL)
        if (synop% p_as % qc /= 0 .and. &
            synop% p_msl% qc == 0 .and. &
            synop% z    % qc == 0 .and. &
            synop% z    % o  <= z_stat  ) then
          synop% p_as   = synop% p_msl
          synop% z_as   = synop% z_msl
          synop% p_used = p_use(i)
          if(prt_data) print *,'p_used: P_MSL (z-gp ok)'
        endif
        if     (synop% p_as % qc                /= 0 .and. &
                synop% p_msl% qc                == 0 .and. &
                synop% ps   % qc                == 0 .and. &
            abs(synop% p_msl% o - synop% ps% o) <= p_stat  ) then
          synop% p_as   = synop% p_msl
          synop% z_as   = synop% z_msl
          synop% p_used = p_use(i)
          if(prt_data) print *,'p_used: P_MSL (p-ps ok)'
        endif
      !-----
      ! none
      !-----
      case (0)
        synop% p_used = 0
      case default
        call finish('check_store_synop',&
                    'p_land, p_ship, p_buoy or p_metar not in 0,1,2,3')
      end select
    end do

    !------------------------------------------
    ! select reference pressure for t, rh, u, v
    !------------------------------------------
    synop% p_ref = synop% ps
    if (synop% p_ref % qc /= 0 .and.     &
        synop% p_msl % qc == 0 .and.     &
                         (ship .or. buoy)) synop% p_ref     = synop% p_msl
    if (synop% p_ref % qc /= 0 .and.     &
        synop% ps    % qc == 0           ) synop% p_ref     = synop% ps
    if (synop% p_ref % qc /= 0)            synop% p_ref % o = use_ps_model

    !--------------------------------------------
    ! insert in list if information is sufficient
    !--------------------------------------------
    if(prt_data) print *,'valid: lon lat z p',&
        spot% col% c% dlon     /= invalid,    &
        spot% col% c% dlat     /= invalid,    &
        synop%        z_as% qc == 0      ,    &
        synop%        p_as% qc == 0

    if (spot% col% c% dlon     /= invalid .and. &
        spot% col% c% dlat     /= invalid       ) then
      if (synop%      z_as% qc == 0       .and. &
          synop%      p_as% qc == 0             ) then
        !-------------------------------------
        ! process individual observation types
        !-------------------------------------
        !-------------------------------------------
        ! geopotential height vs p at station or msl
        !-------------------------------------------
        if (use_h) then
          n = n + 1
          if(prt_data) print *,'z : n =',n
          typ (n) = VN_Z
          lev (n) = synop% p_as% o
          bod (n) = synop% z_as
        endif
        !--------------------------
        ! surface (or msl) pressure
        !--------------------------
        if (use_p) then
          n = n + 1
          if(prt_data) print *,'p : n =',n
          typ (n)          = VN_PS
          lev (n)          = synop% z_as% o
          bod (n)          = synop% p_as
          bod (n)% lev_typ = VN_Z
          bod (n)% plev    = synop% p_as% o
        endif
      endif

      !--------------------------------
      ! check height of sensors (DRIBU)
      !--------------------------------
      if (buoy) then
         if ( synop% h_ts% qc == 0            .and. &
             (synop% h_ts% o  <  buoy_hr_t(1) .or.  &
              synop% h_ts% o  >  buoy_hr_t(2)      )) synop% h_ts% qc = QC_NOUSE
         !-------------------------------------------------------------
         ! Replace missing wind sensor height by height of temp. sensor
         !-------------------------------------------------------------
         if ( synop% h_ws% qc /= 0 .and. synop% h_ts% qc == 0) &
              synop% h_ws = synop% h_ts

         if ( synop% h_ws% qc == 0            .and. &
             (synop% h_ws% o  <  buoy_hr_v(1) .or.  &
              synop% h_ws% o  >  buoy_hr_v(2)      )) synop% h_ws% qc = QC_NOUSE
      end if

      !------------------------------
      ! temperature at station height
      !------------------------------
      if (use_t2) then
        if (synop% t%     qc == QC_OK         .and. &
!           synop% t%     o  /= invalid       .and. &
!          (synop% ps%    qc == QC_OK .or.          &
!           synop% p_msl% qc == QC_OK       ) .and. & ! needed for METAR
            .not. cman                        .and. & ! CMAN lacks sensor height
            spot% z          /= invalid       .and. &
            spot% z          /= empty_spot% z       ) then
          n = n + 1
          if(prt_data) print *,'t : n =',n
          typ (n)             = VN_T2M
          bod (n)             = synop% t
          if (synop% ps%   qc == QC_OK) then
            bod (n)% plev     = synop% ps% o
          else
            bod (n)% plev     = synop% p_ref% o
          end if
          lev (n)             = bod (n)% plev   ! synop% ps% o

          if (land .and. .not. cman .and. version >= 3) then
             bod (n)% lev_typ = VN_HOSAG
             if (synop% h_ts% qc == 0) then
                lev (n)       = synop% h_ts% o  ! use reasonable sensor height
             else
                lev (n)       = 2._wp           ! default for SYNOP land
             end if
             if (lev(n) < hr_t2m(1) .or. lev(n) > hr_t2m(2))                 &
                call decr_use (bod(n)% use, STAT_PASSIVE, check=CHK_BLACKLIST)
          end if
          if (buoy .and. synop% h_ts% qc == 0) then
             bod (n)% lev_typ = VN_HOSAG        ! was VN_HEIGHT previously
             lev (n)          = synop% h_ts% o
          end if
          !---------------------------------------
          ! dewpoint temperature (monitoring only)
          !---------------------------------------
          if (mon_td) then
            if (synop% td % qc == QC_OK .and. &
                synop% t  % qc == QC_OK       ) then
              n = n + 1
              if(prt_data) print *,'td: n =',n
              typ (n)             = VN_TD2M
              bod (n)             = synop% td
              bod (n)% plev       = synop% p_ref% o
              lev (n)             = synop% p_ref% o

              if (land .and. .not. cman .and. version >= 3) then
                 bod (n)% lev_typ = VN_HOSAG
                 if (synop% h_ts% qc == 0) then
                    lev (n)       = synop% h_ts% o  ! use reasonable sensor height
                 else
                    lev (n)       = 2._wp           ! default for SYNOP land
                 end if
                 if (lev(n) < hr_t2m(1) .or. lev(n) > hr_t2m(2))                 &
                    call decr_use (bod(n)% use, STAT_PASSIVE, check=CHK_BLACKLIST)
              end if
              if (buoy .and. synop% h_ts% qc == 0) then
                 bod (n)% lev_typ = VN_HOSAG        ! was VN_HEIGHT previously
                 lev (n)          = synop% h_ts% o
              end if
            endif
          endif
        end if
      endif
      !----------------------------------
      ! rh, if not present derive from td
      !----------------------------------
      if (use_rh) then
        ok = .false.
        if (synop% rh% qc == QC_OK) then
          ok = .true.
        else if (synop% td % qc == 0 .and.     &
                 synop% t  % qc == 0 .and.     &
                 synop% td % o  <= synop% t % o) then
          synop% rh% o   = min (1._wp, max( 0._wp,           &
                           esw_t (real (synop% td% o, wp)) / &
                           esw_t (real (synop% t % o, wp))   ))
          synop% rh% qc  = QC_OK
          synop% rh% src = SRC_DER
          ok = .true.
        end if
        if (ok) then
          n = n + 1
          if(prt_data) print *,'rh: n =',n
          typ (n)             = VN_RH
          bod (n)             = synop% rh
          bod (n)% plev       = synop% p_ref% o
          lev (n)             = synop% p_ref% o

          if (land .and. .not. cman .and. version >= 3) then
             bod (n)% lev_typ = VN_HOSAG
             if (synop% h_ts% qc == 0) then
                lev (n)       = synop% h_ts% o  ! use reasonable sensor height
             else
                lev (n)       = 2._wp           ! default for SYNOP land
             end if
          end if
          if (buoy .and. synop% h_ts% qc == 0) then
             bod (n)% lev_typ = VN_HOSAG        ! was VN_HEIGHT previously
             lev (n)          = synop% h_ts% o
          end if
        endif
      endif
      !----------------------------------------
      ! use wind over sea,
      ! use wind over land only: in the tropics
      !                          below 150 m
      !----------------------------------------
      use_v = .false.
      if (land) then
         if (spot% z /= invalid       .and. &
             spot% z /= empty_spot% z       ) then
            use_v = (     use_vl                               &
                    .and. abs (spot% col% c% dlat) <=  30._wp  &
                    .and.      spot% z             <  150._wp) &
                    .or.  use_vl_gl
            if(prt_data) print *,'land: spot%z, synop%z_ls,qc, synop%zr,qc =',&
                 spot%z, synop%z_ls%o,synop%z_ls%qc,synop%zr%o,synop%zr%qc
         end if
         use_v = use_v  .and. .not. cman        ! CMAN lacks anemometer height
      else
         use_v = use_vs .and. (ship .or. buoy)
      end if
      if (use_v) then

        if      (synop% ff% qc == 0 .and. synop% dd% qc == 0) then
          synop% uu% src = SRC_DER
          synop% uu% qc  = QC_OK
          synop% uu% o   = synop% ff% o * (-sin (d2r * synop% dd% o))
          synop% vv% src = SRC_DER
          synop% vv% qc  = QC_OK
          synop% vv% o   = synop% ff% o * (-cos (d2r * synop% dd% o))

          n = n + 1
          if(prt_data) print *,'u : n =',n
          typ (n)             = VN_U
          bod (n)             = synop% uu
          bod (n)% plev       = synop% p_ref% o
          if (buoy .and. synop% h_ws% qc == 0) then
             bod (n)% lev_typ = VN_HEIGHT
             lev (n)          = synop% h_ws%  o
          else
             lev (n)          = synop% p_ref% o
          end if
          n = n + 1
          if(prt_data) print *,'v : n =',n
          typ (n)             = VN_V
          bod (n)             = synop% vv
          bod (n)% plev       = synop% p_ref% o
          if (buoy .and. synop% h_ws% qc == 0) then
             bod (n)% lev_typ = VN_HEIGHT
             lev (n)          = synop% h_ws%  o
          else
             lev (n)          = synop% p_ref% o
          end if

        end if

        if (synop% ff% qc == 0 .and.                                         &
            ((scat   .and. synop% dd% qc /= 0) .or.                          &
             (monitor_ff .and. synop% ff% qc == 0 .and. synop% dd% qc == 0)))&
          then
          !----------------------------------------------------------
          ! Scatterometer-like wind speed data without wind direction
          ! or FF monitoring
          !----------------------------------------------------------
          n = n + 1
          if (prt_data) print *,'ff: n =',n
          typ (n)             = VN_FF
          bod (n)             = synop% ff
          bod (n)% plev       = synop% p_ref% o
          if (buoy .and. synop% h_ws% qc == 0) then
             bod (n)% lev_typ = VN_HEIGHT
             lev (n)          = synop% h_ws%  o
          else
             lev (n)          = synop% p_ref% o
          end if
        end if

        if (synop% ff% qc == 0 .and. synop% ff% o > 0. .and. &
            synop% dd% qc == 0 .and. monitor_dd              ) then
          n = n + 1
          if (prt_data) print *,'dd: n =',n
          typ (n)             = VN_DD
          bod (n)             = synop% dd
          bod (n)% plev       = synop% p_ref% o
          if (buoy .and. synop% h_ws% qc == 0) then
             bod (n)% lev_typ = VN_HEIGHT
             lev (n)          = synop% h_ws%  o
          else
             lev (n)          = synop% p_ref% o
          end if
        end if

      endif

      if (mon_tsea .and. synop% tsea% qc == QC_OK  &
                   .and. (ship .or. buoy .or. cman)) then
         n = n + 1
         if (prt_data) print *,'tsea: n =',n
         typ (n)             = VN_TSEA
         bod (n)             = synop% tsea
         bod (n)% plev       = synop% p_ref% o
         bod (n)% lev_typ    = VN_DEPTH
         if (synop% h_tsea% qc == QC_OK) then
            lev (n)          = synop% h_tsea% o     ! Depth below sea surface
         else
            lev (n)          = invalid
         end if
      end if

      if (mon_cl) then
        !----------------
        ! cloud variables
        !----------------
        if (synop% n /= -1) then
          n = n + 1
          if (prt_data) print *,'n : n =',n
          typ (n) = VN_N
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          clct = octa_percent (real (synop% n, wp))
          call set_datum (bod (n) ,real (clct, sp) ,spot% corme)
        endif

        if (synop% nh /= rvind) then
          n = n + 1
          if (prt_data) print *,'nh: n =',n
          typ (n) = VN_NH
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% nh ,spot% corme)
        endif

        if (synop% ceil /= rvind) then
          n = n + 1
          if (prt_data) print *,'ceiling: n =',n
          typ (n) = VN_CEIL
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% ceil ,spot% corme)
        endif

        if (synop% vis /= rvind) then
          n = n + 1
          if (prt_data) print *,'visibility: n =',n
          typ (n) = VN_VV
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% vis ,spot% corme)
        endif

!       if (synop% nn >= 0._sp) then
!         n = n + 1
!         if (prt_data) print *,'nn: n =',n
!         typ (n) = VN_N_L
!         lev (n) = synop% p_ref% o
!         bod (n) = inv_datum
!         call set_datum (bod (n) ,synop% nn ,spot% corme)
!       endif

        if (synop% clcl >= 0._sp) then
          n = n + 1
          if (prt_data) print *,'clcl: n =',n
          typ (n) = VN_N_L
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% clcl ,spot% corme)
        endif

        if (synop% clcm >= 0._sp) then
          n = n + 1
          if (prt_data) print *,'clcm: n =',n
          typ (n) = VN_N_M
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% clcm ,spot% corme)
        endif

        if (synop% clch >= 0._sp) then
          n = n + 1
          if (prt_data) print *,'clch: n =',n
          typ (n) = VN_N_H
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% clch ,spot% corme)
        endif
      endif

      if (mon_cwoo) then
        !----------------
        ! cloud / weather related variables without obs operator
        !----------------
        if (synop% ww /= nvind) then
          n = n + 1
          if (prt_data) print *,'present weather: n =',n
          typ (n) = VN_WW
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,real(synop% ww) ,spot% corme)
        endif

        if (synop% gwg /= nvind) then
          n = n + 1
          if (prt_data) print *,'general cloud group: n =',n
          typ (n) = VN_GCLG
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,real(synop% gwg) ,spot% corme)
        endif

        if (synop% icl1 /= nvind) then
          n = n + 1
          if (prt_data) print *,'individual cloud layer 1: n =',n
          typ (n) = VN_ICLG
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,real(synop% icl1) ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_NUM
          if (synop% icl2 /= nvind) then
            n = n + 1
            if (prt_data) print *,'individual cloud layer 2: n =',n
            typ (n) = VN_ICLG
            lev (n) = 2._wp
            bod (n) = inv_datum
            call set_datum (bod (n) ,real(synop% icl2) ,spot% corme)
            bod (n)% plev    = synop% p_ref% o
            bod (n)% lev_typ = VN_NUM
            if (synop% icl3 /= nvind) then
              n = n + 1
              if (prt_data) print *,'individual cloud layer 3: n =',n
              typ (n) = VN_ICLG
              lev (n) = 3._wp
              bod (n) = inv_datum
              call set_datum (bod (n) ,real(synop% icl3) ,spot% corme)
              bod (n)% plev    = synop% p_ref% o
              bod (n)% lev_typ = VN_NUM
              if (synop% icl4 /= nvind) then
                n = n + 1
                if (prt_data) print *,'individual cloud layer 4: n =',n
                typ (n) = VN_ICLG
                lev (n) = 4._wp
                bod (n) = inv_datum
                call set_datum (bod (n) ,real(synop% icl4) ,spot% corme)
                bod (n)% plev    = synop% p_ref% o
                bod (n)% lev_typ = VN_NUM
              endif
            endif
          endif
        endif

      endif

      if (mon_tr) then
        !----------------------
        ! time range variables:
        ! wind gust: 1,3,6 h
        !----------------------
        if (synop% gust1 >= 0._sp) then
          n = n + 1
          typ (n) = VN_GUST
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gust1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% gust3 >= 0._sp) then
          n = n + 1
          typ (n) = VN_GUST
          lev (n) = 3._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gust3 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% gust6 >= 0._sp) then
          n = n + 1
          typ (n) = VN_GUST
          lev (n) = 6._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gust6 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        !--------------------
        ! rain: 1,3,6,12,24 h
        !--------------------
        if (synop% rr1 > -0.11) then
          n = n + 1
          typ (n) = VN_RR
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% rr1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% rr3 > -0.11) then
          n = n + 1
          typ (n) = VN_RR
          lev (n) = 3._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% rr3 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% rr6 > -0.11) then
          n = n + 1
          typ (n) = VN_RR
          lev (n) = 6._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% rr6 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% rr12 > -0.11) then
          n = n + 1
          typ (n) = VN_RR
          lev (n) = 12._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% rr12 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% rr24 > -0.11) then
          n = n + 1
          typ (n) = VN_RR
          lev (n) = 24._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% rr24 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        !--------------------
        ! long wave radiation
        !--------------------
        if (synop% lwr1 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_LW
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% lwr1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        !-----------------
        ! global radiation
        !-----------------
        if (synop% gsr1 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_GL
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gsr1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% gsr3 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_GL
          lev (n) = 3._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gsr3 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% gsr6 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_GL
          lev (n) = 6._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% gsr6 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        !------------------
        ! diffuse radiation
        !------------------
        if (synop% dsr1 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_DF
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% dsr1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% dsr3 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_DF
          lev (n) = 3._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% dsr3 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% dsr6 >= 0._sp) then
          n = n + 1
          typ (n) = VN_RAD_DF
          lev (n) = 6._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% dsr6 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        !--------------------
        ! tmin, tmax: 1h, 12h
        !--------------------
        if (synop% tmin1 /= rvind) then
          n = n + 1
          typ (n) = VN_TMIN
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% tmin1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% tmax1 /= rvind) then
          n = n + 1
          typ (n) = VN_TMAX
          lev (n) = 1._wp
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% tmax1 ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% tmin /= rvind) then
          n = n + 1
          typ (n) = VN_TMIN
          lev (n) = 12._wp
          bod (n) = inv_datum
          if (synop% tmin < 0._sp) synop% tmin = rvind
          call set_datum (bod (n) ,synop% tmin ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif

        if (synop% tmax /= rvind) then
          n = n + 1
          typ (n) = VN_TMAX
          lev (n) = 12._wp
          bod (n) = inv_datum
          if (synop% tmax < 0._sp) synop% tmax = rvind
          call set_datum (bod (n) ,synop% tmax ,spot% corme)
          bod (n)% plev    = synop% p_ref% o
          bod (n)% lev_typ = VN_TRTR
        endif
!       !-----------------------
!       ! pressure tendency: 3 h
!       !-----------------------
!       if (synop% p_tend /= rvind) then
!         n = n + 1
!         typ (n) = VN_PTEND
!         lev (n) = 3._wp
!         bod (n) = inv_datum
!         call set_datum (bod (n) ,synop% p_tend ,spot% corme)
!         bod (n)% plev    = synop% p_ref% o
!         bod (n)% lev_typ = VN_TRTR
!       endif
      endif

      if (mon_snow .and. land) then
        if (synop% sdepth /= rvind .and. synop% sdepth <  40._sp) then
          n = n + 1
          if (prt_data) print *,'sdepth: n =',n
          typ (n) = VN_SDEPTH
          lev (n) = synop% p_ref% o
          bod (n) = inv_datum
          call set_datum (bod (n) ,synop% sdepth ,spot% corme)
        endif
      end if

      !-----------------------------------------
      ! set invalid value for reference pressure
      !-----------------------------------------
    endif

    !-----------------------------------------------------------------
    ! Catch forgotten adjustment of mv if compiled w/o bounds checking
    !-----------------------------------------------------------------
    if (n > mv) call finish('check_store_synop','increase mv')

    !-----------------------------------------------
    ! For BUOYS, derive pcc from quality of location
    !-----------------------------------------------
    if (buoy) then
       select case (synop% qobl)
       case (0)
          spot% pcc = 100
       case (1)
          spot% pcc =  50
       case default
          spot% pcc =   0
       end select
       bod(:n)% pcc = spot% pcc
       if (prt_data) print *,'pcc   =', spot% pcc
    end if

!   if (prt_data) call print_synop (spot, synop)

    !-----------------------------------------
    ! Drop bogus reports at exactly (0N,0E deg.).
    !-----------------------------------------
    if (spot% col% c% dlat == 0._wp .and. &
        spot% col% c% dlon == 0._wp) then
       n = 0
    end if

    !--------------------
    ! store data into OBS
    !--------------------
    if (n>0) then
      lkeep = .true.
      if (present (repl)) then
        s => obs% spot (repl)
      else
        call new_spot (obs, 1, set_id=.true.)
        s => obs% spot (obs% n_spot)
      endif
      id = s% id
      s = spot
!     s% int_type  = ITY_ICOL
      s% id        = id
      s% col% nlev = 1
!     s% cost      = 1._wp
!     s% char      = CHR_ID
!     s% nr        = n
      call new_obs (obs, n, s)
!     call new_int (obs, s, n)
      spot% o      = s% o      ! return observation space info
      obs% varno (s%o%i+1:s%o%i+s%o%n) =      typ (1:n)
      obs% body  (s%o%i+1:s%o%i+s%o%n) =      bod (1:n)
!     obs% t_int (s%i%i+1:s%i%i+s%i%n) =      typ (1:n)
      obs%  olev (s%o%i+1:s%o%i+s%o%n) =      lev (1:n)
!     obs%   lev (s%i%i+1:s%i%i+s%i%n) = log (lev (1:n))
      call store_synop (obs, s, synop)
      call set_xuv (s)
      !--------------------------------------
      ! certain quantities are always passive
      !--------------------------------------
      select case (spot% hd% dbkz)
      case (5, 10005)
        !---------------------------------------------------------------
        ! German national BUFR template used only for RR,GUST monitoring
        !---------------------------------------------------------------
        do i = 1, s% o% n
          select case   (obs% varno(s% o% i+i))
          case (VN_GUST, VN_RR)
            call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
          case default
            call decr_use (obs% body (s% o% i+i)% use, STAT_DISMISS,  check=CHK_NOTUSED)
          end select
        end do
      case default
        do i = 1, s% o% n
          select case   (obs% varno(s% o% i+i))
          case (VN_T, VN_T2M)
            if (.not. land .or. cman .or. version < 3) then
             call decr_use (obs% body(s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
!           else if () then
!            call decr_use (obs% body(s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
            end if
          case (VN_TD, VN_TD2M, VN_TSEA)
            call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
          case (VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND, VN_N, VN_NH, VN_CEIL, VN_N_L&
               ,VN_N_M, VN_N_H, VN_VV, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW, VN_SDEPTH)
!           call decr_use (obs% body (s% o% i+i)% use, STAT_OBS_ONLY, check=CHK_NOTUSED)
            call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
          case (VN_WW, VN_GCLG, VN_ICLG)
            call decr_use (obs% body (s% o% i+i)% use, STAT_OBS_ONLY, check=CHK_NOTUSED)
          case (VN_FF, VN_DD)
            call decr_use (obs% body (s% o% i+i)% use, STAT_PASSIVE,  check=CHK_NOTUSED)
          end select
        end do
      end select
    else
      call decr_rpt_use (spot ,CHK_INSDAT, STAT_DISMISS)
    endif
  end subroutine check_store_synop
!==============================================================================
  pure subroutine construct_synop (s)
  type (t_synop) ,intent(out) :: s

#if defined(__SX__)
    ! Default initialisation does not work with sxf90 rev.360 (SX-6)
!   s = inva_syn
#endif
    s% total %mn =''
    s% z     %mn ='z'
    s% p     %mn ='p'
    s% gp    %mn ='gp'
    s% ps    %mn ='ps'
    s% p_msl %mn ='pmsl'
    s% z_msl %mn ='zmsl'
    s% h_ts  %mn ='h_ts'
    s% t     %mn ='t'
    s% td    %mn ='td'
    s% rh    %mn ='rh'
    s% h_ws  %mn ='h_ws'
    s% ff    %mn ='ff'
    s% dd    %mn ='dd'
    s% uu    %mn ='u'
    s% vv    %mn ='v'
    s% tsea  %mn ='tsea'
    s% h_tsea%mn ='h_tsea'
    s% z_msl %o  = 0.
    s% z_msl %qc = 0_i2

  end subroutine construct_synop
!==============================================================================
  subroutine print_synop (s, d)
  type (t_spot ) ,intent(in) :: s
  type (t_synop) ,intent(in) :: d

    write (6,*)
    write (6,*)' statid               = ',s% statid
    write (6,*)' lat lon              = ',s% col% c% dlat, s% col% c% dlon
    write (6,*)' buftype subtype dbkz = ',s% hd% buf_type,               &
                                          s% hd% buf_subtype, s% hd% dbkz
    write (6,*)' time actual db       = ',chhmm(s% hd% time) ,' ',       &
                                          chhmm(s% actual_time),' ',     &
                                          chhmm(s% hd% db_time)
    write (6,*)' file record          = ',s% hd% source, s% hd% record
    write (6,*)
    write (6,'(a4,4x,a15,3a4)')'name','value','qc','use','src'

    call print (d% total ,'quality control flags')
    call print (d% z_ls  ,'height of land surface     [m]')

    call print (d% z     ,'station height             [m]')
    call print (d% ps    ,'pressure at station height [Pa]')

    call print (d% p_msl ,'pressure mean sea level    [Pa]')
    call print (d% z_msl ,'height   mean sea level    [m]')

    call print (d% p     ,'pressure level             [Pa]')
    call print (d% gp    ,'geopotential               [gpm]')

    call print (d% p_as  ,'pressure     to assimilate [Pa]')
    call print (d% z_as  ,'geopotential to assimilate [gpm]')

    call print (d% p_ref ,'reference pressure         [Pa]')
    call print (d% h_ts  ,'height of temp. sensor     [m]')
    call print (d% t     ,'temperature                [K]')
    call print (d% td    ,'dewpoint temperature       [K]')
    call print (d% rh    ,'relative humidity          [ ]')
    call print (d% h_ws  ,'height of wind sensor      [m]')
    call print (d% ff    ,'wind speed                 [m/s]')
    call print (d% dd    ,'wind direction             [deg]')
    call print (d% uu    ,'wind component             [m/s]')
    call print (d% vv    ,'wind component             [m/s]')

    call print (d% tsea  ,'sea/water temperature      [K]')
    call print (d% h_tsea,'depth of sensor below surf.[m]')
!------------------------------------------------------------------------------
    write (6,'(a,f15.5,9x,a)') 'p_tend  ',d% p_tend,'3 hour pressure change     [Pa]'
    write (6,'(a,f15.5,9x,a)') 'gust1   ',d% gust1 ,'maximum wind speed(gusts)  [m/s] 1 hour'
    write (6,'(a,f15.5,9x,a)') 'gust3   ',d% gust3 ,'maximum wind speed(gusts)  [m/s] 3 hour'
    write (6,'(a,f15.5,9x,a)') 'gust6   ',d% gust6 ,'maximum wind speed(gusts)  [m/s] 6 hour'
    write (6,'(a,f15.5,9x,a)') 'tmax1   ',d% tmax1 ,'maximum temperature        [K]  -1 hour'
    write (6,'(a,f15.5,9x,a)') 'tmin1   ',d% tmin1 ,'minimum temperature        [K]  -1 hour'
    write (6,'(a,f15.5,9x,a)') 'tmax    ',d% tmax  ,'maximum temperature        [K] -12 hour'
    write (6,'(a,f15.5,9x,a)') 'tmin    ',d% tmin  ,'minimum temperature        [K] -12 hour'
    write (6,'(a,i15  ,9x,a)') 'n       ',d% n     ,'cloud cover (total)        [%]'
!   write (6,'(a,i15  ,9x,a)') 'vs      ',d% vs    ,'vertical significance (surface observ.)'
!   write (6,'(a,f15.5,9x,a)') 'nn      ',d% nn    ,'(lowest) cloud amount      [octas]'
    write (6,'(a,f15.5,9x,a)') 'nh      ',d% nh    ,'height of base of cloud    [m]'
    write (6,'(a,f15.5,9x,a)') 'ceiling ',d% ceil  ,'cloud ceiling above ground [m]'
    write (6,'(a,f15.5,9x,a)') 'clcl    ',d% clcl  ,'low    cloud amount        [octas]'
    write (6,'(a,f15.5,9x,a)') 'clcm    ',d% clcm  ,'middle cloud amount        [octas]'
    write (6,'(a,f15.5,9x,a)') 'clch    ',d% clch  ,'high   cloud amount        [octas]'
    write (6,'(a,i15  ,9x,a)') 'ww      ',d% ww    ,'present weather               '
    write (6,'(a,i15  ,9x,a)') 'w       ',d% w     ,'past weather (2) i.e.  6 hour '
    write (6,'(a,f15.5,9x,a)') 'vis     ',d% vis   ,'horizontal visibility      [m]'
    write (6,'(a,i15  ,9x,a)') 'e       ',d% e     ,'state of ground (with or without snow)'
    write (6,'(a,f15.5,9x,a)') 'sdepth  ',d% sdepth,'total snow depth           [m]'
    write (6,'(a,f15.5,9x,a)') 'rr1     ',d% rr1   ,'total precipitation        [kg/m**2]  -1 hour'
    write (6,'(a,f15.5,9x,a)') 'rr3     ',d% rr3   ,'total precipitation        [kg/m**2]  -3 hour'
    write (6,'(a,f15.5,9x,a)') 'rr6     ',d% rr6   ,'total precipitation        [kg/m**2]  -6 hour'
    write (6,'(a,f15.5,9x,a)') 'rr12    ',d% rr12  ,'total precipitation        [kg/m**2] -12 hour'
    write (6,'(a,f15.5,9x,a)') 'rr24    ',d% rr24  ,'total precipitation        [kg/m**2] -24 hour'
    write (6,'(a,f15.5,9x,a)') 'lwr1    ',d% lwr1  ,'long wave radiation        [J/m**2]   -1 hour'
    write (6,'(a,f15.5,9x,a)') 'gsr1    ',d% gsr1  ,'global solar radiation     [J/m**2]   -1 hour'
    write (6,'(a,f15.5,9x,a)') 'gsr3    ',d% gsr3  ,'global solar radiation     [J/m**2]   -3 hour'
    write (6,'(a,f15.5,9x,a)') 'gsr6    ',d% gsr6  ,'global solar radiation     [J/m**2]   -6 hour'
    write (6,'(a,f15.5,9x,a)') 'dsr1    ',d% dsr1  ,'diffuse solar radiation    [J/m**2]   -1 hour'
    write (6,'(a,f15.5,9x,a)') 'dsr3    ',d% dsr3  ,'diffuse solar radiation    [J/m**2]   -3 hour'
    write (6,'(a,f15.5,9x,a)') 'dsr6    ',d% dsr6  ,'diffuse solar radiation    [J/m**2]   -6 hour'
    write (6,'(a,i15  ,9x,a)') 'qobl    ',d% qobl  ,'quality of buoy location      '
    write (6,'(a,i15  ,9x,a)') 'nvs     ',d% nvs   ,'speed of motion            [m/s]'
  end subroutine print_synop
!==============================================================================
  subroutine cmp_synop (s1, s2, diff, m12, m21)
  type (t_synop) ,intent(in)  :: s1, s2
  logical        ,intent(out) :: diff
  logical        ,intent(out) :: m12, m21
    diff = .false.
    m12  = .false.
    m21  = .false.
    call cmp_datum (s1% zr   ,s2% zr   ,diff ,m12 ,m21)
    call cmp_datum (s1% z_ls ,s2% z_ls ,diff ,m12 ,m21)
    call cmp_datum (s1% p    ,s2% p    ,diff ,m12 ,m21)
    call cmp_datum (s1% gp   ,s2% gp   ,diff ,m12 ,m21)
    call cmp_datum (s1% ps   ,s2% ps   ,diff ,m12 ,m21)
    call cmp_datum (s1% p_msl,s2% p_msl,diff ,m12 ,m21)
    call cmp_datum (s1% t    ,s2% t    ,diff ,m12 ,m21)
    call cmp_datum (s1% td   ,s2% td   ,diff ,m12 ,m21)
    call cmp_datum (s1% rh   ,s2% rh   ,diff ,m12 ,m21)
    call cmp_datum (s1% ff   ,s2% ff   ,diff ,m12 ,m21)
    call cmp_datum (s1% dd   ,s2% dd   ,diff ,m12 ,m21)
  end subroutine cmp_synop
!==============================================================================
  subroutine merge_synop (dst, src)
  type (t_synop) ,intent(inout) :: dst
  type (t_synop) ,intent(in)    :: src
    call merge_datum (dst% zr    ,src% zr    )
    call merge_datum (dst% z_ls  ,src% z_ls  )
    call merge_datum (dst% p     ,src% p     )
    call merge_datum (dst% gp    ,src% gp    )
    call merge_datum (dst% ps    ,src% ps    )
    call merge_datum (dst% p_msl ,src% p_msl )
    call merge_datum (dst% t     ,src% t     )
    call merge_datum (dst% td    ,src% td    )
    call merge_datum (dst% rh    ,src% rh    )
    call merge_datum (dst% ff    ,src% ff    )
    call merge_datum (dst% dd    ,src% dd    )
    call merge_datum (dst% h_ts  ,src% h_ts  )
    call merge_datum (dst% h_ws  ,src% h_ws  )
    call merge_datum (dst% tsea  ,src% tsea  )
    call merge_datum (dst% h_tsea,src% h_tsea)
    !--------------------------------------------
    ! Merge simple variables (w/o qc information)
    ! typical for manual correction messages.
    !--------------------------------------------
#define MERGE(x)  if (src% x /= inva_syn% x) dst% x = src% x
    MERGE(gust1)
    MERGE(gust3)
    MERGE(gust6)
    MERGE(tmax1)
    MERGE(tmin1)
    MERGE(tmax)
    MERGE(tmin)
    MERGE(n)
    MERGE(nh)
    MERGE(ceil)
    MERGE(vis)
    MERGE(ww)
!   MERGE(nn)
    MERGE(clcl)
    MERGE(clcm)
    MERGE(clch)
    MERGE(rr1)
    MERGE(rr3)
    MERGE(rr6)
    MERGE(rr12)
    MERGE(rr24)
    MERGE(lwr1)
    MERGE(gsr1)
    MERGE(gsr3)
    MERGE(gsr6)
    MERGE(dsr1)
    MERGE(dsr3)
    MERGE(dsr6)
    MERGE(sdepth)
#undef  MERGE
  end subroutine merge_synop
!==============================================================================
  subroutine count_synop (s1, s2, nsame, ndiff, only1, only2)
    type (t_synop) ,intent(in)  :: s1, s2
    integer        ,intent(out) :: nsame, ndiff
    integer        ,intent(out) :: only1, only2
    !------------------------------------
    ! count useful meteorological content
    !------------------------------------
    nsame = 0
    ndiff = 0
    only1 = 0
    only2 = 0
    call count_datum (s1% zr   ,s2% zr   , nsame, ndiff, only1, only2)
    call count_datum (s1% z_ls ,s2% z_ls , nsame, ndiff, only1, only2)
    call count_datum (s1% p    ,s2% p    , nsame, ndiff, only1, only2)
    call count_datum (s1% gp   ,s2% gp   , nsame, ndiff, only1, only2)
    call count_datum (s1% ps   ,s2% ps   , nsame, ndiff, only1, only2)
    call count_datum (s1% p_msl,s2% p_msl, nsame, ndiff, only1, only2)
    call count_datum (s1% t    ,s2% t    , nsame, ndiff, only1, only2)
    call count_datum (s1% td   ,s2% td   , nsame, ndiff, only1, only2)
    call count_datum (s1% rh   ,s2% rh   , nsame, ndiff, only1, only2)
    call count_datum (s1% ff   ,s2% ff   , nsame, ndiff, only1, only2)
    call count_datum (s1% dd   ,s2% dd   , nsame, ndiff, only1, only2)
    call count_datum (s1% tsea ,s2% tsea , nsame, ndiff, only1, only2)
  end subroutine count_synop
!==============================================================================
  subroutine load_synop (obs, spot, synop)
  type (t_obs)   ,intent(in)  :: obs   ! data of all observations
  type (t_spot)  ,intent(in)  :: spot  ! meta data of this observation
  type (t_synop) ,intent(out) :: synop ! SYNOP information
  !------------------------------------------------------------------------
  ! Load the data from components PAR, OBS of OBS from position provided by
  ! SPOT. Store into SYNOP.
  !------------------------------------------------------------------------
  synop = transfer (obs% par (spot% p% i+1 : spot% p% i + spot% p% n), synop)
  end subroutine load_synop
!------------------------------------------------------------------------------
  subroutine store_synop  (obs, spot, synop)
  type (t_obs)   ,intent(inout) :: obs   ! data of all observations
  type (t_spot)  ,intent(inout) :: spot  ! meta data of this observation
  type (t_synop) ,intent(in)    :: synop ! SYNOP level information
  !-----------------------------------------------------------------------
  ! Store the data from variable SYNOP in the component PAR of
  ! OBS at position provided by SPOT. Allocate memory for PAR if required.
  !-----------------------------------------------------------------------
    integer :: n, par(1)  !+++ work around bug in NEC sxf90/2.0rev360-ia64
    if (synop_int_size==0) synop_int_size = size (transfer (inva_syn,(/0/)))
    n = spot% col% nlev * synop_int_size
    call new_par (obs, n, spot=spot)
    obs % par (spot% p% i+1 : spot% p% i + spot% p% n) = &
      transfer(synop, par)
  end subroutine store_synop
!==============================================================================
  subroutine read_synop_nml
  !-------------------------
  ! read namelist /SYNOP_OBS/
  !-------------------------
    integer :: ierr
    logical :: first = .true.   ! read namelist only once

    if (.not. first) return
    first      = .false.
    !-------------
    ! set defaults
    !-------------
    use_t2     = .false.   ! 2 m temperature
    use_rh     = .false.   ! 2 m humidity
    mon_td     = .true.    ! 2 m dewpoint temperature (monitoring)
    mon_tr     = .true.    ! time range variables (tmin,tmax,gust,rain)
    mon_cl     = .true.    ! cloud amount, etc.
    mon_cwoo   = .false.   ! cloud + weather related, obs only (no operator)
    mon_tsea   = .false.   ! sea temperature monitoring (SHIP,BUOY)
    mon_snow   = .false.   ! monitor snow depth (SYNOP)
    mon_ps     = .false.   ! monitor off-synoptic ps
    mon_t2     = .false.   ! monitor off-synoptic t2m/td2m
    use_vs     = .true.    ! 10 m wind over sea
    use_vl     = .false.   ! 10 m wind over land (tropics, z < 150 m)
    use_vl_gl  = .false.   ! 10 m wind over land (everywhere)
    use_p      = .false.   ! p(h) is observable
    use_h      = .true.    ! h(p) is observable
    ag_rr_1h   = .false.   ! aggregate 1h rain rate
    ag_rr_3h   = .true.    ! aggregate 3h rain rate
    ag_rr_6h   = .true.    ! aggregate 6h rain rate
    ag_rd_1h   = .false.   ! aggregate 1h radiation
    ag_rd_3h   = .true.    ! aggregate 3h radiation
    ag_rd_6h   = .true.    ! aggregate 6h radiation
    ag_gust_3h = .true.    ! aggregate 3h gust
    ag_gust_6h = .true.    ! aggregate 6h gust
    ag_tmaxmin = .true.    ! aggregate tmax/tmin
    ag_combine = .true.    ! combine rr and gust from different BUFRs
    vprc_cloud = 1         ! version of cloud obs pre-processing
    prt_data   = .false.   ! print data
    p_d_msl    =  100      ! max difference between ps and pmsl       [hPa]
    dtdzp      = lapse_cl  ! temperature gradient for p extrapolation [K/m]
    chk_ps     = 0._wp     ! consistency check p_station vs. p_msl
    chk_strict = 0         ! consistency check strictness level
    chk_dbl    = 0         ! extended check for double entries
    chk_dbl_dt =  0._wp    ! tolerance for time difference [min]
    chk_tmaxmin= 0         ! mask for Tmax/Tmin checks
    chk_cldinfo= 0         ! mask for cloud information checks
    chk_ship   = 0         ! extended check for ship reports
    chk_zb_mode= 0         ! consistency check for barometer height
    chk_zb_max =  2._sp    ! barometer height suspicious below 2m
    dzs_zb_max = 20._sp    ! max. diff. |z_station-z_barometer| [m]
    ssdmax_v10m= -1._wp    ! max. SSO standard dev. for using 10m wind [m]
    dz_max_v10m= -1._wp    ! max. height diff. to model surface, v10m [m]
    ssdmax_rh2m= -1._wp    ! max. SSO standard dev. for using 2m humi.[m]
    dz_max_rh2m= -1._wp    ! max. height diff. to model surface, rh2m [m]
    ssdmax_t2m = -1._wp    ! max. SSO standard dev. for using 2m temp.[m]
    dz_max_t2m = -1._wp    ! max. height diff. to model surface, t2m  [m]
    p_land     = (/1,2,3/) ! pressure to use:
    p_ship     = (/3,0,0/) ! 0=none 1=p_s, 2=p_gp, 3=p_msl
    p_buoy     = (/3,0,0/) !
    p_metar    = (/1,3,0/) !
    zls_to_z   = 3         ! set z from land surface height
    zmsl_to_z  = 4         ! set z from mean sea level height
    verbose    = 1         ! Verbosity level of consistency checks
    version    = 0         ! operator version
    v_buoy     = 0         ! operator version for wind, buoys
    buoy_h0_v  = 4.0_wp    ! DRIBU: default height of wind sensor  [m]
    buoy_hr_v  = [ 2, 10 ] ! DRIBU: valid height range of w.sensor [m]
    buoy_hr_t  = [ 2, 10 ] ! DRIBU: valid height range of t.sensor [m]
    hr_t2m     = [ 0, 99 ] ! SYNOP: valid height range of t.sensor [m]
    bc_t2m     = .false.   ! apply b/c for t2m (nonlinear bc)?
    bc_rh2m    = .false.   ! apply b/c for rh2m (nonlinear bc)?
    bc_synop   = .false.   ! apply b/c for SYNOP Land?
    bc_ship    = .false.   ! apply b/c for SYNOP Ship?
    bc_buoy    = .false.   ! apply b/c for DRIBU?
    bc_metar   = .false.   ! apply b/c for METAR?
    chk_zstation=.false.   ! check consistency of reported station height
    z0_bc_synop =  999._wp ! threshold parameter for b/c, SYNOP Land
    z0_bc_ship  =  999._wp ! threshold parameter for b/c, SYNOP Ship
    z0_bc_buoy  =  999._wp ! threshold parameter for b/c, DRIBU
    z0_bc_metar =  999._wp ! threshold parameter for b/c, METAR
    compress_ag = .false.  ! "compress" aggregated data?
    bugfix_t2m  = .true.   ! bugfix for t2m, obs.operator setup (version>=3)
                           ! note: bug was never operational for MEC at DWD
    bugfix_rh2m = .true.   ! bugfix for rh2m, obs.operator setup (version>=3)
                           ! note: bug was never operational for MEC at DWD

    !--------------
    ! read namelist
    !--------------
    if (dace% lpio) then
      write(6,'(a)') repeat('-',79)
      write(6,'()')
      write(6,*) 'Namelist /SYNOP_OBS/:'
      write(6,'()')
      call position_nml ('SYNOP_OBS', status=ierr)
      select case (ierr)
      case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=SYNOP_OBS, iostat=ierr)
        if (ierr/=0) call finish ('read_synop_nml',              &
                                  'ERROR in namelist /SYNOP_OBS/')
#else
        read (nnml ,nml=SYNOP_OBS)
#endif
      case default
        write(6,*) 'Namelist not present, defaults used'
        write(6,'()')
      end select
      !-------------------
      ! consistency checks
      !-------------------
      if (use_vl_gl) use_vl = .true.
      if (biascor_mode <= 0) then
        ag_rr_1h   = .false.; ag_rr_3h   = .false.; ag_rr_6h   = .false.
        ag_rd_1h   = .false.; ag_rd_3h   = .false.; ag_rd_6h   = .false.
        ag_gust_3h = .false.; ag_gust_6h = .false.; ag_tmaxmin = .false.
      endif
      aggregate = (ag_rr_1h   .or. ag_rr_3h   .or. ag_rr_6h .or. &
                   ag_rd_1h   .or. ag_rd_3h   .or. ag_rd_6h .or. &
                   ag_gust_3h .or. ag_gust_6h .or. ag_tmaxmin    )
      chk_dbl_dt = max (chk_dbl_dt, 0._wp)
      !---------
      ! printout
      !---------
      write(6,'(a,l9,a)')    'use_t2       =',use_t2      ,'    use  2 m temperature'
      write(6,'(a,l9,a)')    'use_rh       =',use_rh      ,'    use  2 m humidity'
      write(6,'(a,l9,a)')    'use_vs       =',use_vs      ,'    use 10 m wind over sea'
      write(6,'(a,l9,a)')    'use_vl       =',use_vl      ,'    use 10 m wind over land (TR, z < 150 m)'
      write(6,'(a,l9,a)')    'use_vl_gl    =',use_vl_gl   ,'    use 10 m wind over land (everywhere)'
      write(6,'(a,l9,a)')    'use_p        =',use_p       ,'    p(h) is observable'
      write(6,'(a,l9,a)')    'use_h        =',use_h       ,'    h(p) is observable'
      write(6,'(a,l9,a)')    'mon_td       =',mon_td      ,'    monitor 2 m dewpoint temperature'
      write(6,'(a,l9,a)')    'mon_tr       =',mon_tr      ,'    monitor tmin, tmax, gust, rain'
      write(6,'(a,l9,a)')    'mon_cl       =',mon_cl      ,'    monitor cloud amount, etc.'
      write(6,'(a,l9,a)')    'mon_cwoo     =',mon_cwoo    ,'    monitor obs only cloud + weather var.'
      write(6,'(a,l9,a)')    'mon_tsea     =',mon_tsea    ,'    monitor sea temperature (SHIP,BUOY)'
      write(6,'(a,l9,a)')    'mon_snow     =',mon_snow    ,'    monitor snow depth      (SYNOP)'
      write(6,'(a,l9,a)')    'mon_ps       =',mon_ps      ,'    monitor off-synoptic ps'
      write(6,'(a,l9,a)')    'mon_t2       =',mon_t2      ,'    monitor off-synoptic t2m/td2m'
      write(6,'(a,es9.2,a)') 'dtdzp        =',dtdzp       ,'    temp. gradient for p extrapolation [K/m]'
      write(6,'(a,f9.2,a)')  'p_d_msl      =',p_d_msl     ,'    max difference between ps and pmsl [hPa]'
      write(6,'(a,2f9.2,a)') 'chk_ps       =',chk_ps      ,'    consistency check p_station vs. p_msl'
      write(6,'(a,i9,a)')    'chk_strict   =',chk_strict  ,'    consistency check strictness level (0-3)'
      write(6,'(a,f9.2,a)')  'ssdmax_v10m  =',ssdmax_v10m ,'    max SSO standard dev. for using v10m [m]'
      write(6,'(a,f9.2,a)')  'dz_max_v10m  =',dz_max_v10m ,'    max height diff to model surf., v10m [m]'
      write(6,'(a,f9.2,a)')  'ssdmax_rh2m  =',ssdmax_rh2m ,'    max SSO standard dev. for using rh2m [m]'
      write(6,'(a,f9.2,a)')  'dz_max_rh2m  =',dz_max_rh2m ,'    max height diff to model surf., rh2m [m]'
      write(6,'(a,f9.2,a)')  'ssdmax_t2m   =',ssdmax_t2m  ,'    max SSO standard dev. for using t2m  [m]'
      write(6,'(a,f9.2,a)')  'dz_max_t2m   =',dz_max_t2m  ,'    max height diff to model surf., t2m  [m]'
      write(6,'(a,3i3)')     'p_land       =',pack (p_land ,p_land /=0)
      write(6,'(a,3i3)')     'p_ship       =',pack (p_ship ,p_ship /=0)
      write(6,'(a,3i3)')     'p_buoy       =',pack (p_buoy ,p_buoy /=0)
      write(6,'(a,3i3)')     'p_metar      =',pack (p_metar,p_metar/=0)
      write(6,'(a,i9,a)')    'zls_to_z     =',zls_to_z    ,'    set z from land surface height'
      write(6,'(a,i9,a)')    'zmsl_to_z    =',zmsl_to_z   ,'    set z from mean sea level height'
      write(6,'(a,i9,a)')    'chk_dbl      =',chk_dbl     ,'    extended check for double entries'
      write(6,'(a,f9.2,a)')  'chk_dbl_dt   =',chk_dbl_dt  ,'    tolerance for time difference [min]'
      write(6,'(a,i9,a)')    'chk_tmaxmin  =',chk_tmaxmin ,'    mask for Tmax/Tmin checks (0-7)'
      write(6,'(a,i9,a)')    'chk_cldinfo  =',chk_cldinfo ,'    mask for cloud information checks (0-3)'
      write(6,'(a,i9,a)')    'chk_ship     =',chk_ship    ,'    extended check for ship reports'
      write(6,'(a,i9,a)')    'chk_zb_mode  =',chk_zb_mode ,'    consistency check for barometer height'
      write(6,'(a,f9.2,a)')  'chk_zb_max   =',chk_zb_max  ,'    barometer height suspicious limit [m]'
      write(6,'(a,f9.2,a)')  'dzs_zb_max   =',dzs_zb_max  ,'    max.diff. (z_station-z_barometer) [m]'
      write(6,'(a,l9,a)')    'prt_data     =',prt_data    ,'    print data'
      write(6,'(a,i9,a)')    'verbose      =',verbose     ,'    verbosity'
      write(6,'(a,i9,a)')    'version      =',version     ,'    operator version'
      if (version >= 3) then
        write(6,'(a,l9,a)')  'bugfix_t2m   =',bugfix_t2m  ,'    bugfix for t2m  obs.operator (version>=3)'
        write(6,'(a,l9,a)')  'bugfix_rh2m  =',bugfix_rh2m ,'    bugfix for rh2m obs.operator (version>=3)'
      end if
      write(6,'(a,i9,a)')    'v_buoy       =',v_buoy      ,'    operator version for wind, buoys'
      write(6,'(a,f9.2,a)')  'buoy_h0_v    =',buoy_h0_v   ,'    buoy: default height of wind sensor  [m]'
      write(6,'(a,2f9.2,a)') 'buoy_hr_v    =',buoy_hr_v   ,'    buoy: valid height range of w.sensor [m]'
      write(6,'(a,2f9.2,a)') 'buoy_hr_t    =',buoy_hr_t   ,'    buoy: valid height range of t.sensor [m]'
      write(6,'(a,2f9.2,a)') 'hr_t2m       =',hr_t2m      ,'    land: valid height range of t.sensor [m]'
      write(6,'()')
      write(6,'(a,i9,a)')    'biascor_mode =',biascor_mode,'    bias correction mode'
      write(6,'(a,f9.2,a)')  'reg_param    =',reg_param   ,'    regularization parameter for non-linear bc'
      write(6,'(a,f9.2,a)')  't_decay      =',t_decay     ,'    bc. statistics time decay rate'
      write(6,'(a,i9,a)')    'n_required   =',n_required  ,'    entries required'
      write(6,'(a,l9,a)')    'bc_fallback  =',bc_fallback ,'    create bc.file if missing'
      write(6,'(a,l9,a)')    'ag_rr_1h     =',ag_rr_1h    ,'    aggregate 1h  rain rate observations'
      write(6,'(a,l9,a)')    'ag_rr_3h     =',ag_rr_3h    ,'    aggregate 3h  rain rate observations'
      write(6,'(a,l9,a)')    'ag_rr_6h     =',ag_rr_6h    ,'    aggregate 6h  rain rate observations'
      write(6,'(a,l9,a)')    'ag_rd_1h     =',ag_rd_1h    ,'    aggregate 1h  radiation observations'
      write(6,'(a,l9,a)')    'ag_rd_3h     =',ag_rd_3h    ,'    aggregate 3h  radiation observations'
      write(6,'(a,l9,a)')    'ag_rd_6h     =',ag_rd_6h    ,'    aggregate 6h  radiation observations'
      write(6,'(a,l9,a)')    'ag_gust_3h   =',ag_gust_3h  ,'    aggregate 3h  gust      observations'
      write(6,'(a,l9,a)')    'ag_gust_6h   =',ag_gust_6h  ,'    aggregate 6h  gust      observations'
      write(6,'(a,l9,a)')    'ag_tmaxmin   =',ag_tmaxmin  ,'    aggregate     tmax/tmin observations'
      write(6,'(a,l9,a)')    'ag_combine   =',ag_combine  ,'    combine rr and gust from different BUFRs'
      write(6,'(a,l9,a)')    'aggregate    =',aggregate   ,'    aggregate any           observations'
      write(6,'(a,i9,a)')    'vprc_cloud   =',vprc_cloud  ,'    version of cloud obs pre-processing'
      write(6,'(a,l9,a)')    'bc_t2m       =',bc_t2m      ,'    apply (nonlinear) bias corr. for t2m'
      write(6,'(a,l9,a)')    'bc_rh2m      =',bc_rh2m     ,'    apply (nonlinear) bias corr. for rh2m'
      write(6,'(a,l9,a)')    'bc_synop     =',bc_synop    ,'    apply bias corr. for SYNOP Land'
      write(6,'(a,l9,a)')    'bc_ship      =',bc_ship     ,'    apply bias corr. for SYNOP SHIP'
      write(6,'(a,l9,a)')    'bc_buoy      =',bc_buoy     ,'    apply bias corr. for DRIBU'
      write(6,'(a,l9,a)')    'bc_metar     =',bc_metar    ,'    apply bias corr. for METAR'
      write(6,'(a,l9,a)')    'chk_zstation =',chk_zstation,'    check station height consistency'
      write(6,'(a,f9.2,a)')  'z0_bc_synop  =',z0_bc_synop ,'    bias corr. threshold for SYNOP LAND [m]'
      write(6,'(a,f9.2,a)')  'z0_bc_ship   =',z0_bc_ship  ,'    bias corr. threshold for SYNOP SHIP [m]'
      write(6,'(a,f9.2,a)')  'z0_bc_buoy   =',z0_bc_buoy  ,'    bias corr. threshold for DRIBU      [m]'
      write(6,'(a,f9.2,a)')  'z0_bc_metar  =',z0_bc_metar ,'    bias corr. threshold for METAR      [m]'
      write(6,'()')
      ! Print information if bias correction check is performed for ps
      write(6,'(a,l9,a)') 'do_bccheck_synop_ps =',z0_bc_synop /= 999. ,'    check bc. for PS, SYNOP LAND'
      write(6,'(a,l9,a)') 'do_bccheck_ship_ps  =',z0_bc_ship  /= 999. ,'    check bc. for PS, SYNOP SHIP'
      write(6,'(a,l9,a)') 'do_bccheck_buoy_ps  =',z0_bc_buoy  /= 999. ,'    check bc. for PS, DRIBU'
      write(6,'(a,l9,a)') 'do_bccheck_metar_ps =',z0_bc_metar /= 999. ,'    check bc. for PS, METAR'
      write(6,'()')
    endif
    call p_bcast (use_t2,      dace% pio)
    call p_bcast (mon_td,      dace% pio)
    call p_bcast (mon_tr,      dace% pio)
    call p_bcast (mon_cl,      dace% pio)
    call p_bcast (mon_cwoo,    dace% pio)
    call p_bcast (mon_tsea,    dace% pio)
    call p_bcast (mon_snow,    dace% pio)
    call p_bcast (mon_ps,      dace% pio)
    call p_bcast (mon_t2,      dace% pio)
    call p_bcast (use_rh,      dace% pio)
    call p_bcast (use_vs,      dace% pio)
    call p_bcast (use_vl,      dace% pio)
    call p_bcast (use_vl_gl,   dace% pio)
    call p_bcast (use_p,       dace% pio)
    call p_bcast (use_h,       dace% pio)
    call p_bcast (prt_data,    dace% pio)
    call p_bcast (p_d_msl,     dace% pio)
    call p_bcast (dtdzp,       dace% pio)
    call p_bcast (chk_ps,      dace% pio)
    call p_bcast (chk_strict,  dace% pio)
    call p_bcast (chk_dbl,     dace% pio)
    call p_bcast (chk_dbl_dt,  dace% pio)
    call p_bcast (chk_tmaxmin, dace% pio)
    call p_bcast (chk_cldinfo, dace% pio)
    call p_bcast (chk_ship,    dace% pio)
    call p_bcast (chk_zb_mode, dace% pio)
    call p_bcast (chk_zb_max,  dace% pio)
    call p_bcast (dzs_zb_max,  dace% pio)
    call p_bcast (ssdmax_v10m, dace% pio)
    call p_bcast (dz_max_v10m, dace% pio)
    call p_bcast (ssdmax_rh2m, dace% pio)
    call p_bcast (dz_max_rh2m, dace% pio)
    call p_bcast (ssdmax_t2m,  dace% pio)
    call p_bcast (dz_max_t2m,  dace% pio)
    call p_bcast (p_land,      dace% pio)
    call p_bcast (p_ship,      dace% pio)
    call p_bcast (p_buoy,      dace% pio)
    call p_bcast (p_metar,     dace% pio)
    call p_bcast (zls_to_z,    dace% pio)
    call p_bcast (zmsl_to_z,   dace% pio)
    call p_bcast (verbose,     dace% pio)
    call p_bcast (version,     dace% pio)
    call p_bcast (v_buoy,      dace% pio)
    call p_bcast (buoy_h0_v,   dace% pio)
    call p_bcast (buoy_hr_v,   dace% pio)
    call p_bcast (buoy_hr_t,   dace% pio)
    call p_bcast (hr_t2m,      dace% pio)
    call p_bcast (biascor_mode,dace% pio)
    call p_bcast (reg_param,   dace% pio)
    call p_bcast (t_decay,     dace% pio)
    call p_bcast (n_required,  dace% pio)
    call p_bcast (bc_fallback, dace% pio)
    call p_bcast (ag_rr_1h,    dace% pio)
    call p_bcast (ag_rr_3h,    dace% pio)
    call p_bcast (ag_rr_6h,    dace% pio)
    call p_bcast (ag_rd_1h,    dace% pio)
    call p_bcast (ag_rd_3h,    dace% pio)
    call p_bcast (ag_rd_6h,    dace% pio)
    call p_bcast (ag_gust_3h,  dace% pio)
    call p_bcast (ag_gust_6h,  dace% pio)
    call p_bcast (ag_tmaxmin,  dace% pio)
    call p_bcast (ag_combine,  dace% pio)
    call p_bcast (aggregate,   dace% pio)
    call p_bcast (vprc_cloud,  dace% pio)
    call p_bcast (bc_t2m,      dace% pio)
    call p_bcast (bc_rh2m,     dace% pio)
    call p_bcast (bc_synop,    dace% pio)
    call p_bcast (bc_ship,     dace% pio)
    call p_bcast (bc_buoy ,    dace% pio)
    call p_bcast (bc_metar,    dace% pio)
    call p_bcast (chk_zstation,dace% pio)
    call p_bcast (z0_bc_synop, dace% pio)
    call p_bcast (z0_bc_ship,  dace% pio)
    call p_bcast (z0_bc_buoy,  dace% pio)
    call p_bcast (z0_bc_metar, dace% pio)
    call p_bcast (compress_ag, dace% pio)
    call p_bcast (bugfix_t2m,  dace% pio)
    call p_bcast (bugfix_rh2m, dace% pio)
    p_stat = p_d_msl * 100._wp  ! p_d_msl [Pa]
    z_stat = p_d_msl *   8._wp  ! p_d_msl [gpm]
!   if (.not.dace% lpio) prt_data = .false.
  end subroutine read_synop_nml
!------------------------------------------------------------------------------
  pure function zbs (p,ps,t,zs) ! geopotential height below surface
    !-----------------------------------------------------------
    ! geopotential height below surface (at pressure p)
    ! cf. IFS documentation (CY21R4) , Part II, Paragraph 5.3.2
    !-----------------------------------------------------------
    real(wp) ,intent(in) :: p   ! pressure                   [Pa]
    real(wp) ,intent(in) :: ps  ! model surface pressure     [Pa]
    real(wp) ,intent(in) :: t   ! temperature                [K]
    real(wp) ,intent(in) :: zs  ! model surface geopotential [m^2/s]
    real(wp)             :: zbs ! geopotential

    real(wp) ,parameter :: tx     = 290.5_wp  ! upper temperature bound
    real(wp) ,parameter :: ty     = 255.0_wp  ! lower temperature bound
    real(wp)            :: tstar              ! mean temperature
    real(wp)            :: t0                 ! temperature at msl
    real(wp)            :: gamma              ! adjusted temperature gradient

    tstar   = (t + max (ty, min (tx,t))) / 2._wp
    if (zs > 0._wp .and. p > ps) then
      t0    = tstar + dtdzp * zs / gacc       ! below surface
      t0    = min (t0, max (tx,tstar))
      gamma = R / zs * (t0 - tstar)
    else
      gamma = R * dtdzp / gacc                ! above surface
    endif
    gamma   = max( gamma, 0.01_wp)
    zbs     = zs - R * tstar / gamma * ((p/ps)**gamma-1._wp)
  end function zbs

!==============================================================================
  subroutine read_synop_netcdf (ifile, i_source, obs, rec1, recl, &
                                lkeep, nkeep)
  !====================================================================
  ! Read SYNOP / surface observations from netCDF (converted BUFR data)
  !====================================================================
  integer       ,intent(in)           :: ifile     ! Number of netCDF file read
  integer       ,intent(inout)        :: i_source  ! number of records in source-file
  type (t_obs)  ,intent(inout)        :: obs       ! observations data type to set
  integer       ,intent(in)           :: rec1      ! first record to read
  integer       ,intent(in)           :: recl      ! last  record to read
  logical       ,intent(out)          :: lkeep     ! accept observation ?
  integer       ,intent(out)          :: nkeep     ! number of accepted obsvs.

  type (t_spot)        :: spt0, spti   ! observation meta data
  type (t_spot) ,save  :: empty        !
  type (t_use)         :: use          ! status variable
  type (t_head)        :: head         ! meta information

  integer              :: bufr_type    ! BUFR message type    read
  integer              :: bufr_subtype ! BUFR message subtype read
  integer              :: centre       ! generating centre

  integer              :: obstype      ! observation type
  integer              :: report_subt  ! observation subtype (Datenbankkennz.)
  integer              :: report_subti ! observation subtype index
  type (t_obsid)       :: obsid        ! observation id table entry
  integer              :: stat_time    ! status for time check
  integer              :: idup         ! index of duplicate report

  ! variable for     NetCDF file concerning unlimited dimension, Maximum number of dimensions,
  !                                         number of attributes
  integer              :: status         ! NetCDF status variable
! integer              :: ncdims         ! NetCDF number of dimensions defined in this NetCDF file
! character (LEN=40)   :: yname          ! NetCDF dimension name
  integer              :: ncvars         ! NetCDF number of variables  defined in this NetCDF file
! integer              :: ncatts         ! NetCDF number of attributes defined in this NetCDF file
! integer              :: unlim_dimid    ! NetCDF id for of unlimited dimension defined in this NetCDF file
  integer              :: numDims        ! NetCDF number of dimensions for individual variable in NetCDF variable
  integer              :: numAtts        ! NetCDF number of attributes for individual variable in NetCDF variable
  integer              :: dimid_layer    ! NetCDF dimension id for  Layers
  integer              :: dimid_l_rrr    ! NetCDF dimension id for  precipit.
! integer              :: dimid_report   !                     for  Reports
  integer              :: dimid_l_fxg    !                     for  gusts
! integer              :: dimid_l_nff    !                     for  scat windspeeds

  integer              :: len_layer      ! number of vertical layers in NetCDF-File 1.type of field
  integer              :: len_report     ! number of reports in NetCDF file         1.type of field
  integer              :: len_l_rrr      ! number of precipitations  in NetCDF-File 2.type of field
  integer              :: len_report0    ! number of reports in NetCDF file         2.type of field
  integer              :: len_l_fxg      ! number of gusts           in NetCDF-File 3.type of field
! integer              :: len_l_nff      ! number of scat windspeeds in NetCDF-File 3.type of field
  integer              :: len_report1    ! number of reports in NetCDF file         3.type of field
  integer              :: len_cloud      ! number of vertical layers in NetCDF-File 4.type of field
  integer              :: len_l_rad      ! number of radiances in NetCDF file
  integer              :: j              ! loop index
  integer              :: nc1            ! first  dimension for netcdf 1.type of field getvar in start / count
  integer              :: nc2            ! second dimension for netcdf 1.type of field getvar in start / count
  integer              :: nc3            ! first  dimension for netcdf 2.type of field getvar in start / count
  integer              :: nc4            ! second dimension for netcdf 2.type of field getvar in start / count
  integer              :: nc5            ! first  dimension for netcdf 3.type of field getvar in start / count
  integer              :: nc6            ! second dimension for netcdf 3.type of field getvar in start / count
  integer              :: nc7            ! first  dimension for netcdf 4.type of field getvar in start / count
  integer              :: nc8            ! second dimension for netcdf 4.type of field getvar in start / count
  integer              :: entry1,entry   ! position in source file (subset)

  integer ,allocatable :: ifd_1d  (:)    !
  real    ,allocatable :: rfd_1d  (:)    !
  integer ,allocatable :: ifd_2d1 (:,:)  !
  real    ,allocatable :: rfd_2d1 (:,:)  !
  integer ,allocatable :: ifd_2d2 (:,:)  !
  real    ,allocatable :: rfd_2d2 (:,:)  !
! integer ,allocatable :: ifd_2d3 (:,:)  !
! real    ,allocatable :: rfd_2d3 (:,:)  !
  integer ,allocatable :: ifd_2d4 (:,:)  !
  real    ,allocatable :: rfd_2d4 (:,:)  !

! dummy field
  integer ,allocatable :: ifd_1dm (:)    !
  real    ,allocatable :: rfd_1dm (:)    !
  integer ,allocatable :: ifd_2dm (:,:)  !
  real    ,allocatable :: rfd_2dm (:,:)  !

! real(sp)             :: nhhh_tmp
! integer              :: ivalue

  !----------------
  ! index variables
  !----------------
  integer       :: nreport          ! number of observations (default)

  integer       :: is               ! sub-set    index
! integer       :: ie               ! entry in sub-set
! integer       :: id               ! descriptor index
  integer       :: i                !

  logical       :: lpr_synop        ! extended  printing of synops
  logical       :: lpr_extd         ! extended  printing of synops
  integer       :: npr_extd         ! number of extended printing of aireps
  logical       :: l_buoy           ! report is buoy
  logical       :: l_ship           ! report is synop ship
  logical       :: l_synop          ! report is synop land (but not German national BUFR)
  logical       :: l_nat            ! report is synop land German national BUFR template
  logical       :: l_metar          ! report is metar
  logical       :: l_cman           ! report is cman coastal stations (USA)
  logical       :: l_swis           ! report is road weather station
  logical       :: l_metar_usa      ! report is metar, as synop only for USA
  logical       :: l_seemhem        ! report is seemhem

! logical for meteorological variables(over all reports)
  logical       :: l_mhosnn, l_mhobnn, l_mpn   , l_mppp , l_mpppp,   &
                   l_nhhhn , l_mtdbt , l_mtdnh , l_muuu , l_mhosen,  &
                   l_mvv   , l_mrr24 , l_mggtp1, l_mrrr , l_mn   ,   &
                   l_mvtsu , l_mnh   , l_nh    , l_mcc  , l_mcc0 ,   &
                   l_mcc1  , l_nggtp , l_ndndn , l_nfnfn, l_nggtp0,  &
                   l_nfxgu , l_me    , l_nww   , l_mggtp, l_mhosen5, &
                   l_mw1   , l_mw2   , l_mggtp2, l_mtxtxh,l_mggtp4,  &
                   l_mtntnh, l_nppp  , l_nsss,                       &
                   l_mqobl , l_mphph , l_mggtp7, l_mlwr  ,l_mglsr ,  &
                   l_mdsrh , l_mvtsu0, l_mnh0  , l_nh0   ,l_mhp

! logical for buoys (over all reports)
  logical       :: l_ma1   , l_mwrsa , l_mhawas, l_mhawas1, l_mdrep, l_ndbss, l_mtodb

! logical for metar (over all reports)
  logical       :: l_mdrep2, l_mgwi

! logical for synop sea (buoys, ship)
  logical       :: l_mtn00 , l_nvs

  character(NF90_MAX_NAME)   :: yname_v  ! NetCDF variable name
  character(NF90_MAX_NAME)   :: char_tmp ! NetCDF variable name temporary

  integer            :: i_dft           ! integer default value
  integer, parameter :: qc_dft = -1     ! quality information integer default value
  real               :: r_dft           ! real    default value

  !-----------------------------------------
  ! quantities derived from the BUFR message
  !-----------------------------------------
! integer ,allocatable      :: mds     (:)    ! DIRECTION OF MOTION OF MOV.OBSERV.PLATF.
  integer ,allocatable      :: nvs     (:)    ! SPEED OF MOTION OF MOV.OBSERV.PLATFORM [m/s]
  real    ,allocatable      :: mhosnn  (:)    ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  real    ,allocatable      :: mhp     (:)    ! HEIGHT OF STATION (007001, deprecated)
  real    ,allocatable      :: mhobnn  (:)    ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  real    ,allocatable      :: mpn     (:)    ! PRESSURE (VERT.LOCATION) [Pa]
  real    ,allocatable      :: mppp    (:)    ! PRESSURE [Pa]
  real    ,allocatable      :: mpppp   (:)    ! PRESSURE REDUCED TO MSL [Pa]
  real    ,allocatable      :: mphph   (:)    ! ALTIMETER SETTING (QNH) [Pa]
  real    ,allocatable      :: nppp    (:)    ! 3 HOUR PRESSURE CHANGE [Pa]
! integer ,allocatable      :: na      (:)    ! CHARACTERISTIC OF PRESSURE TENDENCY
! integer ,allocatable      :: np24    (:)    ! 24 HOUR PRESSURE CHANGE [Pa]
  real    ,allocatable      :: nhhhn   (:)    ! GEOPOTENTIAL HEIGHT [gpm]
  real    ,allocatable      :: mhosen  (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  real    ,allocatable      :: mtdbt   (:)    ! TEMPERATURE/DRY BULB TEMPERATURE [K]
  real    ,allocatable      :: mtdnh   (:)    ! DEW-POINT TEMPERATURE            [K]
  real    ,allocatable      :: muuu    (:)    ! RELATIVE HUMIDITY [%]
! real    ,allocatable      :: mhosen0 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  real    ,allocatable      :: mvv     (:)    ! HORIZONTAL VISIBILITY [m]
! real    ,allocatable      :: mhosen1 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  real    ,allocatable      :: mhawas  (:)    ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  real    ,allocatable      :: mhawas1 (:)    ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  real    ,allocatable      :: mrr24   (:)    ! TOTAL PRECIPITATION, PAST 24 HOURS [KG/M**2]
! real    ,allocatable      :: mhosen2 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer ,allocatable      :: mgwi    (:)    ! General weather indicator (TAF/METAR) (020009)
  integer ,allocatable      :: mn      (:)    ! CLOUD COVER (TOTAL) [%]
  integer ,allocatable      :: mvtsu   (:)    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer ,allocatable      :: mnh     (:)    ! CLOUD AMOUNT
  real    ,allocatable      :: nh      (:)    ! HEIGHT OF BASE OF CLOUD [m]
  integer ,allocatable      :: mcc     (:)    ! CLOUD TYPE
  integer ,allocatable      :: mcc0    (:)    ! CLOUD TYPE
  integer ,allocatable      :: mcc1    (:)    ! CLOUD TYPE
  !---------
  ! loop 000
  !---------
  integer ,allocatable      :: mdrep   (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  integer ,allocatable      :: mvtsu0  (:,:)  ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer ,allocatable      :: mnh0    (:,:)  ! CLOUD AMOUNT
! integer ,allocatable      :: mcc2    (:,:)  ! CLOUD TYPE
  real    ,allocatable      :: nh0     (:,:)  ! HEIGHT OF BASE OF CLOUD [m]
  !---------
  ! loop 001
  !---------
! integer ,allocatable      :: mdrep0  (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
! integer ,allocatable      :: mvtsu1  (:,:)  ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! integer ,allocatable      :: mnh1    (:,:)  ! CLOUD AMOUNT
! integer ,allocatable      :: mcc3    (:,:)  ! CLOUD TYPE
! real    ,allocatable      :: nht     (:,:)  ! HEIGHT OF TOP OF CLOUD [m]
! real    ,allocatable      :: mct     (:,:)  ! CLOUD TOP DESCRIPTION
  !---------
  ! loop 002
  !---------
! integer ,allocatable      :: mvtsu2  (:,:)  ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! integer ,allocatable      :: mtdcm   (:,:)  ! TRUE DIR. FROM WHICH CLOUDS ARE MOVING

! integer ,allocatable      :: mvtsu3  (:)    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! real    ,allocatable      :: mda     (:)    ! BEARING OR AZIMUTH [DEGREE_TRUE]
! real    ,allocatable      :: mde     (:)    ! ELEVATION [DEGREE]
! integer ,allocatable      :: mcc4    (:)    ! CLOUD TYPE
! real    ,allocatable      :: mda0    (:)    ! BEARING OR AZIMUTH [DEGREE_TRUE]
! real    ,allocatable      :: mde0    (:)    ! ELEVATION [DEGREE]

  integer ,allocatable      :: me      (:)    ! STATE OF GROUND (WITH OR WITHOUT SNOW)
  real    ,allocatable      :: nsss    (:)    ! TOTAL SNOW DEPTH [m]
! real    ,allocatable      :: mtgtgh  (:)    ! GROUND MINIMUM TEMP., PAST 12 HOURS
! integer ,allocatable      :: mmostm  (:)    ! METHOD OF WATER TEMPERATURE MEASU.
  real    ,allocatable      :: ndbss   (:)    ! DEPTH BELOW SEA/WATER SURFACE [m]
  real    ,allocatable      :: mtn00   (:)    ! SEA/WATER TEMPERATURE [K]
  real    ,allocatable      :: mtn00_  (:,:)  ! SEA/WATER TEMPERATURE [K]         (profiles)
  real    ,allocatable      :: nznzn   (:,:)  ! DEPTH BELOW SEA/WATER SURFACE [m] (profiles)

  integer ,allocatable      :: nww     (:)    ! PRESENT WEATHER
  integer ,allocatable      :: mggtp   (:)    ! TIME PERIOD OR DISPLACEMENT [hour]; buoy: rain
  integer ,allocatable      :: mw1     (:)    ! PAST WEATHER (1)
  integer ,allocatable      :: mw2     (:)    ! PAST WEATHER (2)
  !---------
  ! loop 003
  !---------
! integer ,allocatable      :: mggtp0  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: msssm   (:,:)  ! TOTAL SUNSHINE [minute]
  integer ,allocatable      :: mdrep2  (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !---------
  ! loop 004
  !---------
! real    ,allocatable      :: mhosen3 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer ,allocatable      :: mggtp1  (:,:)  ! TIME PERIOD OR DISPLACEMENT [hour]
  real    ,allocatable      :: mrrr    (:,:)  ! TOTAL PRECIPITATION/TOTAL WATER EQUIV.[kg/m**2]
  !---------
  ! bouy
  !---------
  real    ,allocatable      :: mrrrb   (:)    ! TOTAL PRECIPITATION/TOTAL WATER EQUIV.[kg/m**2]

! real    ,allocatable      :: mhosen4 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer ,allocatable      :: mggtp2  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: mggtp3  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
  real    ,allocatable      :: mtxtxh  (:)    ! MAXIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  integer ,allocatable      :: mggtp4  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: mggtp5  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
  real    ,allocatable      :: mtntnh  (:)    ! MINIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  real    ,allocatable      :: mhosen5 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
! integer ,allocatable      :: niw     (:)    ! TYPE OF INSTRUMENT.FOR WIND MEASUREMENT
! integer ,allocatable      :: mtisi   (:)    ! TIME SIGNIFICANCE
  integer ,allocatable      :: nggtp   (:)    ! TIME PERIOD OR DISPLACEMENT [minute]
  real    ,allocatable      :: ndndn   (:)    ! WIND DIRECTION [degree_true]
  real    ,allocatable      :: ndndn_s (:,:)  ! WIND DIRECTION [degree_true]
  real    ,allocatable      :: nfnfn   (:)    ! WIND SPEED [m/s]
  real    ,allocatable      :: nfnfn_s (:,:)  ! WIND SPEED [m/s]
  !---------
  ! loop 5
  !---------
! integer ,allocatable      :: mtisi0  (:)    ! TIME SIGNIFICANCE
  integer ,allocatable      :: nggtp0  (:,:)  ! TIME PERIOD OR DISPLACEMENT [minute]
! integer ,allocatable      :: nmwgd   (:,:)  ! MAXIMUM WIND GUST DIRECTION [degree_true]
  real    ,allocatable      :: nfxgu   (:,:)  ! MAXIMUM WIND SPEED(GUSTS)   [m/s]
  !---------
  ! buoy
  !---------
  real    ,allocatable      :: nfxgub  (:)    ! MAXIMUM WIND SPEED(GUSTS)   [m/s]

! real    ,allocatable      :: mhosen6 (:)    ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
! integer ,allocatable      :: nggtp1  (:,:)  ! TIME PERIOD OR DISPLACEMENT [minute]
! real    ,allocatable      :: nfxme   (:,:)  ! MAXIMUM WIND SPEED(10 MIN MEAN WIND) [m/s]

! integer ,allocatable      :: mggtp6  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: nie     (:)    ! INSTRUMENT/CROP FOR EVAPO(TRANSPI)RATION
! real    ,allocatable      :: meeev   (:,:)  ! EVAPORATION/EVAPOTRANSPIRATION [KG/M**2]
  !---------
  ! loop 6
  !---------
  integer ,allocatable      :: mggtp7  (:,:)  ! TIME PERIOD OR DISPLACEMENT [hour]
  real    ,allocatable      :: mlwr    (:,:)  ! LONG WAVE RADIAT.,INTEGR./PERIOD SPECIF. [J/M**2]
! real    ,allocatable      :: nriops  (:,:)  ! NET RADIATION,INTEGR.OVER PERIOD SPECIF. [J/M**2]
  real    ,allocatable      :: mglsr   (:,:)  ! GLOBAL SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]
  real    ,allocatable      :: mdsrh   (:,:)  ! DIFFUSE SOLAR RADIAT. (HIGH ACC.) INTEGR [J/M**2]
! real    ,allocatable      :: mdisr   (:,:)  ! DIRECT SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]

! integer ,allocatable      :: mggtp8  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: mggtp9  (:)    ! TIME PERIOD OR DISPLACEMENT [hour]
! integer ,allocatable      :: mtcop   (:)    ! TEMPERATURE CHANGE OVER PERIOD SPECIFIED [K]
! integer ,allocatable      :: mdrep1  (:)    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !------
  ! metar
  !------
  integer ,allocatable      :: mvtsu_  (:,:)  ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer ,allocatable      :: mnh_    (:,:)  ! CLOUD AMOUNT
  real    ,allocatable      :: nh_     (:,:)  ! HEIGHT OF BASE OF CLOUD [m]
  integer ,allocatable      :: mcc_    (:,:)  ! CLOUD TYPE
  !-----
  ! buoy
  !-----
! ! means field present, but no values defined
  integer ,allocatable      :: ma1     (:)    ! WMO REGION NUMBER
  integer ,allocatable      :: mwrsa   (:)    ! WMO REGION SUB-AREA
! integer ,allocatable      :: nix     (:)    ! TYPE OF STATION
! integer ,allocatable      :: nboty   (:)    ! BUOY TYPE
  integer ,allocatable      :: mtodb   (:)    ! TYPE OF DATA BUOY

! integer ,allocatable      :: mdcs    (:)    ! TA COLLECTION AND/OR LOCATION SYSTEM
! real    ,allocatable      :: mdsds   (:)    ! PLATFORM DRIFT SPEED (HIGH PRECISION) [m/s]
! integer ,allocatable      :: mmrvmt  (:)    ! METHOD OF REMOV.VELOCITY+MOTION TABLE...
! integer ,allocatable      :: mqbst   (:)    ! QUALITY OF BUOY SATELLITE TRANSMISSION
  integer ,allocatable      :: mqobl   (:)    ! QUALITY OF BUOY LOCATION
! integer ,allocatable      :: mlqc    (:)    ! LOCATION QUALITY CLASS (..66 CONFIDENCE)
! integer ,allocatable      :: nzd     (:)    ! TOTAL WATER DEPTH [m]
! integer ,allocatable      :: ndw     (:)    ! DIRECTION OF WAVES [DEGREE_TRUE]
! integer ,allocatable      :: mpw     (:)    ! PERIOD OF WAVES [s]
! real    ,allocatable      :: mhw     (:)    ! HEIGHT OF WAVES [m]
! integer ,allocatable      :: ndwdw   (:)    ! DIRECTION OF WIND WAVES [DEGREE_TRUE]
! integer ,allocatable      :: mpwpw   (:)    ! PERIOD OF WIND WAVES [s]
! real    ,allocatable      :: mhwhw   (:)    ! HEIGHT OF WIND WAVES [m]
! integer ,allocatable      :: ndw12   (:)    ! DIRECTION OF SWELL WAVES [DEGREE_TRUE]
! integer ,allocatable      :: mpw12   (:)    ! PERIOD OF SWELL WAVES [s]
! real    ,allocatable      :: mhw12   (:)    ! HEIGHT OF SWELL WAVES [m]
! integer ,allocatable      :: mtyoe   (:)    ! TYPE OF EQUIPMENT
!
! not expanded variables:
! d  defined values
!  * in field list above integrated
!                NBVLR:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE0:long_name = "TYPE OF EQUIPMENT" ;
!                NBVLR0:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE1:long_name = "TYPE OF EQUIPMENT" ;
!                NBVLR1:long_name = "BATTERY VOLTAGE (LARGE RANGE)" ;
!                MTYOE2:long_name = "TYPE OF EQUIPMENT" ;
! d              NDRTY:long_name = "DROGUE TYPE" ;
!                MLDDS:long_name = "LAGRANGIAN DRIFTER DROGUE STATUS" ;
! d              NDRDE:long_name = "DROGUE DEPTH" ;
!                NLDS:long_name = "LAGRANGIAN DRIFTER SUBMERGENCE" ;
! d              MDCI:long_name = "DEPTH CORRECTION INDICATOR" ;
!                NCABL:long_name = "CABLE LENGTH" ;
!                MHPEC:long_name = "HYDROSTATIC PRESS. OF LOWER END OF CABLE" ;
!                MSI:long_name = "ICE DEPOSIT (THICKNESS)" ;
!                MMOSTM:long_name = "METHOD OF SEA-SURFACE TEMPERATURE MEASU." ;
!                NK1:long_name = "INDICATOR FOR DIGITIZATION" ;
! d              NK2:long_name = "METHOD OF SALINITY/DEPTH MEASUREMENT" ;
!                MDREP:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
! d              NZNZN:long_name = "DEPTH BELOW SEA/WATER SURFACE" ;
! d*             MTN00:long_name = "SEA/WATER TEMPERATURE" ;
!                MSNSN:long_name = "SALINITY" ;
!                MMOCM:long_name = "METHOD OF CURRENT MEASUREMENT" ;
!                NK3K4:long_name = "DURATION AND TIME OF CURRENT MEASUREMENT" ;
!
! d              MDREP0:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
!                MHOBNN:long_name = "HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL" ;
!                MTYOE3:long_name = "TYPE OF EQUIPMENT" ;
!                MTINST:long_name = "INSTRUMENT TEMPERATURE" ;
! d*             MPPP:long_name = "PRESSURE" ;
! d*             MPPPP:long_name = "PRESSURE REDUCED TO MSL" ;
! d*             NPPP:long_name = "3 HOUR PRESSURE CHANGE" ;
! d*             NA:long_name = "CHARACTERISTIC OF PRESSURE TENDENCY" ;
!                MTYOE4:long_name = "TYPE OF EQUIPMENT" ;
!                MHOSEN:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MHAWAS:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
! d*             MTDBT:long_name = "TEMPERATURE/DRY BULB TEMPERATURE" ;
!  *             MTDNH:long_name = "DEW-POINT TEMPERATURE" ;
! d*             MUUU:long_name = "RELATIVE HUMIDITY" ;
!                MHOSEN0:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MHAWAS0:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
!                NACH2V:long_name = "ART. CORREC. OF SENSOR H 2 ANOTHER VALUE" ;
! d              MHAWAS1:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
! d              MANTYP:long_name = "ANEMOMETER TYPE" ;
! d              NIW:long_name = "TYPE OF INSTRUMENT.FOR WIND MEASUREMENT" ;
! d*             MTISI1:long_name = "TIME SIGNIFICANCE" ;
! d*             NGGTP:long_name = "TIME PERIOD OR DISPLACEMENT" ;
! d*             NDNDN:long_name = "WIND DIRECTION" ;
! d*             NFNFN:long_name = "WIND SPEED" ;
!
!                MTISI2:long_name = "TIME SIGNIFICANCE" ;
!                NGGTP0:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             NMWGD:long_name = "MAXIMUM WIND GUST DIRECTION" ;
!  *             NFXGU:long_name = "MAXIMUM WIND SPEED(GUSTS)" ;
!                NACH2V0:long_name = "ART. CORREC. OF SENSOR H 2 ANOTHER VALUE" ;
!                MHAWAS2:long_name = "HEIGHT OF SENSOR ABOVE WATER SURFACE" ;
!                MHOSEN1:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
!                MGGTP:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             MRRR:long_name = "TOTAL PRECIPITATION/TOTAL WATER EQUIV." ;
!                MHOSEN2:long_name = "HEIGHT OF SENSOR ABOVE LOCAL GROUND" ;
! d              MTISI3:long_name = "TIME SIGNIFICANCE" ;
!                MGGTP0:long_name = "TIME PERIOD OR DISPLACEMENT" ;
!  *             MSSSS:long_name = "GLOBAL RADIATION INTEGR.O.PERIOD SPECIF." ;
!                MTISI4:long_name = "TIME SIGNIFICANCE" ;
! d              NOMDP:long_name = "OPERATOR OR MANUFACTURER DEF. PARAMETER" ;
!                MDREP1:long_name = "DELAYED DESCRIPTOR REPLICATION FACTOR" ;
!                YSUPL:long_name = "008 CHARACTERS" ;
!                DATE:long_name = "Date as YYYYMMDD" ;
!                TIME:long_name = "Time as HHMMSS" ;
!
!                NORBNU:long_name = "ORBIT NUMBER" ;
!                NCTCN:long_name = "CROSS-TRACK CELL NUMBER" ;
!                NWVCQ:long_name = "ASCAT WIND VECTOR CELL QUALITY" ;
!                MSWCQ:long_name = "SEAWIND WIND VECTOR CELL QUALITY" ;
!                MPROR:long_name = "PROBABILITY OF RAIN" ;
!                NNOVA:long_name = "NUMBER OF VECTOR AMBIGUITIES" ;
!                NISWV:long_name = "INDEX OF SELECTED WIND VECTOR" ;
!                NTNSM:long_name = "TOTAL NUMBER OF SIGMA-0 MEASUREMENTS" ;
!                NFF:long_name = "WIND SPEED AT 10 M" ;
!                NDD:long_name = "WIND DIRECTION AT 10 M" ;
!                MLCFS:long_name = "LIKELIHOOD COMPUTED FOR SOLUTION" ;

! integer ,allocatable      :: lmtodb   (:)    ! TYPE OF DATA BUOY
! integer ,allocatable      :: lmqbst   (:)    ! QUALITY OF BUOY SATELLITE TRANSMISSION
! integer ,allocatable      :: lmqobl   (:)    ! QUALITY OF BUOY LOCATION
! integer ,allocatable      :: lmlqc    (:)    ! LOCATION QUALITY CLASS (..66 CONFIDENCE)

  !---------------------------------------
  ! NetCDF variable IDs for body   section
  !---------------------------------------
  ! variable ID's in NetCDF file for
  !    expansion  in NetCDF file BUFR- data section4
  !
! integer    :: varid_mds      ! DIRECTION OF MOTION OF MOV.OBSERV.PLATF.
  integer    :: varid_nvs      ! SPEED OF MOTION OF MOV.OBSERV.PLATFORM [m/s]
  integer    :: varid_mhosnn   ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  integer    :: varid_mhp      ! HEIGHT OF STATION
  integer    :: varid_mhobnn   ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  integer    :: varid_mpn      ! PRESSURE (VERT.LOCATION) [Pa]
  integer    :: varid_mppp     ! Q-BITS FOR PRESSURE [Pa]
  integer    :: varid_mpppp    ! PRESSURE REDUCED TO MSL [Pa]
  integer    :: varid_mphph    ! ALTIMETER SETTING (QNH) [Pa]
  integer    :: varid_nppp     ! 3 HOUR PRESSURE CHANGE [Pa]
! integer    :: varid_na       ! CHARACTERISTIC OF PRESSURE TENDENCY
! integer    :: varid_np24     ! 24 HOUR PRESSURE CHANGE [Pa]
  integer    :: varid_nhhhn    ! GEOPOTENTIAL HEIGHT [gpm]
  integer    :: varid_mhosen   ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mtdbt    ! TEMPERATURE/DRY BULB TEMPERATURE [K]
  integer    :: varid_mtdnh    ! DEW-POINT TEMPERATURE            [K]
  integer    :: varid_muuu     ! RELATIVE HUMIDITY [%]
! integer    :: varid_mhosen0  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mvv      ! HORIZONTAL VISIBILITY [m]
! integer    :: varid_mhosen1  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mhawas   ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  integer    :: varid_mhawas1  ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  integer    :: varid_mrr24    ! TOTAL PRECIPITATION, PAST 24 HOURS [KG/M**2]
! integer    :: varid_mhosen2  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mgwi     ! General weather indicator (TAF/METAR) (020009)
  integer    :: varid_mn       ! CLOUD COVER (TOTAL) [%]
  integer    :: varid_mvtsu    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer    :: varid_mnh      ! CLOUD AMOUNT
  integer    :: varid_nh       ! HEIGHT OF BASE OF CLOUD [m]
  integer    :: varid_mcc      ! CLOUD TYPE
  integer    :: varid_mcc0     ! CLOUD TYPE
  integer    :: varid_mcc1     ! CLOUD TYPE
  !---------
  ! loop 000
  !---------
  integer    :: varid_mdrep    ! DELAYED DESCRIPTOR REPLICATION FACTOR
  integer    :: varid_mvtsu0   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  integer    :: varid_mnh0     ! CLOUD AMOUNT
! integer    :: varid_mcc2     ! CLOUD TYPE
  integer    :: varid_nh0      ! HEIGHT OF BASE OF CLOUD [m]
  !---------
  ! loop 001
  !---------
! integer    :: varid_mdrep0   ! DELAYED DESCRIPTOR REPLICATION FACTOR
! integer    :: varid_mvtsu1   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! integer    :: varid_mnh1     ! CLOUD AMOUNT
! integer    :: varid_mcc3     ! CLOUD TYPE
! integer    :: varid_nht      ! HEIGHT OF TOP OF CLOUD [m]
! integer    :: varid_mct      ! CLOUD TOP DESCRIPTION
  !---------
  ! loop 002
  !---------
! integer    :: varid_mvtsu2   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! integer    :: varid_mtdcm    ! TRUE DIR. FROM WHICH CLOUDS ARE MOVING

! integer    :: varid_mvtsu3   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
! integer    :: varid_mda      ! BEARING OR AZIMUTH [DEGREE_TRUE]
! integer    :: varid_mde      ! ELEVATION [DEGREE]
! integer    :: varid_mcc4     ! CLOUD TYPE
! integer    :: varid_mda0     ! BEARING OR AZIMUTH [DEGREE_TRUE]
! integer    :: varid_mde0     ! ELEVATION [DEGREE]

  integer    :: varid_me       ! STATE OF GROUND (WITH OR WITHOUT SNOW)
  integer    :: varid_nsss     ! TOTAL SNOW DEPTH [m]
! integer    :: varid_mtgtgh   ! GROUND MINIMUM TEMP., PAST 12 HOURS
! integer    :: varid_mmostm   ! METHOD OF WATER TEMPERATURE MEASU.
  integer    :: varid_ndbss    ! DEPTH BELOW SEA/WATER SURFACE [m]
  integer    :: varid_mtn00    ! SEA/WATER TEMPERATURE [K]

  integer    :: varid_nww      ! PRESENT WEATHER
  integer    :: varid_mggtp    ! TIME PERIOD OR DISPLACEMENT [hour]
  integer    :: varid_mw1      ! PAST WEATHER (1)
  integer    :: varid_mw2      ! PAST WEATHER (2)
  !---------
  ! loop 003
  !---------
! integer    :: varid_mggtp0   ! TIME PERIOD OR DISPLACEMENT [hour]
! integer    :: varid_msssm    ! TOTAL SUNSHINE [minute]
  integer    :: varid_mdrep2   ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !---------
  ! loop 004
  !---------
! integer    :: varid_mhosen3  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mggtp1   ! TIME PERIOD OR DISPLACEMENT [hour]; ship varid_mggtp0
  integer    :: varid_mrrr     ! TOTAL PRECIPITATION/TOTAL WATER EQUIV.[kg/m**2]

! integer    :: varid_mhosen4  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  integer    :: varid_mggtp2   ! TIME PERIOD OR DISPLACEMENT [hour]; ship varid_mggtp1
! integer    :: varid_mggtp3   ! TIME PERIOD OR DISPLACEMENT [hour]; ship varid_mggtp2
  integer    :: varid_mtxtxh   ! MAXIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  integer    :: varid_mggtp4   ! TIME PERIOD OR DISPLACEMENT [hour]; ship varid_mggtp3
! integer    :: varid_mggtp5   ! TIME PERIOD OR DISPLACEMENT [hour]; ship varid_mggtp4
  integer    :: varid_mtntnh   ! MINIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  integer    :: varid_mhosen5  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
! integer    :: varid_niw      ! TYPE OF INSTRUMENT.FOR WIND MEASUREMENT
! integer    :: varid_mtisi    ! TIME SIGNIFICANCE
  integer    :: varid_nggtp    ! TIME PERIOD OR DISPLACEMENT [minute]
  integer    :: varid_ndndn    ! WIND DIRECTION [degree_true]
  integer    :: varid_nfnfn    ! WIND SPEED [m/s]
  !-------
  ! loop 5
  !-------
! integer    :: varid_mtisi0   ! TIME SIGNIFICANCE
  integer    :: varid_nggtp0   ! TIME PERIOD OR DISPLACEMENT [minute]
! integer    :: varid_nmwgd    ! MAXIMUM WIND GUST DIRECTION [degree_true]
  integer    :: varid_nfxgu    ! MAXIMUM WIND SPEED(GUSTS)   [m/s]

! integer    :: varid_mhosen6  ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
! integer    :: varid_nggtp1   ! TIME PERIOD OR DISPLACEMENT [minute]
! integer    :: varid_nfxme    ! MAXIMUM WIND SPEED(10 MIN MEAN WIND) [m/s]

! integer    :: varid_mggtp6   ! TIME PERIOD OR DISPLACEMENT [hour]
! integer    :: varid_nie      ! INSTRUMENT/CROP FOR EVAPO(TRANSPI)RATION
! integer    :: varid_meeev    ! EVAPORATION/EVAPOTRANSPIRATION [KG/M**2]
  !-------
  ! loop 6
  !-------
  integer    :: varid_mggtp7   ! TIME PERIOD OR DISPLACEMENT [hour]
  integer    :: varid_mlwr     ! LONG WAVE RADIAT.,INTEGR./PERIOD SPECIF. [J/M**2]
! integer    :: varid_nriops   ! NET RADIATION,INTEGR.OVER PERIOD SPECIF. [J/M**2]
  integer    :: varid_mglsr    ! GLOBAL SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]
  integer    :: varid_mdsrh    ! DIFFUSE SOLAR RADIAT. (HIGH ACC.) INTEGR [J/M**2]
! integer    :: varid_mdisr    ! DIRECT SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]

! integer    :: varid_mggtp8   ! TIME PERIOD OR DISPLACEMENT [hour]
! integer    :: varid_mggtp9   ! TIME PERIOD OR DISPLACEMENT [hour]
! integer    :: varid_mtcop    ! TEMPERATURE CHANGE OVER PERIOD SPECIFIED [K]
! integer    :: varid_mdrep1   ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !-----
  ! buoy
  !-----
  ! means field present, but no values defined
  integer    :: varid_ma1      ! WMO REGION NUMBER
  integer    :: varid_mwrsa    ! WMO REGION SUB-AREA
! integer    :: varid_nix      ! TYPE OF STATION
! integer    :: varid_nboty    ! BUOY TYPE
  integer    :: varid_mtodb    ! TYPE OF DATA BUOY

! integer    :: varid_mdcs     ! TA COLLECTION AND/OR LOCATION SYSTEM
! integer    :: varid_mdsds    ! PLATFORM DRIFT SPEED (HIGH PRECISION) [m/s]
! integer    :: varid_mmrvmt   ! METHOD OF REMOV.VELOCITY+MOTION TABLE...
! integer    :: varid_mqbst    ! QUALITY OF BUOY SATELLITE TRANSMISSION
  integer    :: varid_mqobl    ! QUALITY OF BUOY LOCATION
! integer    :: varid_mlqc     ! LOCATION QUALITY CLASS (..66 CONFIDENCE)
! integer    :: varid_nzd      ! TOTAL WATER DEPTH [m]
! integer    :: varid_ndw      ! DIRECTION OF WAVES [DEGREE_TRUE]
! integer    :: varid_mpw      ! PERIOD OF WAVES [s]
! integer    :: varid_mhw      ! HEIGHT OF WAVES [m]

    !-----------------------------------------
    ! quantities derived from the BUFR message
    !-----------------------------------------
!   integer            :: icount! Zaehler
!   character(len=8)   :: ymnem ! value decoded from BUFR (string)
    logical            :: diff, m12, m21
    type (t_synop)     :: s1

    type(t_synop)      :: s                ! SYNOP observation read
!   integer            :: yyyy,mo,dd,hh,mi ! actual time read
    character(len=32)  :: com
!   integer            :: ip
    real(wp)           :: dt1, dt2         ! difference to analysis time
    real(wp)           :: dh1, dh2         ! difference to nearest hour
    integer            :: nsame, ndiff     ! data content counters
    integer            :: only1, only2
#if defined (__SX__)
    ! Misc. local variables
    integer              :: n_idx, k, n
    integer, allocatable :: idx(:)         ! Auxiliary arrays, used in
    logical, allocatable :: mask(:)        ! checks for double entries
#endif
  !------------
  ! temporaries
  !------------
  real(wp) :: rhour
  integer  :: ihour
  logical  :: same_time                    ! Same (hourly) time slot
  integer  :: pr1, pr2                     ! priorities
  logical  :: chk
  logical  :: hourly1, hourly2, prefer1
  logical  :: buoy1_as_ship, buoy2_as_ship ! double recoding of buoy via ship?
  real(wp) :: slat, slon                   ! scaled values of lat, lon
  integer  :: len, id
  integer  :: nl, lm1                      ! cloud layers: total, -1
  integer  :: mcc_lmh(3)                   ! cloud type: low,middle,high
  integer  :: dimids (dimids_max)
  integer, allocatable :: mvtsu_tmp(:)
  integer, allocatable :: mnh_tmp  (:)
! integer, allocatable :: mcc_tmp  (:)
  real,    allocatable :: nh_tmp   (:)
!------------------------------------------------------------------------------

    l_mrrr  = .false.
    l_nfnfn = .false.
    l_nfxgu = .false.
    l_mggtp = .false.

    lpr_synop = .false.; if (netcdf_verb > 1) lpr_synop = .true.
    lpr_extd   = .true.
    npr_extd   =   2
!   icount     =   0
!   if( lpr_extd) npr_extd  =   10
    if( lpr_extd) npr_extd  =  100
!   if( lpr_extd) npr_extd  =  300
!   if( lpr_extd) npr_extd  = 1000
!   if( lpr_extd) npr_extd  = 5000

  !------------------------------
  ! get default number of reports
  !------------------------------
  nreport     =  recl - rec1 + 1

  lkeep       = .false.
  nkeep       = 0

  report_subt  = s2ikz(1)
  l_buoy       = any (report_subt == kz_buoy)
  l_ship       = any (report_subt == kz_ship)
  l_synop      = any (report_subt == kz_synop)
  l_nat        = any (report_subt == kz_nat)
  l_metar      = any (report_subt == kz_metar)
  l_cman       = any (report_subt == kz_cman)
  l_swis       = any (report_subt == kz_swis)
  l_metar_usa  = any (report_subt == kz_metar_usa)
  l_seemhem    = any (report_subt == kz_seemhem)
  !---------------------------
  ! get variables ids and name
  !---------------------------
  status = nf90_Inquire (ncid, nVariables=ncvars)
  if ( lpr_synop .and. lpr_extd ) then
    do j = 1 , ncvars
    status = nf90_Inquire_Variable(ncid, j,name=yname_v )
    write (6,'(a,i3,a,a)') 'nf90_Inquire_Variable(', j,          &
                           ')  Variable name(o): ', trim (yname_v)
    end do
  endif
  !========================
  ! get dimension of fields
  !========================
! set some defaults
  len_layer  = 0
  len_report = nreport
  nc1 = 0
  nc2 = 0

  if (l_buoy .or. l_cman) then
! get dimension of water temperature profile
    status = nf90_inq_varid (ncid, 'MTN00',  varid_MTN00)
    if (status == nf90_noerr) then
      status = nf90_Inquire_Variable(ncid, varid_MTN00, ndims=numDims, dimids=dimids, natts=numAtts)
      if(numDims >= 2) then
        dimid_layer = dimids(1)
        status = nf90_Inquire_Dimension(ncid, dimid_layer, len=len_layer)
        nc1 = len_layer
        nc2 = len_report
      endif
    end if
  else
! get dimension of cloud loop
    status = nf90_inq_varid (ncid, 'MNH0' ,  varid_MNH0)
    if (status /= nf90_noerr) then
      nc1       = len_layer
      nc2       = len_report
    else
      status = nf90_Inquire_Variable(ncid, varid_MNH0, ndims=numDims, dimids=dimids, natts=numAtts)
      if(numDims >= 2) then
        dimid_layer  = dimids(1)
!       dimid_report = dimids(2)
        status = nf90_Inquire_Dimension(ncid, dimid_layer, len=len_layer)
        nc1    = len_layer
        nc2    = len_report
      endif
    endif
  end if
  !----------------------
  ! default for buoy case
  !----------------------
  nc3 = 1
  nc4 = 1
  len_l_rrr   = 1
  len_report0 = nreport
  if ( .not.l_buoy) then
  ! get dimension of precipitation loop
    status = nf90_inq_varid (ncid, 'MRRR' ,  varid_MRRR)
    if (status /= nf90_noerr) then
      nc3       = len_l_rrr
      nc4       = len_report0
    else
      status = nf90_Inquire_Variable(ncid, varid_MRRR, ndims=numDims, dimids=dimids, natts=numAtts)

      nc3 = 0
      nc4 = 0

      if(numDims >= 2) then
        dimid_l_rrr  = dimids(1)
!       dimid_report = dimids(2)
        status = nf90_Inquire_Dimension(ncid, dimid_l_rrr,  len=len_l_rrr)
!       status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report0)
        nc3       = len_l_rrr
        nc4       = len_report0
      endif
    endif
  endif

  !----------------------
  ! default for buoy case
  !----------------------
  nc5       = 1
  nc6       = 1
  len_l_fxg   = 1
  len_report1 = nreport
  if ( .not.l_buoy .and. .not. l_metar) then
    ! get dimension of gust loop
    status = nf90_inq_varid (ncid, 'NFXGU' ,  varid_NFXGU)
    if (status /= nf90_noerr) then
      nc5         = len_l_fxg
      nc6         = len_report1
    else
      status = nf90_Inquire_Variable(ncid, varid_NFXGU, ndims=numDims, dimids=dimids, natts=numAtts)

      nc5 = 0
      nc6 = 0

      if(numDims >= 2) then
        dimid_l_fxg  = dimids(1)
!       dimid_report = dimids(2)
        status = nf90_Inquire_Dimension(ncid, dimid_l_fxg,  len=len_l_fxg )
!       status = nf90_Inquire_Dimension(ncid, dimid_report, len=len_report1)
        nc5       = len_l_fxg
        nc6       = len_report1
      endif
    endif
  endif

  nc7 = 0
  nc8 = 1
  if (l_metar) then
    ! get dimension of cloud loop
    status = nf90_inq_varid (ncid, 'MVTSU' ,  varid_mvtsu)
    if (status == nf90_noerr) then
      status = nf90_Inquire_Variable(ncid, varid_mvtsu, ndims=numDims, dimids=dimids, natts=numAtts)
      if(numDims == 2) then
        dimid_l_fxg  = dimids(1)
!       dimid_report = dimids(2)
        status = nf90_Inquire_Dimension(ncid, dimid_l_fxg,  len=len_cloud)
        nc7       = len_cloud
        nc8       = len_report1
      endif
    endif
    if (netcdf_verb > 0) WRITE(6,*) 'NC7/NC8', nc7, nc8
  end if
  len_cloud = nc7

  !===================================
  ! allocate fields for reading netcdf
  !===================================

  allocate ( ifd_1d  (len_report) )
  allocate ( rfd_1d  (len_report) )
  allocate ( ifd_2d1 (nc1,nc2   ) )
  allocate ( rfd_2d1 (nc1,nc2   ) )
  allocate ( ifd_2d2 (nc3,nc4   ) )
  allocate ( rfd_2d2 (nc3,nc4   ) )
! allocate ( ifd_2d3 (nc5,nc6   ) )
! allocate ( rfd_2d3 (nc5,nc6   ) )
  if (l_metar) then
   allocate (ifd_2d4 (nc7,nc8   ) )
   allocate (rfd_2d4 (nc7,nc8   ) )
  end if

  allocate ( ifd_1dm (  1       ) )
  allocate ( rfd_1dm (  1       ) )
  allocate ( ifd_2dm (  1,  1   ) )
  allocate ( rfd_2dm (  1,  1   ) )

  if ( .not.l_buoy .and. lpr_synop ) then
    write (6,'()')
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' BEGINN  MO_SYNOP        !!!!!!!!!!!!!!!!!!!!!!!'
    write (6,'(a,i3,a)' )   'pe=',dace% pe, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    write (6,'(a,i4,/,6x,a,i4,/,6x,a,i4,/,6x,a,10i4)')   'pe=',dace% pe,            &
                            ' varid_MRRR  number of dimensions: ',numDims,          &
                            ' varid_MRRR  number of attributes: ',numAtts,          &
                            ' varid_MRRR  ids    of dimensions: ',dimids(1:numDims)
    write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10) ') &
                            'pe=',dace% pe,  ' Reports in          BUFR File:',len_report ,&
                                         ' Number cloud layer  BUFR File:',len_layer  ,    &
                                         ' nc1=',nc1,' nc2=',nc2

   write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10) ') &
                            'pe=',dace% pe,  ' Reports in          BUFR File:',len_report ,&
                                         ' Rain groups         BUFR File:',len_l_rrr  ,    &
                                         ' nc3=',nc3,' nc4=',nc4

    write (6,'(a,i3,a,i8,/ , &
             & 1x  ,a,i8,/ , &
             & 1x  ,a,i10,a,i10) ') &
                            'pe=',dace% pe,  ' Reports in      BUFR File:',len_report0,&
                                         ' Periods fxgu in BUFR File:',len_l_fxg  ,    &
                                         ' nc5=',nc5,' nc6=',nc6
  endif
  !----------------------------------------
  ! define number of reports in source-file
  !----------------------------------------
  i_source = len_report

  !----------------
  ! allocate fields
  !----------------
  allocate ( ma1     (len_report))     ! WMO REGION NUMBER
  allocate ( mwrsa   (len_report))     ! WMO REGION SUB-AREA
  allocate ( mhosnn  (len_report))     ! HEIGHT OF STATION GROUND ABOVE MEAN SEA
  allocate ( mhobnn  (len_report))     ! HEIGHT OF BAROMETER ABOVE MEAN SEA LEVEL
  allocate ( mhp     (len_report))     ! HEIGHT OF STATION
  allocate ( mpn     (len_report))     ! PRESSURE (VERT.LOCATION) [Pa]
  allocate ( mppp    (len_report))     ! PRESSURE [Pa]
  allocate ( mpppp   (len_report))     ! PRESSURE REDUCED TO MSL [Pa]
  allocate ( mphph   (len_report))     ! ALTIMETER SETTING (QNH) [Pa]
  allocate ( nppp    (len_report))     ! 3 HOUR PRESSURE CHANGE [Pa]
! allocate ( np24    (len_report))     ! 24 HOUR PRESSURE CHANGE [Pa]
  allocate ( nhhhn   (len_report))     ! GEOPOTENTIAL HEIGHT [gpm]
  allocate ( mhosen  (len_report))     ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  allocate ( mhawas  (len_report))     ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  allocate ( mhawas1 (len_report))     ! HEIGHT OF SENSOR ABOVE WATER SURFACE [m]
  allocate ( mtdbt   (len_report))     ! TEMPERATURE/DRY BULB TEMPERATURE [K]
  allocate ( mtdnh   (len_report))     ! DEW-POINT TEMPERATURE            [K]
  allocate ( muuu    (len_report))     ! RELATIVE HUMIDITY [%]
  allocate ( mvv     (len_report))     ! HORIZONTAL VISIBILITY [m]
  allocate ( mrr24   (len_report))     ! TOTAL PRECIPITATION, PAST 24 HOURS [KG/M**2]
  allocate ( mgwi    (len_report))     ! General weather indicator (TAF/METAR)
  allocate ( mn      (len_report))     ! CLOUD COVER (TOTAL) [%]
  allocate ( mvtsu   (len_report))     ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  allocate ( mnh     (len_report))     ! CLOUD AMOUNT
  allocate ( nh      (len_report))     ! HEIGHT OF BASE OF CLOUD [m]
  allocate ( mcc     (len_report))     ! CLOUD TYPE
  allocate ( mcc0    (len_report))     ! CLOUD TYPE
  allocate ( mcc1    (len_report))     ! CLOUD TYPE
  allocate ( nvs     (len_report))     ! SPEED OF MOTION OF MOV.OBSERV.PLATFORM [m/s]
  !---------
  ! loop 000
  !---------
  allocate ( mdrep   (len_report))     ! DELAYED DESCRIPTOR REPLICATION FACTOR
  allocate ( mvtsu0  (len_layer,len_report))   ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
  allocate ( mnh0    (len_layer,len_report))   ! CLOUD AMOUNT
! allocate ( mcc2    (len_layer,len_report))   ! CLOUD TYPE
  allocate ( nh0     (len_layer,len_report))   ! HEIGHT OF BASE OF CLOUD [m]
  !---------
  ! loop 001
  ! loop 002
  !---------
  allocate ( me      (len_report))     ! STATE OF GROUND (WITH OR WITHOUT SNOW)
  allocate ( nsss    (len_report))     ! TOTAL SNOW DEPTH [m]
  allocate ( mtn00   (len_report))     ! SEA/WATER TEMPERATURE [K]
  allocate ( ndbss   (len_report))     ! DEPTH BELOW SEA/WATER SURFACE [m]

  if (l_buoy) then
  allocate ( mtn00_  (len_layer,len_report)) ! SEA/WATER TEMPERATURE [K] (profiles)
  allocate ( nznzn   (len_layer,len_report)) ! DEPTH BELOW SEA/WATER SURFACE [m]
  end if

  if (l_cman) then
  allocate ( mtn00_  (len_layer,len_report)) ! SEA/WATER TEMPERATURE [K] (profiles)
! allocate ( ndbss_  (len_layer,len_report)) ! DEPTH BELOW SEA/WATER SURFACE [m]
  end if

  allocate ( nww     (len_report))     ! PRESENT WEATHER
  allocate ( mggtp   (len_report))     ! TIME PERIOD OR DISPLACEMENT [hour]
  allocate ( mw1     (len_report))     ! PAST WEATHER (1)
  allocate ( mw2     (len_report))     ! PAST WEATHER (2)

  !---------
  ! loop 003
  ! loop 004
  !---------
! allocate ( mhosen3 (len_report))     ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  allocate ( mggtp1  (len_l_rrr,len_report))   ! TIME PERIOD OR DISPLACEMENT [hour]
  allocate ( mrrr    (len_l_rrr,len_report))   ! TOTAL PRECIPITATION/TOTAL WATER EQUIV.[kg/m**2]
! for buoy
  allocate ( mrrrb   (len_report))     ! TOTAL PRECIPITATION/TOTAL WATER EQUIV.[kg/m**2]

! allocate ( mhosen4 (len_report))     ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
  allocate ( mggtp2  (len_report))     ! TIME PERIOD OR DISPLACEMENT [hour]
! allocate ( mggtp3  (len_report))     ! TIME PERIOD OR DISPLACEMENT [hour]
  allocate ( mtxtxh  (len_report))     ! MAXIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  allocate ( mggtp4  (len_report))     ! TIME PERIOD OR DISPLACEMENT [hour]
! allocate ( mggtp5  (len_report))     ! TIME PERIOD OR DISPLACEMENT [hour]
  allocate ( mtntnh  (len_report))     ! MINIMUM TEMP., HEIGHT/PERIOD SPECIFIED [K]

  allocate ( mhosen5 (len_report))     ! HEIGHT OF SENSOR ABOVE LOCAL GROUND [m]
! allocate ( niw     (len_report))     ! TYPE OF INSTRUMENT.FOR WIND MEASUREMENT
! allocate ( mtisi   (len_report))     ! TIME SIGNIFICANCE
  allocate ( nggtp   (len_report))     ! TIME PERIOD OR DISPLACEMENT [minute]
  allocate ( ndndn   (len_report))     ! WIND DIRECTION [degree_true]
  allocate ( ndndn_s (len_l_fxg,len_report))     ! WIND DIRECTION [degree_true]
  allocate ( nfnfn   (len_report))     ! WIND SPEED [m/s]
  allocate ( nfnfn_s (len_l_fxg,len_report))     ! WIND SPEED [m/s]
  allocate ( mdrep2  (len_report))     ! DELAYED DESCRIPTOR REPLICATION FACTOR
  !-------
  ! loop 5
  !-------
! allocate ( mtisi0  (len_report))     ! TIME SIGNIFICANCE
  allocate ( nggtp0  (len_l_fxg,len_report))   ! TIME PERIOD OR DISPLACEMENT [minute]
! allocate ( nmwgd   (len_l_fxg,len_report))   ! MAXIMUM WIND GUST DIRECTION [degree_true]
  allocate ( nfxgu   (len_l_fxg,len_report))   ! MAXIMUM WIND SPEED(GUSTS)   [m/s]
! for buoy
  allocate ( nfxgub  (len_report))     ! MAXIMUM WIND SPEED(GUSTS)   [m/s]

  !---------
  ! loop 6
  !---------
  len_l_rad = 2
  allocate ( mggtp7 (len_l_rad,len_report)) ! TIME PERIOD OR DISPLACEMENT [hour]
  allocate ( mlwr   (len_l_rad,len_report)) ! LONG WAVE RADIAT.,INTEGR./PERIOD SPECIF. [J/M**2]
  allocate ( mglsr  (len_l_rad,len_report)) ! GLOBAL SOLAR RADIATION INTEGR.O.PER.SPEC [J/M**2]
  allocate ( mdsrh  (len_l_rad,len_report)) ! DIFFUSE SOLAR RADIAT. (HIGH ACC.) INTEGR [J/M**2]

  !-----
  ! buoy
  !-----
  if (l_buoy) then
  allocate (mtodb    (len_report))     ! TYPE OF DATA BUOY
! allocate (lmtodb   (len_report))     ! TYPE OF DATA BUOY

! allocate (lmqbst   (len_report))     ! QUALITY OF BUOY SATELLITE TRANSMISSION
! allocate (lmqobl   (len_report))     ! QUALITY OF BUOY LOCATION
  allocate (mqobl    (len_report))     ! QUALITY OF BUOY LOCATION
! allocate (lmlqc    (len_report))     ! LOCATION QUALITY CLASS (..66 CONFIDENCE)
  end if

  if (l_synop .or. l_ship .or. l_metar_usa) then
   len_cloud = 1 + len_layer
   allocate (mvtsu_tmp(len_cloud))
   allocate (mnh_tmp  (len_cloud))
!  allocate (mcc_tmp  (len_cloud))
   allocate (nh_tmp   (len_cloud))
  end if

  if (l_metar) then
   allocate (mvtsu_  (len_cloud,len_report))    ! VERTICAL SIGNIFICANCE (SURFACE OBSERV.)
   allocate (mnh_    (len_cloud,len_report))    ! CLOUD AMOUNT
   allocate (nh_     (len_cloud,len_report))    ! HEIGHT OF BASE OF CLOUD [m]
   allocate (mcc_    (len_cloud,len_report))    ! CLOUD TYPE
  end if
  !=====================
  ! get data from netcdf
  !=====================
  !----------------------------------------
  ! height of station ground above mean sea
  !----------------------------------------
  i_dft = -999
  r_dft = -999.
  call get_dat(ncid, ifd_1dm, mhosnn, ifd_2dm, rfd_2dm, 'MHOSNN',  varid_mhosnn, ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mhosnn ,rec1)
  call get_dat(ncid, ifd_1dm, mhp   , ifd_2dm, rfd_2dm, 'MHP   ',  varid_mhp   , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mhp    ,rec1)
  ! height of barometer above mean sea level
  call get_dat(ncid, ifd_1dm, mhobnn, ifd_2dm, rfd_2dm, 'MHOBNN',  varid_mhobnn, ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mhobnn ,rec1)
  i_dft =    -1
  r_dft = rvind
  !-------------------------------
  ! pressure (vert.location) [Pa]
  ! pressure                 [Pa]
  ! pressure reduced to msl  [Pa]
  ! altimeter setting (QNH)  [Pa]
  ! 3 hour pressure change   [Pa]
  ! geopotential height      [gpm]
  !-------------------------------
  call get_dat(ncid, ifd_1dm, mpn   , ifd_2dm, rfd_2dm, 'MPN'   ,  varid_mpn   , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mpn    ,rec1)
  call get_dat(ncid, ifd_1dm, mppp  , ifd_2dm, rfd_2dm, 'MPPP'  ,  varid_mppp  , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mppp   ,rec1)
  call get_dat(ncid, ifd_1dm, mpppp , ifd_2dm, rfd_2dm, 'MPPPP' ,  varid_mpppp , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mpppp  ,rec1)
  !------------------------------------
  ! Catch case where PMSL=0 is reported
  !------------------------------------
  where (mpppp == 0.) mpppp = r_dft

  call get_dat(ncid, ifd_1dm, mphph , ifd_2dm, rfd_2dm, 'MPHPH' ,  varid_mphph , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mphph  ,rec1)
  call get_dat(ncid, ifd_1dm, nppp  , ifd_2dm, rfd_2dm, 'NPPP'  ,  varid_nppp  , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_nppp   ,rec1)
  call get_dat(ncid, ifd_1dm, nhhhn , ifd_2dm, rfd_2dm, 'NHHHN' ,  varid_nhhhn , ty_f, ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_nhhhn  ,rec1)
  !-------------------------------------
  ! temperature/dry bulb temperature [K]
  ! dew-point temperature            [K]
  ! relative humidity                [%]
  !-------------------------------------

  call get_dat  (ncid, ifd_1dm, mtdbt , ifd_2dm, rfd_2dm, 'MTDBT' ,  varid_mtdbt , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtdbt  ,rec1)
  if (.not. l_mtdbt) then
    call get_dat(ncid, ifd_1dm, mtdbt , ifd_2dm, rfd_2dm, 'MTTT'  ,  varid_mtdbt , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtdbt  ,rec1)
  endif
  call get_dat  (ncid, ifd_1dm, mtdnh , ifd_2dm, rfd_2dm, 'MTDNH' ,  varid_mtdnh , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtdnh  ,rec1)
  if (.not. l_mtdnh) then
    call get_dat(ncid, ifd_1dm, mtdnh , ifd_2dm, rfd_2dm, 'MTDTD' ,  varid_mtdnh , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtdnh  ,rec1)
  endif
  call get_dat  (ncid, ifd_1dm, muuu  , ifd_2dm, rfd_2dm, 'MUUU'  ,  varid_muuu  , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_muuu   ,rec1)

  !---------------------------------------------
  ! height of temp.sensor above local ground [m]
  !---------------------------------------------
  l_mhosen  = .false.
  if (.not. l_buoy)                                                                             &
    call get_dat(ncid, ifd_1dm, mhosen, ifd_2dm, rfd_2dm,'MHOSEN' , varid_mhosen , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mhosen ,rec1)
  !---------------------------------------------
  ! height of wind sensor above local ground [m]
  !---------------------------------------------
  l_mhosen5 = .false.
  if (l_synop .or. l_ship .or. l_cman .or. l_swis .or. l_metar_usa) then
    if (l_swis) then
       char_tmp = 'MHOSEN1'
    else
       char_tmp = 'MHOSEN5'
    end if
    call get_dat(ncid, ifd_1dm,mhosen5, ifd_2dm, rfd_2dm, char_tmp, varid_mhosen5, ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mhosen5,rec1)
  end if
  !--------------------------
  ! sea/water temperature [K]
  !--------------------------
  l_mdrep = .false.
  l_mtn00 = .false.
  l_ndbss = .false.
  l_nvs   = .false.
  ndbss   = rvind
  mtn00   = rvind
  nvs     = -1
  if (l_ship) then
    call get_dat(ncid, nvs    , rfd_1dm,ifd_2dm, rfd_2dm, 'NVS'   ,  varid_nvs   , ty_i, ndim1, &
                       ifd_1d , rfd_1d ,ifd_2d2, rfd_2d2, i_dft, r_dft, l_nvs    ,rec1)
    where (nvs == i_dft) nvs = -1
    call get_dat(ncid, ifd_1dm, mtn00 , ifd_2dm, rfd_2dm, 'MTN00' ,  varid_mtn00 , ty_f, ndim1, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtn00  ,rec1)
  else if (l_buoy) then
    call get_dat(ncid, mdrep  , rfd_1dm,ifd_2dm, rfd_2dm, 'MDREP' ,  varid_mdrep , ty_i, ndim1, &
                       ifd_1d , rfd_1d ,ifd_2d2, rfd_2d2, i_dft, r_dft, l_mdrep  ,rec1)
    call get_dat(ncid, ifd_1dm, rfd_1dm,ifd_2dm, mtn00_ , 'MTN00' ,  varid_mtn00 , ty_f, ndim2, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtn00  ,rec1)
    call get_dat(ncid, ifd_1dm, rfd_1dm,ifd_2dm, nznzn  , 'NZNZN' ,  varid_ndbss , ty_f, ndim2, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_ndbss  ,rec1)
  else if (l_cman) then
    call get_dat(ncid, ifd_1dm, rfd_1dm,ifd_2dm, mtn00_ , 'MTN00' ,  varid_mtn00 , ty_f, ndim2, &
                       ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtn00  ,rec1)
    ! TODO: read NDBSS when available
  end if

  if (l_cman .and. l_mtn00 .and. len_layer >= 1) then
     do i = 1, len_report
        mtn00(i) = mtn00_(1,i)
     end do
  end if

  !---------------------------------------------------
  ! BUOYs: derive water temperature closest to surface
  !---------------------------------------------------
  if (l_buoy .and. l_mdrep .and. l_ndbss) then
     do i = 1, len_report
        if (mdrep(i) > 0) then
           j = minloc (nznzn(:mdrep(i),i), dim=1)
           ndbss(i)              = nznzn (j,i)
           if (l_mtn00) mtn00(i) = mtn00_(j,i)
        end if
     end do
  end if

  !--------------------------
  ! horizontal visibility [m]
  !--------------------------
  call get_dat(ncid, ifd_1dm, mvv    , ifd_2dm, rfd_2dm, 'MVV'   ,  varid_mvv   , ty_f, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mvv  ,rec1)
  !---------------------------------------------
  ! total precipitation, past 24 hours [kg/m**2]
  !---------------------------------------------
  r_dft = -1.
  call get_dat(ncid, ifd_1dm, mrr24 , ifd_2dm, rfd_2dm, 'MRR24' ,  varid_mrr24 , ty_f,  ndim1, &
                     ifd_1d , rfd_1d, ifd_2d1, rfd_2d1, i_dft, r_dft, l_mrr24 ,rec1)
  !-------
  ! loop 4
  !-------
  i_dft = nvind
  r_dft = -99.
  !------------------------------------------------
  ! time period or displacement [hour]
  ! total precipitation/total water equiv.[kg/m**2]
  !------------------------------------------------

  ! name in BUFR depends on template (dbkz)
  char_tmp = ''
  select case (report_subt)
  case (0, 128, 156, 166, 10000, 10128, 10015, 10143, 10150, 10158) ! synop, cman
    char_tmp = 'MGGTP1'
  case (131)                  ! metar USA
    char_tmp = 'MGGTP5'
  case (385, 10385)           ! buoy
    char_tmp = 'MGGTP'
  case (5, 10005)             ! synop, national BUFR template
    char_tmp = 'NGGTM'        ! "DURAT.OF TIME RELAT.TO FOLLOWING VALUE"
  case (170, 10170)           ! road weather data
    char_tmp = 'NGGTP'        ! duration (min)
  case default
    char_tmp = 'MGGTP0'       ! ???
  end select

  mggtp  = i_dft
  mggtp1 = i_dft
  l_mggtp1 = .false.

    if (l_buoy .or. l_nat .or. l_swis .or. l_seemhem) then
      call get_dat(ncid, mggtp  , rfd_1dm, ifd_2dm, rfd_2dm, char_tmp,  varid_mggtp1 , ty_i, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mggtp1 ,rec1)
      call get_dat(ncid, ifd_1dm, mrrrb  , ifd_2dm, rfd_2dm, 'MRRR'  ,  varid_mrrr , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mrrr ,rec1)
      !----------------------------------------------------------
      ! convert time period (min->h) for German national template
      !----------------------------------------------------------
      if (l_nat) then
        where (mggtp /= i_dft .and. modulo (mggtp, 60) == 0)
          mggtp = - mggtp / 60
        elsewhere
          mggtp = i_dft
        endwhere
      else if (l_swis) then
        print *, "WARNING: duration (min) for SWIS not yet properly dealt with, ignored"
        mggtp = i_dft
      endif
    else if (l_synop .or. l_metar .or. l_ship .or. l_metar_usa) then
      call get_dat(ncid, ifd_1dm, rfd_1dm, mggtp1 , rfd_2dm, char_tmp,  varid_mggtp1 , ty_i, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mggtp1 ,rec1)
      call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, mrrr   , 'MRRR'  ,  varid_mrrr , ty_f, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mrrr ,rec1)
    endif

! reset to default
  r_dft = rvind
  i_dft = -1
  !----------------------------------------
  ! cloud cover (total) [%]
  ! vertical significance (surface observ.)
  ! cloud amount
  ! height of base of cloud [m]
  ! cloud type
  ! cloud cover (total) [%]
  !----------------------------------------
  mn        = -1
  mvtsu     = -1
  mnh       = -1
  mcc       = -1
  mcc0      = -1
  mcc1      = -1
  mgwi      = -1
  l_mvtsu   = .false.
  l_mdrep   = .false.
  l_mdrep2  = .false.
  if (l_synop .or. l_ship .or. l_metar_usa) then
   call get_dat(ncid, mn     , rfd_1dm, ifd_2dm, rfd_2dm, 'MN'    , varid_mn    , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mn    ,rec1)
   call get_dat(ncid, mdrep  , rfd_1dm, ifd_2dm, rfd_2dm, 'MDREP' , varid_mdrep , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mdrep ,rec1)
   call get_dat(ncid, mvtsu  , rfd_1dm, ifd_2dm, rfd_2dm, 'MVTSU' , varid_mvtsu , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mvtsu ,rec1)
   call get_dat(ncid, mnh    , rfd_1dm, ifd_2dm, rfd_2dm, 'MNH'   , varid_mnh   , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mnh   ,rec1)
   call get_dat(ncid, ifd_1dm, nh     , ifd_2dm, rfd_2dm, 'NH'    , varid_nh    , ty_f, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nh    ,rec1)
   call get_dat(ncid, mcc    , rfd_1dm, ifd_2dm, rfd_2dm, 'MCC'   , varid_mcc   , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mcc   ,rec1)
   call get_dat(ncid, mcc0   , rfd_1dm, ifd_2dm, rfd_2dm, 'MCC0'  , varid_mcc0  , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mcc0  ,rec1)
   call get_dat(ncid, mcc1   , rfd_1dm, ifd_2dm, rfd_2dm, 'MCC1'  , varid_mcc1  , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mcc1  ,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, mvtsu0 , rfd_2dm, 'MVTSU0', varid_mvtsu0, ty_i, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mvtsu0,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, mnh0   , rfd_2dm, 'MNH0'  , varid_mnh0  , ty_i, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mnh0  ,rec1)
!  call get_dat(ncid, ifd_1dm, rfd_1dm, mcc_   , rfd_2dm, 'MCC2'  , varid_mcc2  , ty_i, ndim2, &
!                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mcc2  ,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, nh0    , 'NH0'   , varid_nh0   , ty_f, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nh0   ,rec1)
  else if (l_metar) then
   mvtsu_    = -1
   mnh_      = -1
   mcc_      = -1
   nh_       = rvind
   call get_dat(ncid, mgwi   , rfd_1dm, ifd_2dm, rfd_2dm, 'MGWI'  , varid_mgwi  , ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mgwi  ,rec1)
   call get_dat(ncid, mdrep2 , rfd_1dm, ifd_2dm, rfd_2dm, 'MDREP2', varid_mdrep2, ty_i, ndim1, &
                      ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mdrep2,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, mvtsu_ , rfd_2dm, 'MVTSU' , varid_mvtsu , ty_i, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d4, rfd_2d4, i_dft, r_dft, l_mvtsu ,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, mnh_   , rfd_2dm, 'MNH'   , varid_mnh   , ty_i, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d4, rfd_2d4, i_dft, r_dft, l_mnh   ,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, mcc_   , rfd_2dm, 'MCC'   , varid_mcc   , ty_i, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d4, rfd_2d4, i_dft, r_dft, l_mcc   ,rec1)
   call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, nh_    , 'NH'    , varid_nh    , ty_f, ndim2, &
                      ifd_1d , rfd_1d , ifd_2d4, rfd_2d4, i_dft, r_dft, l_nh    ,rec1)
   l_mcc0    = .false.
   l_mcc1    = .false.
  end if

  !-------------------
  ! radiances (loop 6)
  !-------------------
  if ((.not. l_metar_usa) .and. (.not. l_seemhem)) then
  call get_dat(ncid, ifd_1dm, rfd_1dm, mggtp7 , rfd_2dm, 'MGGTP7', varid_mggtp7, ty_i, ndim2, &
                     ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_mggtp7 ,rec1)
  call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, mlwr   , 'MLWR'  , varid_mlwr  , ty_f, ndim2, &
                     ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, -99.,  l_mlwr   ,rec1)
  call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, mglsr  , 'MGLSR' , varid_mglsr , ty_f, ndim2, &
                     ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, -99.,  l_mglsr  ,rec1)
  call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, mdsrh  , 'MDSRH' , varid_mdsrh , ty_f, ndim2, &
                     ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, -99.,  l_mdsrh  ,rec1)
  end if
  !-----------------------------------
  ! present weather
  ! time period or displacement [hour]
  ! past weather (1)  past weather (2)
  !-----------------------------------
  nww       = -1
  mw1       = -1
  mw2       = -1
  if (.not. l_buoy .and. .not. l_nat .and. .not. l_swis) mggtp     = -1
  if (.not. l_swis) then
  call get_dat(ncid, nww    , rfd_1dm, ifd_2dm, rfd_2dm, 'NWW'   ,  varid_nww   , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nww   ,rec1)
  if (.not.l_buoy .and. .not. l_nat)                                                           &
  call get_dat(ncid, mggtp  , rfd_1dm, ifd_2dm, rfd_2dm, 'MGGTP' ,  varid_mggtp , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mggtp ,rec1)
  call get_dat(ncid, mw1    , rfd_1dm, ifd_2dm, rfd_2dm, 'MW1'   ,  varid_mw1   , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mw1   ,rec1)
  call get_dat(ncid, mw2    , rfd_1dm, ifd_2dm, rfd_2dm, 'MW2'   ,  varid_mw2   , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mw2   ,rec1)
  end if
  !--------------------------
  ! loop 000 not expanded yet
  ! loop 001 not expanded yet
  ! loop 002 not expanded yet
  ! loop 003 not expanded yet
  !--------------------------

  !-------------------------------------------
  ! time period or displacement [hour]
  ! maximum temp., height/period specified [K]
  ! time period or displacement [hour]
  ! minimum temp., height/period specified [K]
  !-------------------------------------------
  i_dft = nvind
  l_mggtp2 = .false.
  l_mtxtxh = .false.
  l_mggtp4 = .false.
  l_mtntnh = .false.

  ! name in BUFR depends on template
  char_tmp = ''
  if (l_synop) then
    char_tmp = 'MGGTP2'     ! Time range: MGGTP2..MGGTP3 (SYNOP)
  else if (l_ship) then
    char_tmp = 'MGGTP1'     ! Time range: MGGTP1..MGGTP2 (SHIP)
  else if (l_metar_usa) then
    char_tmp = 'MGGTP6'     ! Time range: MGGTP6..MGGTP7 (METAR USA)
  end if
  if (l_synop .or. l_ship .or. l_metar_usa) then
  call get_dat(ncid, mggtp2 , rfd_1dm, ifd_2dm, rfd_2dm, char_tmp, varid_mggtp2 , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mggtp2 ,rec1)
  call get_dat(ncid, ifd_1dm, mtxtxh , ifd_2dm, rfd_2dm, 'MTXTXH', varid_mtxtxh , ty_f, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtxtxh ,rec1)
  end if

  ! name in BUFR depends on template
  char_tmp = ''
  if (l_synop) then
    char_tmp = 'MGGTP4'     ! Time range: MGGTP4..MGGTP5 (SYNOP)
  else if (l_ship) then
    char_tmp = 'MGGTP3'     ! Time range: MGGTP3..MGGTP4 (SHIP)
  else if (l_metar_usa) then
    char_tmp = 'MGGTP8'     ! Time range: MGGTP8..MGGTP9 (METAR USA)
  end if
  if (l_synop .or. l_ship .or. l_metar_usa) then
  call get_dat(ncid, mggtp4 , rfd_1dm, ifd_2dm, rfd_2dm, char_tmp, varid_mggtp4 , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mggtp4 ,rec1)
  call get_dat(ncid, ifd_1dm, mtntnh , ifd_2dm, rfd_2dm, 'MTNTNH', varid_mtntnh , ty_f, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtntnh ,rec1)
  end if
  !-------------------------------------
  ! height of sensor above water surface
  !-------------------------------------
  if (l_buoy) then
    call get_dat (ncid, ifd_1dm, mhawas , ifd_2dm, rfd_2dm, 'MHAWAS' , varid_mhawas , ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft,  r_dft, l_mhawas ,rec1)
    call get_dat (ncid, ifd_1dm, mhawas1, ifd_2dm, rfd_2dm, 'MHAWAS1', varid_mhawas1, ty_f, ndim1, &
                        ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft,  r_dft, l_mhawas1,rec1)
  else
    l_mhawas  = .false.
    l_mhawas1 = .false.
  end if
  !-------------------------------------
  ! time period or displacement [minute]
  ! wind direction [degree_true]
  ! wind speed [m/s]
  !-------------------------------------
  call get_dat  (ncid, nggtp  , rfd_1dm, ifd_2dm, rfd_2dm, 'NGGTP' , varid_nggtp  , ty_i, ndim1, &
                       ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nggtp ,rec1)
    if (l_seemhem) then
    call get_dat  (ncid, ifd_1dm, rfd_1dm, ifd_2dm, ndndn_s, 'NDNDN', varid_ndndn   , ty_f, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_ndndn ,rec1)
    else
    call get_dat  (ncid, ifd_1dm, ndndn  , ifd_2dm, rfd_2dm, 'NDNDN', varid_ndndn   , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_ndndn ,rec1)
    endif
    if (.not. l_ndndn) then
      call get_dat(ncid, ifd_1dm, ndndn  , ifd_2dm, rfd_2dm, 'NDD',   varid_ndndn   , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_ndndn ,rec1)
    endif
    if (l_seemhem) then
    call get_dat  (ncid, ifd_1dm, rfd_1dm, ifd_2dm, nfnfn_s, 'NFNFN', varid_nfnfn   , ty_f, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nfnfn ,rec1)
    else
    call get_dat  (ncid, ifd_1dm, nfnfn  , ifd_2dm, rfd_2dm, 'NFNFN' , varid_nfnfn  , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nfnfn ,rec1)
    endif
    if (.not. l_nfnfn) then
      call get_dat(ncid, ifd_1dm, nfnfn  , ifd_2dm, rfd_2dm, 'NFF' ,   varid_nfnfn  , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nfnfn ,rec1)
    endif
  !-------
  ! loop 5
  !-------
  !-------------------------------------
  ! time period or displacement [minute]
  ! maximum wind speed(gusts)   [m/s]
  !-------------------------------------
!  if (l_buoy) then
!  write(0,*) '***************************************************************'
!  write(0,*) 'read_synop_netcdf: reading NGGTP0'
!  write(0,*)
!  write(0,*) '************************* FIXME!!! ****************************'
!  write(0,*)
!  write(0,*) 'get_dat not called, may crash when reading buoy data!'
!  write(0,*) 'using "invalid values" instead'
!  write(0,*) '***************************************************************'
  nggtp0   = i_dft
  l_nggtp0 = .false.
  r_dft    = -99.

    if (l_buoy .or. l_nat .or. l_swis .or. l_seemhem) then

      select case (report_subt)
      case (5, 10005)            ! German national SYNOP BUFR
        char_tmp = 'NGGTM6'
      case (170, 10170)          ! Road weather stations
        char_tmp = 'NGGTP1'
      case default               ! buoy
        char_tmp = 'NGGTP0'
      end select
      call get_dat(ncid, nggtp  , rfd_1dm, ifd_2dm, rfd_2dm, char_tmp, varid_nggtp  , ty_i, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nggtp0 ,rec1)

      select case (report_subt)
      case (5, 10005)            ! German national SYNOP BUFR
        char_tmp = 'NFXGU0'
      case default               ! buoy
        char_tmp = 'NFXGU'
      end select
      if (l_seemhem) then
      call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm,   nfxgu, char_tmp, varid_nfxgu  , ty_f, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nfxgu  ,rec1)
      else
      call get_dat(ncid, ifd_1dm, nfxgub , ifd_2dm, rfd_2dm, char_tmp, varid_nfxgu  , ty_f, ndim1, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nfxgu  ,rec1)
      endif
      if (l_nat) then
        !--------------------------------------------------------
        ! convert time period (sign) for German national template
        !--------------------------------------------------------
        where (nggtp /= i_dft .and. modulo (nggtp, 60) == 0)
          nggtp = - nggtp
        elsewhere
          nggtp = i_dft
        endwhere
      endif

    else
      call get_dat(ncid, ifd_1dm, rfd_1dm, nggtp0 , rfd_2dm, 'NGGTP0', varid_nggtp0 , ty_i, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nggtp0 ,rec1)
      if (.not. l_metar) then
      call get_dat(ncid, ifd_1dm, rfd_1dm, ifd_2dm, nfxgu  , 'NFXGU' , varid_nfxgu  , ty_f, ndim2, &
                         ifd_1d , rfd_1d , ifd_2d2, rfd_2d2, i_dft, r_dft, l_nfxgu  ,rec1)
      else
         !-----------------------------------
         ! METAR gust data not useful for NWP
         !-----------------------------------
         l_nfxgu = .false.
         nfxgu   = -1
      end if
    end if
  r_dft = rvind

  !---------------------------------------
  ! state of ground (with or without snow)
  !---------------------------------------
  i_dft = -1
  me    = -1
  call get_dat(ncid, me     , rfd_1dm, ifd_2dm, rfd_2dm, 'ME' , varid_me  , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_me    ,rec1)
  !-------------------------------------------------------------------
  ! total snow depth [m]; hint: -.01 rest of snow -.02 patches of snow
  !-------------------------------------------------------------------
  call get_dat(ncid, ifd_1dm, nsss   , ifd_2dm, rfd_2dm, 'NSSS'  , varid_nsss   , ty_f, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_nsss  ,rec1)
  !-----
  ! buoy
  !-----
  l_ma1   = .false.
  l_mwrsa = .false.
  l_mqobl = .false.
  l_mtodb = .false.
  if (l_buoy) then
  call get_dat(ncid, ma1    , rfd_1dm, ifd_2dm, rfd_2dm, 'MA1'  ,  varid_ma1  , ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_ma1  , rec1)
  call get_dat(ncid, mwrsa  , rfd_1dm, ifd_2dm, rfd_2dm, 'MWRSA',  varid_mwrsa, ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mwrsa, rec1)
  call get_dat(ncid, mqobl  , rfd_1dm, ifd_2dm, rfd_2dm, 'MQOBL',  varid_mqobl, ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mqobl, rec1)
  call get_dat(ncid, mtodb  , rfd_1dm, ifd_2dm, rfd_2dm, 'MTODB',  varid_mtodb, ty_i, ndim1, &
                     ifd_1d , rfd_1d , ifd_2d1, rfd_2d1, i_dft, r_dft, l_mtodb, rec1)
  end if

  !-------------------------------------------------------------------------
  ! Fix station IDs of "new buoys" using WMO region number and WMO sub-area:
  ! Fit IDs < 1000 into 5-digit station ID, otherwise extend to 7-digits
  ! Ref.: http://www.wmo.int/pages/prog/amp/mmop/wmo-number-rules.html
  !-------------------------------------------------------------------------
  if (l_buoy .and. l_ma1 .and. l_mwrsa) then
     mwrsa =      min (mwrsa, 9)                ! Note: missing value = -1
     ma1   = max (min (ma1  , 9), 0)
     where (ma1 == 0) ma1 = 7                   ! Antarctic (0) is mapped to 7
     do i = 1, len_report
        if (istidn(i) /= imissing .and. mwrsa(i) >= 0) then
           select case (istidn(i))
           case (      :999)
              istidn(i) = istidn(i) + ma1(i)*  10000 + mwrsa(i)*  1000  ! A1 bw nnn
              write( ystidn(i),'(I5.5)' )  istidn(i)
           case (  1000:9999)
              istidn(i) = istidn(i) + ma1(i)*1000000 + mwrsa(i)*100000  ! A1 bw nnnnn
              write( ystidn(i),'(I7.7)' )  istidn(i)
           case (100000:)
              write( ystidn(i),'(I7.7)' )  istidn(i)
           case default
              ! Station ID is in the range 10000-99999, probably "old buoy"
           end select
        endif
     enddo
  end if

  !----------------------------------
  ! list of defined variables in file
  !----------------------------------
  if (netcdf_verb > 0) then
! meteorological part
    write (6,'(a,i3,1x,a,/,2x, 9l8)') 'pe=',dace% pe , &
     'l_mpn  l_mppp l_mpppp l_nhhhn l_mtdbt l_mtdnh l_muuu   l_mvv  l_mrr24', &
      l_mpn, l_mppp,l_mpppp,l_nhhhn,l_mtdbt,l_mtdnh,l_muuu,  l_mvv, l_mrr24
    write (6,'(a,i3,1x,a,/,3x,10l7)') 'pe=',dace% pe , &
     'l_mrrr l_mn   l_mnh  l_nh  l_mcc l_mcc0 l_mcc1 l_ndndn l_nfnfn l_nfxgu',&
      l_mrrr,l_mn,  l_mnh, l_nh, l_mcc,l_mcc0,l_mcc1,l_ndndn,l_nfnfn,l_nfxgu
    write (6,'(a,i3,1x,a,/,2x, 9l8)')  'pe=',dace% pe , &
     'l_me    l_nww   l_mw1   l_mw2 l_mtxtxh l_mtntnh l_nsss l_mtn00 l_ndbss',&
      l_me,   l_nww,  l_mw1,  l_mw2,l_mtxtxh,l_mtntnh,l_nsss,l_mtn00,l_ndbss
    write (6,'(a,i3,1x,a,/,2x, 8l9)')  'pe=',dace% pe , &
     'l_mphph l_mhobnn l_mhosnn   l_nvs   l_mhp',&
      l_mphph,l_mhobnn,l_mhosnn,  l_nvs,  l_mhp

! technological / cloud part
    write (6,'(a,i3,1x,a,/,3x, 9l9)') 'pe=',dace% pe , &
     'l_mggtp1 l_nggtp  l_nggtp0 l_mggtp  l_mggtp2 l_mggtp4 l_mqobl  l_ma1    l_mwrsa',&
      l_mggtp1,l_nggtp, l_nggtp0,l_mggtp, l_mggtp2,l_mggtp4,l_mqobl, l_ma1,   l_mwrsa
    write (6,'(a,i3,1x,a,/,3x, 9l9)') 'pe=',dace% pe , &
     'l_mhosen l_mhosen5 l_mhawas l_mhawas1 l_mtodb',&
      l_mhosen,l_mhosen5,l_mhawas,l_mhawas1,l_mtodb
  end if
  if (prt_data) then
    write( *, '("### Cloud -                 reported     derived  "           &
         &     ," derived (incl individual cloud layers)")' )
    write( *, '("### Cloud --          general cloud group   |     "           &
         &     ," and   saved   as   observations")' )
    write( *, '("### Cloud ---          20003 oct 8002 [m] index   "           &
         &     ," cld-cov-octa   [m]   [m]    [m]")' )
    write( *, '("### Cloud -sta ID   UTC   WW MNH vsig CBH  fog    "           &
         &     ," N_L _M _H  N   CBH  CEIL    VIS")' )
  endif

  !-------------------------------
  ! preset total number of reports
  !-------------------------------
  entry   = sum (source(1:ifile-1)% entries) + rec1 -1

  !------------------
  ! loop over reports
  !------------------
! netcdf  : is number of report
! do is = 1,min( 2,len_report)
! do is = 1,min(10,len_report)
  do is = 1, len_report

  entry1  = entry    + 1
  entry   = entry    + 1

  !---------------
  ! initialize use
  !---------------
  use = use_0

  !--------------------
  ! define head section
  !--------------------
  !=======================================
  ! derive observation type specifications
  !=======================================
  report_subt  = s2ikz(is)
  bufr_type    = s1cat(is)
  bufr_subtype = s1catls(is)
  centre       = s1cent(is)

  if (report_subt >= 0) then
    !---------------------------------
    ! derive information from DWD dbkz
    !---------------------------------
    obsid      = obstype_dbkz (report_subt)
  else
    !---------------------------------------------------------------
    ! or from bufr_type, bufr_subtype specified by generating center
    !---------------------------------------------------------------
    obsid      = obstype_bufr (bufr_type, bufr_subtype, centre)
    !------------------------
    ! optionally set DWD dbkz
    !------------------------
    if (derive_dbkz) report_subt   = obsid% dbkz
  endif
  !-------------------------------------------------------
  ! set CMA obstype, BUFR type, BUFR subtype, if not given
  !-------------------------------------------------------
  if (bufr_type   <0) bufr_type    = obsid% bufrtype
  if (bufr_subtype<1 .and. obsid% centre == WMO0_ECMWF) &
                      bufr_subtype = obsid% subtype
  obstype                          = obsid% obstype
  if (obstype < 0) cycle
  report_subti      = idb_dbk (report_subt, obstype)
  !---------------------------------------------------------------------------
  ! Fix bufr_type, bufr_subtype missing values for otherwise known "good" data
  !---------------------------------------------------------------------------
  if (obsid% centre == WMO0_DWD) then
     if (bufr_type    == 255) bufr_type    = obsid% bufrtype
     if (bufr_subtype == 255) bufr_subtype = obsid% subtype
  end if

  head% obstype     = obstype
  head% dbkz        = report_subt
! head% modtype     = rept_char(obstype)% mod
  head% modtype     = SYNOP
  head% buf_type    = bufr_type
  head% buf_subtype = bufr_subtype
  head% codetype    = obsid% codetype
  head% time        = stime(is)
  head% db_time     = db_time(is)
  head% idbk        = report_subti
  head% source      = ifile
  head% record      = is + rec1 - 1
  head% id          = entry1
  head% center      = s1cent(is)
  head% subcenter   = s1cents(is)

  if ( lpr_synop .and. is < npr_extd ) then
    write (6,'()')
    write (6,'( 8(a16, i6  ,/),   &
              & 2(a16, a   ,/),   &
              & 6(a16, i6  ,/) )' )                      &
      'pe='         ,dace% pe,                           &
      'head is='    ,is,                                 &
      'obstype='    , head% obstype  ,                   &
      'dbkz='       , head% dbkz     ,                   &
      'modtype='    , head% modtype  ,                   &
      'buf_type='   , head% buf_type ,                   &
      'buf_subtype=', head% buf_subtype,                 &
      'codetype='   , head% codetype   ,                 &
      'time='       , cyyyymmddhhmmss (head% time)   ,   &
      'db_time='    , cyyyymmddhhmmss (head% db_time),   &
      'dbk='        , head% idbk     ,                   &
      'source='     , head% source   ,                   &
      'record='     , head% record   ,                   &
      'id='         , head% id       ,                   &
      'center='     , head% center    ,                  &
      'subcenter='  , head% subcenter
  endif

  !--------------------------------------------
  ! perform simple generic check on report type
  !--------------------------------------------
  call check_report_0 (use, head, 1)
  if (use% state <= STAT_DISMISS) then
    if (lpr_synop .and. is < npr_extd) write (6,*) '  check_report_0: state =', use% state
    lkeep = .false.
    cycle
  endif
  !------------------
  ! create new report
  !------------------
  spt0             = empty
  spt0% use        = use
  spt0% hd         = head
  spti             = spt0
! spti% hd% subset = is

  lkeep = .false.
  call construct_synop (s)

  !------------------------------
  ! process the following entries
  !------------------------------
  spti% corme        = max ( s1updat(is), 0)
  spti% col% c% dlat = mlah(is)
  spti% col% c% dlon = mloh(is)
  spti% z            = mhobnn(is)
  spti% actual_time  = obs_time(is)
  spti% statid       = ystidn(is)
  spti% ident        = istidn(is)
  !--------------------------------------------------------------------
  ! Map station id to integer "hash value" for vectorizable comparisons
  ! Be careful about reports with invalid/non-standard WIGOS station ID
  ! that may report redundantly using traditional station identifiers.
  !--------------------------------------------------------------------
  if (lwsi(is)) then
     spti% wsi       = wsi    (is)
  end if
  if (lwsi(is) .and. spti% statid(1:1) == "_") then
     spti% stat_hash = wsihash(is)
  else
     spti% stat_hash = transfer (spti% statid, spti% stat_hash)
  end if

  !------------------------------------------------------------
  ! Check if automatic ship report is actually from buoy or rig
  ! guess by looking at the station identifier
  !------------------------------------------------------------
  if (l_ship .and. spti% hd% codetype == 24 .and. chk_ship > 0) then
     len = len_trim (spti% statid)
     chk = ( len == 5 .or. len == 7 ) .and.                                  &
           (spti% statid(1:1) >= "1") .and. (spti% statid(1:1) <= "7") .and. &
           (spti% statid(2:2) >= "1") .and. (spti% statid(2:2) <= "9")
     if (chk) then
        chk = verify (spti% statid(1:len),"0123456789") == 0
        if (chk) then
           id = -1
           read (spti% statid(1:len),*,iostat=status) id
           chk = (status == 0)
           if (chk) then
              chk = ((len == 5) .and. (id < 74999  )) .or. &
                    ((len == 7) .and. (id < 7499999))
              if (chk) spti% ident = id
              if (chk .and. verbose > 0) &
                   write(6,'(a,a,i6,2f10.3)') "DRIBU as SHIP?  ident,dbkz,lat,lon = ",&
                   spti% statid, spti% hd% dbkz, spti% col% c% dlat, spti% col% c% dlon
           end if
        end if
     end if
  end if

  if ( lpr_synop  .and. is < npr_extd ) then
    write (6,'()')
    write (6,'(   a21, i6  ,a, /, &
              &   a21, i6  ,   /, &
              & 2(a21,f8.3 ,   /),&
              &   a21,2i8  ,   / ,&
              &   a21, a   ,   / ,&
              &   a21, i5      / )' )                    &
          'pe= ',dace% pe,'  spti',                      &
          'spti% corme        = ', spti% corme        ,  &
          'spti% col% c% dlat = ', spti% col% c% dlat ,  &
          'spti% col% c% dlon = ', spti% col% c% dlon ,  &
          'spti% actual_time  = ', spti% actual_time  ,  &
          'spti% statid       = ', spti% statid       ,  &
          'spti% ident        = ', spti% ident
  endif

  !--------------------
  ! location (vertical)
  !--------------------
! default: -999._sp.  Bad values (9999) have been seen...
! height of barometer above mean sea level
  if ( mhobnn(is) /= -999._sp .and. mhobnn(is) /= 9999._sp ) then
    call set_datum (s% zr ,mhobnn(is) ,spti% corme)
  endif

! for coastal stations set default station height to 0 when values missing
  if (l_cman .and. mhobnn(is) == -999._sp .and. mhosnn(is) == -999._sp) then
     mhosnn(is) = 0._sp
  end if

! height of station ground above mean sea  [m]
  if ( mhosnn(is) /= -999._sp .and. mhosnn(is) /= 9999._sp ) then
    call set_datum (s% z_ls ,mhosnn(is) ,spti% corme)
  else if ( mhp(is) /= -999._sp .and. mhp(is) /= 9999._sp ) then
    call set_datum (s% z_ls ,mhp(is) ,spti% corme)
  endif

  !-------------------------------
  ! vertical elements and pressure
  !-------------------------------
  ! pressure (vert.location)                [Pa]
  call set_datum (s% p     ,mpn(is)   ,spti% corme)
  ! pressure                                [Pa]
  call set_datum (s% ps    ,mppp(is)  ,spti% corme)
  ! q-BITS for pressure                     [code table]
!!if ( mlahq(is) /= qc_dft ) call set_qbits (s% ps , mlahq(is))
  ! pressure reduced to mean sea level      [Pa]
  call set_datum (s% p_msl ,mpppp(is) ,spti% corme)
  ! q-bits for pressure reduced to mean sea level
!!if ( mlahq(is) /= qc_dft ) call set_qbits (s% p_msl ,mlahq(is))

  ! 3 hour pressure change                  [Pa]
  s% p_tend = nppp(is)

  ! altimeter setting (QNH) (METAR)         [Pa]
  if (l_mphph) &
  call set_datum (s% p_msl ,mphph(is) ,spti% corme)

  ! geopotential   [gpm] corresponds to pressure (vert.location) mpn
  call set_datum (s% gp    ,nhhhn(is) ,spti% corme)

  if (l_ship)     s% nvs = nvs(is)

  !-----------
  ! wind, gust
  !-----------
  if (l_buoy) then
    call set_datum (s% h_ts ,mhawas (is) ,spti% corme)  ! Height of temp.sensor
    call set_datum (s% h_ws ,mhawas1(is) ,spti% corme)  ! Height of wind sensor
  end if
  !-----------------------------
  ! wind direction [degree_true]
  !-----------------------------
  if ( l_seemhem ) then
    call set_datum (s% dd ,ndndn_s(1,is) ,spti% corme)
  else
    call set_datum (s% dd ,ndndn(is)     ,spti% corme)
  endif
  !-----------------
  ! wind speed [m/s]
  !-----------------
  if ( l_seemhem ) then
    call set_datum (s% dd ,nfnfn_s(1,is) ,spti% corme)
  else
    call set_datum (s% ff ,nfnfn(is)     ,spti% corme)
  endif
  !---------------------------------------------------------
  ! time period or displacement [minute]  e [-10,-60, -360 ]
  ! verification needs only -60, -360
  ! maximum wind speed(gusts)   [m/s]
  !---------------------------------------------------------
! default

    if (l_buoy .or. l_nat .or. l_swis .or. l_seemhem) then
      select case (nggtp  (is))
      case ( -60 )
       if (l_seemhem) then
        s% gust1 = nfxgu (1,is)
       else
        s% gust1 = nfxgub (is)
       endif
      case ( -180 )
        s% gust3 = nfxgub (is)
      case ( -360 )
        s% gust6 = nfxgub (is)
      case default
        continue
      end select
    else
      do  j = 1, size (nggtp0(:,is))
        select case (nggtp0(j,is))
        case ( -60 )
          s% gust1 = nfxgu (j,is)
        case ( -180 )
          s% gust3 = nfxgu (j,is)
        case ( -360 )
          s% gust6 = nfxgu (j,is)
        case default
          continue
        end select
      end do
    endif

  !----------------------------------
  ! temperature, dewpoint-temperature
  ! relative humidity
  !----------------------------------
  ! temperature/dry bulb temperature [K]
  ! dew-point temperature            [K]
  ! relative humidity                [ ]
  call set_datum (s% t , mtdbt(is) ,spti% corme)
  call set_datum (s% td, mtdnh(is) ,spti% corme)
  call set_datum (s% rh, muuu (is) ,spti% corme)
  s% rh% o = 0.01 * s% rh% o       ! [%] -> [ ]

  if (l_mhosen ) call set_datum (s% h_ts, mhosen (is) ,spti% corme)
  if (l_mhosen5) call set_datum (s% h_ws, mhosen5(is) ,spti% corme)

  !-------------------------------------------
  ! maximum temp., height/period specified [K]
  ! minimum temp., height/period specified [K]
  !-------------------------------------------
  ! verification needs -12h, -1h (aggregation)
  !-------------------------------------------
  if (l_mggtp2) then
    select case (mggtp2(is))
    case (-1)
      s% tmax1 = mtxtxh(is)
    case (-12)
      s% tmax  = mtxtxh(is)
    end select
  end if
  if (l_mggtp4) then
    select case (mggtp4(is))
    case (-1)
      s% tmin1 = mtntnh(is)
    case (-12)
      s% tmin  = mtntnh(is)
    end select
  end if

  !----------------------
  ! Sea/water temperature
  !----------------------
  if (l_mtn00) call set_datum (s% tsea,   mtn00(is) ,spti% corme)
  if (l_ndbss) call set_datum (s% h_tsea, ndbss(is), spti% corme)

  !-------------------------------------------
  ! present weather
  ! time period or displacement [hour]
  ! past weather (1)  past weather (2) i.e. 1,6 hour
  ! horizontal visibility [m]
  !-------------------------------------------
  if (nww(is) /= -1) s% ww  = nww(is)
! if ( mw1(is) /= -1 ) then
!   if ( mw2(is) /= -1 ) then
!     ivalue = max( mw1(is), mw2(is) )
!   else
!     ivalue =      mw1(is)
!   endif
! else if ( mw2(is) /= -1 ) then
!   ivalue =      mw2(is)
! else
!   ivalue =      -1
! endif
  ! verification uses only past weather (2)
  if (mw2(is) /= -1) s% w   = mw2(is)
  s% vis = mvv(is)

  !-------------------------------------------
  ! cloud cover (total) [%]
  ! vertical significance (surface observ.)
  ! cloud amount
  ! height of base of cloud [m]
  ! cloud type low medium high
  !-------------------------------------------
  if (l_synop .or. l_ship .or. l_metar_usa) then
    s% n   = mn    (is)
    lm1    = mdrep (is)
    nl     = lm1 + 1
    if (nl > 0) then
      mvtsu_tmp(:) = [ mvtsu (is), mvtsu0 (:,is) ]
      mnh_tmp  (:) = [ mnh   (is), mnh0   (:,is) ]
!     mcc_tmp  (:) = [ mcc   (is), mcc2   (:,is) ]
      nh_tmp   (:) = [ nh    (is), nh0    (:,is) ]
      mcc_lmh(1:3) = [ mcc(is), mcc0(is), mcc1(is) ]
      if (vprc_cloud == 1) then
        call cloudcover_synop (s, s% n, nl, mvtsu_tmp, mnh_tmp, mcc_lmh, nh_tmp)
      else
        call cloudcover_synop_lam (s, s% n, nl, mvtsu_tmp, mnh_tmp, mcc_lmh     &
                                  ,nh_tmp, s% vis, s% ww                        &
                                  , spti% statid, spti% actual_time% secs )
!                                 , spti% statid )
      endif
      if (prt_data) print *, "### is,statid,n,clcl,mdrep,mvtsu,mnh,nh,mcc =",    &
           is, spti% statid, s% n, s% clcl, mdrep(is), mvtsu_tmp(1), mnh_tmp(1), &
           nh_tmp(1), mcc_lmh
    end if
  else if (l_metar) then
    call cloudcover_metar (s, mgwi(is), mdrep2(is), mvtsu_(:,is), &
                           mnh_(:,is), mcc_(:,is), nh_(:,is))
    if (prt_data) print *, "### is,statid,clcl,mgwi,mdrep2,mnh,nh =",         &
         is, spti% statid, s% clcl, mgwi(is), mdrep2(is), mnh_(1,is), nh_(1,is)
  end if

  !---------------------------------------
  ! state of ground (with or without snow)
  ! total snow depth [m]
  !---------------------------------------
  if (l_me  ) s% e      = me   (is)
  if (l_nsss) s% sdepth = nsss (is)
  !----------------------------------------------------------------------
  ! set special values to zero: -0.01 rest of snow, -0.02 patches of snow
  !----------------------------------------------------------------------
  if (l_nsss .and. s% sdepth > -0.021) s% sdepth = max (0., s% sdepth)

  !------------------------------------------------
  ! total precipitation/total water equiv.[kg/m**2]
  ! values of -0.1 allowed(not measurable prec.)
  ! time period or displacement [hour] e [-1, -6, -9, -12, -18,  -24]
  ! verification needs only [-1, -6, -12, -24]
  !------------------------------------------------
  ! default: -1.0
    if (l_buoy .or. l_nat .or. l_swis)  then
      select case (mggtp (is))
      case ( -1  )
        s% rr1  = mrrrb (is)
      case ( -3  )
        s% rr3  = mrrrb (is)
      case ( -6   )
        s% rr6  = mrrrb (is)
      case ( -12  )
        s% rr12 = mrrrb (is)
      case ( -24  )
        s% rr24 = mrrrb (is)
      case default
        continue
      end select
    else
      do  j = 1, size (mggtp1(:,is))
        select case (mggtp1(j,is))
        case ( -1  )
          s% rr1  = mrrr  (j,is)
        case ( -3  )
          s% rr3  = mrrr  (j,is)
        case ( -6   )
          s% rr6  = mrrr  (j,is)
        case ( -12  )
          s% rr12 = mrrr  (j,is)
        case ( -24  )
          s% rr24 = mrrr  (j,is)
        case default
          continue
        end select
      end do
      if (s% rr24 < -0.11 .and. mrr24(is) >= -0.11 ) s% rr24 = mrr24(is)
    endif

  !-------------------------------
  ! set special value -0.1 to zero
  !-------------------------------
  if (s% rr1  > -0.11) s% rr1  = max (0., s% rr1 )
  if (s% rr3  > -0.11) s% rr3  = max (0., s% rr3 )
  if (s% rr6  > -0.11) s% rr6  = max (0., s% rr6 )
  if (s% rr12 > -0.11) s% rr12 = max (0., s% rr12)
  if (s% rr24 > -0.11) s% rr24 = max (0., s% rr24)

  !-----------------------------------
  ! radiation (loop 6)
  ! 1 and 24 h present, we only use 1h
  !-----------------------------------
  if (l_mggtp7) then
    do  j = 1, len_l_rad
      select case (mggtp7(j,is))
      case ( -1 )
        if (l_mlwr ) s% lwr1 = mlwr  (j,is)
        if (l_mglsr) s% gsr1 = mglsr (j,is)
        if (l_mglsr) s% dsr1 = mdsrh (j,is)
      case default
      end select
    end do
  end if

  !--------------------------------------------------
  ! write placeholder (invalid value) for aggregation
  !--------------------------------------------------
  if (l_synop .or. l_nat .or. l_metar_usa) then
    rhour = hours (spti% actual_time)
    ihour = nint (rhour)
!   if ((rhour - ihour) == 0._wp) then
!   if (abs (rhour - ihour) < 0.16_wp) then     ! up to ~9 minutes (US BUFR)
    if (abs (rhour - ihour) < 0.34_wp) then     ! up to ~20 minutes (Hungary)

      if (s% rr1  > -0.11 .or. &
          s% rr3  > -0.11 .or. &
          s% rr6  > -0.11 .or. &
          s% rr12 > -0.11      ) then
        if (s% rr1  < -0.11 .and. ag_rr_1h                          ) s% rr1 = rvind
        if (s% rr3  < -0.11 .and. ag_rr_3h .and. modulo (ihour,3)==0) s% rr3 = rvind
        if (s% rr6  < -0.11 .and. ag_rr_6h .and. modulo (ihour,6)==0) s% rr6 = rvind
      endif

      if   (s% gsr1 >= 0.      ) then
        if (s% gsr3 <  0. .and. ag_rd_3h .and. modulo (ihour,3)==0) s% gsr3 = rvind
        if (s% gsr6 <  0. .and. ag_rd_6h .and. modulo (ihour,6)==0) s% gsr6 = rvind
      endif

      if   (s% dsr1 >= 0.      ) then
        if (s% dsr3 <  0. .and. ag_rd_3h .and. modulo (ihour,3)==0) s% dsr3 = rvind
        if (s% dsr6 <  0. .and. ag_rd_6h .and. modulo (ihour,6)==0) s% dsr6 = rvind
      endif

      if (s% gust1  >= 0. .or. &
          s% gust3  >= 0.      ) then
        if (s% gust3  < 0. .and. ag_gust_3h .and. modulo (ihour,3)==0) s% gust3 = rvind
        if (s% gust6  < 0. .and. ag_gust_6h .and. modulo (ihour,6)==0) s% gust6 = rvind
      endif

      !-------------------------------------------
      ! Aggregate hourly tmax/tmin to 12h every 6h
      !-------------------------------------------
      if (ag_tmaxmin) then
        if (s% tmax1 /= rvind .and.                    &
            s% tmax  == rvind .and. modulo (ihour,6)==0) s% tmax = -99.
        if (s% tmin1 /= rvind .and.                    &
            s% tmin  == rvind .and. modulo (ihour,6)==0) s% tmin = -99.
      end if

    endif
  endif
  !--------------------------
  ! Quality of buoy location:
  !--------------------------
  if (l_buoy .and. l_mqobl) s% qobl = mqobl(is)

  !-------------------
  ! Type of data buoy:
  !-------------------
  if (l_buoy .and. l_mtodb) spti% sttyp = mtodb(is)

  !-----------------------------------------------
  ! if station name is not set, use station number
  !-----------------------------------------------
  if (spti% statid == '' .and. spti% ident > 0) &
    write(spti% statid,'(i5.5)') spti% ident

  s% z = s% zr                  ! Preset used height from reported height
  select case (spti% hd% dbkz)
  case (1)
    !-----------------------------------------------------------
    ! metar: reported pressure is truncated, add fractional part
    !+++ disabled as height/bias correction is the way to go +++
    !-----------------------------------------------------------
!   if (s% p_msl% qc == 0) then
!     ip = int (s% p_msl% o)
!     if (mod(ip,100)==0) then
!       s% p_msl% o = s% p_msl% o + 50._sp    ! German METARs.+0.5hPa
!     else
!       s% p_msl% o = s% p_msl% o + 16.925_sp ! Internat.METARs+.16..
!     endif
!   endif
  end select
  !----------------
  ! standard checks
  !----------------
  lkeep = .true.
! call init_time (spti% actual_time, yyyy, mo, dd, hh, mi)
  call check_report_1 (spti)


  if (spti% use% state <= STAT_DISMISS) lkeep = .false.

  if (l_synop .or. l_nat .or. (l_ship .and. mon_tsea)) then
    stat_time = STAT_PASSIVE
  else
    stat_time = STAT_DISMISS
  endif
  idup = 0

  if (lkeep) then
    !-------------------------
    ! check for double entries
    !-------------------------
    select case (spti% hd% obstype)
    case (OT_SYNOP, OT_DRIBU)

#if defined (__SX__)
    !--------------------------------------------------------------------
    ! Map station id to integer "hash value" for vectorizable comparisons
    !--------------------------------------------------------------------
!   spti% stat_hash = transfer (spti% statid, spti% stat_hash)
    !------------------------------------------------------
    ! Vectorized preselection of observations to compare to
    !------------------------------------------------------
    n = obs% n_spot
    k = 0; if (allocated (mask)) k = size (mask)
    if (n > k) then
       if (allocated (mask)) deallocate (mask)
       allocate (mask(max (n, 2*k)))
    end if
    n_idx = 0
    do i = 1, n
       mask(i) = obs% spot(i)% hd% obstype  == spti% hd% obstype  .and. (&
                 obs% spot(i)% stat_hash    == spti% stat_hash .or.      &
                (obs% spot(i)% ident        /= 0           .and.         &
                 obs% spot(i)% ident        == spti% ident)    .or.      &
                (obs% spot(i)% col% c% dlat == spti% col% c% dlat .and.  &
                 obs% spot(i)% col% c% dlon == spti% col% c% dlon)      )
       if (mask(i)) n_idx = n_idx + 1
    end do
    if (n_idx > 0) then
       k = 0; if (allocated (idx)) k = size (idx)
       if (n_idx > k) then
          if (allocated (idx)) deallocate (idx)
          allocate (idx(n_idx))
       end if
       idx(1:n_idx) = pack ( (/ (k, k=1,n) /), mask(1:n))
    end if
!if (n_idx > 1) write (0,*) "n_idx =", n_idx
    !------------------------------------------
    ! Now we walk only through the list of hits
    !------------------------------------------
    do k = 1, n_idx
       i = idx(k)
#else
    !-------------------------------
    ! Full scan through observations
    !-------------------------------
    do i=1,obs% n_spot
!   do i = obs% n_spot, 1, -1    ! for sorted input, reverse scan may be faster
#endif
     !---------------------------------
     ! ignore reports already dismissed
     !---------------------------------
     if (obs% spot(i)% use% state <= STAT_DISMISS) cycle
     !----------------------
     ! same observation type
     !----------------------
     chk = (obs% spot(i)% hd% obstype == spti% hd% obstype)
     if (.not. chk .and. chk_ship >= 2) then
        !-----------------------------------------------------
        ! Special duplicate check for wrongly encoded as ship?
        !-----------------------------------------------------
        chk =  obs% spot(i)% ident > 0   .and.       spti% ident > 0  .and.  &
              (obs% spot(i)% hd% codetype == 24 .or. spti% hd% codetype == 24)
     end if
     !----------------------------------------------------------------------
     ! "Same station" check for stations with valid, hashed WIGOS station ID
     !----------------------------------------------------------------------
     if (chk .and. obs% spot(i)% wsi% valid .and. spti% statid(1:1) == "_") then
        chk = obs% spot(i)% wsi == spti% wsi
     end if
     if (chk) then
      !----------------
      ! same station id
      !----------------
      if ( obs% spot(i)% hd% buf_type == spti% hd% buf_type .and. &
          (obs% spot(i)% stat_hash    == spti% stat_hash .and.    &
           obs% spot(i)% statid       /= "SHIP"          .and.    &
           obs% spot(i)% statid       == spti% statid    )   .or. &
          (obs% spot(i)% ident        /= 0      .and.             &
           obs% spot(i)% ident        == spti% ident )            ) then
       !-----------------------
       ! same observation time?
       ! (allow tolerance for
       ! differing report path)
       !-----------------------
       dt1 = minutes(obs% spot(i)% actual_time - ana_time)
       dt2 = minutes(     spti   % actual_time - ana_time)
       same_time =              dt1 == dt2          .or.     &
                   (abs (dt1 - dt2) <= chk_dbl_dt  .and.     &
                    spti% hd% dbkz  /= obs% spot(i)% hd% dbkz)
       !---------------------------------------------------
       ! SYNOP land reporting practices (c.f. WMO B/C1.1.1)
       !---------------------------------------------------
       if (l_synop) then
          same_time = same_time .or. (obs% spot(i)% hd% time == spti% hd% time)
       end if
       if (same_time) then
        !------------
        ! same 'dbkz'
        !------------
        if(obs% spot(i)% hd% dbkz     == spti% hd% dbkz           ) then
         !-------------------------
         ! same station id and time
         !-------------------------
         call load_synop (obs, obs% spot(i), s1)
         call cmp_synop (s1, s, diff, m12, m21)
!print *, "cmp_synop:", spti% statid, spti% hd% dbkz, spti% actual_time, diff, m12, m21
!        if (obs%spot(i)% corme>0 .or. spti% corme>0) then
         !-------------------------------------------------------
         ! Correction level (duplicates possible in DWD database)
         !-------------------------------------------------------
         if (obs%spot(i)% corme    /=  spti% corme  ) then
          !-------------------
          ! correction reports
          !-------------------
          if      (spti% corme > obs%spot(i)% corme) then       ! s corrects s1
            if (prt_data) print *, "CORR 1<2:", obs%spot(i)% corme, spti% corme
            write (com,'("corme = ",i0)') spti% corme
            spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
            call decr_rpt_use (obs% spot(i), CHK_CORR, comment=com)
            call merge_synop (s1, s)
            call check_store_synop (s1, spti, obs, lkeep, repl=i)
            lkeep = .false.
            s     = s1         ! Retain saved data for diagnostic printout below
          else if (spti% corme < obs%spot(i)% corme) then       ! s1 corrects s
            if (prt_data) print *, "CORR 1>2:", obs%spot(i)% corme, spti% corme
            spti% corme = obs%spot(i)% corme
            write (com,'("corme = ",i0)') spti% corme
            spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
            call decr_rpt_use (obs% spot(i), CHK_CORR, comment=com)
            call merge_synop (s, s1)
            call check_store_synop (s , spti, obs, lkeep, repl=i)
            lkeep = .false.
          else
            if (prt_data) print *, "CORR 1==2:", obs%spot(i)% corme, spti% corme
            write (com,'("corme = ",i0)') obs%spot(i)% corme
            call decr_rpt_use (obs% spot(i) ,CHK_CORRERR ,STAT_DISMISS, com)
            write (com,'("corme = ",i0)') spti% corme
            call decr_rpt_use (spti         ,CHK_CORRERR ,STAT_DISMISS, com)
            lkeep = .false.
          endif
         else ! corme
          !----------------------
          ! no correction reports
          !----------------------
          if (chk_dbl > 1 .and. spti% hd% obstype == OT_SYNOP .and. dt1 /= dt2) then
            !-----------------------------------------------------------
            ! multiple / sub-hourly(?) reports of same station for same
            ! synoptic date and time.  Keep report closest to full hour.
            !-----------------------------------------------------------
            dh1 = hours (obs% spot(i)% actual_time)
            dh2 = hours (     spti   % actual_time)
            dh1 = dh1 - nint (dh1)
            dh2 = dh2 - nint (dh2)
            if (abs (dh1) <= abs (dh2)) then
              write (com,'("synop.time, keep: ",i0)') obs% spot(i)% hd% id
              call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
              lkeep = .false.
            else
              if (ag_combine .and. .not. l_swis) then
                call combine_ag_into (s, s1)
              end if
              write (com,'("synop.time, use: ",i0)') spti% hd% id
              spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
              call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
              call check_store_synop (s, spti, obs, lkeep, repl=i)
              lkeep = .false.
            end if
          elseif (chk_dbl > 1 .and. spti% hd% obstype == OT_DRIBU) then
            !----------------------------------------------
            ! conflicting double occurence of buoys:
            ! check for possible recoding via ship template
            ! reports with location apparently rounded
            ! to 0.1 degree are assigned lower priority
            !----------------------------------------------
            slon = obs% spot(i)% col% c% dlon * 10._wp
            slat = obs% spot(i)% col% c% dlat * 10._wp
            buoy1_as_ship = (abs (slat - nint (slat)) + &
                             abs (slon - nint (slon))   ) < 2.e-4_wp
            if (prt_data) print *, "chk_dbl:buoy1/ship?", buoy1_as_ship,slat,slon
            slon = spti% col% c% dlon         * 10._wp
            slat = spti% col% c% dlat         * 10._wp
            buoy2_as_ship = (abs (slat - nint (slat)) + &
                             abs (slon - nint (slon))   ) < 2.e-4_wp
            if (prt_data) print *, "chk_dbl:buoy2/ship?", buoy2_as_ship,slat,slon
            if      (buoy1_as_ship .and. .not. buoy2_as_ship) then
              !--------------------------------
              ! drop first message, keep second
              !--------------------------------
              write (com,'("buoy location, use: ",i0)') spti% hd% id
              spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
              call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
              call check_store_synop (s, spti, obs, lkeep, repl=i)
              lkeep = .false.
            else if (buoy2_as_ship .and. .not. buoy1_as_ship) then
              !--------------------------------
              ! drop second message, keep first
              !--------------------------------
              write (com,'("buoy location, keep: ",i0)') obs% spot(i)% hd% id
              call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
              lkeep = .false.
            else if (diff) then
              write (com,'("conflict: ",i0)') spti% hd% id
              call decr_rpt_use (obs% spot(i) ,CHK_DBLERR ,STAT_DISMISS ,comment=com)
              write (com,'("conflict: ",i0)') obs% spot(i)% hd% id
              call decr_rpt_use (spti         ,CHK_DBLERR ,STAT_DISMISS ,comment=com)
              lkeep = .false.
            else
              if (prt_data) print *, "chk_dbl(BUOY):diff,m12,m21:",diff,m12,m21
              ! (handle chk_dbl=2 regression for split buoy reports)
              if (m21 .and. chk_dbl > 2) then
                !------------------------------
                ! non conflicting messages,
                ! the latter extends the former
                !------------------------------
                write (com,'("nc, use: ",i0)') spti% hd% id
                spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
                call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
                call check_store_synop (s, spti, obs, lkeep, repl=i)
                lkeep = .false.
              else if (m12 .and. chk_dbl > 2) then
                !---------------------------------
                ! non conflicting messages,
                ! the former may extend the latter
                !---------------------------------
                write (com,'("nc, keep: ",i0)') obs% spot(i)% hd% id
                call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
                lkeep = .false.
              !---------------------------------
              ! non conflicting messages, buoys:
              ! compare quality of buoy location
              !---------------------------------
!------------------------------------------------------------------
! Quality of buoy location (0 33 023):
!  0 Reliable     (location was made over two satellite passes)
!  1 Latest known (no location over the corresponding pass)
!  2 Dubious      (location made over one pass only; a second
!                  solution is possible in 5 per cent of the cases)
!  3 Missing value
!------------------------------------------------------------------
              else if (s% qobl == s1% qobl) then
                write (com,'("prev: ",i0)') obs% spot(i)% hd% id
                call decr_rpt_use (spti ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
                lkeep = .false.
              else if (abs (s% qobl) > abs (s1% qobl)) then
                write (com,'("qual,keep: ",i0)') obs% spot(i)% hd% id
                call decr_rpt_use (spti ,CHK_QI        ,STAT_DISMISS ,comment=com)
                lkeep = .false.
              else
                spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
                call decr_rpt_use (obs% spot(i) ,CHK_QI ,STAT_DISMISS)
                call check_store_synop (s, spti, obs, lkeep, repl=i)
                lkeep = .false.
              end if
            end if
            !--------------
            ! chk_dbl <= 1:
            !--------------
          elseif (diff) then
            !-----------------------------
            ! conflicting double occurence
            !-----------------------------
            write (com,'("conflict: ",i0)') spti% hd% id
            call decr_rpt_use (obs% spot(i) ,CHK_DBLERR ,STAT_DISMISS ,comment=com)
            write (com,'("conflict: ",i0)') obs% spot(i)% hd% id
            call decr_rpt_use (spti         ,CHK_DBLERR ,STAT_DISMISS ,comment=com)
            lkeep = .false.
          elseif (m21) then
            !------------------------------
            ! non conflicting messages,
            ! the latter extends the former
            !------------------------------
            write (com,'("nc, use: ",i0)') spti% hd% id
            spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
            call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
            call check_store_synop (s, spti, obs, lkeep, repl=i)
            lkeep = .false.
          elseif (m12 .or. spti% hd% obstype /= OT_DRIBU) then
            !---------------------------------
            ! non conflicting messages,
            ! the former may extend the latter
            !---------------------------------
            write (com,'("nc, keep: ",i0)') obs% spot(i)% hd% id
            call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
            lkeep = .false.
          else  ! .not.(m12.or.m21) .and. obstype == OT_DRIBU
            !---------------------------------
            ! non conflicting messages, buoys:
            ! compare quality of buoy location
            !---------------------------------
!------------------------------------------------------------------
! Quality of buoy location (0 33 023):
!  0 Reliable     (location was made over two satellite passes)
!  1 Latest known (no location over the corresponding pass)
!  2 Dubious      (location made over one pass only; a second
!                  solution is possible in 5 per cent of the cases)
!  3 Missing value
!------------------------------------------------------------------
            if (s% qobl == s1% qobl) then
              write (com,'("prev: ",i0)') obs% spot(i)% hd% id
              call decr_rpt_use (spti ,CHK_REDUNDANT ,STAT_DISMISS ,comment=com)
              lkeep = .false.
            else if (abs (s% qobl) > abs (s1% qobl)) then
              call decr_rpt_use (spti ,CHK_QI        ,STAT_DISMISS)
              lkeep = .false.
            else
              spti% use% state = min (spti% use% state, obs% spot(i)% use% state)
              call decr_rpt_use (obs% spot(i) ,CHK_QI ,STAT_DISMISS)
              call check_store_synop (s, spti, obs, lkeep, repl=i)
              lkeep = .false.
            end if
          endif
         endif   ! no corme
        endif    ! same dbkz
        !-------------------------------------------------
        ! same time and station id, prefer new BUFR format
        ! keep rain and gust information from both reports
        ! special duplicate check for autom. ship/non-ship
        !-------------------------------------------------
        if (lkeep) then
         if (prt_data) print *, "DBKZ:", obs%spot(i)% hd% dbkz, spti% hd% dbkz
         call load_synop (obs, obs% spot(i), s1)

         if (                   chk_ship >= 2                          .and. &
              obs% spot(i)% hd% codetype /=         spti% hd% codetype .and. &
             (obs% spot(i)% hd% codetype == 24 .or. spti% hd% codetype == 24)) then
           if (obs% spot(i)% hd% codetype == 24) then
              write (com,'("coloc. SHIP/other, use: ",i0)') spti% hd% id
              call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,&
                                 comment=com                                )
           else
              write (com,'("coloc. SHIP/other, use: ",i0)') obs% spot(i)% hd% id
              call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,&
                                 comment=com                                )
              lkeep = .false.
           end if

         else if (any (obs% spot(i)% hd% dbkz == kz_nat) .or. &
                      (obs% spot(i)% hd% dbkz <  10000  .and. &
                       spti%         hd% dbkz >= 10000       ))then
          if (prt_data) print *, "chk_dbl:dbkz:1,2=", obs% spot(i)% hd% dbkz, &
                                                      spti%         hd% dbkz
          if (ag_combine .and. .not. l_swis) then
            call combine_ag_into (s, s1)
          endif
          write (com,'("coloc. SYNOP/BUFR, use: ",i0)') spti% hd% id
          call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,&
                             comment=com                                )
!                            comment='collocated SYNOP/SYNOP BUFR'      )
!         call check_store_synop (s, spti, obs, lkeep, repl=i)

         else if (any (spti%         hd% dbkz == kz_nat) .or. &
                      (spti%         hd% dbkz <  10000  .and. &
                       obs% spot(i)% hd% dbkz >= 10000       )) then
          if (prt_data) print *, "chk_dbl:dbkz:2,1=", spti%         hd% dbkz, &
                                                      obs% spot(i)% hd% dbkz
          if (ag_combine .and. .not. l_swis) then
            call combine_ag_into (s1, s)
            call check_store_synop (s1, obs% spot(i), obs, lkeep, repl=i)
          endif
          write (com,'("coloc. SYNOP/BUFR, use: ",i0)') obs% spot(i)% hd% id
          call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,&
                             comment=com                                )
!                            comment='collocated SYNOP/SYNOP BUFR'      )
          lkeep = .false.
         else
          !----------------------------------------------
          ! same reporting technology,
          ! choice based on useful meteorological content
          !----------------------------------------------
          call load_synop  (obs, obs% spot(i), s1)
          call count_synop (s1, s, nsame, ndiff, only1, only2)
          if (prt_data) print *, "chk_dbl:only1,only2,ndiff=", only1, only2, ndiff
          if (only1 >= only2) then
            write (com,'("coloc. SYNOP/BUFR, use: ",i0)') obs% spot(i)% hd% id
            call decr_rpt_use (spti         ,CHK_REDUNDANT ,STAT_DISMISS ,&
                               comment=com                                )
!                              comment='collocated SYNOP/SYNOP BUFR'      )
            lkeep = .false.
          else
            write (com,'("coloc. SYNOP/BUFR, use: ",i0)') spti% hd% id
            call decr_rpt_use (obs% spot(i) ,CHK_REDUNDANT ,STAT_DISMISS ,&
                               comment=com                                )
!                              comment='collocated SYNOP/SYNOP BUFR'      )
!           call check_store_synop (s, spti, obs, lkeep, repl=i)

          end if
         endif
        endif ! lkeep

       else if (chk_dbl > 0 .and. obs% spot(i)% use% state > stat_time) then
        !-----------------------------------
        ! same station id but different time
        !-----------------------------------
!       dt1 = minutes(obs% spot(i)% actual_time - ana_time)
!       dt2 = minutes(     spti   % actual_time - ana_time)
        if (prt_data) print *, "chk_dbl:dt1,dt2=", nint (dt1), nint (dt2)

        if (l_ship .and. stat_time > STAT_DISMISS) then
          call load_synop (obs, obs% spot(i), s1)
          !--------------------------------------------------------
          ! Keep all tsea observations if ship is apparently moving
          ! Try to guess if speed of motion not reported.
          !--------------------------------------------------------
          if (     (s% nvs <= 0 .or. s1% nvs <= 0)                             .and. &
               abs (obs% spot(i)% col% c% dlat - spti% col% c% dlat) < 0.02_wp .and. &
               abs (obs% spot(i)% col% c% dlon - spti% col% c% dlon) < 0.02_wp) then
             stat_time = STAT_DISMISS
          end if
          if (prt_data) print *,"SHIP,v1,v2,state= ",spti%statid,s%nvs,s1%nvs,stat_time
        end if

        idup = 0
        if (abs(dt1) < abs(dt2)) then
          write (com,'("chk_dbl:time, keep: ",i0)') obs% spot(i)% hd% id
          call decr_rpt_use (spti         ,CHK_TIME, stat_time ,comment=com)
          if (stat_time > STAT_DISMISS) idup =  -1

        else if (abs(dt1) == abs(dt2)) then

          if (dt1 <= dt2) then
            write (com,'("chk_dbl:time, keep: ",i0)') obs% spot(i)% hd% id
            call decr_rpt_use (spti         ,CHK_TIME, stat_time ,comment=com)
            if (stat_time > STAT_DISMISS) idup =  -1
          else
            write (com,'("chk_dbl:time, use: ",i0)') spti% hd% id
            call decr_rpt_use (obs% spot(i) ,CHK_TIME, stat_time ,comment=com)
            if (stat_time > STAT_DISMISS) idup =  i
          endif

        else
          write (com,'("chk_dbl:time, use: ",i0)') spti% hd% id
          call decr_rpt_use (obs% spot(i) ,CHK_TIME, stat_time ,comment=com)
          if (stat_time > STAT_DISMISS) idup =  i

        endif

        lkeep = spti% use% state > STAT_DISMISS
        if (prt_data) print *, spti% statid, trim (com), lkeep

       endif    ! time
      endif     ! statid

      if (.not.lkeep) exit
      !-------------------------------------------
      ! two stations still active at same location
      !-------------------------------------------
      if (chk_dbl > 0                                      .and. &
          obs% spot(i)% use% state   >= STAT_ACTIVE        .and. &
          spti%         use% state   >= STAT_ACTIVE        .and. &
          obs% spot(i)% col% c% dlat == spti% col% c% dlat .and. &
          obs% spot(i)% col% c% dlon == spti% col% c% dlon       ) then
        idup = 0
        dt1 = abs(minutes(obs% spot(i)% actual_time - ana_time))
        dt2 = abs(minutes(     spti   % actual_time - ana_time))
        pr1 = get_prio (obs% spot(i)% hd)
        pr2 = get_prio (spti        % hd)
        if (prt_data) print *,"chk_dbl:loc:dt1,dt2,pr1,pr2=",dt1,dt2,pr1,pr2
        !----------------------------------------------------
        ! Handle observations reporting at non-synoptic times
        ! prefer higher quality (SYNOP>METAR), fixed location
        !----------------------------------------------------
        hourly1 = mod (dt1, 60._wp) < 1._wp
        hourly2 = mod (dt2, 60._wp) < 1._wp
        if (hourly1 .and. hourly2) then
           if (pr1 < 0 .or. pr2 < 0) then
              prefer1 = (pr1 >= pr2)
           else
              prefer1 = (dt1 < dt2 .or. (dt1 == dt2 .and. pr1 >= pr2))
           end if
        else if (hourly1 .and. pr1 >= 0) then
           prefer1 = .true.
        else if (hourly2 .and. pr2 >= 0) then
           prefer1 = .false.
        else
           prefer1 = (dt1 < dt2 .or. (dt1 == dt2 .and. pr1 >= pr2))
        end if
        if (prefer1) then
          write (com,'("chk_dbl:loc, keep: ",i0)') obs% spot(i)% hd% id
          call decr_rpt_use (spti         ,CHK_THIN ,STAT_DISMISS ,comment=com)
          lkeep = .false.
        else
          write (com,'("chk_dbl:loc, use: ",i0)') spti% hd% id
          call decr_rpt_use (obs% spot(i) ,CHK_THIN ,STAT_DISMISS ,comment=com)
!         call check_store_synop (s, spti, obs, lkeep, repl=i)
        endif
        if (prt_data) print *, spti% statid, trim (com), lkeep
      endif     ! active
     endif      ! obstype
     if (.not.lkeep) exit
    end do
    end select

    !--------------------
    ! no double occurence
    !--------------------
    if (lkeep) then
      call check_store_synop (s, spti, obs, lkeep)
      if (lkeep) then
         nkeep = nkeep + 1
         if (idup == -1) idup = obs% n_spot
         !--------------------------------------------
         ! Decrease status of non-aggregated variables
         !--------------------------------------------
         if (idup > 0 .and. stat_time > STAT_DISMISS) then
            call decr_use_nonag (obs% spot(idup), obs, STAT_DISMISS, CHK_TIME)
         end if
      end if
    endif
  endif

  if (prt_data) then
    print *, 'pe=',dace% pe,' mo_synop: station-id lat lon typ subtyp lkeep'
    write(*,'(1x,a,2x,2f15.5,2i8,l6)') spti% statid, &
         spti% col% c% dlat, spti% col% c% dlon,     &
         spti% hd% buf_type, spti% hd% buf_subtype, lkeep
    print *, 'pe=',dace% pe,' mo_synop: station-time  actual time  db-time'
    call print (spti% hd% time)
    call print (spti%     actual_time)
    call print (spti% hd% db_time)
!   print *, spti% hd% dbkz

    call print_synop (spti, s)
  endif
  !-------------------------
  ! end of loop over reports
  !-------------------------
  end do

  !-----------
  ! deallocate
  !-----------
 if (allocated ( ifd_1d  )) deallocate ( ifd_1d  )
 if (allocated ( rfd_1d  )) deallocate ( rfd_1d  )
 if (allocated ( ifd_2d1 )) deallocate ( ifd_2d1 )
 if (allocated ( rfd_2d1 )) deallocate ( rfd_2d1 )
 if (allocated ( ifd_2d2 )) deallocate ( ifd_2d2 )
 if (allocated ( rfd_2d2 )) deallocate ( rfd_2d2 )
!if (allocated ( ifd_2d3 )) deallocate ( ifd_2d3 )
!if (allocated ( rfd_2d3 )) deallocate ( rfd_2d3 )
 if (allocated ( ifd_2d4 )) deallocate ( ifd_2d4 )
 if (allocated ( rfd_2d4 )) deallocate ( rfd_2d4 )

 if (allocated ( ifd_1dm )) deallocate ( ifd_1dm )
 if (allocated ( rfd_1dm )) deallocate ( rfd_1dm )
 if (allocated ( ifd_2dm )) deallocate ( ifd_2dm )
 if (allocated ( rfd_2dm )) deallocate ( rfd_2dm )

 if (allocated ( ma1     )) deallocate ( ma1     )
 if (allocated ( mwrsa   )) deallocate ( mwrsa   )
 if (allocated ( mhosnn  )) deallocate ( mhosnn  )
 if (allocated ( mhp     )) deallocate ( mhp     )
 if (allocated ( mhobnn  )) deallocate ( mhobnn  )
 if (allocated ( mpn     )) deallocate ( mpn     )
 if (allocated ( mppp    )) deallocate ( mppp    )
 if (allocated ( mpppp   )) deallocate ( mpppp   )
 if (allocated ( mphph   )) deallocate ( mphph   )
 if (allocated ( nppp    )) deallocate ( nppp    )
!if (allocated ( np24    )) deallocate ( np24    )
 if (allocated ( nhhhn   )) deallocate ( nhhhn   )
 if (allocated ( mhosen  )) deallocate ( mhosen  )
 if (allocated ( mhawas  )) deallocate ( mhawas  )
 if (allocated ( mhawas1 )) deallocate ( mhawas1 )
 if (allocated ( mtdbt   )) deallocate ( mtdbt   )
 if (allocated ( mtdnh   )) deallocate ( mtdnh   )
 if (allocated ( muuu    )) deallocate ( muuu    )
 if (allocated ( mvv     )) deallocate ( mvv     )
 if (allocated ( mrr24   )) deallocate ( mrr24   )
 if (allocated ( mn      )) deallocate ( mn      )
 if (allocated ( mvtsu   )) deallocate ( mvtsu   )
 if (allocated ( mnh     )) deallocate ( mnh     )
 if (allocated ( nh      )) deallocate ( nh      )
 if (allocated ( mcc     )) deallocate ( mcc     )
 if (allocated ( mcc0    )) deallocate ( mcc0    )
 if (allocated ( mcc1    )) deallocate ( mcc1    )
! loop 000
 if (allocated ( mdrep   )) deallocate ( mdrep   )
 if (allocated ( mvtsu0  )) deallocate ( mvtsu0  )
 if (allocated ( mnh0    )) deallocate ( mnh0    )
!if (allocated ( mcc2    )) deallocate ( mcc2    )
 if (allocated ( nh0     )) deallocate ( nh0     )
! loop 001
! loop 002
 if (allocated ( me      )) deallocate ( me      )
 if (allocated ( nsss    )) deallocate ( nsss    )
 if (allocated ( mtn00   )) deallocate ( mtn00   )

 if (allocated ( nww     )) deallocate ( nww     )
 if (allocated ( mggtp   )) deallocate ( mggtp   )
 if (allocated ( mw1     )) deallocate ( mw1     )
 if (allocated ( mw2     )) deallocate ( mw2     )

! loop 003
 if (allocated (mdrep2   )) deallocate ( mdrep2  )
! loop 004
!if (allocated ( mhosen3 )) deallocate ( mhosen3 )
 if (allocated ( mggtp1  )) deallocate ( mggtp1  )
 if (allocated ( mrrr    )) deallocate ( mrrr    )
 if (allocated ( mrrrb   )) deallocate ( mrrrb   )

!if (allocated ( mhosen4 )) deallocate ( mhosen4 )
 if (allocated ( mggtp2  )) deallocate ( mggtp2  )
!if (allocated ( mggtp3  )) deallocate ( mggtp3  )
 if (allocated ( mtxtxh  )) deallocate ( mtxtxh  )

 if (allocated ( mggtp4  )) deallocate ( mggtp4  )
!if (allocated ( mggtp5  )) deallocate ( mggtp5  )
 if (allocated ( mtntnh  )) deallocate ( mtntnh  )

 if (allocated ( mhosen5 )) deallocate ( mhosen5 )
!if (allocated ( niw     )) deallocate ( niw     )
!if (allocated ( mtisi   )) deallocate ( mtisi   )
 if (allocated ( nggtp   )) deallocate ( nggtp   )
 if (allocated ( ndndn   )) deallocate ( ndndn   )
 if (allocated ( ndndn_s )) deallocate ( ndndn_s )
 if (allocated ( nfnfn   )) deallocate ( nfnfn   )
 if (allocated ( nfnfn_s )) deallocate ( nfnfn_s )

! loop 5
!if (allocated ( mtisi0  )) deallocate ( mtisi0  )
 if (allocated ( nggtp0  )) deallocate ( nggtp0  )
!if (allocated ( nmwgd   )) deallocate ( nmwgd   )
 if (allocated ( nfxgu   )) deallocate ( nfxgu   )

! loop 6
 if (allocated ( mggtp7  )) deallocate ( mggtp7  )
 if (allocated ( mlwr    )) deallocate ( mlwr    )
 if (allocated ( mglsr   )) deallocate ( mglsr   )
 if (allocated ( mdsrh   )) deallocate ( mdsrh   )

! buoy
!
!if (allocated (lmtodb   )) deallocate ( lmtodb  )
 if (allocated (mtodb    )) deallocate ( mtodb   )

!if (allocated (lmqbst   )) deallocate ( lmqbst  )
!if (allocated (lmqobl   )) deallocate ( lmqobl  )
 if (allocated (mqobl    )) deallocate ( mqobl   )
!if (allocated (lmlqc    )) deallocate ( lmlqc   )

  contains
    function valid (entry)
    real(sp) ,intent(in) :: entry
    logical              :: valid
      valid = entry >= 0. .and. entry /= rvind
    end function valid

    subroutine combine_ag_into (s, s1)
      type(t_synop), intent(inout) :: s
      type(t_synop), intent(in)    :: s1
      !----------------------------------------------
      ! Combine aggregated values for same valid date
      !----------------------------------------------
      if (.not. valid(s% rr1  ) .and. s1% rr1   >= 0.) s% rr1   = s1% rr1
      if (.not. valid(s% rr3  ) .and. s1% rr3   >= 0.) s% rr3   = s1% rr3
      if (.not. valid(s% rr6  ) .and. s1% rr6   >= 0.) s% rr6   = s1% rr6
      if (.not. valid(s% rr12 ) .and. s1% rr12  >= 0.) s% rr12  = s1% rr12
      if (.not. valid(s% rr24 ) .and. s1% rr24  >= 0.) s% rr24  = s1% rr24
      if (.not. valid(s% lwr1 ) .and. s1% lwr1  >= 0.) s% lwr1  = s1% lwr1
      if (.not. valid(s% gsr1 ) .and. s1% gsr1  >= 0.) s% gsr1  = s1% gsr1
      if (.not. valid(s% gsr3 ) .and. s1% gsr3  >= 0.) s% gsr3  = s1% gsr3
      if (.not. valid(s% gsr6 ) .and. s1% gsr6  >= 0.) s% gsr6  = s1% gsr6
      if (.not. valid(s% dsr1 ) .and. s1% dsr1  >= 0.) s% dsr1  = s1% dsr1
      if (.not. valid(s% dsr3 ) .and. s1% dsr3  >= 0.) s% dsr1  = s1% dsr3
      if (.not. valid(s% dsr6 ) .and. s1% dsr6  >= 0.) s% dsr1  = s1% dsr6
      if (.not. valid(s% gust1) .and. s1% gust1 >= 0.) s% gust1 = s1% gust1
      if (.not. valid(s% gust3) .and. s1% gust3 >= 0.) s% gust3 = s1% gust3
      if (.not. valid(s% gust6) .and. s1% gust6 >= 0.) s% gust6 = s1% gust6
      if (.not. valid(s% tmax1) .and. s1% tmax1 >= 0.) s% tmax1 = s1% tmax1
      if (.not. valid(s% tmin1) .and. s1% tmin1 >= 0.) s% tmin1 = s1% tmin1
      ! potential issues with 12-hourly data between TAC and BUFR
      if (.not. valid(s% tmax ) .and. valid(s1% tmax)) s% tmax  = s1% tmax
      if (.not. valid(s% tmin ) .and. valid(s1% tmin)) s% tmin  = s1% tmin
    end subroutine combine_ag_into

  end subroutine read_synop_netcdf
!------------------------------------------------------------------------------
  !----------------------------------------
  ! Derive priority from header information
  !----------------------------------------
  elemental function get_prio (head) result (prio)
    type(t_head), intent(in) :: head
    integer                  :: prio
    select case (head% dbkz)
    !--------------------------------------
    ! Observations at mostly synoptic times
    !--------------------------------------
    case (10000:)               ! BUFR SYNOP, SYNOP national; SHIP, BUOY
       prio = 2
    case (0,5,128,385)          ! TAC  SYNOP, SYNOP national; BUOY
       prio = 1
    !-----------------------------------------
    ! Observations at possibly asynoptic times
    !-----------------------------------------
    case (131)                  ! METAR USA (recoded as SYNOP)
       prio = 0
    case (1)                    ! METAR (report completeness, accuracy)
       prio = -1
    case (256,384)              ! TAC  SHIP (varying location)
       prio = -2
    case default                ! unknown
       prio = -3
    end select
  end function get_prio
!------------------------------------------------------------------------------
!   subroutine get_dat(ncid, ito1  , rto1  , ito2  , rto2  , charid, varid, ty_fto, ty_ffr, &
!                      ndim, ifd_1d, rfd_1d, ifd_2d, rfd_2d, i_dft, r_dft, l_var,           &
!                      start)
  subroutine get_dat(ncid, ito1  , rto1  , ito2  , rto2  , charid, varid, ty_fto, &
                     ndim, ifd_1d, rfd_1d, ifd_2d, rfd_2d, i_dft, r_dft, l_var,   &
                     start)

  integer     ,intent(in)                    :: ncid    ! NetCDF file handle

  integer     ,intent(inout), dimension(:)   :: ito1    ! 1dim integer field (to)
  real        ,intent(inout), dimension(:)   :: rto1    ! 1dim real    field (to)
  integer     ,intent(inout), dimension(:,:) :: ito2    ! 2dim integer field (to)
  real        ,intent(inout), dimension(:,:) :: rto2    ! 2dim real    field (to)
  character(*),intent(in   )                 :: charid  ! netCDF variable name
  integer     ,intent(  out)                 :: varid   ! netCDF variable id
  integer     ,intent(in   )                 :: ty_fto  ! type of field (to) ; 1 integer ; 2 real
!  integer     ,intent(in   )                 :: ty_ffr  ! type of field (from)
  integer     ,intent(in   )                 :: ndim    ! dimension of field (to)
  integer     ,intent(inout), dimension(:)   :: ifd_1d  ! 1dim integer field (from)
  real        ,intent(inout), dimension(:)   :: rfd_1d  ! 1dim integer field (from)
  integer     ,intent(inout), dimension(:,:) :: ifd_2d  ! 2dim integer field (from)
  real        ,intent(inout), dimension(:,:) :: rfd_2d  ! 2dim integer field (from)
  integer     ,intent(in   )                 :: i_dft   ! default value for integer field (to)
  real        ,intent(in   )                 :: r_dft   ! default value for real    field (to)
  logical     ,intent(  out)                 :: l_var   ! logical for variable exists and any values defined
  integer     ,intent(in   )                 :: start   ! parameter for nf90_get_var


  integer              :: status         ! NetCDF status variable
  integer              :: itype          ! NetCDF type of variable
  integer              :: ndim_aux       ! NetCDF number of dimensions
  integer              :: ty_ffr         ! type of field
  character(len=120)   :: msg = ''

  l_var = .FALSE.

  status = nf90_inq_varid (ncid,trim( charid ),  varid )
  if (status == nf90_noerr) then
    status = nf90_inquire_variable (ncid, varid, xtype=itype, ndims=ndim_aux)
  end if

  if (status == nf90_noerr) then
    l_var = .true.
    if (ndim /= ndim_aux) then
      write(msg, '(A,I1,A,I1,A)') 'unexpected number of dimensions for variable '//&
           trim(charid)//': ',ndim,' expected, but ',ndim_aux,' in NetCDF file.'
      call finish('get_dat', msg)
    end if
    if (itype==nf90_int .or. itype==nf90_byte .or. itype==nf90_short) then
      ty_ffr = TY_I
    else if (itype==nf90_float .or. itype==nf90_double) then
      ty_ffr = TY_F
    end if

    !--------------
    ! 1 dimensional
    !--------------
    if ( ndim == 1 ) then
      if ( ty_ffr == TY_I ) then
        status = nf90_get_var (ncid, varid , ifd_1d ,start=(/start/))
        l_var  = any( ifd_1d /= imissing )
        if ( ty_fto == TY_I ) then
          where ( ifd_1d == imissing )
            ito1 = i_dft
          elsewhere
            ito1 = ifd_1d
          end where
        elseif ( ty_fto == TY_F ) then
          where ( ifd_1d == imissing )
            rto1 = r_dft
          elsewhere
            rto1 = ifd_1d
          end where
        end if
      else if ( ty_fto == TY_F ) then
        status = nf90_get_var (ncid, varid , rfd_1d ,start=(/start/))
        l_var  = any( rfd_1d /= rmissing )
        where ( rfd_1d == rmissing )
          rto1 = r_dft
        elsewhere
          rto1 = rfd_1d
        end where
      end if
    !--------------
    ! 2 dimensional
    !--------------
    elseif ( ndim == 2 ) then
      if ( ty_ffr == TY_I ) then
        status = nf90_get_var (ncid, varid , ifd_2d ,start=(/1,start/))
        l_var  = any( ifd_2d /= imissing )
        if ( ty_fto == TY_I ) then
          where ( ifd_2d == imissing )
            ito2 = i_dft
          elsewhere
            ito2 = ifd_2d
          end where
        elseif ( ty_fto == TY_F ) then
          where ( ifd_2d == imissing )
            rto2 = r_dft
          elsewhere
            rto2 = ifd_2d
          end where
        end if
      else if ( ty_fto == TY_F ) then
        status = nf90_get_var (ncid, varid , rfd_2d ,start=(/1,start/))
        l_var  = any( rfd_2d /= rmissing )
        where ( rfd_2d == rmissing )
          rto2 = r_dft
        elsewhere
          rto2 = rfd_2d
        end where
      end if
    end if
  else
  ! if not defined , preset field with default
    if ( ndim == 1 ) then
      if ( ty_fto == TY_I ) then
        ito1 = i_dft
      else if ( ty_fto == TY_F ) then
        rto1 = r_dft
      endif
    else if ( ndim == 2 ) then
      if ( ty_fto == TY_I ) then
        ito2 = i_dft
      else if ( ty_fto == TY_F ) then
        rto2 = r_dft
      endif
    endif
  endif
  end subroutine get_dat
!------------------------------------------------------------------------------
  !--------------------------------------------
  ! Decrease status of non-aggregated variables
  !--------------------------------------------
  subroutine decr_use_nonag (spot, obs, status, check)
    type(t_spot) ,intent(in)    :: spot   ! report header meta data
    type(t_obs)  ,intent(inout) :: obs    ! data of all observations
    integer      ,intent(in)    :: status ! new state
    integer      ,intent(in)    :: check  ! check causing new state

    integer :: i

    do i = 1, spot% o% n
       select case (obs% varno(spot% o% i+i))
       case (VN_P, VN_PS)
          if (.not. mon_ps) &
             call decr_use (obs% body (spot% o% i+i)% use, status, check=check)
       case (VN_T, VN_TD, VN_T2M, VN_TD2M)
          if (.not. mon_t2) &
             call decr_use (obs% body (spot% o% i+i)% use, status, check=check)
       case (VN_Z, VN_RH, VN_RH2M, VN_U, VN_V, VN_U10M, VN_V10M, &
             VN_FF, VN_DD, VN_N, VN_NH, VN_CEIL, VN_N_L, VN_N_M, VN_N_H, VN_VV,&
             VN_WW, VN_GCLG, VN_ICLG, VN_SDEPTH)
          call decr_use (obs% body (spot% o% i+i)% use, status, check=check)
       end select
    end do
  end subroutine decr_use_nonag
!------------------------------------------------------------------------------
  subroutine cloudcover_synop (s, tcc, nlay, mvtsu, mnh, mcc, nh)
    type(t_synop), intent(inout) :: s           ! SYNOP derived type
    integer,       intent(in)    :: tcc         ! Total cloud cover [%] (020010)
    integer,       intent(in)    :: nlay        ! Number of cloud layers
    integer,       intent(in)    :: mvtsu(:)    ! Vertical signif. (008002)
    integer,       intent(in)    :: mnh  (:)    ! Cloud amount     (020011)
    integer,       intent(in)    :: mcc  (:)    ! Cloud type       (020012)
    real,          intent(in)    :: nh   (:)    ! Height of cloud base [m]
    !-------------------------------------------------------------------------
    ! Interpret cloud data in SYNOP reports
    !-------------------------------------------------------------------------
    integer             :: i, k
    real                :: h                    ! Cloud height above msl
    real                :: octa(4)              ! Low/mid/high/total cloud cover
    integer             :: clcx(4)              ! Low/mid/high/total cloud cover
    integer             :: clc (nlay)           ! Temporary cloud cover
    integer             :: ix  (nlay)           ! Auxiliary index
    integer             :: vs                   ! Vertical significance
    real(sp), parameter :: H_LOW = 2000._sp     ! Limit between low and mid
    real(sp), parameter :: H_MID = 7000._sp     ! Limit between mid and high
    real(sp), parameter :: N_FEW = 1.5_sp       ! Average amount for FEW (1-2)
    real(sp), parameter :: N_SCT = 3.5_sp       ! Average amount for SCT (3-4)
    real(sp), parameter :: N_BKN = 6.0_sp       ! Average amount for BKN (5-7)
    real(sp), parameter :: N_OVC = 8.0_sp       ! Overcast

    if (nlay > size(mvtsu)) return              ! Should never happen

    !---------------------------------------------------------------------------
    ! Sky obscured by fog and/or other meteorological phenomena (020010, note 5)
    !---------------------------------------------------------------------------
    if (tcc == 113) then
       s% n    = 100
       s% clcl = N_OVC
       return
    end if

    s% n = tcc
    octa = -1
    if (tcc >= 0) octa(4) = octa_percent (real (tcc, wp))

    ! if total cloud cover 'N' is zero, then cloudless (MNH, NH not reported)
    if (octa(4) == 0) then
       s% clcl = 0._sp
       s% clcm = 0._sp
       s% clch = 0._sp
       return
    end if

    ! clear sky (vert. signif. = value not applicable) and total cloud undefined
    ! B/C1.4.4.2(e) requires this, but too many bad data, so we don't handle it.
!   if (tcc < 0 .and. mvtsu(1) == 62) then
!      s% clcl = -1._sp
!      s% clcm = -1._sp
!      s% clch = -1._sp
!      return
!   end if

    !-------------------------------------------------------
    ! Single layer: treat clear sky and fog/height undefined
    !-------------------------------------------------------
    if (nlay == 1 .and. nh(1) == rvind) then
       if      (mnh(1) == 0) then
          s% clcl = 0._sp
          return
       else if (mnh(1) == 9) then
          s% clcl = N_OVC
          return
       end if
    end if

    clcx = -1
    clc  = -1
    ix   = -1

    !-----------------------
    ! Check cloud type first
    ! B/C1.4.4.5 and 020012
    !-----------------------
    if (btest (chk_cldinfo, 1)) then
       if (mcc(1) == 30) clcx(1) = 0   ! No C_L clouds
       if (mcc(2) == 20) clcx(2) = 0   ! No C_M clouds
       if (mcc(3) == 10) clcx(3) = 0   ! No C_H clouds
    end if

    !-------------------
    ! Check cloud height
    !-------------------
    do i = 1, nlay
       vs = mvtsu(i)
       if (vs == -1) then
          if (i == 1 .and. btest (chk_cldinfo, 0)) vs = 0
       end if
       select case (vs)
       case (7)           ! Low cloud
          ix(i) = 1
       case (8)           ! Middle cloud
          ix(i) = 2
       case (9)           ! High cloud
          ix(i) = 3
       case (0,1:5,21:24) ! 1st, 2nd, ... cloud layer
          h = nh(i)       ! Height of base of cloud above surface (0 20 013)
                          ! Should we check height ordering (B/C1.4.5.1.3)?
!         if (h /= rvind .and. s% zr% o /= rvind) h = h + s% zr% o
          if      (h <  H_LOW) then
             ix(i) = 1
          else if (h <  H_MID) then
             ix(i) = 2
          else if (h /= rvind) then
             ix(i) = 3
          else
             ix(i) = -1
          end if
!      case (20)  ! No clouds detected by system (c.f. B/C1.4.5.2.1h)
!                 ! * Do not use because of broken Swedish automatic stations!
!         if (mnh(i) == 0) clcx(4) = 0
!      case (0)   ! FM12/FM13 rules apply
!      case (6)   ! Clouds not detected below following height(s)
!      case (62)  ! Value not applicable
!      case (:-1) ! Missing value
       end select
    end do

    !---------------------------------------------
    ! Check cloud amount for classes we can handle
    !---------------------------------------------
    do i = 1, nlay
       select case (mnh(i))
       case (0)
          if (nh(i) == rvind) then
             clc(i) = 0     ! see B/C1.4.4.4.2
          else
             ix (i) = -1    ! Dismiss 0 octas for defined cloud base height!
          end if
       case (1:9,11:13)
          clc(i) = mnh(i)
       case default
          ix (i) = -1
       end select
    end do

    if (prt_data) then
       write(*,*) "### cloudcover_synop: n=", s% n
       write(*,*) "###  nh =",  nh(1:nlay)
       write(*,*) "###  ix =", ix (1:nlay)
       write(*,*) "### mnh =", mnh(1:nlay)
       write(*,*) "### clc =", clc(1:nlay)
    end if

    !--------------------------------
    ! Merge information within layers
    ! keep highest amount (020011)
    !--------------------------------
    do i = 1, nlay
       k = ix(i)
       if (k > 0) then
          select case (clc(i))
          case (0:9)
             clcx(k) = max (clcx(k), clc(i))
             if (clcx(k) == 9) clcx(4) = 8      ! Sky obscured by fog etc.
          case (13)     ! Few       (1-2 octas)
             if (clcx(k) < 2) then
                clcx(k) = 13
             end if
          case (11)     ! Scattered (3-4 octas)
             if (clcx(k) < 4 .or. clcx(k) == 13) then
                clcx(k) = 11
             end if
          case (12)     ! Broken    (5-7 octas)
             if (clcx(k) < 6 .or. clcx(k) == 11 .or. clcx(k) == 13) then
                clcx(k) = 12
             end if
          end select
       end if
    end do

    !------------------------------------------
    ! Canonicalize for categorical verification
    ! convert to (representative) octas
    !------------------------------------------
    do i = 1, 3
       select case (clcx(i))
       case (0:8)
          octa(i) = real (clcx(i), sp)
       case (9)             ! Sky obscured by fog etc. (8 octas)
          octa(i) = N_OVC
       case (11)            ! Scattered (3-4 octas)
          octa(i) = N_SCT
       case (12)            ! Broken    (5-7 octas)
          octa(i) = N_BKN
       case (13)            ! Few       (1-2 octas)
          octa(i) = N_FEW
       end select
    end do
    if (clcx(4) >= 0) then
       octa(4) = real (clcx(4), sp)
    end if

    if (prt_data) then
       write(*,*) "### clcx:", clcx(:)
       write(*,*) "### octa:", octa(:), "#1"
    end if

    !------------------------------------
    ! Recover low/medium/high "clear sky"
    ! when no low cloud observed/reported
    !------------------------------------
    if   (octa(1) <= 0._sp) then
      if (octa(4) == 0._sp .and. octa(3) < 0._sp) octa(3) = 0._sp
      if (octa(3) >= 0._sp .and. octa(2) < 0._sp) octa(2) = 0._sp
      if (octa(2) >= 0._sp .and. octa(1) < 0._sp) octa(1) = 0._sp
    end if

    !-------------------------------------
    ! Plausibility checks for cloud layers
    ! (Not to be used for manned stations;
    ! we also do not enforce B/C1.4.5.1.1)
    !-------------------------------------
!   if ((octa(2) >= 0._sp .and.  octa(2) <  octa(1)) .or. &
!       (octa(2) >  0._sp .and.  octa(2) == octa(1))      ) octa(2) = -1._sp
!   if ( octa(3) >= 0._sp .and. (octa(3) <  octa(1)  .or. &
!                                octa(3) <  octa(2)      )) octa(3) = -1._sp
!   if ( octa(3) >  0._sp .and. (octa(3) == octa(1)  .or. &
!                                octa(3) == octa(2)      )) octa(3) = -1._sp

    if (prt_data) then
       write(*,*) "### octa:", octa(:), "#2"
    end if

    !--------------------------------------
    ! Cross-check against total cloud cover
    !--------------------------------------
    if (octa(4) >= 0._sp) octa(4) = maxval (octa)

    if (octa(1) >= 0._sp) s% clcl = octa(1)
    if (octa(2) >= 0._sp) s% clcm = octa(2)
    if (octa(3) >= 0._sp) s% clch = octa(3)
    if (octa(4) >= 0._sp) s% n    = max (s% n, nint (octa(4) * 12.5_sp))

  end subroutine cloudcover_synop
!------------------------------------------------------------------------------
  subroutine cloudcover_synop_lam (s, tcc, nlay_in, mvtsu, mnh_in, mcc, nh    &
                                  ,vis, ww, statid, hhmmss )
    type(t_synop), intent(inout) :: s          ! SYNOP derived type
    integer,       intent(in)    :: tcc        ! total cloud cover [%] [020010]
    integer,       intent(in)    :: nlay_in    ! number of cloud layers (>=1)
                                               !-> following arrays: indices:
                                               !    1     : general cloud group
                                               !    2-nlay: indiv. cl. layer gr.
                                               !   missing val.: -1 resp. rvind:
    integer,       intent(in)    :: mvtsu (:)  ! vertical significance [008002]
    integer,       intent(in)    :: mnh_in(:)  ! cloud amount          [020011]
    integer,       intent(in)    :: mcc   (:)  ! cloud type            [020012]
    real,          intent(in)    :: nh    (:)  ! height of cloud base AGL   [m]
    real,          intent(in)    :: vis        ! visibility                 [m]
    integer,       intent(in)    :: ww         ! present weather       [020003]
    character(len=*), intent(in) :: statid     ! station id (for printing)
    integer,       intent(in)    :: hhmmss     ! hour, minute, sec. of obs
    !-------------------------------------------------------------------------
    ! Interpretation of cloud data in SYNOP reports
    ! ---------------------------------------------
    ! - this mostly agrees well with the COSMO obs operators (for KENDA)
    ! - code tables given by a '[6-digit number]' are from:
    !     WMO manual on Codes, Vol. I.2, Annex II, Part B - Binary Codes:
    !     FM 94 BUFR - BUFR/CREX Table B
    ! - code tables '[A- followed by a 4-digit number]' are from:
    !     WMO manual on Codes, Vol. I.1, Annex II, Part A - Alphanumeric Codes:
    !     Section C, Code tables
    !-------------------------------------------------------------------------
    integer             :: nlay                ! number of reported cloud layers
    integer             :: inh  (nlay_in)      ! nh    , adjusted, miss.val.=-1
    integer             :: knh  (nlay_in)      ! nh    , reported
    integer             :: mnh  (nlay_in)      ! mnh_in, reported
    integer             :: kvs  (nlay_in)      ! mvtsu , missing val = MIS_SG
    real(sp)            :: r_mnh(nlay_in)      ! accurate mnh    , miss.val.=-1
    real(sp)            :: clcd (4)            ! low/mid/high/total cloud amount
                                               !  derived from general cld group
    real(sp)            :: clcx (4)            ! low/mid/high/total cloud amount
                                               !  derived from indiv. cld group
    real(sp)            :: clcml               ! cld amount of lowest mid cloud
    integer             :: ceilx               ! cloud ceiling [m] above ground
    integer             :: ceilg               ! cloud ceiling [m] from GCG
    integer             :: cbhx                ! cloud base height [m] AGL
    integer             :: kfog                ! fog (probability) indicator:
                                               !   = 3: yes   : reported fog
                                               !   = 2: likely: low visibility
                                               !   = 1: doubt : low vis + precip
                                               !   = 0: no    : (no indication)
                                               !   =-1: dust, sand --> no fog
                                               !        even if vis is low
    integer             :: i, k                ! indices
    integer             :: ivsg, icbh, icil ,& ! for printing
                           ivis, icbx, icix ,& !
                           ivsgg               !
    logical             :: lcil                ! ceiling reported
    logical             :: lcc1, lcc2          ! results of consistency checks
    integer             :: igcg                ! general cloud group
    integer             :: iclg (nlay_in)      ! individual cloud layers
    integer             :: kcc  (3)            ! mcc [A-0509, A-0513, A-0515]
    integer             :: nhkey               ! cloud base height  [A-1600]
    integer             :: nclct               ! total cloud amount [A-2700]
    integer             :: hhmm                ! hour, minute of obs
    integer , parameter :: H_LOW   = 1999      ! Limit between low and mid
    integer , parameter :: H_MID   = 6999      ! Limit between mid and high
    real(sp), parameter :: N_FEW   = 1.5_sp    ! Average amount for FEW (1-2)
    real(sp), parameter :: N_SCT   = 3.5_sp    ! Average amount for SCT (3-4)
    real(sp), parameter :: N_BKN   = 6.0_sp    ! Average amount for BKN (5-7)
    real(sp), parameter :: N_OVC   = 8.0_sp    ! Overcast
    real(sp), parameter :: TOL     = 0.001_sp  ! tolerance
    real    , parameter :: VFOGLIM = 500._wp   ! visibility lim. [m] below which
    !                      low cloud (fog) is assumed to exist in case of precip
    !   bit positions in general cloud group (NBP_SG from COSMO operators):
    integer , parameter :: NBP_CTH =  0        ! type of high cloud
    integer , parameter :: NBP_CTM =  4        ! type of mid-level cloud
    integer , parameter :: NBP_CTL = 12        ! type of low cloud
    integer , parameter :: NBP_HK  =  8        ! cloud base height
    integer , parameter :: NBP_NH  = 16        ! amount of low (else mid) cloud
    integer , parameter :: NBP_N   = 20        ! total cloud amount
    integer , parameter :: NBP_SG  = 24        ! vertical significance
    !   bit positions in individual cloud layer group (differs from COSMO op.):
    integer , parameter :: IBP_SG  =  0        ! vertical significance
    integer , parameter :: IBP_N   =  6        ! cloud amount
    integer , parameter :: IBP_H   = 10        ! cloud base height
    !   no. of bits occupied in cloud groups:
!   integer , parameter :: NBO_CT  =  4        ! cloud type    [A-0509 - A-0515]
!   integer , parameter :: NBO_N_  =  4        ! cloud amount           [A-2700]
!   integer , parameter :: NBO_HK  =  4        ! cloud base height      [A-1600]
!   integer , parameter :: NBO_SG  =  6        ! vertical significance  [008002]
!   integer , parameter :: NBO_H   = 11        ! cloud base height        [10*m]
    !   all bits set to indicate missing value:
    integer , parameter :: MIS_CT  = 15        ! 2**NBO_CT - 1
    integer , parameter :: MIS_N_  = 15        ! 2**NBO_N_ - 1
    integer , parameter :: MIS_HK  = 15        ! 2**NBO_HK - 1
    integer , parameter :: MIS_H   = 2047      ! 2**NBO_H  - 1
    integer , parameter :: MIS_SG  = 63        ! 2**NBO_SG - 1
    integer , parameter :: MIS_GCG = 1073741823! 2**(3*NBO_CT + NBO_HK +
                                               !     2*NBO_N_ + NBO_SG) - 1
    integer , parameter :: MIS_ICG = 2097151   ! 2**(NBO_SG + NBO_N_ + NBO_H) -1
    !-------------------------------------------------------------------------

    if (nlay_in > size(mnh_in)) return            ! should never happen
    nlay = nlay_in

    ! pre-set flag for fog based on visibility and present weather [020003]
    ! ---------------------------------------------------------------------
    kfog = 0
    if ((vis < VFOGLIM) .and. (vis >= 0._wp))  kfog = 2
    !   present weather is fog (excluding sky visible)
    if (         ((ww >=  43) .and. (ww <=  49) .and. (mod( ww,2 ) == 1))      &
            .or. ((ww >= 132) .and. (ww <= 135)) .or. (ww == 130)              &
            .or. ((ww >= 247) .and. (ww <= 249))) then
       kfog = 3
    !   dust(storm), sand(storm), smoke, (blowing) snow, haze, volcanic ash, ...
    !       04-09 haze, dust,sand, smoke
    !       30-39 duststorm, sandstorm, drifting or blowing snow
    !       104, 105, 127-129 haze, smoke, dust, snow
    !       204, 206-208, 210 volcanic ash, dust haze, blowing spray, snow haze
    !       230, 239 duststorm, sandstorm, blowing snow
    elseif (     ((ww >=   4) .and. (ww <=   9))                               &
            .or. ((ww >=  30) .and. (ww <=  39))                               &
            .or. ((ww >= 104) .and. (ww <= 105))                               &
            .or. ((ww >= 127) .and. (ww <= 129)) .or. (ww == 230)              &
            .or. ((ww >= 204) .and. (ww <= 210)) .or. (ww == 239)) then
       kfog = -1
    !   in case of heavy precip, it is unclear whether
    !   (moderately) low visibility is due to precip or fog
    elseif (      (kfog == 2) .and. (vis > 100._wp)                            &
            .and. (     ((ww >=  50) .and. (ww <=  99))                        &
                   .or. ((ww >= 150) .and. (ww <= 199)))) then
       select case (ww)
       case (54, 55, 57, 59, 64, 65, 67, 69, 74, 75, 81, 82, 84, 86, 88, 90, 97&
            ,98, 99,153,156,158,163,166,168,173,176,183,184,187,189,195,196)
          kfog = 1
       end select
    endif

    ! cloud base height 'nh': convert into integer 'inh' for convenience
    ! ------------------------------------------------------------------
    ! (nh: precision = 10[m], 11 bit in BUFR --> largest possible value = 20470)
    inh (1:nlay) = -1
    do i = 1 , nlay
       if (      (nh(i) > -0.1_wp) .and. (nh(i) <= 20470.1_wp)                 &
           .and. (nh(i) /= rvind))  inh(i) = min( nint(nh(i)) , CBH_CLR )
    enddo
    !   after this, missing value is -1 for 'inh' ('nh' is not used any more)

    ! cloud amounts 'mnh': convert to real 'r_mnh' and canonicalize
    ! -------------------------------------------------------------
    !    (cloud amounts have to be reals (only) due to N_SCT, N_FEW)
    r_mnh (1:nlay) = -1._sp
    do i = 1, nlay
       if ((mnh_in(i) >= 0) .and. (mnh_in(i) <= 13))  r_mnh(i) = real(mnh_in(i))
       !   FM 15 Metar, regulation 15.9.1.1, for [020011]:
       if (mnh_in(i) == 12)  r_mnh(i) = N_BKN  ! 12: 'broken'    --> 5 - 7 octas
       if (mnh_in(i) == 11)  r_mnh(i) = N_SCT  ! 11: 'scattered' --> 3 - 4 octas
       if (mnh_in(i) == 13)  r_mnh(i) = N_FEW  ! 13: 'few'       --> 1 - 2 octas
       if (mnh_in(i) == 10)  r_mnh(i) = -1._sp ! 10: 'part. obscured (e.g. fog)'
    enddo
    !   note: after this point, -1. <= r_mnh <= 9.

    ! total cloud cover
    ! -----------------
    !   conversion to octa
    clcd = -1._sp
    if (tcc == 113)  clcd(4) = 9._sp
    if ((tcc >= 0) .and. (tcc <= 100))  clcd(4) = octa_percent (real (tcc, wp))

    !   preparation for packed cloud info
    knh(:) = inh(:)
    mnh(:) = nint( r_mnh(:) )
    nclct  = nint( clcd(4) )

    !   zero total cloud: consistency checks
    if (nint(clcd(4)) == 0) then
       !   inconsistency if cloud base height (CBH) < CBH_CLR, except:
       !   - [A-1600]: CBH = 2500 m also for zero cloud, hence not inconsistent
       !   - [008002]: mvtsu = 6 = 'clouds not detected below following height'
       !               (-> France)
       if (     ((inh(1) > 0) .and. (inh(1) < 2500))                           &
           .or. (      (inh(1) > 2500) .and. (inh(1) < CBH_CLR)                &
                 .and. (mvtsu(1) /= 6))) then
          clcd (4)  = -1._sp
          inh  (1)  = -1
       !  CBH not reported or zero or >= 2500 m: assume cloud-free is correct
       else
          inh  (1)  = CBH_CLR
       endif
       !  inconsistency:  non-zero MNH(1) (= low or mid-level cloud) reported
       if (r_mnh(1) >= 1._sp) then
          clcd (4)  = -1._sp
          r_mnh(1)  = -1._sp
          inh  (1)  = -1
       endif
    endif

    ! -----------------------------------------------------------------------
    ! derive low, mid-level, high cloud amount from general cloud group (GCG)
    ! -----------------------------------------------------------------------
    !   WMO Manual on Codes, Vol I.1, Annex II, Alphanumeric codes: 12.2.7.3:
    !       The coding of CL, CM and CH clouds (for 'MHN' + cloud types) shall
    !       be as specified in the Internatl. Cloud Atlas (WMO-No. 407), Vol I.
    !       --> This is according to cloud type, not height (above ground).
    !       --> CBH is given preference over vert. significance (or cloud type).
    !   In ICON, layer of low  cloud is at p > 0.8 ps,
    !            layer of high cloud is at p < 0.4 ps, mid-level cloud in betw.
    !   In standard atmosphere, p = 0.8 ps is at ~ 1940 m,
    !                           p = 0.4 ps is at ~ 7180 m;
    !   hence height limits of <= H_LOW = 1999 m and <= H_MID = 6999 m are used.
    !   Keep in mind that for FM12, cloud base height = 2500 m means >= 2500 m
    !     or no cloud (see [A-1600]).
    !   Note: Only 1 cloud layer (low, mid, or high) derived from GCG can be >0.
    !         Total cloud is used here only if it indicates zero cloud or fog,
    !         or for high cloud if low + mid-level cloud 'mnh(1)' is zero.
    clcd (1)  = r_mnh(1)
    clcd (2)  = r_mnh(1)
    clcd (3) = -1._sp
    !   zero total cloud: (mnh(1) > 0 excluded, see above --> mnh(1) = 0 or -1)
    if (nint(clcd(4)) == 0) then
       clcd (1)  = 0._sp
       clcd (2)  = 0._sp
       clcd (3)  = 0._sp
    !   non-zero total cloud and zero mnh(1) --> high cloud = total cloud
    !   (note: this may lead to cloud layers slightly inconsistent with model
    !          equivalents because reported low / mid-level cloud is based on
    !          cloud type, not on cloud base / top height as in the model)
    elseif (      (nint(clcd(4)) >= 1) .and. (nint(clcd(4)) /= 9)              &
            .and. (nint(r_mnh(1)) == 0)                                        &
            .and. ((inh(1) >= 2500) .or. (inh(1) == -1))) then
       clcd (3)  = clcd(4)
    !   cloud amount = 9 : 'sky obscured by fog or other meteorol. phenomena'
    !   assume fog / low cloud with overcast sky (8 octas) if both
    !   - sky obscured ('MNH'=9, or 'N'=9 (then 'MNH', 'NH' not reported)), and
    !   - fog flag is set or vertical significance = ceiling
    elseif (      ((nint(r_mnh(1)) == 9) .or. (nint(clcd(4)) == 9))            &
            .and. ((kfog >= 3) .or. ((mvtsu(1) == 5) .and. (kfog /= -1)))) then
       clcd (4)  = N_OVC
       clcd (1)  = -1._sp
       clcd (2)  = -1._sp
       if (((inh(1) >= 0) .and. (inh(1) <= H_LOW)) .or. (kfog >= 2)) then
         clcd (1)  = N_OVC
       elseif ((inh(1) >  H_LOW) .and. (inh(1) <= H_MID)) then
         clcd (1)  = 0._sp
         clcd (2)  = N_OVC
       endif
    !   otherwise if sky obscured then do not trust CBH or vert. significance
    !     (and hope to get info from individual cloud layer group)
    elseif (nint(r_mnh(1)) == 9) then
       clcd (1)  = -1._sp
       clcd (2)  = -1._sp
    !   assign cloud layer from cloud base height (CBH)
    elseif (inh(1) /= -1) then
       !   inconsistency: ('MNH' = 0 .and. CBH < 2500 m)
       !                 (note: CBH = 2500 m in [A-1600] for zero cloud is ok)
       if ((inh(1) < 2500) .and. (nint(r_mnh(1)) == 0)) then
          clcd (1)  = -1._sp
          clcd (2)  = -1._sp
          inh  (1)  = -1          ! discard CBH if inconsistent with MNH(1)?
       else
          if (inh(1) <= H_LOW)  clcd (2) = -1._sp  ! --> low cloud
          if (inh(1) >  H_LOW)  clcd (1) =  0._sp  ! --> no low cloud
          if (inh(1) >  H_MID)  clcd (2) =  0._sp  ! --> no mid-level cloud, but
          if (inh(1) >  H_MID)  clcd (3) = clcd(4) !     ... high cloud
       endif
    !   vertical significance ('MVTSU'): low, middle, or high cloud
    !   -> very few reports exist with missing CBH and existing 'MHN' +'MVTSU'
    !      from which cloud layer could be derived (if trusting 'MVTSU')
    elseif ((mvtsu(1) >= 7) .and. (mvtsu(1) <= 9)) then
       if (mvtsu(1) == 7)  clcd(2) = -1._sp       ! low cloud (WMO Table 008002)
       if (mvtsu(1) >= 8)  clcd(1) =  0._sp       ! mid-level or high cloud
       if (mvtsu(1) == 9)  clcd(2) =  0._sp       ! high cloud
       if (mvtsu(1) == 9)  clcd(3) = clcd(4)      ! high cloud
       !  (mvtsu = 62 = 'value not applicable': in many countries (DE, CH, ...)
       !   this denotes 'clear sky (no clouds)'; however e.g in DK (mvtsu == 62)
       !   .and.('N' > 0) occurs --> do not infer (zero) low/mid/high cloud
    !   cloud type check (B/C 1.4.4.5; 020012)     (not done in COSMO operators)
    elseif ((btest (chk_cldinfo, 1)) .and. (     (mcc(1) == 30)                &
                                            .or. (mcc(2) == 20))) then
       if (mcc(1) == 30)   clcd(1) =  0._sp       ! no C_L clouds
       if (mcc(2) == 20)   clcd(2) =  0._sp       ! no C_M clouds
    !   it cannot be determined whether 'MNH' is low or mid-level cloud
    elseif (min( clcd(1), clcd(2) ) > TOL) then
       clcd (1) = -1._sp
       clcd (2) = -1._sp
    elseif ((btest (chk_cldinfo, 1)) .and. (mcc(3) == 10)) then
       clcd (3) =  0._sp        ! no C_H clouds
    endif
    !   cloud and derived high cloud: obscured or undefined, or missing
    if (nint(clcd(4)) == 9)  clcd (4) = -1._sp
    if (nint(clcd(3)) == 9)  clcd (3) = -1._sp
    ! (hereafter -1. <= clcd(:) <= 8.)

    ! consistency checks:
    ! - discard cloud levels if vert. significance inconsistent with
    !   CBH (emphasis on low cloud, with tolerances) or with 'MNH'
    !   (note: tolerances for low / middle cloud bounds are asymmetric because
    !          middle cloud (type) with CBH < H_LOW can reach up into mid-level
    !          cloud layer, whereas low cloud (type) with CBH >= H_LOW cannot
    !          reach down into low cloud layer)
    !   (note: for high cloud (ncsig == 9), reported CBH can be 2500 m)
    lcc1 =      ((mvtsu(1) == 7) .and. (inh(1) > H_LOW+ 400))                  &
           .or. ((mvtsu(1) == 8) .and. (inh(1) < H_LOW- 600))                  &
           .or. ((mvtsu(1) == 9) .and. (inh(1) < H_MID-2000)                   &
                                 .and. (inh(1) /= 2500))
    lcc2 =      ((mvtsu(1) == 7) .and. (nint(r_mnh(1)) == 0))                  &
           .or. ((mvtsu(1) == 8) .and. (nint(r_mnh(1)) == 0))                  &
           .or. ((mvtsu(1) == 9) .and. (nint(r_mnh(1)) >= 1))
    if ((lcc1) .or. (lcc2)) then
       clcd (1) = -1._sp
       clcd (2) = -1._sp
       clcd (3) = -1._sp
       if (lcc1)  inh  (1) = -1    ! discard CBH if inconsistent with mvtsu(1)?
       if (lcc2)  r_mnh(1) = -1._sp
    endif
    ! - total cloud cover < low or mid-level cloud cover:
    !   - allow tolerance of 1 octa and adjust total cloud cover
    if (nint(clcd(4)) /= -1) then
       if (clcd(4) >= clcd(1)-1._sp-TOL)  clcd (4)  = max( clcd(4), clcd(1) )
       if (clcd(4) >= clcd(2)-1._sp-TOL)  clcd (4)  = max( clcd(4), clcd(2) )
    !   - otherwise reject cloud amounts
       if (clcd(4) < max( clcd(1), clcd(2) )-TOL) then
          clcd (1)  = -1._sp
          clcd (2)  = -1._sp
          clcd (4)  = -1._sp
       endif
    endif
    !   if high cloud cover, adjust (TAC) reported CBH of 2500 m to high CBH
    !                        to make it consistent
    if ((nint(clcd(3)) >= 1) .and. (inh(1) == 2500))  inh(1) = H_MID + 501
!   if (prt_data) then
!      write(*,*) "### cloudcover_synop: GCG: tcc=", tcc, ", fog index=", kfog
!      write(*,*) "###  nh  =", inh  (1:nlay)
!      write(*,*) "### mnh  =", r_mnh(1:nlay)
!      write(*,*) "### mvtsu=", mvtsu(1:nlay)
!      write(*,*) "### clcd l/m/h/t =", clcd (1:4)
!   endif

    ! Ceiling from GCG (note: only either low .xor. mid-level .xor. high cloud
    ! ----------------        may be non-zero according to the above derivation)
    !   (general question is whether GCG info is reliable enough to construct
    !    a ceiling observation, or not?)
    ceilg = -1
    if ((maxval( clcd(1:3) ) > 4.1_sp ) .and. (inh(1) /= -1)) then
      ceilg = inh(1)
    !    if total cloud cover is > 0 (but <= 4 octa), CBH is requested to be
    !    reported here to avoid the situation that stations never reporting
    !    CBH can only yield ceiling = cbh_clr (which leads to biased obs)
    !    (if total cloud cover == 0, CBH is not requested to be reported)
    elseif (       (clcd(4) >= -TOL) .and. (clcd(4) <= 4.1_sp)                 &
            .and. ((clcd(4) <=  TOL) .or. (inh(1) /= -1))) then
      ceilg = CBH_CLR
    endif

    !  discard CBH unless non-zero cloud cover is oberved or zero total cloud
    if ((maxval( clcd(1:3) ) <= TOL) .and. (clcd(4) < -TOL))  inh(1) = -1

    ! --------------------------------------------------------------------------
    ! derive / update low, mid-level, high cloud amount, cloud base height (CBH)
    ! and ceiling from individual cloud layer group (ICLG)
    ! --------------------------------------------------------------------------
    ! - CBH, ceiling: use ICLG to replace info from GCG (general cloud group)
    !   ~~~~~~~~~~~~              =======          ... where possible, because
    !   - CBH in GCG does no need to be exact but often indicates some value
    !     within the correct layer according to [A-1600]
    !   - e.g. in case of several low cloud layers and cloud cover <= 4 octas
    !     for the lowest layer and >= 5 octas for another low layer, the use of
    !     GCG info will result in either (both usages of 'MHN' are found):
    !     - too little low cloud cover if 'MHN' refers to the lowest layer only,
    !     - ceiling being erroneously assigned the base of the lowest layer if
    !       'MNH' refers to all low cloud layers
    ! - cloud amounts: use ICLG to cross-check and complement info from GCG
    !   ~~~~~~~~~~~~~              ===========     ==========

    !   remove 'cloud' layers with layer index > 1 (i.e. i>2) and zero CBH
    !     (reported e.g. by some Belgian stations)
    k = nlay
    do i = k, 3, -1
       if (inh(i) == 0)  nlay = max( nlay , i-1 )
    enddo
    !   cloud amount = 9: 'sky obscured by fog or other meteorol. phenomena'
    !   --> overcast cloud only if lowest level .and.
    !       ((reported fog).or.((ceiling or low cloud base height) but no dust))
    if (nlay >= 2) then
      lcil = (mvtsu(2) == 5) .or. ((inh(2) < 100) .and. (inh(2) /= -1))
      if ((nint(r_mnh(2)) == 9) .and.(     ((kfog /= -1) .and. (lcil))         &
                                      .or. ( kfog >=  3)))  r_mnh(2) = N_OVC
    endif
    !   --> otherwise doubtful and discarded
    do i = 2, nlay
       if (nint(r_mnh(i)) == 9)  r_mnh(i) = -1._sp
    enddo
    ! (hereafter -1. <= r_mnh(:) <= 8.)
    !   (note: cross-checking with vert. significance is not done for individual
    !          cloud layers as values 7 - 9 (low, middle, high) are not reported
    !          (only values 1 - 4, 21 - 24, 5 (and 0, 63) are found)

    ! derive low / mid-level / high cloud amounts, cloud base height, ceiling
    ! based on individual cloud layer group info (i.e. i >= 2) only
    ! -----------------------------------------------------------------------
    clcx(:)  =  -1._sp ! cloud amount (low, mid-level, high, total)
    clcml    =  15._sp ! lowest mid-level cloud layer
    cbhx     =  -1     ! CBH
    ceilx    =  -1     ! ceiling
    do i = nlay , 2 , -1
       !   use only layers with both valid + non-zero cloud amount and valid CBH
       if ((nint(r_mnh(i)) >= 1) .and. (inh(i) /= -1)) then
          !   derive values in top down loop
          !   - cloud base height, ceiling
                                  cbhx   = inh(i)
          if (r_mnh(i) > 4.1_sp)  ceilx  = inh(i)
          !   - cloud amount of low / mid-level / high cloud layers
          if     (inh(i) <= H_LOW) then
             clcx (1)  = max( clcx(1), r_mnh(i) )
          elseif (inh(i) <= H_MID) then
             clcx (2)  = max( clcx(2), r_mnh(i) )
             !     lowest mid-level layer
             clcml     =               r_mnh(i)
          else
             clcx (3)  = max( clcx(3), r_mnh(i) )
          endif
       endif
    enddo
    !   complement cloud-free layers, if cloud layers exist only further above
    !   (missing mid-level cloud remains missing if low + high cloud report > 0)
    if ((clcx(3) > TOL) .and. (clcx(1) <= TOL)) clcx(2)  = max( clcx(2), 0._sp )
    if (max( clcx(3), clcx(2) ) > TOL)          clcx(1)  = max( clcx(1), 0._sp )

    ! consistency check on ceiling from GCG:
    ! -------------------------------------
    !   ceiling based on 'MNH'-based cloud cover >= 5 octa assumed suspect
    !   if ICLG-derived cloud layer further above is < 5 octa
    if ((ceilg < CBH_CLR) .and. (ceilg /= -1)) then
      if (     ((ceilg <= H_MID) .and. (clcx(3) < 4.1_sp)                       &
                                 .and. (clcx(3) >= -TOL))                       &
          .or. ((ceilg <= H_LOW) .and. (clcx(2) < 4.1_sp)                       &
                                 .and. (clcx(2) >= -TOL)))  ceilg = -1
    endif

    !   for some British stations, mid-level or high cloud < low cloud cover
    !   --> false reporting of mid-level resp. high cloud cover
    if ((min( clcx(1), clcx(2) ) >= -TOL) .and. (clcx(2) < clcx(1)-TOL))       &
      clcx (2) = -1._sp
    if ((min( clcx(1), clcx(3) ) >= -TOL) .and. (clcx(3) < clcx(1)-TOL))       &
      clcx (3) = -1._sp
    if ((min( clcx(2), clcx(3) ) >= -TOL) .and. (clcx(3) < clcx(2)-TOL))       &
      clcx (3) = -1._sp

    ! consistency checks on cloud cover from GCG:
    ! ------------------------------------------
    !   discard ('MNH'-derived) cloud amounts if suspect (overestimated)
    !   (i.e. if contradicting ICLG which is given preference;
    !    after these checks, amounts from ICLG or from GCG may still be
    !    imcomplete, i.e. underestimated, since not all cloud layers are
    !    always reported resp. accounted)
    !   high cloud:
    !    - 'MNH' sometimes reports zero by wrongly assigning (mid-level) cloud
    !      (mainly above 6000m) -> high cloud derived from total cloud suspect
    !      - unless its cloud amount is clcd(3) > clcx(2): in this case assume
    !        that 'clcd(3)' was based on a different (higher) cloud layer); or
    !      - if mid-level (or low) cloud with base below ~ 5000 m exists
    !    - high cloud derived from total cloud is suspect if low cloud exists
    if (      (clcd(3) > TOL)                                                  &
        .and. (     (clcd(3) <= clcx(2)+TOL)                                   &
               .or. ((cbhx < 5000) .and. (cbhx /= -1))                         &
               .or. (clcx(1) > TOL)))                                          &
       clcd (3)  = -1._sp
    !   mid-level cloud:
    !    - high cloud being misplaced as mid-level cloud in 'MNH' is assumed
    !      if individual cloud layer group contains only high cloud;
    !    - low cloud with base height <= 1999 in individual cloud layer group
    !      is occasionally misplaced as mid-level cloud in 'MNH' (this must be
    !      the lowest cloud layer - if an even lower cloud layer existed, it
    !      would have been interpreted as low cloud in 'MNH')
    if (      (clcd(2) > TOL)                                                  &
        .and. (     ((max( clcx(1), clcx(2) ) <= TOL) .and. (clcx(3) >  TOL))  &
               .or. (clcx(1) > TOL)))                                          &
       clcd (2)  = -1._sp
    !   low cloud: suspect
    !    - if individual cloud group reports cloud layers, but no low cloud; or
    !    - if 'MNH' (>)== clcml (= cloud amount of lowest mid-level layer)
    !      ('MNH'-derived low cloud occasionally includes 1 cloud layer with
    !       2000 m <= CBH < 2500 m (above another low cloud layer))
    if (      (clcd(1) > TOL)                                                  &
        .and. (     ((max( clcx(3), clcx(2) ) >  TOL) .and. (clcx(1) <= TOL))  &
               .or. (clcd(1) >= clcml-TOL)))                                   &
       clcd (1)  = -1._sp

!   if (prt_data) then
!      write(*,*) "### clcx l/m/h/t =", clcx (1:4)
!      write(*,*) "### clcd l/m/h/t =", clcd (1:4)
!      write(*,*) "### ceiling, cbh=", ceilx, cbhx
!   endif

    ! -----------------------------------
    ! assign values to SYNOP derived type
    ! -----------------------------------
    ! cloud amounts: take max of general and individual cloud group info
    ! ~~~~~~~~~~~~~  (info from any cloud group may still be incomplete and
    !                 not include all layersclcd. of low cloud)
    !                (zero amount is also allowed; missing is < 0)
    clcd (1)  = max( clcd(1), clcx(1) )
    clcd (2)  = max( clcd(2), clcx(2) )
    clcd (3)  = max( clcd(3), clcx(3) )
    if (clcd(1) >= -TOL)  s% clcl = real( clcd(1) , wp )
    if (clcd(2) >= -TOL)  s% clcm = real( clcd(2) , wp )
    if (clcd(3) >= -TOL)  s% clch = real( clcd(3) , wp )
    if (clcd(4) >= -TOL)  s% n    = max( s% n, nint (clcd(4) * 12.5_sp) )
    !                     ------- subr output

    ! cloud base height CBH, ceiling:  prefer individual layer cloud group info
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    icbx = cbhx
    icix = ceilx
    if ((cbhx  == -1) .and. (inh(1) /= -1))  cbhx  = inh(1)
    if ((ceilx == -1) .and. (ceilg  /= -1))  ceilx = ceilg
    if (cbhx  /= -1)  s% nh   = real( cbhx  , wp )
    if (ceilx /= -1)  s% ceil = real( ceilx , wp )
    !                 ------- subr output

    !   control output
    if (prt_data) then
      hhmm = nint( hhmmss / 60._wp )
      hhmm = (hhmm / 60) *100 + mod( hhmm, 60 )
      icbh = -1
      icil = -1
      ivis = -1
      if (s% ceil /= rvind)  icil = nint( s% ceil )
      if (s% nh   /= rvind)  icbh = nint( s% nh   )
      if (s% vis  /= rvind)  ivis = nint( s% vis  )
      ivsgg = mvtsu(1)
      if (ivsgg == -1)  ivsgg = 63
      do i = 2 , nlay
        ivsg = min( mvtsu(i), 63 )
        if (ivsg == -1)  ivsg = 63
        write( *,'("### CLX ",a8,i5  ,2i3,3f4.0,i6,i3,2i6)' ) statid, hhmm     &
             , i-1, ivsg, clcx(1), clcx(2), clcx(3), inh(i), mnh(i), icbx, icix
      enddo
      write( *,'("### Cloud ",a8,i2,":",i2.2,i5,i4,i3,i6,i3," saved:",4i3,2i6  &
           &    ,i7)')                     statid, hhmm/100, mod(hhmm,100)     &
           , ww, mod(  mnh(1)   +100,100), ivsgg, inh(1), kfog                 &
           , abs(nint( s% clcl )), abs(nint( s% clcm )), abs(nint( s% clch ))  &
           , mod(nint( clcd(4) )+100,100), icbh, icil, ivis
    endif

    !   store (nearly) originally reported cloud info in packed form
    !   ------------------------------------------------------------
    !   (these bit patterns differ from the COSMO operators)
    kcc(:) = MIS_CT
    !   type of low cloud       -> [A-0513]
    if ((mcc(1) >= 30) .and. (mcc(1) <= 39)) kcc(1) = mcc(1) - 30
    if ((mcc(1) == 59) .or.  (mcc(1) == 62)) kcc(1) = 10
    !   type of mid-level cloud -> [A-0515]
    if ((mcc(2) >= 20) .and. (mcc(2) <= 29)) kcc(2) = mcc(2) - 20
    if ((mcc(2) == 59) .or.  (mcc(2) == 62)) kcc(2) = 10
    !   type of high cloud      -> [A-0509]
    if ((mcc(3) >= 10) .and. (mcc(3) <= 19)) kcc(3) = mcc(3) - 10
    if ((mcc(3) == 59) .or.  (mcc(3) == 62)) kcc(3) = 10
    !   height key: [A-1600], lower limits of bins: 0,   50,  100,  200,  300,
    !                                             600, 1000, 1500, 2000, 2500 m)
    if (knh(1) == -1) then
      nhkey =  MIS_HK
    elseif (knh(1) <=   99) then
      nhkey =  knh(1)        /  50
    elseif (knh(1) <=  399) then
      nhkey =  knh(1)        / 100 + 1
    elseif (knh(1) <=  999) then
      nhkey = (knh(1) - 200) / 400 + 1
    elseif (knh(1) <= 2499) then
      nhkey =  knh(1)        / 500 + 1
    else
      nhkey =  9
    endif
    if (nclct == -1)  nclct = MIS_N_
    do i = 1 , nlay
      kvs (i) = mvtsu(i)
      if ((kvs(i) < 0) .or. (kvs(i) > MIS_SG))  kvs(i) = MIS_SG
      if (mnh(i) == -1)  mnh(i) = MIS_N_
      if (knh(i) == -1)  knh(i) = MIS_H * 10
      if (i >= 2)        knh(i) = knh(i) / 10
    enddo

    igcg  =  MIS_GCG
    igcg  =  ireplace( igcg, NBP_CTH, MIS_CT, kcc(3) )
    igcg  =  ireplace( igcg, NBP_CTM, MIS_CT, kcc(2) )
    igcg  =  ireplace( igcg, NBP_HK , MIS_HK, nhkey  )
    igcg  =  ireplace( igcg, NBP_CTL, MIS_CT, kcc(1) )
    igcg  =  ireplace( igcg, NBP_NH , MIS_N_, mnh(1) )
    igcg  =  ireplace( igcg, NBP_N  , MIS_N_, nclct  )
    igcg  =  ireplace( igcg, NBP_SG , MIS_SG, kvs(1) )
    if (igcg /= MIS_GCG)  s% gwg  =  igcg
    !                     ------ subr output
    iclg (:)  =  MIS_ICG
    do i = 1 , min( nlay-1 , 4 )
      iclg (i)  =  ireplace( iclg(i), IBP_SG, MIS_SG, kvs(i+1) )
      iclg (i)  =  ireplace( iclg(i), IBP_N , MIS_N_, mnh(i+1) )
      iclg (i)  =  ireplace( iclg(i), IBP_H , MIS_H , knh(i+1) )
      if (iclg(i) /= MIS_ICG) then
        if (i == 1)  s% icl1  =  iclg(i)
        if (i == 2)  s% icl2  =  iclg(i)
        if (i == 3)  s% icl3  =  iclg(i)
        if (i == 4)  s% icl4  =  iclg(i)
        !            ------- subr output
      endif
    enddo

  end subroutine cloudcover_synop_lam
!------------------------------------------------------------------------------
  subroutine cloudcover_metar (s, mgwi, nlay, mvtsu, mnh, mcc, nh)
    type(t_synop), intent(inout) :: s           ! SYNOP derived type
    integer,       intent(in)    :: mgwi        ! Gen.weather ind. (020009)
    integer,       intent(in)    :: nlay        ! Number of cloud layers
    integer,       intent(in)    :: mvtsu(:)    ! Vertical signif. (008002)
    integer,       intent(in)    :: mnh  (:)    ! Cloud amount     (020011)
    integer,       intent(in)    :: mcc  (:)    ! Cloud type       (020012)
    real,          intent(in)    :: nh   (:)    ! Height of cloud base [m]
    !-------------------------------------------------------------------------
    ! Interpret cloud data in METAR reports (c.f. ICAO annex 3, section 4.5.4)
    ! (Currently only low cloud cover is kept)
    !-------------------------------------------------------------------------
    integer             :: i, k
    real                :: h                    ! Cloud height above msl
    real                :: octa(4)              ! Low/mid/high/total cloud cover
    integer             :: clcx(4)              ! Low/mid/high/total cloud cover
    integer             :: clc (nlay)           ! Temporary cloud cover
    integer             :: ix  (nlay)           ! Auxiliary index
    real(sp), parameter :: H_LOW = 2000._sp     ! Limit between low and mid
    real(sp), parameter :: H_MID = 7000._sp     ! Limit between mid and high
    real(sp), parameter :: N_FEW = 1.5_sp       ! Average amount for FEW (1-2)
    real(sp), parameter :: N_SCT = 3.5_sp       ! Average amount for SCT (3-4)
    real(sp), parameter :: N_BKN = 6.0_sp       ! Average amount for BKN (5-7)
    real(sp), parameter :: N_OVC = 8.0_sp       ! Overcast

    !--------------------------------------------
    ! No significant cloud detected or sky clear?
    !--------------------------------------------
    if (mgwi >= 1 .and. mgwi <= 3) then
       s% clcl = 0._sp
       return
    end if

    !-------------------------------
    ! No cloud information in report
    !-------------------------------
    if (nlay == 0) then
       s% clcl = inva_syn% clcl
       return
    end if

    !-------------------------------------------------------
    ! Single layer: treat clear sky and fog/height undefined
    !-------------------------------------------------------
    if (nlay == 1 .and. nh(1) == rvind) then
       if      (mnh(1) == 0) then
          s% clcl = 0._sp
          return
       else if (mnh(1) == 9) then
          s% clcl = N_OVC
          return
       end if
    end if

    clcx = -1
    clc  = -1
    ix   = -1

    !-------------------
    ! Check cloud height
    !-------------------
    do i = 1, nlay
       select case (mvtsu(i))
       case (7)
          ix(i) = 1
       case (8)
          ix(i) = 2
       case (9)
          ix(i) = 3
       case (1:5)
          h = nh(i)       ! Height of base of cloud above local ground
!         if (h /= rvind .and. s% zr% o /= rvind) h = h + s% zr% o
          if      (h <  H_LOW) then
             ix(i) = 1
          else if (h <  H_MID) then
             ix(i) = 2
          else if (h /= rvind) then
             ix(i) = 3
          else
             ix(i) = -1
          end if
       end select
    end do

    !---------------------------------------------
    ! Check cloud amount for classes we can handle
    !---------------------------------------------
    do i = 1, nlay
       select case (mnh(i))
       case (0)
          if (nh(i) == rvind) then
             clc(i) = 0     ! analogous to B/C1.4.4.4.2 ?
          else
             ix (i) = -1    ! Dismiss 0 octas for defined cloud base height!
          end if
       case (1:9,11:13)
          clc(i) = mnh(i)
       case default
          ix (i) = -1
       end select
    end do

    if (prt_data) then
       write(*,*) "### cloudcover_metar: n=", s% n
       write(*,*) "###  nh =",  nh(1:nlay)
       write(*,*) "###  ix =", ix (1:nlay)
       write(*,*) "### mnh =", mnh(1:nlay)
       write(*,*) "### clc =", clc(1:nlay)
    end if

    !--------------------------------
    ! Merge information within layers
    ! keep highest amount (020011)
    !--------------------------------
    do i = 1, nlay
       k = ix(i)
       if (k > 0) then
          select case (clc(i))
          case (0:9)
             clcx(k) = max (clcx(k), clc(i))
             if (clcx(k) == 9) clcx(4) = 8      ! Sky obscured by fog etc.
          case (13)     ! Few       (1-2 octas)
             if (clcx(k) < 2) then
                clcx(k) = 13
             end if
          case (11)     ! Scattered (3-4 octas)
             if (clcx(k) < 4 .or. clcx(k) == 13) then
                clcx(k) = 11
             end if
          case (12)     ! Broken    (5-7 octas)
             if (clcx(k) < 6 .or. clcx(k) == 11 .or. clcx(k) == 13) then
                clcx(k) = 12
             end if
          end select
       end if
    end do

    !------------------------------------------
    ! Canonicalize for categorical verification
    ! convert to (representative) octas
    !------------------------------------------
    octa = -1._sp
    do i = 1, 3
       select case (clcx(i))
       case (0:8)
          octa(i) = real (clcx(i), sp)
       case (9)             ! Sky obscured by fog etc. (8 octas)
          octa(i) = N_OVC
       case (11)            ! Scattered (3-4 octas)
          octa(i) = N_SCT
       case (12)            ! Broken    (5-7 octas)
          octa(i) = N_BKN
       case (13)            ! Few       (1-2 octas)
          octa(i) = N_FEW
       end select
    end do

    if (prt_data) then
       write(*,*) "### clcx:", clcx(:)
       write(*,*) "### octa:", octa(:), "#1"
    end if

    !-------------------------------
    ! Recover low/medium "clear sky"
    ! when no low cloud reported
    !-------------------------------
    if   (octa(1) <= 0._sp) then
      if (octa(3) >= 0._sp .and. octa(2) < 0._sp) octa(2) = 0._sp
      if (octa(2) >= 0._sp .and. octa(1) < 0._sp) octa(1) = 0._sp
    end if

    !-------------------------------------
    ! Plausibility checks for cloud layers
    !-------------------------------------
!   if ((octa(2) >= 0._sp .and.  octa(2) <  octa(1)) .or. &
!       (octa(2) >  0._sp .and.  octa(2) == octa(1))      ) octa(2) = -1._sp
!   if ( octa(3) >= 0._sp .and. (octa(3) <  octa(1)  .or. &
!                                octa(3) <  octa(2)      )) octa(3) = -1._sp
!   if ( octa(3) >  0._sp .and. (octa(3) == octa(1)  .or. &
!                                octa(3) == octa(2)      )) octa(3) = -1._sp

    if (prt_data) then
       write(*,*) "### octa:", octa(:), "#2"
    end if

    if (octa(1) >= 0._sp) s% clcl = octa(1)
    if (octa(2) >= 0._sp) s% clcm = octa(2)

  end subroutine cloudcover_metar
!------------------------------------------------------------------------------
  elemental integer function ireplace ( invar, ipos, ibocset, irepl )
    integer,       intent(in)    :: invar, ipos, ibocset, irepl
  !------------------------------------------------------------------------
  ! replaces 'iboc' bits starting at bit position 'ipos' of integer 'invar'
  !   by the 'iboc' bits from integer word 'irepl',
  !   assuming that:  ibocset = 2**(iboc) - 1
  !------------------------------------------------------------------------
    ireplace = ior( iand( invar, not( ishft( ibocset, ipos ) ) )              &
                  , ishft( iand( irepl, ibocset ), ipos ) )

  end function ireplace
!------------------------------------------------------------------------------
end module mo_synop
