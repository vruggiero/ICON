!
!+ Interface 3D-Var/LETKF <-> COSMO conventional observation operators
!
MODULE mo_cosmo_conv
!
!-------------------------------------------------------------------------------
! Description:
!   Interface 3D-Var/LETKF <-> COSMO conventional observation operators
!   This module contains the driver routines for
!     - the reading of conventional obs by the COSMO routines
!     - the COSMO observation operators and quality control (first guess and
!       multi-level checks) for conventional non-gridded data
!
!   The module contains the following public procedure:
!     - process_cosmo_conv   : driver routine for all tasks
!   and the following private procedures:
!     - run_operator         : driver routine for running the COSMO obs operator
!                              and quality control for one report
!     - cosmo_to_obs         : converts COSMO obs data structure 'ODR' for
!                              into 3dvar/LETKF obs data structure 't_obs'
!     - check_store_conv     : converts observation body of 'ODR' into 't_obs'
!
! Written by        : DWD, Christoph Schraff and Andreas Rhodin
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_22        2013-02-13 Christoph Schraff
!  Interface to COSMO conventional observation operators
! V1_28        2014/02/26 Andreas Rhodin
!  new interface to new_int
! V1_35        2014-11-07 Andreas Rhodin
!  adaptions for MEC
! V1_42        2015-06-08 Andreas Rhodin
!  implement temporal interpolation for MEC
! V1_43        2015-08-19 Andreas Rhodin
!  changes for COSMO MEC
! V1_44        2015-09-30 Andreas Rhodin
!  update shared modules to COSMO 5.03-beta
! V1_45        2015-12-15 Andreas Rhodin
!  fix bug pointed out by Josue Gehring (meteoswiss) in level determination
! V1_47        2016-06-06 Andreas Rhodin
!  changes for COSMO-MEC
! V1_48        2016-10-06 Andreas Rhodin
!  move subroutine obs_trange to mo_obs_trange, rename to apply_trang
!  namelist /observations/: select pressure bounds for lapse rate estimation
! V1_50        2017-01-09 Harald Anlauf
!  MEC: bugfix for total cloud cover (octas!)
! V1_51        2017-02-24 Andreas Rhodin
!  compile with COSMO observation operators
!              2018-05-15 Christoph Schraff
!  processor element index 'p_pe' replaced by 'dace% pe'.
! 2022-08-01 Christoph Schraff
!   Cloud ceiling introduced as (derived) observation.
!   Visibility collected from model to prepare model equivalent.
!   Std. deviation of SSO filled in 'spt' for station selection of 10-m wind.
! 2023-05-24 Christoph Schraff
!   Processing of tower obs added, including call of 'shift_profile' to allow
!   (alternatively) for use of obs at correct height of sensor above ground.
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!===============================================================================

!=============
! Modules used
!=============
  !-------------------------
  ! general purpose routines
  !-------------------------
  use mo_exception,   only: finish        ! abort routine
  use mo_kind,        only: wp,          &! working precision
                            sp,          &! single precision (used in 't_datum')
                            i8            ! 8-byte integer
  use mo_mpi_dace,    only: dace          ! MPI group info
  use mo_time,        only: t_time,      &! Data type to hold time
                            zero_time,   &! zero time = fill value
                            init_time,   &! set time variable
                            operator (+)  ! operator adding times
  !-----------------------------
  ! access atmospheric data type
  !-----------------------------
  use mo_atm_state, only: t_atm            ! atm. state data type
  use mo_atm_grid,  only: t_grid           ! atm. grid  data type
  use mo_grid_intpol,only: alloc_imcol,   &! (re)allocate t_imcol
                           Grid_Indices,  &! get indices of nearest gridpoint(s)
                           add_index       ! set model column indices
  use mo_physics,   only: gacc,           &! gravity acceleration
                          rdv => RdRd,    &! R/Rd
                          r_d => R,       &! gas constant of dry air [J/(kg*K)]
                          r_v => Rd,      &! gas constant of water vapour
                          b1,             &! \  constants for computing
                          b2w,            &!  \ the saturation vapour pressure
                          b3,             &!  / over water
                          b4w              ! /  by the Magnus formula
  use mo_t_col,     only: t_cols,         &! model columns data type
                          COL_UV,         &!              hor.wind comp.
                          COL_T, COL_Q,   &!              temp., spec.hum.
                          COL_GEO,        &!              geop.h. full levs.
                          COL_P,          &!              press.  full levs.
                          COL_CLC,        &!              cloud cover
                          COL_RANGE,      &! variables with time range
                          nt,             &! number of time ranges,
                          tr               ! time ranges in minutes
  !-----------------------------
  ! access observation data type
  !-----------------------------
  use mo_cosmo_obs, only:  init_cosmo_obs, &! initialise mo_cosmo
                           read_cosmo_obs   ! read COSMO data
  use mo_obs_set,   only:  t_obs_block      ! obs data type
  use mo_t_obs,     only:  t_obs,          &!
                           t_spot,         &!
                           t_head,         &! observation data type
                           ptop_lapse,     &!
                           pbot_lapse,     &!
                           read_cdfin,     &! flag to read COSMO observations
                           new_spot,       &! reserve memory
                           new_obs,        &! reserve memory
                           new_int,        &! reserve memory for interp.space
                           set_xuv,        &! set unit vectors, zenith angle
                           TSK_INIT,       &! FLAGS: initialize module
                           TSK_READ,       &!  read observations
                           TSK_SET_CHR,    &!  set observation characteristics
                           TSK_SETUP_COLS, &!  determine model columns required
                           TSK_SETUP_FULL, &!  setup description of PSAS-space
                           TSK_SETUP_FUL0, &!  setup size of PSAS-space
                           TSK_SHRINK,     &!  release unused obs. in report
                           TSK_Y,          &!          evaluate nonlinear oper.
                           TSK_R,          &!   setup observational errors
                           CHR_NONL,       &! H is nonlinear
                           CHR_EXP,        &! H is 'expensive'
!                          OBS_TV,         &! interp. observation type: Tv
!                          OBS_RH,         &! interp. observation type: rh
!                          OBS_HS,         &! interp. observation type: geop.
!                          OBS_U,          &! interp. observation type: u
!                          OBS_V,          &! interp. observation type: v
                           COSMO,          &! module type
                           ITY_MCOLS        ! interpolation type: column
  use mo_t_use,      only: STAT_ACTIVE,    &! active
                           STAT_REJECTED,  &! rejected
                           STAT_DISMISS,   &! dismissed
                           STAT_OBS_ONLY,  &!
!                          STAT_FORGET,    &!
                           CHK_NONE,       &! valid data        check id
!                          CHK_DATASET,    &!
                           CHK_FG,         &! first guess       check id
                           CHK_CONSIST,    &! consistency check flag
                           CHK_DOMAIN,     &! domain check flag
                           CHK_NO_OBS,     &! flag for no real obs
                           t_use,          &! data type to hold state
                           decr_use,       &! decrease state of datum
                           use_0,          &! default values of type use
                           change_code,    &! convert code (3dvar to feedback)
                           change_bits,    &! convert bit-field
                           reverse_code,   &! convert code back (fdbk to 3dvar)
                           reverse_bits,   &! convert bit-field back
                           checks => chk,  &! table of 'checks'
                           stats            ! table of 'states'
  use mo_t_datum,    only: t_datum,        &! date+time derived type
                           rvind            ! invalid value indicator
  use mo_t_table,    only: bit1_pos_nr,    &! search first bit=1 in given order
                           name_value       ! derive name from (table, value)
  use mo_obs_tables, only: decr_rpt_use,   &! change use-flags of report
                           check_report_1, &! basic checks on reports
                           check_report_0   ! basic checks on reports
  use mo_obs_trange, only: apply_trange     ! apply time range observation operator
  use mo_fdbk_tables,only: flags          ,&! table for check flags
                           init_fdbk_tables,&! initialise the tables
                           obstype        ,&! observationtype table
                           ST_OBS_ONLY    ,&! no model equivalent available
                           ST_ACTIVE      ,&! used in the assimilation
                           ST_PASSIVE     ,&! not used, only monitored
                           ST_REJECTED    ,&! not used, suspicious quality
!                          ST_PAS_REJ     ,&! passive and rejected
                           FL_BLACKLIST   ,&! blacklist (or not whitelist)
                           FL_HEIGHT      ,&! '' not in valid height range
                           FL_PRACTICE    ,&! bad reporting practice/insf.data
                           FL_DATASET     ,&! dataset quality flags
                           FL_MERGE       ,&! merged reports (e.g. TEMP ABCD)
                           FL_GROSS       ,&! gross error flag
                           FL_THIN        ,&! thinning
                           FL_FG          ,&! obs-fg check
                           FL_FG_LBC      ,&! obs - lateral BC field check
                           FL_NONE        ,&! no flag set
                           OT_AIREP       ,&! aircraft      obs type identifier
                           OT_SATOB       ,&! SATOBS        obs type identifier
                           OT_SYNOP       ,&! SYNOP + SHIP  obs type identifier
                           OT_DRIBU       ,&! drifting buoy obs type identifier
                           OT_TEMP        ,&! TEMP          obs type identifier
                           OT_PILOT       ,&! PILOT         obs type identifier
                           OT_SCATT       ,&! scatterometer obs type identifier
                           OT_GPSGB       ,&! GPS ground-b. obs type identifier
                           OC_WP_EU       ,&! European  wind profiler  code type
                           OC_RA_EU       ,&! European sodar/rass rep. code type
                           OC_PR_US       ,&! wind/profiler/rass (USA) code type
                           OC_RAVAD       ,&! radar VAD wind profile   code type
                           OC_TOWER       ,&! tower profile            code type
                           OC_ICOS        ,&! ICOS tower profile       code type
                           OC_SYTEMP      ,&! TEMP-derived surface rep code type
                           LS_SURFACE     ,&! surface level
                           LS_STANDARD      ! standard level
!                          LS_SIGN          ! significant level
  use mo_fdbk_tables,only: VN_P           ,&! pressure variable indicator
                           VN_PS          ,&! surface pressure  indicator
                           VN_U , VN_U10M ,&! zonal  wind (10m) var. indicator
                           VN_V , VN_V10M ,&! merid. wind (10m) var. indicator
                           VN_T , VN_T2M  ,&! temperature (2m)  var. indicator
                           VN_RH, VN_RH2M ,&! rel. hum (2m) variable indicator
                           VN_HEIGHT      ,&! height
                           VN_PTEND       ,&! pressure tendency
                           VN_W           ,&! vertical wind speed
                           VN_VGUST       ,&! vertical gust (aircraft)
                           VN_GUST        ,&! wind gust
                           VN_NH          ,&! cloud base height
                           VN_CEIL        ,&! cloud ceiling
                           VN_N           ,&! total cloud amount        (020011)
                           VN_N_L         ,&! low cloud amount
                           VN_N_M         ,&! middle cloud amount       (020011)
                           VN_N_H         ,&! high   cloud amount       (020011)
                           VN_VV          ,&! visibility
                           VN_WW          ,&! present weather
                           VN_GCLG        ,&! general cloud group          Table
                           VN_ICLG        ,&! individual cloud layer group Table
                           VN_PWC         ,&! precipitable water content kg/m**2
                           VN_ZPD         ,&! zenith path delay
                           VN_ZWD         ,&! zenith wet delay
!                          VN_SPD         ,&! slant path delay
                           VN_TMAX        ,&! maximum temperature        (K)
                           VN_TMIN        ,&! minimum Temperature        (K)
                           VN_TURB        ,&! degree of turbulence      (011031)
                           VN_SDEPTH      ,&! snow depth
                           VN_RR          ,&! precipitation amount      (kg/m^2)
                           VN_TRTR        ,&! time period of information (h)
                           VN_FLEV        ,&! nominal flight level       (m)
                           VN_NUM         ,&! ordinal (channel) number   (  )
                           VN_RAD_GL      ,&!
                           VN_RAD_DF      ,&!
                           VN_RAD_LW      ,&!
                           VN_TD2M        ,&!
                           VN_DD          ,&!
                           VN_FF
  !------------------------
  ! access matrix data type
  !------------------------
  use mo_dec_matrix,only: t_vector_segm    ! vector segment data type

!#undef NO_COSMO_OBS
#ifndef NO_COSMO_OBS
  !------------------------------------------------------
  ! access observation data type of COSMO data structures
  !------------------------------------------------------
  use data_obs_lib_cosmo , only :   &
    rmdi       ,& ! =-1.E31_ireals : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_ireals : check value in ODR for missing data
    epsy       ,& ! = 1.E-8_ireals : commonly used very small value > 0
    nupr       ,& ! unit number for printout (nupr<0: no printout)
    lwonl         ! .true. for the node (sub-domain) at which file with
                  !        the unit number 'nupr' is open
                  ! (i.e. where grid pt. (ionl ,jonl ) lies)

! use data_obs_cdfin     , only :   &
!   rprlim        ! = 1.1 : vertical distance limit for replacing 'missing data
                  ! [hPa]   in multi-level reports or removing colocated levels

  use data_obs_qc_limits , only :   &
    nqclev     ,& ! number of levels in the quality control threshold tables
    tabqclp    ,& ! ln(tabqcp(15))
    tabqcp     ,& ! levels of the error / threshold tables
    qczcorl    ,& ! radiosonde height error correlation matrix
    nzcorl        ! radiosonde height error correlation matrix [*1000]

  use data_obs_record    , only :   &
                  ! --> ODR header format parameters -->
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
!   mxghed     ,& ! header length of GPS reports
!   mxthed     ,& ! header length of satellite retrieval reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! time of observat. in forecast hours
!   nhsurf     ,& ! height of model grid pt. to which obs. is assigned
    nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
    nhtddb     ,& ! data base decoding time in forecast hours
    nhsolz     ,& ! solar zenith angle [deg]
    nhssos     ,& ! standard devation of sub-grid scale model orogrophy [m]
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
!   mxghdf     ,& ! header length of GPS reports
!   mxthdf     ,& ! header length of satellite retrieval reports
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhflag     ,& ! report flags (obs type, surf., alt., sta ID) (see 1.2.1)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
!   nhqcfw     ,& ! QC flags for surface pressure increments from TEMPs
!   nhcorr     ,& ! update sequence number (station correction indicator)
    nhcat      ,& ! data     category (from BUFR Section 1)
    nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre  +  (1000* sub-centre)
    nhstid     ,& ! station identity number
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhstyp     ,& ! sing-lv obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
    nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
!   ilstid     ,& ! character length of the station identity
!   ilstidp    ,& ! char. length used for printing the station ID
!   nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvscbp     ,& ! bit pos. for station correction indicator         "
!   nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
!   nvaaoc     ,& ! no. of bits occ. by aircraft roll angle           "
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
    nvapoc        ! no. of bits occ. by (extended) phase of flight    "
  use data_obs_record    , only :   &
                  ! --> ODR body format parameters for multi-level reports -->
    mxrbdy     ,& ! body length of multi-level reports
    nbtu       ,& ! u wind component [m/s]
    nbtv       ,& ! v wind component [m/s]
    nbtt       ,& ! temperature [K]
    nbtrh      ,& ! relative humidity [/]
    nbtp       ,& ! pressure [Pa]
    nbtz       ,& ! height [m]
    nbtuer     ,& ! error of observed wind component
    nbtter     ,& ! error of observed temperature
    nbtqer     ,& ! error of observed rel. humidity
    nbtzer     ,& ! error of observed height
    nbtlop     ,& ! LOG( pressure )
    nbtdrh     ,& ! bias correction for relative humidity [/] or RASS Tv
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! signal to noise ratio
    nbtuac     ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]
    mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')
  use data_obs_record    , only :   &
                  ! --> ODR body format parameters for single-level reports -->
    mxsbdy     ,& ! body length of single-level reports
    nbsu       ,& ! u wind component                                   [m/s]
    nbsv       ,& ! v wind component                                   [m/s]
    nbst       ,& ! temperature                                        [K]
    nbsrh      ,& ! relative humidity                                  [/]
    nbsp       ,& ! pressure                                           [Pa]
    nbsz       ,& ! height                                             [m]
    nbsuer     ,& ! error of observed wind component
    nbster     ,& ! error of observed temperature
    nbsqer     ,& ! error of observed relative humidity
    nbszer     ,& ! error of observed height
    nbspst     ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    nbscbs     ,& ! (lowest) cloud base height CBH above surface (AGL) [m]
    nbscil     ,& ! ceiling (CBH of lowest cloud layer >= 5 octas)     [m]
    nbscl      ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm      ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch      ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct      ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis     ,& ! (horizontal) visibility                            [m]
    nbsff      ,& ! wind speed                                         [m/s]
    nbsdd      ,& ! wind direction                                     [deg]
    nbstd      ,& ! dewpoint temperature                               [K]
    nbsrr1     ,& ! precipitation amount over 1  hour                  [mm]
    nbsrr3     ,& ! precipitation amount over 3  hours                 [mm]
    nbsrr6     ,& ! precipitation amount over 6  hours                 [mm]
    nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
!   nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3     ,& ! max. wind speed of gusts over 3 hour               [m/s]
    nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn      ,& ! minimum temperature (at 2m during past 12 hrs)     [K]
    nbstx      ,& ! maximum temperature (at 2m during past 12 hrs)     [K]
    nbsrad     ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
    nbsrdd     ,& ! diffuse   solar    radiation, sum over 1 hour      [J/m2]
    nbsrdl     ,& ! long-wave downward radiation, sum over 1 hour      [J/m2]
!   nbshsw     ,& ! total snow depth                                   [m]
    nbsdrh     ,& ! bias correction for relative humidity              [/]
    nbsvip     ,& ! model surface pressure (for LETKF localisation)    [Pa]
    mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
!   nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
    nbswwe     ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
    nbstur     ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
    nbsclg     ,& ! general           cloud       group (code)
    nbscl1     ,& ! first  individual cloud layer group (code)
    nbscl2     ,& ! second individual cloud layer group (code)
    nbscl3     ,& ! third  individual cloud layer group (code)
    nbscl4     ,& ! forth  individual cloud layer group (code)
    nbsttr     ,& ! time periods    (bit pattern of values, see below: nbsttr)
    nbsqal        ! quality info: per cent confidence (for AMV)        [%]
  use data_obs_record    , only :   &
                  ! --> ODR body format parameters for ground-based GPS rep. -->
    mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay [mm]
!   nbgiwv     ,& ! integrated water vapour [mm]
    nbgp       ,& ! pressure [Pa]
    nbgt       ,& ! temperature [K]
    nbgrh      ,& ! relative humidity [/]
    nbgbia     ,& ! bias correction to integrated water vapour [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour [mm]
    nbgz       ,& ! height [m]
    nbgzer     ,& ! error of observed height
    nbgqer     ,& ! error of observed relative humidity
!   nbgdrh     ,& ! bias correction for observed relative humidity
    nbgter     ,& ! error of observed temperature
    mxgbdf     ,& ! body length of GPS reports
    nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbgqcf        ! threshold quality control flags      (see below: 'nb?qcf')
  use data_obs_record    , only :   &
                  ! --> bit pattern parameters in ODR body format -->
    nvru       ,& ! bit pos. for status/QC flags for horiz. wind  nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature        "
    nvrq       ,& ! bit pos. for status/QC flags for humidity           "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height    "
    nvrw       ,& ! bit pos. for status/QC flags for vertical wind      "
    nvriwv     ,& ! bit pos. for status/QC flags for IWV                "
    nvrzpd     ,& ! bit pos. for status/QC flags for zenith path delay  "
!   nvrspd     ,& ! bit pos. for status/QC flags for slant path delay   "
    nvrct      ,& ! bit pos. for status/QC flags for (total) cloud      "
    nvrcl      ,& ! bit pos. for status/QC flags for low     cloud      "
    nvrcm      ,& ! bit pos. for status/QC flags for middle  cloud      "
    nvrch      ,& ! bit pos. for status/QC flags for high    cloud      "
    nvrcbs     ,& ! bit pos. for status/QC flags for cloud base height  "
    nvrcil     ,& ! bit pos. for status/QC flags for cloud ceiling      "
    nvrvis     ,& ! bit pos. for status/QC flags for visibility         "
!   nvrzbc     ,& ! bit pos. for temporary flag: QC ag. LBC pressure    "
!   nvrqbc     ,& ! bit pos. for temporary flag: QC against LBC IWV     "
    nvfubp     ,& ! bit pos. for main flag on wind                    nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature               "
    nvfqbp     ,& ! bit pos. for main flag on humidity                  "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.        "
    nvfgbp     ,& ! bit pos. for main flag on integr. water vapour      "
!   nvfaoc     ,& ! no. of bits occ. by each main flag                  "
    nvfbps     ,& ! bit pattern for main flags:                         "
    nvfboc     ,& ! no. of bits occ. for main flags                     "
    nvflbp     ,& ! bit pos. for level flag: level below surface        "
!   nvfloc     ,& ! no. of bits occ. by level flag                      "
    nvffbp     ,& ! bit pos. for flag: tower obs operator for wind    "
                  !                    disagrees with 'dhosag'
    nvfmbp     ,& ! bit pos. for flag: tower obs operator for mass    "
                  !                    (T, qv) disagrees with 'dhosag'
    nvlidp     ,& ! level id. bit pattern                             nb?lid
!   nvlido     ,& ! no. bits occ. by each indicator in level id.        "
    nvw0bp     ,& ! bit position for ww (present wea.) code (WMO Table 020003)
    nvw0oc     ,& ! no. bits occupied for ww code           (WMO Table 020003)
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
!   nxclbp     ,& ! bit position for cloud amount      code (WMO Table 020011)
!   nxctbp     ,& ! bit position for cloud type        code (WMO Table 020012)
!   nxbsbp     ,& ! bit position for cloud base height                 [m]
    nxsibp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxcloc     ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc     ,& ! no. bits occupied for cloud type   code (WMO Table 020012)
    nxbsoc     ,& ! no. bits occupied for cloud base height            [m]
    nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
    ntxbp      ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp      ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc         ! no. of bits occ. by each time period              "
  use data_obs_record    , only :   &
                  ! --> other parameters -->
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    nibits     ,& ! masking constants
                  ! --> ODR arrays and counter variables -->
    omlbdy     ,& ! body of multi-level ODR
    momlbd     ,& ! body of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    osgbdy     ,& ! body of single-level ODR
    mosgbd     ,& ! body of single-level ODR
    osghed     ,& ! header of single-level ODR
    mosghd     ,& ! header of single-level ODR
    ogpbdy     ,& ! body of GPS ODR
    mogpbd     ,& ! body of GPS ODR
    ogphed     ,& ! header of GPS ODR
    mogphd     ,& ! header of GPS ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd     ,& ! header of GPS ODR
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    ntotgp        ! tot. number of stored GPS reports
  use data_obs_record    , only :   &
                  ! format of SOR (Simulated Observation Record)
    mxsoml     ,& ! SOR body length for multi-level reports
    mxsosg     ,& ! SOR body length for single-level reports
    mxsogp     ,& ! SOR body length for GPS reports
!   mxsops     ,& ! SOR header length for multi-level reports
    nso_u      ,& ! u wind component                               [m/s]
    nso_v      ,& ! v wind component                               [m/s]
    nso_t      ,& ! temperature                                    [K]
    nso_rh     ,& ! relative humidity                              [%]
    nso_p      ,& ! pressure (sfc.) | geopotential (upper-air)  [Pa] | [m2/s2]
    nso_ct     ,& ! total cloud cover                              [octas]
    nso_cl     ,& ! low cloud cover                                [octas]
    nso_cm     ,& ! mid-level cloud cover                          [octas]
    nso_ch     ,& ! high cloud cover                               [octas]
    nso_cbs    ,& ! cloud base height CBH above surface (AGL)      [m]
    nso_cil    ,& ! ceiling (CBH of lowest cloud layer >= 5 octas) [m]
    nso_vis    ,& ! visibility                                     [m]
    nso_iq     ,& ! integrated water vapour (increment)            [mm]
!   nso_ps     ,& ! pressure                                       [Pa]
    nso_ff     ,&
    nso_dd     ,&
    nso_td

  !---------------------------------------------------------------------
  ! access namelist parameters: quality control, use of obs
  !---------------------------------------------------------------------
  use mo_cosmo_obs_data , only :   &
!      quality control for multi-level data
    qcc        ,& !  0.,500: constant parts of the quality control thresholds
                  !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf       ,& !  5., 1.: multiplication factor to the vertically varying
                  ! 10., 0.  part of the QCT (as def. in 'data_nudge_local')
!      quality control for surface-level data
    qccsu      ,& ! 12.,500: constant parts of the quality control thresholds
                  ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
!      QC for integrated water vapour (IWV) derived from GPS or radiosonde
!   qcciq      ,& !  1.   constant part of QC threshold for IWV
!   qcsiq      ,& !  .15  IWV QC threshold, as a fraction of IWV of saturated
                  !       model temperature profile
!      use of observations
    doromx     ,& ! SYNOP obs. with height differences betw. model orography
                  !   and station height >'doromx' are set passive
    fplev_ps   ,& !  1.   factor to 'plevel' (for localisation in LETKF)
                  !         for surface pressure obs from surface stations
    ilocv_sfc     !  0    mode for setting 'plevel' in surface reports
                  !        = 0 : observed station pressure only if present
                  !              (else set to missing value)
                  !        = 1 : set to model surface pressure instead
                  !              (only) if station pressure not observed
                  !        = 2 : always set to model surface pressure

  !---------------------------------------------------------------------
  ! access observation operator, quality control and auxilliary routines
  !---------------------------------------------------------------------
  use src_obs_cdfin_util , only :   &
    ncfeedobs_status   ! function: returns report status in fdbk format

  use src_obs_operator_conv , only :   &
    q2rh_col            ,& ! computes relative humidity from specific humidity
    tvirt_col           ,& ! computes model columns of virtual temperature
    hhl_col             ,& ! computes model height on model half levels
    lhyphl_col          ,& ! computes LOG of hydrostatic pressure on half levels
    shift_profile       ,& ! shifts model levels for height of sensor a. ground
    sing_obs_operator   ,& ! forward obs operator for upper-air single-level obs
    surf_obs_operator   ,& ! forward obs operator for surface-level obs
    cloud_obs_operator  ,& ! forward obs operator for cloud obs
    ps_obs_operator_sing,& ! obs opr for surface pressure from a surface station
    ps_obs_operator_mult,& ! obs opr for surf pressure from a multi-level report
    mult_obs_operator   ,& ! forward observation operator for multi-level obs
    mult_obs_2_modlev   ,& ! inverse observation operator for multi-level obs
    mult_obs_operator_z    ! supplementary obs operator for multi-level T, z
!   frac2octas             ! conversion of fraction into octas

  use src_obs_qc_conv       , only :   &
    sing_quality_cntl   ,& ! quality control of individual single-level obs
    ps_quality_cntl     ,& ! quality control of (near-) surface pressure obs
    mult_obs_qc_fg      ,& ! quality control of multi-level obs
!   mult_quality_cntl   ,& ! quality control of 1 level of a multi-level report
    mult_obs_qc_dz         ! height / thickness QC for multi-level temperature
#endif

  implicit none

  !================
  ! Public entities
  !================
  private
  public :: process_cosmo_conv ! specific tasks for COSMO observation operators
  public :: present_clc        ! flag indicating if CLC is present
  public :: lqc                ! flag indicating fg-check application
  public :: octa_percent       ! convert cloud cover (%) to octa

  !-----------------
  ! module variables
  !-----------------
  logical :: present_clc = .false. ! set to indicate if CLC is present
  logical :: lqc = .true.          ! set to apply fg-check

!===============================================================================
contains
!===============================================================================
  subroutine process_cosmo_conv (task, spot, obs, atm, cols, xi, y, Jo, Jo_atm, &
                                 state)
  integer            ,intent(in)             :: task    ! what to do
  type(t_spot)       ,intent(inout),optional :: spot    ! SPOT observations
  type(t_obs_block)  ,intent(inout),optional :: obs     ! observation data type
  type(t_atm)        ,intent(in)             :: atm     ! atmospheric state
  type(t_cols)       ,intent(in)   ,optional :: cols    ! model columns
  type(t_vector_segm),intent(in)   ,optional :: xi      ! interpolated values
  type(t_vector_segm),intent(inout),optional :: y       ! model equivalents
  real(wp)           ,intent(inout),optional :: Jo      ! obs. cost funct. Jo
  type(t_atm)        ,intent(inout),optional :: Jo_atm  ! gradient:d Jo/d atm
  integer            ,intent(in)   ,optional :: state   ! status flag
  !-----------------------------------------------------------------------
  ! This subroutine is called from various points in the assimilation code
  ! in order to perform specific tasks for the COSMO observation
  ! operators ported from COSMO to DACE.
  !------------------------------------------------------------------------
    integer                  :: tsk          ! task, local copy
    integer                  :: nimcol       ! # of model columns required
    type (t_grid)   ,pointer :: grid         ! atmospheric grid
    integer                  :: idx (1,4)    ! model indices & processor
    integer                  :: ixmcol(1)    ! indices to imcol
    integer                  :: ix, iy       ! horizontal model indices
    integer                  :: i, j, k, n   ! indices
    integer                  :: natm         ! number of model columns required
    integer(i8)              :: iatm         ! model columns reqrd. (bitfield)
    real(wp)        ,pointer :: yr (:)       ! model equivalents for a report
    type(t_datum)   ,pointer :: bd (:)       ! pointer to body data
!   logical                  :: change       ! argument to shrink_report
!   integer                  :: nlev         ! argument to shrink_report
    integer                  :: idx_(4,4)    ! Grid indices (fallback)
    real(wp)                 :: w   (4)      ! Weights
    integer                  :: np           ! number of neighbors

    !==============================
    ! observation non_specific part (called once)
    !==============================
    tsk = task
!   IF (MOD( dace% pe,10 ) == 0)  &
!   write(6,*) 'ZZD ', dace% pe, ' ZZD1', tsk, TSK_INIT, TSK_READ, TSK_SETUP_COLS, TSK_SET_CHR &
!                                            , TSK_SETUP_FUL0, TSK_SETUP_FULL, TSK_Y, TSK_R
    !------------------------------------------
    ! tsk == TSK_INIT:
    ! module initialisation, namelist read etc.
    !------------------------------------------
    if (iand (TSK_INIT,tsk) /= 0) then
      call init_cosmo_obs
      tsk=tsk-TSK_INIT
    endif
    if (tsk==0) return

    !-----------------------------
    ! tsk == TSK_READ:
    ! read COSMO observation files
    !-----------------------------
    if (iand (TSK_READ,tsk) /= 0) then
      if (read_cdfin) then
        !--------------
        ! read the file
        !--------------
        call init_fdbk_tables
        call read_cosmo_obs (atm% grid, atm)
#ifndef NO_COSMO_OBS
        !----------------------------------------------------
        ! store in observation data types
        ! seperatly for: muli-level, single-level, GPSGB data
        ! +++ currently multi level data only +++
        !----------------------------------------------------
!   write(6,*) 'ZZ4 ', dace% pe, ntotml, ntotsg, ntotgp
        call cosmo_to_obs (obs% o, 1, ntotml, omlhed, momlhd, yomlhd ) ! ml
        call cosmo_to_obs (obs% o, 2, ntotsg, osghed, mosghd, yosghd ) ! sl
        call cosmo_to_obs (obs% o, 3, ntotgp, ogphed, mogphd, yogphd ) ! gpsgb
#else
  call finish ('cosmo_to_obs','not in source code')
#endif
      endif
      tsk=tsk-TSK_READ
    endif
    if (tsk==0) return

    !========================== -------------------------
    ! observation specific part (called for every report)
    !========================== -------------------------

    !---------------------------------------------------------------
    ! tsk == SETUP_COLS:
    ! specify the input date (model column) required by the operator
    !---------------------------------------------------------------
    if (iand (TSK_SETUP_COLS,tsk) /= 0) then

!     write(0,*) dace% pe, ' ZZ6 ', spot% statid, spot% o% n, spot% col% nlev  &
!                        , obs% o% varno(spot% o% i + 1)                       &
!                        , obs% o% body (spot% o% i + 1)% use% flags
      grid => atm% grid
      select case (spot% hd% obstype)
      !----------------------------------------------
      ! possibly observation type specific processing
      !----------------------------------------------
      case (OT_TEMP,  OT_PILOT,           &! multilevel data
            OT_AIREP, OT_SATOB,           &! upper-air single level data
            OT_SYNOP, OT_SCATT, OT_DRIBU, &! single level data
            OT_GPSGB                      )! ground based GPS
        !----------------------------
        ! specify parameters required
        !----------------------------
        select case (spot% hd% obstype)
        case (OT_SYNOP)
          iatm        = COL_T + COL_Q   + &! parameters required
                        COL_P + COL_GEO + &!
                        COL_RANGE          ! (does not count in natm)
          natm        = 4                  ! number of parameters
        case default
          iatm        = COL_T + COL_Q   + &! parameters required
                        COL_UV          + &!
                        COL_P + COL_GEO    !
          natm        = 6                  ! number of parameters
        end select
        if (present_clc) then
          iatm      = iatm + COL_CLC
          natm      = natm + 1
        endif
        !----------------------------------
        ! specify the model column required
        !----------------------------------
        spot% n_spt = 1
        spot% mke   = grid% nz
        ix          = spot% col% h% ijdp(1)
        iy          = spot% col% h% ijdp(2)
        if (ix > 0 .and. iy > 0) then
          idx (1,1) = grid% marr (2,ix,iy,1)   ! ix
          idx (1,2) = grid% marr (3,ix,iy,1)   ! iy
          idx (1,3) =                     1    ! diamond
          idx (1,4) = grid% marr (1,ix,iy,1)   ! pe
          np        = 1
        else
          !----------------------------------------------
          ! Fallback if index_x, index_y not properly set
          !----------------------------------------------
          call Grid_Indices             &
               (spot% col%c% dlon,      & ! longitude
                spot% col%c% dlat,      & ! latitude
                             grid,      & ! grid data type
                             idx_(:,:), & ! Grid point indices
                             w   (:),   & ! Weight
                             np,        & ! # of coefficients
                             order = 1  ) ! nearest neighbor
          select case (np)
          case (1)
             idx(1,:) = idx_(1,:)
          case (0)
             call decr_rpt_use (spot, CHK_DOMAIN)
          case default
             write(0,*) "process_cosmo_conv: lat,lon,np =", &
                  spot% col%c% dlat, spot% col%c% dlon, np
             call finish ("process_cosmo_conv","np/=1")
          end select
        end if
        !------------------------------------------------------
        ! store model column info in the 3dvar/LETKF structures
        !------------------------------------------------------
        if (np == 1) then
           nimcol = 0
           call alloc_imcol (spot% imcol, 1)
           call add_index (idx, 1, obs% o, iatm, natm,        &
                           spot% imcol, nimcol, ixmcol, grid, &
                           spot% i_time, spot% w_time         )
        end if
      !-----------------------
      ! unknown obstype: abort
      !-----------------------
      case default
        write(0,*)  'process_cosmo_conv',                        &
                    'TSK_SETUP_COLS not implemented for obstype',&
                     spot% hd% obstype
        call finish('process_cosmo_conv',                          &
                    'TSK_SETUP_COLS not implemented for obstype '//&
                     name_value (obstype, spot% hd% obstype)       )
      end select
      tsk=tsk-TSK_SETUP_COLS
    endif
    if (tsk==0) return

    !----------------------------------------------------------------
    ! tsk == TSK_SET_CHR: set characteristics of observation operator
    !----------------------------------------------------------------
    if (iand (TSK_SET_CHR,tsk) /= 0) then
      select case (spot% hd% obstype)
      !----------------
      ! multilevel data
      !----------------
      case (OT_TEMP, OT_PILOT)
        spot% int_type  = ITY_MCOLS        ! requires model columns as input
        spot% cost      = spot% col% nlev  ! estimate of cost
        spot% nr        = spot% o%n        ! diagonal R so far
        spot% char      = CHR_NONL+CHR_EXP ! characteristics of operator
!       if(count(obs%o% varno (spot%o%i+1:spot%o%i+spot%o%n) == VN_Z) > 1) &
!       spot% nr      = spot% o%n * spot% o%n    ! geop. correlations only
      case (OT_SYNOP, OT_AIREP, OT_SATOB, OT_SCATT, OT_GPSGB, OT_DRIBU)
        spot% int_type  = ITY_MCOLS        ! requires model columns as input
        spot% cost      = 1._wp            ! estimate of cost
        spot% nr        = spot% o%n        ! diagonal R so far
        spot% char      = CHR_NONL         ! characteristics of operator
      !-----------------------
      ! unknown obstype: abort
      !-----------------------
      case default
        write(0,*)  'process_cosmo_conv',                       &
                    'TSK_SET_CHR not implemented for obstype',  &
                     spot% hd% obstype
        call finish('process_cosmo_conv',                       &
                    'TSK_SET_CHR not implemented for obstype'// &
                     name_value (obstype, spot% hd% obstype)    )
      end select
      tsk=tsk-TSK_SET_CHR
    endif
    if (tsk == 0) return

    !----------------------------------------------------------
    ! tsk == TSK_SETUP_FUL0: specify size of PSAS-space
    ! (currently not used. required for 3dvar, not verification
    !----------------------------------------------------------
    if (iand (TSK_SETUP_FUL0,tsk) /= 0) then
      call new_int (obs% o, spot, 0)
      tsk=tsk-TSK_SETUP_FUL0
    endif
    if (tsk == 0) return

    !----------------------------------------------------------
    ! tsk == TSK_SETUP_FULL: setup description of PSAS-space
    ! ++++ TO BE ADAPTED FOR COSMO VERTICAL GRID ++++
    ! (currently not used. required for 3dvar, not verification
    !----------------------------------------------------------
    if (iand (TSK_SETUP_FULL,tsk) /= 0) then
      tsk=tsk-TSK_SETUP_FULL
    endif
    if(tsk==0) return

    !------------------------------------------
    ! tsk == TSK_SHRINK:
    ! does not work yet, ignore so far
    !------------------------------------------
    if (iand (TSK_SHRINK,tsk) /= 0) then
!     call shrink_report (spot, obs%o, state, change, nl=nlev)
!     if (change) spot% col% nlev = nlev
      tsk=tsk-TSK_SHRINK
    endif
    if (tsk==0) return

    !--------------------------------------------
    ! tsk == TSK_Y: evaluate observation operator
    !--------------------------------------------
    if (iand (TSK_Y,tsk) /= 0) then

      !-------------------------------------
      ! set all model equivalents to missing
      !-------------------------------------
      yr => y% x         (spot% o% i + 1 : spot% o% i + spot% o% n)
      bd => obs% o% body (spot% o% i + 1 : spot% o% i + spot% o% n)
      yr = rvind

#ifndef NO_COSMO_OBS
      !------------------------------
      ! call the observation operator
      !------------------------------
      call run_operator (spot, obs, cols, yr, spot% mke)
#else
      call finish ('process_cosmo_conv','TSK_Y: not in source code')
#endif

      !------------------------------------------------------------
      ! set NO_OBS status if observation operator is not applicable
      !------------------------------------------------------------
      do i = 1, spot% o% n
        if (yr (i) == rvind) &
          call decr_use (bd(i)% use, STAT_OBS_ONLY, check=CHK_NONE)
      end do

      tsk=tsk-TSK_Y

! write(6,*) 'ZZ9 ',spot% statid, spot% hd% obstype, spot% col% nlev

!call finish('process_cosmo_conv','TSK_Y')

    endif
    if(tsk==0) return


    !-----------------------------------------------
    ! tsk == TSK_R
    ! set up R (observation error covariance matrix)
    ! take over values read and stored into body
    !-----------------------------------------------
    if (iand (TSK_R,tsk) /= 0) then
      if (obs% o% pe == dace% pe) then
        n = spot% o% n
        i = spot% o% i
        k = obs% R% ia (i + 1)
        do j=1,n
          obs% R% ia (i + j) = k
          obs% R% packed (k) = obs% o% body(i+j)% eo ** 2
          obs% R% ja (k) = i + j
          k = k + 1
        end do
        obs% R% ia (i + n + 1) = k
      endif
      tsk = tsk - TSK_R
    endif
    if (tsk == 0) return


    write(0,*)'process_cosmo_conv:','TSK =',tsk
    call finish('process_cosmo_conv','TSK_?')

  end subroutine process_cosmo_conv

!===============================================================================

  elemental function octa_percent (clc) result (octa)
  real(wp) ,intent(in) :: clc  ! cloud cover
  real(wp)             :: octa ! cloud amount
  !----------------------------------------------------------------------
  ! convert cloud cover (%) to cloud amount (octa)
  !
  ! WMO table 020011: CLOUD AMOUNT
  !
  ! code  meaning (synoptic stations)     meaning (climatologic stations)
  ! ----  -------  --------                        ------------
  !    0  0                               0
  !    1  1 OKTA OR LESS, BUT NOT ZERO    1/10 OR LESS, BUT NOT ZERO
  !    2  2 OKTAS                         2/10 - 3/10
  !    3  3 OKTAS                         4/10
  !    4  4 OKTAS                         5/10
  !    5  5 OKTAS                         6/10
  !    6  6 OKTAS                         7/10 - 8/10
  !    7  7 OKTAS OR MORE,BUT NOT 8 OKTAS 9/10 OR MORE, BUT NOT 10/10
  !    8  8 OKTAS                         10/10
  !    9  Sky obscured by fog and or other meteorological phenomena
  !   10  Sky partially obscured by fog and/or other meteorological phenomena
  !   11  Scattered
  !   12  Broken
  !   13  Few
  !   14  Reserved
  !   15  Cloud cover is indiscernible for reasons other than fog or other
  !       meteorological phenomena, or observation is not made
  !----------------------------------------------------------------------
  ! The following choices were made:
  ! - use the meaning of table 020011 for synoptic stations
  !   i.e. prefer (equal-spaced) octa stepping instead of mostly 10% steps
  ! - handle rounding issues near 0% or 100% (threshold: 0.5%)
  ! - assume valid input value >= 0 %, <= 100 %
  !----------------------------------------------------------------------
    if (     clc >= 99.5_wp) then
      octa = 8._wp
    else if (clc <   0.0_wp) then
      octa = rvind
    else if (clc <=  0.5_wp) then
      octa = 0._wp
    else
      octa = real (max (1, min (7, nint (0.08_wp * clc))), wp)
    endif
  end function octa_percent

!===============================================================================
#ifndef NO_COSMO_OBS

  subroutine run_operator (spot, obs, cols, y, ke)
  !--------------------------------------------------
  ! run the COSMO observation operator for one report
  !--------------------------------------------------
  type(t_spot)      ,intent(inout) :: spot    ! report header meta data
  type(t_obs_block) ,intent(inout) :: obs     ! observation data
  type(t_cols)      ,intent(in)    :: cols    ! model columns
  real(wp)          ,intent(inout) :: y (:)   ! model equivalent to observation
  integer           ,intent(in)    :: ke      ! number of model levels

    !-----------
    ! parameters
    !-----------
    integer  , parameter  :: ndqc       = 3   ! max. number of QC control
                                              !   messages per report
    integer  , parameter  :: mxvimo     = 5   ! dimension of 'vimtoob'
    logical  , parameter  :: lscadj (4) = (/ .true., .true., .true., .false. /)
                                  ! method of vertical interpolation:
                                  !   .true.  --> vertical scale adjustment
                                  !   .false. --> linear interpol. (in log( p ))
    !----------------------
    ! 3dvar/LETKF meta data
    !----------------------
    type(t_datum) ,pointer :: bd (:)  ! pointer to body data
    integer                :: nobs    ! number of observations in the report
    integer                :: nlev    ! number of vertical levels in report
    integer                :: ic      ! pointer to model column
    integer                :: ic2     ! pointer to model column (2nd time slot)
    integer                :: icn     ! pointer to model column (next time slot)
    real(wp)               :: w1, w2  ! time interpolation weights
    integer                :: kobtyp  ! CMA observation type
    integer                :: kcdtyp  ! CMA observation code type
    real(wp)               :: zstalt  ! station altitude          [m]
    real(wp)               :: zoblat  ! latitude                  [deg]
    real(wp)               :: zoblon  ! longitude                 [deg]
    character(len=8)       :: ystid   ! station identity
    !-------------------------------------------------
    ! model fields:  input to the observation operator
    !-------------------------------------------------
    real(wp)     :: col_t       (ke)  ! temperature                [  K  ]
    real(wp)     :: col_qv      (ke)  ! specific humidity          [kg/kg]
    real(wp)     :: col_u       (ke)  ! zonal  u-wind component    [ m/s ]
    real(wp)     :: col_v       (ke)  ! merid. v-wind component    [ m/s ]
    real(wp)     :: col_lnp     (ke)  ! ln( pressure ) on main levels
    real(wp)     :: col_p       (ke)  ! pressure                   [ Pa  ]
    real(wp)     :: col_p2      (ke,1)! pressure (2-d array)       [ Pa  ]
    real(wp)     :: col_z       (ke)  ! geopotential height        [  m  ]
    real(wp)     :: col_clc     (ke)  ! cloud fraction stratiform + convection
    real(wp)     :: ps                ! surface pressure           [ Pa  ]
    real(wp)     :: u_10m             ! 10 m wind component        [ m/s ]
    real(wp)     :: v_10m             ! 10 m wind component        [ m/s ]
    real(wp)     :: t_2m              !  2 m temperature           [  K  ]
    real(wp)     :: td_2m             !  2 m dewpoint temperature  [  K  ]
    real(wp)     :: hsurf             ! model orography            [  m  ]
    !---------------------------------------
    ! variables valid for a given time range
    !---------------------------------------
    real(wp)     :: vmax_10m    (nt)  ! max. 10m wind speed
    real(wp)     :: tmin_2m     (nt)  ! min.  2m temperature
    real(wp)     :: tmax_2m     (nt)  ! max.  2m temperature
    real(wp)     :: tot_prec    (nt)  ! total precipitation
    real(wp)     :: aswdir_s    (nt)  ! downw.direct SW rad
    real(wp)     :: aswdifd_s   (nt)  ! down. diffusive r.
    !!CS!!!!!--------------------------------
    ! model fields:  input not handled so far
    !!CS!!!!!--------------------------------
    real(wp)     :: col_qc     (ke)   ! specific cloud water content     [kg/kg]
    real(wp)     :: col_qrs    (ke)   ! specific content of hydrometeors
                                      !   excluding cloud water          [kg/kg]
    real(wp)     :: col_hhl    (ke+1) ! geometrical height of half levels    [m]
!   real(wp)     :: col_clc_sgs  (ke) ! subgrid-scale stratiform cloud cover [ ]
!   real(wp)     :: col_clc_con  (ke) ! cloud fraction due to convection     [ ]
    !----------------------------------
    ! model fields:  derived quantities
    !----------------------------------
    real(wp)     :: col_tv     (ke)   ! virtual temperature
    real(wp)     :: col_rh     (ke)   ! relative humidity
    real(wp)     :: col_t_og   (ke)   ! observed variable related to temperature
                                      !   (temperature or virtual temperature)
                                      !   (used only for stability-dep. humidity
                                      !    QC thresholds)
    real(wp)     :: col_hyphl  (ke+1) ! hydrostatic pressure at half levels
    !--------------------
    ! temporary variables
    !--------------------
    integer      :: kfmt         ! type of report (0: upper-air single-level;
                                 !           1: multi-level; 2: surface; 3: GPS)
    integer      :: nbxqcf       ! ODR: threshold quality control flags
    integer      :: nbxerr       ! ODR: pre-processing status flags
    integer      :: nbxflg       ! ODR: main flag word
    integer      :: mxxbdy       ! ODR: body length, real part
    integer      :: mxxbdf       ! ODR: body length, integer part
    integer      :: nbxx         ! ODR: index for variable
    integer      :: nbxxer       ! ODR: index for observation error
    integer      :: nvrx         ! ODR: index for status/QC flags
    integer      :: nvfxbp       ! ODR: index for main flag
    integer      :: mxsoxx       ! SOR: body length
    integer      :: nso_x        ! SOR: index for variable
    integer      :: istatus (1)  ! feedback file observation status
    integer      :: kflags  (1)  ! feedback file observation flag word
    integer      :: kcheck  (1)  ! feedback file observation check status
    integer      :: ilevsig      ! feedback file level significance
    integer      :: mexi    (3)  ! active (+1) or passive (-1) obs exists
    integer      :: mflgqc       ! QC bit flag (not used any further)
    integer      :: nuprin       ! file unit 'nupr', if < 0 then no printing
    integer      :: kbotlev      ! number of obs levels below lowest model level
    integer      :: ktoplev      ! number of obs levels above top model level
    integer      :: mbotlv  (3)  ! number of model levels below lowest obs level
    integer      :: mtoplv  (3)  ! number of model levels below lowest obs level
    integer      :: ilvvip  (2)  ! 1: index of lowest obs level with p-z obs
                                 ! 2: index of (lower) obs level used for intpl.
    integer      :: ilvpr  (ke)  ! obs level used for vertical interpolation
    integer      :: kb           ! index of model level immediat.below obs level
    integer      :: iobt         ! obs index in 3dvar data structure
    integer      :: iobs         ! loop index over observations within this rep.
    integer      :: ilev         ! loop index over vertical obs levels
    integer      :: ivar         ! loop index over variables
    integer      :: ivrs         ! loop index over variables
    integer      :: kk           ! loop index over vertical model levels
    integer      :: nqc          ! loop index over control messages
    integer      :: ntqc         ! index of control message
    integer      :: mv           ! auxilliary variable
    integer      :: kml850       ! model level at about 850 hPa
    integer      :: kml700       ! model level at about 700 hPa
    integer      :: nlqc         ! counter for QC control messages in 'zyqc'
    real(wp)     :: zpsvob  (1)  ! obs pressure interpol. to lowest model level
    real(wp)     :: zpsvim  (1)  ! model pressure interpolated to station height
    real(wp)     :: zpsdz        ! scaled height distance between level 'ke'
                                 !   and nearest obs.
    real(wp)     :: vip     (4)  ! model values interpolated to obs level
    real(wp)     :: zlopf        ! weight for vertical interpol. to obs level
    real(wp)     :: zpp          ! pressure at observation level
    real(wp)     :: zzz          ! height     at observation level
    integer      :: ilt          ! level type at observation level
    real(wp)     :: zobdps       ! pressure increment at lowest model level
    real(wp)     :: zf850        ! 0.85 * surface pressure
    real(wp)     :: zf700        ! 0.70 * surface pressure
    real(wp)     :: tabdif       ! difference obs time - ref. time     [hour]
    real(wp)     :: tobtim       ! obs time                            [hour]
    logical      :: lev_new      ! true, if new vertical level (if kfmt <= 1)
    logical      :: lgetso       ! true, if variable is needed for simulated obs
    logical      :: lveridat     ! true, if simulated obs to be computed
    logical      :: lobinc  (3)  ! true, if interpolate obs incr., not obs val.
    logical      :: lvirt        ! true, if 'virtual' temperature observed
    logical      :: lqcdz        ! true, if height and thicknees QC to be done
    logical      :: lwrqc        ! true, if info to be written to 'zyqc(ps)' for
                                 !   QC-rejected obs
    logical      :: lprqcm  (2)  ! true, if control messages about QC
    logical      :: lprfrs  (2)  ! true, if control messages to file unit 'nupr'
    logical, save:: lfirst=.true.! true, if routine is called the first time
    logical      :: lsoz         ! compute simulated z-obs from multi-level rep.
    character(len=1) :: yeq      ! auxilliary variable for control output only
    real(wp)     :: clct         ! total cloud cover
    real(wp)     :: clcl         ! low   cloud cover
    real(wp)     :: clcm         ! mid-level cloud cover
    real(wp)     :: clch         ! high  cloud cover
    real(wp)     :: ceil         ! ceiling (cloud base height)
    real(wp)     :: vis          ! visibility
    logical      :: l_uv         ! requires u,v (in the atmosphere)
    integer      :: nupr_        ! nupr on pe where lwonl=T
    !-----------------------------
    ! temporary allocatable arrays
    !-----------------------------
    integer  , allocatable :: mzobbd (:,:)   ! ODR: obs body, integer part
    real(wp) , allocatable :: zobbdy (:,:)   ! ODR: obs body, real part
    real(wp) , allocatable :: zsobdy (:,:)   ! SOR body (simulated obs record)
    integer  , allocatable :: ixb_odr  (:)   ! index for ODR level and variable
    integer  , allocatable :: ixv_odr  (:)   ! index for ODR status / QC flag
    integer  , allocatable :: ixs_odr  (:)   ! index for SOR variable
    real(wp) , allocatable :: vimtoob(:,:)   ! model interpolated to obs levels
    real(wp) , allocatable :: viobtom(:,:)   ! obs interpolated to model levels
    real(wp) , allocatable :: zt_o     (:)   ! obs-derived (dry bulb) temperat.
    real(wp) , allocatable :: zdzob    (:)   ! profile of height obs increments
    real(wp) , allocatable :: zyqc   (:,:)   ! info for QC control messages
    real(wp) , allocatable :: zyqcps   (:)   ! info for QC control messages

!   logical      :: lwonl =.true.!
!   integer      :: nupr         ! file unit number for control output
!                                !   (if < 0 then no control output is written)
!
!
    character(len=8) :: yvar  ! part of control message (variable type)
    real(wp):: acthr = 0._wp  ! actual hour (model time)
    integer :: iyqc   ! loop index over control messages
    integer :: myqcvar! type of control message
    integer :: maxqcp ! max.  number of rejected data to be printed per timestep
    integer :: ntotqc ! total number of rejected data to be printed per timestep
    integer         , allocatable ::  myqc (:,:)  ! info on data rejected by QC
    real(wp)        , allocatable ::  oyqc (:,:)  ! info on data rejected by QC
    character(len=8), allocatable ::  yyqc   (:)  ! info on data rejected by QC
!
!
    maxqcp      = spot% col% nlev * 4  +  2
!   maxqcp      = (maxsgo  +  maxmlo *(maxmlv + ke)  +  maxgpo/2) /20  +  1
    ntotqc      = 0
    allocate ( oyqc    (maxqcp, 11) )
    allocate ( myqc    (maxqcp,  2) )
    allocate ( yyqc    (maxqcp    ) )
!
!CS: to do: end open issues !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !--------------------------------------------------
    ! select header and body meta data from 3dvar/LETKF
    !--------------------------------------------------
    bd     => obs% o% body (spot% o% i + 1 : spot% o% i + spot% o% n)
    nobs   = spot% o% n
    kobtyp = spot% hd% obstype
    kcdtyp = spot% hd% codetype
    nlev   = spot% col% nlev
    zstalt = spot% z
    zoblat = spot% col% c% dlat
    zoblon = spot% col% c% dlon
    ystid  = spot% statid

    if (nlev < 1) write(0,*) 'ZZnlev ', nlev, ystid, kobtyp, kcdtyp

    !   determine type of report
    if ((nlev >= 2) .and. ((kobtyp == OT_TEMP) .or. (kobtyp == OT_PILOT))) then
      !   multi-level obs report
      kfmt   = 1
    elseif ((kobtyp == OT_AIREP) .or. (kobtyp == OT_SATOB)) then
      !   upper-air single-level obs report
      kfmt   = 0
    elseif  (kobtyp == OT_GPSGB) then
      !   ground-based GPS (IWV) obs report
      kfmt   = 3
    elseif ((kobtyp == OT_SYNOP) .or. (kobtyp == OT_DRIBU)                     &
                                 .or. (kobtyp == OT_SCATT)) then
      !   surface-level obs report
      kfmt   = 2
    elseif ((kobtyp == OT_TEMP ) .or. (kobtyp == OT_PILOT)) then
      !   from TEMP and PILOT baloon reports, surface reports are derived
      !   in the COSMO obs reading routines and have to be distinguished
      !   from regular upper-air reports with only 1 level
      ilevsig = bd( 1  )% lev_sig
      mv      = bd(nobs)% lev_sig
      if (     (kcdtyp == OC_RA_EU) .or. (kcdtyp == OC_WP_EU)                  &
          .or. (kcdtyp == OC_RAVAD) .or. (kcdtyp == OC_PR_US)) then
        kfmt = 0
      elseif (     (ibit1( ilevsig , LS_SURFACE ) == 1)                        &
              .or. (ibit1( mv      , LS_SURFACE ) == 1)) then
        !   for convenience, level significance of only the first and last obs
        !   (instead of all obs) are checked for surface flag
        kfmt = 2
      else
        kfmt = 0
      endif
    else
      kfmt   = - 1
    endif
! write(0,*) 'ZZ8x ', ystid, ' ', nlev, kobtyp, kcdtyp, kfmt
    if (kfmt == -1) then
      !   obs type unknown or not allowed for current obs operators
                                                                          return
    endif
    !   set indices and initialise arrays related to ODR / SOR (of COSMO)
    if (kfmt == 1) then
      nbxqcf = nbtqcf
      nbxerr = nbterr
      nbxflg = nbtflg
      mxxbdy = mxrbdy
      mxxbdf = mxrbdf
      mxsoxx = mxsoml
    elseif ((kfmt == 0) .or. (kfmt == 2)) then
      nbxqcf = nbsqcf
      nbxerr = nbserr
      nbxflg = nbsflg
      mxxbdy = mxsbdy
      mxxbdf = mxsbdf
      mxsoxx = mxsosg
    elseif (kfmt == 3) then
      nbxqcf = nbgqcf
      nbxerr = nbgerr
      nbxflg = nbgflg
      mxxbdy = mxgbdy
      mxxbdf = mxgbdf
      mxsoxx = mxsogp
    endif
    select case (kfmt)
    case (2, 3)
      l_uv = .false.
    case default
      l_uv = .true.
    end select
    allocate ( zobbdy (nlev,mxxbdy) )
    allocate ( mzobbd (nlev,mxxbdf) )
    allocate ( zsobdy (mxsoxx,nlev) )
    allocate ( ixb_odr (nobs) )
    allocate ( ixv_odr (nobs) )
    allocate ( ixs_odr (nobs) )
    zobbdy  (:,:) = rmdi
    mzobbd  (:,:) = 0
    zsobdy  (:,:) = rmdi
    ixb_odr   (:) = 0

    !-------------------------------------
    ! gather model columns for this report
    !-------------------------------------
    col_clc     = 0._wp
    ic          = spot% imcol(1)% imc(1)
    if (spot% w_time == 0._wp) then
         col_t     = cols% col(ic)%    t
         col_qv    = cols% col(ic)%    q
if(l_uv) col_u     = cols% col(ic)%    u
if(l_uv) col_v     = cols% col(ic)%    v
         col_z     = cols% col(ic)%    geo / gacc
         col_lnp   = cols% col(ic)%    p
         ps        = cols% col(ic)% s% ps
         t_2m      = cols% col(ic)% s% t2m
         u_10m     = cols% col(ic)% s% u10m
         v_10m     = cols% col(ic)% s% v10m
         td_2m     = cols% col(ic)% s% td2m
         hsurf     = cols% col(ic)% s% geosp / gacc
!        sso_std   = cols% col(ic)% s% ssd
         clct      = cols% col(ic)% s% clct
         clcl      = cols% col(ic)% s% clcl
         clcm      = cols% col(ic)% s% clcm
         clch      = cols% col(ic)% s% clch
         ceil      = cols% col(ic)% s% ceiling
         vis       = cols% col(ic)% s% vis
         vmax_10m  = cols% col(ic)% s% vmax_10m
         tmin_2m   = cols% col(ic)% s% tmin_2m
         tmax_2m   = cols% col(ic)% s% tmax_2m
         tot_prec  = cols% col(ic)% s% tot_prec
         aswdir_s  = cols% col(ic)% s% aswdir_s
         aswdifd_s = cols% col(ic)% s% aswdifd_s
         if (present_clc) then
           col_clc =    cols% col(ic)%    clc
         endif
    else
         !-----------------------
         ! temporal interpolation
         !-----------------------
         ic2     = spot% imcol(1)% imc(2)
         w2      = spot% w_time
         w1      = 1._wp - spot% w_time
         icn     = ic2; if (w2 < 0.5_wp) icn = ic
         col_t   = w1 * cols% col(ic )%    t           &
                 + w2 * cols% col(ic2)%    t
         col_qv  = w1 * cols% col(ic )%    q           &
                 + w2 * cols% col(ic2)%    q
if(l_uv) col_u   = w1 * cols% col(ic )%    u           &
                 + w2 * cols% col(ic2)%    u
if(l_uv) col_v   = w1 * cols% col(ic )%    v           &
                 + w2 * cols% col(ic2)%    v
         col_z   =(w1 * cols% col(ic )%    geo         &
                 + w2 * cols% col(ic2)%    geo) / gacc
         col_lnp = w1 * cols% col(ic )%    p           &
                 + w2 * cols% col(ic2)%    p
         ps      = w1 * cols% col(ic )% s% ps          &
                 + w2 * cols% col(ic2)% s% ps
         t_2m    = w1 * cols% col(ic )% s% t2m         &
                 + w2 * cols% col(ic2)% s% t2m
         u_10m   = w1 * cols% col(ic )% s% u10m        &
                 + w2 * cols% col(ic2)% s% u10m
         v_10m   = w1 * cols% col(ic )% s% v10m        &
                 + w2 * cols% col(ic2)% s% v10m
         td_2m   = w1 * cols% col(ic )% s% td2m        &
                 + w2 * cols% col(ic2)% s% td2m
         hsurf   =(w1 * cols% col(ic )% s% geosp       &
                 + w2 * cols% col(ic2)% s% geosp) / gacc
         clct    = w1 * cols% col(ic )% s% clct        &
                 + w2 * cols% col(ic2)% s% clct
         clcl    = w1 * cols% col(ic )% s% clcl        &
                 + w2 * cols% col(ic2)% s% clcl
         clcm    = w1 * cols% col(ic )% s% clcm        &
                 + w2 * cols% col(ic2)% s% clcm
         clch    = w1 * cols% col(ic )% s% clch        &
                 + w2 * cols% col(ic2)% s% clch
         ! no interpolation for ceiling, visibility, CBH, ...
!        ceil    = w1 * cols% col(ic )% s% ceiling     &
!                + w2 * cols% col(ic2)% s% ceiling
         ceil         = cols% col(icn)% s% ceiling
         vis          = cols% col(icn)% s% vis
         vmax_10m     = cols% col(icn)% s% vmax_10m
         tmin_2m      = cols% col(icn)% s% tmin_2m
         tmax_2m      = cols% col(icn)% s% tmax_2m
         tot_prec     = cols% col(icn)% s% tot_prec
         aswdir_s     = cols% col(icn)% s% aswdir_s
         aswdifd_s    = cols% col(icn)% s% aswdifd_s
         if (present_clc) then
           col_clc = w1 * cols% col(ic )%    clc       &
                   + w2 * cols% col(ic2)%    clc
      !--------------------
      ! keep invalid values (values: see atm2col in mo_t_col.f90)
      !--------------------
         endif
         if (cols% col(ic2)% s% t2m     <= 0._wp) t_2m    =   0._wp
         if (cols% col(ic2)% s% td2m    <= 0._wp) td_2m   =   0._wp
         if (cols% col(ic2)% s% clct    <  0._wp) clct    = -99._wp
         if (cols% col(ic2)% s% clcl    <  0._wp) clcl    = -99._wp
         if (cols% col(ic2)% s% clch    <  0._wp) clch    = -99._wp
         if (cols% col(ic2)% s% clcm    <  0._wp) clcm    = -99._wp
         if (cols% col(ic2)% s% ceiling <  0._wp) ceil    = -99._wp
         if (cols% col(ic2)% s% vis     <  0._wp) vis     = -99._wp
    endif
    if (cols% col(ic )% s% t2m     <= 0._wp) t_2m    =   0._wp
    if (cols% col(ic )% s% td2m    <= 0._wp) td_2m   =   0._wp
    if (cols% col(ic )% s% clct    <  0._wp) clct    = -99._wp
    if (cols% col(ic )% s% clcl    <  0._wp) clcl    = -99._wp
    if (cols% col(ic )% s% clch    <  0._wp) clch    = -99._wp
    if (cols% col(ic )% s% clcm    <  0._wp) clcm    = -99._wp
    if (cols% col(ic )% s% ceiling <  0._wp) ceil    = -99._wp
    if (cols% col(ic )% s% vis     <  0._wp) vis     = -99._wp
    col_p   = exp (col_lnp)
    col_p2 (:,1)  =  col_p(:)

    !---------------------------------------------------------------------
    ! Missing values in the model boundary region produce repeated values
    ! of pressure within a column and do not allow vertical interpolation.
    ! Catch these and return without further calculation.
    !---------------------------------------------------------------------
    if (any (col_lnp(1:ke-1) == col_lnp(2:ke))) return

! if (ystid(1:5) == '10394') &
! write(0   ,*) 'ZZ1a ', ystid, ' ', nlev, iobs, kcdtyp, col_p(ke), col_u(ke)   &
!                                                      , col_v(ke), col_t(ke)
! if (ystid(1:5) == '10394') &
! write(0   ,*) 'ZZ1b ', ystid, ' ', nlev, iobs, kcdtyp, col_p(ke-1), col_u(ke-1)   &
!                                                      , col_v(ke-1), col_t(ke-1)
! if (ystid(1:5) == '10394') &
! write(0   ,*) 'ZZ1c ', ystid, ' ', nlev, iobs, kcdtyp, col_p(1), col_u(1)   &
!                                                      , col_v(1), col_t(1)
!CS: to do
    !   temporary setting
    col_qc  =      0._wp
    col_qrs =      0._wp
!CS: to do end

    !-------------------------------------------------------------------------
    ! For tower profiles (in which the standard-level bit is set (for all
    ! observations)), model equivalents should be computed at the same
    ! 'height of sensor above ground' rather than at the same geometric height
    ! or pressure as the observations.
    ! For this purpose, the height and pressure levels of the model column are
    ! shifted vertically (as a first step of the forward operator) to account
    ! for the difference between station height and model orography, and
    ! a height correction is applied to temperature (and to specific humidity
    ! by preserving relative humidity).
    ! (Subsequently, the model equivalent can be computed like for multi-level
    ! observations by vertical interpolation to the (apparently) observed
    ! pressure level.)
    !-------------------------------------------------------------------------
    if ((kcdtyp == OC_TOWER) .or. (kcdtyp == OC_ICOS)) then
      if (btest( bd(1)% lev_sig , LS_STANDARD )) then

        call shift_profile ( zstalt, hsurf, ke, gacc, r_d, r_v, b1,b2w, b3,b4w &
                           , col_qc + col_qrs                                  &
                           , col_lnp, col_p, col_z, col_t, col_qv )
      ! ==================
        col_p2 (:,1)  =  col_p(:)
        !  'hsurf' is also shifted as it may be used later to compute 'col_hhl'
        hsurf  =  zstalt
      endif
    endif

    !----------------------------------------------------------------
    ! reconstruct ODR obs body 'zobbdy', 'mzobbd' as far as necessary
    !----------------------------------------------------------------
    ilt = 0
    if (kfmt <= 1)  lev_new = .true.
    if (kfmt == 1)  ilev = 0
    if (kfmt /= 1)  ilev = 1
!   write(6,*) 'ZZE ', dace% pe, ' ZZE1', nobs, kfmt, kobtyp, kcdtyp, nlev, ystid

    loop_over_obs:  do iobs = 1 , nobs
!     bd     => obs% o% body (spot% o% i + 1 : spot% o% i + spot% o% n)
      iobt  =  spot% o% i  +  iobs

!     write(0,*) dace% pe, ' ZZO3 ', spot% statid, ' ', kobtyp, ilev, iobs     &
!              , obs% o% varno(iobt), obs% o% olev(iobt), bd(iobs)% use% flags

      !   determine whether obs variable is needed for simulated obs or QC
      !   (the following check should be consistent with table 'i**_varno'!)
      lgetso = .false.
      if (kfmt <= 1) then
        lgetso =      (obs% o% varno(iobt) == VN_U      )                      &
                 .or. (obs% o% varno(iobt) == VN_V      )                      &
                 .or. (obs% o% varno(iobt) == VN_T      )                      &
                 .or. (obs% o% varno(iobt) == VN_RH     )                      &
                 .or. (obs% o% varno(iobt) == VN_HEIGHT )
!                .or. (obs% o% varno(iobt) == VN_W      )
      elseif (kfmt == 2) then
        lgetso =      (obs% o% varno(iobt) == VN_U10M   )                      &
                 .or. (obs% o% varno(iobt) == VN_V10M   )                      &
                 .or. (obs% o% varno(iobt) == VN_FF     )                      &
                 .or. (obs% o% varno(iobt) == VN_DD     )                      &
                 .or. (obs% o% varno(iobt) == VN_T2M    )                      &
                 .or. (obs% o% varno(iobt) == VN_RH2M   )                      &
                 .or. (obs% o% varno(iobt) == VN_TD2M   )                      &
                 .or. (obs% o% varno(iobt) == VN_PS     )                      &
                 .or. (obs% o% varno(iobt) == VN_PTEND  )                      &
                 .or. (obs% o% varno(iobt) == VN_NH     )                      &
                 .or. (obs% o% varno(iobt) == VN_CEIL   )                      &
                 .or. (obs% o% varno(iobt) == VN_VV     )                      &
                 .or. (obs% o% varno(iobt) == VN_N_L    )                      &
                 .or. (obs% o% varno(iobt) == VN_N_M    )                      &
                 .or. (obs% o% varno(iobt) == VN_N_H    )                      &
                 .or. (obs% o% varno(iobt) == VN_N      )
      elseif (kfmt == 3) then
        lgetso = .false.
!       lgetso =      (obs% o% varno(iobt) == VN_PWC    )
      endif
      !   hereafter, only obs of type 'real' (rather then integer) and
      !              with level_typ VN_P, VN_FLEV, or VN_HEIGHT are processed
      if (      (bd(iobs)% lev_typ /= VN_P     )                               &
          .and. (bd(iobs)% lev_typ /= VN_FLEV  )                               &
          .and. (bd(iobs)% lev_typ /= VN_HEIGHT)) then
!       PRINT *, 'run_operator: ERROR lev_typ ', bd(iobs)% lev_typ, kobtyp
                                                             cycle loop_over_obs
      endif

      !   upper-air obs: determine - pressure level
      !                            - whether current obs is at a new obs level
      if (kfmt <= 1) then
        if     (bd(iobs)% lev_typ == VN_P     ) then
          zpp  =  obs% o% olev(iobt)
          zzz  =  rmdi
        elseif (bd(iobs)% lev_typ == VN_FLEV  ) then
          zpp  =  bd(iobs)% plev
          zzz  =  rmdi
        elseif (bd(iobs)% lev_typ == VN_HEIGHT) then
          zpp  =  bd(iobs)% plev
          zzz  =  obs% o% olev(iobt)
        endif
        if (kfmt == 1) then
          lev_new = (ilev == 0)
          if (ilev >= 1) then
            if (ilt /= bd(iobs)% lev_typ) then
              lev_new = .true.
            else
              if (      (bd(iobs)% lev_typ == VN_HEIGHT)                         &
                  .and. (zobbdy(ilev,nbtz) > rmdich)) then
                lev_new  =  (zzz > zobbdy(ilev,nbtz) + 0.5_wp)
              else
! AR            lev_new  =  (zpp < zobbdy(ilev,nbtp) - rprlim*100._wp)
                lev_new  =  (zpp < zobbdy(ilev,nbtp)                 )
              endif
            endif
          endif
          if (lev_new)  ilev = ilev + 1
        endif
! if ((ystid(1:5) == '10394') .and. (nlev >= 51)) &
! write(0   ,*) 'ZZ2 ', ystid, ' ', nlev, iobs, kcdtyp, bd(iobs)% lev_typ, VN_P  &
!                     , VN_HEIGHT, lev_new, zpp, zzz, obs% o% varno(iobt), ilev
!       if ((lev_new) .or. (kfmt /= 1)) then
        if (lev_new) then
          ilt = bd(iobs)% lev_typ
          if (kfmt == 1)  zobbdy (ilev,nbtp) = zpp
          if (kfmt == 0)  zobbdy (ilev,nbsp) = zpp
          !   for upper-air obs, 'zobbdy(:,nbtz)' is filled only further below,
          !   when 'varno'=VN_HEIGHT;  for 'varno'=VN_HEIGHT,
          !   'lev_typ' must be VN_P and cannot be VN_HEIGHT or VN_FLEV
        endif
      endif

      !   note: this cycling has to be done after the pressure level is computed
      !         for multi-level obs, otherwise obs operator routines would crash
      if (.not. lgetso)                                      cycle loop_over_obs

      !   surface-level obs: get height level for 'surface' pressure obs
      !   (for other 'surface' obs, height level is 'station altitude')
      if ((kfmt == 2) .and. (obs% o% varno(iobt) == VN_PS   ))                 &
        zobbdy (ilev,nbsz) = obs% o% olev(iobt)

! if (kfmt == 2)  &
!   write(6,*) 'ZZ7a ', iobt, ystid, obs% o% varno(iobt), VN_PS, VN_P, obs% o% olev(iobt)

      !   fill observations and errors in COSMO ODR body
      nbxx   = 0
      nbxxer = 0
      nvrx   = -1
      nvfxbp = -1
      if (kfmt <= 2) then
        if     (     (obs% o% varno(iobt) == VN_U     )                        &
                .or. (obs% o% varno(iobt) == VN_U10M  )) then
          nbxx   = nbsu
!         nbxxer = nbsuer
          nvfxbp = nvfubp
          nvrx   = nvru
          nso_x  = nso_u
          if (kfmt == 1)  nbxx   = nbtu
!         if (kfmt == 1)  nbxxer = nbtuer
        elseif (     (obs% o% varno(iobt) == VN_V     )                        &
                .or. (obs% o% varno(iobt) == VN_V10M  )) then
          nbxx   = nbsv
!         nbxxer = nbsuer
          nvfxbp = nvfubp
          nvrx   = nvru
          nso_x  = nso_v
          if (kfmt == 1)  nbxx   = nbtv
!         if (kfmt == 1)  nbxxer = nbtuer
        elseif (     (obs% o% varno(iobt) == VN_T     )                        &
                .or. (obs% o% varno(iobt) == VN_T2M   )) then
          nbxx   = nbst
!         nbxxer = nbster
          nvfxbp = nvftbp
          nvrx   = nvrt
          nso_x  = nso_t
          if (kfmt == 1)  nbxx   = nbtt
!         if (kfmt == 1)  nbxxer = nbtter
        elseif (     (obs% o% varno(iobt) == VN_RH    )                        &
                .or. (obs% o% varno(iobt) == VN_RH2M  )) then
          nbxx   = nbsrh
          nbxxer = nbsqer
          nvfxbp = nvfqbp
          nvrx   = nvrq
          nso_x  = nso_rh
          if (kfmt == 1)  nbxx   = nbtrh
          if (kfmt == 1)  nbxxer = nbtqer
        elseif (     (obs% o% varno(iobt) == VN_HEIGHT)                        &
                .or. (obs% o% varno(iobt) == VN_PS    )) then
          nbxx   = nbsz
!         nbxxer = nbszer
          nvfxbp = nvfzbp
          nvrx   = nvrz
          nso_x  = nso_p
          if (kfmt == 2)  nbxx   = nbsp
          if (kfmt == 1)  nbxx   = nbtz
!         if (kfmt == 1)  nbxxer = nbtzer
        elseif (kfmt == 2) then
          if (obs% o% varno(iobt) == VN_PTEND )  nbxx   = nbspst
          if (obs% o% varno(iobt) == VN_NH    )  nbxx   = nbscbs
          if (obs% o% varno(iobt) == VN_CEIL  )  nbxx   = nbscil
          if (obs% o% varno(iobt) == VN_VV    )  nbxx   = nbsvis
          if (obs% o% varno(iobt) == VN_N_L   )  nbxx   = nbscl
          if (obs% o% varno(iobt) == VN_N_M   )  nbxx   = nbscm
          if (obs% o% varno(iobt) == VN_N_H   )  nbxx   = nbsch
          if (obs% o% varno(iobt) == VN_N     )  nbxx   = nbsct
          if (obs% o% varno(iobt) == VN_FF    )  nbxx   = nbsff
          if (obs% o% varno(iobt) == VN_DD    )  nbxx   = nbsdd
          if (obs% o% varno(iobt) == VN_TD2M  )  nbxx   = nbstd
          if (obs% o% varno(iobt) == VN_NH    )  nso_x  = nso_cbs
          if (obs% o% varno(iobt) == VN_CEIL  )  nso_x  = nso_cil
          if (obs% o% varno(iobt) == VN_VV    )  nso_x  = nso_vis
          if (obs% o% varno(iobt) == VN_N_L   )  nso_x  = nso_cl
          if (obs% o% varno(iobt) == VN_N_M   )  nso_x  = nso_cm
          if (obs% o% varno(iobt) == VN_N_H   )  nso_x  = nso_ch
          if (obs% o% varno(iobt) == VN_N     )  nso_x  = nso_ct
          if (obs% o% varno(iobt) == VN_FF    )  nso_x  = nso_ff
          if (obs% o% varno(iobt) == VN_DD    )  nso_x  = nso_dd
          if (obs% o% varno(iobt) == VN_TD2M  )  nso_x  = nso_td
        endif
      elseif (kfmt == 3) then
        if           (obs% o% varno(iobt) == VN_PWC   )  then
          nbxx   = nbgiwa
          nvfxbp = nvfgbp
          nvrx   = nvriwv
          nso_x  = nso_iq
        endif
      endif
      if (nbxx  == 0)                                        cycle loop_over_obs
      if (nbxx   > 0)  zobbdy (ilev,nbxx  )  =  bd(iobs)% o
      if (nbxxer > 0)  zobbdy (ilev,nbxxer)  =  bd(iobs)% eo
!     if (nbxx   > 0)  ix_bd  (ilev,nbxx  )  =  iobs
      ixb_odr (iobs)  =  (ilev - 1)* mxxbdy  +  nbxx
      ixv_odr (iobs)  =  nvrx
      ixs_odr (iobs)  =  nso_x

      if ((nvrx >= 0) .or. (nvfxbp >= 0)) then

        !   get observation flags and status
        istatus (1)  =  bd(iobs)% use% state
        kflags  (1)  =  bd(iobs)% use% flags
        kcheck  (1)  =  bd(iobs)% use% check

        !   convert into feedback file format
        call change_code (istatus, stats % code)
        call change_bits (kflags , checks% code)
        call change_code (kcheck , checks% code)

        !   fill required flags in ODR
        if (     (ibit1( kflags(1), FL_FG     ) == 1)                          &
            .or. (ibit1( kflags(1), FL_FG_LBC ) == 1))                         &
          mzobbd (ilev,nbxqcf) = ireplace( mzobbd(ilev,nbxqcf), nvrx, 1, 1, 0 )
        if (istatus(1) <= ST_ACTIVE)                                           &
          mzobbd (ilev,nbxerr) = ireplace( mzobbd(ilev,nbxerr), nvrx, 1, 1, 0 )
        !   required flags: - height flag
        !                   - multi-level obs   : blacklist
        !                   - multi-level height: no temperature obs available
        !     (see 'check_store_conv': pressure converted from height is not an
        !      obs in the feedback file --> no need to set corresponding flag)
!       if (ibit1( kflags(1), FL_DATASET  ) == 1)                              &
!         mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
!                                        , nvfxbp + nvfbps(1), 1, 1, 0 )
        if (ibit1( kflags(1), FL_BLACKLIST) == 1)                              &
          mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
                                         , nvfxbp + nvfbps(2), 1, 1, 0 )
!       if (ibit1( kflags(1), FL_GROSS    ) == 1)                              &
!         mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
!                                        , nvfxbp + nvfbps(3), 1, 1, 0 )
        if (ibit1( kflags(1), FL_HEIGHT   ) == 1)                              &
          mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
                                         , nvfxbp + nvfbps(4), 1, 1, 0 )
        if (      (ibit1( kflags(1), FL_PRACTICE) == 1) .and. (kfmt == 1)      &
            .and. (nvfxbp == nvfzbp))                                          &
          mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
                                         , nvfxbp + nvfbps(5), 1, 1, 0 )
!       if (      (ibit1( kflags(1), FL_NO_OBS  ) == 1) .and. (kfmt == 1)      &
!           .and. (nvfxbp == nvfzbp))                                          &
!         mzobbd (ilev,nbxflg) = ireplace( mzobbd(ilev,nbxflg)                 &
!                                        , nvfxbp + nvfbps(6), 1, 1, 0 )
        !   level significance: surface level of a multi-level report
        ilevsig = bd(iobs)% lev_sig
        if ((ibit1( ilevsig, LS_SURFACE ) == 1) .and. (kfmt == 1))             &
          mzobbd (ilev,nbtlid) = ishft( 1 , nvlidp(7) )
      endif

!     if (kfmt == 1)  lev_new = .false.
      if (kfmt <= 1)  lev_new = .false.
    enddo  loop_over_obs

    !-------------------------------------------------------------------------
    ! check consistency between:
    !   nlev: number of different levels specified in observation derived type
    !   ilev: number of different recognised here by comparing level values
    !-------------------------------------------------------------------------
    if (kfmt == 1 .and. ilev /= nlev) then
      write(0,*) 'run_operator: ilev /= nlev; kfmt,ilev,nlev,ystid=',kfmt, ilev, nlev, ystid
      call finish('run_operator','ilev/=nlev')
    endif

    if (kfmt == 1) then
      do ilev = 1 , nlev
        zobbdy (ilev,nbtlop)  =  LOG( zobbdy(ilev,nbtp) )
      enddo
      mexi = 0
      do ilev = 1 , nlev
        if (zobbdy(ilev,nbtu ) > rmdich)  mexi (1) = -1
        if (zobbdy(ilev,nbtt ) > rmdich)  mexi (2) = -1
        if (zobbdy(ilev,nbtrh) > rmdich)  mexi (3) = -1
      enddo
      do ilev = 1 , nlev
        if (ibit1( mzobbd(ilev,nbxerr), nvru ) == 1)  mexi (1) = 1
        if (ibit1( mzobbd(ilev,nbxerr), nvrt ) == 1)  mexi (2) = 1
        if (ibit1( mzobbd(ilev,nbxerr), nvrq ) == 1)  mexi (3) = 1
      enddo
    endif

    !----------------------------------------------------------------
    ! apply observation operators and quality control
    !----------------------------------------------------------------
    lveridat = .true.

!CS: to do: fill correct values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lsoz   = .true.   ! compute simulated z-obs from multi-level reports
    tabdif = 0._wp
    tobtim = 0._wp
!CS: to do: fill correct values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   if ((lqc) .and. (lfirst)) then
    if (lfirst) then
      do kk = 1 , nqclev
        tabqclp (kk) = LOG( tabqcp(kk) )
      enddo
      do kk = 1 , nqclev
        do mv = 1 , nqclev
          qczcorl(kk,mv) = nzcorl(kk,mv) * 0.001_wp
        enddo
      enddo
    endif

    if (lqc) then
      allocate ( zyqc   ( nlev*4 +2 , 11 )   )
      allocate ( zyqcps (             11 )   )
      zyqc (:,:)  = 0._wp
      zyqcps (:)  = 0._wp
      nlqc        = 0
    endif

    if (kfmt == 1) then

      !   forward observation operator for surface pressure from multi-level obs
      !   ----------------------------------------------------------------------

      call ps_obs_operator_mult ( ke, 1, col_p2(:,1:1), col_z, col_t, col_qv   &
                                , col_qc, col_qrs, r_d, gacc, rdv, doromx(2)   &
                                , nlev, zobbdy, mzobbd, .false., lveridat      &
                                , zsobdy , zpsvob, zpsdz, zpsvim, ilvvip )
!     =========================

      if (lwonl .and. nupr >=0) then
        write( nupr,'(A,", level with lowest z obs",I3                         &
                    &,", height diff ",F5.0,", oro",F6.0,", ps",F8.0)' )       &
               ystid, ilvvip(1), zpsdz, col_z(ke), col_p2(ke,1)
        if (zpsdz < -epsy)                                                     &
          write( nupr,'("hs-obs < hs-mod : ",A ,",lower level of obs.",I3      &
                      &,",dh-min",F5.0,",ipol. ps-obs",F8.0)' )                &
                 ystid, ilvvip(2), zpsdz, zpsvob(1)
        if (zpsdz >  epsy)                                                     &
          write( nupr,'("hs-obs > hs-mod : ",A ,",dh",F5.0                     &
                      &,",p-obs orig.",F7.0,",ps: mod/obs",2F8.0)' )           &
                 ystid, zpsdz, zobbdy(ilvvip(1),nbtp), col_p2(ke,1), zpsvob(1)
      endif

      !   f.g. check for surface pressure
      !   -------------------------------
      if ((lqc) .and. (zpsvob(1) > rmdich)) then
        lwrqc       = (kobtyp == OT_TEMP)

        call ps_quality_cntl ( 1, zobbdy(ilvvip(1):ilvvip(1),nbtp)             &
                             , zpsvim(1:1), .false., kobtyp, tabdif            &
                             , qcc(2), 1._wp, rmdi, lwrqc, mxrbdf              &
                             , mzobbd(ilvvip(1),1:mxrbdf), mflgqc , zyqcps )
!       ====================
      endif

      !   forward observation operator for multi-level wind, T(v), RH
      !   -----------------------------------------------------------
      lvirt  =  (     (kcdtyp == OC_RA_EU) .or. (kcdtyp == OC_WP_EU)           &
                 .or. (kcdtyp == OC_RAVAD) .or. (kcdtyp == OC_PR_US))

      lwrqc      = (lqc)
      lprqcm (1) = (lqc) .and. (qcvf(4) > epsy) .and. (lwonl)
!     lprqcm (2) = (ystid(1:5) =='11520') .or. ((io == ionl) .and. (jo == jonl))
!     lprfrs (1) = (lfirst) .and. (lwonl)
!     lprfrs (2) = (lprfrs(1)) .and. (io == ionl) .and. (jo == jonl)
      lprqcm (2) = (ystid(1:5) =='11520')
      lprfrs (1) = (lfirst) .and. (lwonl)
      lprfrs (2) = (lprfrs(1))
      nuprin = -1
      if (lprfrs(2))  nuprin = nupr
      !---------------------------------------------------
      ! pass nupr to subroutines only on pes where lwonl=T
      !---------------------------------------------------
      nupr_  = -1; if (lwonl) nupr_ = nupr

      allocate ( vimtoob (nlev,mxvimo) )
      allocate ( viobtom (ke  ,mxvimo) )
      allocate ( zt_o    (nlev)        )
      allocate ( zdzob   (nlev)        )

      call q2rh_col ( ke, col_t, col_qv, col_qc, col_p                         &
                    , rdv, b1, b2w, b3, b4w, nupr_, col_rh )
!     =============

      if ((lvirt) .or. (lqc) .or. (lsoz))                                      &

        call tvirt_col ( ke, col_t, col_qv, col_qc, col_qrs, rdv , col_tv )
!       ==============
      if (      lvirt)  col_t_og (:)  =  col_tv(:)
      if (.not. lvirt)  col_t_og (:)  =  col_t (:)

      call mult_obs_operator ( ke, col_u, col_v, col_t_og, col_rh, col_p,col_z &
                             , nlev, zobbdy, mzobbd, kobtyp, lveridat, lqc     &
                             , nuprin , zsobdy , vimtoob                       &
                             , lvirt, col_t , zt_o, kbotlev, ktoplev )
!     ======================

      !   f.g. check for multi-level wind, T, RH
      !   --------------------------------------
      if (lqc) then
        !   by setting ktoplev = nlev, multi-level check is skipped for towers
        if ((kcdtyp == OC_TOWER) .or. (kcdtyp == OC_ICOS)) ktoplev = nlev

        call mult_obs_qc_fg ( nlev, zobbdy, zt_o, kbotlev, ktoplev, kobtyp     &
                            , tabdif, mexi, zoblat, qcvf, r_d, gacc, lqc       &
                            , lwrqc, lprqcm, nupr_, ystid , mzobbd, vimtoob    &
                            , nlqc, zyqc(1:4*nlev,:), qcc )
!       ===================

        if (nupr >= 0) then
        if (lwonl) then
          if (lprfrs(1)) write( nupr,'("Sta. ",A ,3X,I4,": kbotlev, ktoplev "  &
                                     &,2I4)' )  ystid, kbotlev, ktoplev
        endif
        endif
      endif

      !   forward observation operator for multi-level height
      !   ---------------------------------------------------
      lqcdz = .false.
      if ((lqc) .or. (lsoz)) then

        call hhl_col    ( ke, col_z, hsurf , col_hhl )
!       ============
        call lhyphl_col ( ke, col_p, col_tv, col_hhl, r_d, gacc , col_hyphl )
!       ===============

        !  required: inverse obs operator for T, RH

        call mult_obs_2_modlev ( ke, col_u, col_v, col_t_og, col_rh, col_p     &
                               , col_hyphl, nlev, zobbdy, mzobbd, vimtoob      &
                               , kobtyp, nuprin, 2, 3, lscadj(3:4)             &
                               , viobtom(1:ke,3:4), mbotlv(2:3), mtoplv(2:3)   &
                               , lobinc(2:3), ilvpr )
!       ======================

        !   control output
        if (nupr >= 0) then
        if (lwonl) then
          if (lprfrs(1)) then
            write( nupr,'("sta. ",A,2X,": mbotlv, lobinc",4I4,3(2X,L1))' )       &
                   ystid, (mbotlv(ivrs),ivrs=2,3)                                &
                        , (mtoplv(ivrs),ivrs=2,3), (lobinc(ivrs),ivrs=1,3)
            mv = ke - min( mbotlv(2), mbotlv(3) )
            do kk = mv , 1 , -1
              write( nupr,'(16X,F8.1,F6.2,F8.0,I4)' )                            &
                    (viobtom(kk,ivar), ivar=3,4), col_p(kk), ilvpr(kk)
            enddo
          endif
        endif
        endif

        zobdps    =  0._wp
        if (zpsvob(1) > rmdich)  zobdps  =  ABS( zpsvob(1) ) - col_p(ke)

        call mult_obs_operator_z ( ke, col_t, col_p, col_z, col_tv, col_hhl    &
                                 , nlev, zobbdy, mzobbd, zobdps                &
                                 , viobtom(1:ke,3:3), mbotlv(2), mtoplv(2)     &
                                 , kobtyp, r_d, gacc, lveridat                 &
                                 , zsobdy, zdzob , lqcdz )
      ! ========================

      !   f.g. check of height and thickness for multi-level T and z
      !   ----------------------------------------------------------
        if ((lqc) .and. (lqcdz)) then
        !   (zdzob = zobbdy(.,nbtz) - zsobdy(nso_p) for non-missing values)

          call mult_obs_qc_dz ( nlev, zobbdy, zdzob, tabdif, qcvf(2), r_d, gacc  &
                              , mzobbd, zyqc(nlqc+1:nlqc+1,:) )
      ! ===================
          if (zyqc(nlqc+1,1) >= epsy)  nlqc = nlqc + 1
        endif
      endif

      deallocate ( viobtom )
      deallocate ( vimtoob )
      deallocate ( zt_o    )
      deallocate ( zdzob   )


    elseif (kfmt == 0) then

      !   forward observation operator for upper-air single-level wind, T, RH
      !   -------------------------------------------------------------------

      call sing_obs_operator ( ke, col_u, col_v, col_t, col_qv, col_qc, col_p  &
                             , zobbdy(1,:), lveridat, rdv, b1, b2w, b3, b4w    &
                             , zsobdy(:,1) , vip, kb, zlopf )
    ! ======================

      !   f.g. check for upper-air single-level wind, T, RH
      !   -------------------------------------------------
      lwrqc      = (lqc)

      if (lqc)  call sing_quality_cntl ( zobbdy(1,:), vip, ps, kb, kobtyp      &
                                       , tabdif, qcvf, qcc, ndqc, lwrqc        &
                                       , mzobbd(1,:) , nlqc, zyqc(1:ndqc,:) )
      !         ======================

    elseif (kfmt == 2) then

      !   forward observation operator for surface pressure from surface station
      !   ----------------------------------------------------------------------

      if ((zobbdy(1,nbsz) > rmdich) .and. (zobbdy(1,nbsp) > rmdich)) then

! write(0,*) 'ZZ80 ', ystid, ' ', kobtyp, kcdtyp, zobbdy(1,nbsz), col_z(ke), zobbdy(1,nbsp)

        call ps_obs_operator_sing ( ke, 1, col_p2(:,1:1), col_z, col_t, col_qv &
                                  , col_qc, col_qrs, r_d, gacc, rdv            &
                                  , zobbdy(1,:), lveridat                      &
                                  , zsobdy(:,1) , zpsvob, zpsdz, zpsvim )
      ! =========================

        if (nupr >= 0) then
        if (lwonl) then
          yeq = '='
          if (zpsdz < -epsy)  yeq = '<'
          if (zpsdz >  epsy)  yeq = '>'
          write( nupr,'("surf: hs-obs ",A," hs-mod:",A,", height diff obs-mod" &
                      &,F6.0,", ps mod/obs",2F8.0)' )                          &
                 yeq, ystid, zpsdz, col_p2(ke,1), zpsvob(1)
        endif
        endif

        !   f.g. check for surface pressure
        !   -------------------------------
        if ((lqc) .and. (zpsvob(1) > rmdich)) then
          lwrqc       = (kobtyp == OT_TEMP)

          call ps_quality_cntl ( 1, zobbdy(1:1,nbsp), zpsvim(1:1), .true.      &
                               , kobtyp, tabdif, qccsu(2), 1._wp               &
                               , zobbdy(1,nbspst), lwrqc , mxsbdf              &
                               , mzobbd(1,1:mxsbdf), mflgqc , zyqcps )
        ! ====================
        endif
      endif

      !   forward observation operator for screen-level wind, T, RH
      !   ---------------------------------------------------------
      zf850      = pbot_lapse / 1000._wp * ps
      zf700      = ptop_lapse / 1000._wp * ps
      kml850 = ke
      kml700 = ke - 1
      do kk = 1 , ke - 1
        if ((col_p(kk) <= zf850) .and. (col_p(kk+1) > zf850))  kml850 = kk
        if ((col_p(kk) <= zf700) .and. (col_p(kk+1) > zf700))  kml700 = kk
      enddo
      kml850 = max( kml850 , kml700 + 1 )

      call surf_obs_operator ( ke, u_10m, v_10m, t_2m, td_2m, col_t, col_z, ps &
                             , hsurf, zobbdy(1,1:mxsbdy), zstalt               &
                             , kml700, kml850, lveridat, rdv, b1, b2w, b3, b4w &
                             , zsobdy(:,1) , vip )
    ! ======================

!CS: to be done !!!
      !   forward observation operator for cloud cover
      !   --------------------------------------------
      if (.true.) then
!     if (lrad) then
!       col_clc (:)  =  col_clc_sgs(:) + col_clc_con(:) *(c1 - col_clc_sgs(:))

        call hhl_col  ( ke, col_z, hsurf , col_hhl )
      ! ============

        clct = octa_percent (clct)
        clcl = octa_percent (clcl)
        clcm = octa_percent (clcm)
        clch = octa_percent (clch)
        call cloud_obs_operator ( ke, col_clc,  clct, clcl, clcm, clch         &
                                , col_hhl, zobbdy(1,:) , zsobdy(:,1)           &
                                , ceil, vis )
      ! =======================
      endif

      !   f.g. check for screen-level (10-m) wind, (2-m) T and RH
      !   -------------------------------------------------------
      lwrqc      = (lqc)

      if (lqc)  call sing_quality_cntl ( zobbdy(1,:), vip, ps, -1, kobtyp      &
                                       , tabdif, qcvf, qccsu, ndqc, lwrqc      &
                                       , mzobbd(1,:) , nlqc, zyqc(1:ndqc,:) )
      !         ======================

    endif

    !----------------------------------------------------------
    ! observation operator for variables valid for a time range
    ! (for verification with MEC)
    !----------------------------------------------------------
    call apply_trange (spot, obs, y, vmax_10m, tmin_2m, tmax_2m,   &
                                     tot_prec, aswdir_s, aswdifd_s )

    !   fill record for later printing for control of QC
    !   ------------------------------------------------
    if (lqc) then
      if (((kfmt == 1) .or. (kfmt == 2)) .and. (zyqcps(1) >= epsy)) then
        nlqc = nlqc + 1
        zyqc (nlqc,:)  =  zyqcps(:)
      endif
      !   adjust 'nlqc' if space in arrays 'oyqc' is insufficient
      nlqc    = MIN( maxqcp - ntotqc , nlqc )
      do nqc = 1 , nlqc
        ntqc  =  ntotqc + nqc
        yyqc (ntqc   ) = ystid
        myqc (ntqc, 1) = kcdtyp
        myqc (ntqc, 2) = NINT( zyqc(nqc,1) )
        oyqc (ntqc, 1) = tobtim
        oyqc (ntqc, 2) = zyqc (nqc, 2)
        oyqc (ntqc, 3) = zoblat
        oyqc (ntqc, 4) = zoblon
        oyqc (ntqc, 5) = zyqc (nqc, 5)
        oyqc (ntqc, 6) = zyqc (nqc, 6)
        oyqc (ntqc, 7) = zyqc (nqc, 7)
        oyqc (ntqc, 8) = zyqc (nqc, 8)
        oyqc (ntqc, 9) = zyqc (nqc, 9)
        oyqc (ntqc,10) = zyqc (nqc,10)
        oyqc (ntqc,11) = zyqc (nqc,11)
      enddo
      ntotqc  =  ntotqc  +  nlqc

      deallocate ( zyqc   )
      deallocate ( zyqcps )
    endif


    !----------------------------------------------------------------
    ! convey simulated observations and f.g.check rejection flag
    ! from ODR to 3DVAR data structure
    !----------------------------------------------------------------
    do iobs = 1 , nobs
      if (ixb_odr(iobs) > 0) then
        ilev   =  (ixb_odr(iobs) - 1) /mxxbdy  +  1
        nbxx   =   ixb_odr(iobs)  -  (ilev - 1)* mxxbdy
        nvrx   =   ixv_odr(iobs)
        nso_x  =   ixs_odr(iobs)
        iobt   =  spot% o% i + iobs
        !------------------------------------
        ! skip if operator was not applicable
        !------------------------------------
        if (obs% o% varno(iobt) == VN_T2M  .and. t_2m  <= 0._wp) cycle
        if (obs% o% varno(iobt) == VN_RH2M .and. t_2m  <= 0._wp) cycle
        if (obs% o% varno(iobt) == VN_RH2M .and. td_2m <= 0._wp) cycle
        if (obs% o% varno(iobt) == VN_TD2M .and. t_2m  <= 0._wp) cycle
        if (obs% o% varno(iobt) == VN_TD2M .and. td_2m <= 0._wp) cycle
        !-----------------------------------------------
        ! finally set model equivalent and fg check flag
        !-----------------------------------------------
        if ((nso_x > 0) .and. (zsobdy(max(nso_x,1),ilev) > rmdich))            &
          !   if simulated obs exists, pass it to 't_obs'
          !   (note: y => y% x (spot% o% i + 1 : spot% o% i + spot% o% n)
          !          --> use index 'iobs' rather than 'iobt' for 'y')
          y (iobs)  =  zsobdy(nso_x,ilev)
        if (lqc) then
          !-------------------------------------------
          ! handle some passively monitored quantities
          !-------------------------------------------
          select case (obs% o% varno(iobt))
          case (VN_FF, VN_DD)
            nvrx = nvru
          case (VN_TD2M)
            nvrx = nvrq
          end select
          if (nvrx >= 0) then
            !----------------------------------------------------
            ! if f.g.check flag in ODR is set, pass it to 't_obs'
            !----------------------------------------------------
            if (ibit1( mzobbd(ilev,nbxqcf), nvrx ) == 1)                         &
              call decr_use ( bd(iobs)% use, STAT_REJECTED, check=CHK_FG )
          endif
        endif
      endif
    enddo
    ystid  = spot% statid

!
! debug printout
!
!if (lfirst) &
!print *,'## OBS ##',col_t(1),col_t(ke),col_qv(1),col_qv(ke),col_p(1),col_p(ke),col_z(1),col_z(ke),ps,t_2m,u_10m,td_2m,hsurf
!lfirst = .false.
    deallocate ( zobbdy  )
    deallocate ( mzobbd  )
    deallocate ( zsobdy  )
    deallocate ( ixb_odr )
    deallocate ( ixv_odr )
    deallocate ( ixs_odr )
    lfirst = .false.

!CS: to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   nupr = 6
    if (nupr >= 0) then
    if (lwonl) then
      do iyqc = 1 , ntotqc
        myqcvar = ABS( myqc(iyqc,2) )
        if (myqcvar ==  1) yvar = 'uv    : '
        if (myqcvar ==  2) yvar = 'p-TEMP: '
        if (myqcvar ==  3) yvar = 'T     : '
        if (myqcvar ==  4) yvar = 'RH    : '
        if (myqcvar ==  5) yvar = 'uv-10m: '
        if (myqcvar ==  6) yvar = 'ps    : '
        if (myqcvar ==  7) yvar = 'T -2m : '
        if (myqcvar ==  8) yvar = 'RH-2m : '
        if (myqcvar ==  9) yvar = 'z     : '
        if (myqcvar == 10) yvar = 'dz    : '
        if (myqcvar == 11) yvar = 'V-mult: '
        if (myqcvar == 12) yvar = 'T-mult: '
        if (myqcvar == 13) yvar = 'q-mult: '
        if (myqcvar == 14) yvar = 'z-mult: '
        if (myqcvar == 16) yvar = 'IWV   : '
        if (myqc(iyqc,1) == 137) cycle   ! no print for RADAR VAD
! horizontal wind
        if ((myqcvar == 1) .OR. (myqcvar == 5)) then
          write( nupr,'(A,''-QC: '',A ,I4,'' Obs/Mod/Thr'',F5.1,3F6.1,F5.1     &
                      &,'', Time O/M'',2F5.1  ,'', P '',F5.0)' )               &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1)                      &
               , oyqc(iyqc,6), oyqc(iyqc,8), oyqc(iyqc,7), oyqc(iyqc,9)        &
               , oyqc(iyqc,5), oyqc(iyqc,1), acthr       , oyqc(iyqc,2)
! 'surface' pressure
        elseif ((myqcvar == 2) .OR. (myqcvar == 6)) then
          write( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',2F8.2,F6.2      &
                      &,'', Time Obs/Mod'',2F6.1)' )                           &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)        &
               , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr
! temperature
        elseif ((myqcvar == 3) .OR. (myqcvar == 7)) then
          write( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.1           &
                      &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )           &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)        &
               , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
! relative humidity
        elseif ((myqcvar == 4) .OR. (myqcvar == 8)) then
          write( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.2           &
                      &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )           &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)        &
               , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
! hydrostatic height or thickness
        elseif ((myqcvar == 9) .OR. (myqcvar == 10)) then
          write( nupr,'(A ,''-QC: '',A ,I4,'' Differ./Thresh'',2F6.1           &
                      &,'', Time Obs/Mod'',2F6.1,'', P '',F5.0,''-'',F5.0)')   &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,7)        &
               , oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2), oyqc(iyqc,6)
! multi-level check
        elseif ((myqcvar >= 11) .AND. (myqcvar <= 14)) then
          write( nupr,'(A ,''-QC: '',A ,I4, 27X                                &
                      &,'', Time Obs/Mod'',2F6.1,'', P '',F5.0,''-'',F5.0)')   &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,1)        &
               , acthr       , oyqc(iyqc,2), oyqc(iyqc,5)
! integrated water vapour
        elseif (myqcvar == 16) then
          write( nupr,'(A ,''-QC: '',A ,I4,'' Obs/Mod/Thresh'',3F6.2           &
                      &,'', Time Obs/Mod'',2F6.1  ,'', P '',F5.0)' )           &
                 yvar(1:6)   , yyqc(iyqc)  , myqc(iyqc,1), oyqc(iyqc,6)        &
               , oyqc(iyqc,7), oyqc(iyqc,5), oyqc(iyqc,1), acthr, oyqc(iyqc,2)
        endif
      enddo
    endif
    endif
    deallocate ( oyqc    )
    deallocate ( myqc    )
    deallocate ( yyqc    )
!CS: to do !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine run_operator

!===============================================================================

  subroutine cosmo_to_obs  ( obs, kfmt, nrep, odrhed, modrhd, yodrhd )
  !-------------------------------------------------------------------
  type (t_obs)     ,intent(inout) :: obs          ! observations data type
  integer          ,intent(in)    :: kfmt         ! report type:
                                                  !   = 1 : multi-level
                                                  !   = 2 : single-level
                                                  !   = 3 : ground-based GPS
  integer          ,intent(in)    :: nrep         ! number of reports
  real(wp)         ,intent(in)    :: odrhed (:,:) ! COSMO ODR header (real part)
  integer          ,intent(in)    :: modrhd (:,:) ! COSMO ODR header (int  part)
  character(len=*) ,intent(in)    :: yodrhd (:)   ! COSMO ODR header (char part)
  !---------------------------------------------------------------------
  ! convert COSMO obs data structure 'ODR' for conventional observations
  ! to 3dvar/LETKF obs data structure 't_obs'
  !---------------------------------------------------------------------
    integer            :: ir          ! report index
    integer            :: ss          ! time difference (decoding - obs) [s]
    integer            :: kphase      ! aircraft flight phase
    integer            :: ktyp        ! report type (similar to kfmt)
    integer            :: r_state (1) ! report status
    integer            :: r_flags (1) ! report flags
    integer            :: r_check (1) ! report flag number changing the status
    integer            :: nlev        ! number of vertical levels
    type(t_spot)       :: spt         ! report meta data variable
    type(t_spot)       :: empty       ! default initialised
    type(t_head)       :: head        ! data usually stored in BUFR header
    type(t_use)        :: use         ! state of the report
    type(t_time)       :: t_obtime    ! actual observation time
    type(t_time)       :: t_diff      ! decoding time - obs time
!   integer ,parameter :: imiss = 2147483647  ! missing value
    !
    !----------------- end of subroutine header --------------------------------

    !------------------
    ! loop over reports
    !------------------
    do ir = 1, nrep
      !--------------------------
      ! convert basic header data
      !--------------------------
      head% modtype     = COSMO                ! module to handle this report
      head% obstype     = modrhd(ir,nhobtp)    ! CMA observation type
      head% codetype    = modrhd(ir,nhcode)    ! CMA obs code type
      head% buf_type    = modrhd(ir,nhcat )    ! BUFR type (data category)
      head% buf_subtype = modrhd(ir,nhcats)    ! BUFR type (data sub-category)
      head% dbkz        = modrhd(ir,nhkz  )    ! DWD internal class (optional)
      if (modrhd(ir,nhkz) == imdi)  head% dbkz = -1  !  fill value
      head% center      = MOD( modrhd(ir,nhcent) , 1000)   ! originating centre
      head% subcenter   =      modrhd(ir,nhcent) / 1000    ! origin. sub-centre

      ! absolute actual observation time: 't_obtime'
      call init_time ( t_obtime   , yyyymmdd = modrhd(ir,nhdate)               &
                                  , hhmmss   = modrhd(ir,nhhrmn)*100 )
      ! for absolute synoptic obs time, use 't_obtime' and time difference
      ss = NINT( 3600._wp* (odrhed(ir,nhsynt) - odrhed(ir,nhtime)) )
      call init_time ( t_diff     , ss = ss )
      head% time        = t_obtime + t_diff    ! synoptic obs decoding time
      ! for absolute data base decoding time, use 't_obtime' and time difference
      !   (note: data base decoding time is optional, may have fill value !)
      if (odrhed(ir,nhtddb) > rmdich) then
        ss = NINT( 3600._wp* (odrhed(ir,nhtddb) - odrhed(ir,nhtime)) )
        call init_time ( t_diff     , ss = ss )
        head% db_time     = t_obtime + t_diff    ! data base decoding time
      else
        head% db_time     = zero_time
      endif
!     ! for year of synoptic time, add century of obs time
!     iyyyyob = modrhd(ir,nhdate) / 10000
!     iyyyysy = modrhd(ir,nhsyhr) / 1000000 + MOD( iyyyyob, 100 ) * 100
!     ! correct century, if synoptic time and obs time are not in same century
!     if (ABS( MOD( iyyyyob, 100 ) - MOD( iyyyysy, 100 ) ) == 99) then
!       if (MOD( iyyyysy, 100 ) == 99)  iyyyysy = iyyyysy - 100
!       if (MOD( iyyyysy, 100 ) ==  0)  iyyyysy = iyyyysy + 100
!     endif
!     call init_time ( head% time , yyyymmdd =   modrhd(ir,nhsyhr) /100        &
!                                              + 1000000* (iyyyysy /100)
!                                 , hhmmss   = MOD(  modrhd(ir,nhsyhr), 100 )  &
!                                                  * 10000 )

!!    head% source      =                      ! source file number (optional)
!!    head% record      =                      ! record in file     (optional)
      !-------------------------
      ! perform some basic tests
      !-------------------------
      use = use_0
      call check_report_0 (use, head, 1)

!     IF (MOD( dace% pe,40 ) == 0)  &
!       write(6,*) 'ZZX4 ', dace% pe, ir, modrhd(ir,nhobtp), modrhd(ir,nhcode) &
!                                       , odrhed(ir,nhtime), modrhd(ir,nhpass) &
!                                       , modrhd(ir,nhflag), modrhd(ir,nhkz  ) &
!                                       , use% state, STAT_DISMISS
      if (use% state <= STAT_DISMISS) cycle
      !------------------------------
      ! convert remaining header data
      !------------------------------
      spt = empty
      spt% use                 = use
      spt% hd                  = head
!!    spt% hd% satid           =         ! satellite Id (optional)
      spt% statid       = yodrhd(ir)           ! station id (alphanumeric)
!     if (modrhd(ir,nhstid) /= imdi)                                           &
        spt% ident      = modrhd(ir,nhstid)    ! station id (numeric)
!     spt% stlsf        =                      ! station land / sea flag
!     spt% stclf        =                      ! station cloud flag
!     spt% stzen        =                      ! satellite zenith angle
      spt% sozen        = odrhed(ir,nhsolz)    ! solar zenith angle
      spt% ssd_bg       = odrhed(ir,nhssos)    ! std dev. of model SSO
      if (odrhed(ir,nhalt) > rmdich)                                           &
        spt% z          = odrhed(ir,nhalt )    ! station altitude
      if (kfmt == 1) then
        spt% sttyp      = modrhd(ir,nhrtyp)    ! radiosonde type
        spt% tracking   = modrhd(ir,nhtrac)    ! tracking technique
        spt% stret      = modrhd(ir,nhrad )    ! radiation correction
        spt% meas_type  = modrhd(ir,nhna4 )    ! measuring equipment (002003)
      elseif (kfmt == 2) then
        if (head% obstype == OT_SATOB) then
          spt% stret    = modrhd(ir,nhstyp)    ! wind computation method
        else
          spt% sttyp    = modrhd(ir,nhstyp)    ! station type
        endif
      endif
!!!!  spt% ps           = ?????????????????    ! pressure at station height ???
      if (head% obstype == OT_AIREP) then
        kphase  =  ibits( modrhd(ir,nhschr), nvapbp, nvapoc-1 )
        if ((kphase /= nibits(nvapoc-1)) .and. (kphase /= 0)) then
          spt% phase    = kphase               ! aircraft flight phase
          if (ibit1( modrhd(ir,nhschr), nvapbp+nvapoc-1 ) == 1)                &
            spt% tracking  = 18
        endif
      endif
      spt% actual_time  = t_obtime             ! actual observation time
!     if (modrhd(ir,nhcorr) /= imdi)                                           &
!       spt% corme      = MIN( 1 , modrhd(ir,nhcorr)         ! correction report
      spt% corme        = ibit1( modrhd(ir,nhschr), nvscbp ) ! correction report
      !------------------------------------------------------
      !    state, flags, check: first convert to fdbk format,
      !                         then to 't_obs'
      !------------------------------------------------------
      r_flags (1)       = modrhd(ir,nhflag)
      if ((modrhd(ir,nhpass) == 1) .and. (ibit1( r_flags(1), FL_MERGE ) == 0)) &
        r_flags (1)     = ior( r_flags(1) , ishft( 1 , FL_MERGE ) )
      r_check (1)       = bit1_pos_nr      ( r_flags(1) , flags% n             &
                                           , flags% e(1:flags% n)% value )
      r_state (1)       = ncfeedobs_status ( modrhd(ir,nhpass)                 &
                                           , modrhd(ir,nhflag) )
      call reverse_code ( r_state , stats )
      call reverse_bits ( r_flags , checks , unknown = CHK_CONSIST )
      call reverse_code ( r_check , checks , unknown = CHK_CONSIST )
!     !---------------------------------------------------------------
!     ! overwrite state/check/flags from DA system with COSMO settings
!     !---------------------------------------------------------------
!     spt% use% state      = r_state(1)
!     spt% use% flags      = r_flags(1)
!     spt% use% check      = r_check(1)
!!    spt% use% check      = CHK_DATASET
      !-------------------------------------------------
      ! merge state/check/flags from COSMO and DA system
      !-------------------------------------------------
      spt%   use% flags = ior (spt% use% flags, r_flags(1))
      if (r_state(1) <= spt% use% state) then
        spt% use% state = r_state(1)
        spt% use% check = r_check(1)
 !      spt% use% check = CHK_DATASET
      endif
      !------------------------------------------------
      ! station coordinates and associated model column
      !------------------------------------------------
      spt% col% c% dlat    = odrhed(ir,nhjlat)    ! station latitude (degree)
      spt% col% c% dlon    = odrhed(ir,nhilon)    ! station longitude (degree)
      spt% col% h% ijdp(1) = modrhd(ir,nhitot)    ! \ global indices (lon /lat)
      spt% col% h% ijdp(2) = modrhd(ir,nhjtot)    ! / of grid pt assigned to obs
      spt% col% h% ijdp(3) = 1                    ! diamond index (icosah. grid)
      if (kfmt == 1)  spt% col% nlev   = modrhd(ir,nhnlev)    ! number of levels
      if (kfmt /= 1)  spt% col% nlev   = 1                    ! number of levels
      !   the following 2 assignments are done in 'mo_psasutil.f90'
!     spt% gp_bg     = gacc* odrhed(ir,nhsurf)    ! model surface geopotential
!     call frac_mdlsfc ( spt, ibit1( modrhd(ir,nhschr), nvsebp ) )
      call set_xuv (spt)                 ! derive quantities from coordinates
      !-------------------------------------------
      ! convert COSMO missing value to LETKF/3DVAR
      !-------------------------------------------
      if (spt% ident     == imdi ) spt% ident     = empty% ident
!     if (spt% stret     == imiss) spt% stret     = empty% stret
      if (spt% stret     == imdi ) spt% stret     = empty% stret
      if (spt% sttyp     == imdi ) spt% sttyp     = empty% sttyp
      if (spt% tracking  == imdi ) spt% tracking  = empty% tracking
      if (spt% meas_type == imdi ) spt% meas_type = empty% meas_type
      if (spt% ssd_bg    == imdi ) spt% ssd_bg    = empty% ssd_bg
      !------------------
      ! convert body data
      !------------------
! write(6,*) 'ZZ5 ',spt% statid, ir, kfmt, spt% hd% obstype, spt% col% nlev
      if (kfmt == 1) then
        nlev = spt% col% nlev
        call check_store_conv (spt, obs, kfmt, nlev, omlbdy(ir,1:nlev,:)       &
                                                   , momlbd(ir,1:nlev,:) )
      elseif (kfmt == 2) then
        ktyp = kfmt
        if (     (head% obstype == OT_AIREP)                                   &
            .or. (head% obstype == OT_SATOB))  ktyp = 0
        call check_store_conv (spt, obs, ktyp, 1   , osgbdy(ir:ir,    :)       &
                                                   , mosgbd(ir:ir,    :) )
      elseif (kfmt == 3) then
        call check_store_conv (spt, obs, kfmt, 1   , ogpbdy(ir:ir,    :)       &
                                                   , mogpbd(ir:ir,    :) )
      endif
!     IF ((MOD( dace% pe,40 ) == 0) .OR. (spt% use% state >= STAT_ACTIVE))     &
!       write(6,*) 'ZZX6 ', dace% pe, spt% statid, spt% hd% obstype            &
!                         , head% codetype, spt% col% nlev                     &
!                         , r_state(1), spt% use% state                        &
!                         , r_flags(1), spt% use% flags
    end do

  end subroutine cosmo_to_obs

!===============================================================================

  subroutine check_store_conv  ( spot, obs, kfmt, nlev, zobbdy, mzobbd )
  !---------------------------------------------------------------
  type(t_spot)  ,intent(inout)  :: spot         ! header of 1 obs in 't_obs'
  type(t_obs)   ,intent(inout)  :: obs
  integer       ,intent(in)     :: kfmt         ! report type:
                                                !   = 0 : upper-air single-level
                                                !   = 1 : multi-level
                                                !   = 2 : synoptic surface level
                                                !   = 3 : ground-based GPS
  integer       ,intent(in)     :: nlev         ! number of vertical levels
  real(wp)      ,intent(in)     :: zobbdy (:,:) ! COSMO ODR body (real part)
  integer       ,intent(in)     :: mzobbd (:,:) ! COSMO ODR body (integer part)
  !--------------------------------------------------------------------
  ! convert body of COSMO obs data structure 'ODR' for conventional obs
  ! to 3dvar/LETKF obs data structure 't_obs';   perform some checks
  !--------------------------------------------------------------------

    type(t_spot) ,pointer :: spt    ! pointer to report entry
!   type(t_datum)    :: bod         ! temporary body derived type
    integer          :: id          ! temporary for observation id
    integer          :: ib          ! loop index over observations in report
    integer          :: ilev        ! loop index over vertical levels
    integer          :: ivar        ! loop index over variables (entries) in ODR
    integer          :: ibt         ! body index
    integer          :: nobs        ! number of new entries in body
    integer          :: mxodrr      ! number of variables (entries) in ODR
!   integer          :: mflgqcf     ! ODR quality control flag word
    integer          :: mflgerf     ! ODR status flag word
    integer          :: mflgmfw     ! ODR main flag word
    integer          :: nvrx        ! ODR index for status / QC flags
    integer          :: nvfxbp      ! ODR index for main flag
    integer          :: istatus (1) ! observation status
    integer          :: kflags  (1) ! observation flags
    integer          :: kcheck  (1) ! obs flag number changing the status
    integer          :: kflg        ! observation flags
    integer          :: kfbit       ! temporary buffer
    integer          :: icl         ! temporary buffer
!   integer          :: ibuf    (3) ! temporary buffer
    logical          :: lnewobs     ! new obs to be filled into 't_obs'
    logical          :: lraso       ! multi-level report is radiosonde
    integer   :: kobs (max( nlev*mxrbdy, mxsbdy+mxsbdf ) ) ! ODR variable index

    ! relation between ODR body and observation variable in fdbk format:
    ! 'ibx_varno' = 0 : quantity is not to be written as obs to 't_datum'
    !             > 0 : obs variable in fdbk format, to be written to 't_datum'
    !             < 0 : minus( obs variable in fdbk format ), to be written to
    !                   to 't_datum' with state=ST_OBS_ONLY (no simulated obs)
    !             (note: fdbk variable number 0 = NUM is never an observation)
    integer , parameter  ::  &
      ibm_varno (mxrbdy) = (/ VN_U, VN_V, VN_T, VN_RH, 0  , VN_HEIGHT         ,&
                              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0           ,&
                              VN_W, 0, 0 /)                                   ,&
      ibu_varno (mxsbdy) = (/ VN_U, VN_V, VN_T, VN_RH   , 0 , 0 , 0 , 0 , 0   ,&
                              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0   ,&
                              0 , 0 , 0 , 0 , 0 ,-VN_VGUST  , 0 , 0 , 0 , 0   ,&
                              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 /),&
      ibs_varno (mxsbdy) = (/ VN_U10M, VN_V10M, VN_T2M, VN_RH2M, VN_PS, 0     ,&
                              0 , 0 , 0 , 0 , -VN_PTEND, VN_NH, VN_CEIL       ,&
                              VN_N_L, VN_N_M, VN_N_H, VN_N, VN_VV             ,&
                              VN_FF, VN_DD, VN_TD2M                           ,&
                              -VN_RR, -VN_RR, -VN_RR, -VN_RR, -VN_RR          ,&
                              -VN_GUST, -VN_GUST, -VN_GUST, -VN_TMIN, -VN_TMAX,&
                              -VN_RAD_GL, -VN_RAD_DF, -VN_RAD_LW, -VN_SDEPTH  ,&
                              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 /)                ,&
      ibg_varno (mxgbdy) = (/ 0, -VN_ZPD, -VN_ZWD, 0, -VN_PS, -VN_T2M,-VN_RH2M,&
                              0 , VN_PWC, 0 , 0 , 0 , 0 , 0 , 0 /)            ,&
      iis_varno (mxsbdf) = (/ 0 , 0 , 0 , 0 , 0 , -VN_WW, -VN_GCLG            ,&
                              -VN_ICLG, -VN_ICLG, -VN_ICLG, -VN_ICLG, 0, 0 /) ,&
      iiu_varno (mxsbdf) = (/ 0 , 0 , 0 , 0 , 0 , -VN_TURB, 0, 0, 0, 0, 0, 0,0/)
    !
    !----------------- end of subroutine header --------------------------------
    !------------------------
    ! perform some more tests
    !------------------------
    call check_report_1 (spot)
    if (spot% use% state > STAT_DISMISS) then
      !----------------------------
      ! get number of obs in report
      !----------------------------
!     nlev  =  spt% col% nlev
      if (kfmt == 1)  mxodrr = mxrbdy
      if (kfmt == 0)  mxodrr = mxsbdy
      if (kfmt == 2)  mxodrr = mxsbdy
      if (kfmt == 3)  mxodrr = mxgbdy
      nobs  =  0
      do ilev = 1, nlev
        !   all types of reports: check real elements
        do ivar = 1, mxodrr
          if (zobbdy(ilev,ivar) > rmdich) then
            lnewobs = .false.
            if (kfmt == 1) then
              lnewobs = (ibm_varno(ivar) /= 0)
              if (ibm_varno(ivar) == VN_HEIGHT) then
                if (spot% hd% obstype == OT_AIREP) then
                  !   aircraft flight level is not an observation
                  lnewobs = (ibit1( mzobbd(ilev,nbtflg), nvfzbp+nvfbps(5)) /= 1)
                else
                  !   if bit nvfbps(6) is set then profiler height is not an obs
                  lnewobs = (ibit1( mzobbd(ilev,nbtflg), nvfzbp+nvfbps(6)) /= 1)
                endif
              endif
            elseif (kfmt == 0) then
              lnewobs = (ibu_varno(ivar) /= 0)
!             if (      (ibu_varno(ivar) == VN_HEIGHT)                         &
!                 .and. (spot% hd% obstype == OT_AIREP))                       &
!               !   aircraft flight level is not an observation
!               lnewobs = (ibit1( mzobbd(ilev,nbsflg), nvfzbp+nvfbps(5)) /= 1)
            elseif (kfmt == 2) then
              lnewobs = (ibs_varno(ivar) /= 0)
            elseif (kfmt == 3) then
              lnewobs = (ibg_varno(ivar) /= 0)
            endif
!           if ((kfmt <= 1) .and. (lnewobs) .and. () then
            if (lnewobs) then
              nobs = nobs + 1
              kobs (nobs)  =  (ilev - 1)* mxodrr  +  ivar
            endif
          endif
        enddo
        !   surface-level reports: integer obs (WW, cloud layers (not TEMPs))
        if ((kfmt == 2) .and. (spot% hd% obstype /= OT_TEMP)) then
          do ivar = 1, mxsbdf
            if ((mzobbd(ilev,ivar) /= imdi) .and. (iis_varno(ivar) /= 0)) then
              if (ivar /= nbswwe) then
                nobs = nobs + 1
                kobs (nobs)  =  - ((ilev - 1)* mxsbdf  +  ivar)
              elseif (   ibits( mzobbd(ilev,ivar), nvw0bp, nvw0oc )            &
                      /= nibits( nvw0oc )) then
                nobs = nobs + 1
                !    for integer obs, set 'kobs' to negative value
                kobs (nobs)  =  - ((ilev - 1)* mxsbdf  +  ivar)
              endif
            endif
          enddo
        !   upper-air single-level reports: integer obs (turbulence)
        elseif (kfmt == 0) then
          do ivar = 1, mxsbdf
            if ((mzobbd(ilev,ivar) /= imdi) .and. (iiu_varno(ivar) /= 0)) then
              if ((ivar == nbstur) .and. (mzobbd(ilev,ivar) /= 15)) then
                nobs = nobs + 1
                !    for integer obs, set 'kobs' to negative value
                kobs (nobs)  =  - ((ilev - 1)* mxsbdf  +  ivar)
              endif
            endif
          enddo
        endif
      enddo
      ! check if multi-level report is radiosonde
      if (kfmt == 1) then
        lraso = (     (       spot% hd%  obstype == OT_TEMP )                  &
                 .or. (      (spot% hd%  obstype == OT_PILOT)                  &
                       .and. (spot% hd% codetype /= OC_RA_EU)                  &
                       .and. (spot% hd% codetype /= OC_WP_EU)                  &
                       .and. (spot% hd% codetype /= OC_RAVAD)                  &
                       .and. (spot% hd% codetype /= OC_PR_US)))
      endif
      !---------------------------------------
      ! insert new header entry into the table
      !---------------------------------------
      call new_spot (obs, 1, set_id=.true.)
      spt => obs% spot (obs% n_spot)
      id      = spt% id
      spt     = spot
      spt% id = id
      !-----------------------------------
      ! insert body entries into the table
      !-----------------------------------
      call new_obs (obs, nobs, spot=spt)
      !------------------------
      ! set common body entries
      !------------------------
!     bod % use % state = STAT_ACTIVE
!     bod % use % check = CHK_NONE
      !-----------------
      ! set body entries
      !-----------------
      do ib = 1, nobs
!!      obs % varno (ib)      =
!!      obs %  olev (ib)      =
!       obs %  body (ib)      = bod
!!      obs %  body (ib)% o   =
!!      obs %  body (ib)% eo  =
        ibt = spt% o% i + ib
!       obs %  body (ibt)     = bod
        if (kobs(ib) > 0) then
          ilev        =  ( kobs(ib) - 1) /mxodrr  +  1
          ivar        =    kobs(ib)  -  (ilev - 1)* mxodrr
        else
          ilev        =  (-kobs(ib) - 1) /mxsbdf  +  1
          ivar        =   -kobs(ib)  -  (ilev - 1)* mxsbdf
        endif
        istatus (1)              = ST_ACTIVE
        !   pre-set missing value for level significance
        obs% body (ibt)% lev_sig = -1
        if (kfmt == 1) then
          !   multi-level obs
!         mflgqcf   = mzobbd(ilev,nbtqcf)
          mflgerf   = mzobbd(ilev,nbterr)
          mflgmfw   = mzobbd(ilev,nbtflg)
          if (ibm_varno(ivar) < 0)  istatus (1) = ST_OBS_ONLY
          obs% varno(ibt)          = ABS( ibm_varno (ivar) )
!  PRINT *, 'ZZ4 ', dace% pe, ilev, ivar, ib, ibt, ibm_varno (ivar), obs% varno(ibt), kobs(ib)
          obs% olev (ibt)          = zobbdy(ilev,nbtp)
          obs% body (ibt)% plev    = zobbdy(ilev,nbtp)
          obs% body (ibt)% lev_typ = VN_P
          if (zobbdy(ilev,nbtsnr) > rmdich)                                    &
            obs% body (ibt)% pcc     = nint( zobbdy(ilev,nbtsnr) )
          if (lraso)                                                           &
            obs% body (ibt)% lev_sig = mzobbd(ilev,nbtlsg)
          if (spot% hd% obstype == OT_AIREP) then
            obs% body (ibt)% lev_sig = spt% phase
            if (      (zobbdy(ilev,nbtz) > rmdich)                             &
                .and. (ibit1( mflgmfw, nvfzbp+nvfbps(5) ) == 1)) then
              !   aircraft reports: pressure level is converted from reported
              !                     flight level using ICAO standard atmosphere
              obs% body (ibt)% lev_typ = VN_FLEV
              obs% olev (ibt)          = zobbdy(ilev,nbtz)
              obs% body (ibt)% plev    = zobbdy(ilev,nbtp)
            endif
          elseif (ibit1( mflgmfw, nvfzbp+nvfbps(6) ) == 1) then
            !   pressure level is converted from reported height using model
            !   atmosphere, therefore height should be used as level
            !   because level should be independent from model (ensemble member)
            obs% body (ibt)% lev_typ = VN_HEIGHT
            obs% olev (ibt)          = zobbdy(ilev,nbtz)
            obs% body (ibt)% plev    = zobbdy(ilev,nbtp)
          endif
        elseif ((kfmt == 0) .and. (kobs(ib) > 0)) then
          !   upper-air single-level obs
!         mflgqcf   = mzobbd(ilev,nbsqcf)
          mflgerf   = mzobbd(ilev,nbserr)
          mflgmfw   = mzobbd(ilev,nbsflg)
          if (ibu_varno(ivar) < 0)  istatus (1) = ST_OBS_ONLY
          obs% varno(ibt)          = ABS( ibu_varno (ivar) )
          obs% olev (ibt)          = zobbdy(ilev,nbsp  )
          obs% body (ibt)% lev_typ = VN_P
!         obs% body (ibt)% lev_sig = ishft( 1 , LS_SIGN )
          if (spot% hd% obstype == OT_AIREP) then
            obs% body (ibt)% lev_sig = spt% phase
            if (      (zobbdy(ilev,nbsz) > rmdich)                             &
                .and. (ibit1( mflgmfw, nvfzbp+nvfbps(5) ) == 1)) then
              !   aircraft reports: pressure level is converted from reported
              !                     flight level using ICAO standard atmosphere
              obs% body (ibt)% lev_typ = VN_FLEV
              obs% olev (ibt)          = zobbdy(ilev,nbsz)
              obs% body (ibt)% plev    = zobbdy(ilev,nbsp)
            endif
          elseif (ibit1( mflgmfw, nvfzbp+nvfbps(6) ) == 1) then
            !   pressure level is converted from reported height using model
            !   atmosphere, therefore height should be used as level
            !   because level should be independent from model (ensemble member)
            obs% body (ibt)% lev_typ = VN_HEIGHT
            obs% olev (ibt)          = zobbdy(ilev,nbsz)
            obs% body (ibt)% plev    = zobbdy(ilev,nbsp)
          endif
        elseif ((kfmt == 2) .and. (kobs(ib) > 0)) then
          !   surface synoptic obs
!         mflgqcf   = mzobbd(ilev,nbsqcf)
          mflgerf   = mzobbd(ilev,nbserr)
          mflgmfw   = mzobbd(ilev,nbsflg)
          if (ibs_varno(ivar) < 0)  istatus (1) = ST_OBS_ONLY
          obs% varno(ibt)          = ABS( ibs_varno (ivar) )
          obs% olev (ibt)          = spt% z
          obs% body (ibt)% lev_typ = VN_HEIGHT
            !   ('plev', 'lev_sig' from TEMP-/tower-derived surface reports
            !    are treated / overruled further below)
            !   (.(nbslid) == 1 means that Synop station pressure is observed)
          if (      (zobbdy(ilev,nbsp) > rmdich) .and. (ilocv_sfc <= 1)         &
              .and. (     (mzobbd(ilev,nbslid) == 1)                            &
                     .or. (spot% hd% obstype == OT_DRIBU))) then
            !   set plev to observed station pressure
            obs% body (ibt)% plev    = zobbdy(ilev,nbsp)
          elseif ((zobbdy(ilev,nbsvip) > rmdich) .and. (ilocv_sfc >= 1)) then
            !   note that this makes 'plev' inconsistent with 'level'
            obs% body (ibt)% plev    = zobbdy(ilev,nbsvip)
          endif
          if (spot% hd% obstype == OT_TEMP)  obs% olev (ibt) = zobbdy(ilev,nbsz)

          if (      (spot% hd% obstype == OT_SYNOP)                             &
              .and. (mzobbd(ilev,nbslid) /= ibset( 0 , nvlidp(10) ))) then
            !   SYNOP: pressure code
            if (ivar == nbsp)  obs% body (ibt)% lev_sig = mzobbd(ilev,nbslid)
!         elseif (ibit1( mzobbd(ilev,nbslid), nvlidp(7) ) == 1) then
!           if (     (ivar == nbsu ) .or. (ivar == nbsv ) .or. (ivar == nbst ) &
!               .or. (ivar == nbsrh) .or. (ivar == nbsp ))                     &
!             obs% body (ibt)% lev_sig = ishft( 1 , LS_SURFACE )
          endif
!         obs% body (ibt)% pcc       = -1
        elseif (((kfmt == 2) .or. (kfmt == 0)) .and. (kobs(ib) < 0)) then
          !   integer obs
!         mflgqcf   = mzobbd(ilev,nbsqcf)
          mflgerf   = mzobbd(ilev,nbserr)
          mflgmfw   = mzobbd(ilev,nbsflg)
          if ((kfmt == 2) .and. (iis_varno(ivar) < 0))  istatus (1) = ST_OBS_ONLY
          if ((kfmt == 0) .and. (iiu_varno(ivar) < 0))  istatus (1) = ST_OBS_ONLY
          if  (kfmt == 2)   obs% varno(ibt)   = ABS( iis_varno (ivar) )
          if  (kfmt == 0)   obs% varno(ibt)   = ABS( iiu_varno (ivar) )
          obs% olev (ibt)          = 0._wp
        elseif (kfmt == 3) then
          !   ground-based GPS (IWV) obs
!         mflgqcf   = mzobbd(ilev,nbgqcf)
          mflgerf   = mzobbd(ilev,nbgerr)
          mflgmfw   = mzobbd(ilev,nbgflg)
          if (ibg_varno(ivar) < 0)  istatus (1) = ST_OBS_ONLY
          obs% varno(ibt)          = ABS( ibg_varno (ivar) )
          obs% olev (ibt)          = spt% z
          obs% body (ibt)% lev_typ = VN_HEIGHT
!         obs% body (ibt)% pcc     = -1
        endif
        if (kobs(ib) > 0)   obs% body (ibt)% o   =       zobbdy(ilev,ivar)
        if (kobs(ib) < 0)   obs% body (ibt)% o   = REAL( mzobbd(ilev,ivar), wp )
        obs% body (ibt)% bc        =  0._sp
!       obs% body (ibt)% eo        =  0._sp
!       obs% body (ibt)% ac        = -1._sp

        !   set variables that depend on observed quantity
        !   note: for all variables with 'i**_varno > 0', 'nvrx' must be >= 0
        nvrx   = -1
!       nvrx2  = -1
        nvfxbp = -1
        if (kobs(ib) > 0) then
          if (      ((kfmt == 1) .and. ((ivar == nbtu) .or. (ivar == nbtv)))   &
              .or.  ((abs( kfmt-1 ) == 1)                                      &
                               .and. ((ivar == nbsu) .or. (ivar == nbsv)))) then
            !   horizontal wind
            if (kfmt == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbtuer)
            if (kfmt /= 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbsuer)
            if ((kfmt == 1) .and. (zobbdy(ilev,nbtuac) > rmdich)) then
              if (     (spot% hd% codetype == OC_TOWER)                        &
                  .or. (spot% hd% codetype == OC_ICOS )) then
!               obs% body (ibt)% azimuth    = zobbdy(ilev,nbtuac)
                obs% body (ibt)% obs_par(1) = zobbdy(ilev,nbtuac)
              else
                obs% body (ibt)% ac         = zobbdy(ilev,nbtuac)
              endif
            endif
            nvfxbp = nvfubp
            nvrx   = nvru
            if (spot% hd% obstype == OT_AIREP)                                 &
              obs% body (ibt)% pcc     = ibits( mflgmfw, nvfubp+nvfbps(1), 2 )
            if (spot% hd% obstype == OT_SATOB)                                 &
              obs% body (ibt)% pcc     = mzobbd(ilev,nbsqal)
          elseif (      ((     kfmt     == 1) .and. (ivar == nbtt))            &
                  .or.  ((abs( kfmt-1 ) == 1) .and. (ivar == nbst))            &
                  .or.  ((     kfmt     == 3) .and. (ivar == nbgt))) then
            !   temperature
            if (     kfmt     == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbtter)
            if (abs( kfmt-1 ) == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbster)
            if (     kfmt     == 3)  obs% body (ibt)% eo   = zobbdy(ilev,nbgter)
            nvfxbp = nvftbp
            nvrx   = nvrt
          elseif (      ((     kfmt     == 1) .and. (ivar == nbtrh))           &
                  .or.  ((abs( kfmt-1 ) == 1) .and. (ivar == nbsrh))           &
                  .or.  ((     kfmt     == 3) .and. (ivar == nbgrh))) then
            !   humidity
            if (     kfmt     == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbtqer)
            if (abs( kfmt-1 ) == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbsqer)
            if (     kfmt     == 3)  obs% body (ibt)% eo   = zobbdy(ilev,nbgqer)
            if (     kfmt     == 1)  obs% body (ibt)% bc   = zobbdy(ilev,nbtdrh)
            if (abs( kfmt-1 ) == 1)  obs% body (ibt)% bc   = zobbdy(ilev,nbsdrh)
            nvfxbp = nvfqbp
            nvrx   = nvrq
!           nvrx2  = nvrqbc
          elseif (     ((kfmt == 1) .and. ((ivar == nbtp) .or.(ivar == nbtz))) &
                  .or. ((kfmt == 3) .and. ((ivar == nbgp) .or.(ivar == nbgz))) &
                  .or. ((abs( kfmt-1 ) == 1)                                   &
                               .and. ((ivar == nbsp) .or. (ivar == nbsz)))) then
            !   pressure / height  ( --> use obs height, not station height)
            if (     kfmt     == 2)  obs% olev (ibt)       = zobbdy(ilev,nbsz )
            if (     kfmt     == 3)  obs% olev (ibt)       = zobbdy(ilev,nbgz )
            if (     kfmt     == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbtzer)
            if (abs( kfmt-1 ) == 1)  obs% body (ibt)% eo   = zobbdy(ilev,nbszer)
            if (     kfmt     == 3)  obs% body (ibt)% eo   = zobbdy(ilev,nbgzer)
            nvfxbp = nvfzbp
            nvrx   = nvrz
!           nvrx2  = nvrzbc
          endif
          !   RASS (virtual) temperature bias correction if present
          if ((     (spot% hd% codetype == OC_RA_EU)                           &
               .or. (spot% hd% codetype == OC_WP_EU)) .and. (kfmt == 1)        &
                                                      .and. (ivar == nbtt)) then
            if (zobbdy(ilev,nbtdrh) > rmdich)                                  &
                                     obs% body (ibt)% bc   = zobbdy(ilev,nbtdrh)
          endif
        endif
        !   variables specific to upper-air reports:
        if ((kfmt == 1) .and. (kobs(ib) > 0) .and. (ivar == nbtw  )) then
          !   vertical velocity
          nvrx   = nvrw
        !   variables specific to ground-based GPS reports
        elseif ((kfmt == 3) .and. (kobs(ib) > 0)) then
          if (ivar == nbgiwa) then
            !   IWV
            if (zobbdy(ilev,nbgbia) > rmdich)                                  &
              obs% body (ibt)% bc   = zobbdy(ilev,nbgbia)
            !   compute observation error for IWV !
            if (zobbdy(ilev,nbgzwd) > 0.5_wp) then
              !   if ZWD exists, not close to 0 : IWV-err = ZWD-err * IWV / ZWD
              obs% body (ibt)% eo   = zobbdy(ilev,nbgtze) *zobbdy(ilev,nbgiwa) &
                                                          /zobbdy(ilev,nbgzwd)
            else
              !   0.17 = approx. ratio IWV[mm] / ZWD[mm]
              obs% body (ibt)% eo   = zobbdy(ilev,nbgtze) * 0.17_sp
            endif
            nvfxbp = nvfgbp
            nvrx   = nvriwv
!           nvrx2  = nvrqbc
          elseif (ivar == nbgzwd) then
            !   ZWD
            obs% body (ibt)% eo     = zobbdy(ilev,nbgtze)
            nvfxbp = nvfgbp
            nvrx   = nvrzpd
!           nvrx2  = nvrqbc
          elseif (ivar == nbgzpd) then
            !   ZPD
            obs% body (ibt)% eo     = zobbdy(ilev,nbgtze)
            nvfxbp = nvfgbp
            nvrx   = nvrzpd
!           nvrx2  = nvrqbc
          endif
        !   variables specific to surface reports
        elseif ((kfmt == 2) .and. (kobs(ib) > 0)) then
          if (ivar == nbsct ) then
            !   total cloud
            obs% body (ibt)% lev_sig   =  0
            icl  =  ibits( mzobbd(ilev,nbsclg), nxsgbp, nxsgoc )
            if (icl /= 63)  obs% body (ibt)% lev_sig   = icl
            nvrx   = nvrct
          elseif (ivar == nbscl ) then
            !   low cloud
            nvrx   = nvrcl
            obs% body (ibt)% lev_sig   =  7
          elseif (ivar == nbscm ) then
            !   mid-level cloud
            nvrx   = nvrcm
            obs% body (ibt)% lev_sig   =  8
          elseif (ivar == nbsch ) then
            !   high cloud
            nvrx   = nvrch
            obs% body (ibt)% lev_sig   =  9
          elseif (ivar == nbscbs) then
            !   cloud base height
            nvrx   = nvrcbs
!           obs% body (ibt)% lev_sig   =  0
          elseif (ivar == nbscil) then
            !   cloud ceiling
            nvrx   = nvrcil
!           obs% body (ibt)% lev_sig   =  0
          elseif (ivar == nbsvis) then
            !   visibility
            nvrx   = nvrvis
!           obs% body (ibt)% lev_sig   =  0
          elseif ((ivar == nbsff ) .or. (ivar == nbsdd )) then
            !   wind speed
            nvfxbp = nvfubp
            nvrx   = nvru
          elseif (ivar == nbstd ) then
            !   2-m dewpoint temperature
            nvfxbp = nvfqbp
            nvrx   = nvrq
          !   temporally non-local observations: pressure tendency, precip,
          !                                      gusts, min./max temperature
          elseif (ivar == nbspst) then
            obs% olev (ibt)            =  3._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif ((ivar == nbsrr1) .or. (ivar == nbsfg1)) then
            obs% olev (ibt)            =  1._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif ((ivar == nbsrr3) .or. (ivar == nbsfg3)) then
            obs% olev (ibt)            =  3._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif ((ivar == nbsrr6) .or. (ivar == nbsfg6)) then
            obs% olev (ibt)            =  6._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif  (ivar == nbsr12) then
            obs% olev (ibt)            =  12._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif  (ivar == nbsr24) then
            obs% olev (ibt)            =  24._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif  (ivar == nbstx) then
            icl  =  IBITS( mzobbd  (ilev, nbsttr), ntxbp, ntxoc )
            obs% olev (ibt)            =  icl
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif  (ivar == nbstn) then
            icl  =  IBITS( mzobbd  (ilev, nbsttr), ntnbp, ntxoc )
            obs% olev (ibt)            =  icl
            obs% body (ibt)% lev_typ   =  VN_TRTR
          elseif ((ivar == nbsrad) .or. (ivar == nbsrdd)                     &
                                   .or. (ivar == nbsrdl)) then
            obs% olev (ibt)            =  1._wp
            obs% body (ibt)% lev_typ   =  VN_TRTR
          endif
        elseif ((kfmt == 2) .and. (kobs(ib) < 0)) then
          if (ivar == nbswwe) then
            obs% olev (ibt)            =  0._wp
            obs% body (ibt)% lev_typ   =  VN_NUM
          elseif (ivar == nbsclg) then
            obs% olev (ibt)     =  0._wp
            obs% body (ibt)% o  =  real( iand( mzobbd(ilev,ivar)               &
                                             , nibits( nxcloc+3*nxctoc ) ) , wp)
            obs% body (ibt)% lev_typ   =  VN_NUM
            obs% body (ibt)% lev_sig   =  ibits( mzobbd(ilev,ivar) &
                                               , nxsgbp, nxsgoc    )
          elseif (     (ivar == nbscl1) .or. (ivar == nbscl2)                  &
                  .or. (ivar == nbscl3) .or. (ivar == nbscl4)) then
            if (ivar == nbscl1)   obs% olev (ibt)   =  1._wp
            if (ivar == nbscl2)   obs% olev (ibt)   =  2._wp
            if (ivar == nbscl3)   obs% olev (ibt)   =  3._wp
            if (ivar == nbscl4)   obs% olev (ibt)   =  4._wp
            obs% body (ibt)% o  = real( iand( mzobbd(ilev,ivar)                &
                                            , nibits(nxcloc+nxctoc+nxbsoc) ),wp)
            obs% body (ibt)% lev_typ   =  VN_NUM
            obs% body (ibt)% lev_sig   =  ibits( mzobbd(ilev,ivar) &
                                               , nxsibp, nxsgoc    )
          endif
        elseif ((kfmt == 0) .and. (kobs(ib) < 0)) then
          if (ivar == nbstur) then
            obs% olev (ibt)            =  0._wp
            obs% body (ibt)% lev_typ   =  VN_NUM
          endif
        endif
        if (kfmt == 2) then
          !   special treatment of 'surface reports' derived from TEMP or towers
          if (mzobbd(ilev,nbslid) == ibset( 0 , nvlidp(10) )) then
            !   this is to cope with checks in 'read_fdbk_temp' in LETKF
            if (obs% body (ibt)% lev_typ == VN_HEIGHT) then
              obs% body (ibt)% lev_sig = ishft( 1 , LS_SURFACE )
            else
              obs% body (ibt)% lev_sig = 0
            endif
            !   pressure obs can be ~100 m above the surface (HPB, OXK)
            !   --> for plev: - use observed pressure for pressure obs
            !                   or (if ilocv_sfc == 0) also for other variables
            !                 - use model surface pressure for other variables
            if ((zobbdy(ilev,nbsp) > rmdich) .and. (nvrx == nvrz)) then
              obs% body (ibt)% plev    = zobbdy(ilev,nbsp)
            elseif ((zobbdy(ilev,nbsvip) > rmdich) .and. (ilocv_sfc >= 1)) then
              obs% body (ibt)% plev    = zobbdy(ilev,nbsvip)
            elseif ((zobbdy(ilev,nbsp  ) > rmdich) .and. (ilocv_sfc <= 1)) then
              obs% body (ibt)% plev    = zobbdy(ilev,nbsp)
            endif
          endif
          !   correct plev by 'fplev_ps' for all surface pressure obs
          !    (missing value of plev is -1., see mo_t_datum)
          if ((nvrx == nvrz) .and. (obs% body (ibt)% plev > 0._wp))            &
            !   note that this makes 'plev' inconsistent with 'level'
            obs% body (ibt)% plev  =  obs% body (ibt)% plev  *  fplev_ps
        endif

        ! state, flags, check: first convert to fdbk format,
        !                      then to 't_obs'
        ! --> get obs status in fdbk format:
        !   default status is 'obs_only' for obs with 'i**_varno < 0',
        !                     'active'   for obs with 'i**_varno > 0'
        !                     (note: if 'i**_varno > 0' then 'nvrx >= 0')
        if ((istatus(1) == ST_ACTIVE) .and. (nvrx >= 0)) then
          if (ibit1( mflgerf, nvrx ) == 0) then
            !   only if (nvfxbp >= 0) there is a specific reason (check) to
            !   reject an obs --> default is 'rejected', otherwise 'passive'
            if (nvfxbp >=  0)  istatus (1) = ST_REJECTED
            if (nvfxbp == -1)  istatus (1) = ST_PASSIVE
          endif
          !   threshold quality control:  no flags set yet at this point
!         if (ibit1( mflgqcf, nvrx ) == 1) then
!           if (istatus(1) == ST_ACTIVE )  istatus (1) = ST_REJECTED
!           if (istatus(1) == ST_PASSIVE)  istatus (1) = ST_PAS_REJ
!         endif
        endif
        ! --> get obs flags in fdbk format:
        kflg = 0
        if (nvfxbp >= 0) then
          !   set dataset flag if 2-bit ODR flag = 1
          if (ibits( mflgmfw, nvfxbp + nvfbps(1), nvfboc(1) ) == 1)            &
            kflg = ishft( 1 , FL_DATASET )
!           kflg = insert( kflg, 1, FL_DATASET )
          kflg = ireplace( kflg, FL_BLACKLIST, 1, mflgmfw, nvfxbp + nvfbps(2) )
          kflg = ireplace( kflg, FL_PRACTICE , 1, mflgmfw, nvfxbp + nvfbps(5) )
          kfbit = ibit1( mflgmfw, nvfxbp + nvfbps(3) )
          if ((kfmt <= 1) .and. ((nvfxbp == nvfubp) .or. (nvfxbp == nvfubp)))  &
            kfbit = MAX( kfbit , ibit1( mflgmfw, nvfxbp + nvfbps(6) ) )
          kflg = ireplace( kflg, FL_GROSS    , 1, kfbit  , 0 )
          kfbit = MAX( ibit1( mflgmfw, nvfxbp + nvfbps(4) )                    &
                     , ibit1( mflgmfw, nvflbp             ) )
          kflg = ireplace( kflg, FL_HEIGHT   , 1, kfbit  , 0 )
        endif
!       if (nvrx  >= 0)  kflg = ireplace( kflg, FL_FG    , 1, mflgqcf, nvrx  )
!       if (nvrx2 >= 0)  kflg = ireplace( kflg, FL_FG_LBC, 1, mflgqcf, nvrx2 )
        if (     (kfmt == 1)                                                   &
            .and.(     ((btest( mflgmfw, nvffbp )).and.(  nvrx == nvru))       &
                  .or. ((btest( mflgmfw, nvfmbp )).and.( (nvrx == nvrt) .or.   &
                                                         (nvrx == nvrq))))) then
          if (istatus(1) == ST_ACTIVE)  istatus (1) =  ST_PASSIVE
          kflg = ibset( kflg, FL_THIN )
        endif
        kflags (1)       = kflg
        ! --> get obs check in fdbk format:
        kcheck (1)       = bit1_pos_nr ( kflags(1) , flags% n                  &
                                       , flags% e(1:flags% n)% value )
        if (kflags(1) == 0)  kcheck (1) = FL_NONE
!       ibuf(:) = (/ istatus(1), kcheck(1), kflags(1) /)

        ! --> state, flags, check: first convert fdbk format into 't_obs' format
        call reverse_code ( istatus , stats )
        call reverse_bits ( kflags  , checks , unknown = CHK_CONSIST )
        call reverse_code ( kcheck  , checks , unknown = CHK_CONSIST )
        obs% body (ibt)% use% state   =  istatus(1)
        obs% body (ibt)% use% flags   =  kflags (1)
        obs% body (ibt)% use% check   =  kcheck (1)

        !   cancel upper-air height obs if pressure was converted from height
        if (      (      obs% varno(ibt)          == VN_HEIGHT )               &
            .and. (     (obs% body (ibt)% lev_typ == VN_HEIGHT)                &
                   .or. (obs% body (ibt)% lev_typ == VN_FLEV )))               &
          call decr_use (obs% body (ibt)% use, STAT_DISMISS, check=CHK_NO_OBS )

!       if (obs% body (ibt)% use% state < STAT_ACTIVE) &
!           obs% body (ibt)% use% check = CHK_DATASET

!       if (ivar >= 12)                                                        &
!         write( 0,'("FLAGS",I3," varno=",I3," state=",2I3," chk=",2I3         &
!                  &," flg=",I8,I11," obs=",F11.2)' )  ivar, obs% varno(ibt)   &
!                   , ibuf(1), obs% body (ibt)% use% state                     &
!                   , ibuf(2), obs% body (ibt)% use% check                     &
!                   , ibuf(3), obs% body (ibt)% use% flags                     &
!                   , obs% body (ibt)% o

      enddo
    endif

  end subroutine check_store_conv

!===============================================================================
  pure &
  integer function ireplace   ( invar, ipos, iboc, irepl, ipsr )
  !-------------------------------------------------------------
  integer ,intent(in) :: invar, ipos, iboc, irepl, ipsr
  !-----------------------------------------------------------------------------
  ! replaces 'iboc' bits starting at bit position 'ipos' of integer word 'inver'
  ! by the 'iboc' bits starting at bit position 'ipsr' from integer word 'irepl'
  !-----------------------------------------------------------------------------
  !
  ireplace = ior( iand( invar, not( ishft( NIBITS(iboc), ipos ) ) )            &
                , ishft( iand( ishft( irepl,-ipsr ), NIBITS(iboc) ), ipos ) )
  !
  end function ireplace

!===============================================================================
  pure &
  integer function ibits   ( invar, ibp, iboc )
  !--------------------------------------------
  integer ,intent(in) :: invar, ibp, iboc
  !---------------------------------------------------------------------------
  ! returns 'iboc' bits starting at bit position 'ibp' of integer word 'inver'
  !---------------------------------------------------------------------------
  !
  ibits = iand( ishft( invar,-ibp ), NIBITS(iboc) )
  !
  end function ibits

!===============================================================================
  pure &
  integer function ibit1   ( invar, ibp )
  !--------------------------------------------
  integer ,intent(in) :: invar, ibp
  !---------------------------------------------------------------------------
  ! returns 1 bit at bit position 'ibp' of integer word 'inver'
  !---------------------------------------------------------------------------
  !
  ibit1 = iand( ishft( invar,-ibp ), 1 )
  !
  end function ibit1

#endif
!===============================================================================
end module mo_cosmo_conv
