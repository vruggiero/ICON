!
!+ 3DVAR/COSMO feedback file tables and table entry values
!
!------------------------------------------------------------------------------
!
MODULE mo_fdbk_tables
!
!------------------------------------------------------------------------------
! Description:
!   3DVAR/COSMO feedback file tables and table entry values,
!   specifying the valid values of the variables in the feedback file.
!   This module is jointly used by the COSMO model and 3DVAR program packages
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Harald Anlauf                    DWD, Christoph Schraff
!    phone: +49 69 8062 4941               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de           email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_2         2008/12/04 Andreas Rhodin
!  tables added for: general cloud group, individual cloud group
! V1_5         2009/05/25 Andreas Rhodin
!  new table surf_char; modified tables ind_cg, gen_cg
!  change varno=RR from rain to precipitation amount
! V1_8         2009/12/09 Andreas Rhodin
!  define: VE_VQC_WEIGHT (variational quality control weight)
! V1_9         2010/04/20 Andreas Rhodin
!  for radar volume observations define: OT_GPSGB, OT_RADAR, VN_RADVEL
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  define codetype OC_WP_JP = 134 (Japanese  wind profiler)
!  define n_vn (number of variable numbers)
!  public :: VN_FLEV
!  new parameter: VN_FLEV (nominal flight level)
!  new flag ST_OBS_ONLY (obs only, no model); new body entry: accuracy (from data provider)
!  define ensemble mean in observation space (VE_ENS_MEAN_OBS)
!  make public: n_ot  ! number of observation types
!  put table "flags" into correct order
!  new module variable "init" : flag if module was already initialised
!  extend table flag_entries
!  new flags FL_NO_BIASCOR, FL_NO_OBS (no bias correction available, no observations in report)
!  modifications required for COSMO to write feedobs files (C.Schraff, M.Lazanowicz)
!  make public: OT_RAD
!  define VE_MEMBER = -6 ! generic value for ensemble member
! V1_15        2011/12/06 Andreas Rhodin
!  define new flag "OPERATOR" (observation operator not applicable)
! V1_17        2011/12/21 Andreas Rhodin
!  new flag values: ST_DISMISS, VN_REFL, VN_RADIANCE
! V1_19        2012-04-16 Andreas Rhodin
!  define VT_LIN_ANA linear operator on analysis (Y^a)
! V1_20        2012-06-18 Andreas Rhodin
!  new flags FG_LBC, VN_CTH, VN_TRH
! V1_22        2013-02-13 Andreas Rhodin
!  new parameter VN_VGUST (vertical gust); fix mnemonic for VN_CTH (CTH)
! V1_23        2013-03-26 Andreas Rhodin
!  new variable ct_nwc (Cloud Type according to NWC SAF)
! V1_27        2013-11-08 Andreas Rhodin
!  new entry VE_BIASCOR for variational bias correction
! V1_29        2014/04/02 Andreas Rhodin
!  define constants for wind/solar power observation operator:
!  OT_POWER, OC_PWIND, OC_PWSOL, VN_PWIND, VN_PWSOL
! V1_31        2014-08-21 Christoph Schraff
!  add OC_ACARS; make OC_WP_JP public
! V1_35        2014-11-07 Andreas Rhodin
!  Add codetype=146 for MODE-S in FF tables;
!  changes for GPSGB (body variable azimuth, level type elevation)
! V1_37        2014-12-23 Andreas Rhodin
!  change OC_ACARS from 244 to 145;
!  homogenize feedback file interface with COSMO V5.1
! V1_42        2015-06-08 Andreas Rhodin
!  add codetype OC_GPSRO = 250, OC_GPSGB = 251
! V1_43        2015-08-19 Andreas Rhodin
!  feedback file: new entry 'veri_operator_flag'
! V1_44        2015-09-30 Christoph Schraff
!  new definitions for COSMO SYNOP: VN_RAD_GL, VN_RAD_DF, VN_RAD_LW
! V1_47        2016-06-06 Andreas Rhodin
!  define 'varno' VN_PRH=17 (pseudo relative humidity); rename code JJ to TMAX
! V1_48        2016-10-06 Roland Potthast
!  fix spelling
! V1_50        2017-01-09 Harald Anlauf
!  add LS_SUPEROBS to table level_sig_entries
! V1_51        2017-02-24 Andreas Rhodin
!  new table entries for COMET: OT_SOIL OC_ASCWS VN_PRH2M VN_Q2M VN_LWC
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO main programs.!!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2007  original source to be used in 3DVAR/LETKF/COSMO
!------------------------------------------------------------------------------
!-------------
! Modules used
!-------------
use mo_t_table, only :      t_table,   &! derived type definition for tables
                       e => t_entry,   &!   and table entries
                            init_table  ! initialise a table
implicit none

!================
! public entities
!================
private
!-------
! Tables
!-------
public :: status            ! report and observation status flags
public :: flags             ! report and observation check  flags
public :: obstype           ! observation types
public :: obstype_oce       ! oceanographic observation types
public :: codetype          ! code types
public :: varno             ! variables
public :: runtype           ! run type
public :: runclass          ! run class
public :: satsens           ! satellite sensors
public :: rsondtype         ! radiosonde type
public :: trackteqn         ! tracking technique
public :: meas_equip        ! measuring equipment used
public :: radiation_corr    ! solar & infrared radiation correction
public :: surftype          ! surface type
public :: soiltype          ! soil type
public :: surf_char         ! model surface characteristics
public :: ct_nwc            ! cloud type according to NWC SAF
!public :: flg_1dvar         ! 1dvar processing flags
public :: flg_cld           ! 1dvar cloud flags
public :: level_sig         ! vertical level significance
public :: phase             ! aircraft phase, radiances fov, GPSRO PCD
public :: rollangle         ! aircraft roll angle
public :: retrtype          ! AMV retrieval type
public :: ensmem            ! generalised ensemble member flag
public :: oflag             ! observation operator processing flag
public :: tovsflag          ! TOVS specific flags

!------------
! subroutines
!------------
public :: init_fdbk_tables  ! initialise the tables (fill in table entries)
public :: clean_fdbk_tables ! deallocate the tables

!----------
! run class
!----------
public :: RC_HAUPT, RC_VOR, RC_ASS, RC_TEST

!-----------------------------------------------
! constants: observation or report status values
!-----------------------------------------------
public :: ST_ACCEPTED, ST_ACTIVE, ST_MERGED, ST_PASSIVE, ST_REJECTED, &
          ST_PAS_REJ, ST_OBS_ONLY, ST_DISMISS

!----------------------------------
! observation or report check flags
!----------------------------------
public :: FL_OBSTYPE, FL_BLACKLIST, FL_SUSP_LOCT, FL_TIME, FL_AREA,      &
          FL_HEIGHT, FL_SURF, FL_CLOUD, FL_PRACTICE, FL_DATASET,         &
          FL_REDUNDANT, FL_FLIGHTTRACK, FL_MERGE, FL_THIN, FL_RULE,      &
          FL_OBS_ERR, FL_GROSS, FL_NO_BIASCOR, FL_FG, FL_FG_LBC, FL_NONE,&
          FL_NO_OBS, FL_OPERATOR

!------------------
! observation types
!------------------
public :: OT_SYNOP,  OT_AIREP, OT_SATOB, OT_DRIBU,  OT_TEMP,  OT_PILOT, &
          OT_SATEM,  OT_PAOB,  OT_SCATT, OT_RAD,    OT_GPSRO, OT_GPSGB, &
          OT_RADAR,  OT_POWER, OT_SOIL,  OT_OBJECT, OT_LIGHTN,          &
          OT_WLIDAR, OT_MWR,   n_ot

!-----------------------
! observation code types
!-----------------------
public :: OC_SRSCD,   OC_ATSCD, OC_AHSCD, OC_ATSHS, OC_AIRCD, OC_CODAR, &
          OC_AMDAR,   OC_STBCD, OC_DRBCD, OC_TESAC, OC_LDTCD, OC_SHTCD, &
          OC_TDROP,   OC_LDPCD, OC_SHPCD, OC_ATOVS, OC_WP_EU, OC_RA_EU, &
          OC_PR_US,   OC_RAVAD, OC_SEVIR, OC_ASCAT, OC_QSCAT, OC_TMPMB, &
          OC_PLTMB,   OC_AIRS,  OC_IASI,  OC_CLPRD, OC_PWIND, OC_PWSOL, &
          OC_WP_JP,   OC_ACARS, OC_MODES, OC_GPSRO, OC_GPS  , OC_GPSGB, &
          OC_ASCWS,   OC_BTEMP, OC_BSHIP, OC_BDROP, OC_TEMPD, OC_METAR, &
          OC_REFLOBJ, OC_CMAN,  OC_TOWER, OC_SWS,   OC_CARS,  OC_ALSD,  &
          OC_STATIST, OC_WL_GB, OC_ICOS,  OC_SCADA, OC_SYNTST,OC_CWS,   &
          OC_GBLIGHT, OC_MWR,   OC_RAMAN, OC_TAMDAR,OC_SSNOW ,OC_SYTEMP,&
          OC_SYTOWR


!------------------------------------
! observation code types (rof format)
!------------------------------------
public :: OC_SRSCD_ROF,OC_CWS_ROF,OC_LDPCD_ROF

!--------------------------------
! oceanographic observation types
!--------------------------------
public :: OT_OCE_PROF, OT_OCE_TRAJ, OT_OCE_MOOR, OT_OCE_SURF, OT_OCE_SSH

!-------------------------------------
! oceanographic observation code types
!-------------------------------------
public :: OC_OCE_OSTIA_SST, OC_OCE_SSH,    OC_OCE_SMOS_L2, OC_OCE_SMOS_L3, &
          OC_OCE_ARGO,      OC_OCE_GLIDER, OC_OCE_ANIBOS,  OC_OCE_MOOREDB, &
          OC_OCE_DRIBU,     OC_OCE_SEAICE

!--------------------------
! for preliminary use only:
!--------------------------
public :: OC_SYNOPDUM, OC_AIREPDUM, OC_TEMPDUM, OC_PILOTDUM


!-----------------
! variable numbers
!-----------------
public :: VN_U, VN_V, VN_Z, VN_DZ, VN_RH, VN_RH2M, VN_PRH, VN_PRH2M, VN_T, VN_TD,&
          VN_T2M, VN_TD2M, VN_TS, VN_TSEA, VN_PTEND, VN_W1, VN_WW, VN_VV, VN_CH, &
          VN_CM, VN_CL, VN_NH, VN_N_L, VN_C, VN_NS, VN_SDEPTH, VN_E, VN_TRTR,    &
          VN_RR, VN_TMAX, VN_GCLG, VN_N, VN_SFALL, VN_PS, VN_DD, VN_FF,          &
          VN_RAWBT, VN_U10M, VN_V10M, VN_Q, VN_Q2M, VN_VT, VN_HEIGHT, VN_BENDANG,&
          VN_IMPPAR, VN_REFR, VN_ZPD, VN_ZWD, VN_SPD, VN_GUST, VN_P,             &
          VN_TMIN, VN_PRED, VN_N_M, VN_N_H, VN_CEIL, VN_W, VN_TURB, VN_NFXME,    &
          VN_ICLG, VN_PWC, VN_NUM, VN_RADVEL, VN_FLEV, VN_ELEV, VN_REFL, VN_HLOS,&
          VN_RADIANCE, VN_RREFL, VN_CTH, VN_TRH, VN_VGUST,                       &
          VN_PWIND, VN_PWSOL, VN_RAD_GL, VN_RAD_DI, VN_RAD_DF, VN_RAD_LW,        &
          VN_LWC, VN_OBJ_LAT, VN_OBJ_LON, VN_OBJ_Z, VN_OBJ_AREA, VN_OBJ_CVIL,    &
          VN_OBJ_NUM, VN_LIGH_FLR, VN_DEPTH, n_vn, VN_HOSAG, VN_NSOILM, VN_SOILM,&
          VN_MIXR, VN_FR_ICE

!------------------
! Greenhouse gases:
!------------------
public :: VN_CO2, VN_CH4, VN_N2O

!-------
! Ocean:
!-------
public :: VN_SWPSAL, VN_SWPT, VN_SWT, VN_SSH, VN_SLA

!--------------------------
! for preliminary use only:
!--------------------------
public :: VN_DUMMY1, VN_DUMMY2, VN_DUMMY3

!--------------------------------
! RTTOV instrument IDs (sensors):
!--------------------------------
public :: SI_HIRS, SI_AMSU_A, SI_AMSU_B, SI_AIRS, SI_IASI, SI_SEVIRI, SI_MWR

!----------------
! radiosonde type
!----------------
public :: RS_GRAW, RS_BASORA, RS_RU_A_MRZ, RS_RU_MET1, RS_80, RS_VIZ_M2,  &
          RS_DC_MODEM, RS_RU_A_BAR, RS_90_DIG12, RS_RU_ARMA, RS_92_DIG12, &
          RS_92_DIG3, RS_92_AUTO, RS_RU_V_MRZ, RS_RU_V_BAR, RS_SA_DAT4G,  &
          RS_MISS

!-------------------
! tracking technique
!-------------------
public :: TT_NOWIND, TT_AUX_OPTIC, TT_AUX_RANGE, TT_LORANC, TT_SATNAV, &
          TT_NOTSPEC, TT_NORMAL, TT_MISS, TT_AIR_PHASE

!----------------------------
! type of measuring equipment
!----------------------------
public :: TME_PRESS, TME_OPTTHEO, TME_RADTHEO, TME_RADAR, TME_VLFOMEGA, &
          TME_LORANC, TME_WINDPROF, TME_SATNAV, TME_RASS, TME_SODAR, TME_MISS

!--------------------------------------
! solar & infrared radiation correction
!--------------------------------------
public :: RC_NO, RC_CS_CI, RC_CS_IN, RC_CS, RC_SO_IN_AUTO, RC_SO_AUTO, &
          RC_SO_IN_CNTRY, RC_SO_CNTRY, RC_MISS

!---------------------------------------
! ct_nwc cloud type according to NWC SAF
!---------------------------------------
public :: CT_NOPR,        CT_LAND_FREE,   CT_SEA_FREE,     CT_LAND_SNOW,  &
          CT_SEA_ICE,     CT_CUM_VLOW,    CT_STRAT_VLOW,   CT_LOW_CUM,    &
          CT_LOW_STRAT,   CT_MED_CUM,     CT_MED_STRAT,    CT_HI_OP_CUM,  &
          CT_HI_OP_STRAT, CT_VHI_OP_CUM,  CT_VHI_OP_STRAT, CT_HI_ST_THIN, &
          CT_HI_ST_MEAN,  CT_HI_ST_THICK, CT_HI_ST_ABOVE,  CT_FRAC,       &
          CT_UNDEF

!-----------------------------------------------------
! surftype (for radiances)
!-----------------------------------------------------
public :: SUR_SEA,      SUR_BLK_PP,   SUR_LAND,     SUR_HIGHLAND, &
          SUR_MISMATCH, SUR_ICE,      SUR_NOICE,    SUR_MISSING,  &
          SUR_MWSURF,   SUR_SNOW,     SUR_NOSNOW

!public :: D1_DATA, D1_MIN, D1_SUR, D1_CLD
!public :: CL_CLEAR, CL_IR_CLOUDY, CL_MW_CLEAR, CL_MW_CLOUDY

!----------
! soiltype
!----------
public :: SOT_ICE,      SOT_ROCK, SOT_SAND, SOT_SANDLOAM, SOT_LOAM, &
          SOT_CLAYLOAM, SOT_CLAY, SOT_PEAT, SOT_SEAWATER, SOT_SEAICE

!-----------------------------------------------------
! surf_char (model surface characteristics)
!-----------------------------------------------------
public :: MS_LAND, MS_SEA, MS_ICE, MS_NO_ICE, MS_SNOW, MS_NO_SNOW

!-------------------
! level significance
!-------------------
public :: LS_SURFACE, LS_STANDARD, LS_TROPO, LS_MAX, LS_SIGN, LS_SUPEROBS

!------------------------------------------------
! phase of aircraft flight and roll angle quality
!------------------------------------------------
public :: PH_UNS, PH_LVR, PH_LVW, PH_ASC, PH_DES, PH_MIS
public :: RA_GOOD, RA_BAD, RA_MIS

!------------------------------------------------------------------
! retrieval type: for AMV satellite derived wind computation method
!------------------------------------------------------------------
public :: RT_IR,  RT_IR1,  RT_IR2,  RT_IR3,  &
          RT_VIS, RT_VIS1, RT_VIS2, RT_VIS3, &
          RT_WV,  RT_WV1,  RT_WV2,  RT_WV3

!-----------------------------------------------------
! specification of verification data (runtype, ensmem)
!-----------------------------------------------------
public :: VT_FORECAST, VT_FIRSTGUESS, VT_PREL_ANA, VT_ANALYSIS, VT_INIT_ANA, &
          VT_LIN_ANA,  VT_FC_SENS
public :: VE_ENS_MEAN, VE_DETERM, VE_ENS_SPREAD, VE_BG_ERROR, VE_TALAGRAND,  &
          VE_VQC_WEIGHT, VE_MEMBER, VE_ENS_MEAN_OBS, VE_BIASCOR

!--------------------------
! observation operator flag
!--------------------------
public :: OF_MISSING, OF_RAD_CLEAR_SKY, OF_RAD_CLOUDY, OF_BT_CLEAR_SKY

!--------------------------
! TOVS flag
!--------------------------
public :: TF_EMIS, TF_EMIS_FASTEM, TF_EMIS_DYNRET, TF_EMIS_GRODY, TF_EMIS_TLSM, &
          TF_EMIS_CNRM, TF_EMIS_ATLS, TF_EMIS_SEA_MOD, TF_EMIS_CAMEL07,         &
          TF_EMIS_CAMELCL, TF_EMIS_UWIR, TF_REFL_BRDF, TF_EMIS_FAILED,          &
          TF_SURF_INFL, TF_SURF_RETR, TF_SURF_TYPE, TF_SURF_MODEL, TF_CLOUD,    &
          TF_AEROSOL, TF_DESERT_DUST, TF_VOLCANIC_ASH, TF_TRACE_GAS, TF_TSKIN,  &
          TF_TSKIN_DYNRET, TF_TSKIN_FAILED, TF_BC_FAILED

!-------------------------------------
! optional TOVS specific variables
!-------------------------------------
public :: n_optv
public :: c_optv

!==============================================================================
! Module variables
!==============================================================================
save
logical :: init = .false.   ! flag if module was already initialised

!==============================================================================
! Run Class
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: RC_HAUPT = 0 ,&
                      RC_VOR   = 1 ,&
                      RC_ASS   = 2 ,&
                      RC_TEST  = 3
!--------------
! table entries
!--------------
type(e) ,target :: runclass_entries (4) = (/  &
e(RC_HAUPT ,'HAUPT','','main forecast cycle'),&
e(RC_VOR   ,'VOR  ','','pre-run            '),&
e(RC_ASS   ,'ASS  ','','assimilation cycle '),&
e(RC_TEST  ,'TEST ','','test (offline)     ')/)

!------
! table
!------
type(t_table) ,pointer :: runclass

!==============================================================================
! Observation or Report Status Flags
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: ST_ACCEPTED =  0 ,&
                      ST_ACTIVE   =  1 ,&
                      ST_MERGED   =  3 ,&
                      ST_PASSIVE  =  5 ,&
                      ST_REJECTED =  7 ,&
                      ST_PAS_REJ  =  9 ,&
                      ST_OBS_ONLY = 11 ,&
                      ST_DISMISS  = 13
!--------------
! table entries
!--------------
type (e) ,target :: status_entries (8) = (/                                     &
e(ST_ACCEPTED ,'ACCEPTED','','active and VQC accepted (used in 3D-Var only)  '),&
e(ST_ACTIVE   ,'ACTIVE  ','','used in the assimilation                       '),&
e(ST_MERGED   ,'MERGED  ','','not used, merged into multilevel report        '),&
e(ST_PASSIVE  ,'PASSIVE ','','not used, only monitored                       '),&
e(ST_REJECTED ,'REJECTED','','not used due to suspicious quality             '),&
e(ST_PAS_REJ  ,'PAS_REJ ','','passive and rejected                           '),&
e(ST_OBS_ONLY ,'OBS_ONLY','','observation only, no model equivalent available'),&
e(ST_DISMISS  ,'DISMISS' ,'','dismiss observation, should not appear in file ')/)

!------
! table
!------
type(t_table) ,pointer :: status

!==============================================================================
! Observation or Report Check Flags
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &
      FL_OBSTYPE     =  0 ,&! passive report type (at obs. location)
      FL_BLACKLIST   =  1 ,&! blacklist (or not whitelist)
      FL_SUSP_LOCT   =  2 ,&! suspicious location or date/time
      FL_TIME        =  3 ,&! time     not in valid range
      FL_AREA        =  4 ,&! location not in valid area
      FL_HEIGHT      =  5 ,&! location not in valid height range
      FL_SURF        =  6 ,&! incorrect surface (land,ice,etc)
      FL_CLOUD       =  7 ,&! cloud check
      FL_PRACTICE    =  8 ,&! bad reporting practice / insufficient data
      FL_DATASET     =  9 ,&! dataset quality flags
      FL_REDUNDANT   = 10 ,&! redundant report
      FL_FLIGHTTRACK = 11 ,&! flight track error flag
      FL_MERGE       = 12 ,&! merged reports (e.g. TEMP ABCD)
      FL_THIN        = 13 ,&! thinning
      FL_RULE        = 14 ,&! complex rule
      FL_OBS_ERR     = 15 ,&! observation error too large
      FL_GROSS       = 16 ,&! gross error flag
      FL_NO_BIASCOR  = 17 ,&! no bias correction available
      FL_FG          = 18 ,&! observation - first guess check
      FL_NO_OBS      = 19 ,&! no observations in report
      FL_OPERATOR    = 20 ,&! observation operator not applicable
      FL_FG_LBC      = 21, &! obs - lateral boundary condition check
      FL_NONE        = 32   ! no flag set
!--------------
! table entries
!--------------
type (e) ,target :: flag_entries (23) = (/                                    &
e(FL_SUSP_LOCT  ,'SUSP_LOCT  ','','suspicious location or date/time     ~ 1'),&
e(FL_TIME       ,'TIME       ','','time     not in valid range          ~ 2'),&
e(FL_AREA       ,'AREA       ','','location not in valid area           ~ 3'),&
e(FL_PRACTICE   ,'PRACTICE   ','','bad reporting practice/insuff. data  ~ 4'),&
e(FL_DATASET    ,'DATASET    ','','dataset quality flags                ~ 5'),&
e(FL_BLACKLIST  ,'BLACKLIST  ','','blacklist (or not whitelist)         ~ 6'),&
e(FL_HEIGHT     ,'HEIGHT     ','','location not in valid height range   ~ 7'),&
e(FL_SURF       ,'SURF       ','','incorrect surface (land,ice,etc)     ~ 8'),&
e(FL_CLOUD      ,'CLOUD      ','','cloud check                          ~ 9'),&
e(FL_GROSS      ,'GROSS      ','','gross error flag                     ~10'),&
e(FL_OBSTYPE    ,'OBSTYPE    ','','passive report type (at obs.location)~11'),&
e(FL_REDUNDANT  ,'REDUNDANT  ','','redundant report                     ~12'),&
e(FL_FLIGHTTRACK,'FLIGHTTRACK','','flight track error flag              ~13'),&
e(FL_MERGE      ,'MERGE      ','','merged reports (e.g. TEMP ABCD)      ~14'),&
e(FL_THIN       ,'THIN       ','','thinning                             ~15'),&
e(FL_RULE       ,'RULE       ','','complex rule                         ~16'),&
e(FL_NO_BIASCOR ,'NO_BIASCOR ','','no bias correction available         ~17'),&
e(FL_OBS_ERR    ,'OBS_ERR    ','','observation error too large          ~18'),&
e(FL_NO_OBS     ,'NO_OBS     ','','no observations in report            ~19'),&
e(FL_FG         ,'FG         ','','observation - first guess check      ~20'),&
e(FL_FG_LBC     ,'FG_LB      ','','obs- lateral boundary condition check~21'),&
e(FL_OPERATOR   ,'OPERATOR   ','','observation operator not applicable  ~22'),&
e(FL_NONE       ,'NONE       ','','no flag set                             ')/)

!------
! table
!------
type(t_table) ,pointer :: flags

!==============================================================================
! Observation Types
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: OT_SYNOP  =  1 ,&
                      OT_AIREP  =  2 ,&
                      OT_SATOB  =  3 ,&
                      OT_DRIBU  =  4 ,&
                      OT_TEMP   =  5 ,&
                      OT_PILOT  =  6 ,&
                      OT_SATEM  =  7 ,&
                      OT_PAOB   =  8 ,&
                      OT_SCATT  =  9 ,&
                      OT_RAD    = 10 ,&
                      OT_GPSRO  = 11 ,&
                      OT_GPSGB  = 12 ,&
                      OT_RADAR  = 13 ,&
                      OT_POWER  = 14 ,&
                      OT_SOIL   = 15 ,&
                      OT_OBJECT = 16 ,&
                      OT_LIGHTN = 17 ,&
                      OT_WLIDAR = 18, &
                      OT_MWR    = 19
integer ,parameter :: n_ot      = 19   ! number of observation types

!-----------------------------------------------
! oceanographic obstypes (starting at 100, tbc)
!-----------------------------------------------
integer ,parameter :: OT_OCE_PROF = 101,  &! in-situ vertical profile, e.g. ARGO + CTD
                      OT_OCE_TRAJ = 102,  &! trajectories, e.g. glider and AniBOS
                      OT_OCE_MOOR = 103,  &! moored
                      OT_OCE_SURF = 104,  &! e.g. OSTIA, SMOS
                      OT_OCE_SSH  = 105    ! sea surface height from satellite
integer ,parameter :: n_ot_oce    = 5      ! number of observation types

!--------------
! table entries
!--------------
type (e) ,target :: obstype_entries (n_ot) = (/                            &
e(OT_SYNOP ,'SYNOP' ,'','SYNOP report,               ~(ECMWF convention)'),&
e(OT_AIREP ,'AIREP' ,'','AIREP report,               ~(ECMWF convention)'),&
e(OT_SATOB ,'SATOB' ,'','SATOB report (AMV),         ~(ECMWF convention)'),&
e(OT_DRIBU ,'DRIBU' ,'','DRIBU report,               ~(ECMWF convention)'),&
e(OT_TEMP  ,'TEMP ' ,'','TEMP  report,               ~(ECMWF convention)'),&
e(OT_PILOT ,'PILOT' ,'','PILOT report,               ~(ECMWF convention)'),&
e(OT_SATEM ,'SATEM' ,'','SATEM report,               ~(ECMWF convention)'),&
e(OT_PAOB  ,'PAOB ' ,'','PAOB  report.               ~(ECMWF convention)'),&
e(OT_SCATT ,'SCATT' ,'','Scatterometer report,       ~(ECMWF convention)'),&
e(OT_RAD   ,'RAD'   ,'','Radiances                   ~(ECMWF convention)'),&
e(OT_GPSRO ,'GPSRO' ,'','GPS Radio occultations,       ~(DWD convention)'),&
e(OT_GPSGB ,'GPSGB' ,'','GPS ground based observations ~(DWD convention)'),&
e(OT_RADAR ,'RADAR' ,'','RADAR (volume data)           ~(DWD convention)'),&
e(OT_POWER ,'POWER' ,'','POWER (win, solar) data       ~(DWD convention)'),&
e(OT_SOIL  ,'SOIL'  ,'','Soil retrieval                ~(DWD convention)'),&
e(OT_OBJECT,'OBJECT','','Objects                       ~(DWD convention)'),&
e(OT_LIGHTN,'LIGHTN','','Lightning                     ~(DWD convention)'),&
e(OT_WLIDAR,'WLIDAR','','Atmospheric Lidar             ~(DWD convention)'),&
e(OT_MWR   ,'MWR'   ,'','Microwave radiometer          ~(DWD convention)')/)

type (e) ,target :: obstype_oce_entries (n_ot_oce) = (/                        &
e(OT_OCE_PROF,'OCE_PROF','','Ocean profile                 ~(DWD convention)'),&
e(OT_OCE_TRAJ,'OCE_TRAJ','','Ocean trajectory              ~(DWD convention)'),&
e(OT_OCE_MOOR,'OCE_MOOR','','Ocean moored report           ~(DWD convention)'),&
e(OT_OCE_SURF,'OCE_SURF','','Ocean surface report          ~(DWD convention)'),&
e(OT_OCE_SSH ,'OCE_SSH' ,'','Satellite sea surface height  ~(DWD convention)')/)

!------
! table
!------
type(t_table) ,pointer :: obstype
type(t_table) ,pointer :: obstype_oce

!==============================================================================
! Observation Code Types
!==============================================================================
! An attempt has been made to align observation code types with ECMWF's ODB,
! when a corresponding code type has been defined there but not at DWD.
! c.f. https://codes.ecmwf.int/odb/codetype/
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: OC_SRSCD   =  11 ,&! synop surface report
                      OC_ATSCD   =  14 ,&! automatic synop surface report
                      OC_SWS     =  15 ,&! Road-Weather-Station data, SWS/SWIS
                      OC_SSNOW   =  16 ,&! Snow observations from SYNOP
                      OC_ALSD    =  17 ,&! Additional Land Surface Data/METAR USA
                      OC_CARS    =  18 ,&! car data
                      OC_CWS     =  19 ,&! Citizen Weather Station (Netatmo)
                      OC_CMAN    =  20 ,&! Coastal Marine Automated Network
                      OC_AHSCD   =  21 ,&! ship synop report
                  !!! OC_ABSCD   =  22 ,&! ship synop abbreviated report
                  !!! OC_SHRED   =  23 ,&! shred reportKp synop report
                      OC_ATSHS   =  24 ,&! automatic ship synop report
                      OC_GPS     = 110 ,&! GPS zenith delay
                      OC_GPSGB   = 251 ,&! GPS slant delay
                      OC_METAR   = 140 ,&! METAR

                      OC_AIRCD   = 141 ,&! AIREP report
                      OC_CODAR   =  41 ,&! CODAR report
                  !!! OC_COLBA   = 241 ,&! COLBA report
                      OC_AMDAR   = 144 ,&! AMDAR report
                      OC_ACARS   = 145 ,&! ACARS report
                      OC_MODES   = 146 ,&! MODE-S report          (ECMWF: 147)
                      OC_TAMDAR  = 147 ,&! TAMDAR or AFIRS report (ECMWF: 148)

                      OC_CLPRD   =  87 ,&! Cloud (height) product
                      OC_STBCD   =  88 ,&! SATOB report
                  !!!               89 ,&! High Resolution VIS wind
                      OC_AMV     =  90 ,&! AMV
                  !!! OC_SST     = 188 ,&! SST report

                      OC_DRBCD   = 165 ,&! DRIBU report
                  !!! OC_BATHY   =  63 ,&! BATHY report
                  !!!              160     ERS as DRIBU
                      OC_TESAC   =  64 ,&! scatterometer

                      OC_LDTCD   =  35 ,&! land TEMP report
                      OC_SHTCD   =  36 ,&! ship TEMP report
                      OC_TMPMB   =  37 ,&! TEMP  mobile
                      OC_PLTMB   =  38 ,&! PILOT mobile
                      OC_TDROP   = 135 ,&! temp-drop report
                      OC_BTEMP   = 109 ,&! LAND TEMP report (high res. BUFR)
                      OC_BSHIP   = 111 ,&! SHIP TEMP report (high res. BUFR)
                      OC_BDROP   = 230 ,&! TEMP DROP report (high res. BUFR)
                      OC_TEMPD   = 231 ,&! TEMP DESCENT report (high res. BUFR)
                  !!! OC_ROCOB   =  39 ,&! ROCOB rep
                  !!! OC_ROCSH   =  40 ,&! ROCOB ship report

                      OC_LDPCD   =  32 ,&! land pilot report
                      OC_SHPCD   =  33 ,&! ship pilot report

                      OC_WP_EU   = 132 ,&! European  wind profiler
                      OC_RA_EU   = 133 ,&! European sodar/rass report
                      OC_WP_JP   = 134 ,&! Japanese  wind profiler
                      OC_PR_US   = 136 ,&! wind/profiler/rass report (USA)
                      OC_RAVAD   = 137 ,&! radar VAD wind profile report
                      OC_TOWER   = 139 ,&! tower profile data
                      OC_ICOS    = 159 ,&! tower profile data (ICOS)
                      OC_SCADA   = 157 ,&! SCADA report
                      OC_WL_GB   = 187 ,&! wind lidar (ground-based)
                      OC_RAMAN   = 190 ,&! Raman lidar
                      OC_DIAL    = 191 ,&! Differential absorption lidar (DIAL)
                      OC_PWIND   = 150 ,&! wind  power data
                      OC_PWSOL   = 151 ,&! solar power data
                  !!!               34     American wind profiler (ECMWF)

                  !!!        8 122 210     Scatterometer

                      OC_ATOVS   = 210 ,&! ATOVS satellite data (radiances)
                  !!! OC_STMCD   =  86 ,&! satem report
                  !!! OC_STOVS   = 186 ,&! high resolution ATOVS satellite data
                  !!! OC_SMSG1   =  71 ,&! MSG_1 satellite retrieval
                  !!! OC_NOA15   = 206 ,&! NOAA15 satellite retrieval
                  !!! OC_NOA16   = 207 ,&! NOAA16 satellite retrieval
                  !!! OC_NOA17   = 208 ,&! NOAA17 satellite retrieval
                  !!! OC_NOA18   = 209 ,&! NOAA18 satellite retrieval
                  !!! OC_GPGFZ   =  96 ,&! GPS report processed by GFZ
                  !!! OC_1DVAR   = 999 ,&! satellite retrieval (1dvar)
                      OC_SEVIR   = 218 ,&! SEVIRI
                      OC_ASCAT   = 123 ,&! ASCAT scatterometer
                      OC_QSCAT   = 122 ,&! QSCAT scatterometer
                      OC_AIRS    = 216 ,&! AIRS
                      OC_IASI    = 217 ,&! IASI
                      OC_GPSRO   = 250 ,&! GPS Radio Occultation
                      OC_ASCWS   = 305 ,&! ASCAT soil moisture retrieval
                      OC_REFLOBJ = 400 ,&! Reflectivity object
                      OC_STATIST = 401 ,&! Statistical object
                      OC_LGT_IMG = 200 ,&! Satellite lightning imager
                      OC_GBLIGHT = 201 ,&! Ground-based lightning network
                      OC_MWR     = 501 ,&! MW radiometer data (ground based)
!--------------------------------------------------------------------------
! 600 <= codetype <= 699 : reserved for oceanographic codetypes,
!--------------------------------------------------------------------------
                      OC_OCE_OSTIA_SST = 600 ,&! OSTIA SST product
                      OC_OCE_SSH       = 601 ,&! SSH from satellite
                      OC_OCE_SMOS_L2   = 602 ,&! Level 2 SMOS product
                      OC_OCE_SMOS_L3   = 603 ,&! Level 3 SMOS product
                      OC_OCE_SEAICE    = 604 ,&! sea ice product

                      OC_OCE_ARGO      = 650 ,&! generic vertical ocean profile esp. ARGO
                      OC_OCE_GLIDER    = 651 ,&! ocean glider
                      OC_OCE_ANIBOS    = 652 ,&! ocean obs from animal
                      OC_OCE_MOOREDB   = 653 ,&! ocean obs from moored buoy
                      OC_OCE_DRIBU     = 654 ,&! ocean obs from DRIBU
!--------------------------------------------------------------------------
! 800 <= codetype <= 899 : reserved for experimental observation codetypes,
!                          intended for non-operational use, data exchange.
!--------------------------------------------------------------------------
                      OC_SYNTST    = 811, &! test synop surface report
                      OC_SYTEMP    = 835 ,&! surface report derived from TEMP
                      OC_SYTOWR    = 839 ,&! surface report derived from tower
                      OC_SRSCD_ROF = 850, &! rof SYNOP surface report
                      OC_CWS_ROF   = 851, &! rof CWS report
                      OC_LDPCD_ROF = 852, &! rof land PILOT report
                      OC_IAGOS_ROF = 853   ! rof AIREP report
!-----------------------------------------------------------------------
! 900 <= codetype <= 989 : reserved for local observation codetypes,
!                          long-term definitions, non-exchangeable data.
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! temporary observation codetypes, not to be used for operational code!
!----------------------------------------------------------------------
integer ,parameter :: OC_SYNOPDUM = 991,&! temporary SYNOP codetype
                      OC_AIREPDUM = 992,&! temporary AIREP codetype
                      OC_TEMPDUM  = 995,&! temporary TEMP  codetype
                      OC_PILOTDUM = 996  ! temporary PILOT codetype

!--------------
! table entries
!--------------
type (e) ,target :: codetype_entries (80) = (/                 &
e(OC_SRSCD   ,'SRSCD'   ,'','synop surface report           '),&
e(OC_ATSCD   ,'ATSCD'   ,'','automatic synop surface report '),&
e(OC_CMAN    ,'CMAN'    ,'','Coastal Marine Automated Network'),&
e(OC_CWS     ,'CWS'     ,'','Citizen Weather Station (Netatmo)'),&
e(OC_SWS     ,'SWIS'    ,'','Road weather station data (SWIS)'),&
e(OC_SSNOW   ,'SSNOW'   ,'','Snow observations from SYNOP   '),&
e(OC_ALSD    ,'ALSD'    ,'','Additional Land Surface Data   '),&
e(OC_CARS    ,'CARS'    ,'','Car data                       '),&
e(OC_AHSCD   ,'AHSCD'   ,'','ship synop report              '),&
e(OC_ATSHS   ,'ATSHS'   ,'','automatic ship synop report    '),&
e(OC_METAR   ,'METAR'   ,'','METAR                          '),&
e(OC_GPS     ,'GPS'     ,'','GPS zenith delay               '),&
e(OC_AIRCD   ,'AIRCD'   ,'','airep report                   '),&
e(OC_CODAR   ,'CODAR'   ,'','codar report                   '),&
e(OC_AMDAR   ,'AMDAR'   ,'','amdar report                   '),&
e(OC_MODES   ,'MODES'   ,'','mode-s report                  '),&
e(OC_ACARS   ,'ACARS'   ,'','acars report                   '),&
e(OC_TAMDAR  ,'TAMDAR'  ,'','TAMDAR or AFIRS report         '),&
e(OC_CLPRD   ,'CLPRD'   ,'','cloud (height) product         '),&
e(OC_STBCD   ,'STBCD'   ,'','satob report                   '),&
e(OC_AMV     ,'AMV'     ,'','AMV                            '),&
e(OC_DRBCD   ,'DRBCD'   ,'','dribu report                   '),&
e(OC_TESAC   ,'TESAC'   ,'','scatterometer                  '),&
e(OC_LDTCD   ,'LDTCD'   ,'','land temp report               '),&
e(OC_SHTCD   ,'SHTCD'   ,'','temp ship report               '),&
e(OC_TDROP   ,'TDROP'   ,'','temp-drop report               '),&
e(OC_TMPMB   ,'TMPMB'   ,'','temp mobile                    '),&
e(OC_BTEMP   ,'BTEMP'   ,'','BUFR TEMP land (high resol.)   '),&
e(OC_BSHIP   ,'BSHIP'   ,'','BUFR TEMP ship (high resol.)   '),&
e(OC_BDROP   ,'BDROP'   ,'','BUFR TEMP drop (high resol.)   '),&
e(OC_TEMPD   ,'TEMPD'   ,'','BUFR TEMP descent (high resol.)'),&
e(OC_LDPCD   ,'LDPCD'   ,'','land pilot report              '),&
e(OC_SHPCD   ,'SHPCD'   ,'','ship pilot report              '),&
e(OC_PLTMB   ,'PLTMB'   ,'','pilot mobile                   '),&
e(OC_ATOVS   ,'ATOVS'   ,'','ATOVS satellite data (1dvar)   '),&
e(OC_WP_EU   ,'WP_EU'   ,'','European  wind profiler        '),&
e(OC_RA_EU   ,'RA_EU'   ,'','European sodar/rass report     '),&
e(OC_WP_JP   ,'WP_JP'   ,'','Japanese  wind profiler        '),&
e(OC_PR_US   ,'PR_US'   ,'','wind/profiler/rass report (USA)'),&
e(OC_RAVAD   ,'RAVAD'   ,'','radar VAD wind profile report  '),&
e(OC_TOWER   ,'TOWER'   ,'','tower profile data             '),&
e(OC_ICOS    ,'ICOS '   ,'','tower profile data (ICOS)      '),&
e(OC_SCADA   ,'SCADA'   ,'','SCADA report                   '),&
e(OC_WL_GB   ,'WL_GB'   ,'','wind lidar (ground-based)      '),&
e(OC_RAMAN   ,'RAMAN'   ,'','Raman lidar                    '),&
e(OC_DIAL    ,'DIAL'    ,'','Differential absorption lidar  '),&
e(OC_PWIND   ,'PWIND'   ,'','wind power data                '),&
e(OC_PWSOL   ,'PWSOL'   ,'','solar power data               '),&
e(OC_SEVIR   ,'SEVIR'   ,'','SEVIRI                         '),&
e(OC_ASCAT   ,'ASCAT'   ,'','ASCAT scatterometer            '),&
e(OC_QSCAT   ,'QSCAT'   ,'','QSCAT scatterometer            '),&
e(OC_AIRS    ,'AIRS'    ,'','AIRS                           '),&
e(OC_IASI    ,'IASI'    ,'','IASI                           '),&
e(OC_GPSRO   ,'GPSRO'   ,'','GPS Radio Occultation          '),&
e(OC_GPSGB   ,'GPSGB'   ,'','GPS slant delay                '),&
e(OC_ASCWS   ,'ASCWS'   ,'','ASCAT soil moisture retrieval  '),&
e(OC_REFLOBJ ,'REFLOBJ' ,'','Reflectivity object            '),&
e(OC_STATIST ,'STATIST' ,'','Statistical object             '),&
e(OC_GBLIGHT ,'GBLIGHT' ,'','Ground-based lightning data    '),&
e(OC_SYNTST  ,'SYNTST ' ,'','test synop surface report      '),&
e(OC_SYTEMP  ,'SYTEMP ' ,'','TEMP-derived surface report    '),&
e(OC_SYTOWR  ,'SYTOWR ' ,'','tower-derived surface report   '),&
e(OC_MWR     ,'MWR'     ,'','MW radiom. data (ground-based) '),&
e(OC_SRSCD_ROF ,'SRSCD_ROF','','SYNOP surface (rof)         '),&
e(OC_CWS_ROF   ,'CWS_ROF  ','','SYNOP CWS (rof)             '),&
e(OC_LDPCD_ROF ,'LDPCD_ROF','','land PILOT (rof)            '),&
!-----------------------------------------------------------------------------
e(OC_OCE_OSTIA_SST ,'OCE_OSTIA_SST','','OSTIA SST product              '),&
e(OC_OCE_SSH       ,'OCE_SSH'      ,'','SSH from satellite             '),&
e(OC_OCE_SMOS_L2   ,'OCE_SMOS_L2'  ,'','Level 2 SMOS product           '),&
e(OC_OCE_SMOS_L3   ,'OCE_SMOS_L3'  ,'','Level 3 SMOS product           '),&
e(OC_OCE_SEAICE    ,'OCE_SEAICE'   ,'','sea ice product                '),&
e(OC_OCE_ARGO      ,'OCE_ARGO'     ,'','generic vertical ocean profile '),&
e(OC_OCE_GLIDER    ,'OCE_GLIDER'   ,'','ocean glider                   '),&
e(OC_OCE_ANIBOS    ,'OCE_ANIBOS'   ,'','ocean obs from animal          '),&
e(OC_OCE_MOOREDB   ,'OCE_MOOREDB'  ,'','ocean obs from moored buoy     '),&
e(OC_OCE_DRIBU     ,'OCE_DRIBU'    ,'','ocean obs from DRIBU           '),&
!-----------------------------------------------------------------------------
e(OC_SYNOPDUM,'SYNOPDUM','','dummy SYNOP                    '),&
e(OC_AIREPDUM,'AIREPDUM','','dummy AIREP                    '),&
e(OC_TEMPDUM ,'TEMPDUM ','','dummy TEMP                     '),&
e(OC_PILOTDUM,'PILOTDUM','','dummy PILOT                    ')/)

!------
! table
!------
type(t_table) ,pointer :: codetype

!==============================================================================
! Variable Numbers
!==============================================================================
! An initial attempt was made to align variable numbers with ECMWF's ODB,
! c.f. https://codes.ecmwf.int/odb/varno/
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: VN_U       =   3 ,&! u-component of wind
                      VN_V       =   4 ,&! v-component of wind
                      VN_Z       =   1 ,&! geopotential
                      VN_DZ      =  57 ,&! thickness
                      VN_PWC     =   9 ,&! precipitable water content kg/m**2
                      VN_TRH     =  28 ,&! transformed relative humidity
                      VN_RH      =  29 ,&! relative humidity
                      VN_RH2M    =  58 ,&! 2 metre relative humidity
                      VN_PRH     =  17 ,&! pseudo relative humidity
                      VN_PRH2M   =  46 ,&! 2 metre pseudo relative humidity
                      VN_T       =   2 ,&! upper air temperature
                      VN_TD      =  59 ,&! upper air dew point
                      VN_T2M     =  39 ,&! 2 metre temperature
                      VN_LWC     =  10 ,&! liquid water content kg/m**2
                      VN_TD2M    =  40 ,&! 2 metre dew point
                      VN_TS      =  11 ,&! surface temperature
                      VN_TSEA    =  12 ,&! sea/water temperature
                      VN_PTEND   =  30 ,&! pressure tendency
                      VN_W1      =  60 ,&! past weather
                      VN_WW      =  61 ,&! present weather
                      VN_VV      =  62 ,&! visibility
                      VN_CH      =  63 ,&! type of high clouds
                      VN_CM      =  64 ,&! type of middle clouds
                      VN_CL      =  65 ,&! type of low clouds
                      VN_NH      =  66   ! cloud base height

integer ,parameter :: VN_N_L     =  67 ,&! low cloud amount
                  !!! VN_HSHS    =  68 ,&! additional cloud group height
                      VN_C       =  69 ,&! additional cloud group type
                      VN_NS      =  70 ,&! additional cloud group amount
                      VN_SDEPTH  =  71 ,&! snow depth
                      VN_E       =  72 ,&! state of ground
                  !!! VN_TGTG    =  73 ,&! ground temperature
                  !!! VN_SPSP1   =  74 ,&! special phenomena 1

                  !!! VN_SPSP2   =  75 ,&! special phenomena 2
                  !!! VN_RS      =  76 ,&! ice code type

                  !!! VN_ESES    =  77 ,&! ice thickness              (1751)
                  !!! VN_IS      =  78 ,&! ice                        (1751)
                      VN_TRTR    =  79 ,&! time period of information (h)
                      VN_RR      =  80 ,&! precipitation amount       (kg/m^2)
                      VN_TMAX    =  81   ! maximum temperature        (K)
                  !!! VN_VS      =  82 ,&! ship speed                 (m/s)
                  !!! VN_DS      =  83 ,&! ship direction             (degree)
                  !!! VN_HWHW    =  84 ,&! wave height                (m)
                  !!! VN_PWPW    =  85 ,&! wave period                (s)
                  !!! VN_DWDW    =  86   ! wave direction             (degree)

integer ,parameter :: VN_GCLG    =  87 ,&! general cloud group         Table
                  !!! VN_RHLC    =  88 ,&! relative humidity from low clouds
                  !!! VN_RHMC    =  89 ,&! relative humidity from middle clouds
                  !!! VN_RHHC    =  90 ,&! relative humidity from high clouds
                      VN_N       =  91 ,&! total cloud amount         (20011)
                      VN_SFALL   =  92 ,&! 6h snow fall               (m)
                      VN_PS      = 110 ,&! surface pressure           (Pa)
                      VN_DD      = 111 ,&! wind direction             (degree)
                      VN_FF      = 112 ,&! wind speed                 (m/s)

                      VN_REFL    = 118 ,&! reflectance                (0..1)
                      VN_RAWBT   = 119 ,&! brightness temperature     (K)
                      VN_RADIANCE= 120 ,&! radiance                   (W/sr/m^3)

                  !!! VN_SATCL   = 121 ,&! cloud amount from satellite(%)
                  !!! VN_SCATSS  = 122 ,&! backscatter                (dB)
                  !!! VN_DU      =   5 ,&! wind shear u-component     (1/s)
                  !!! VN_DV      =   6 ,&! wind shear v-component     (1/s)
                      VN_U10M    =  41 ,&! 10m u-component of wind    (m/s)
                      VN_V10M    =  42 ,&! 10m v-component of wind    (m/s)
                  !!! VN_RHLAY   =  19 ,&! layer relative humidity
                  !!! VN_AUXIL   = 200 ,&! auxiliary variable

                  !!! VN_CLLQW   = 123 ,&! cloud liquid water         (kg/kg)
                  !!! VN_SCATDD  = 124 ,&! ambiguous v-component      (m/s)
                  !!! VN_SCATFF  = 125 ,&! ambiguous u-component      (m/s)
                      VN_Q       =   7 ,&! specific humidity          (kg/kg)
                      VN_Q2M     =  45 ,&! 2 metre specific humidity  (kg/kg)
                  !!! VN_SCATWD  = 126 ,&! ambiguous wind direction   (degree)
                  !!! VN_SCARWS  = 127 ,&! ambiguous wind speed       (m/s)
                      VN_W       =   8 ,&! vertical speed             (m/s)
                      VN_VT      =  56 ,&! virtual temperature        (K)
                  !!! VN_O3LAY   = 130 ,&! ozone layer density        (kg/m^2)
                      VN_DEPTH   = 154 ,&! depth below surface        (m)
                      VN_CTH     = 155 ,&! cloud top height           (m)
                      VN_HEIGHT  = 156 ,&! height                     (m)
                      VN_HOSAG   = 153 ,&! height of sensor above surface (m)
                      VN_FLEV    = 157   ! nominal flight level       (m)
                  !!! VN_XXXX    = 215   ! SSM/I multilevel variable
                  !!! --------- new --------------------
                  !!!             206     ozone                       (Dopson)
                  !!!             160     past weather                 numeric
                  !!!             130     pressure tendency charact.   numeric
integer ,parameter :: VN_TURB    = 244 ,&! degree of turbulence         WMO table 011031
                      VN_NFXME   = 249 ,&! max wind speed (10min mean)(m/s)
                      VN_RREFL   = 192 ,&! radar reflectivity         (Db)
                      VN_RADVEL  = 193 ,&! radial velocity            (m/s)
                      VN_HLOS    = 194 ,&! horizontal line of sight wind (m/s)
                      VN_PDELAY  = 128 ,&! atmospheric path delay     (m)
                      VN_BENDANG = 162 ,&! radio occ. bending angle   (Rad)
                      VN_ICLG    =  95 ,&! individual cloud layer group Table
                      VN_N_M     =  93 ,&! middle cloud amount         WMO table 020011
                      VN_N_H     =  94 ,&! high   cloud amount         WMO table 020011
                      VN_CEIL    =  96 ,&! cloud ceiling a. orography (m)
                  !!! ---------------------------------------------
                      VN_ELEV    = 158 ,&! elevation                  (degree)
                      VN_NSOILM  = 179 ,&! normalized soil moisture   (0..1)
                      VN_SOILM   = 180 ,&! volumetric soil moisture   (m^3/m^3)
                      VN_MIXR    = 226 ,&! water vapour mixing ratio  (kg/kg)
                      VN_PWIND   = 230 ,&! wind  power data
                      VN_PWSOL   = 231 ,&! solar power data
                      VN_RAD_DI  = 236 ,&! acc.direct solar radiation (J/m2)
                      VN_RAD_GL  = 237 ,&! acc.global solar radiation (J/m2)
                      VN_RAD_DF  = 238 ,&! acc.diffuse solar radiat.  (J/m2)
                      VN_RAD_LW  = 239 ,&! acc.long-wave radiation    (J/m2)
                      VN_IMPPAR  = 252 ,&! impact parameter           (m)
                      VN_REFR    = 248 ,&! refractivity
                      VN_ZPD     = 245 ,&! zenith path delay
                      VN_ZWD     = 246 ,&! zenith wet delay
                      VN_SPD     = 247 ,&! slant path delay
                      VN_VGUST   = 240 ,&! vertical gust (aircrafts)  (m/s)
                      VN_GUST    = 242 ,&! wind gust                  (m/s)
                      VN_P       = 251 ,&! pressure                   (Pa)
                      VN_TMIN    = 243 ,&! minimum temperature        (K)
                      VN_PRED    = 241 ,&! reduced pressure           (Pa)
                      VN_FR_ICE  = 253 ,&! sea ice area fraction      (0..1)
                      VN_NUM     =   0   ! ordinal (channel) number   (  )

!------------------
! Greenhouse gases:
!------------------
integer ,parameter :: VN_CO2     = 186 ,&! carbon dioxide
                      VN_CH4     = 188 ,&! methane
                      VN_N2O     = 189   ! nitrous oxide ("laughing gas")

!-------
! Ocean:
!-------
integer ,parameter :: VN_SWPSAL  = 224 ,&! ocean practical salinity   (PSU)
                      VN_SWPT    = 254 ,&! sea water pot. temperature (K)
                      VN_SWT     = 255 ,&! sea water temperature      (K)
                      VN_SSH     = 273 ,&! sea surface height         (m)
                      VN_SLA     = 287   ! sea level anomaly          (m)

!---------------
! "Objects" etc.
!---------------
integer ,parameter :: VN_OBJ_LAT = 500 ,&! centroid latitude          (degree)
                      VN_OBJ_LON = 501 ,&! centroid longitude         (degree)
                      VN_OBJ_Z   = 502 ,&! centroid height            (m)
                      VN_OBJ_AREA= 503 ,&! area of projected polygon  (m**2)
                      VN_OBJ_CVIL= 504 ,&! cell based VIL             (kg/m^2)
                      VN_OBJ_NUM = 505 ,&! number of objects          (1)
                      VN_LIGH_FLR= 600   ! lightning flash rate       (/m**2/s)
                  !!! ---------------------------------------------

!-----------------------------------------------------------------
! temporary variable numbers, not to be used for operational code!
!-----------------------------------------------------------------
integer ,parameter :: VN_DUMMY1  = 991, &! dummy variable numbers..
                      VN_DUMMY2  = 992, &! ..
                      VN_DUMMY3  = 993   ! .. for testing only

integer, parameter :: n_vn       = 104   ! number of variable numbers

!--------------
! table entries
!--------------
type (e) ,target :: varno_entries (n_vn) = (/                                 &
e(VN_NUM     ,'NUM'     ,''          ,'ordinal (channel) number         ~L '),&
e(VN_U       ,'U'       ,'m/s'       ,'u-component of wind                 '),&
e(VN_V       ,'V'       ,'m/s'       ,'v-component of wind                 '),&
e(VN_W       ,'W'       ,'m/s'       ,'vertical velocity                   '),&
e(VN_Z       ,'Z'       ,'(m/s)**2'  ,'geopotential               ~(also L)'),&
e(VN_DZ      ,'DZ'      ,'(m/s)**2'  ,'thickness                           '),&
e(VN_PWC     ,'PWC'     ,'kg/m**2'   ,'precipitable water content          '),&
e(VN_TRH     ,'TRH'     ,'0..1'      ,'transformed relative humidity       '),&
e(VN_RH      ,'RH'      ,'0..1'      ,'relative humidity                   '),&
e(VN_RH2M    ,'RH2M'    ,'0..1'      ,'2 metre relative humidity           '),&
e(VN_PRH     ,'PRH'     ,'0..1'      ,'pseudo relative humidity            '),&
e(VN_PRH2M   ,'PRH2M'   ,'0..1'      ,'2 metre pseudo relative humidity    '),&
e(VN_LWC     ,'LWC'     ,'kg/m**2'   ,'liquid water content                '),&
e(VN_T       ,'T'       ,'K'         ,'upper air temperature               '),&
e(VN_TD      ,'TD'      ,'K'         ,'upper air dew point                 '),&
e(VN_T2M     ,'T2M'     ,'K'         ,'2 metre temperature                 '),&
e(VN_TD2M    ,'TD2M'    ,'K'         ,'2 metre dew point                   '),&
e(VN_TS      ,'TS'      ,'K'         ,'surface temperature                 '),&
e(VN_TSEA    ,'TSEA'    ,'K'         ,'sea/water temperature               '),&
e(VN_PTEND   ,'PTEND'   ,'Pa/3h'     ,'pressure tendency                   '),&
e(VN_W1      ,'W1'      ,'WMO 020004','past weather                        '),&
e(VN_WW      ,'WW'      ,'WMO 020003','present weather                     '),&
e(VN_VV      ,'VV'      ,'m'         ,'visibility                          '),&
e(VN_CH      ,'CH'      ,'WMO 020012','type of high clouds                 '),&
e(VN_CM      ,'CM'      ,'WMO 020012','type of middle clouds               '),&
e(VN_CL      ,'CL'      ,'WMO 020012','type of low clouds                  '),&
e(VN_NH      ,'NH'      ,'m'         ,'cloud base height (cover >= 1/8)    '),&
e(VN_CEIL    ,'CEIL'    ,'m'         ,'cloud ceiling AGL (cover >= 5/8)    '),&
e(VN_N_L     ,'N_L'     ,'WMO 020011','low cloud amount                    '),&
e(VN_N_M     ,'N_M'     ,'WMO 020011','medium cloud amount                 '),&
e(VN_N_H     ,'N_H'     ,'WMO 020011','high cloud amount                   '),&
e(VN_C       ,'C'       ,'WMO  500'  ,'additional cloud group type      ~L '),&
e(VN_NS      ,'NS'      ,'WMO 2700'  ,'additional cloud group amount       '),&
e(VN_SDEPTH  ,'SDEPTH'  ,'m'         ,'snow depth                          '),&
e(VN_E       ,'E'       ,'WMO 020062','state of ground                     '),&
e(VN_TRTR    ,'TRTR'    ,'h'         ,'time period of information       ~L '),&
e(VN_RR      ,'RR'      ,'kg/m**2'   ,'precipitation amount                '),&
e(VN_TMAX    ,'TMAX'    ,'K'         ,'maximum temperature                 '),&
e(VN_GCLG    ,'GCLG'    ,'Table 6'   ,'general cloud group                 '),&
e(VN_N       ,'N'       ,'WMO 020011','total cloud amount                  '),&
e(VN_SFALL   ,'SFALL'   ,'m'         ,'6h snow fall                        '),&
e(VN_PS      ,'PS'      ,'Pa'        ,'surface (station) pressure          '),&
e(VN_DD      ,'DD'      ,'degree'    ,'wind direction                      '),&
e(VN_FF      ,'FF'      ,'m/s'       ,'wind speed                          '),&
e(VN_REFL    ,'REFL'    ,'0..1'      ,'reflectance                         '),&
e(VN_RAWBT   ,'RAWBT'   ,'K'         ,'brightness temperature              '),&
e(VN_RADIANCE,'RADIANCE','W/sr/m**3' ,'radiance                            '),&
e(VN_U10M    ,'U10M'    ,'m/s'       ,'10m u-component of wind             '),&
e(VN_V10M    ,'V10M'    ,'m/s'       ,'10m v-component of wind             '),&
e(VN_Q       ,'Q'       ,'kg/kg'     ,'specific humidity                   '),&
e(VN_Q2M     ,'Q2M'     ,'kg/kg'     ,'2 metre specific humidity           '),&
e(VN_VT      ,'VT'      ,'K'         ,'virtual temperature                 '),&
e(VN_DEPTH   ,'DEPTH'   ,'m'         ,'depth below surface                 '),&
e(VN_CTH     ,'CTH'     ,'m'         ,'cloud top height                    '),&
e(VN_HEIGHT  ,'HEIGHT'  ,'m'         ,'height                           ~L '),&
e(VN_HOSAG   ,'HOSAG'   ,'m'         ,'height of sensor above ground    ~L '),&
e(VN_FLEV    ,'FLEV'    ,'m'         ,'nominal flight level             ~L '),&
e(VN_ELEV    ,'ELEV'    ,'degree'    ,'elevation                        ~L '),&
e(VN_NSOILM,  'NSOILM'  ,'0..1'      ,'normalized soil moisture            '),&
e(VN_SOILM,   'SOILM'   ,'m**3/m**3' ,'volumetric soil moisture            '),&
e(VN_MIXR    ,'MIXR'    ,'kg/kg'     ,'water vapour mixing ratio           '),&
e(VN_PWIND   ,'PWIND'   ,'W'         ,'wind power data                     '),&
e(VN_PWSOL   ,'PWSOL'   ,'W'         ,'solar power data                    '),&
e(VN_RREFL   ,'RREFL'   ,'Db'        ,'radar reflectivity                  '),&
e(VN_RADVEL  ,'RADVEL'  ,'m/s'       ,'radial velocity                     '),&
e(VN_HLOS    ,'HLOS'    ,'m/s'       ,'horizontal line of sight wind       '),&
e(VN_PDELAY  ,'PDELAY'  ,'m'         ,'atmospheric path delay              '),&
e(VN_BENDANG ,'BENDANG' ,'rad'       ,'bending angle                       '),&
e(VN_IMPPAR  ,'IMPPAR'  ,'m'         ,'impact parameter                 ~L '),&
e(VN_REFR    ,'REFR'    ,''          ,'refractivity                        '),&
e(VN_ZPD     ,'ZPD'     ,''          ,'zenith path delay                   '),&
e(VN_ZWD     ,'ZWD'     ,''          ,'zenith wet delay                    '),&
e(VN_SPD     ,'SPD'     ,''          ,'slant path delay                    '),&
e(VN_VGUST   ,'VGUST'   ,'m/s'       ,'vertical gust (aircrafts)           '),&
e(VN_GUST    ,'GUST'    ,'m/s'       ,'wind gust                           '),&
e(VN_P       ,'P'       ,'Pa'        ,'pressure                         ~L '),&
e(VN_TMIN    ,'TMIN'    ,'K'         ,'minimum temperature                 '),&
e(VN_RAD_DI  ,'RAD_DI'  ,'J/m**2'    ,'direct solar radiation              '),&
e(VN_RAD_GL  ,'RAD_GL'  ,'J/m**2'    ,'global solar radiation              '),&
e(VN_RAD_DF  ,'RAD_DF'  ,'J/m**2'    ,'diffuse solar radiation             '),&
e(VN_RAD_LW  ,'RAD_LW'  ,'J/m**2'    ,'long-wave (downward) radiation      '),&
e(VN_PRED    ,'PRED'    ,'Pa'        ,'reduced pressure                    '),&
e(VN_TURB    ,'TURB'    ,'WMO 011031','degree of turbulence                '),&
e(VN_NFXME   ,'NFXME'   ,'m/s'       ,'max wind speed (10min mean)         '),&
e(VN_ICLG    ,'ICLG'    ,'Table 7'   ,'individual cloud layer group        '),&
e(VN_CO2     ,'CO2'     ,'kg/kg'     ,'carbon dioxide                      '),&
e(VN_CH4     ,'CH4'     ,'kg/kg'     ,'methane                             '),&
e(VN_N2O     ,'N2O'     ,'kg/kg'     ,'nitrous oxide                       '),&
e(VN_FR_ICE  ,'FR_ICE'  ,'0..1'      ,'sea ice area fraction               '),&
e(VN_SWPSAL  ,'SALINITY','PSU'       ,'ocean practical salinity            '),&
e(VN_SWPT    ,'SWPT'    ,'K'         ,'sea water potential temperature     '),&
e(VN_SWT     ,'SWT'     ,'K'         ,'sea water temperature               '),&
e(VN_SSH     ,'SSH'     ,'m'         ,'sea surface height                  '),&
e(VN_SLA     ,'SLA'     ,'m'         ,'sea level anomaly                   '),&
e(VN_OBJ_LAT ,'OBJ_LAT' ,'degree'    ,'centroid latitude                   '),&
e(VN_OBJ_LON ,'OBJ_LON' ,'degree'    ,'centroid longitude                  '),&
e(VN_OBJ_Z   ,'OBJ_Z'   ,'m'         ,'centroid height                     '),&
e(VN_OBJ_AREA,'OBJ_AREA','m**2'      ,'area of projected polygon           '),&
e(VN_OBJ_CVIL,'OBJ_CVIL','kg/m**2'   ,'cell-based vertical integrated liq. '),&
e(VN_OBJ_NUM ,'OBJ_NUM' ,''          ,'number of objects                   '),&
e(VN_LIGH_FLR,'LIGH_FLR','/m**2/s'   ,'flash rate                          '),&
e(VN_DUMMY1,  'DUMMY1  ',''          ,'dummy varno for preliminary use     '),&
e(VN_DUMMY2,  'DUMMY2  ',''          ,'dummy varno for preliminary use     '),&
e(VN_DUMMY3,  'DUMMY3  ',''          ,'dummy varno for preliminary use     ')/)

!------
! table
!------
type(t_table) ,pointer :: varno

!==============================================================================
! General cloud group
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &!
!       GC_SGBP =  0 ,&! bit position for vert.signif.      (WMO 008002)
        GC_CLBP =  0 ,&! bit position for cloud amount      (WMO 020011)
        GC_LCBP =  4 ,&! bit position for low cloud type    (WMO 020012)
        GC_MCBP = 10 ,&! bit position for middle cloud type (WMO 020012)
        GC_HCBP = 16 ,&! bit position for high cloud type   (WMO 020012)
!       GC_SGOC =  6 ,&! no.bits used for vert.signif.      (WMO 008002)
        GC_CLOC =  4 ,&! no.bits used for cloud amount      (WMO 020011)
        GC_LCOC =  6 ,&! no.bits used for low cloud type    (WMO 020012)
        GC_MCOC =  6 ,&! no.bits used for middle cloud type (WMO 020012)
        GC_HCOC =  6   ! no.bits used for high cloud type   (WMO 020012)

!--------------
! table entries
!--------------
type (e) ,target :: gen_cg_entries (8) = (/                              &
!e(GC_SGBP,'SGBP'   ,'WMO 008002' ,'bit position for vert.signif.'     ),&
 e(GC_CLBP,'CLBP'   ,'WMO 020011' ,'bit position for cloud amount'     ),&
 e(GC_LCBP,'LCBP'   ,'WMO 020012' ,'bit position for low cloud type'   ),&
 e(GC_MCBP,'MCBP'   ,'WMO 020012' ,'bit position for middle cloud type'),&
 e(GC_HCBP,'HCBP'   ,'WMO 020012' ,'bit position for high type'        ),&
!e(GC_SGOC,'SGOC'   ,'WMO 008002' ,'no.bits used for vert.signif.'     ),&
 e(GC_CLOC,'CLOC'   ,'WMO 020011' ,'no.bits used for cloud amount'     ),&
 e(GC_LCOC,'LCOC'   ,'WMO 020012' ,'no.bits used for low cloud type'   ),&
 e(GC_MCOC,'MCOC'   ,'WMO 020012' ,'no.bits used for middle cloud type'),&
 e(GC_HCOC,'HCOC'   ,'WMO 020012' ,'no.bits used for high type'        )/)

!------
! table
!------
type(t_table) ,pointer :: gen_cg

!==============================================================================
! Individual cloud group
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: &!
!       IC_SGBP =  0 ,&! bit position for vert.signif. (WMO 008002)
        IC_CLBP =  0 ,&! bit position for cloud amount (WMO 020011)
        IC_CTBP =  4 ,&! bit position for cloud type   (WMO 020012)
        IC_BSBP = 10 ,&! bit position for cloud base height     (m)
!       IC_SGOC =  6 ,&! no.bits used for vert.signif. (WMO 008002)
        IC_CLOC =  4 ,&! no.bits used for cloud amount (WMO 020011)
        IC_CTOC =  6 ,&! no.bits used for cloud type   (WMO 020012)
        IC_BSOC = 14   ! no.bits used for cloud base height     (m)
!--------------
! table entries
!--------------
type (e) ,target :: ind_cg_entries (6) = (/                              &
!e(IC_SGBP,'SGBP'   ,'WMO 008002' ,'bit position for vert.signif.'     ),&
 e(IC_CLBP,'CLBP'   ,'WMO 020011' ,'bit position for cloud amount'     ),&
 e(IC_CTBP,'CTBP'   ,'WMO 020012' ,'bit position for cloud type'       ),&
 e(IC_BSBP,'BSBP'   ,'m'          ,'bit position for cloud base height'),&
!e(IC_SGOC,'SGOC'   ,'WMO 008002' ,'no.bits used for vert.signif.'     ),&
 e(IC_CLOC,'CLOC'   ,'WMO 020011' ,'no.bits used for cloud amount'     ),&
 e(IC_CTOC,'CTOC'   ,'WMO 020012' ,'no.bits used for cloud type'       ),&
 e(IC_BSOC,'BSOC'   ,'m'          ,'no.bits used for cloud base height')/)

!------
! table
!------
type(t_table) ,pointer :: ind_cg

!==============================================================================
! Satellite Instruments (Sensors)
! (RTTOV Instrument ID codes, see inst_id_* in rttov*/rttov_const.f90)
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: SI_HIRS      =  0 ,&! RTTOV8 channels numbers : 1-  19
                  !!! SI_MSU       =  1 ,&!                           1-   4
                  !!! SI_SSU       =  2 ,&!                           1-   3
                      SI_AMSU_A    =  3 ,&!                           1-  15
                      SI_AMSU_B    =  4 ,&!                           1-   5
                  !!! SI_AVHRR     =  5 ,&!                           1-   3
                  !!! SI_SSMI      =  6 ,&!                           1-   4
                  !!! SI_VTPR1     =  7 ,&!                           1-   8
                  !!! SI_VTPR2     =  8 ,&!                           1-   8
                  !!! SI_TMI       =  9 ,&!                           1-   9
                  !!! SI_SSMIS     = 10 ,&!                           1-  21
                      SI_AIRS      = 11 ,&!                           1-2378
                  !!! SI_HSB       = 12 ,&!                           1-   4
                  !!! SI_MODIS     = 13 ,&!                           1-  17
                  !!! SI_ATSR      = 14 ,&!                           1-   3
                  !!! SI_MHS       = 15 ,&!                           1-   5
                      SI_IASI      = 16 ,&!                           1-8461
                  !!! SI_AMSR      = 17 ,&!                           1-   7
                  !!! SI_MVIRI     = 20 ,&!                           1-   2
                      SI_SEVIRI    = 21 ,&!                           1-   8
                  !!! SI_GOES_IM   = 22 ,&! GOES Imager               1-   4
                  !!! SI_GOES_SND  = 23 ,&! GOES Sounder              1-  18
                  !!! SI_GMS_IM    = 24 ,&! GMS/MTSAT imager          1-   4
                  !!! SI_FY2_VISSR = 25 ,&!                           1-   2
                  !!! SI_FY1_VISSR = 26 ,&!                           1-   3
                  !!! SI_CRIS      = 27 ,&!
                  !!! SI_CMSS      = 28 ,&!
                  !!! SI_VIIRS     = 29 ,&!
                  !!! SI_WINDSAT   = 30   !                           1-   5
                      SI_MWR       = 31   ! MWR data RTTOV_GB         1-  22
!--------------
! table entries
!--------------
type (e) ,target :: satsens_entries (7) = (/&
e(SI_HIRS   ,'HIRS'   ,'' ,''),&
e(SI_AMSU_A ,'AMSU_A' ,'' ,''),&
e(SI_AMSU_B ,'AMSU_B' ,'' ,''),&
e(SI_AIRS   ,'AIRS'   ,'' ,''),&
e(SI_IASI   ,'IASI'   ,'' ,''),&
e(SI_SEVIRI ,'SEVIRI' ,'' ,''),&
e(SI_MWR    ,'MWR'    ,'' ,'')/)

!------
! table
!------
type(t_table) ,pointer :: satsens


!==============================================================================
! radiosonde type
!==============================================================================
integer ,parameter :: RS_GRAW     =  17 ,&! Graw   (D)
                      RS_BASORA   =  26 ,&! Basora (CH)
                      RS_RU_A_MRZ =  27 ,&! AVK-MRZ (Russia)
                      RS_RU_MET1  =  28 ,&! Meteorit Marz2-1 (Russia)
                      RS_80       =  37 ,&! Vaisala RS 80
                      RS_VIZ_M2   =  49 ,&! VIZ MARK II (USA)
                      RS_DC_MODEM =  57 ,&! M2K2-DC Modem (France)
                      RS_RU_A_BAR =  58 ,&! AVK-BAR (Russia)
                      RS_90_DIG12 =  71 ,&! Vaisala RS 90 Digicora I,II or Marvin
                      RS_RU_ARMA  =  75 ,&! AVK-MRZ-ARMA (Russia)
                      RS_92_DIG12 =  79 ,&! Vaisala RS 92 Digicora I,II or Marvin
                      RS_92_DIG3  =  80 ,&! Vaisala RS 92 Digicora III
                      RS_92_AUTO  =  81 ,&! Vaisala RS 92 Autosonde
                      RS_RU_V_MRZ =  88 ,&! MARL-A or Vektor-M-MRZ (Russia)
                      RS_RU_V_BAR =  89 ,&! MARL-A or Vektor-M-BAR (Russia)
                      RS_SA_DAT4G =  99 ,&! BAT-4G (South Africa)
                      RS_MISS     = 255   ! missing value

!--------------
! table entries
!--------------
type (e) ,target :: rsondtype_entries (17) = (/                          &
e(RS_GRAW    ,'GRAW'       ,'' ,'Graw   (D)')                           ,&
e(RS_BASORA  ,'BASORA'     ,'' ,'Basora (CH)')                          ,&
e(RS_RU_A_MRZ,'RS_RU_A_MRZ','' ,'AVK-MRZ (Russia)')                     ,&
e(RS_RU_MET1 ,'RS_RU_MET1' ,'' ,'Meteorit Marz2-1 (Russia)')            ,&
e(RS_80      ,'RS_80'      ,'' ,'Vaisala RS 80')                        ,&
e(RS_VIZ_M2  ,'RS_VIZ_M2'  ,'' ,'VIZ MARK II (USA)')                    ,&
e(RS_DC_MODEM,'RS_DC_MODEM','' ,'M2K2-DC Modem (France)')               ,&
e(RS_RU_A_BAR,'RS_RU_A_BAR','' ,'AVK-BAR (Russia)')                     ,&
e(RS_90_DIG12,'RS_90_DIG12','' ,'Vaisala RS 90 Digicora I,II or Marvin'),&
e(RS_RU_ARMA ,'RS_RU_ARMA' ,'' ,'AVK-MRZ-ARMA (Russia)')                ,&
e(RS_92_DIG12,'RS_92_DIG12','' ,'Vaisala RS 92 Digicora I,II or Marvin'),&
e(RS_92_DIG3 ,'RS_92_DIG3' ,'' ,'Vaisala RS 92 Digicora III')           ,&
e(RS_92_AUTO ,'RS_92_AUTO' ,'' ,'Vaisala RS 92 Autosonde')              ,&
e(RS_RU_V_MRZ,'RS_RU_V_MRZ','' ,'MARL-A or Vektor-M-MRZ (Russia)')      ,&
e(RS_RU_V_BAR,'RS_RU_V_BAR','' ,'MARL-A or Vektor-M-BAR (Russia)')      ,&
e(RS_SA_DAT4G,'RS_SA_DAT4G','' ,'BAT-4G (South Africa)')                ,&
e(RS_MISS    ,'MISS'       ,'' ,'missing value')/)

!------
! table
!------
type(t_table) ,pointer :: rsondtype

!==============================================================================
! tracking technique
!==============================================================================
integer ,parameter :: TT_NOWIND    =   0 ,&! no windfinding
                      TT_AUX_OPTIC =   2 ,&! automatic with aux. optical direction finding
                      TT_AUX_RANGE =   3 ,&! automatic with auxiliary ranging
                      TT_LORANC    =   6 ,&! automatic cross chain Loran-C
                      TT_SATNAV    =   8 ,&! automatic satellite navigation
                      TT_NOTSPEC   =  19 ,&! tracking technique not specified
                      TT_NORMAL    =  70 ,&! all systems in normal operation
                      TT_MISS      = 127 ,&! missing value
                      TT_AIR_PHASE =  18   ! aircraft obs: flight phase
                                           !   determined in data assimilation
!--------------
! table entries
!--------------
type (e) ,target :: trackteqn_entries (9) = (/                                  &
e(TT_NOWIND   ,'NOWIND'   ,'' ,'no windfinding')                               ,&
e(TT_AUX_OPTIC,'AUX_OPTIC','' ,'automatic with aux. optical direction finding'),&
e(TT_AUX_RANGE,'AUX_RANGE','' ,'automatic with auxiliary ranging')             ,&
e(TT_LORANC   ,'LORANC'   ,'' ,'automatic cross chain Loran-C')                ,&
e(TT_SATNAV   ,'SATNAV'   ,'' ,'automatic satellite navigation')               ,&
e(TT_NOTSPEC  ,'NOTSPEC'  ,'' ,'tracking technique not specified')             ,&
e(TT_NORMAL   ,'NORMAL'   ,'' ,'all systems in normal operation')              ,&
e(TT_MISS     ,'MISS'     ,'' ,'missing value')                                ,&
e(TT_AIR_PHASE,'AIR_PHASE','' ,'aircraft obs: flight phase from data assimil.')/)

!------
! table
!------
type(t_table) ,pointer :: trackteqn

!==============================================================================
! type of measuring equipment used
!==============================================================================
integer ,parameter :: &
  TME_PRESS    = 0 ,&! pressure instrument associated with wind measuring equipment
  TME_OPTTHEO  = 1 ,&! optical theodolite
  TME_RADTHEO  = 2 ,&! radio theodolite
  TME_RADAR    = 3 ,&! radar
  TME_VLFOMEGA = 4 ,&! VLF-Omega
  TME_LORANC   = 5 ,&! Loran-C
  TME_WINDPROF = 6 ,&! wind profiler
  TME_SATNAV   = 7 ,&! satellite navigation
  TME_RASS     = 8 ,&! radio acoustic sounding system (RASS)
  TME_SODAR    = 9 ,&! SODAR
  TME_MISS     = 15  ! missing

!--------------
! table entries
!--------------
type (e) ,target :: meas_equip_entries (11) = (/     &
e(TME_PRESS   ,'PRESS'   ,'','pressure instrument associated with wind measuring equipment'),&
e(TME_OPTTHEO ,'OPTTHEO' ,'','optical theodolite'  ),&
e(TME_RADTHEO ,'RADTHEO' ,'','radio theodolite'    ),&
e(TME_RADAR   ,'RADAR'   ,'','radar'               ),&
e(TME_VLFOMEGA,'VLFOMEGA','','VLF-Omega'           ),&
e(TME_WINDPROF,'WINDPROF','','wind profiler'       ),&
e(TME_LORANC  ,'LORANC'  ,'','Loran-C'             ),&
e(TME_SATNAV  ,'SATNAV'  ,'','satellite navigation'),&
e(TME_RASS    ,'RASS'    ,'','radio acoustic sounding system (RASS)'),&
e(TME_SODAR   ,'SODAR'   ,'','SODAR'               ),&
e(TME_MISS    ,'MISS'    ,'','missing'             )/)

!------
! table
!------
type(t_table) ,pointer :: meas_equip

!==============================================================================
! solar & infrared radiation correction (NSR, 002013)
!==============================================================================
integer ,parameter ::  &
  RC_NO          =  0 ,&! no correction
  RC_CS_CI       =  1 ,&! CIMO solar + CIMO infrared corrected
  RC_CS_IN       =  2 ,&! CIMO solar + infrared corrected
  RC_CS          =  3 ,&! CIMO solar corrected only
  RC_SO_IN_AUTO  =  4 ,&! solar + infrared corr., automatic. by rsond. system
  RC_SO_AUTO     =  5 ,&! solar corrected automatically by radiosonde system
  RC_SO_IN_CNTRY =  6 ,&! solar + infrared corr. as specified by country
  RC_SO_CNTRY    =  7 ,&! solar corrected by country
  RC_MISS        = 15   ! missing

!--------------
! table entries
!--------------
type (e) ,target :: radiation_corr_entries (9) = (/&
e(RC_NO         ,'NO'         ,'','no correction'),&
e(RC_CS_CI      ,'CS_CI'      ,'','CIMO solar + CIMO infrared corrected'),&
e(RC_CS_IN      ,'CS_IN'      ,'','CIMO solar + infrared corrected'),&
e(RC_CS         ,'CS'         ,'','CIMO solar corrected only'),&
e(RC_SO_IN_AUTO ,'SO_IN_AUTO' ,'','solar + infrared corr., automatic. by rsond. system'),&
e(RC_SO_AUTO    ,'SO_AUTO'    ,'','solar corrected automatically by radiosonde system'),&
e(RC_SO_IN_CNTRY,'SO_IN_CNTRY','','solar + infrared corr. as specified by country'),&
e(RC_SO_CNTRY   ,'SO_CNTRY'   ,'','solar corrected by country'),&
e(RC_MISS       ,'MISS'       ,'','missing')/)

!------
! table
!------
type(t_table) ,pointer :: radiation_corr

!==============================================================================
! surf_char: model surface characteristics (bit pattern)
!==============================================================================

integer ,parameter :: MS_LAND    = 0  ! ( 1)
integer ,parameter :: MS_SEA     = 1  ! ( 2)
integer ,parameter :: MS_ICE     = 2  ! ( 4)
integer ,parameter :: MS_NO_ICE  = 3  ! ( 8)
integer ,parameter :: MS_SNOW    = 4  ! (16)
integer ,parameter :: MS_NO_SNOW = 5  ! (32)

type (e) ,target :: surf_char_entries (6) = (/                              &
e(MS_LAND   ,'LAND'   ,'','set if some fraction is covered by land'       ),&
e(MS_SEA    ,'SEA'    ,'','set if some fraction is covered by sea'        ),&
e(MS_ICE    ,'ICE'    ,'','set if some fraction is covered by sea-ice'    ),&
e(MS_NO_ICE ,'NO_ICE' ,'','set if some fraction is not covered by sea-ice'),&
e(MS_SNOW   ,'SNOW'   ,'','set if some fraction is covered by snow'       ),&
e(MS_NO_SNOW,'NO_SNOW','','set if some fraction is not covered by snow'   )/)

type(t_table) ,pointer :: surf_char

!=======================================
! ct_nwc cloud type according to NWC SAF
!=======================================

integer ,parameter :: CT_NOPR         =  0 ! non-processed containing no data or corrupted data
integer ,parameter :: CT_LAND_FREE    =  1 ! cloud free land
integer ,parameter :: CT_SEA_FREE     =  2 ! cloud free sea
integer ,parameter :: CT_LAND_SNOW    =  3 ! land contaminated by snow
integer ,parameter :: CT_SEA_ICE      =  4 ! sea contaminated by snow/ice
integer ,parameter :: CT_CUM_VLOW     =  5 ! very low and cumuliform clouds
integer ,parameter :: CT_STRAT_VLOW   =  6 ! very low and stratiform clouds
integer ,parameter :: CT_LOW_CUM      =  7 ! low and cumuliform clouds
integer ,parameter :: CT_LOW_STRAT    =  8 ! low and stratiform clouds
integer ,parameter :: CT_MED_CUM      =  9 ! medium and cumuliform clouds
integer ,parameter :: CT_MED_STRAT    = 10 ! medium and stratiform clouds
integer ,parameter :: CT_HI_OP_CUM    = 11 ! high opaque and cumuliform clouds
integer ,parameter :: CT_HI_OP_STRAT  = 12 ! high opaque and stratiform clouds
integer ,parameter :: CT_VHI_OP_CUM   = 13 ! very high opaque and cumuliform clouds
integer ,parameter :: CT_VHI_OP_STRAT = 14 ! very high opaque and stratiform clouds
integer ,parameter :: CT_HI_ST_THIN   = 15 ! high semitransparent thin clouds
integer ,parameter :: CT_HI_ST_MEAN   = 16 ! high semitransparent meanly thick clouds
integer ,parameter :: CT_HI_ST_THICK  = 17 ! high semitransparent thick clouds
integer ,parameter :: CT_HI_ST_ABOVE  = 18 ! high semitransparent above low or medium clouds
integer ,parameter :: CT_FRAC         = 19 ! fractional clouds (sub-pixel water clouds)
integer ,parameter :: CT_UNDEF        = 20 ! undefined (undefined by CMa)

type (e) ,target :: ct_nwc_entries (21) = (/                                              &
e(CT_NOPR        ,'NOPR'        ,'','non-processed containing no data or corrupted data'),&
e(CT_LAND_FREE   ,'LAND_FREE'   ,'','cloud free land'                                   ),&
e(CT_SEA_FREE    ,'SEA_FREE'    ,'','cloud free sea'                                    ),&
e(CT_LAND_SNOW   ,'LAND_SNOW'   ,'','land contaminated by snow'                         ),&
e(CT_SEA_ICE     ,'SEA_ICE'     ,'','sea contaminated by snow/ice'                      ),&
e(CT_CUM_VLOW    ,'CUM_VLOW'    ,'','very low and cumuliform clouds'                    ),&
e(CT_STRAT_VLOW  ,'STRAT_VLOW'  ,'','very low and stratiform clouds'                    ),&
e(CT_LOW_CUM     ,'LOW_CUM'     ,'','low and cumuliform clouds'                         ),&
e(CT_LOW_STRAT   ,'LOW_STRAT'   ,'','low and stratiform clouds'                         ),&
e(CT_MED_CUM     ,'MED_CUM'     ,'','medium and cumuliform clouds'                      ),&
e(CT_MED_STRAT   ,'MED_STRAT'   ,'','medium and stratiform clouds'                      ),&
e(CT_HI_OP_CUM   ,'HI_OP_CUM'   ,'','high opaque and cumuliform clouds'                 ),&
e(CT_HI_OP_STRAT ,'HI_OP_STRAT' ,'','high opaque and stratiform clouds'                 ),&
e(CT_VHI_OP_CUM  ,'VHI_OP_CUM'  ,'','very high opaque and cumuliform clouds'            ),&
e(CT_VHI_OP_STRAT,'VHI_OP_STRAT','','very high opaque and stratiform clouds'            ),&
e(CT_HI_ST_THIN  ,'HI_ST_THIN'  ,'','high semitransparent thin clouds'                  ),&
e(CT_HI_ST_MEAN  ,'HI_ST_MEAN'  ,'','high semitransparent meanly thick clouds'          ),&
e(CT_HI_ST_THICK ,'HI_ST_THICK' ,'','high semitransparent thick clouds'                 ),&
e(CT_HI_ST_ABOVE ,'HI_ST_ABOVE' ,'','high semitransparent above low or medium clouds'   ),&
e(CT_FRAC        ,'FRAC'        ,'','fractional clouds (sub-pixel water clouds)'        ),&
e(CT_UNDEF       ,'UNDEF'       ,'','undefined (undefined by CMa)'                      )/)

type(t_table) ,pointer :: ct_nwc

!==============================================================================
! surftype (for radiances)
!==============================================================================

! Just as MS_* (mdlsfc) bits:
integer, parameter :: SUR_LAND     =  0 !    (1)  Land
integer, parameter :: SUR_SEA      =  1 !    (2)  Sea
integer, parameter :: SUR_ICE      =  2 !    (4)  Ice/snow detected
integer, parameter :: SUR_NOICE    =  3 !    (8)  ice-free sea
integer, parameter :: SUR_SNOW     =  4 !   (16)  snow detected
integer, parameter :: SUR_NOSNOW   =  5 !   (32)  snow-free land
! Additional stuff:
integer, parameter :: SUR_HIGHLAND =  6 !   (64)  High land
integer, parameter :: SUR_MWSURF   =  7 !  (128)  MW-surface test
integer, parameter :: SUR_BLK_PP   =  8 !  (256)  Blacklisted by sat_pp
integer, parameter :: SUR_MISMATCH =  9 !  (512)  Surface mismatch
integer, parameter :: SUR_MISSING  = 10 ! (1024)  missing surface classification

type (e) ,target :: surftype_entries (11) = (/&
e(SUR_LAND     ,'LAND'    ,'' , 'land'                                      ),&
e(SUR_SEA      ,'SEA'     ,'' , 'sea'                                       ),&
e(SUR_ICE      ,'ICE'     ,'' , 'seaice'                                    ),&
e(SUR_NOICE    ,'NOICE'   ,'' , 'ice-free sea'                              ),&
e(SUR_SNOW     ,'SNOW'    ,'' , 'snow'                                      ),&
e(SUR_NOSNOW   ,'NOSNOW'  ,'' , 'snow-free land'                            ),&
e(SUR_HIGHLAND ,'HIGHLAND','' , 'highland'                                  ),&
e(SUR_MWSURF   ,'MWSURF'  ,'' , 'bad surface type detected by MW sounder'   ),&
e(SUR_BLK_PP   ,'BLK_PP'  ,'' , 'surface blacklisted by sat_pp'             ),&
e(SUR_MISMATCH ,'MISMATCH','' , 'mismatch of surface types'                 ),&
e(SUR_MISSING  ,'MISSING' ,'' , 'no surftype available'                     )/)

! OLD (pre DACE-2.13):
! integer, parameter :: SUR_SEA      =  0 !    (1)  Sea
! integer, parameter :: SUR_BLK_PP   =  1 !    (2)  Blacklisted by sat_pp
! integer, parameter :: SUR_LAND     =  2 !    (4)  Land
! integer, parameter :: SUR_HIGHLAND =  3 !    (8)  High land
! integer, parameter :: SUR_MISMATCH =  4 !   (16)  Surface mismatch
! integer, parameter :: SUR_ICE      =  5 !   (32)  Ice/snow detected
! integer, parameter :: SUR_NOICE    =  6 !   (64)  ice-free sea
! integer, parameter :: SUR_MISSING  =  7 !  (128)  missing surface classification
! integer, parameter :: SUR_MWSURF   =  8 !  (256)  MW-surface test

! type (e) ,target :: surftype_entries (9) = (/&
! e(SUR_SEA      ,'SEA'     ,'' , 'sea'                                       ),&
! e(SUR_BLK_PP   ,'BLK_PP'  ,'' , 'surface blacklisted by sat_pp'             ),&
! e(SUR_LAND     ,'LAND'    ,'' , 'land'                                      ),&
! e(SUR_HIGHLAND ,'HIGHLAND','' , 'highland'                                  ),&
! e(SUR_MISMATCH ,'MISMATCH','' , 'mismatch of surface types'                 ),&
! e(SUR_ICE      ,'ICE'     ,'' , 'seaice'                                    ),&
! e(SUR_NOICE    ,'NOICE'   ,'' , 'ice-free sea'                              ),&
! e(SUR_MISSING  ,'MISSING' ,'' , 'no surftype available'                     ),&
! e(SUR_MWSURF   ,'MWSURF'  ,'' , 'bad surface type detected by MW sounder'   )/)

type(t_table) ,pointer :: surftype

!==============================================================================
! soiltype
!==============================================================================

integer, parameter :: SOT_ICE      =  1   ! ice
integer, parameter :: SOT_ROCK     =  2   ! rock
integer, parameter :: SOT_SAND     =  3   ! sand
integer, parameter :: SOT_SANDLOAM =  4   ! sandy loam
integer, parameter :: SOT_LOAM     =  5   ! loam
integer, parameter :: SOT_CLAYLOAM =  6   ! clay loam
integer, parameter :: SOT_CLAY     =  7   ! clay
integer, parameter :: SOT_PEAT     =  8   ! peat
integer, parameter :: SOT_SEAWATER =  9   ! sea water
integer, parameter :: SOT_SEAICE   = 10   ! sea ice

type (e) ,target :: soiltype_entries (10) = (/&
e(SOT_ICE       ,'ICE'      ,'' ,''),&
e(SOT_ROCK      ,'ROCK'     ,'' ,''),&
e(SOT_SAND      ,'SAND'     ,'' ,''),&
e(SOT_SANDLOAM  ,'SANDLOAM' ,'' ,''),&
e(SOT_LOAM      ,'LOAM'     ,'' ,''),&
e(SOT_CLAYLOAM  ,'CLAYLOAM' ,'' ,''),&
e(SOT_CLAY      ,'CLAY'     ,'' ,''),&
e(SOT_PEAT      ,'PEAT'     ,'' ,''),&
e(SOT_SEAWATER  ,'SEAWATER' ,'' ,''),&
e(SOT_SEAICE    ,'SEAICE'   ,'' ,'')/)

type(t_table) ,pointer :: soiltype

!------------------------------------------------------------------------------
!
!integer ,parameter :: D1_DATA = 0  !
!integer ,parameter :: D1_MIN  = 1  ! 1dvar minimisation failed
!integer ,parameter :: D1_SUR  = 2  ! wrong surface type
!integer ,parameter :: D1_CLD  = 3  ! cloud flag
!
!type (e) ,target :: flg_1dvar_entries (4) = (/      &
!e(D1_DATA ,'DATA' ,'' ,'                         '),&
!e(D1_MIN  ,'MIN'  ,'' ,'1dvar minimisation failed'),&
!e(D1_SUR  ,'SUR'  ,'' ,'wrong surface type       '),&
!e(D1_CLD  ,'CLD'  ,'' ,'cloud flag               ')/)
!
!type(t_table) ,pointer :: flg_1dvar
!
!------------------------------------------------------------------------------

integer ,parameter :: CL_CLEAR     = 0
integer ,parameter :: CL_IR_CLOUDY = 1
integer ,parameter :: CL_MW_CLEAR  = 2
integer ,parameter :: CL_MW_CLOUDY = 3

type (e) ,target :: flg_cld_entries (4) = (/&
e(CL_CLEAR     ,'CLEAR'     ,'' ,''),&
e(CL_IR_CLOUDY ,'IR_CLOUDY' ,'' ,''),&
e(CL_MW_CLEAR  ,'MW_CLEAR'  ,'' ,''),&
e(CL_MW_CLOUDY ,'MW_CLOUDY' ,'' ,'')/)

type(t_table) ,pointer :: flg_cld

!==============================================================================
! level significance
!==============================================================================
!----------
! constants
!----------
integer ,parameter :: LS_SURFACE       = 0 ! surface
integer ,parameter :: LS_STANDARD      = 1 ! standard level
integer ,parameter :: LS_TROPO         = 2 ! tropopause level
integer ,parameter :: LS_MAX           = 3 ! maximum wind level
integer ,parameter :: LS_SIGN          = 4 ! significant level
integer ,parameter :: LS_SUPEROBS      = 5 ! superobservation (layer average)

!--------------
! table entries
!--------------
type (e) ,target :: level_sig_entries (6) = (/                               &
e(LS_SURFACE       ,'SURFACE'       ,'' ,'surface                         '),&
e(LS_STANDARD      ,'STANDARD'      ,'' ,'standard level                  '),&
e(LS_TROPO         ,'TROPO'         ,'' ,'tropopause level                '),&
e(LS_MAX           ,'MAX'           ,'' ,'maximum wind level              '),&
e(LS_SIGN          ,'SIGN'          ,'' ,'significant level               '),&
e(LS_SUPEROBS      ,'SUPEROBS'      ,'' ,'superobservation (layer average)')/)

!------
! table
!------
type(t_table) ,pointer :: level_sig

!==============================================================================
! phase of aircraft flight and roll angle quality
!==============================================================================
!----------
! constants
!----------
!                               0,1 ! Reserved
integer ,parameter :: PH_UNS  = 2   ! Unsteady
integer ,parameter :: PH_LVR  = 3   ! Level flight, routine observation
integer ,parameter :: PH_LVW  = 4   ! Level flight, highest wind encountered
integer ,parameter :: PH_ASC  = 5   ! Ascending
integer ,parameter :: PH_DES  = 6   ! Descending
integer ,parameter :: PH_MIS  = 7   ! Missing value

integer ,parameter :: RA_GOOD = 0   ! Good
integer ,parameter :: RA_BAD  = 1   ! Bad
!                               2   ! Reserved
integer ,parameter :: RA_MIS  = 3   ! Missing value

!--------------
! table entries
!--------------
type (e) ,target :: phase_entries (6) = (/                     &
e(PH_UNS  ,'UNS' ,'','Unsteady                              '),&
e(PH_LVR  ,'LVR' ,'','Level flight, routine observation     '),&
e(PH_LVW  ,'LVW' ,'','Level flight, highest wind encountered'),&
e(PH_ASC  ,'ASC' ,'','Ascending                             '),&
e(PH_DES  ,'DES' ,'','Descending                            '),&
e(PH_MIS  ,'MIS' ,'','Missing                               ')/)

type (e) ,target :: rollangle_entries (3) = (/&
e(RA_GOOD ,'GOOD','','Good                 '),&
e(RA_BAD  ,'BAD' ,'','Bad                  '),&
e(RA_MIS  ,'MIS' ,'','Missing value        ')/)

!-------
! tables
!-------
type(t_table) ,pointer :: phase
type(t_table) ,pointer :: rollangle

!==============================================================================
! retrtype: for AMV satellite derived wind computation method
!==============================================================================
!----------
! constants
!----------
!Wind derived from cloud motion observed in the ..
integer ,parameter :: RT_IR   =   1 ! infrared channel
integer ,parameter :: RT_VIS  =   2 ! visible channel
integer ,parameter :: RT_WV   =   3 ! water vapour channel
!                                 4 ! combination of spectral chans
!                                 5 ! water vapour ch. in clear air
!                                 6 ! ozon channel
!                                 7 ! water vap.ch., cloud or clear not spec.
integer ,parameter :: RT_IR1  = 101 !  8.7 um Meteos8-9, 10.7 um GOES 10-12
integer ,parameter :: RT_IR2  = 201 !  9.7 um Meteos8-9,  3.9 um GOES 10-12
integer ,parameter :: RT_IR3  = 301 ! 10.8 um Meteos8-9
integer ,parameter :: RT_VIS1 = 102 !  0.6 um Meteos8-9,  0.65um GOES 10-12
integer ,parameter :: RT_VIS2 = 302 !  0.75um Meteos8-9
integer ,parameter :: RT_VIS3 = 202 !  0.8 um Meteos8-9
integer ,parameter :: RT_WV1  = 103 !  6.2 um Meteos8-9,  7.4 um GOES 10-12
integer ,parameter :: RT_WV2  = 203 !  7.3 um Meteos8-9,  7.0 um GOES 10-12
integer ,parameter :: RT_WV3  = 303 !                     6.8 um GOES 10-12

!--------------
! table entries
!--------------
type (e) ,target :: retrtype_entries (12) = (/                        &
e(RT_IR  ,'IR     ','','infrared channel                           '),&
e(RT_VIS ,'VIS    ','','visible channel                            '),&
e(RT_WV  ,'WV     ','','water vapour channel                       '),&
e(RT_IR1 ,'IR1    ','',' 8.7 um Meteosat 8-11, ~ 10.7 um GOES 10-17'),&
e(RT_IR2 ,'IR2    ','',' 9.7 um Meteosat 8-11, ~  3.9 um GOES 10-17'),&
e(RT_IR3 ,'IR3    ','','10.8 um Meteosat 8-11                      '),&
e(RT_VIS1,'VIS1   ','',' 0.6 um Meteosat 8-11, ~  0.65um GOES 10-17, 0.65um Himawari 8-9'),&
e(RT_VIS3,'VIS3   ','',' 0.8 um Meteosat 8-11                      '),&
e(RT_VIS2,'VIS2-HR','',' 0.75um Meteosat 8-11                      '),&
e(RT_WV1 ,'WV1    ','',' 6.2 um Meteosat 8-11, ~  7.4 um GOES 10-17, 6.25um Himawari 8-9'),&
e(RT_WV2 ,'WV2    ','',' 7.3 um Meteosat 8-11, ~  7.0 um GOES 10-17, 6.95um Himawari 8-9'),&
e(RT_WV3 ,'WV3    ','','                    ~ 6.8-6.2 um GOES 10-17, 7.35um Himawari 8-9')/)

!-------
! tables
!-------
type(t_table) ,pointer :: retrtype

!==============================================================================
! specification of verification data (runtype, ensmem)
!==============================================================================

integer ,parameter :: VT_FORECAST    = 0 ! forecast
integer ,parameter :: VT_FIRSTGUESS  = 1 ! first guess
integer ,parameter :: VT_PREL_ANA    = 2 ! preliminary analysis
integer ,parameter :: VT_ANALYSIS    = 3 ! analysis
integer ,parameter :: VT_INIT_ANA    = 4 ! initialised analysis
integer ,parameter :: VT_LIN_ANA     = 5 ! linear operator on analysis
integer ,parameter :: VT_FC_SENS     = 6 ! forecast sensitivity

type (e) ,target :: runtype_entries (7) = (/                                   &
e(VT_FORECAST   ,'FORECAST   ','','forecast                                 '),&
e(VT_FIRSTGUESS ,'FIRSTGUESS ','','first guess                              '),&
e(VT_PREL_ANA   ,'PREL_ANA   ','','preliminary analysis in observation space'),&
e(VT_ANALYSIS   ,'ANALYSIS   ','','analysis                                 '),&
e(VT_INIT_ANA   ,'INIT_ANA   ','','initialised analysis                     '),&
e(VT_LIN_ANA    ,'LIN_ANA    ','','linear operator on analysis (Y_a)        '),&
e(VT_FC_SENS    ,'FC_SENS    ','','forecast sensitivity                     ')/)

type(t_table) ,pointer :: runtype

!------------------------------------------------------------------------------
integer ,parameter :: VE_ENS_MEAN     =  0 ! ensemble mean
integer ,parameter :: VE_DETERM       = -1 ! deterministic model run
integer ,parameter :: VE_ENS_SPREAD   = -2 ! ensemble spread
integer ,parameter :: VE_BG_ERROR     = -3 ! 3dvar background error
integer ,parameter :: VE_TALAGRAND    = -4 ! Talagrand index
integer ,parameter :: VE_VQC_WEIGHT   = -5 ! variational quality control weight
integer ,parameter :: VE_MEMBER       = -6 ! generic value for ensemble member
integer ,parameter :: VE_ENS_MEAN_OBS = -7 ! ensemble mean in observation space
integer ,parameter :: VE_BIASCOR      = -8 ! bias correction applied

type (e) ,target :: ensmem_entries (9) = (/                               &
e(VE_ENS_MEAN    ,'ENS_MEAN    ','','ensemble mean                     '),&
e(VE_DETERM      ,'DETERM      ','','deterministic model run           '),&
e(VE_ENS_SPREAD  ,'ENS_SPREAD  ','','ensemble spread                   '),&
e(VE_BG_ERROR    ,'BG_ERROR    ','','3dvar background error            '),&
e(VE_TALAGRAND   ,'TALAGRAND   ','','Talagrand index                   '),&
e(VE_VQC_WEIGHT  ,'VQC_WEIGHT  ','','variational quality control weight'),&
e(VE_MEMBER      ,'MEMBER      ','','generic value for ensemble member '),&
e(VE_ENS_MEAN_OBS,'ENS_MEAN_OBS','','ensemble mean in observation space'),&
e(VE_BIASCOR     ,'BIASCOR     ','','bias correction                   ')/)
type(t_table) ,pointer :: ensmem


!==============================================================================
! observation operator flag (observation type dependent)
!==============================================================================
integer ,parameter :: OF_MISSING       = 0 ! missing (not present in old files)
integer ,parameter :: OF_RAD_CLEAR_SKY = 1 ! clear sky radiances
integer ,parameter :: OF_RAD_CLOUDY    = 2 ! cloudy radiances
integer ,parameter :: OF_BT_CLEAR_SKY  = 3 ! clear sky radiances
type (e) ,target :: oflag_entries (4) = (/                   &
e(OF_MISSING      ,'MISSING      ','','missing value      '),&
e(OF_RAD_CLEAR_SKY,'RAD_CLEAR_SKY','','clear sky radiances'),&
e(OF_RAD_CLOUDY   ,'RAD_CLOUDY   ','','cloudy radiances   '),&
e(OF_BT_CLEAR_SKY ,'BT_CLEAR_SKY ','','clear sky br.temp. ')/)
type(t_table) ,pointer :: oflag

!==============================================================================
! TOVS specific flags
!==============================================================================
integer :: i
integer ,parameter :: TF_EMIS(0:7)     =  (/ (i, i=0,7) /) ! emissivity flags
integer ,parameter :: TF_EMIS_DYNRET   =  TF_EMIS(1) ! emissivity calculated with dynamical retrieval
integer ,parameter :: TF_EMIS_ATLS     =  TF_EMIS(5) ! emissivity from other atlas
integer ,parameter :: TF_EMIS_PC_RETR  =  TF_EMIS(6) ! emissivity from PC retrievel
integer ,parameter :: TF_EMIS_FAILED   =  TF_EMIS(7) ! emissivity calc. failed
                                                     ! MW emissivity calc. methods
integer ,parameter :: TF_EMIS_FASTEM   =  TF_EMIS(0) ! emissivity calculated with FASTEM
integer ,parameter :: TF_EMIS_GRODY    =  TF_EMIS(2) ! emissivity calculated with Grody method
integer ,parameter :: TF_EMIS_TLSM     =  TF_EMIS(3) ! emissivity from TELSEM atlas
integer ,parameter :: TF_EMIS_CNRM     =  TF_EMIS(4) ! emissivity from CNRM atlas
                                                     ! IR emissivity calc. methods
integer ,parameter :: TF_EMIS_SEA_MOD  =  TF_EMIS(0) ! emissivity calculated with IR sea model
integer ,parameter :: TF_EMIS_CAMEL07  =  TF_EMIS(2) ! emissivity from CAMEL2007 atlas
integer ,parameter :: TF_EMIS_CAMELCL  =  TF_EMIS(3) ! emissivity from CAMEL climatology atlas
integer ,parameter :: TF_EMIS_UWIR     =  TF_EMIS(4) ! emissivity from UWIR atlas
                                                     ! VIS refl. calc. method
integer ,parameter :: TF_REFL_BRDF     =  TF_EMIS(5) ! BRDF atlas
integer ,parameter :: TF_SURF_INFL     =  8 ! surface influence too large
integer ,parameter :: TF_SURF_RETR     =  9 ! wrong surface type detected in cloud check
integer ,parameter :: TF_SURF_TYPE     = 10 ! wrong surface type (not activated in namelist)
integer ,parameter :: TF_SURF_MODEL    = 11 ! surface properties (given by model) not adequate
integer ,parameter :: TF_CLOUD(6)      = (/ (11+i, i=1, 6) /)
integer ,parameter :: TF_AEROSOL       = 18
integer ,parameter :: TF_DESERT_DUST   = 19
integer ,parameter :: TF_VOLCANIC_ASH  = 20
integer ,parameter :: TF_TRACE_GAS     = 21
integer ,parameter :: TF_BC_FAILED     = 22
integer ,parameter :: TF_TSKIN(0:1)    = (/23,24/)
integer ,parameter :: TF_TSKIN_DYNRET  = TF_TSKIN(0)
integer ,parameter :: TF_TSKIN_FAILED  = TF_TSKIN(1)
type (e) ,target :: tovsflag_entries (25) = (/                                           &
e(TF_EMIS(1)      ,'EMIS_SEA_MOD ','','emiss. calc. with sea mod. (FASTEM,ISEM,...)     '),&
e(TF_EMIS_DYNRET  ,'EMIS_DYNRET  ','','emiss. calc. with dynamical retrieval            '),&
e(TF_EMIS(2)      ,'EMIS_2       ','','emiss. from Grody(MW), CAMEL2007(IR)             '),&
e(TF_EMIS(3)      ,'EMIS_3       ','','emiss. from TELSEM(MW), CAMEL clim.(IR) atlas    '),&
e(TF_EMIS(4)      ,'EMIS_4       ','','emiss. from CNRM(MW), UWIR(IR) atlas             '),&
e(TF_EMIS(5)      ,'EMIS_5       ','','emiss. from BRDF(VIS) or other atlas             '),&
e(TF_EMIS_PC_RETR ,'EMIS_PC_RETR ','','emiss. from PC retrieval                         '),&
e(TF_EMIS_FAILED  ,'EMIS_FAILED  ','','emiss./reflectance calc. failed                  '),&
e(TF_SURF_INFL    ,'SURF_INFL    ','','surface influence too large                      '),&
e(TF_SURF_RETR    ,'SURF_RETR    ','','wrong surface type retrieved                     '),&
e(TF_SURF_TYPE    ,'SURF_TYPE    ','','wrong surface type (not activated in namelist)   '),&
e(TF_SURF_MODEL   ,'SURF_MODEL   ','','surface properties (given by model) not adequate '),&
e(TF_CLOUD(1)     ,'CLOUD1       ','','cloud check 1                                    '),&
e(TF_CLOUD(2)     ,'CLOUD2       ','','cloud check 2                                    '),&
e(TF_CLOUD(3)     ,'CLOUD3       ','','cloud check 3                                    '),&
e(TF_CLOUD(4)     ,'CLOUD4       ','','cloud check 4                                    '),&
e(TF_CLOUD(5)     ,'CLOUD5       ','','cloud check 5                                    '),&
e(TF_CLOUD(6)     ,'CLOUD6       ','','cloud check 6                                    '),&
e(TF_AEROSOL      ,'AEROSOL      ','','aerosol affected                                 '),&
e(TF_DESERT_DUST  ,'DESERT_DUST  ','','desert dust affected                             '),&
e(TF_VOLCANIC_ASH ,'VOLCANIC_ASH ','','volcanic ash affected                            '),&
e(TF_TRACE_GAS    ,'TRACE_GAS    ','','trace gas affected                               '),&
e(TF_TSKIN_DYNRET ,'TSKIN_DYNRET ','','skin temperature calc. with dynamical retrieval  '),&
e(TF_TSKIN_FAILED ,'TSKIN_FAILED ','','skin temperature calc. failed                    '),&
e(TF_BC_FAILED    ,'BC_FAILED    ','','biascorrection failed                            ')/)
type(t_table) ,pointer :: tovsflag


!==============================================================================
! optional TOVS specific variables
!==============================================================================

! Please, make sure that this is consistent with the OPTV_* parameters in
! mo_rad.f90
integer, parameter :: n_optv           = 6
character(len=7)   :: c_optv(0:5)
data c_optv(0)/'L2C'    /
data c_optv(1)/'NWC_FLG'/
data c_optv(2)/'ORB_PH' /
data c_optv(3)/'INS_TMP'/
data c_optv(4)/'CLD_FRC'/
data c_optv(5)/'CLD_FLG'/

!==============================================================================
contains
!==============================================================================
  subroutine clean_fdbk_tables
    if (.not.init) return
    init = .false.
    deallocate (status)
    deallocate (flags)
    deallocate (obstype)
    deallocate (obstype_oce)
    deallocate (codetype)
    deallocate (varno)
    deallocate (satsens)
    deallocate (rsondtype)
    deallocate (trackteqn)
    deallocate (meas_equip)
    deallocate (radiation_corr)
    deallocate (surftype)
    deallocate (soiltype)
!   deallocate (flg_1dvar)
    deallocate (surf_char)
    deallocate (flg_cld)
    deallocate (level_sig)
    deallocate (phase)
    deallocate (rollangle)
    deallocate (retrtype)
    deallocate (runtype)
    deallocate (runclass)
    deallocate (ensmem)
    deallocate (ind_cg)
    deallocate (gen_cg)
    deallocate (oflag)
  end subroutine clean_fdbk_tables
!------------------------------------------------------------------------------
  subroutine init_fdbk_tables (latex, csv)
  logical, intent(in), optional :: latex
  logical, intent(in), optional :: csv
  !-------------------------------------------------------------------
  ! join table entries and table meta data (name and caption).
  !   optionally (if 'latex' is given and true)
  !   write LaTeX file with tables for inclusion in the documentation.
  !-------------------------------------------------------------------

    if (init) return
    init = .true.

    call init_table (status                                ,&
                     status_entries                        ,&
                     'status'                              ,&
                     'Observation or report status values' ,&
                     .false.                               ,&
                     latex                                 ,&
                     csv)

    call init_table (flags                         ,&
                     flag_entries                  ,&
                     'flags'                       ,&
                     'Report or Observation Flags' ,&
                     .true.                        ,&
                     latex                         ,&
                     csv)

    call init_table (obstype             ,&
                     obstype_entries     ,&
                     'obstype'           ,&
                     'Observation Types' ,&
                     .false.             ,&
                     latex               ,&
                     csv)

    call init_table (obstype_oce                       ,&
                     obstype_oce_entries               ,&
                     'obstype_oce'                     ,&
                     'Oceanographic Observation Types' ,&
                     .false.                           ,&
                     latex                             ,&
                     csv)

    call init_table (codetype                 ,&
                     codetype_entries         ,&
                     'codetype'               ,&
                     'Observation Code Types' ,&
                     .false.                  ,&
                     latex                    ,&
                     csv)

    call init_table (varno                                    ,&
                     varno_entries                            ,&
                     'varno'                                  ,&
                     'Observation and level variable numbers' ,&
                     .false.                                  ,&
                     latex                                    ,&
                     csv)

    call init_table (satsens                       ,&
                     satsens_entries               ,&
                     'sat-sensor'                  ,&
                     'RTTOV instrument ID (sensor)',&
                     .false.                       ,&
                     latex                         ,&
                     csv)

    call init_table (rsondtype                                          ,&
                     rsondtype_entries                                  ,&
                     'rsondtype'                                        ,&
'radiosonde type (NRARA), for other values see WMO common code table C2',&
                     .false.                                            ,&
                     latex                                              ,&
                     csv)

    call init_table (trackteqn                                             ,&
                     trackteqn_entries                                     ,&
                     'trackteqn'                                           ,&
'tracking technique (NSASA), for other values see WMO common code table C7',&
                     .false.                                               ,&
                     latex                                                 ,&
                     csv)

    call init_table (meas_equip                                      ,&
                     meas_equip_entries                              ,&
                     'measequip'                                     ,&
                     'type of measuring equipment used (NA4, 002003)',&
                      .false.                                        ,&
                      latex                                          ,&
                      csv)

    call init_table (radiation_corr                                         ,&
                     radiation_corr_entries                                 ,&
                     'radiationcorr'                                        ,&
                     'solar and infrared radiation correction (NSR, 002013)',&
                     .false.                                                ,&
                      latex                                                 ,&
                      csv)

    call init_table (surftype                               ,&
                     surftype_entries                       ,&
                     'surftype'                             ,&
                     'surface types consistent with 1d-Var' ,&
                     .true.                                 ,&
                     latex                                  ,&
                     csv)

    call init_table (soiltype                               ,&
                     soiltype_entries                       ,&
                     'soiltype'                             ,&
                     'type of soil in the retrieval point'  ,&
                     .false.                                ,&
                     latex                                  ,&
                     csv)

!   call init_table (flg_1dvar               ,&
!                    flg_1dvar_entries       ,&
!                    'flg_1dvar'             ,&
!                    '1dvar processing flag' ,&
!                    .true.                  ,&
!                    latex                   ,&
!                    csv)

    call init_table (surf_char                      ,&
                     surf_char_entries              ,&
                     'surf_char'                    ,&
                     'model surface characteristics',&
                     .true.                         ,&
                     latex                          ,&
                     csv)

    call init_table (flg_cld               ,&
                     flg_cld_entries       ,&
                     'flg_cld'             ,&
                     'cloud flag'          ,&
                     .true.                ,&
                     latex                 ,&
                     csv)

    call init_table (level_sig                           ,&
                     level_sig_entries                   ,&
                     'level_sig'                         ,&
                     'level significance for TEMP/PILOT.',&
                     .true.                              ,&
                     latex                               ,&
                     csv)

    call init_table (phase                                     ,&
                     phase_entries                             ,&
                     'phase'                                   ,&
                     'aircraft phase, radiances fov, gpsro pcd',&
                     .false.                                   ,&
                     latex                                     ,&
                     csv)

    call init_table (rollangle                        ,&
                     rollangle_entries                ,&
                     'rollangle'                      ,&
                     'aircraft roll angle'            ,&
                     .false.                          ,&
                     latex                            ,&
                     csv)

    call init_table (retrtype                                    ,&
                     retrtype_entries                            ,&
                     'retrtype'                                  ,&
                     'satellite derived wind computation method' ,&
                     .false.                                     ,&
                     latex                                       ,&
                     csv)

    call init_table (runtype                    ,&
                     runtype_entries            ,&
                     'runtype'                  ,&
                     'type of verification run' ,&
                     .false.                    ,&
                     latex                      ,&
                     csv)

    call init_table (runclass                    ,&
                     runclass_entries            ,&
                     'runclass'                  ,&
                     'class of verification run' ,&
                     .false.                     ,&
                     latex                       ,&
                     csv)

    call init_table (ensmem                                   ,&
                     ensmem_entries                           ,&
                     'ensmem'                                 ,&
                     'specification of the verification data' ,&
                     .false.                                  ,&
                     latex                                    ,&
                     csv)

    call init_table (ind_cg                                 ,&
                     ind_cg_entries                         ,&
                     'ind_cg'                               ,&
                     'individual cloud group bit positions ',&
                     .false.                                ,&
                     latex                                  ,&
                     csv)

    call init_table (gen_cg                             ,&
                     gen_cg_entries                     ,&
                     'gen_cg'                           ,&
                     'general cloud group bit positions',&
                     .false.                            ,&
                     latex                              ,&
                     csv)

    call init_table (oflag                                                    ,&
                     oflag_entries                                            ,&
                     'oflag'                                                  ,&
                     'Observation operator flag (observation type dependent)' ,&
                     .false.                                                  ,&
                     latex                                                    ,&
                     csv)

    call init_table (tovsflag             ,&
                     tovsflag_entries     ,&
                     'tovsflag'           ,&
                     'TOVS specific flag' ,&
                     .false.              ,&
                     latex                ,&
                     csv)

    call init_table (ct_nwc                            ,&
                     ct_nwc_entries                    ,&
                     'ct_nwc'                          ,&
                     'Cloud type according to NWC SAF' ,&
                     .false.                           ,&
                     latex                             ,&
                     csv)

  end subroutine init_fdbk_tables
!==============================================================================
end module mo_fdbk_tables
