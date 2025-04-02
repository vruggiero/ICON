!
!+ observation error covariance model
!
MODULE mo_obs_err
!
! Description:
!   Observation error covariance model.
!   Provides observational errors for the observed quantities.
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
! V1_4         2009/03/26 Andreas Rhodin
!  errors for ASCAT/QSCAT; changes for feedback file input
! V1_8         2009/12/09 Harald Anlauf
!  read_nml: print diagnostics if observation error not present in given table
!            implement "extern" tables of observation errors; minor cleanups
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Changes for Oceansat-2; properly use obs.error for TEMP and PILOT
! V1_15        2011/12/06 Andreas Rhodin
!  option to specify codetypes for observation error tables
! V1_16        2011/12/09 Andreas Rhodin
!  changed comment lines
! V1_19        2012-04-16 Alexander Cress
!  reduce observation error for OSCAT
! V1_22        2013-02-13 Andreas Rhodin
!  changes for COSMO observation operators
! V1_43        2015-08-19 Andreas Rhodin
!  remove dependency from external coefficient file 'rszcoef_f'
! V1_44        2015-09-30 Harald Anlauf
!  Preparations for Jason-2
! V1_45        2015-12-15 Andreas Rhodin
!  generalize observation error tables to all 'obstype's and 'varno's
! V1_46        2016-02-05 Harald Anlauf
!  bug fix for SHIP obs.error; more reasonable defaults for SYNOP/DRIBU t2m
! V1_47        2016-06-06 Andreas Rhodin
!  specific observation error table to be used for adaptive localisation
! V1_48        2016-10-06 Andreas Rhodin
!  handle TD2M,FF,DD,cloud,..
! V1_50        2017-01-09 Andreas Rhodin
!  changes to run COMET observations thru 3dvar
! V1_51        2017-02-24 Andreas Rhodin
!  allow VN_U10M,VN_V10M for DRIBUs (for COMET);
!  pass observations from COSMO fof-files through the 3dvar
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! DWD OI observation error covariance model
!
! Authors:
! Andreas Rhodin  DWD  2002-2008  original code
!==============================================================================
  !=============
  ! Modules used
  !=============
  use mo_kind,      only: wp, sp, i8     ! kind parameters
  use mo_exception, only: finish,       &! abort routine
                          message        ! issue a warning
  use mo_mpi_dace,  only: dace,         &! MPI group info
                          p_bcast        ! broadcast routine
  use mo_namelist,  only: position_nml, &! routine to position nml group
                          nnml,         &! namelist fortran unit number
                          POSITIONED     ! position_nml: OK  return flag
  use mo_dwd,       only: rs             ! radiosonde height observ.error coeff.
  use mo_physics,   only: t0c,          &! 0 degree celsius
                          gacc,         &! gravity acceleration
                          r2d            ! factor radiand -> degree
  use mo_t_obs,     only: t_spot         !   component of t_obs
  use mo_t_use,     only: STAT_DISMISS   ! status flag value: dismissed
  use mo_fdbk_tables,only:varno,        &! variable number table
                          VN_U,         &!          wind component code
                          VN_V,         &!                         code
                          VN_U10M,      &!          wind component code
                          VN_V10M,      &!                         code
                          VN_FF,        &!          wind speed     code
                          VN_DD,        &!          wind direction code
                          VN_RH,VN_RH2M,&!          rh             code
                          VN_T,         &!          T              code
                          VN_T2M,       &!          T2m            code
                          VN_TD2M,      &!          2m dewpoint    code
                          VN_Z,         &!          geopotential   code
                          VN_P,         &!          pressure       code
                          VN_PS,        &!          surface pressure
                          VN_GUST,      &!          gust     (passive)
                          VN_RR,        &!          rain     (passive)
                          VN_TMIN,      &!          tmin     (passive)
                          VN_TMAX,      &!          tmax     (passive)
                          VN_PTEND,     &! pressure tendency (passive)
                          VN_N, VN_NH,  &!
                          VN_N_L,VN_N_M,&!
                          VN_N_H,       &!
                          VN_SDEPTH,    &!
                          VN_TSEA,      &!
                          VN_CEIL,      &! cloud ceiling       (passive)
                          VN_VV,        &! visibility          (passive)
                          VN_WW,        &! present weather code    code
                          VN_GCLG,      &! general cloud group     code
                          VN_ICLG,      &! individual cloud layer group
                          VN_RAD_GL,    &! global    radiation (passive)
                          VN_RAD_DF,    &! diffusive radiation (passive)
                          VN_RAD_LW,    &! longwave  radiation (passive)
                          OT_SYNOP,     &!
                          OT_AIREP,     &!
                          OT_SATOB,     &!
                          OT_DRIBU,     &!
                          OT_TEMP,      &!
                          OT_PILOT,     &!
                          OT_PAOB,      &!
                          OT_SCATT,     &!
                          init_fdbk_tables! initialise the tables
  use mo_t_table,   only: value_name,   &! find table entry value (name given)
                          name_value,   &! find name of table entry
                          INVALID_VALUE  ! returned by 'value_name'
  use mo_obs_set,   only: t_obs_block    !
  use mo_obs_tables,only: obstyp,       &!
                          rept_char      ! report type characteristics
  use mo_dace_string,only:toupper,      &! convert to upper case
                          char3          ! convert integer -> char(len=3)
  use mo_endian,    only: little,       &! .true. for little endian
                          flip           ! flip integer field
  use mo_wigos,     only: t_wsi,        &! WIGOS id data type
                          wsi_mode,     &! WIGOS station id mode
                          wsi_prefix,   &! Prefix of "fake station id"
                          wsi_decode,   &! Derive WSI from given string
                          wsi_encode,   &! statid,stathash from t_wsi
                          wsi_to_text,  &! Convert t_wsi to text format
                          operator(==)   ! Same WIGOS station ID?
  implicit none
  !================
  ! public entities
  !================
  private
  public :: temp_obs_err       ! set TEMP observation error covariances
  public :: synop_obs_err      ! set SYNOP, DRIBU, PAOB observation errors
  public :: scatt_obs_err      ! set SCATTerometer observation errors
  public :: init_obs_err       ! initialize module mo_obs_err
  public :: obs_err            ! function to get observation error from table
  public :: n_adp              ! number of entries in table adaptloc

  !=======
  ! tables
  !=======

  real(wp) ,parameter :: e_rh = 0.1_wp                        ! RH error ?

  !------------------------------------------------------------------------------
  !
  ! Options for using namelist /OBSERR/ :
  !
  !  1) Set specific observation errors for TEMP or SYNOP observations.
  !
  !  2) Set height (pressure, cf. array 'pt' below) dependent tables of
  !     observation errors for atmospheric in-situ data, characterised by
  !     observation type (TEMP, SYNOP, AMV,..), code type (optional) and
  !     variable type (t,q,...). The namelist group may be used
  !     repeatedly. Tables may be either explicitly specified in the namelist
  !     via the variable 'err' (use the key table='extern' for this case) or
  !     taken from predefined tables (use the key table='OI' or
  !     table='IFS'). In any case the values of a given table may be
  !     multiplied by the factor provided in 'scale'.
  !
  !  3) A table of observation errors to be used only in the estimation of
  !     the adaptive localisation scale may be given by specifying the key
  !     table='adaploc'.
  !
  !------------------------------------------------------------------------------
  integer  ,parameter :: ntab = 30               ! max. number of tables
  integer  ,parameter :: ncde = 10               ! max. number of codetypes
  integer  ,parameter :: nlev = 15               ! number of pressure levels
  real(wp) ,parameter :: pt (nlev) = (/100000., &! pressure levels (Pa)
                                        85000., &
                                        70000., &
                                        50000., &
                                        40000., &
                                        30000., &
                                        25000., &
                                        20000., &
                                        15000., &
                                        10000., &
                                         7000., &
                                         5000., &
                                         3000., &
                                         2000., &
                                         1000. /)
  !---------------------------------
  ! namelist parameters (SYNOP land)
  !---------------------------------
  real(wp)         :: obserr_h_reduction  =  0.01_wp ! h extrapolation error(1%)
  real(wp)         :: obserr_synopland_h  =  7._wp   ! geop.h.  SYNOP land
  real(wp)         :: obserr_synopauto_h  = -1._wp   ! geop.h.  SYNOP automatic
  real(wp)         :: obserr_coastal_h    =  7._wp   ! geop.h.  coastal stations
  real(wp)         :: obserr_synopland_rh = 0.40_wp  ! rel.hum  SYNOP land
  real(wp)         :: obserr_synophigh_rh = -9._wp   ! rel.hum  SYNOP 5000 m
  real(wp)         :: obserr_synopland_v  = 2.52_wp  ! windsp.  SYNOP land
  real(wp)         :: obserr_synopland_t  = 3.00_wp  ! 2m temp. SYNOP land
  real(wp)         :: obserr_synopland_dt = 0.00_wp  ! extrapolation error (K/m)
  !---------------------------------
  ! namelist parameters (SYNOP ship)
  !---------------------------------
  real(wp)         :: obserr_synopship_h  =  9._wp   ! geop.h.  SYNOP ship
  real(wp)         :: obserr_synopship_rh = 0.40_wp  ! rel.hum  SYNOP ship
  real(wp)         :: obserr_synopship_v  = 2.52_wp  ! windsp.  SYNOP ship
  real(wp)         :: obserr_synopship_t  = 2.00_wp  ! 2m temp. SYNOP ship
  !----------------------------
  ! namelist parameters (METAR)
  !----------------------------
  real(wp)         :: obserr_metar_h      = 10.0_wp  ! geop.h.  METAR
  real(wp)         :: obserr_metar_rh     = 0.15_wp  ! rel.hum  METAR
  real(wp)         :: obserr_metar_v      = 2.50_wp  ! windsp.  METAR
  real(wp)         :: obserr_metar_t      = 2.00_wp  ! 2m temp. METAR
  !---------------------------------
  ! namelist parameters (BUOY, PAOB)
  !---------------------------------
  real(wp)         :: obserr_paob_h       = 12._wp   ! geop.h.  PAOB
  real(wp)         :: obserr_dribu_h      =  7._wp   ! geop.h.  DRIBU
  real(wp)         :: obserr_dribu_rh     = 0.40_wp  ! rel.hum  DRIBU
  real(wp)         :: obserr_dribu_v      = 5.40_wp  ! windsp.  DRIBU
  real(wp)         :: obserr_dribu_t      = 2.00_wp  ! 2m temp. DRIBU
  !----------------------------
  ! namelist parameters (TEMPs)
  !----------------------------
  real(wp)         :: obserr_temp_rh      = 0.15_wp  ! rel.hum  TEMP
  real(wp)         :: obserr_temp_rh_dry  = 0.15_wp  ! below 20% rel.hum.
  real(wp)         :: obserr_temp_rh_cold = 0.20_wp  ! below -40 degree C
  !--------------------------------------
  ! namelist parameters (free atmosphere)
  !--------------------------------------
  character(len=8) :: obstype         ! observation type ('TEMP','PILOT',etc.)
  integer          :: codetype (ncde) ! code type        (numeric)
  character(len=8) :: quantity        ! quantity         ('uv','tv', 'rh')
  character(len=8) :: table           ! 'OI' or 'IFS'
  real(wp)         :: scale(nlev)     ! rescaling factor
  real(wp)         :: err  (nlev)     ! externally specified obs. errors
  real(wp)         :: r_err(nlev)     ! externally specified relative error
  !------------------
  ! tuning parameters
  !------------------
  real(wp)         :: a   != 0.5_wp               ! asathc in oi
  real(wp)         :: b   != 0.963621748665073_wp ! czrcor in oi: rszcoef

  namelist /OBSERR/ a, b,                          &! observation error model
                    obstype, codetype,             &! specifies table to use
                    quantity, table,               &! specifies table to use
                    scale, err, r_err,             &! error specification
                    obserr_synopship_h,            &!
                    obserr_synopland_h,            &!
                    obserr_synopauto_h,            &!
                    obserr_coastal_h,              &!
                    obserr_h_reduction,            &!
                    obserr_metar_h,                &!
                    obserr_paob_h, obserr_dribu_h, &!
                    obserr_synopship_rh,           &!
                    obserr_synopland_rh,           &!
                    obserr_synophigh_rh,           &!
                    obserr_metar_rh,               &!
                    obserr_dribu_rh,               &!
                    obserr_synopland_v,            &!
                    obserr_synopship_v,            &!
                    obserr_metar_v,                &!
                    obserr_dribu_v,                &!
                    obserr_synopland_t,            &!
                    obserr_synopland_dt,           &!
                    obserr_synopship_t,            &!
                    obserr_metar_t,                &!
                    obserr_dribu_t,                &!
                    obserr_temp_rh,                &!
                    obserr_temp_rh_dry,            &!
                    obserr_temp_rh_cold

  !--------------------------------------
  ! Observation error data type definiton
  !--------------------------------------
  type t_obserr
    integer  :: obstype         =  0     ! valid for observation types
    integer  :: codetype (ncde) = -1     ! ..for code types (-1: not specified)
    integer  :: varno           =  0     ! valid for parameter (varno)
    real(wp) :: err   (nlev)    =  0._wp ! observation error table
    real(wp) :: r_err (nlev)    =  0._wp ! relative observation error table
  end type t_obserr

  !---------------------------------------------
  ! Station-ID based observation error data type
  !---------------------------------------------
  type t_stn_obserr
    integer          :: code(4) = -1     ! code types (-1: not specified)
    character(len=8) :: statid  = ''     ! station id (wildcard)
    integer(i8)      :: mask    = 0      ! mask (for wildcards)
    integer(i8)      :: id      = 0      ! station id (masked)
    type(t_wsi)      :: wsi              ! WIGOS station ID
    integer          :: len     = 0      ! len_trim (station id)
    integer          :: varno   =  0     ! parameter (varno)
    real(wp)         :: err     =  0._wp ! observation error
!   real(wp)         :: par(2)  =  0._wp ! auxiliary parameter(s), reseverved
  end type t_stn_obserr

  logical ,save               :: bigend = .true. ! true for big endian
  integer ,save               :: nentry = 0      ! no. active individual entries
  type(t_stn_obserr), pointer :: stn_obserr(:) => NULL()

  !----------------------------------------------------------
  ! IFS observation errors (IFS documentation Part I (CY21R4)
  !----------------------------------------------------------
!#if defined (__SX__)
!type (t_obserr), parameter :: &
!IFS_99=t_obserr(OT_SYNOP, -1, VN_U,  99.99)
!type (t_obserr) ,save ,target :: ifs_obserr (16) = IFS_99
!#else

  type (t_obserr) ,save ,target :: ifs_obserr (16) = (/ &
  t_obserr(OT_SYNOP ,-1, VN_U  ,&
  (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.4_wp,  3.6_wp, &
     3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
     2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /),&
     0._wp),&
  t_obserr(OT_AIREP ,-1, VN_U  ,&
  (/ 2.5_wp,  2.5_wp,  3.0_wp,  3.5_wp,  4.0_wp, &
     4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp, &
     4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp /),&
     0._wp),&
  t_obserr(OT_SATOB ,-1, VN_U  ,&
  (/ 2.0_wp,  2.0_wp,  2.0_wp,  3.5_wp,  4.3_wp, &
     5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp, &
     5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp,  5.7_wp /),&
     0._wp),&
  t_obserr(OT_DRIBU ,-1, VN_U  ,&
  (/ 2.4_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
     3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
     2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /),&
     0._wp),&
  t_obserr(OT_TEMP  ,-1, VN_U  ,&
  (/ 2.3_wp,  2.3_wp,  2.5_wp,  3.0_wp,  3.5_wp, &
     3.7_wp,  3.5_wp,  3.5_wp,  3.4_wp,  3.3_wp, &
     3.2_wp,  3.2_wp,  3.3_wp,  3.6_wp,  4.5_wp /),&
     0._wp),&
  t_obserr(OT_PILOT ,-1, VN_U  ,&
  (/ 2.3_wp,  2.3_wp,  2.5_wp,  3.0_wp,  3.5_wp, &
     3.7_wp,  3.5_wp,  3.5_wp,  3.4_wp,  3.3_wp, &
     3.2_wp,  3.2_wp,  3.3_wp,  3.6_wp,  4.5_wp /),&
     0._wp),&
  t_obserr(OT_SCATT ,-1, VN_U  ,&
  (/ 2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp, &
     2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp, &
     2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp /),&
     0._wp),&
  t_obserr(OT_SYNOP ,-1, VN_Z ,&
  (/ 7.0_wp,  8.0_wp,  8.6_wp, 12.1_wp, 14.9_wp, &
    18.8_wp, 25.4_wp, 27.7_wp, 32.4_wp, 39.4_wp, &
    50.3_wp, 59.3_wp, 69.8_wp, 96.0_wp,114.2_wp /),&
     0._wp),&
  t_obserr(OT_DRIBU ,-1, VN_Z ,&
  (/11.5_wp, 11.5_wp, 17.2_wp, 24.2_wp, 29.8_wp, &
    37.6_wp, 50.8_wp, 55.4_wp, 64.8_wp, 78.8_wp, &
   100.6_wp,118.6_wp,139.6_wp,192.0_wp,228.4_wp /),&
     0._wp),&
  t_obserr(OT_TEMP  ,-1, VN_Z  ,&
  (/ 4.3_wp,  4.4_wp,  5.2_wp,  8.4_wp,  9.8_wp, &
    10.7_wp, 11.8_wp, 13.2_wp, 15.2_wp, 18.1_wp, &
    19.5_wp, 22.5_wp, 25.0_wp, 32.0_wp, 40.0_wp /),&
     0._wp),&
  t_obserr(OT_PILOT ,-1, VN_Z  ,&
  (/ 4.3_wp,  4.4_wp,  5.2_wp,  8.4_wp,  9.8_wp, &
    10.7_wp, 11.8_wp, 13.2_wp, 15.2_wp, 18.1_wp, &
    19.5_wp, 22.5_wp, 25.0_wp, 32.0_wp, 40.0_wp /),&
     0._wp),&
  t_obserr(OT_PAOB  ,-1, VN_Z ,&
  (/24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, &
    24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, &
    24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp /),&
     0._wp),&
  t_obserr(OT_SYNOP ,-1, VN_T ,&
  (/ 2.0_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.3_wp, &
     1.5_wp,  1.8_wp,  1.8_wp,  1.9_wp,  2.0_wp, &
     2.2_wp,  2.4_wp,  2.5_wp,  2.5_wp,  2.5_wp /),&
     0._wp),&
  t_obserr(OT_AIREP ,-1, VN_T ,&
  (/ 1.4_wp,  1.3_wp,  1.2_wp,  1.2_wp,  1.2_wp, &
     1.3_wp,  1.3_wp,  1.4_wp,  1.4_wp,  1.4_wp, &
     1.5_wp,  1.6_wp,  1.8_wp,  2.0_wp,  2.2_wp /),&
     0._wp),&
  t_obserr(OT_DRIBU ,-1, VN_T ,&
  (/ 1.8_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.3_wp, &
     1.5_wp,  1.8_wp,  1.8_wp,  1.9_wp,  2.0_wp, &
     2.2_wp,  2.4_wp,  2.5_wp,  2.5_wp,  2.5_wp /),&
     0._wp),&
  t_obserr(OT_TEMP  ,-1, VN_T ,&
  (/ 1.7_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.2_wp, &
     1.4_wp,  1.5_wp,  1.5_wp,  1.6_wp,  1.7_wp, &
     1.8_wp,  1.9_wp,  2.0_wp,  2.2_wp,  2.5_wp /),&
     0._wp)/)

!#endif

  !----------------------
  ! OI observation errors
  !----------------------
!#if defined (__SX__)
! type (t_obserr) ,save ,target :: dwd_obserr (6) = IFS_99
!#else
  type (t_obserr) ,save ,target :: dwd_obserr (6) = (/ &
  t_obserr(OT_PILOT ,-1, VN_U  , &
           (/ 2.0_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
              3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
              2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /),&
            0._wp),&
  t_obserr(OT_AIREP ,-1, VN_U  , &
           (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.0_wp,  3.5_wp, &
              4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp, &
              4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp /),&
            0._wp),&
  t_obserr(OT_AIREP ,-1, VN_RH , &
           (/ .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp, &
              .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp, &
              .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp /),&
            0._wp),&
  t_obserr(OT_SATOB ,-1, VN_U  , &
           (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.0_wp,  6.0_wp, &
              6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp, &
              6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp /),&
            0._wp),&
  t_obserr(OT_TEMP,  -1, VN_U  , &
           (/ 2.0_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
              3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
              2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /),&
            0._wp),&
  t_obserr(OT_TEMP  ,-1, VN_Z  , &
           (/ 6.5_wp,  6.6_wp,  7.8_wp, 12.6_wp, 14.7_wp, &
             16.1_wp, 17.7_wp, 19.8_wp, 22.8_wp, 27.2_wp, &
             29.3_wp, 33.8_wp, 37.5_wp, 48.0_wp, 60.0_wp /),&
            0._wp) /)
!#endif

  !-----------------------
  ! used observation errors
  !-----------------------
  type (t_obserr) ,save, target :: use_obserr (ntab)
  type (t_obserr) ,save, target :: adaptloc   (ntab)
  integer :: n_use = 0
  integer :: n_adp = 0

  !----------------------------------------
  ! Temp categories (derived from location)
  !----------------------------------------
  type t_cat
    real(wp) :: lat1
    real(wp) :: lat2
    real(wp) :: lon1
    real(wp) :: lon2
    real(wp) :: f
  end type t_cat

  type (t_cat) ,parameter :: cat (4) = &
     !         lat1    lat2   lon1    lon2    factor
     (/ t_cat( 30._wp, 80._wp,190._wp,310._wp,0.66_wp),& ! cat1 North America
        t_cat(-60._wp, 10._wp,270._wp,330._wp,1.33_wp),& ! cat3 South America
        t_cat(-40._wp, 30._wp,340._wp,360._wp,1.33_wp),& ! cat3 Africa,Arabia
        t_cat(-40._wp, 30._wp,  0._wp, 60._wp,1.33_wp)/) ! cat3 Africa,Arabia
  !----------------
  ! local variables
  !----------------
  logical  :: not_init = .true. ! true if module is not yet initialized

contains
!==============================================================================
  subroutine init_obs_err

!#if defined (__SX__)
!!------------------------------------------------------------
!! default initialisation didnt work on NEC SX: set explicitly
!!------------------------------------------------------------
!ifs_obserr(01)=t_obserr(OT_SYNOP,-1, VN_U,&
!  (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.4_wp,  3.6_wp, &
!     3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
!     2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /))
!ifs_obserr(02)=t_obserr(OT_AIREP,-1, VN_U,&
!  (/ 2.5_wp,  2.5_wp,  3.0_wp,  3.5_wp,  4.0_wp, &
!     4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp, &
!     4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp /))
!ifs_obserr(03)=t_obserr(OT_SATOB,-1, VN_U,&
!  (/ 2.0_wp,  2.0_wp,  2.0_wp,  3.5_wp,  4.3_wp, &
!     5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp, &
!     5.0_wp,  5.0_wp,  5.0_wp,  5.0_wp,  5.7_wp /))
!ifs_obserr(04)=t_obserr(OT_DRIBU,-1, VN_U,&
!  (/ 2.4_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
!     3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
!     2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /))
!ifs_obserr(05)=t_obserr(OT_TEMP, -1, VN_U,&
!  (/ 2.3_wp,  2.3_wp,  2.5_wp,  3.0_wp,  3.5_wp, &
!     3.7_wp,  3.5_wp,  3.5_wp,  3.4_wp,  3.3_wp, &
!     3.2_wp,  3.2_wp,  3.3_wp,  3.6_wp,  4.5_wp /))
!ifs_obserr(06)=t_obserr(OT_PILOT,-1, VN_U,&
!  (/ 2.3_wp,  2.3_wp,  2.5_wp,  3.0_wp,  3.5_wp, &
!     3.7_wp,  3.5_wp,  3.5_wp,  3.4_wp,  3.3_wp, &
!     3.2_wp,  3.2_wp,  3.3_wp,  3.6_wp,  4.5_wp /))
!ifs_obserr(07)=t_obserr(OT_SCATT,-1, VN_U,&
!  (/ 2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp, &
!     2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp, &
!     2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp,  2.0_wp /))
!ifs_obserr(08)=t_obserr(OT_SYNOP,-1, VN_Z,&
!  (/ 7.0_wp,  8.0_wp,  8.6_wp, 12.1_wp, 14.9_wp, &
!    18.8_wp, 25.4_wp, 27.7_wp, 32.4_wp, 39.4_wp, &
!    50.3_wp, 59.3_wp, 69.8_wp, 96.0_wp,114.2_wp /))
!ifs_obserr(09)=t_obserr(OT_DRIBU,-1, VN_Z,&
!  (/11.5_wp, 11.5_wp, 17.2_wp, 24.2_wp, 29.8_wp, &
!    37.6_wp, 50.8_wp, 55.4_wp, 64.8_wp, 78.8_wp, &
!   100.6_wp,118.6_wp,139.6_wp,192.0_wp,228.4_wp /))
!ifs_obserr(10)=t_obserr(OT_TEMP,-1, VN_Z,&
!  (/ 4.3_wp,  4.4_wp,  5.2_wp,  8.4_wp,  9.8_wp, &
!    10.7_wp, 11.8_wp, 13.2_wp, 15.2_wp, 18.1_wp, &
!    19.5_wp, 22.5_wp, 25.0_wp, 32.0_wp, 40.0_wp /))
!ifs_obserr(11)=t_obserr(OT_PILOT,-1, VN_Z,&
!  (/ 4.3_wp,  4.4_wp,  5.2_wp,  8.4_wp,  9.8_wp, &
!    10.7_wp, 11.8_wp, 13.2_wp, 15.2_wp, 18.1_wp, &
!    19.5_wp, 22.5_wp, 25.0_wp, 32.0_wp, 40.0_wp /))
!ifs_obserr(12)=t_obserr(OT_PAOB,-1, VN_Z,&
!  (/24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, &
!    24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, &
!    24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp, 24.0_wp /))
!ifs_obserr(13)=t_obserr(OT_SYNOP,-1, VN_T,&
!  (/ 2.0_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.3_wp, &
!     1.5_wp,  1.8_wp,  1.8_wp,  1.9_wp,  2.0_wp, &
!     2.2_wp,  2.4_wp,  2.5_wp,  2.5_wp,  2.5_wp /))
!ifs_obserr(14)=t_obserr(OT_AIREP,-1, VN_T,&
!  (/ 1.4_wp,  1.3_wp,  1.2_wp,  1.2_wp,  1.2_wp, &
!     1.3_wp,  1.3_wp,  1.4_wp,  1.4_wp,  1.4_wp, &
!     1.5_wp,  1.6_wp,  1.8_wp,  2.0_wp,  2.2_wp /))
!ifs_obserr(15)=t_obserr(OT_DRIBU,-1, VN_T,&
!  (/ 1.8_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.3_wp, &
!     1.5_wp,  1.8_wp,  1.8_wp,  1.9_wp,  2.0_wp, &
!     2.2_wp,  2.4_wp,  2.5_wp,  2.5_wp,  2.5_wp /))
!ifs_obserr(16)=t_obserr(OT_TEMP, -1, VN_T,&
!  (/ 1.7_wp,  1.5_wp,  1.3_wp,  1.2_wp,  1.2_wp, &
!     1.4_wp,  1.5_wp,  1.5_wp,  1.6_wp,  1.7_wp, &
!     1.8_wp,  1.9_wp,  2.0_wp,  2.2_wp,  2.5_wp /))
!
!dwd_obserr(01)=t_obserr(OT_PILOT,-1, VN_U  , &
!           (/ 2.0_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
!              3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
!              2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /))
!dwd_obserr(02)=t_obserr(OT_AIREP,-1, VN_U  , &
!           (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.0_wp,  3.5_wp, &
!              4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp, &
!              4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp,  4.0_wp /))
!dwd_obserr(03)=t_obserr(OT_AIREP,-1, VN_RH , &
!           (/ .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp, &
!              .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp, &
!              .15_wp,  .15_wp,  .15_wp,  .15_wp,  .15_wp /))
!dwd_obserr(04)=t_obserr(OT_SATOB,-1, VN_U  , &
!           (/ 3.0_wp,  3.0_wp,  3.0_wp,  3.0_wp,  6.0_wp, &
!              6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp, &
!              6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp,  6.0_wp /))
!dwd_obserr(05)=t_obserr(OT_TEMP, -1, VN_U  , &
!           (/ 2.0_wp,  2.4_wp,  2.5_wp,  3.4_wp,  3.6_wp, &
!              3.8_wp,  3.2_wp,  3.2_wp,  2.4_wp,  2.2_wp, &
!              2.0_wp,  2.0_wp,  2.0_wp,  2.5_wp,  3.0_wp /))
!dwd_obserr(06)=t_obserr(OT_TEMP ,-1, VN_Z  , &
!           (/ 6.5_wp,  6.6_wp,  7.8_wp, 12.6_wp, 14.7_wp, &
!             16.1_wp, 17.7_wp, 19.8_wp, 22.8_wp, 27.2_wp, &
!             29.3_wp, 33.8_wp, 37.5_wp, 48.0_wp, 60.0_wp /))
!
!#endif

    !-------------------------
    ! preprocess lookup tables
    !-------------------------
    if (.not. not_init) return
    not_init = .false.
    call read_nml
  end subroutine init_obs_err
!------------------------------------------------------------------------------
  subroutine read_nml
    integer                  :: ierr, i, k, n
    logical                  :: first
    integer                  :: obst            ! observation type
    integer                  :: qnt             ! quantity (varno)
    type (t_obserr) ,pointer :: src  (:)        ! pointer to source table
    type (t_obserr) ,pointer :: dest (:)        ! pointer to destination
    type (t_obserr) ,target  :: extern(1)       ! external table entry
    logical                  :: prhead = .true.
#if defined(__ibm__)
    integer                  :: ios
#endif

    !------------------------------------
    ! make sure feedback tables are valid
    !------------------------------------
    call init_fdbk_tables      ! initialize tables
    !-------------
    ! set defaults
    !-------------
    a   = 0.8_wp
    b   = 1.0_wp
    !--------------
    ! read namelist
    !--------------
    first = .true.
    do
      if (dace% lpio) then
        obstype  = ''
        codetype = -1
        quantity = ''
        table    = ''
        scale    = -1._wp
        scale(1) =  1._wp
        err      = -HUGE (0._wp)
        r_err    = -HUGE (0._wp)
        call position_nml ('OBSERR', lrewind=first, status=ierr)
        first = .false.
        select case (ierr)
        case (POSITIONED)
#if defined(__ibm__)
        read (nnml ,nml=OBSERR ,iostat=ios)
        if (ios/=0) call finish ('read_nml','ERROR in namelist /OBSERR/')
#else
        read (nnml ,nml=OBSERR)
#endif
        end select
      endif
      !-------------------------------------------
      ! exit if no further namelist group is found
      !-------------------------------------------
      call p_bcast (ierr, dace% pio)
      if (ierr /= POSITIONED) exit
      !-----------------------------
      ! broadcast namelist variables
      !-----------------------------
      call p_bcast (obstype  ,dace% pio)
      call p_bcast (codetype ,dace% pio)
      call p_bcast (quantity ,dace% pio)
      call p_bcast (table    ,dace% pio)
      call p_bcast (scale    ,dace% pio)
      call p_bcast (err      ,dace% pio)
      call p_bcast (r_err    ,dace% pio)
      do i=2,size(scale)
        if (scale(i)<0._wp) scale(i) = scale(i-1)
      end do
      !---------------------------------------------
      ! if 'obstype' is given, set up table for the
      ! vertical dependency of the observation error
      !---------------------------------------------
      if (obstype/='') then
        !-----------------------------------------
        ! select pre-defined tables (IFS/OI/3DVAR)
        ! or externally specified profile
        !-----------------------------------------
        select case (table)
        case ('oi'    ,'OI'    )
          dest => use_obserr
          src  => dwd_obserr
        case ('ifs'   ,'IFS'   )
          dest => use_obserr
          src  => ifs_obserr
        case ('extern','EXTERN')
          dest => use_obserr
          src  => extern
        case ('adaptloc','ADAPTLOC')
          dest => adaptloc
          src  => extern
        case default
          call finish('read_nml(mo_obs_err)','invalid table = '//table)
        end select
        !---------------------------------
        ! search a free index in the table
        !---------------------------------
        k = 0
        do i=1,size(dest)
          if (dest(i)% obstype == 0) then
            k = i
            exit
          endif
        end do
        if (k==0) call finish('read_nml(mo_obs_err)','dest too small')
        !------------------------------------------------------------
        ! convert observation type from character string to index key
        !------------------------------------------------------------
        obst = 0
        obstype = toupper(obstype)
        do i=1,size(obstyp)
          if (obstyp(i)% name == obstype) then
            obst = obstyp(i)% key
            exit
          endif
        end do
        if (obst == 0) &
          call finish('read_nml(mo_obs_err)','invalid obstype: '//obstype)
        !--------------------------------------
        ! convert quantity to index key (varno)
        !--------------------------------------
        qnt = 0
        if (quantity/='') qnt = value_name (varno, toupper (quantity))
        if (qnt == 0) call finish ("read_nml(mo_obs_err)",   &
                                   "invalid varno "//quantity)
        !-----------------------------------
        ! prepare externally set table entry
        !-----------------------------------
        if (associated (src, extern)) then
          if (qnt == 0) then
             call finish ("read_nml(mo_obs_err)", &
                          "must specify quantity for external table entry!")
          end if
          r_err     = max (r_err, 0._wp)
          extern(1) = t_obserr (obst, codetype, qnt, err, r_err)
        end if
        !-----------------------------------------
        ! copy specified entries from source table
        !-----------------------------------------
        n=0
        do i=1,size(src)
          if (obst == src(i)% obstype .and.       &
              (qnt == src(i)% varno   .or. qnt==0)) then
            dest(k)      = src(i)
            dest(k)% err = src(i)% err * scale
            if (dace% lpio) then
              if (prhead) then
                 write(6,'()')
                 write(6,'(a)')      ' Namelist /OBSERR/:'
                 write(6,'()')
                 prhead = .false.
              end if
              write(6,'(5A)')        ' observation error ', obstype,      &
                                       name_value (varno, src(i)% varno), &
                                     '  taken from ', trim (table)
              write(6,'(A,20I5   )') '  codetype =', pack (dest(k)% codetype,  &
                                                           dest(k)% codetype>=0)
              write(6,'(A,20F10.6)') '     scale =', scale
              write(6,'(A,20F10.6)') '       err =', dest(k)% err
              write(6,'(A,20F10.6)') '     r_err =', dest(k)% r_err
              write(6,'()')
            end if
            k=k+1
            n=n+1
            if(k>size(dest)) &
              call finish('read_nml(mo_obs_err)','dest too small')
          endif
        end do
        !----------------------------
        ! abort if no entry was found
        !----------------------------
        if (n == 0) then
           if (dace% lpio) then
              write(6,'(5A)') "warning: observation error ", obstype, &
                   quantity, " not found in table ", trim (table)
              write(6,'()')
              write(0,'(5A)') "warning: observation error ", obstype, &
                   quantity, " not found in table ", trim (table)
              write(0,'()')
           end if
        end if
        nullify (src)
      endif
    end do
    !-----------------------------
    ! broadcast namelist variables
    !-----------------------------
    call p_bcast (a                   ,dace% pio)
    call p_bcast (b                   ,dace% pio)
    call p_bcast (obserr_synopship_h  ,dace% pio)
    call p_bcast (obserr_synopland_h  ,dace% pio)
    call p_bcast (obserr_synopauto_h  ,dace% pio)
    call p_bcast (obserr_coastal_h    ,dace% pio)
    call p_bcast (obserr_paob_h       ,dace% pio)
    call p_bcast (obserr_metar_h      ,dace% pio)
    call p_bcast (obserr_dribu_h      ,dace% pio)
    call p_bcast (obserr_h_reduction  ,dace% pio)
    call p_bcast (obserr_synopship_rh ,dace% pio)
    call p_bcast (obserr_synopland_rh ,dace% pio)
    call p_bcast (obserr_synophigh_rh ,dace% pio)
    call p_bcast (obserr_metar_rh     ,dace% pio)
    call p_bcast (obserr_dribu_rh     ,dace% pio)
    call p_bcast (obserr_synopland_v  ,dace% pio)
    call p_bcast (obserr_synopship_v  ,dace% pio)
    call p_bcast (obserr_metar_v      ,dace% pio)
    call p_bcast (obserr_dribu_v      ,dace% pio)
    call p_bcast (obserr_synopland_t  ,dace% pio)
    call p_bcast (obserr_synopland_dt ,dace% pio)
    call p_bcast (obserr_synopship_t  ,dace% pio)
    call p_bcast (obserr_metar_t      ,dace% pio)
    call p_bcast (obserr_dribu_t      ,dace% pio)
    call p_bcast (obserr_temp_rh      ,dace% pio)
    call p_bcast (obserr_temp_rh_dry  ,dace% pio)
    call p_bcast (obserr_temp_rh_cold ,dace% pio)
    !-------------------
    ! derived quantities
    !-------------------
    n_use = count (use_obserr% obstype /= 0)
    n_adp = count (adaptloc  % obstype /= 0)
    if (obserr_synopauto_h  < 0._wp) obserr_synopauto_h  = obserr_synopland_h
    if (obserr_synophigh_rh < 0._wp) obserr_synophigh_rh = obserr_synopland_rh
    !-----------------------------------------------
    ! Print parameters of surface observation errors
    !-----------------------------------------------
    if (dace% lpio) then
       if (prhead) then
          write(6,'()')
          write(6,'(a)')   ' Namelist /OBSERR/:'
          write(6,'()')
          prhead = .false.
       end if
       write(6,'(/A/)')     ' Surface observations:'
       write(6,'(A,F8.3)')  '  obserr_h_reduction  =',obserr_h_reduction
       write(6,'(A,F8.3,A)')'  obserr_synopship_h  =',obserr_synopship_h," gpm"
       write(6,'(A,F8.3,A)')'  obserr_synopland_h  =',obserr_synopland_h," gpm"
       write(6,'(A,F8.3,A)')'  obserr_synopauto_h  =',obserr_synopauto_h," gpm"
       write(6,'(A,F8.3,A)')'  obserr_coastal_h    =',obserr_coastal_h  ," gpm"
       write(6,'(A,F8.3,A)')'  obserr_metar_h      =',obserr_metar_h    ," gpm"
       write(6,'(A,F8.3,A)')'  obserr_paob_h       =',obserr_paob_h     ," gpm"
       write(6,'(A,F8.3,A)')'  obserr_dribu_h      =',obserr_dribu_h    ," gpm"
       write(6,'(A,F8.3)')  '  obserr_synopship_rh =',obserr_synopship_rh
       write(6,'(A,F8.3)')  '  obserr_synopland_rh =',obserr_synopland_rh
       write(6,'(A,F8.3)')  '  obserr_synophigh_rh =',obserr_synophigh_rh
       write(6,'(A,F8.3)')  '  obserr_metar_rh     =',obserr_metar_rh
       write(6,'(A,F8.3)')  '  obserr_dribu_rh     =',obserr_dribu_rh
       write(6,'(A,F8.3,A)')'  obserr_synopship_v  =',obserr_synopship_v," m/s"
       write(6,'(A,F8.3,A)')'  obserr_synopland_v  =',obserr_synopland_v," m/s"
       write(6,'(A,F8.3,A)')'  obserr_metar_v      =',obserr_metar_v    ," m/s"
       write(6,'(A,F8.3,A)')'  obserr_dribu_v      =',obserr_dribu_v    ," m/s"
       write(6,'(A,F8.3,A)')'  obserr_synopland_t  =',obserr_synopland_t, " K"
       write(6,'(A,F8.3,A)')'  obserr_synopland_dt =',obserr_synopland_dt," K/m"
       write(6,'(A,F8.3,A)')'  obserr_synopship_t  =',obserr_synopship_t," K"
       write(6,'(A,F8.3,A)')'  obserr_metar_t      =',obserr_metar_t    ," K"
       write(6,'(A,F8.3,A)')'  obserr_dribu_t      =',obserr_dribu_t    ," K"
       write(6,'(/A/)')     ' TEMP humidity observations:'
       write(6,'(A,F8.3)')  '  obserr_temp_rh      =',obserr_temp_rh
       write(6,'(A,F8.3)')  '  obserr_temp_rh_dry  =',obserr_temp_rh_dry
       write(6,'(A,F8.3)')  '  obserr_temp_rh_cold =',obserr_temp_rh_cold
       write(6,'()')
    end if
    !-----------------------------------------
    ! Read station specific observation errors
    !-----------------------------------------
    call read_stn_obserr ()
  end subroutine read_nml
!------------------------------------------------------------------------------
  subroutine synop_obs_err (obs, spot, it)
  type(t_obs_block),intent(inout) :: obs    ! observation data type
  type(t_spot)     ,intent(in)    :: spot   ! SPOT observations
  integer          ,intent(in)    :: it (:) ! observation data type
    integer  :: i,k
    real(wp) :: dh2
    real(wp) :: ff, eff
    real(wp) :: dh, err, var    ! height-difference, adjusted error, variance
    integer  :: iii  !+++ workaround for NEC SX6 bug
    integer, parameter :: me = 4     ! max.# if station obserr entries
    integer            :: nentry     ! number of entries found
    type(t_stn_obserr) :: entry(me)  ! station obserr entries
    real(wp)           :: std        ! error standard deviation

    !--------------------------------------------------------------
    ! Safeguard: set default null values for dismissed observations
    !--------------------------------------------------------------
    if (spot% use% state <= STAT_DISMISS) then
       call zero_obs_err (obs, spot)
       return
    end if

    k = obs% R% ia (spot% o% i+1)
    select case (spot% hd% obstype)
    case (OT_SYNOP)
      nentry = obserr_entry (spot% statid, spot% wsi, spot% hd% codetype, entry)
      ff = 0._wp; eff = 0._wp
      do i=1,spot% o% n
        iii = spot% o% i + i   !+++ workaround for NEC SX6 bug
        obs% R% ia (iii) = k
        select case (it(i))

        case (VN_U, VN_V, VN_U10M, VN_V10M, VN_FF)
          if (it(i)==VN_FF) ff = obs% o% body(iii)% o
          select case (spot% hd% codetype)
          case (11,14)  ! Manual, automatic land station
             eff = obserr_synopland_v
          case (17)     ! Metar USA (Synop Bufr)
             eff = obserr_synopland_v
          case (15)     ! SWIS (Road weather stations)
             eff = obserr_synopland_v
          case (20)     ! Coastal marine automatic station
             eff = obserr_synopland_v
          case (21:24)  ! SHIP
             eff = obserr_synopship_v
          case (140)    ! METAR
             eff = obserr_metar_v
          case default
             call finish('synop_obs_err','SYNOP: invalid codetype, wind obsv.')
          end select
          !----------------------------------------
          ! Use station specific observation error?
          !----------------------------------------
          if (nentry > 0) then
             std = get_err (VN_U)
             if (std > 0._wp) eff = std
          end if
          var = eff ** 2

        case (VN_DD)
          if (ff > 0._wp) then
            eff = r2d * eff / ff
          else
            eff = 360._wp
          endif
          eff = min (eff, 360._wp)
          var = eff ** 2

        case (VN_Z)
          select case (spot% hd% dbkz)
          case default
            dh2 = (obs% o% body(iii)% o - spot% gp_bg / gacc) ** 2 &
                + (obs% o% body(iii)% o - spot% z           ) ** 2
            var = obserr_synopland_h ** 2 &     ! SYNOP height error
                + obserr_h_reduction ** 2 * dh2
          case (128, 10128, 10143, 10150)  ! SYNOP automatic
            dh2 = (obs% o% body(iii)% o - spot% gp_bg / gacc) ** 2 &
                + (obs% o% body(iii)% o - spot% z           ) ** 2
            var = obserr_synopauto_h ** 2 &     ! SYNOP height error
                + obserr_h_reduction ** 2 * dh2
          case (10158)                  ! CMAN stations
            dh2 = (obs% o% body(iii)% o - spot% gp_bg / gacc) ** 2 &
                + (obs% o% body(iii)% o - spot% z           ) ** 2
            var = obserr_coastal_h   ** 2 &     ! CMAN height error
                + obserr_h_reduction ** 2 * dh2
          case (256, 10256)             ! manual SHIP
            var = obserr_synopship_h**2 ! manualSHIP height error
          case (1)                      ! METAR
            var = obserr_metar_h    **2 ! METAR      height error
          end select
          !----------------------------------------
          ! Use station specific observation error?
          !----------------------------------------
          if (nentry > 0) then
             std = get_err (VN_Z)
             if (std > 0._wp) var = std ** 2
          end if

        case (VN_PS, VN_P)
          select case (spot% hd% dbkz)
          case default
            dh2 = (obs% o% olev (iii) - spot% gp_bg / gacc) ** 2 &
                + (obs% o% olev (iii) - spot% z           ) ** 2
            var = ( obserr_synopland_h ** 2 &      ! SYNOP height error
                  + obserr_h_reduction ** 2 * dh2) * spot% pz_bg ** 2
          case (128, 10128, 10143, 10150)   ! SYNOP automatic
            dh2 = (obs% o% olev (iii) - spot% gp_bg / gacc) ** 2 &
                + (obs% o% olev (iii) - spot% z           ) ** 2
            var = ( obserr_synopauto_h ** 2 &      ! SYNOP height error
                  + obserr_h_reduction ** 2 * dh2) * spot% pz_bg ** 2
          case (10158)                  ! CMAN stations
            dh2 = (obs% o% olev (iii) - spot% gp_bg / gacc) ** 2 &
                + (obs% o% olev (iii) - spot% z           ) ** 2
            var = ( obserr_coastal_h   ** 2 &      ! CMAN height error
                  + obserr_h_reduction ** 2 * dh2) * spot% pz_bg ** 2
          case (256, 10256)             ! manual SHIP
            var = (obserr_synopship_h * spot% pz_bg) **2 ! manualSHIP height error
          case (1)                      ! METAR
            var = (obserr_metar_h     * spot% pz_bg) **2 ! METAR      height error
          end select
          !----------------------------------------
          ! Use station specific observation error?
          !----------------------------------------
          if (nentry > 0) then
             std = get_err (VN_Z)
             if (std > 0._wp) var = (std * spot% pz_bg) ** 2
          end if

        case (VN_RH, VN_RH2M)
          select case (spot% hd% dbkz)
          case default                  ! SYNOP
            var = ( obserr_synopland_rh * (1._wp - spot% z/5000._wp)     &
                +   obserr_synophigh_rh * (        spot% z/5000._wp)) ** 2
          case (256, 384, 10256, 10384) ! SHIP
            var = obserr_synopship_rh ** 2
          case (1)                      ! METAR
            var = obserr_metar_rh     ** 2
          end select
          !----------------------------------------
          ! Use station specific observation error?
          !----------------------------------------
          if (nentry > 0) then
             std = get_err (VN_RH)
             if (std > 0._wp) var = std ** 2
          end if

        case (VN_T2M, VN_TD2M, VN_T)
          !--------------------------------------------
          ! dewpoint temperature is used only passively
          ! currently dummy error (same as temperature)
          !--------------------------------------------
          dh = abs (spot% gp_bg / gacc - spot% z)
          select case (spot% hd% codetype)
          case (11,14)  ! Manual, automatic land station
             err = obserr_synopland_t + dh * obserr_synopland_dt
          case (17)     ! Metar USA (Synop Bufr)
             err = obserr_synopland_t + dh * obserr_synopland_dt
          case (15)     ! SWIS (Road weather stations)
             err = obserr_synopland_t + dh * obserr_synopland_dt
          case (16)     ! SNOW (SYNOP)
             err = obserr_synopland_t + dh * obserr_synopland_dt
          case (21:24)  ! SHIP
             err = obserr_synopship_t
          case (140)    ! METAR
             err = obserr_metar_t     + dh * obserr_synopland_dt
          case (20)     ! CMAN has sensors at different height (not 2m!)
             select case (it(i))
             case (VN_T)
                err = obserr_synopland_t + dh * obserr_synopland_dt
             case default
                call finish('synop_obs_err','CMAN: invalid t2m/td2m obsv.')
             end select
          case default
             call finish('synop_obs_err','SYNOP: invalid codetype, t2m obsv.')
          end select
          !----------------------------------------
          ! Use station specific observation error?
          !----------------------------------------
          if (nentry > 0) then
             std = get_err (VN_T)
             if (std > 0._wp) err = std
          end if
          var = err ** 2

        case (VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND, &
              VN_N, VN_NH, VN_N_L, VN_N_M, VN_N_H, VN_VV, &
              VN_CEIL, VN_RAD_GL, VN_RAD_DF, VN_RAD_LW,   &
              VN_SDEPTH, VN_TSEA, VN_WW, VN_GCLG, VN_ICLG )
          var = 0._wp
        case default
          call finish('synop_obs_err','SYNOP: invalid varno: '//char3(it(i)))
        end select
        obs% R% ja    (k)     = iii
        obs% R% packed(k)     = var
!       obs% o% body(iii)% eo = sqrt (obs% R% packed(k))
        k = k + 1
      end do
    case (OT_DRIBU)
      ff = 0._wp; eff = 0._wp
      do i=1,spot% o% n
        iii = spot% o% i + i   !+++ workaround for NEC SX6 bug
        obs% R% ia (iii) = k
        select case (it(i))

        case (VN_U, VN_V, VN_U10M, VN_V10M, VN_FF)
          if (it(i)==VN_FF) ff = obs% o% body(iii)% o
          select case (spot% hd% dbkz)
          case (1697, 1698, 1699) ! QSCAT ASCAT
            eff = 2.0_wp          ! SCATT wind error (as IFS)
          case default
            eff = obserr_dribu_v  ! DRIBU wind error
          end select
          var = eff ** 2

        case (VN_DD)
          if (ff > 0._wp) then
            eff = r2d * eff / ff
          else
            eff = 360._wp
          endif
          eff = min (eff, 360._wp)
          var = eff ** 2

        case (VN_Z)
          var = obserr_dribu_h **2    ! DRIBU height error
        case (VN_PS)
          var = (obserr_dribu_h * spot% pz_bg) **2
        case (VN_RH, VN_RH2M)
          var = obserr_dribu_rh **2   ! DRIBU rel.hum. error
        case (VN_T2M, VN_TD2M)
          !--------------------------------------------
          ! dewpoint temperature is used only passively
          ! currently dummy error (same as temperature)
          !--------------------------------------------
          select case (spot% hd% codetype)
          case (165)                        ! DRIBU
             var = obserr_dribu_t **2
          case default
          write(0,*)' invalid code type:', spot% hd% codetype
             call finish('synop_obs_err','DRIBU: invalid codetype, t2m obsv.')
          end select
        case (VN_GUST, VN_RR, VN_TMIN, VN_TMAX, VN_PTEND, &
              VN_N, VN_NH, VN_N_L, VN_N_M, VN_N_H, VN_VV, &
              VN_SDEPTH, VN_TSEA                          )
          var = 0._wp
        case default
          write(0,*)' invalid varno:',it(i)
          call finish('synop_obs_err','DRIBU: invalid varno')
        end select
        obs% R% ja    (k)     = iii
        obs% R% packed(k)     = var
!       obs% o% body(iii)% eo = sqrt (obs% R% packed(k))
        k = k + 1
      end do
    case (OT_PAOB)
      do i=1,spot% o% n
        iii = spot% o% i + i   !+++ workaround for NEC SX6 bug
        obs% R% ia (iii) = k
        select case (it(i))
        case (VN_Z)
          obs% R% packed(k) = obserr_paob_h **2  ! PAOB height error
        case (VN_PS)
          obs% R% packed(k) = (obserr_paob_h * spot% pz_bg) **2
        case default
          write(0,*)' invalid varno:',it(i)
          call finish('synop_obs_err','PAOB: invalid varno')
        end select
        obs% R% ja (k) = iii
        k = k + 1
      end do
    case default
      call finish('synop_obs_err','invalid observation type ')
    end select
    obs% R% ia (spot% o% i + spot% o% n + 1) = k

  contains

    function get_err (varno) result(err)
      integer, intent(in)  :: varno
      real(wp)             :: err
      integer :: i
      err = 0._wp
      do i = 1, nentry
         if (entry(i)% varno == varno) then
            err = entry(i)% err
            exit
         end if
      end do
    end function get_err

  end subroutine synop_obs_err
!------------------------------------------------------------------------------
  subroutine scatt_obs_err (obs, spot, it)
    type(t_obs_block),intent(inout) :: obs    ! observation data type
    type(t_spot)     ,intent(in)    :: spot   ! SPOT observations
    integer          ,intent(in)    :: it (:) ! observation data type

    integer  :: i,k
    real(wp) :: ff, eff, var
    integer  :: iii  !+++ workaround for NEC SX6 bug

    k = obs% R% ia (spot% o% i+1)
    select case (spot% hd% obstype)
    case (OT_SCATT)
      ff = 0._wp; eff = 0._wp
      do i=1,spot% o% n
        iii = spot% o% i + i   !+++ workaround for NEC SX6 bug
        obs% R% ia (iii) = k
        select case (it(i))

        case (VN_U, VN_U10M, VN_V, VN_V10M, VN_FF)
         if (it(i)==VN_FF) ff = obs% o% body(iii)% o
!        select case (spot% hd% dbkz)
!          case (1701,1702)
!            eff = 1.5_wp   ! JASON, SARAL wind error
!          case (1699,1700)
!            eff = 1.5_wp   ! ASCAT, OSCAT wind error
!          case (1697)
!            eff = 2.0_wp   ! QSCAT wind error
!          case default
             !---------------------------------------------------------
             ! no DBKZ given, probably third party report from fof-file
             ! keep observation error if already set
             !---------------------------------------------------------
!            if (obs% o% body(iii)% eo > 0._sp) then
!              eff = obs% o% body(iii)% eo
!            else
!              write(0,*)'SCATT: invalid KZ type',spot% hd% dbkz
!              call finish('scatt_obs_err','SCATT: invalid KZ type' )
!            endif
!         end select

          select case (spot% ident)
           case (261,262,441)
             eff = 1.5_wp   ! JASON, SARAL wind error
           case (61,65)
             eff = 2.0_wp   ! Sentinel-3 wind error
           case (3,4,5,421)
             eff = 1.5_wp   ! ASCAT, OSCAT wind error
           case (281,422)
             eff = 2.0_wp   ! QSCAT/ScatSat wind error
           case (503,504,505)
             eff = 2.0_wp   ! HSCAT
           case (801)
             eff = 1.5_wp   ! RSCAT wind error
           case default
             !---------------------------------------------------------
             ! no DBKZ given, probably third party report from fof-file
             ! keep observation error if already set
             !---------------------------------------------------------
             if (obs% o% body(iii)% eo > 0._sp) then
               eff = obs% o% body(iii)% eo
             else
               write(0,*)'SCATT: invalid SatID:',spot% ident
               call finish('scatt_obs_err','SCATT: invalid SatID' )
             endif
          end select

          var = eff ** 2

        case (VN_DD)
          if (ff > 0._wp) then
            eff = r2d * eff / ff
          else
            eff = 360._wp
          endif
          eff = min (eff, 360._wp)
          var = eff ** 2

        case default
          write(0,*)' invalid varno:',it(i)
          call finish('scatt_obs_err','SCATT: invalid varno')
        end select
        obs% R% ja    (k)     = iii
        obs% R% packed(k)     = var
!       obs% o% body(iii)% eo = sqrt (obs% R% packed(k))
        k = k + 1
      end do
    case default
      call finish('scatt_obs_err','invalid observation type ')
    end select
    obs% R% ia (spot% o% i + spot% o% n + 1) = k

  end subroutine scatt_obs_err
!------------------------------------------------------------------------------
  subroutine zero_obs_err (obs, spot)
    !-------------------------------------
    ! Set observation error to null values
    !-------------------------------------
    type(t_obs_block) ,intent(inout) :: obs     ! observation data type
    type(t_spot)      ,intent(in)    :: spot    ! SPOT observations

    integer :: i, k

    k = obs% R% ia (spot% o% i+1)
    do i = spot% o% i + 1, spot% o% i + spot% o% n
       obs% R% ia    (i)     = k
       obs% R% ja    (k)     = i
       obs% R% packed(k)     = 0._sp
       obs% o% body  (i)% eo = 0._sp
       k = k + 1
    end do
    obs% R% ia (spot% o% i + spot% o% n + 1) = k

  end subroutine zero_obs_err
!------------------------------------------------------------------------------
  subroutine temp_obs_err (obs, spot, e, it, z, o, pl, tl, lcov)

  !-------------------------------------------------------------
  ! calculate observation error correlation or covariance matrix
  ! for TEMPs (DWD OI model)
  !-------------------------------------------------------------
  type(t_obs_block) ,intent(inout) :: obs    ! observation data type
  type(t_spot)      ,intent(in)    :: spot   ! SPOT observations
  real(wp)          ,intent(out)   :: e  (:) ! observation error
  integer           ,intent(in)    :: it (:) ! observation data type
  real(wp)          ,intent(in)    :: z  (:) ! vertical coordinate p[hPa]
  real(sp)          ,intent(in)    :: o  (:) ! observed quantities
  real(sp)          ,intent(in)    :: pl (:) ! pressure    on TEMP levels
  real(sp)          ,intent(in)    :: tl (:) ! temperature on TEMP levels
  logical           ,intent(in)    :: lcov   ! .true.: calculate covariances

    !----------------
    ! local variables
    !----------------
    integer  :: n            ! number of observations
    real(wp),allocatable :: Rnew(:) !
    real(wp) :: f            ! height error factor, depending on category
    integer  :: k            ! level index
    integer  :: ierr         ! return code of obs_err
    !----------------------------------------
    ! variables valid for current row, column
    !----------------------------------------
    integer  :: i,j,l,l1     ! indices
    real(wp) :: pi, pj       ! vertical coordinate    (p)
    real(wp) :: zj           ! vertical coordinate (ln(p))
    real(wp) :: xi, xj       ! transformed vertical coordinate
    real(wp) :: ehj          ! height error
    real(wp) :: ff           ! wind speed
    integer  :: iti, itj     ! observation type
    !------------------------------------------------------
    ! arrays to keep quantities for inner loop calculations
    !------------------------------------------------------
    real(wp) :: x (size(z))  ! transformed coordinate
    !--------------------------------------------------
    ! temporaries for final calculation of correlations
    !--------------------------------------------------
    real(wp) :: c            ! height error correlation
    real(wp) :: dx, dx2      ! vertical coordinate difference, squared

    if (not_init) call finish ('temp_obs_err','module not initialized')
    n = size (z)
    allocate (Rnew(n*n))
    Rnew = 0._wp
    !--------------------------------------------
    ! category dependent factor for height errors
    !--------------------------------------------
    f = 1._wp
    do i=1,size(cat)
      if (spot% col% c% dlon >= cat(i)% lon1 .and. &
          spot% col% c% dlon <= cat(i)% lon2 .and. &
          spot% col% c% dlat >= cat(i)% lat1 .and. &
          spot% col% c% dlat <= cat(i)% lat2) then
        f = cat(i)% f
        exit
      endif
    end do
    !-----------------------------
    ! outer loop over observations
    !-----------------------------
    pj = huge (pj)
    do j = 1, n
      l = (j - 1) * n + j
      Rnew(l)  = 1._wp
      !------------
      ! new level ?
      !------------
      if (pj/=z(j)) then
        pj = z(j)
        zj = log(pj)
        ff = 0._wp
        !----------------------------------
        ! polynominal expansion for x and e
        !----------------------------------
        ehj  = 0._wp
        xj   = 0._wp
        do i = rs% nfpars(1), 1, -1
          ehj = ehj * zj + rs% cz (i)
        end do
        do i = rs% nfpars(2), 1, -1
          xj = xj * zj + rs% cx (i)
        end do
        !-----------------------------------------
        ! rescale height error for category 1 or 3
        !-----------------------------------------
        ehj = ehj * f
        !---------------------------------
        ! store temporaries for inner loop
        !---------------------------------
        x (j) = xj
      endif
      !-------------
      ! store errors
      !-------------
      itj = it(j)
!print *,zj,itj,ti,o(j)
      select case (itj)
      case (VN_U, VN_V, VN_FF)
        e(j) = obs_err (spot% hd% obstype,                 &
                        spot% hd% codetype, VN_U, pj, 0._sp)
        if (itj==VN_FF) ff = o(j)
      case (VN_DD)
        if (ff > 0._wp) then
          e(j) = obs_err (spot% hd% obstype,                  &
                          spot% hd% codetype, VN_U, pj, 0._sp)&
               * r2d / ff
        else
          e(j) = 360._wp
        endif
        e(j) = min (e(j), 360._wp)
      case (VN_Z)
        e(j) = ehj
      case (VN_T)                                  ! Temperature
        e(j) = obs_err (OT_TEMP,                          &
                        spot% hd% codetype, VN_T, pj, o(j))
      case (VN_RH)
        e(j) = obs_err (OT_TEMP,                                      &
                        spot% hd% codetype, VN_RH, pj, o(j), ierr=ierr)
        if (ierr /= 0) then
          e(j)                              = obserr_temp_rh
          if (o(j) <= 0.2_wp ) e(j)         = obserr_temp_rh_dry
          do k=1,size(pl)
            if (abs(pl(k)-pj) < 0.1) then
              if (tl(k) <= t0c-40._wp) e(j) = obserr_temp_rh_cold
              exit
            endif
            if (k==size(pl)) call finish('temp_obs_err','no matching p level')
          end do
        end if
      case default
        write(0,*) dace% pe,'temp_obs_err: unknown observation type =',itj
        call finish ('temp_obs_err','unknown observation type.')
      end select
      !-----------------------------
      ! inner loop over observations
      !-----------------------------
      pi = huge (pi)
      do i=1,j
        l=(j - 1) * n + i
        !------------
        ! new level ?
        !------------
        if (pi /= z(i)) then
          pi  = z  (i)
          !---------------------------------
          ! load temporaries from outer loop
          !---------------------------------
          xi  = x  (i)
          !-----------------------------------------------------
          ! height correlation function, 1st, and 2nd derivative
          !-----------------------------------------------------
          dx    = xi-xj
          dx2   = dx**2
          c     = a * exp ( -b * dx2 )
        endif
        !-----------------------
        ! calculate correlations
        !-----------------------
        iti = it(i)
        !-------------
        ! geopotential
        !-------------
        if     (iti==VN_Z  .and. itj==VN_Z ) then
          Rnew(l) = c
        endif
        !-----------------------------------
        ! account for nugget effect if a < 1
        !-----------------------------------
        if (i==j) then
          Rnew(l) = max (Rnew(l), 1._wp)
        endif
        !-----------------------------------
        ! return variances, not correlations
        !-----------------------------------
        if (lcov) Rnew(l) = Rnew(l) * e(i) * e(j)
        !--------------------------------
        ! mirror lower triangle of matrix
        !--------------------------------
        l1 = (i - 1) * n + j
        Rnew(l1) = Rnew(l)
      end do
    end do
    deallocate (Rnew)
  end subroutine temp_obs_err
!------------------------------------------------------------------------------
  function obs_err (obstype, codetype, qty, p, o, adl, ierr) result (err)
  !----------------------------------------------------------------------------
  ! Find the appropriate observation error table
  ! based on obstype/codetype.
  ! If a  matching codetype is specified in the table the entry is preferred.
  ! If no matching codetype is found the table with unspecified codetype (-1)
  ! is taken.
  ! Observation error is interpolated linearily in height (pressure).
  ! If 'ierr' is present and no entry is found the function returns with ierr=1
  ! If 'ierr' is not present and no entry is found the function aborts
  !----------------------------------------------------------------------------
  real(wp)                             :: err      ! observation error
  integer        ,intent(in)           :: obstype  ! CMA obstype
  integer        ,intent(in)           :: codetype ! CMA codetype
  integer        ,intent(in)           :: qty      ! observed quantity (varno)
  real(wp)       ,intent(in)           :: p        ! height (Pa)
  real(sp)       ,intent(in)           :: o        ! observed value
  logical        ,intent(in) ,optional :: adl      ! use table adaptloc
  integer        ,intent(out),optional :: ierr     ! error return variable

    integer  :: i, k, n
    integer  :: io     ! table index for matching obstype / any codetype
    integer  :: ioc    ! table index for matching obstype / codetype
    real(wp) :: w1, w2 ! interpolation weights
    type(t_obserr) ,pointer :: t (:)


    if (not_init) call finish ('obs_err','module not initialized')
    if (present (ierr)) ierr = 0
    t => use_obserr
    n =  n_use
    if (present(adl)) then
      if (adl) then
        t => adaptloc
        n =  n_adp
      endif
    endif

    !--------------------------
    ! find matching table entry
    !--------------------------
    io  = -1
    ioc = -1
    if (codetype >= 0) then
      !-----------------------------------
      ! codetype specified for observation
      !-----------------------------------
      do i = 1, n
        if   (    t(i)% obstype     == obstype &
            .and. t(i)% varno       == qty      ) then
          if (    t(i)% codetype(1) == -1       ) io  = i ! c.type in table
          if (any(t(i)% codetype    == codetype)) ioc = i !    not in table
        endif
      end do
    else
      !---------------------------------------
      ! codetype not specified for observation
      !---------------------------------------
      do i = 1, n
        if (t(i)% obstype     == obstype .and. &
            t(i)% varno       == qty     .and. &
            t(i)% codetype(1) == -1            ) io  = i
      end do
    endif
    i = io; if (ioc > 0) i = ioc
    !------------------------
    ! no matching table found
    !------------------------
    if (i < 1) then
      if (present (ierr)) then
        ierr =  1
        err  = -1._wp
        return
      else
        write (0,*) 'obs_err: not in table: obstype, codetype, qty =',&
                                            obstype, codetype, qty
        write (0,*) '                      ',rept_char(obstype)% mnem, qty
        call finish ('obs_err','not in table: '//rept_char(obstype)% mnem)
      endif
    endif
    !-------------------------------------------
    ! table found, interpolate in between levels
    !-------------------------------------------
    if (p>=pt(1)) then
      err = t(i)% err ( 1)
    elseif (p<=pt(nlev)) then
      err = t(i)% err (nlev)
    else
      do k = 2, nlev
        if (pt(k)<p) then
          w1 = (p-pt(k))/(pt(k-1)-pt(k))
          w2 = 1._wp-w1
          err = max (w2 * t(i)%   err(k)          + w1 * t(i)%   err(k-1),        &
                     w2 * t(i)% r_err(k) * abs(o) + w1 * t(i)% r_err(k-1) * abs(o))
          exit
        endif
      end do
    endif

  end function obs_err
!==============================================================================
  subroutine read_stn_obserr ()
    !-----------------------------------------
    ! Read station specific observation errors
    !-----------------------------------------
    character(len=8)  :: statid
    character(len=32) :: wsi
    character(len=10) :: quantity        ! quantity         ('uv','tv', 'rh')
    integer           :: codetype(4)     ! CMA codetype
    real(wp)          :: err             ! observation error

    namelist /OBSERR_SYNOP/ statid, wsi, codetype, quantity, err

    logical           :: first   ! flag for first namelist group
    integer           :: ierr    ! namelist position status variable
    integer           :: qnt     ! quantity (varno)
    integer           :: i       ! loop index
    type(t_wsi)       :: tmpwsi  ! Decoded WIGOS station ID
    integer(i8)       :: wsihash ! WIGOS station ID (internal hash)
    character(16)     :: c16     ! Temporary

    integer, parameter :: mentry = 10000 ! estimated max no. individual entries

    type(t_stn_obserr), pointer :: list(:) => NULL()

    nentry = 0

    bigend = .not. little()

    if (dace% lpio) then
       allocate (list(mentry))
       first = .true.
       do  ! read namelist group
          call position_nml ('OBSERR_SYNOP', lrewind=first, status=ierr)
          first = .false.

          statid   = ""
          wsi      = ""
          codetype = -1
          quantity = ""
          err      = 0._wp
          if (ierr == POSITIONED) then
             !------------------------------------------------------------
             ! namelist group positioned: read path/filename from namelist
             !------------------------------------------------------------
#if defined(__ibm__)
             read (nnml ,nml=OBSERR_SYNOP, iostat=ierr)
             if (ierr/=0) call finish ('read_stn_obserr',                &
                                       'ERROR in namelist /OBSERR_SYNOP/')
#else
             read (nnml ,nml=OBSERR_SYNOP)
#endif
             if (statid == "" .and. wsi == "" .and. all (codetype <= 0)) &
                  call finish ("read_stn_obserr","missing station ID pattern")
             !--------------------------------------
             ! convert quantity to index key (varno)
             !--------------------------------------
             qnt = 0
             if (quantity /= "") qnt = value_name (varno, toupper (quantity))
             if (qnt == 0) call finish ("read_stn_obserr",        &
                                        "invalid varno "//quantity)
             if (err <= 0) call finish ("read_stn_obserr", "err <= 0")
             nentry = nentry + 1
             if (nentry > mentry) call finish ("read_stn_obserr","n > mentry")

             !-------------
             ! handle WSI's
             !-------------
             select case (wsi(1:2))
             case ("0-")
                if (wsi_mode == 0) call finish ("read_stn_obserr","wsi_mode=0")
                call wsi_decode (wsi, tmpwsi, ierr)
                if (ierr /= 0) &
                     call finish ("read_stn_obserr","bad WSI: " // trim (wsi))
                list(nentry)% wsi   = tmpwsi

                if (statid == "") then
                   ! Derive statid from WSI
                   call wsi_encode (tmpwsi, statid, wsihash)
                else if (statid(1:1) /= wsi_prefix) then
                   call finish ("read_stn_obserr","bad prefix: "//trim (statid))
                else if (statid(2:8) /= "......." .and. &
                         statid(5:8) /= "...."          ) then
                   call finish ("read_stn_obserr","unsupported: "//trim (statid))
                   ! Assume statid has a valid pattern
                end if
             case ("  ")
                ! No WSI
             case default
                call finish ("read_stn_obserr","bad WSI: " // trim (wsi))
             end select

             list(nentry)% statid   = statid
             call set_mask (list(nentry))
             list(nentry)% varno    = qnt
             list(nentry)% code     = codetype
             list(nentry)% err      = err
          else
             exit
          endif
       end do

       if (nentry == 0) then
          write(6,*) 'Namelist /OBSERR_SYNOP/ not present'
       else
          allocate (stn_obserr(nentry))
          stn_obserr(:) = list(1:nentry)
          write(6,*) 'Namelist /OBSERR_SYNOP/:'
          write(6,'()')
          write(6,'(a)') "  statid      wsi                             &
               &        codetype  quantity       err"
          do i = 1, nentry
             statid   = list(i)% statid
             if (list(i)% wsi% valid) then
                call wsi_to_text (list(i)% wsi, wsi)
             else
                wsi   = ""
             end if
             quantity = name_value (varno, list(i)% varno)
             write(c16,'(4i4)') pack (list(i)% code, list(i)% code >=0)
             write(6,'(2x,a,4x,a,a16,2x,a,f8.3)')               &
                  statid, wsi, trim (c16), quantity, list(i)% err
          end do
       end if
       write(6,'()')
       deallocate (list)
    end if
    call bcast_stn_obserr (stn_obserr, dace% pio)
  end subroutine read_stn_obserr
!------------------------------------------------------------------------------
  subroutine bcast_stn_obserr (list, source)
    type(t_stn_obserr), pointer    :: list(:)
    integer           , intent(in) :: source

#if (defined (__GFORTRAN__) && (__GNUC__ >= 10))
    !----------------------------------------------
    ! include interfaces for external bcast routine
    !----------------------------------------------
    interface p_bcast_derivedtype
       subroutine p_bcast_derivedtype (buffer, count, source, comm)
         type(*) ,intent(inout)     :: buffer(*)    ! variable to bcast
         integer ,intent(in)        :: count        ! len(byte) of variable
         integer ,intent(in)        :: source       ! source processor index
         integer ,intent(in)        :: comm         ! communicator
       end subroutine p_bcast_derivedtype
    end interface p_bcast_derivedtype
#endif

    type(t_stn_obserr) :: d_err ! Dummy data type
    integer            :: s_err ! Size of data type in bytes

    call p_bcast (nentry, dace% pio)
    if (.not. dace% lpio) then
       allocate (list(nentry))
    end if
    if (nentry > 0) then
       s_err = size (transfer (d_err, ["*"]))
       call p_bcast_derivedtype (list, nentry*s_err, source, dace% comm)
    end if
  end subroutine bcast_stn_obserr
!------------------------------------------------------------------------------
  subroutine set_mask (e)
    type (t_stn_obserr), intent (inout) :: e
    !------------------------------------------
    ! derive integer representation of wildcard
    !------------------------------------------
    integer(i8) :: m     ! mask variable
    integer     :: i     ! index variable
    e% len    = len_trim  (e% statid)
    e% id     = transfer  (e% statid, e% id)
    if (bigend) call flip (e% id)
    e% mask   = not (  0_i8)
    m         = not (255_i8)
    do i=1,8
      if(e% statid(i:i)=='.') then
        e% id   = iand (e% id,   m)
        e% mask = iand (e% mask, m)
      endif
      m = ior (ishft(m, 8), 255_i8)
    end do
  end subroutine set_mask
!------------------------------------------------------------------------------
  function obserr_entry (statid, wsi, code, entry) result (found)
    character(len=8)   ,intent(in)  :: statid    ! station identifier
    type(t_wsi)        ,intent(in)  :: wsi       ! WIGOS station identifier
    integer            ,intent(in)  :: code      ! observation code type
    type(t_stn_obserr) ,intent(out) :: entry(:)  ! station obserr entry
    integer                         :: found     ! number of matching entries
    !----------------
    ! local variables
    !----------------
    integer                     :: j, len, me
    integer(i8)                 :: id
    type(t_stn_obserr), pointer :: e
    !----------------------
    ! executable statements
    !----------------------
    found = 0
    if (nentry <= 0) return
    me  = size (entry)
    len = len_trim (statid)
    id  = transfer (statid, id)
    if (bigend) call flip (id)
    do j = 1, nentry
       e => stn_obserr(j)
       if (e% code(1) >= 0 .and. code >= 0) then
          if (all (e% code /= code)) cycle
       end if
       if ((e% len  == 0)                 .or. &
           (e% len  == len     .and.           &
            e% id   == iand (e% mask, id))     ) then
          if (e% wsi% valid) then
             if (.not. (e% wsi == wsi)) cycle
          end if
          found = found + 1
          if (found <= me) entry(found) = e
       endif
    end do
    !----------------------------------------------
    ! check for sufficient size of output parameter
    !----------------------------------------------
    if (found > me) then
       call message ('obserr_entry',"increase max.number of entries 'me'")
       found = - found
    endif
  end function obserr_entry
!==============================================================================
end module mo_obs_err
