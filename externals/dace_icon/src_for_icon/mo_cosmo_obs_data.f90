!
!+ Variables used by mo_cosmo_obs for computing local information
!
MODULE mo_cosmo_obs_data
!
! Description:
!  This module contains variables (scalars and arrays) which need to be
!  available in the observation reading part or for quality control
!   These are
!    - NAMELIST variables controlling the reading of observation
!    - some parameters (constants) and associated variables
!    - some I/O device numbers and file names for nudging
!    - surface analysis limits and allocatable arrays
!    - some other variables, e.g. the switch for the observation processing
!
!  variables are taken from data_nudge_all.f90 by C. Schraff
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_11        2010/06/16 Tanja Weusthoff
!  Variables used by mo_cosmo_obs for computing local information
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_22        2013-02-13 Christoph Schraff
!  adapt to COSMO multilevel observation operator
! V1_43        2015-08-19 Andreas Rhodin
!  COSMO MEC: increase defaults for 'maxmlo', 'maxsgo'.
!   add 'verification_end' to list of namelist variables.
! V1_44        2015-09-30 Andreas Rhodin
!  update shared modules to COSMO 5.03-beta
! V1_47        2016-06-06 Andreas Rhodin
!  add namelist variables; increase defaults maxmlo:->1000, maxsgo:->10000
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
!-------------------------------------------------------------------------------
!
! Modules used:
!
!-------------------------------------------------------------------------------

use mo_kind, only : wp       ! working precision kind parameter

use mo_run_params, only :   &
    obsinput,               &! observation input directory
    path_file                ! concatenate: path / file . suffix

use mo_namelist,  only:     &
    position_nml,           &! routine to position nml group
    nnml,                   &! namelist fortran unit number
    POSITIONED               ! position_nml: OK return flag

use mo_mpi_dace,  only:     &
    dace,                   &! MPI group info
    p_bcast                  ! generic broadcast routine

!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: read_nml_cosmo_obs, lverpas, lcloud_ice, icdfdirlen, mxav, mxgpc,   &
          qcc, qccsu, qcvf, qcciq, qcsiq,                                     &
          obnlat, obslat, obwlon, obelon, exnlat, exslat, exwlon, exelon,     &
          doromx, altopsu, zlimv10, dhosag, av_levs, av_incr, av_reso,        &
          thairh, rhtsat, itim_wp, icdt_tws, icdt_rss, ilocv_sfc, fplev_ps,   &
          rtmlrs, rtmlsy, rtmlair, rtmltow, rtmlrsy, rtmlim,                  &
          mxtwex, ntwex,  htwex,  ivtwex, ytwex, ilstidtw,                    &
          mxbcrr, nbcrr,  ybcrr,  bcrrt,  bcrrhl, bcrrhu, lredn_repro,        &
          lsytac, maxmlo, maxsgo, maxgpo, maxmlv, nolbc, mqcorr92,            &
          lsynop, laircf, lsatob, ldribu, ltemp, lpilot, lsatem, lgps, lscatt,&
          lcd011, lcd014, lcd021, lcd024, lcd140, lcd811, lcd835, lcd839,     &
          lcd041, lcd141, lcd144, lcd146, lcd244,                             &
          lcd088, lcd090,         lcd064, lcd165,         lcd039, lcd040,     &
          lcd035, lcd036, lcd037, lcd135, lcd109, lcd111, lcd230, lcd231,     &
          lcd032, lcd033, lcd038, lcd132, lcd133, lcd136, lcd137, lcd139,     &
          lcd159, lcd187,                                                     &
          lcd086, lcd186,         lcd122, lcd123,         lcd096, igpscen,    &
          ycdfdir, verification_start, verification_end, ionl, jonl

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 1:  General parameters and related variables
!-------------------------------------------------------------------------------

  INTEGER        , PARAMETER  :: &
    icdfdirlen = 250 ,& ! max. length of name of directory where
                        !   NetCDF observation input files reside
    mxav       =  15 ,& ! max. length of level definition list for
                        !   superobbing of high-resolution radisondes
    mxgpc      =  20 ,& ! max. number of GPS processing centres used
    mxtwex     =  20 ,& ! max. number of exceptions for tower processing
    mxbcrr     =  20 ,& ! max. number of RASS bias correction rules
    ilstidtw   =  10    ! max. length of tower station ID's

  INTEGER        ::       &
    ntwex            ,& ! number of exceptions for tower processing
    nbcrr               ! number of RASS bias correction rules


!-------------------------------------------------------------------------------
! Section 2:  Namelist variables controlling the data assimilation
!-------------------------------------------------------------------------------

!      0.    General steering switches
!      -------------------------------

  LOGICAL        ::       &
    lverpas      ,& ! .t.: on - off switch for verif.also of passive reports
    lcloud_ice      ! .t.: on - off switch for cloud_ice in grid-scale
                    !      precipitaiton physics (.t. if itype_gscp > 2)

!      1.    Time window for reading of observations
!      ---------------------------------------------

  INTEGER        ::       &
    verification_start ,& ! - 29 : start of time window (minutes before time_verif)
    verification_end      !    0 :  end  of time window (minutes before time_verif)

!      7.    Threshold quality control
!      -------------------------------

  REAL (wp)      ::       &
!      for upper-air data
    qcc      (4) ,& !  0.,500: constant part of the 'quality control thresholds'
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf     (4) ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (for height/ thickness instead of
                    !          pressure ps, and not available for humidity RH),
                    !          given on following pressure levels:
                    !          level          ,hPa:1000, 850, 700, 500, 400, 300
                    !               ,250, 200, 150, 100,  70,  50,  30,  20,  10
                    !          TEMP  wind     ,m/s: 2.3, 2.3, 2.5, 3.0, 3.5, 3.7
                    !               ,3.5, 3.5, 3.4, 3.3, 3.2, 3.2, 3.3, 3.6, 4.5
                    !          AIREP wind     ,m/s: 2.5, 2.5, 3.0, 3.5, 4.0, 4.0
                    !               ,4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0
                    !          TEMP  temperature,K: 1.2, 1.0, 0,7, 0.4, 0.4, 0.5
                    !               ,0.5, 0.6, 0,7, 0.8, 0.8, 0,9, 0.9, 1.0, 1.2
                    !          AIREP temperature,K: 1.2 ,1.0, 0.7, 0.5, 0.5 ,0.6
                    !               ,0.6 ,0.7, 0.8, 0.9, 1.0, 1.1, 1.1, 1.2, 1,4
                    !          (cf. 'data_nudge_local' for tables, thickness QC,
                    !           and other factors used for thresholds, e.g. an
                    !           additional factor of 0.8 used for AIREP wind)
!      for surface-level data
    qccsu    (4) ,& ! 12.,500: constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
!      for integrated water vapour (IWV) derived from GPS or radiosonde data
    qcciq        ,& !  1.   constant part of QC threshold for IWV
    qcsiq           !  .15  IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile

!!      8.    Observation processing
!-------------------------------------------------------------------------------

!      8.0   Reading of observation reports
!      ------------------------------------

  CHARACTER (LEN=icdfdirlen) ::       &
    ycdfdir         ! './'   : directory where NetCDF observation input files
                    !            and the blacklist file reside


!      8.1   Use of stations / reports
!      -------------------------------

  REAL (wp)      :: &
    obnlat       ,& !   90. : northern boundary of observation area
    obslat       ,& !  -90. : southern boundary of observation area
    obwlon       ,& ! -180. : western boundary of observation area
    obelon       ,& !  180. : eastern boundary of observation area
    exnlat       ,& !   90. : northern boundary for exclusion area
    exslat       ,& !  -90. : southern boundary for exclusion area
    exwlon       ,& ! -180. : western boundary for exclusion area
    exelon       ,& !  180. : eastern boundary for exclusion area
    doromx   (4) ,& !  100.,: vertical extrapolation cut-off and gaussian
                    !  150.,  radius of height differences between model
                    !  150.,  orography and surface station height for a factor
                    !  150.,  contributing to the quality weight factor as part
                    !         of the nudging weights
                    !         (height diff of SYNOP/ GPS with (z-obs > z-model,
                    !          --> interpolation instead of extrapolation) are
                    !          divided by 4 (fdoro) for surf pressure/ IWV obs)
    altopsu (4)  ,& !  100.,: SYNOP obs. above height 'altopsu' are not assimi-
                    ! 3*5000. lated. If (altopsu == 0.) then SYNOP / surf. TEMP
                    !         assigned to land grid pts. are not assimilated
    zlimv10 (3)  ,& !  400.,: additional limits for the use of 10-m wind obs:
                    !  800.,  positive resp. negative scaled Laplacian of
                    !    5.1  orography [m], surface roughness length [m]
                    !         reasonable values would be: 75., 200., 0.7
    dhosag       ,& !    5.: for active use of tower obs profiles:
                    !        limit in terms of height of sensor above ground [m]
                    !          below which model equivalents are computed at the
                    !                observed height of sensor above ground, and
                    !          above which model equivalents are computed at the
                    !                observed height (altitude)
    thairh       ,& !    0. : maximum horizontal distance [km] between the
                    !         lowest report and any single level report that
                    !         is added to a multi-level AIRCRAFT report
    rhtsat  (5)  ,& ! relative humidity threshold above which obs is set =100%
                    !   (1: all, 2: surface obs, 3: TEMP, 4: PILOT, 5: AIREP;
                    !    non-default values 2 - 5 overrule value at index 1)
    rtmlrs       ,& ! redundancy time limit for radiosondes              [hrs]
    rtmlsy       ,& ! redundancy time limit for SYNOP                    [hrs]
    rtmlair      ,& ! redundancy time limit for AIREP                    [hrs]
    rtmltow      ,& ! redundancy time limit for tower                    [hrs]
    rtmlrsy      ,& ! redundancy time limit for raso sfc. level vs Synop [hrs]
    rtmlim       ,& ! redundancy time limit for other obs                [hrs]
    fplev_ps        !    1.  : factor to 'plevel' (for localisation in LETKF)
                    !            for surface pressure obs from surface stations

  REAL (wp)      , TARGET :: &
    av_levs(mxav),& !        : level definition list       \  for superobbing
                    ! (1075., 755., 710., 90., 75., 5.)     \ layers of
    av_incr(mxav),& !        : level increment  list        / high-resolution
                    ! (  10.,  15.,  20., 15., 10., 0.)    /  radiosonde reports
    av_reso         !    3.  : apply superobbing if the averaged resolution of
                    !          the observed profile exceeds 'av_reso' times the
                    !          model resolution

  LOGICAL        ::       &
    lredn_repro  ,& ! .f. ensure reproducibility of redundancy check irrespect.
                    !       of domain decomposition by allowing for redundancy
                    !       only between reports assigned to the same grid point
    lsytac          ! .t.   : if .t.: TAC preferred over BUFR Synop reports
                    !         if .f.: vice versa; this uses DWD DB KZ

  INTEGER        :: &
    maxmlo      ,& !  600   : max. number of multi-level reports
    maxsgo      ,& ! 4000   : max. number of (surface-level and upper-air)
                   !                         single-level reports
    maxuso      ,& !  900   : max. number of upper-air single-level reports
    maxgpo      ,& !  200   : max. number of GPS reports within total domain
    maxmlv      ,& !  100   : max. number of observation levels in multi-level
                   !                         reports
!   maxtvo      ,& !  200   : max. number of sat retrievals within totaldomain
    nolbc       ,& !    5   : number of grid rows at lateral boundaries
                   !          where obs are neglected
    itim_wp     ,& !    0   : mode of correction of obs time for wind profiler
                   !          = 0 : no correction
                   !          = 1 : correct by half of the obs averaging period
                   !          = 2 : correct to the end of the obs aver. period
                   !          = 3 : correct to the beginning of the obs period
    icdt_tws    ,& !    0   : mode for obs/code type of tower surface-level rep.
                   !          = 0 : no correction (remains PILOT/ ICOS tower)
                   !          = +/- 1 : automatic Synop
                   !          = +/- 2 : Synop code type 839
                   !          > 0: obstype status active if Synop to be active
                   !          < 0: active only if both Synop + tower active
    icdt_rss    ,& !    0   : mode for radiosonde-derived surface-level report
                   !          =-1 : no single-level surface report (SR) created
                   !          = 0 : SR: type + status like parent TEMP
                   !          = 1 : SR: obs type SYNOP, code 835, active only if
                   !                parent TEMP code type active and lcd835 true
                   !          = 2 : SR: obs type SYNOP, code 835, status lcd835
    ilocv_sfc   ,& !    0   : mode for setting 'plevel' in surface reports
                   !          = 0 : set to observed station pressure only if
                   !                present (else set missing value
                   !                 -> US Std Atmos used in LETKF for v-loc.)
                   !          = 1 : set to model surface pressure instead
                   !                if station pressure not observed
                   !          = 2 : always set to model surface pressure
    mqcorr92       !    0   : switch for bias correction for Vaisala RS92
                   !            radiosonde humidity
                   !          = 0 : no correction for humidity
                   !          = 1 : correct only solar radiation bias
                   !          = 2 : correct total bias (incl. nighttime bias)

  INTEGER        , TARGET   :: &
    igpscen (mxgpc) ! X* -1  : array of processing centres of GPS reports used
                    !          actively (order of centres determines preference
                    !          in redundancy check; '-1' means no active centre)

!      8.1.1  List of 'exception levels' of single tower stations for which
!             the active use of data shall deviate from that defined by 'dhosag'
!      -------------------------------------------------------------------------
  REAL (wp)      , TARGET   :: &
    htwex  (mxtwex) ! X* -1. : heights of sensor above ground [m] of levels;
                    !            if = 0 then apply to all levels

  INTEGER        , TARGET   :: &
    ivtwex (mxtwex) ! X* 0   : variable indicators related to exception levels:
                    !          O: none; 1: wind; 2:temperature + humidity; 3: all

  CHARACTER (LEN=ilstidtw) , TARGET :: &
    ytwex  (mxtwex) ! X*' '  : station ID's related to exception levels

!      8.1.2  List of 'bias correction rules' for RASS reports
!      -------------------------------------------------------
!             (for defaults, see correction values set in Section 1 for
!              stations 10266, 10394, 10678; in order to switch off the
!              bias correction, set bcrrt(1)=0., or ybcrr(1)=' ')
  REAL (wp)      , TARGET   :: &
    bcrrt  (mxbcrr) ,& !     : bias correction of RASS (virtual) temperature [K]
    bcrrhl (mxbcrr) ,& !     : lower \ height limit [m] in RASS profile to apply
    bcrrhu (mxbcrr)    !     : upper / bias correct. 'bcrrt' for station 'ybcrr'

  CHARACTER (LEN=ilstidtw) , TARGET :: &
    ybcrr  (mxbcrr)    !     : station ID's for RASS bias correction rules

!      8.2   Use of observation and code types
!      ---------------------------------------

  LOGICAL        :: &
    lsynop       ,& ! .t.    : .t. if SYNOP data is used
    laircf       ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob       ,& ! .false.: .t. if SATOB data is used
    ldribu       ,& ! .t.    : .t. if BUOY  data is used (drifting buoy)
    ltemp        ,& ! .t.    : .t. if TEMP  data is used
    lpilot       ,& ! .t.    : .t. if PILOT data is used
    lsatem       ,& ! .false.: .t. if SATEM data is used
    lgps         ,& ! .false.: .t. if GPS   data is used
!   lgnssstd     ,& ! .false.: .t. if GNSS-STD data is used
    lscatt          ! .t.    : .t. if SCATT data is used (scatterometer)

  LOGICAL        :: &
    lcd011       ,& ! .t.    : synop code  11 data is used (land synop)
    lcd014       ,& ! .t.    : synop code  14 data is used (automatic)
    lcd021       ,& ! .t.    : synop code  21 data is used (ship)
!   lcd022       ,& ! .t.    : synop code  22 data is used (ship abbrev.)
!   lcd023       ,& ! .t.    : synop code  23 data is used (shred)
    lcd024       ,& ! .t.    : synop code  24 data is used (autom. ship)
    lcd140       ,& ! .t.    : synop code 140 data is used (metar)
    lcd811       ,& ! .t.    : synop code 811 data is used (synop test)
    lcd835       ,& ! .t.    : synop code 835 data is used (surface TEMP)
    lcd839       ,& ! .t.    : synop code 839 data is used (surface tower)
    lcd041       ,& ! .t.    : airep code  41 data is used (codar)
    lcd141       ,& ! .t.    : airep code 141 data is used (airep)
!   lcd241       ,& ! .t.    : airep code 241 data is used (colba)
    lcd144       ,& ! .t.    : airep code 144 data is used (amdar)
    lcd146       ,& ! .t.    : airep code 146 data is used (mode-s)
    lcd244       ,& ! .t.    : airep code 244 data is used (acars)
    lcd088       ,& ! .t.    : satob code  88 data is used (satob)
    lcd090       ,& ! .t.    : satob code  90 data is used (amv)
!   lcd188       ,& ! .f.    : satob code 188 data is used (sst)
!   lcd063       ,& ! .t.    : dribu code  63 data is used (bathy)
    lcd064       ,& ! .t.    : dribu code  64 data is used (tesac)
    lcd165       ,& ! .t.    : dribu code 165 data is used (drift. buoy)
    lcd035       ,& ! .t.    : temp  code  35 data is used (land temp)
    lcd036       ,& ! .t.    : temp  code  36 data is used (temp ship)
    lcd037       ,& ! .t.    : temp  code  37 data is used (mobile)
    lcd135       ,& ! .t.    : temp  code 135 data is used (dropsonde)
    lcd109       ,& ! .t.    : temp  code 109 data is used (land temp hi-res)
    lcd111       ,& ! .t.    : temp  code 111 data is used (ship temp hi-res)
    lcd230       ,& ! .t.    : temp  code 230 data is used (drop temp hi-res)
    lcd231       ,& ! .t.    : temp  code 231 data is used (desc temp hi-res)
    lcd039       ,& ! .t.    : temp  code  39 data is used (rocob)
    lcd040       ,& ! .t.    : temp  code  40 data is used (rocob ship)
    lcd032       ,& ! .t.    : pilot code  32 data is used (land pilot)
    lcd033       ,& ! .t.    : pilot code  33 data is used (pilot ship)
    lcd038       ,& ! .t.    : pilot code  38 data is used (mobile)
    lcd132       ,& ! .t.    : pilot code 132 data is used (win-prof eu)
    lcd133       ,& ! .t.    : pilot code 133 data is used (sod/rass eu)
    lcd136       ,& ! .t.    : pilot code 136 data is used (pro/rass us)
    lcd137       ,& ! .t.    : pilot code 137 data is used (Radar VAD)
    lcd139       ,& ! .t.    : pilot code 139 data is used (tower)
    lcd159       ,& ! .f.    : pilot code 159 data is used (icos tower)
    lcd187       ,& ! .f.    : pilot code 187 data is used (wind lidar)
    lcd086       ,& ! .t.    : satem code  86 data is used (satem)
    lcd186       ,& ! .t.    : atovs code 186 data is used (hi-res ATOVS)
    lcd122       ,& ! .t.    : scatt code 122 data is used (QuickScat)
    lcd123       ,& ! .t.    : scatt code 123 data is used (ASCAT)
    lcd096          ! .t.    : gps data from COST ASCII file is used

! INTEGER        :: &
!   mcdmsg1      ,& !  processing / use of MSG1   code  71 data
!   mcdmsg2      ,& !  processing / use of MSG2   code  72 data
!   mcdno15      ,& !  processing / use of NOAA15 code 206 data
!   mcdno16      ,& !  processing / use of NOAA16 code 207 data
!   mcdno17      ,& !  processing / use of NOAA17 code 208 data
!   mcdno18         !  processing / use of NOAA18 code 209 data

!      10.   Diagnostic output
!      -----------------------
  INTEGER        :: &
    ionl         ,& ! 167    : / grid point coordinates
    jonl            ! 103    : \ for standard output on nudging


!-------------------------------------------------------------------------------


CONTAINS

!===============================================================================
!+ Internal procedure in "data_nl_cosmo_obs" for NAMELIST input
!-------------------------------------------------------------------------------

SUBROUTINE read_nml_cosmo_obs

! Local variables:
! ---------------

!      0.    General steering switches
!      -------------------------------

  LOGICAL        :: &
    lverpas_d     ,&  ! .t. : on - off switch for verif. also of passive report
    lcloud_ice_d      ! .t. : on - off switch for cloud_ice in grid-scale
                      !       precipitaiton physics (.t. if itype_gscp > 2)

!      1.    Time window for reading of observations
!      ---------------------------------------------

  INTEGER        :: &
    verification_end_d ,& !  0  : end of time window (minutes before time_verif)
    verification_start_d  ! -29 : start of time window (min. before time_verif)

!      7.1   (Threshold) Quality control
!      ---------------------------------

  REAL (wp)      :: &
!      for upper-air data
    qcc_d    (4) ,& !  0.,500: constant parts of the quality control thresholds
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf_d   (4) ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (as def. in 'data_nudge_local')
!      for surface-level data
    qccsu_d  (4) ,& ! 12.,500: constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
!      for integrated water vapour (IWV) derived from GPS or radiosonde data
    qcciq_d      ,& !  1.   constant part of QC threshold for IWV
    qcsiq_d         !  .15  IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile

  CHARACTER (LEN=128)       ::       &
    ycdfdir_d         !     : directory where NetCDF observation input files
                      !       and the blacklist file reside

!      8.1   Use of stations / reports
!      -------------------------------
  REAL (wp)      :: &
    obnlat_d     ,& !   90.  : northern boundary of observation area
    obslat_d     ,& !  -90.  : southern boundary of observation area
    obwlon_d     ,& ! -180.  : western boundary of observation area
    obelon_d     ,& !  180.  : eastern boundary of observation area
    exnlat_d     ,& !   90.  : northern boundary for exclusion area
    exslat_d     ,& !  -90.  : southern boundary for exclusion area
    exwlon_d     ,& ! -180.  : western boundary for exclusion area
    exelon_d     ,& !  180.  : eastern boundary for exclusion area
    doromx_d (4) ,& !  100., : cut-off and gaussian radius of height differences
                    !  400.,   between model orography and station height for a
                    !  160.,   factor contributing to the quality weight factor
                    !  160.,   as part of the nudging weights
                    !          (height diff. of SYNOPs with (z-obs < z-model)
                    !           are multiplied by 4 for surface pressure obs.)
    altopsu_d(4) ,& !  100., : SYNOP obs. above height 'altopsu' are not assimi-
                    ! 3*5000.  lated. If (altopsu == 0.) then SYNOP / surf. TEMP
    zlimv10_d(3) ,& !  400., : additional limits for the use of 10-m wind obs:
                    !  800.,   positive resp. negative scaled Laplacian of
                    !    5.1   orography [m], surface roughness length [m]
                    !          reasonable values would be: 75., 200., 0.7
    dhosag_d     ,& !    5.: for active use of tower obs profiles:
                    !        limit in terms of height of sensor above ground [m]
                    !          below which model equivalents are computed at the
                    !                observed height of sensor above ground, and
                    !          above which model equivalents are computed at the
                    !                observed height (altitude)
    thairh_d     ,& !    0.  : maximum horizontal distance [km] between the
                    !          lowest report and any single level report that
                    !          is ad
    rhtsat_d (5) ,& !   0.96 : relative humidity threshold above which a
                    !            (relative) humidity obs is set to saturation
                    !            (1: all, 2: sfc, 3: TEMP, 4: PILOT, 5: AIREP;
                    !            non-default values at 2-5 overrule value at 1)
    rtmlrs_d     ,& !   0.751: redundancy time limit for radiosondes       [hrs]
    rtmlsy_d     ,& !   0.35 : redundancy time limit for SYNOP             [hrs]
    rtmlair_d    ,& !   0.25 : redundancy time limit for AIREP             [hrs]
    rtmltow_d    ,& !   0.51 : redundancy time limit for tower             [hrs]
    rtmlrsy_d    ,& !   0.99 : red time limit for raso sfc. level vs Synop [hrs]
    rtmlim_d     ,& !   0.15 : redundancy time limit for other obs         [hrs]
    fplev_ps_d   ,& !    1.  : factor to 'plevel' (for localisation in LETKF)
                    !            for surface pressure obs from surface stations
    av_levs_d(mxav),& !        : level definition list     \  for superobbing
                      ! (1075., 755., 710., 90., 75., 5.)   \ layers of
    av_incr_d(mxav),& !        : level increment  list      / high-resolution
                      ! (  10.,  15.,  20., 15., 10., 0.)  /  radiosonde reports
    av_reso_d      ,& !    3.  : apply superobbing if the averaged resolution of
                      !          the observed profile exceeds 'av_reso' times
                      !          the model resolution
    htwex_d (mxtwex),&! X* -1. : sensor heights above ground [m] of exception
                      !          levels; if = 0 then apply to all levels
    bcrrt_d (mxbcrr),&!      : bias correction of RASS (virtual) temperature [K]
    bcrrhl_d(mxbcrr),&!      : lower \ height limit [m] in RASS profile to apply
    bcrrhu_d(mxbcrr)  !      : upper / bias correct. 'bcrrt' for station 'ybcrr'

  LOGICAL        ::       &
    lredn_repro_d,& ! .f. : if true then reproducibility of redundancy check
                    !       ensured irrespective of domain decomposition by
                    !       allowing for redundancy only between reports
                    !       assigned to the same grid point
    lsytac_d        ! .t.   : if .t.: TAC preferred over BUFR Synop reports
                    !         if .f.: vice versa; this uses DWD DB KZ

  INTEGER        :: &
!   size def. of the 'ODR' (obs. data record) for internal storage of all obs.
!   within the max. time window given by 'wtuk??a', wtuk??e', 'tip??mx'
    maxmlo_d     ,& !  600   : max. number of multi-level reports in the ODR
    maxsgo_d     ,& ! 4000   : max. number of (surface-level and upper-air)
                    !                         single-level reports in the ODR
    maxuso_d     ,& !  900   : max. number of upper-air single-level rep. in ODR
    maxgpo_d     ,& !  200   : max. number of GPS reports in ODR on total domain
    maxmlv_d     ,& !  100   : max. number of observation levels in multi-level
                    !                         reports
!   maxtvo_d     ,&
    nolbc_d      ,& !    5   : number of grid rows at lateral boundaries
                    !          where obs are neglected
    itim_wp_d    ,& !    0   : mode of correction of obs time for wind profiler
                    !          = 0 : no correction
                    !          = 1 : correct by half of the obs averaging period
                    !          = 2 : correct to the end of the obs aver. period
                    !          = 3 : correct to the beginning of the obs period
    icdt_tws_d   ,& !    0   : mode for obs/code type of tower surface-level rep
                    !          = 0 : no correction (remains PILOT/ ICOS tower)
                    !          = +/- 1 : automatic Synop     \ > 0: 
                    !          = +/- 2 : Synop code type 839 /
                    !          > 0: obstype status active if Synop to be active
                    !          < 0: active only if both Synop + tower active
    icdt_rss_d   ,& !    0   : mode for radiosonde-derived surface-level report
                    !          =-1 : no single-level surface report (SR) created
                    !          = 0 : SR: type + status like parent TEMP
                    !          = 1 : SR: obs type SYNOP, code 835, active onlyif
                    !                parent TEMP codetype active and lcd835 true
                    !          = 2 : SR: obs type SYNOP, code 835, status lcd835
    ilocv_sfc_d  ,& !    0   : mode for setting 'plevel' in surface reports
                    !          = 0 : set to observed station pressure only if
                    !                present (else set missing value
                    !                 -> US Std Atmos used in LETKF for v-loc.)
                    !          = 1 : set to model surface pressure instead
                    !                if station pressure not observed
                    !          = 2 : always set to model surface pressure
    mqcorr92_d   ,& !    0   : switch for bias correction for Vaisala RS92
                    !            radiosonde humidity
                    !          = 0 : no correction for humidity
                    !          = 1 : correct only solar radiation bias
                    !          = 2 : correct total bias (incl. nighttime bias)
    igpscen_d (mxgpc),& ! X* -1  : array of used GPS processing centres
    ivtwex_d (mxtwex)   ! X* 0   : variable indicators rel. to exception levels:
                        !          O: none; 1: wind; 2: T + humidity; 3: all

  CHARACTER (LEN=ilstidtw) :: &
    ytwex_d  (mxtwex),& ! X*' '  : station ID's related to exception levels
    ybcrr_d  (mxbcrr)   !        : station ID's for RASS bias correction rules

!      8.2   Use of observation and code types
!      ---------------------------------------

  LOGICAL        :: &
    lsynop_d     ,& ! .t.    : .t. if SYNOP data is used
    laircf_d     ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob_d     ,& ! .t.    : .t. if SATOB data is used
    ldribu_d     ,& ! .t.    : .t. if DRIBU data is used (drifting buoy)
    ltemp_d      ,& ! .t.    : .t. if TEMP  data is used
    lpilot_d     ,& ! .t.    : .t. if PILOT data is used
    lsatem_d     ,& ! .false.: .t. if SATEM data is used
    lgps_d       ,& ! .false.: .t. if GPS   data is used
    lscatt_d        ! .t.    : .t. if SCATT data is used

  LOGICAL        :: &
    lcd011_d     ,& ! .t.    : synop code  11 data is used (land synop)
    lcd014_d     ,& ! .t.    : synop code  14 data is used (automatic)
    lcd021_d     ,& ! .t.    : synop code  21 data is used (ship)
!   lcd022_d     ,& ! .t.    : synop code  22 data is used (ship abbrev.)
!   lcd023_d     ,& ! .t.    : synop code  23 data is used (shred)
    lcd024_d     ,& ! .t.    : synop code  24 data is used (autom. ship)
    lcd140_d     ,& ! .t.    : synop code 140 data is used (metar)
    lcd811_d     ,& ! .t.    : synop code 811 data is used (synop test)
    lcd835_d     ,& ! .t.    : synop code 835 data is used (surface TEMP)
    lcd839_d     ,& ! .t.    : synop code 839 data is used (surface tower)
    lcd041_d     ,& ! .t.    : airep code  41 data is used (codar)
    lcd141_d     ,& ! .t.    : airep code 141 data is used (airep)
!   lcd241_d     ,& ! .t.    : airep code 241 data is used (colba)
    lcd144_d     ,& ! .t.    : airep code 144 data is used (amdar)
    lcd244_d     ,& ! .t.    : airep code 244 data is used (acars)
    lcd146_d     ,& ! .t.    : airep code 146 data is used (mode-s)
    lcd088_d     ,& ! .t.    : satob code  88 data is used (satob)
    lcd090_d     ,& ! .t.    : satob code  90 data is used (amv)
!   lcd188_d     ,& ! .false.: satob code 188 data is used (sst)
!   lcd063_d     ,& ! .t.    : dribu code  63 data is used (bathy)
    lcd064_d     ,& ! .t.    : dribu code  64 data is used (tesac)
    lcd165_d     ,& ! .t.    : dribu code 165 data is used (drift. buoy)
    lcd035_d     ,& ! .t.    : temp  code  35 data is used (land temp)
    lcd036_d     ,& ! .t.    : temp  code  36 data is used (temp ship)
    lcd037_d     ,& ! .t.    : temp  code  37 data is used (mobile)
    lcd135_d     ,& ! .t.    : temp  code 135 data is used (dropsonde)
    lcd109_d     ,& ! .t.    : temp  code 109 data is used (land temp hi-res)
    lcd111_d     ,& ! .t.    : temp  code 111 data is used (ship temp hi-res)
    lcd230_d     ,& ! .t.    : temp  code 230 data is used (drop temp hi-res)
    lcd231_d     ,& ! .t.    : temp  code 231 data is used (desc temp hi-res)
    lcd039_d     ,& ! .t.    : temp  code  39 data is used (rocob)
    lcd040_d     ,& ! .t.    : temp  code  40 data is used (rocob ship)
    lcd032_d     ,& ! .t.    : pilot code  32 data is used (land pilot)
    lcd033_d     ,& ! .t.    : pilot code  33 data is used (pilot ship)
    lcd038_d     ,& ! .t.    : pilot code  38 data is used (mobile)
    lcd132_d     ,& ! .t.    : pilot code 132 data is used (win-prof eu)
    lcd133_d     ,& ! .t.    : pilot code 133 data is used (sod/rass eu)
    lcd136_d     ,& ! .t.    : pilot code 136 data is used (pro/rass us)
    lcd137_d     ,& ! .t.    : pilot code 137 data is used (Radar VAD)
    lcd139_d     ,& ! .t.    : pilot code 139 data is used (tower)
    lcd159_d     ,& ! .f.    : pilot code 159 data is used (icos tower)
    lcd187_d     ,& ! .f.    : pilot code 187 data is used (wind lidar)
    lcd086_d     ,& ! .false.: satem code  86 data is used (satem)
    lcd186_d     ,& ! .false.: atovs code 186 data is used (hi-res ATOVS)
    lcd122_d     ,& ! .t.    : scatt code 122 data is used (QuickScat)
    lcd123_d     ,& ! .t.    : scatt code 123 data is used (ASCAT)
    lcd096_d        ! .t.    : gps data from COST ASCII file is used

! INTEGER        :: &
!   mcdmsg1_d    ,& !  processing / use of MSG1   code  71 data
!   mcdmsg2_d    ,& !  processing / use of MSG2   code  72 data
!   mcdno15_d    ,& !  processing / use of NOAA15 code 206 data
!   mcdno16_d    ,& !  processing / use of NOAA16 code 207 data
!   mcdno17_d    ,& !  processing / use of NOAA17 code 208 data
!   mcdno18_d       !  processing / use of NOAA18 code 209 data

!      10.   Diagnostic output
!      -----------------------
  INTEGER        :: &
    ionl_d         ,& ! 100    : / grid point coordinates
    jonl_d            ! 1      : \ for standard output on nudging

  INTEGER        :: &
    nav            ,& ! number of definition levels for superobbing
    ii             ,& ! loop index
    ierr              ! iostat return variable for namelist input

! Define the namelist group
! -------------------------

 namelist /COSMO_OBS/ ycdfdir, verification_start, verification_end,           &
                      qcc, qccsu, qcvf, doromx, altopsu, zlimv10, dhosag,      &
                      thairh, rhtsat, itim_wp, icdt_tws, icdt_rss, ilocv_sfc,  &
                      fplev_ps,                                                &
                      rtmlrs, rtmlsy, rtmlair, rtmltow, rtmlrsy, rtmlim,       &
                      av_levs, av_incr, av_reso, lredn_repro, lsytac,          &
                      htwex, ivtwex, ytwex, ybcrr, bcrrt, bcrrhl, bcrrhu,      &
                      maxmlo, maxsgo, maxgpo, maxmlv, mqcorr92,                &
                      lsynop, laircf, lsatob, ldribu, ltemp, lpilot, lscatt,   &
                      lcd144, lcd146,         lcd811, lcd835, lcd839,          &
                      lcd035, lcd036, lcd109, lcd111, lcd231,                  &
                      lcd032, lcd132, lcd133, lcd136, lcd137, lcd139, lcd159,  &
                      lcd187,         ionl, jonl, nolbc

!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE read_nml_cosmo_obs
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!-------------------------------------------------------------------------------
 lverpas_d    = .TRUE.
 lcloud_ice_d = .TRUE.

 verification_start_d = - 29
 verification_end_d   =    0

 qcc_d    = (/ 0.0_wp, 500.0_wp,  0.0_wp, 0.7_wp/)
 qccsu_d  = (/12.0_wp, 500.0_wp, 12.0_wp, 0.7_wp/)
 qcvf_d   = (/ 5.0_wp,   1.0_wp, 10.0_wp, 0.0_wp/)
 qcciq_d  = 1.0_wp
 qcsiq_d  = .15_wp

 ycdfdir_d = obsinput

 maxmlo_d =  3000  ! max. number of multi-level reports in the ODR
 maxsgo_d = 20000  ! max. number of (surface-level and upper-air)
 maxuso_d =   900  ! max. number of upper-air single-level rep. in ODR
 maxgpo_d =   200  ! max. number of GPS reports in ODR on total domain
 maxmlv_d =   100  ! max. number of observation levels in multi-level reports
!maxtvo_d =   200  ! max. number of sat retrievals within totaldomain

 nolbc_d = 5
 mqcorr92_d = 0
 igpscen_d = (/30,23,26,24,29,33,34,37,32,0,21,35,25,-1,-1,-1,-1,-1,-1,-1/)

 obnlat_d =   90._wp
 obslat_d =  -90._wp
 obwlon_d = -180._wp
 obelon_d =  180._wp
 exnlat_d =   90._wp
 exslat_d =  -90._wp
 exwlon_d = -180._wp
 exelon_d =  180._wp

 doromx_d   = (/100._wp, 150._wp, 150._wp, 150._wp/)
 altopsu_d  = (/100._wp, 5000._wp, 5000._wp, 5000._wp/)
 zlimv10_d  = (/400.0_wp, 800.0_wp, 5.1_wp/)
!zlimv10_d  = (/75.0_wp, 200.0_wp, 0.7_wp/)
 dhosag_d   = 5._wp
 thairh_d   = 0._wp
 rhtsat_d   = (/0.96_wp, 0.96_wp, 0.96_wp, 0.96_wp, 0.96_wp/)
 itim_wp_d  = 0
 icdt_tws_d = 0
 icdt_rss_d = 0
 ilocv_sfc_d= 0
 rtmlrs_d   = 0.751_wp
 rtmlsy_d   = 0.35_wp
 rtmlair_d  = 0.25_wp
 rtmltow_d  = 0.51_wp
 rtmlrsy_d  = 0.99_wp
 rtmlim_d   = 0.15_wp
 fplev_ps_d = 1.0_wp
 av_levs_d  = 0.0_wp
 av_incr_d  = 0.0_wp
 av_levs_d(1:6) = (/1075._wp, 755._wp, 710._wp, 90._wp, 75._wp, 5._wp/)
 av_incr_d(1:6) = (/  10._wp,  15._wp,  20._wp, 15._wp, 10._wp, 0._wp/)
 av_reso_d      = 3.0_wp
 lredn_repro_d  = .FALSE.
 lsytac_d       = .TRUE.
 htwex_d (1:mxtwex) = -1._wp
 ivtwex_d(1:mxtwex) =  0
 ytwex_d (1:mxtwex) = ' '
 ybcrr_d (1:mxbcrr) = ' '
 bcrrt_d (1:mxbcrr) = 0._wp
 bcrrhl_d(1:mxbcrr) = 0._wp
 bcrrhu_d(1:mxbcrr) = 0._wp
 ybcrr_d (1:4) = (/ '10394', '10394', '10678', '10266'/)
 bcrrt_d (1:4) = (/ -0.2_wp, -0.5_wp, -0.6_wp, -0.6_wp/)
 bcrrhl_d(1:4) = (/   0._wp, 600._wp,   0._wp,   0._wp/)
 bcrrhu_d(1:4) = (/ 600._wp,9000._wp,9000._wp,9000._wp/)

 lsynop_d = .TRUE.
 laircf_d = .TRUE.
 lsatob_d = .TRUE.
 ldribu_d = .TRUE.
 ltemp_d  = .TRUE.
 lpilot_d = .TRUE.
 lsatem_d = .FALSE.
 lgps_d   = .FALSE.
 lscatt_d = .TRUE.

 lcd011_d = .TRUE.
 lcd014_d = .TRUE.
 lcd021_d = .TRUE.
!lcd022_d = .TRUE.
!lcd023_d = .TRUE.
 lcd024_d = .TRUE.
 lcd140_d = .TRUE.
 lcd811_d = .TRUE.
 lcd835_d = .TRUE.
 lcd839_d = .TRUE.
 lcd035_d = .TRUE.
 lcd041_d = .TRUE.
 lcd141_d = .TRUE.
!lcd241_d = .TRUE.
 lcd144_d = .TRUE.
 lcd146_d = .TRUE.
 lcd244_d = .TRUE.
 lcd088_d = .TRUE.
 lcd090_d = .TRUE.
!lcd188_d = .FALSE.
!lcd063_d = .TRUE.
 lcd064_d = .TRUE.
 lcd165_d = .TRUE.
 lcd035_d = .TRUE.
 lcd036_d = .TRUE.
 lcd037_d = .TRUE.
 lcd135_d = .TRUE.
 lcd109_d = .TRUE.
 lcd111_d = .TRUE.
 lcd230_d = .TRUE.
 lcd231_d = .TRUE.
 lcd039_d = .TRUE.
 lcd040_d = .TRUE.
 lcd032_d = .TRUE.
 lcd033_d = .TRUE.
 lcd038_d = .TRUE.
 lcd132_d = .TRUE.
 lcd133_d = .TRUE.
 lcd136_d = .FALSE.
 lcd137_d = .FALSE.
 lcd139_d = .FALSE.
 lcd159_d = .FALSE.
 lcd187_d = .FALSE.
 lcd086_d = .FALSE.
 lcd186_d = .FALSE.
 lcd122_d = .TRUE.
 lcd123_d = .TRUE.
 lcd096_d = .TRUE.
!mcdmsg1_d = 0
!mcdmsg2_d = 0
!mcdno15_d = 0
!mcdno16_d = 0
!mcdno17_d = 0
!mcdno18_d = 0

!ionl_d = 100
 ionl_d =  -1   ! Default: do not open yuprint
 jonl_d = 1

!-------------------------------------------------------------------------------



!- Section 2: Initialize variables with defaults
!-------------------------------------------------------------------------------
  lverpas     = lverpas_d
  lcloud_ice  = lcloud_ice_d

  verification_start = verification_start_d
  verification_end   = verification_end_d

  qcc         = qcc_d
  qccsu       = qccsu_d
  qcvf        = qcvf_d
  qcciq       = qcciq_d
  qcsiq       = qcsiq_d

  ycdfdir     = ycdfdir_d

  doromx      = doromx_d
  altopsu     = altopsu_d
  zlimv10     = zlimv10_d
  dhosag      = dhosag_d
  thairh      = thairh_d
  rhtsat      = rhtsat_d
  itim_wp     = itim_wp_d
  icdt_tws    = icdt_tws_d
  icdt_rss    = icdt_rss_d
  ilocv_sfc   = ilocv_sfc_d
  rtmlrs      = rtmlrs_d
  rtmlsy      = rtmlsy_d
  rtmlair     = rtmlair_d
  rtmltow     = rtmltow_d
  rtmlrsy     = rtmlrsy_d
  rtmlim      = rtmlim_d
  fplev_ps    = fplev_ps_d
  av_levs     = av_levs_d
  av_incr     = av_incr_d
  av_reso     = av_reso_d
  lredn_repro = lredn_repro_d
  lsytac      = lsytac_d
  htwex       = htwex_d
  ivtwex      = ivtwex_d
  ytwex       = ytwex_d
  ybcrr       = ybcrr_d
  bcrrt       = bcrrt_d
  bcrrhl      = bcrrhl_d
  bcrrhu      = bcrrhu_d

  maxmlo      = maxmlo_d
  maxsgo      = maxsgo_d
  maxuso      = maxuso_d
  maxgpo      = maxgpo_d
  maxmlv      = maxmlv_d
  nolbc       = nolbc_d
  mqcorr92    = mqcorr92_d
  igpscen     = igpscen_d

  obnlat      = obnlat_d
  obslat      = obslat_d
  obwlon      = obwlon_d
  obelon      = obelon_d
  exnlat      = exnlat_d
  exslat      = exslat_d
  exwlon      = exwlon_d
  exelon      = exelon_d

  lsynop      = lsynop_d
  laircf      = laircf_d
  lsatob      = lsatob_d
  ldribu      = ldribu_d
  ltemp       = ltemp_d
  lpilot      = lpilot_d
  lsatem      = lsatem_d
  lgps        = lgps_d
  lscatt      = lscatt_d

  lcd011      = lcd011_d
  lcd014      = lcd014_d
  lcd021      = lcd021_d
! lcd022      = lcd022_d
! lcd023      = lcd023_d
  lcd024      = lcd024_d
  lcd140      = lcd140_d
  lcd811      = lcd811_d
  lcd835      = lcd835_d
  lcd839      = lcd839_d
  lcd041      = lcd041_d
  lcd141      = lcd141_d
! lcd241      = lcd241_d
  lcd144      = lcd144_d
  lcd146      = lcd146_d
  lcd244      = lcd244_d
  lcd088      = lcd088_d
  lcd090      = lcd090_d
! lcd188      = lcd188_d
! lcd063      = lcd063_d
  lcd064      = lcd064_d
  lcd165      = lcd165_d
  lcd035      = lcd035_d
  lcd036      = lcd036_d
  lcd037      = lcd037_d
  lcd135      = lcd135_d
  lcd109      = lcd109_d
  lcd111      = lcd111_d
  lcd230      = lcd230_d
  lcd231      = lcd231_d
  lcd039      = lcd039_d
  lcd040      = lcd040_d
  lcd032      = lcd032_d
  lcd033      = lcd033_d
  lcd038      = lcd038_d
  lcd132      = lcd132_d
  lcd133      = lcd133_d
  lcd136      = lcd136_d
  lcd137      = lcd137_d
  lcd139      = lcd139_d
  lcd159      = lcd159_d
  lcd187      = lcd187_d
  lcd086      = lcd086_d
  lcd186      = lcd186_d
  lcd122      = lcd122_d
  lcd123      = lcd123_d
  lcd096      = lcd096_d
! mcdmsg1     = mcdmsg1_d
! mcdmsg2     = mcdmsg2_d
! mcdno15     = mcdno15_d
! mcdno16     = mcdno16_d
! mcdno17     = mcdno17_d
! mcdno18     = mcdno18_d

  ionl        = ionl_d
  jonl        = jonl_d

!- Section 3: Read namelist variables and replace default variables if necessar
!------------------------------------------------------------------------------

  !----------------------------
  ! read on processor p_io only
  !----------------------------
  if (dace% lpio) then
    !----------------------------------
    ! search for namelist group in file
    !----------------------------------
    call position_nml ('COSMO_OBS', status=ierr)
    select case (ierr)
    case (POSITIONED)
    !----------------------------------------------
    ! read namelist, special error handling for IBM
    !----------------------------------------------
#if defined(__ibm__)
      read (nnml ,nml=COSMO_OBS, iostat=ierr)
      if (ierr/=0) call finish ('read_psas_nml','ERROR in namelist /COSMO_OBS/')
#else
      read (nnml ,nml=COSMO_OBS)
#endif
    end select
    !---------------------------------
    ! check values                 ...
    !---------------------------------
    nav = count( av_levs > 0._wp )
    if (     (nav <= 1) .or. (any( av_levs(1:nav-1) <= av_levs(2:nav) ))       &
        .or. (any( av_levs(1:nav)   <= 0._wp ))                                &
        .or. (any( av_incr(1:nav-1) <= 0._wp ))) then
      ! if values are invalid, set them to default
      av_levs = av_levs_d
      av_incr = av_incr_d
    endif
    ! hPa --> Pa  (values in diagnostic print below are in hPa)
    av_levs = av_levs_d * 100._wp
    av_incr = av_incr_d * 100._wp
    if (rhtsat(2) == rhtsat_d(2))  rhtsat(2) = rhtsat(1)
    if (rhtsat(3) == rhtsat_d(3))  rhtsat(3) = rhtsat(1)
    if (rhtsat(4) == rhtsat_d(4))  rhtsat(4) = rhtsat(1)
    if (rhtsat(5) == rhtsat_d(5))  rhtsat(5) = rhtsat(1)
    ! 'exception levels' of single tower stations
    ivtwex = max( 0, min( 3, ivtwex ) )
    ntwex  = 0
    tower_exception: do ii = 1, mxtwex
      if (     (htwex(ii) < -0.1_wp) .or. (ivtwex(ii) <= 0)                    &
          .or. (ytwex(ii)(1:1) == ' '))                     exit tower_exception
      ntwex = ii
    enddo tower_exception
    if (ntwex == 0) then
      htwex  = htwex_d
      ivtwex = ivtwex_d
      ytwex  = ytwex_d
    endif
    ! RASS bias correction rules
    nbcrr = 0
    rass_bias_correction_rules: do ii = 1, mxbcrr
      if (     (min( bcrrhl(ii), bcrrhu(ii) ) < -400.1_wp)                     &
          .or. (bcrrhl(ii) > bcrrhu(ii)) .or. (abs( bcrrt(ii) ) >  3.01_wp)    &
                                         .or. (abs( bcrrt(ii) ) <= 0.01_wp)    &
          .or. (ybcrr(ii)(1:1) == ' '))          exit rass_bias_correction_rules
      nbcrr = ii
    enddo rass_bias_correction_rules
    do ii = nbcrr+1 , mxbcrr
      ybcrr  (ii) = ' '
      bcrrt  (ii) = 0._wp
      bcrrhl (ii) = 0._wp
      bcrrhu (ii) = 0._wp
    enddo

    !---------------------------------
    ! consistency checks, printout ...
    !---------------------------------
    ycdfdir = path_file (obsinput, ycdfdir)
    write(6,*)
    write(6,*) 'namelist /COSMO_OBS/ read'
    write(6,*)
    write(6,*) '  ycdfdir            = ',trim(ycdfdir)
    write(6,*) '  verification_start = ',verification_start
    write(6,*) '  verification_end   = ',verification_end
    write(6,*) '  qcc                = ',qcc
    write(6,*) '  qccsu              = ',qccsu
    write(6,*) '  qcvf               = ',qcvf
    write(6,*) '  doromx             = ',doromx
    write(6,*) '  altopsu            = ',altopsu
    write(6,*) '  zlimv10            = ',zlimv10
    write(6,*) '  dhosag             = ',dhosag
    write(6,*) '  thairh             = ',thairh
    write(6,*) '  rhtsat             = ',rhtsat
    write(6,*) '  itim_wp            = ',itim_wp
    write(6,*) '  icdt_tws           = ',icdt_tws
    write(6,*) '  icdt_rss           = ',icdt_rss
    write(6,*) '  ilocv_sfc          = ',ilocv_sfc
    write(6,*) '  rtmlrs             = ',rtmlrs
    write(6,*) '  rtmlsy             = ',rtmlsy
    write(6,*) '  rtmlair            = ',rtmlair
    write(6,*) '  rtmltow            = ',rtmltow
    write(6,*) '  rtmlrsy            = ',rtmlrsy
    write(6,*) '  rtmlim             = ',rtmlim
    write(6,*) '  fplev_ps           = ',fplev_ps
    write(6,*) '  av_levs            = ',pack (av_levs *0.01_wp, av_levs > 0._wp)
    write(6,*) '  av_incr            = ',pack (av_incr *0.01_wp, av_incr > 0._wp)
    write(6,*) '  av_reso            = ',av_reso
    write(6,*) '  lredn_repro        = ',lredn_repro
    write(6,*) '  lsytac             = ',lsytac
    write(6,*) '  ytwex              = ',(ytwex (ii), ii=1,max(ntwex,1))
    write(6,*) '  htwex              = ',(htwex (ii), ii=1,max(ntwex,1))
    write(6,*) '  ivtwex             = ',(ivtwex(ii), ii=1,max(ntwex,1))
    write(6,*) '  ybcrr              = ',(ybcrr (ii), ii=1,max(nbcrr,1))
    write(6,*) '  bcrrt              = ',(bcrrt (ii), ii=1,max(nbcrr,1))
    write(6,*) '  bcrrhl             = ',(bcrrhl(ii), ii=1,max(nbcrr,1))
    write(6,*) '  bcrrhu             = ',(bcrrhu(ii), ii=1,max(nbcrr,1))
    write(6,*) '  maxmlo             = ',maxmlo
    write(6,*) '  maxsgo             = ',maxsgo
    write(6,*) '  maxgpo             = ',maxgpo
    write(6,*) '  maxmlv             = ',maxmlv
    write(6,*) '  mqcorr92           = ',mqcorr92
    write(6,*) '  lsynop             = ',lsynop
    write(6,*) '  laircf             = ',laircf
    write(6,*) '  lsatob             = ',lsatob
    write(6,*) '  ldribu             = ',ldribu
    write(6,*) '  ltemp              = ',ltemp
    write(6,*) '  lpilot             = ',lpilot
    write(6,*) '  lscatt             = ',lscatt
    write(6,*) '  lcd144             = ',lcd144
    write(6,*) '  lcd146             = ',lcd146
    write(6,*) '  lcd035             = ',lcd035
    write(6,*) '  lcd036             = ',lcd036
    write(6,*) '  lcd109             = ',lcd109
    write(6,*) '  lcd111             = ',lcd111
    write(6,*) '  lcd231             = ',lcd231
    write(6,*) '  lcd032             = ',lcd032
    write(6,*) '  lcd132             = ',lcd132
    write(6,*) '  lcd133             = ',lcd133
    write(6,*) '  lcd136             = ',lcd136
    write(6,*) '  lcd137             = ',lcd137
    write(6,*) '  lcd139             = ',lcd139
    write(6,*) '  lcd159             = ',lcd159
    write(6,*) '  lcd187             = ',lcd187
    write(6,*) '  lcd811             = ',lcd811
    write(6,*) '  lcd835             = ',lcd835
    write(6,*) '  lcd839             = ',lcd839
    write(6,*) '  ionl               = ',ionl
    write(6,*) '  jonl               = ',jonl
    write(6,*) '  nolbc              = ',nolbc
    write(6,*)
  endif
  !----------------------------
  ! broadcast to all processors
  !----------------------------
  call p_bcast (ycdfdir            ,dace% pio)
  call p_bcast (verification_start ,dace% pio)
  call p_bcast (verification_end   ,dace% pio)
  call p_bcast (qcc                ,dace% pio)
  call p_bcast (qccsu              ,dace% pio)
  call p_bcast (qcvf               ,dace% pio)
  call p_bcast (doromx             ,dace% pio)
  call p_bcast (altopsu            ,dace% pio)
  call p_bcast (zlimv10            ,dace% pio)
  call p_bcast (dhosag             ,dace% pio)
  call p_bcast (thairh             ,dace% pio)
  call p_bcast (rhtsat             ,dace% pio)
  call p_bcast (itim_wp            ,dace% pio)
  call p_bcast (icdt_tws           ,dace% pio)
  call p_bcast (icdt_rss           ,dace% pio)
  call p_bcast (ilocv_sfc          ,dace% pio)
  call p_bcast (rtmlrs             ,dace% pio)
  call p_bcast (rtmlsy             ,dace% pio)
  call p_bcast (rtmlair            ,dace% pio)
  call p_bcast (rtmltow            ,dace% pio)
  call p_bcast (rtmlrsy            ,dace% pio)
  call p_bcast (rtmlim             ,dace% pio)
  call p_bcast (fplev_ps           ,dace% pio)
  call p_bcast (av_levs            ,dace% pio)
  call p_bcast (av_incr            ,dace% pio)
  call p_bcast (av_reso            ,dace% pio)
  call p_bcast (lredn_repro        ,dace% pio)
  call p_bcast (lsytac             ,dace% pio)
  call p_bcast (maxmlo             ,dace% pio)
  call p_bcast (maxsgo             ,dace% pio)
  call p_bcast (maxgpo             ,dace% pio)
  call p_bcast (maxmlv             ,dace% pio)
  call p_bcast (mqcorr92           ,dace% pio)
  call p_bcast (lsynop             ,dace% pio)
  call p_bcast (laircf             ,dace% pio)
  call p_bcast (lsatob             ,dace% pio)
  call p_bcast (ldribu             ,dace% pio)
  call p_bcast (ltemp              ,dace% pio)
  call p_bcast (lpilot             ,dace% pio)
  call p_bcast (lscatt             ,dace% pio)
  call p_bcast (lcd144             ,dace% pio)
  call p_bcast (lcd146             ,dace% pio)
  call p_bcast (lcd035             ,dace% pio)
  call p_bcast (lcd036             ,dace% pio)
  call p_bcast (lcd109             ,dace% pio)
  call p_bcast (lcd111             ,dace% pio)
  call p_bcast (lcd231             ,dace% pio)
  call p_bcast (lcd032             ,dace% pio)
  call p_bcast (lcd132             ,dace% pio)
  call p_bcast (lcd133             ,dace% pio)
  call p_bcast (lcd136             ,dace% pio)
  call p_bcast (lcd137             ,dace% pio)
  call p_bcast (lcd139             ,dace% pio)
  call p_bcast (lcd159             ,dace% pio)
  call p_bcast (lcd187             ,dace% pio)
  call p_bcast (lcd811             ,dace% pio)
  call p_bcast (lcd835             ,dace% pio)
  call p_bcast (lcd839             ,dace% pio)
  call p_bcast (ionl               ,dace% pio)
  call p_bcast (jonl               ,dace% pio)
  call p_bcast (nolbc              ,dace% pio)
  call p_bcast (ntwex              ,dace% pio)
  if (ntwex > 0) then
    call p_bcast (ytwex (1:ntwex)  ,dace% pio)
    call p_bcast (htwex (1:ntwex)  ,dace% pio)
    call p_bcast (ivtwex(1:ntwex)  ,dace% pio)
  endif
  call p_bcast (nbcrr              ,dace% pio)
  if (nbcrr > 0) then
    call p_bcast (ybcrr (1:nbcrr)  ,dace% pio)
    call p_bcast (bcrrt (1:nbcrr)  ,dace% pio)
    call p_bcast (bcrrhl(1:nbcrr)  ,dace% pio)
    call p_bcast (bcrrhu(1:nbcrr)  ,dace% pio)
  endif

!print *, 'ZZ ', dace% pe, lcd146, lcd133, lcd136, lcd137, lcd139, lcd159

END SUBROUTINE read_nml_cosmo_obs


END MODULE mo_cosmo_obs_data
