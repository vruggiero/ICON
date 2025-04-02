!+ Source module for special observation processing of aircraft observations
!-------------------------------------------------------------------------------

MODULE src_obs_proc_air

!-------------------------------------------------------------------------------
! Description:
!   This module performs observation processing steps that are applied to
!   aircraft observations only. This includes:
!    - production of vertical profiles from single-level aircraft reports.
!    - flight track check and thinning
!    - adaptation of vertical correlation scales for dense sequences of reports
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_air_org_mult  : organizing production of multi-level aircraft reports
!    - obs_air_make_mult : production of multi-level aircraft reports
!    - obs_air_correl    : shorten vertical correlation scale for dense reports
!    - obs_air_list      : lists of aircraft reports with specified requirements
!    - obs_air_interface : interface to model environment if AOF-read
!   Driver routine 'obs_air_org_mult' is called in 'obs_org_cdf_proc' of module
!   src_obs_proc_cdf.f90. Other procedures are called in 'obs_air_org_mult'.
!   'obs_air_interface' is called in module 'src_obs_processing.f90' only if
!   observation data are read from AOF file.
!
!   It uses from:
!    - src_obs_cdfin_util:    - obs_assign_sort_node
!    - parallel_utilities:    - global_values
!                             - gather_all
!                             - reports_all2all
!    - environment:           - model_abort
!                           \ - comm_barrier
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin
!    - data_obs_record
!    - mo_fdbk_tables
!
!   Note: This module contains MPI-based communication between PE's.
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.12       2004/09/15 Christoph Schraff
!  Initial release
! 3.18       2006/03/03 Christoph Schraff
!  Caution messages on insufficient ODR size to specific file unit 'nucautn'.
! V4_4         2008/07/16 Ulrich Schaettler
!  Adapted interface of get_timings
! V4_5         2008/09/10 Christoph Schraff
!  Adaptions to read observations from NetCDF files (changed report number and
!  index names, additional ODR elements)
! V4_8         2009/02/16 Oliver Fuhrer
!  Use global values only if num_compute greater 1
! V4_9         2009/07/16 Ulrich Schaettler
!  Replaced some quotes with double quotes (problems with IBM preprocessor)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Modified interface so that the only data modules used are 'data_parameters'
!    'data_obs_cdfin', 'data_obs_lib_cosmo', and 'data_obs_record', whereas
!    'data_nudge_all', 'data_modelconfig', 'data_runcontrol', 'data_parallel,
!    and 'data_obs_process' are not used any more. The additional interface
!    routine 'obs_air_interface' is used only when the obs are read from AOF.
!    'i_subpos' is included in the argument list of routine 'obs_air_org_mult'.
!  - 'get_timings', writing rejection messages, and writing of messages on
!    number of reports exceeding array size is shifted to outside this module.
!  - Some array sizes are re-defined (and more safe to avoid violation of array
!    bounds).
!  - More consistent setting of ODR report flags and updating of counters for
!    statistics.
! V4_28        2013/07/12 Christoph Schraff
!  - Adaptation for adding the reading from feedobs files: Reduce 'ntotmlo' only
!    by cancelled reports that have really been read at previous timesteps.
!    (This is achieved by using the new variable 'ndifo'.)
!  - Direct calls of mpi routines (mpi_alltoallv) replaced by calls of a new
!    routine in module 'parallel_utilities' ('reports_all2all').
!  - Statement functions replaced by intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Processing of Mode-S aircraft obs introduced (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Removal of (all dependencies on) the AOF interface.
!  For VOF: if blacklist flag set, set passive bit in station characteristics.
! V5_4d        2016-12-12 Christoph Schraff
!  For Mode-S data (for which the data with the same 'station id' are not
!  usually produced by one flight), the determination of the flight phase(s),
!  the flight track check, the production of multi-level reports, and the
!  thinning are skipped.
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Editorial changes only.
! V5_5a        2018-06-22 Ulrich Schaettler
!  Use kind_parameters instead of data_parameters (US)
!  Removed all iintegers (US)
!  Removed irealgrib, which is not used here
! V5_6b        2019-10-16 Christoph Schraff + Andreas Rhodin
!  Entries 'nbtplv', 'nbsplv' added.
!  Set subroutines (obs_air_org_mult) used by other modules explicitly to public
!  (required by ICON Compiler directives when modules are ported to ICON).
! @VERSION@  @DATE@     Christoph Schraff
!  - Stricter criteria for vertical colocation to account for active use of
!    (correct) short-distance flights below 350 hPa throughout.
!  - Meaning of entry 'nbtzio' changed from grid pt. to geographic coordinates.
!  - 16.05.2022: 'nhssos' added.
!
! ==============================================================================
! DACE history:
! ------------
!  1.51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
! ==============================================================================
! CAUTION: This module is used by both the DACE and COSMO model codes.       !!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :   &

! 1. General parameters
! ---------------------

    c0          ,& ! standard real constant 0.0
    c1          ,& ! standard real constant 1.0
    c2          ,& ! standard real constant 2.0
    c05         ,& ! standard real constant 0.5
!   c3600       ,& ! standard real constant 3600.0
    rmdi        ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich      ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy        ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
!                                   or by calling 'obs_air_interface'
! -------------------------------------------------------------------

    ! variables related to parallelisation / domain decomposition
    num_compute ,& ! number of compute PEs
    nboundlines ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id  ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart  ,& ! communicator for the virtual cartesian topology
    imp_reals   ,& ! REAL      type used for MPI
    imp_integers,& ! INTEGER   type used for MPI

    ! report dimensions of ODR
    maxmlv      ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    maxmll      ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl      ,& ! size (report dimension) of the single-level (s-l)  ODR

    ! constants for the horizontal rotated grid and related variables
    r_degrad    ,& ! factor for transforming degree to rad

    ! constants, variables related to namelist parameters, and other variables
    acthr       ,& ! actual model time [hours] with respect to 'yydate_ref'

    ! switches related to namelist parameters and other
    lwonl       ,& ! .true. for the node (sub-domain) at which file with
                   !        the unit number 'nupr' is open

! 3. pressure dependent scales and geometry of horizontal correlations
! --------------------------------------------------------------------

    ncolev      ,& ! number of levels in the correlation scale tables
    tabcolp     ,& ! LOG( tabcop )
    rhvsond     ,& ! upper-air wind horizontal correlation scales
                   ! (pressure dependent part)
    rhtsond     ,& ! upper-air temperature horiz. correlation scales
    rhqsond     ,& ! upper-air humidity horiz. correlation scales

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    nupr           ! all the remaining information

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nairep     ,& ! aircraft reports
    naircd     ,& ! aircraft report
    ncodar     ,& ! codar report
    namdar     ,& ! AMDAR report
    nacar      ,& ! ACAR  report
    nmodes     ,& ! Mode-S report

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

!   t_cmatyp   ,& ! data type for information on CMA observation and code types
    cma        ,& ! array of meta data on CMA observation and code types

! 7. Functions
! ------------

    i_cma         ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         3.1.1   Report event counter array format
!                 ---------------------------------
    neblak     ,& ! blacklisted ship
    netrac     ,& ! (flight) track suspicious
    nethin     ,& ! thinning of aircraft (flight) track
    neslml     ,& ! one multi-level report made from  single-level reports
    neslps     ,& ! single-level report put in multi-level rep. and set passive
    nenoml     ,& ! multi-levl report not made due to ODR array size

!         3.2     Event counter arrays
!                 --------------------
    neventr    ,& ! counter array of report events

!         5       Different kinds of limits
!                 -------------------------
!   rtmlair    ,& !  time limit for reports of obs type 'airep' [hrs]
                  !  (time of lowest level of multi-level ODR)
!   rhzlair    ,& !  horizont. distance limit for airep reports  [km]
    rvtlair       !  vertical  distance limit for airep reports [hpa]

USE data_obs_cdfin, ONLY :  &

!         6       Variables used for producing aircraft multi-level reports
!                 ---------------------------------------------------------

    minslml    ,& !  min. no. of levels in multi-level aircraft report
    thairt     ,& !  maximum temporal distance [h] between the lowest
                  !       report and any single level report that is added to
                  !       a multi-level AIRCRAFT report
    thairv     ,& !  maximum vertical distance [Pa] between two successive
                  !       levels within a multi-level AIRCRAFT report
    maxaid     ,& ! max. number of single-level aircraft reports with the
                  ! same station id.
    nairls     ,& ! length of list 'iairls' (see below)
    mzmlbd     ,& ! body of temporary aircraft multi-level ODR
    mzmlhd     ,& ! header of temporary aircraft multi-level ODR
    mzsgbd     ,& ! body of temporary aircraft single-level ODR
    mzsghd     ,& ! header of temporary aircraft single-level ODR
    mzslon     ,& ! list for single-level reports, as part of the temporary
                  ! single-level ODR, containing the 'long' list indices
                  ! (see below) of the single-level reports
    mzmlon     ,& ! list for multi -level reports, as part of the temporary
                  ! multi -level ODR, containing the 'long' list indices
                  ! (see below) of all single-level reports, which the
                  ! multi-level reports are made of
    iarlls     ,& ! 'long' list, containing all aircraft single-level
                  ! reports, that are processed at the current timestep
    iairls     ,& ! other lists with aircraft single-level report meeting
                  ! various requirements
    ismls      ,& ! list for multi-level reports, that contains the report
                  ! index in the temporary multi-level report array of the
                  ! reports with the current station id
    zmlbdy     ,& ! body of temporary aircraft multi-level ODR
    zmlhed     ,& ! header of temporary aircraft multi-level ODR
    zsgbdy     ,& ! body of temporary aircraft single-level ODR
    zsghed     ,& ! header of temporary aircraft single-level ODR
    rarlls     ,& ! model layer thickness at reports in 'long' list
    rairls     ,& ! model layer thickness at reports in other list 1
    yarlls     ,& ! station id's in 'long' list
    yairls     ,& ! station id's in other lists

!         7       To report rejection of data: Output buffer, size and formats
!                 ------------------------------------------------------------
    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& !  actual number of records stored in the output buffer
    nmxoln     ,& !  maximum length of output buffer
    istrej     ,& !  length of strings (station id) in output buffer
    nfmt17     ,& ! thinning of aircraft reports
    nfmt18     ,& ! exaggerated flight colocation
    nfmt19     ,& ! flight track error
    nfmt26        ! suspicious aircraft identity

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.1    ODR header format
!       ------------------------

!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------
    mxrhed,     & ! header length of multi-level reports
    mxshed,     & ! header length of single-level reports
    nhilon,     & ! longitude of observing station
    nhjlat,     & ! latitude  of observing station
    nhalt ,     & ! station altitude [m]
    nhtime,     & ! (exact) time of observation in forecast hours
    nhsurf,     & ! height of model grid pt. to which obs. is assigned
    nhzio ,     & ! longitude of obs. station (or lowest datum) in grid pt. unit
    nhzjo ,     & ! latitude  of obs. station in grid pt. units
    nhsynt,     & ! nominal (synoptic) time of observation in forecast hours
    nhtddb,     & ! data base decoding time in forecast hours
    nhsolz,     & ! solar zenith angle [deg]
    nhssos,     & ! standard devation of sub-grid scale model orogrophy [m]
    nhvcbu,     & ! correction factor to vertical correlation scale for wind
                  ! at base of report
    nhvcbt,     & ! as 'nhvcbu', but for temperature
    nhvcbq,     & ! as 'nhvcbu', but for humidity
    nhvctu,     & ! correction factor to vertical correlation scale for wind
                  ! at top of report
    nhvctt,     & ! as 'nhvctu', but for temperature
    nhvctq,     & ! as 'nhvctu', but for humidity
    nhtvip,     & ! obs multi-level pressure interpol. to the lowest model level
    nhtviz        ! vertical distance to nearest observation

USE data_obs_record, ONLY :   &

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------
    mxrhdf,     & ! header length of multi-level reports
    mxshdf,     & ! header length of single-level reports
    nhio  ,     & ! (local) x-coord. of grid pt. assigned to obs
    nhjo  ,     & ! (local) y-coord. of grid pt. assigned to obs
    nhitot,     & ! global x-coord. of grid pt. assigned to obs
    nhjtot,     & ! global y-coord. of grid pt. assigned to obs
    nhobtp,     & ! observation type
    nhcode,     & ! code type
    nhschr,     & ! station characteristics (packed as in VOF, see below)
    nhpass,     & ! flag for report being set to 'passive' (as in VOF)
    nhqcfw,     & ! threshold quality control flags for pressure, status of
                  ! verification output
    nhflag,     & ! report flags (obs type, surface, altitude, station ID)
    nhcorr,     & ! update sequence number (station correction indicator)
    nhcat ,     & ! data     category (from BUFR Section 1)
    nhcats,     & ! data sub-category (from BUFR Section 1)
    nhkz  ,     & ! DWD internal classification number (observation type)
    nhcent,     & ! originating centre
    nhstid,     & ! station identity number
    nhdate,     & ! absolute exact observation date [yyyymmdd]
    nhhrmn,     & ! absolute exact observation time [hhmm]
    nhsyhr,     & ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev,     & ! number of obs. levels (for multi-level reports)
    nhuexi,     & ! flag for existence of wind data        in multi-level report
    nhaexi,     & ! flag for exist. of wind or temperature in multi-level report
    nhtexi,     & ! flag for existence of temperature data in multi-level report
    nhqexi,     & ! flag for existence of humidity data    in multi-level report
    nhrtyp,     & ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhtrac,     & ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad ,     & ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4 ,     & ! instrument type                   (NA4, BUFR Table 002003)
    nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
    nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]


!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------

    ilstid,     & ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
!   nvpabp,     & ! bit pos. for report set passive since it is     nvschr
                  !              used in a multi-level pseudo report
    nvpsbp,     & ! bit pos. for report set passive since 1 of next   "
                  !              5 flags or flight track flag applies
    nvbkbp,     & ! bit pos. for flag: 'blacklisted station (ship)'   "
    nvhtbp,     & ! bit pos. for flight track error flag              "
!   nvhtoc,     & ! no. of bits occ. by flight track error flag       "
    nvhhbp,     & ! bit pos. for flight thinning flag                 "
!   nvhhoc,     & ! no. of bits occ. by flight thinning flag          "
    nvapbp,     & ! bit pos. for phase of flight (aircraft)           "
    nvapoc        ! no. of bits occ. by phase of flight               "
!   nvaabp,     & ! bit pos. for aircraft roll angle (code)           "
!   nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd'
!              ------------------------------------------------------
!   maxrsl,     & ! max. number of levels in multi-level ODR
    maxarl,     & ! max. number of levels for multi-level aircraft reports

!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
    mxrbdy,     & ! body length of multi-level reports
    nbtu  ,     & ! u wind component [m/s]
    nbtv  ,     & ! v wind component [m/s]
    nbtt  ,     & ! temperature [K]
    nbtrh ,     & ! relative humidity [/]
    nbtp  ,     & ! pressure [Pa]
    nbtz  ,     & ! height [m]
    nbtuer,     & ! error of observed wind component
    nbtter,     & ! error of observed temperature
    nbtqer,     & ! error of observed rel. humidity
    nbtzer,     & ! error of observed height
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control
    nbtzio,     & ! longitude in geographical coordinates
    nbtzjo,     & ! latitude  in geographical coordinates
    nbttim,     & ! observation time relative to report (header) time
    nbtlop,     & ! LOG( pressure )
    nbtplv,     & ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbtdrh,     & ! bias correction for relative humidity [/]
    nbtw  ,     & ! vertical velocity [m/s]
    nbtsnr,     & ! signal to noise ratio
    nbtuac,     & ! accuracy (std dev from data provider) of horiz. wind [m/s]

!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
    mxrbdf,     & ! body length of multi-level reports
    nbtflg,     & ! main flag word (bit pattern as in VOF body, word 'nvbmfw')
    nbterr,     & ! status flag word        (bit pattern, see below: 'nb?err')
    nbtqcf,     & ! threshold quality control flags
    nbtlsg,     & ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level id (bit pattern, as in VOF)

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
    mxsbdy,     & ! body length of single-level reports
    nbsu  ,     & ! u wind component                                   [m/s]
    nbsv  ,     & ! v wind component                                   [m/s]
    nbst  ,     & ! temperature                                        [K]
    nbsrh ,     & ! relative humidity                                  [/]
    nbsp  ,     & ! pressure                                           [Pa]
    nbsz  ,     & ! height                                             [m]
    nbsuer,     & ! error of observed wind component
    nbster,     & ! error of observed temperature
    nbsqer,     & ! error of observed relative humidity
    nbszer,     & ! error of observed height
    nbsplv,     & ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbsdrh,     & ! bias correction for relative humidity [/]

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
    mxsbdf,     & ! body length of single-level reports
    nbsflg,     & ! main flag word (bit pattern as in VOF body, word 'nvbmfw')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbsqcf,     & ! threshold quality control flags (as in VOF)
    nbslid,     & ! pressure code (SYNOP) or (else) level id. (bit pattern, VOF)
!   nbstur,     & ! degree of turbulence (aircraft, WMO Table 011031)

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status 'active' for horiz. wind    nb?err
    nvrt       ,& ! bit pos. for status 'active' for temperature      "
    nvrq          ! bit pos. for status 'active' for humidity         "
!   nvrz       ,& ! bit pos. for status 'active' for pressure/height  "
!   nvrw          ! bit pos. for status 'active' for vertical wind    "

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi  ,     & ! missing data indicator for ODR integers (2^31-1)
    ntotml,     & ! tot. number of stored multi-level reports
    ntotsg,     & ! tot. number of stored single-level reports

!       2.     Observation data records (ODR)
!       -------------------------------------

    omlbdy,     & ! body   of multi-level ODR
    omlhed,     & ! header of multi-level ODR
    osgbdy,     & ! body   of single-level ODR
    osghed,     & ! header of single-level ODR
    momlbd,     & ! body   of multi-level ODR
    momlhd,     & ! header of multi-level ODR
    mosgbd,     & ! body   of single-level ODR
    mosghd,     & ! header of single-level ODR
    yomlhd,     & ! header of multi-level ODR
    yosghd        ! header of single-level ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_BLACKLIST    ,& ! blacklist (or not whitelist)  (in particular:
                       !   bad (hard blacklisted) aircraft station identity)
    FL_HEIGHT       ,& ! location not in valid height range
    FL_REDUNDANT    ,& ! redundant report
    FL_FLIGHTTRACK  ,& ! flight track error flag
    FL_MERGE        ,& ! (report used for) merged reports (only)
                       !   (e.g. single-level aircraft used for multi-level rep)
    FL_THIN         ,& ! thinning
    LS_SIGN            ! significant level

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
!   comm_barrier,    & ! explicit synchronization point
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values,   & ! computes global values by operating on local arrays
    gather_all,      & ! gathers at all nodes a set of arrays from all nodes
    reports_all2all    ! Individual (amount of) info sent/received by each node

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_sort_node   ! assign node to reports and sort them accordingly

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_air_org_mult

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "obs_air_org_mult.incf"
! INCLUDE "obs_air_make_mult.incf"
! INCLUDE "obs_air_correl.incf"
! INCLUDE "obs_air_list.incf"


!-------------------------------------------------------------------------------
!+ Module procedure for organizing production of multi-level aircraft reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_air_org_mult ( lsvcorl , ntotsgo , ntotmlo , thairh             &
                            , lfirst  , i_subpos , nexceair                    &
                            , odp , vscale_in , tscale_in , hscale_c_in        &
                            , hscale_v_in )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure organizes and performs the production of multi-level
!   aircraft reports from single-level reports with the same flight number.
!
! Method:
!   1. Making of a 'long' list containing the ODR index, code type, source node
!      station id, and model thickness at the obs. location of (globally) all
!      single-level aircraft reports to be processed, i.e. active reports that
!      have been read at the current timestep and reports with passive flag less
!      than two that have been read earlier but have a station id identical to
!      that of a report read at the current timestep.
!      Making of a list of the station id's that occur in the 'long' list.
!      Removal of all multi-level reports with a station id in that list.
!      (Note, that both lists are global and known by each node.)
!   2. Determination of the target (host) nodes of the elements in the 'long'
!      list so that all reports with one identical station id are hosted by one
!      target node.
!      Selective scattering and gathering of the denoted single-level reports
!      by 'mpi_alltoallv' (by call of reports_all2all) to keep communication and
!      memory demands as small as possible.
!   3. Production of multi-level reports (see subr. 'OBS_AIR_MAKE_MULT') by the
!      target nodes for each station id separately. Preparation to set single-
!      level reports to passive which are used in the multi-level reports,
!      or which are rejected by the flight track check.
!   4. Check if the number of multi-level reports exceeded the size of the
!      local ODR and probably the global observation increment arrays.
!      If required, conversion of the multi-level reports with the smallest
!      number of vertical levels back to single-level reports.
!   5. Computation of (reductional) factors to the vertical correlation scales
!      for closely located reports (see subr. 'OBS_AIR_CORREL').
!   6. Selective scattering of these factors and of the passive bit of the
!      single-level reports by the host nodes to the source nodes by using
!      'mpi_alltoallv' (by call of reports_all2all).
!   7. Selective scattering of the complete produced multi-level aircraft
!      reports by the host nodes to the nodes which cover the local domain
!      containing the observation location of the multi-level reports.
!   Note, that only steps 3 and 5 do not contain any communication between PE's.
!   Note, that for passive single-level reports, the 2nd variable of the lists
!         used here is not set equal to the code type, but set to:
!         -1 , if the report is outside the user-specified observation area
!              ==> no change to passive flag in ODR, but used in flight tracking
!         -2 , if the report is rejected by the flight track check
!              ==> set passive flag in ODR to 2, and set flight track check flag
!         After 'obs_air_make_mult', it is set to 0 (in the long list 'iarlls'),
!                if the report is used in a multi-level report
!                ==> set passive flag in ODR to 1
!         Redundant reports (or with inconsistency to orography) are not put to
!         any list.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  Revised subroutine argument list to 'gather_all'. ANSI violations removed.
! 1.27       1999/03/29 Christoph Schraff
!  32-bit missing data indicator in the ODR.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO
! 1.31       1999/07/01 Christoph Schraff
!  Quantities related to MPI communicator 'icomm_world' replaced; input for
!  calls of 'global_values' adjusted.
! 1.34       1999/12/10 Ulrich Schaettler
!  Changed interfaces to timing routines
! 1.39       2000/05/03 Ulrich Schaettler
!  Included dt to call of timing routines.
! 2.5        2001/06/01 Christoph Schraff
!  Handling and setting to passive of reports rejected by the flight track check
! 2.6        2001/06/12 Christoph Schraff
!  Bug correction at report events counting.
! 2.12       2001/11/07 Christoph Schraff
!  Bug correction: 'ibufls' to 'ybufls' deallocated also without aircraft data.
! 2.13       2002/01/18 Christoph Schraff
!  Flight track check and thinning flags moved to station characteristcs word.
! 2.14       2002/02/15 Christoph Schraff
!  Cancellation of suspicious (default) aircraft id's.
! 3.6        2003/12/11 Christoph Schraff
!  Ensuring error message at model_abort.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments:
! --------------------

  LOGICAL                 , INTENT (IN)         ::  &
    lsvcorl          ,& ! adjustment of vertical correlation scales
    lfirst              ! .TRUE when obs processing is called the first time

  INTEGER                 , INTENT (IN)         ::  &
    ntotsgo          ,& ! number of single-level reports before having read new
                        ! data at current timestep
    i_subpos (0:,:)     ! positions of the subdomains in the total domain
                        ! (i-, j-indices of the lower left + upper right grid pt
                        !  in the order: i_ll, j_ll, i_ur, j_ur; only the domain
                        !  interior is considered, not the boundary lines)
                        ! (dimension: (0:npe-1,4) , npe = number of processors)

  INTEGER                 , INTENT (INOUT)      ::  &
    ntotmlo             ! number of multi-level reports before having read new
                        ! data at current timestep
                        !   (needed only for aircraft reports and if (lsvcorl) ,

  INTEGER                 , INTENT (OUT)        ::  &
    nexceair     (2)    ! number of multi-level reports in excess of ODR array

  REAL    (KIND=wp   ),     INTENT (IN)         ::       &
    thairh              ! maximum horizontal distance for combining single level
                        ! aircraft reports to a multi-level report

!   the following variables are used only if (lsvcorl) ,
!   to compute an adapted short vertical correlation scale

  REAL    (KIND=wp   ),     INTENT (IN)   , OPTIONAL  ::  &
    vscale_in    (4) ,& ! scale (in [ln(p)] units) of the static Gaussian
                        ! vertical correlation function
                        ! --> usual and reasonable values: 0.577,0.577, 0.2, 0.2
    tscale_in        ,& ! scale [hrs] of the temporal Gaussian weight
    hscale_c_in  (4) ,& ! vertically invariant part [km]  \  of the scale of
    hscale_v_in  (4) ,& ! fraction of vertically variable  > the autoregressive
                        ! part given by tables 'rh?sond'  /  horizontal weight
                        ! ! where the total Gaussian weight (product of temporal
                        ! ! and horizontal weight) is given to the adapted
                        ! ! short correlation scale (which adapts to the
                        ! ! vertical distance between aircraft reports)
                        ! ! (1 minus the Gaussian weight is given to the large
                        ! !  static correlation scale)
                        ! --> usual and reasonable values for:
                        !     tscale_in   = 0.25
                        !     hscale_c_in = 0., 70., 0.  , 0.
                        !     hscale_v_in = 1.,  0., 0.83, 0.83
    odp (MAX(ntotsg,1)) ! approx. model layer thickness at the obs level
                        !   (to compute the minimum vertical correlation scale)

! Local parameters: None
! ----------------

  REAL    (KIND=wp   ),     PARAMETER           ::       &
    c100  =   100.0_wp  ! 100.0

! Local scalars:
! -------------

  INTEGER         ::  &
    nairll           ,& ! length of 'long' list
    ntotnsl ,nglobnsl,& ! number of new single-level (sing-l.) aircraft reports
    ntotasl ,nglobasl,& ! number of all sing-l. reports with station id in list
    ntotlon          ,& ! number of all sing-l. rep. to be processed:'long' list
    ntotlsl ,ntotlsr ,& ! number of local sing-l. reports to be processed
    ntottsl          ,& ! number of locally hosted reports to be processed
    ntotssl          ,& ! number of sing-l. reports with current station id
    ntotaml          ,& ! number of locally produced (hosted) multi-lev. reports
    ntotsml          ,& ! number of multi-level reports with current station id
    nstaid           ,& ! number of station id's in list of new aircraft reports
    ntarget          ,& ! number of all multi-level reports in the target node
    nglobaml,kglobaml,& ! global number of produced (hosted) multi-lev. reports
    nall             ,& ! total global number of elements to be gathered
    ncandls          ,& ! starting / ending index on candidates list
    ncanda  ,ncande  ,& ! index denoting the candidates list
!   nexcloc          ,& ! number of reports in excess of local array size
    irtot            ,& ! total global number of reports in 'long' list
    ndiffid          ,& ! number of rejected station id's
    maxaidp          ,& ! max. allowed number of single-level aircraft reports
                        !   with the same station ID - if this is exceeded,
                        !   the station ID is assumed to be incorrect
    ntotacd (5)      ,& ! number of reports of a certain aircraft code type
    izerror             ! error status variable

  INTEGER         ::  &
    nonl_cart_id     ,& ! 'my_cart_id' of node (sub-domain) for which lwonl=.t.
    irep, irps, inls ,& ! loop indices over reports
    ista             ,& ! loop index over stations
    inode            ,& ! loop index over nodes
    mvar             ,& ! loop index over variables
    ilsq             ,& ! loop index over list elements
    ilev             ,& ! loop index over vertical observation levels
    nsga    ,nmla    ,& ! sing-l./ multi-l. report index in ODR or temporary ODR
    kobtyp , kcdtyp  ,& ! observation type , observation code type
    ndiff            ,& ! number of removed multi-level aircraft reports
    ndifo            ,& ! number of removed old multi-level aircraft reports
    narest           ,& ! ODR station index used to reset elements to missing
    nrep             ,& ! number of reports in some lists / arrays
    istoffs          ,& ! station number offset, so that most reports at 'lwonl'
    ntsl    ,mxtsl   ,& ! quantities used to determine 'istoffs'
    nrepnm1          ,& ! number of reports up to inode-1
    nlev             ,& ! number of levels in multi-level report
    nlon             ,& ! report index in 'long' list
!   irnode           ,& ! node with subdomain containing the observ. location
    minlev           ,& ! smallest number of vertical levels in multi-l. reports
    irepass          ,& ! index of multi-level report to be cancelled
    npassiv          ,& ! indicator to set single-level report to passive
    icma             ,& ! indices for statistics ('cma', events)
    ilen   , icl     ,& ! length of output record
    iante            ,& ! number of elements already used
    ndsizel          ,& ! dimension (size) of local buffers
    ndsizet          ,& ! dimension (size) of hosted buffers, temporary arrays
    npej             ,& ! number of PEs j receiving data from PE i
    nft              ,& ! number of printed lines used for PE i
    iftn             ,& ! number of PEs j on current line
    ifta    ,ifte    ,& ! interval of PEs j on current line
    ift              ,& ! loop index over PEs j
    inpe             ,& ! loop index over nodes
    izmplcode        ,& ! error variable for MPI error code / opening files
    nzerr  ,ierr,jerr   ! indicator of status for (de-)allocation of fields

  INTEGER         ::  &
    mxobs            ,& ! max. number of observ. available at local process
    mxibls           ,& ! no. of integer elements to be gathered per obs & lev
    mxrbls           ,& ! number of real elements to be gathered per obs & lev
    mxybls           ,& ! number of charac. elements to be gathered per obs.
    ichlen           ,& ! length of characters elements to be gathered
    nelemi           ,& ! number of integer elements \  in each report
    nelemr              ! number of real    elements /  (for reports_all2all)

  REAL    (KIND=wp   )     ::  &
    tscale       =    0.25_wp                                 ,& ! see: 'tscale_in'
    vscale   (4) = (/ 0.577_wp, 0.577_wp, 0.2_wp , 0.2_wp  /) ,& ! see: 'vscale_in'
    hscale_c (4) = (/ 0.0_wp  ,70.0_wp  , 0.0_wp , 0.0_wp  /) ,& ! see: 'hscale_c_in'
    hscale_v (4) = (/ 1.0_wp  , 0.0_wp  , 0.83_wp, 0.83_wp /) ,& ! see: 'hscale_v_in'
    zodp(MAX(ntotsg,1)) ! approx. model layer thickness at the obs level

  LOGICAL                  ::  &
    lfound           ,& ! current station id found in list of station id's
    lremove             ! current multi-level report to be removed

  CHARACTER (LEN=75) yerrmsg,yerr ! error message
  CHARACTER (LEN=25) yroutine     ! name of this subroutine
  CHARACTER (LEN= 9) yformat      ! format
  CHARACTER (LEN=23) yformatl     ! long format
  CHARACTER (LEN= 3)       , PARAMETER  :: &
    yfh (12) = (/ 'PEi', ' | ', 'snd', 'rcv', 'src', 'srd'  &   ! specification
                , 'sic', 'sid', 'rrc', 'rrd', 'ric', 'rid' /)   ! of the output
  CHARACTER (LEN= 1)       , PARAMETER  :: &
    yfhs     =    ':'             ! for output

! Local arrays:
! ------------

  INTEGER              ::  &
    iscount  (num_compute),& ! number of rep. sent by local node to other nodes
    ircount  (num_compute),& ! number rep. received by local node from other n.
    isrcount (num_compute),& ! number of real    elements sent     by local node
    isicount (num_compute),& ! number of integer elements sent     by local node
    irrcount (num_compute),& ! number of real    elements received by local node
    iricount (num_compute),& ! number of integer elements received by local node
    isrdispl (num_compute),& ! displacement of real elements sent by local node
    isidispl (num_compute),& ! displacement of int. elements sent by local node
    irrdispl (num_compute),& ! displacement of real elements received by loc. n.
    iridispl (num_compute),& ! displacement of int. elements received by loc. n.
    ntotoml  (num_compute),& ! number of old multi-level reports in each node
    ntargsl  (num_compute),& ! number of all single-level reports in each
                             !                                   (target) node
    ntargaml (num_compute),& ! number of all multi-level  reports in each
                             !                                   (target) node
    naltoal  (num_compute,2) ! number of reports sent to PE j, number of PE j

  INTEGER, ALLOCATABLE  ::       &
    ibufls     (:,:) ,& ! for intermediate storage for gathering
    ibufls2    (:,:) ,& ! for intermediate storage for gathering
    isendbuf   (:,:) ,& ! sending   buffer for all-to-all communication
    irecvbuf   (:,:) ,& ! receiving buffer for all-to-all communication
    imlls      (:,:) ,& ! list of produced multi-level reports
    alltoall     (:) ,& ! number of reports sent by PE 'i' to PE 'j'
    irprcs       (:) ,& ! indices of reports to be processed
    irnode       (:) ,& ! nodes to which the reports will be distributed
    irsort       (:)    ! report indices sorted according to 'irnode'

  REAL    (KIND=wp)       , ALLOCATABLE  ::       &
    rbufls     (:,:) ,& ! for intermediate storage
    rbufls2    (:,:) ,& ! for intermediate storage
    rsendbuf   (:,:) ,& ! sending   buffer for all-to-all communication
    rrecvbuf   (:,:)    ! receiving buffer for all-to-all communication

  CHARACTER (LEN=ilstid ) , ALLOCATABLE  ::       &
    ybufls     (:,:) ,& ! for intermediate storage
    ybufls2    (:,:)    ! for intermediate storage
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_air_org_mult
!-------------------------------------------------------------------------------

  izerror = 0
  yroutine = 'obs_air_org_mult'

  nexceair(1) = 0
  nexceair(2) = 0

  IF ((PRESENT(odp)) .AND. (lsvcorl)) THEN
    DO irep = 1 , ntotsg
      zodp (irep) = odp(irep)
    ENDDO
  ELSE
    DO irep = 1 , ntotsg
      zodp (irep) = epsy
    ENDDO
  ENDIF

  IF ((PRESENT(tscale_in)) .AND. (lsvcorl))  tscale = tscale_in
  IF ((PRESENT(vscale_in)) .AND. (lsvcorl)) THEN
    DO mvar = 1 , 4
      vscale   (mvar) = vscale_in  (mvar)
    ENDDO
  ENDIF
  IF ((PRESENT(hscale_c_in)) .AND. (PRESENT(hscale_v_in)) .AND. (lsvcorl)) THEN
    DO mvar = 1 , 4
      hscale_c (mvar) = hscale_c_in(mvar)
      hscale_v (mvar) = hscale_v_in(mvar)
    ENDDO
  ENDIF

  nonl_cart_id = 0
  IF (lwonl)  nonl_cart_id = my_cart_id
  IF (num_compute > 1) THEN
    CALL global_values ( nonl_cart_id, 1, 'MAX', imp_integers, icomm_cart, -1  &
                       , yerr, ierr )
!   ==================
  ENDIF

  ntotacd = 0
  DO irep = 1 , ntotsg
    IF (mosghd(irep,nhcode) == namdar)  ntotacd (1) = ntotacd(1) + 1
    IF (mosghd(irep,nhcode) == nacar )  ntotacd (2) = ntotacd(2) + 1
    IF (mosghd(irep,nhcode) == naircd)  ntotacd (3) = ntotacd(3) + 1
    IF (mosghd(irep,nhcode) == ncodar)  ntotacd (4) = ntotacd(4) + 1
    IF (mosghd(irep,nhcode) == nmodes)  ntotacd (5) = ntotacd(5) + 1
  ENDDO
  IF (num_compute > 1) THEN
    CALL global_values ( ntotacd, 5, 'SUM', imp_integers, icomm_cart, -1       &
                       , yerr, ierr )
!   ==================
  ENDIF
  maxaidp  =  MAX( maxaid , MAXVAL( ntotacd ) / 2 )

!-------------------------------------------------------------------------------
!  Section 1: Preliminaries, including
!             - the making of a list of station id's from all aircraft reports
!               that have been read from the AOF at the current timestep
!             - the making of a complete list of single-level reports with these
!               station id's
!             - the removal of all multi-level reports with these station id's
!-------------------------------------------------------------------------------

  nairls = MAX( maxsgl , maxaidp , 1 )

  ALLOCATE ( iairls   (nairls , 2 , 3)  , STAT=nzerr )
  ALLOCATE ( yairls   (nairls     , 3)  , STAT=nzerr )

!-------------------------------------------------------------------------------
!  Section 1.1: Making of a ('long') list with all active single-level aircraft
!               reports (from the total ('global') model domain) that have been
!               read at the current timestep (from AOF file or NetCDF files)
!-------------------------------------------------------------------------------

  ALLOCATE ( rairls   (nairls        )  , STAT=nzerr )

  ncanda  = ntotsgo + 1
  ntotnsl = 0

  CALL obs_air_list ( 9, ncanda, ntotsg , 0, 0, 0 , 1 , ntotnsl )
! =================

! model layer thickness at report (approximated; for minimum vertical
! -------------------------------                correlation scale)

  DO irep = 1 , ntotnsl
    nsga = iairls (irep,1,1)
    rairls (irep)  =  zodp (nsga)
!   io   = mosghd (nsga,nhio)
!   jo   = mosghd (nsga,nhjo)
!   zpob = osgbdy (nsga,nbsp)
!   rairls (irep)  =  p2dp ( zpob , io , jo )
!                     ====
  ENDDO

! gather the local complements from all nodes and add them to the list
! --------------------------------------------------------------------

  IF (num_compute > 1) THEN

    nglobnsl = ntotnsl

    CALL global_values ( nglobnsl, 1,'SUM', imp_integers, icomm_cart, -1       &
                       , yerr, ierr )
!   ==================

!   mxobs   = nairll
    mxobs   = MAX( nglobnsl , 1 )
    mxibls  = 4
    mxrbls  = 1
    mxybls  = 1
    ichlen  = ilstid

    IF (nzerr == 0) ALLOCATE (ibufls (mxibls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufls (mxrbls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufls (mxybls, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'("* Allocation of xbufls failed *",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 11704, yerrmsg, yroutine)
    ENDIF

! write 'local' information to buffers
    DO   irep = 1 , ntotnsl
      ibufls (1,irep)  =  iairls (irep,1,1)
      ibufls (2,irep)  =  iairls (irep,2,1)
      ibufls (3,irep)  =  my_cart_id
      ibufls (4,irep)  =  - 1
      ybufls (1,irep)  =  yairls (irep  ,1)
      rbufls (1,irep)  =  rairls (irep    )
    ENDDO

! gather information globally

    CALL gather_all ( mxibls, mxrbls, mxybls, mxobs  , ichlen                  &
                    , ibufls, rbufls, ybufls, ntotnsl, nall , yerrmsg, nzerr)
!   ===============

!   (this should never happen)
    IF ((my_cart_id == 0) .AND. (nall > mxobs))                                &
      PRINT       '("CAUTION !!!!! t=",F6.3,":",I5," NEW SINGLE-LEVEL AIR"     &
                  &,"CRAFT REPS., ARRAY SIZE ",I5)' , acthr, nall, ntotnsl
!   IF ((lwonl ) .AND. (nall > mxobs))                                         &
!     WRITE( nupr,'("CAUTION !!!!! t=",F6.3,":",I5," NEW SINGLE-LEVEL AIR"     &
!                 &,"CRAFT REPS., ARRAY SIZE ",I5)' ) acthr, nall, ntotnsl

! write the global information from the buffers into the 'long' list
    nairll = MAX( ntotnsl , 1 )
    ALLOCATE ( iarlls   (nairll , 4)      , STAT=nzerr )
    ALLOCATE ( yarlls   (nairll    )      , STAT=nzerr )
    ALLOCATE ( rarlls   (nairll    )      , STAT=nzerr )

    DO   irep = 1 , ntotnsl
      iarlls (irep,1)  =  ibufls (1,irep)
      iarlls (irep,2)  =  ibufls (2,irep)
      iarlls (irep,3)  =  ibufls (3,irep)
      iarlls (irep,4)  =  ibufls (4,irep)
      yarlls (irep  )  =  ybufls (1,irep)
      rarlls (irep  )  =  rbufls (1,irep)
    ENDDO
    DEALLOCATE ( ibufls, STAT=nzerr )
    DEALLOCATE ( rbufls, STAT=nzerr )
    DEALLOCATE ( ybufls, STAT=nzerr )

  ELSE    !  num_compute = 1

! write to the 'long' list
!   (if (num_compute == 1) then the long list can never be longer than 'maxsgl')
    nairll = maxsgl
    ALLOCATE ( iarlls   (nairll , 4)      , STAT=nzerr )
    ALLOCATE ( yarlls   (nairll    )      , STAT=nzerr )
    ALLOCATE ( rarlls   (nairll    )      , STAT=nzerr )

    DO irep = 1 , ntotnsl
      iarlls (irep,1)  =  iairls (irep,1,1)
      iarlls (irep,2)  =  iairls (irep,2,1)
      iarlls (irep,3)  =  my_cart_id
      iarlls (irep,4)  =  - 1
      yarlls (irep  )  =  yairls (irep  ,1)
      rarlls (irep  )  =  rairls (irep    )
    ENDDO

  ENDIF

!-------------------------------------------------------------------------------
!  Section 1.2: Numbering and listing of the station id's
!               (which exist in the long list)
!-------------------------------------------------------------------------------

  nstaid = 0
  DO   irep = 1 , ntotnsl
    lfound = .FALSE.
    innerlist: DO inls = 1 , irep-1
      IF (yarlls(inls) (1:ilstidp) == yarlls(irep) (1:ilstidp)) THEN
        lfound = .TRUE.
        EXIT innerlist
      ENDIF
    ENDDO innerlist
    IF (.NOT. lfound) THEN
      nstaid = nstaid + 1
! here, as an exception, 'iairls(:,1,2)' contains the index 'irep' of list
! 'iairls(irep,:,1)' instead of an ODR index
      iairls (nstaid,1,2) = irep
      iairls (nstaid,2,2) = iarlls(irep,2)
      yairls (nstaid  ,2) = yarlls(irep  )
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
!  Section 1.3: Complementation of the ('long') list with old single-level
!               reports with the same station id. and passive flag less than 2
!               (or passive flag = 2 if the reason is not thinning, redundancy,
!                or invalid height)
!-------------------------------------------------------------------------------

! list with old local reports (local complement)
! ----------------------------------------------

  ntotasl = 0

  IF (nstaid > 0)                                                              &

    CALL obs_air_list ( 9, 1, ntotsgo , 2, 1, nstaid, 1 , ntotasl )
!   =================

! model layer thickness at report (approximated; for minimum vertical
! -------------------------------                correlation scale)

  DO irep = 1 , ntotasl
    nsga = iairls (irep,1,1)
    rairls (irep)  =  zodp (nsga)
!   io   = mosghd (nsga,nhio)
!   jo   = mosghd (nsga,nhjo)
!   zpob = osgbdy (nsga,nbsp)
!   rairls (irep)  =  p2dp ( zpob , io , jo )
!                     ====
  ENDDO

! gather the local complements from all nodes and add them to the list
! --------------------------------------------------------------------

  IF ( (num_compute > 1)  .AND. (nstaid > 0)) THEN

!   at first, copy hitherto long list into buffer arrays
!   (to allow for adaptive memory re-allocation of long list)

    mxobs   = MAX( ntotnsl , 1 )
    mxibls  = 4
    mxrbls  = 1
    mxybls  = 1
    IF (nzerr == 0) ALLOCATE (ibufls2 (mxibls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufls2 (mxrbls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufls2 (mxybls, mxobs), STAT=nzerr)

    DO   irep = 1 , ntotnsl
      ibufls2 (1,irep)  =  iarlls (irep,1)
      ibufls2 (2,irep)  =  iarlls (irep,2)
      ibufls2 (3,irep)  =  iarlls (irep,3)
      ibufls2 (4,irep)  =  iarlls (irep,4)
      ybufls2 (1,irep)  =  yarlls (irep  )
      rbufls2 (1,irep)  =  rarlls (irep  )
    ENDDO

    DEALLOCATE ( iarlls  , STAT=nzerr )
    DEALLOCATE ( yarlls  , STAT=nzerr )
    DEALLOCATE ( rarlls  , STAT=nzerr )

!  now gather the supplementary info and add it to the list

    nglobasl = ntotasl

    CALL global_values ( nglobasl, 1,'SUM', imp_integers, icomm_cart, -1       &
                       , yerr, ierr )
!   ==================

!   mxobs   = nairll
    mxobs   = MAX( nglobasl , 1 )
    mxibls  = 4
    mxrbls  = 1
    mxybls  = 1
    ichlen  = ilstid

    IF (nzerr == 0) ALLOCATE (ibufls (mxibls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufls (mxrbls, mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufls (mxybls, mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'("* Allocation of xbufls failed *",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 11704, yerrmsg, yroutine)
    ENDIF

! write 'local' information to buffers
    DO   irep = 1 , ntotasl
      ibufls (1,irep)  =  iairls (irep,1,1)
      ibufls (2,irep)  =  iairls (irep,2,1)
      ibufls (3,irep)  =  my_cart_id
      ibufls (4,irep)  =  - 1
      ybufls (1,irep)  =  yairls (irep  ,1)
      rbufls (1,irep)  =  rairls (irep    )
    ENDDO

! gather information globally

    CALL gather_all ( mxibls, mxrbls, mxybls, mxobs  , ichlen                  &
                    , ibufls, rbufls, ybufls, ntotasl, nall , yerrmsg, nzerr)
!   ===============

!   (this should never happen)
    IF ((my_cart_id == 0) .AND. (nall > mxobs))                                &
      PRINT       '("CAUTION !!!!! t=",F6.3,":",I5," AIRCRAFT REPS. TO BE"     &
                  &," PROCESSED , ARRAY SIZE ",I5)' , acthr, nall, ntotasl
!   IF ((lwonl ) .AND. (nall > mxobs))                                         &
!     WRITE( nupr,'("CAUTION !!!!! t=",F6.3,":",I5," AIRCRAFT REPS. TO BE"     &
!                 &," PROCESSED , ARRAY SIZE ",I5)' ) acthr, nall, ntotasl

! write the global information from the buffers into the 'long' list
    nairll = MAX( ntotnsl  + ntotasl , 1 )
!   (by construction, these arrays are large enough
!    to accommodate the complete long list)
    ALLOCATE ( iarlls   (nairll , 4)      , STAT=nzerr )
    ALLOCATE ( yarlls   (nairll    )      , STAT=nzerr )
    ALLOCATE ( rarlls   (nairll    )      , STAT=nzerr )

    DO   irep = 1 , ntotnsl
      iarlls (irep,1)  =  ibufls2 (1,irep)
      iarlls (irep,2)  =  ibufls2 (2,irep)
      iarlls (irep,3)  =  ibufls2 (3,irep)
      iarlls (irep,4)  =  ibufls2 (4,irep)
      yarlls (irep  )  =  ybufls2 (1,irep)
      rarlls (irep  )  =  rbufls2 (1,irep)
    ENDDO

    DO   irep = 1 , ntotasl
      irtot  =  ntotnsl + irep
      iarlls (irtot,1)  =  ibufls (1,irep)
      iarlls (irtot,2)  =  ibufls (2,irep)
      iarlls (irtot,3)  =  ibufls (3,irep)
      iarlls (irtot,4)  =  ibufls (4,irep)
      yarlls (irtot  )  =  ybufls (1,irep)
      rarlls (irtot  )  =  rbufls (1,irep)
    ENDDO

! read global information from buffers into the 'long' list
    DEALLOCATE ( ibufls , STAT=nzerr )
    DEALLOCATE ( rbufls , STAT=nzerr )
    DEALLOCATE ( ybufls , STAT=nzerr )
    DEALLOCATE ( ibufls2, STAT=nzerr )
    DEALLOCATE ( rbufls2, STAT=nzerr )
    DEALLOCATE ( ybufls2, STAT=nzerr )

  ELSE    !  num_compute = 1

! complement the 'long' list
    DO   irep = 1 , ntotasl
      irtot  =  ntotnsl + irep
      iarlls (irtot,1)  =  iairls (irep,1,1)
      iarlls (irtot,2)  =  iairls (irep,2,1)
      iarlls (irtot,3)  =  my_cart_id
      iarlls (irtot,4)  =  - 1
      yarlls (irtot  )  =  yairls (irep  ,1)
      rarlls (irtot  )  =  rairls (irep    )
    ENDDO

  ENDIF

  ntotlon = ntotasl + ntotnsl

  DEALLOCATE ( rairls  , STAT=nzerr )

!-------------------------------------------------------------------------------
!  Section 1.4: Exclusion of suspicious (default) aircraft id
!               (i.e. id which appears too often since it is probably used as
!                default id by the data base if something is wrong with the id)
!-------------------------------------------------------------------------------

  ndiffid = 0
  DO ista = 1 , nstaid
    nrep = 0
    DO irep = 1 , ntotlon
      IF (yarlls(irep) == yairls(ista,2))  nrep = nrep + 1
    ENDDO
    IF (nrep >= maxaidp) THEN
!   IF (nrep >= nairls) THEN
! cancel reports on long list which have the current bad station id,
! and set such reports to passive in local ODR
! (Note: local list 1 does not have to be adjusted since it is not used any more
!        before being re-made in section 2)
      ndiff   = 0
      DO   irep = 1 , ntotlon
        IF (yarlls(irep) == yairls(ista,2)) THEN
          ndiff = ndiff + 1
          IF (iarlls(irep,3) == my_cart_id) THEN
            ilen = 4 + istrej
            IF (nacout+ilen <= nmxoln) THEN
              outbuf(nacout+1) = ilen
              outbuf(nacout+2) = nfmt26
              outbuf(nacout+3) = NINT( osghed(iarlls(irep,1),nhtime) * c100 )
              outbuf(nacout+4) = NINT( osgbdy(iarlls(irep,1),nbsp  )        )
              DO icl = 1 , istrej
                outbuf(nacout+6+icl) = ICHAR( yarlls(irep) (icl:icl))
              ENDDO
              nacout  = nacout + ilen
            ENDIF
!           IF (lwonl)                                                         &
!             WRITE( nupr,'(A20,": CAUTION: suspicious aircraft station ID:"   &
!                         &,A ," ==> problem in data base / makeaof ???")' )   &
!                    yroutine, yairls(ista,2)
            PRINT         '(A20,": CAUTION: suspicious aircraft station ID:"   &
                          &,A ," ==> problem in data base / makeaof ???")' ,   &
                   yroutine, yairls(ista,2)
            IF (mosghd(iarlls(irep,1),nhpass) < 2) THEN
              kcdtyp = mosghd(iarlls(irep,1),nhcode)
              kobtyp = mosghd(iarlls(irep,1),nhobtp)
              icma  = i_cma ( kobtyp , kcdtyp )
              cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
              cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
              neventr (neblak,icma) = neventr(neblak,icma) + 1
              irps = iarlls(irep,1)
              mosghd (irps,nhpass) = 2
              mosghd (irps,nhschr) = IBSET ( mosghd(irps,nhschr), nvbkbp )
              mosghd (irps,nhschr) = IBSET ( mosghd(irps,nhschr), nvpsbp )
              mosghd (irps,nhschr) = IBSET ( mosghd(irps,nhschr), FL_BLACKLIST )
            ENDIF
          ENDIF
        ELSE
          iarlls (irep-ndiff,1)  =  iarlls(irep,1)
          iarlls (irep-ndiff,2)  =  iarlls(irep,2)
          iarlls (irep-ndiff,3)  =  iarlls(irep,3)
          iarlls (irep-ndiff,4)  =  iarlls(irep,4)
          yarlls (irep-ndiff  )  =  yarlls(irep  )
          rarlls (irep-ndiff  )  =  rarlls(irep  )
        ENDIF
      ENDDO
      ntotlon = ntotlon - ndiff
! cancel bad station id on station id list
      ndiffid = ndiffid + 1
    ELSE
      iairls (ista-ndiffid,1,2) = iairls(ista,1,2)
      iairls (ista-ndiffid,2,2) = iairls(ista,1,2)
      yairls (ista-ndiffid  ,2) = yairls(ista  ,2)
    ENDIF
  ENDDO
  nstaid = nstaid - ndiffid

!-------------------------------------------------------------------------------
!  Section 1.5: Removal of all multi-level reports with station id in the list
!-------------------------------------------------------------------------------

  ndiff  = 0
  ndifo  = 0
  DO nmla  = 1 , ntotml
    lremove = .FALSE.
    IF (      (momlhd(nmla,nhobtp) == nairep)                                  &
        .AND. (momlhd(nmla,nhcode) /= ncodar)) THEN
      listctl: DO ista = 1 , nstaid
        IF (yomlhd(nmla) (1:ilstidp) == yairls(ista,2) (1:ilstidp)) THEN
          lremove = .TRUE.
          ndiff = ndiff + 1
          IF (nmla <= ntotmlo)  ndifo = ndifo + 1
          kcdtyp = momlhd(nmla,nhcode)
          kobtyp = momlhd(nmla,nhobtp)
          icma  = i_cma ( kobtyp , kcdtyp )
          neventr (neslml,icma) = neventr (neslml,icma) - 1
          EXIT listctl
        ENDIF
      ENDDO listctl
    ENDIF
    IF ((.NOT. lremove) .AND. (ndiff > 0)) THEN
      omlbdy (nmla-ndiff,1:maxmlv,1:mxrbdy) = omlbdy (nmla,1:maxmlv,1:mxrbdy)
      momlbd (nmla-ndiff,1:maxmlv,1:mxrbdf) = momlbd (nmla,1:maxmlv,1:mxrbdf)
      omlhed (nmla-ndiff,1:mxrhed)          = omlhed (nmla,1:mxrhed)
      momlhd (nmla-ndiff,1:mxrhdf)          = momlhd (nmla,1:mxrhdf)
      yomlhd (nmla-ndiff)                   = yomlhd (nmla)
    ENDIF
  ENDDO

  IF (ndiff > 0) THEN
    narest = ntotml - ndiff + 1
    omlbdy (narest:ntotml,1:maxmlv,1:mxrbdy) = rmdi
    momlbd (narest:ntotml,1:maxmlv,1:mxrbdf) = imdi
    omlhed (narest:ntotml,1:mxrhed)          = rmdi
    momlhd (narest:ntotml,1:mxrhdf)          = imdi
    yomlhd (narest:ntotml)                   = '     '
  ENDIF

  ntotml  =  ntotml  - ndiff
  ntotmlo =  ntotmlo - ndifo

!-------------------------------------------------------------------------------
!  Section 1.6: No further processing in the absence of new aircraft reports
!-------------------------------------------------------------------------------

  IF (nstaid == 0) THEN
    DEALLOCATE ( iairls  , STAT=nzerr )
    DEALLOCATE ( yairls  , STAT=nzerr )
    DEALLOCATE ( iarlls  , STAT=nzerr )
    DEALLOCATE ( rarlls  , STAT=nzerr )
    DEALLOCATE ( yarlls  , STAT=nzerr )
  ENDIF


IF (nstaid > 0) THEN


!-------------------------------------------------------------------------------
!  Section 2  : Scattering and gathering of the required single-level reports
!-------------------------------------------------------------------------------

  IF (num_compute > 1) THEN

! determine the target nodes of the single-level reports in the 'long' list
! for the production of multi-level reports
! -------------------------------------------------------------------------

! get maximum number of reports on any target host
    ntargsl (:) = 0
    DO ista = 1 , nstaid
      DO irep = 1 , ntotlon
        IF (yarlls(irep) == yairls(ista,2))                                    &
          ntargsl (MOD(ista,num_compute)+1) =ntargsl(MOD(ista,num_compute)+1) +1
      ENDDO
    ENDDO
    mxtsl = -1
    DO inode = 1 , num_compute
      IF (ntargsl(inode) > mxtsl) THEN
        ntsl  = inode
        mxtsl = ntargsl (inode)
      ENDIF
    ENDDO
! get offset, so that if (lwonl) then the number of hosted reports is maximum
    istoffs  =  num_compute +  nonl_cart_id  -  (ntsl - 1)

    DO ista = 1 , nstaid
      DO irep = 1 , ntotlon
        IF (yarlls(irep) == yairls(ista,2))                                    &
          iarlls (irep,4) = MOD( ista + istoffs , num_compute )
      ENDDO
!     PRINT '("LOCAL:",3I4,1X,A ,I7,I4,2I6)'                                   &
!           , my_cart_id, nstaid, ista, yairls(ista,2)                         &
!           , ntargsl(MOD(ista,num_compute)+1), istoffs, nairll, nairls
    ENDDO

! re-make a list with the local single-level reports to be processed
! ------------------------------------------------------------------

    ntotlsl = 0
    DO irep = 1 , ntotlon
      IF (iarlls(irep,3) == my_cart_id) THEN
        ntotlsl = ntotlsl + 1
        iairls (ntotlsl,1,1) = irep
        iairls (ntotlsl,2,1) = iarlls (irep,2)
        yairls (ntotlsl,  1) = yarlls (irep  )
      ENDIF
    ENDDO

! Determine dimensions of buffers and temporary arrays
! ----------------------------------------------------

    ndsizel  =  ntotlsl

    CALL global_values ( ndsizel, 1, 'MAX', imp_integers, icomm_cart, -1       &
                       , yerrmsg, nzerr )
!   ==================

    IF (lwonl) WRITE( nupr,'(''ndsizel lsl '',I4)' ) ndsizel

! Store the local single-level reports (including the 'long' list index) in a
! sending buffer, sorted according to the target nodes of the reports.
! Determine the number of local reports to be sent to the different target nodes
! ------------------------------------------------------------------------------

    nelemr = mxshed + mxsbdy
    nelemi = mxshdf + mxsbdf + 1
    ALLOCATE ( rsendbuf (nelemr , ndsizel)            , STAT=nzerr )
    ALLOCATE ( isendbuf (nelemi , ndsizel)            , STAT=nzerr )
    ALLOCATE ( alltoall (num_compute * num_compute)   , STAT=nzerr )

    rsendbuf (:,:) = c0
    isendbuf (:,:) = 0
    alltoall   (:) = 0
    nrep           = 0

    DO inode  = 0 , num_compute-1
      nrepnm1 = nrep
      DO irep = 1 , ntotlsl
        nlon = iairls (irep,1,1)
        IF (inode == iarlls(nlon,4)) THEN
          nrep = nrep + 1
          nsga = iarlls (nlon,1)
          DO mvar = 1 , mxshed
            rsendbuf (mvar,nrep) = osghed (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxshdf
            isendbuf (mvar,nrep) = mosghd (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxsbdy
            rsendbuf (mxshed+mvar,nrep) = osgbdy (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxsbdf
            isendbuf (mxshdf+mvar,nrep) = mosgbd (nsga,mvar)
          ENDDO
          isendbuf (mxshdf+mxsbdf+1,nrep) = nlon
        ENDIF
      ENDDO
      alltoall (my_cart_id*num_compute+inode+1)  =  nrep - nrepnm1
    ENDDO

! Get the numbers of local reports to be sent to the different target nodes,
! and the numbers of reports to be received from the different source nodes.
! --------------------------------------------------------------------------

    CALL global_values ( alltoall, num_compute*num_compute, 'MAX'              &
                       , imp_integers, icomm_cart, -1, yerrmsg, nzerr )
!   ==================

    ndsizet = 0
    DO inode = 1, num_compute
      iscount (inode) =  alltoall (my_cart_id*num_compute +inode       )
      ircount (inode) =  alltoall ((inode-1) *num_compute +my_cart_id+1)
      ndsizet  =  ndsizet + ircount(inode)
    ENDDO

! printing for control
    IF (lwonl) THEN
      WRITE( nupr, '(''ALL_TO_ALLV of single-level aircraft reports: ''        &
                   &,''# reports sent by PE i to PE j'')' )
      WRITE( yformatl,'(''(3X,'',I3,''I3)           '')' ) num_compute
      DO inode = 0 , num_compute-1
        IF (num_compute <= 25) THEN
          WRITE(nupr,yformatl) (alltoall(ift),ift=inode*num_compute+1          &
                                                 ,(inode+1)*num_compute)
        ELSE
          npej = 0
          DO inpe = 1 , num_compute
            IF (alltoall(inode*num_compute+inpe) > 0) THEN
              npej = npej + 1
              naltoal (npej,1) = alltoall(inode*num_compute+inpe)
              naltoal (npej,2) = inpe - 1
            ENDIF
          ENDDO
          nft = MIN( INT( (npej-1) / 8 + 1 ) , npej )
          DO inpe = 1 , nft
            ifta = (inpe-1)*8 + 1
            ifte = MIN( INT( inpe*8 ) , npej )
            iftn = ifte - ifta + 1
            WRITE( yformatl,'(''(A ,I4,'',I2,''(A ,I3,A1,I3))'')' ) iftn
            WRITE( nupr, yformatl ) yfh(1), inode                              &
                 , (yfh(2),naltoal(ift,2),yfhs,naltoal(ift,1), ift=ifta,ifte)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    DEALLOCATE (alltoall  , STAT=nzerr)

! Determine dimensions of receiving buffers and temporary arrays
! --------------------------------------------------------------

    CALL global_values ( ndsizet, 1, 'MAX', imp_integers, icomm_cart, -1       &
                       , yerrmsg, nzerr)
!   ==================

    IF (lwonl) WRITE( nupr,'(''ndsizet tsl '',I4)' ) ndsizet

! Print for control output
! ------------------------

    IF ((lwonl) .AND. (lfirst)) THEN
      DO inode = 1, num_compute
        IF (inode >= 2) THEN
          isrdispl (inode)  =  isrdispl(inode-1) + isrcount(inode-1)
          isidispl (inode)  =  isidispl(inode-1) + isicount(inode-1)
          irrdispl (inode)  =  irrdispl(inode-1) + irrcount(inode-1)
          iridispl (inode)  =  iridispl(inode-1) + iricount(inode-1)
        ELSE
          isrdispl (inode)  =  0
          isidispl (inode)  =  0
          irrdispl (inode)  =  0
          iridispl (inode)  =  0
        ENDIF
        isrcount (inode)  =  iscount(inode) * nelemr
        isicount (inode)  =  iscount(inode) * nelemi
        irrcount (inode)  =  ircount(inode) * nelemr
        iricount (inode)  =  ircount(inode) * nelemi
      ENDDO
      nft = (num_compute-1) / 20 + 1
      DO inode = 1 , nft
        ifta = (inode-1)*20 + 1
        ifte = MIN( INT( inode*20 ) , num_compute )
        iftn = ifte - ifta + 1
        WRITE( yformat,'(''(A ,'',I2,''I4)'')' ) iftn
        WRITE( nupr, '(''PE '',I2,'': nodes'',I4,'' to'',I4)' )                &
               my_cart_id, ifta, ifte
        WRITE( nupr, yformat ) yfh( 3), (iscount (ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 4), (ircount (ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 5), (isrcount(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 6), (isrdispl(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 7), (isicount(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 8), (isidispl(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh( 9), (irrcount(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh(10), (irrdispl(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh(11), (iricount(ift),ift=ifta,ifte)
        WRITE( nupr, yformat ) yfh(12), (iridispl(ift),ift=ifta,ifte)
      ENDDO
    ENDIF

! Scatter and gather the single-level reports
! -------------------------------------------

    ALLOCATE ( rrecvbuf (nelemr , ndsizet)     , STAT=nzerr )
    ALLOCATE ( irecvbuf (nelemi , ndsizet)     , STAT=nzerr )

    CALL reports_all2all ( iscount , ircount , num_compute                     &
                         , ndsizel , isendbuf, rsendbuf , nelemi , nelemr      &
                         , ndsizet , irecvbuf, rrecvbuf , yerrmsg, izmplcode )
  ! ====================

    IF (izmplcode /= 0)                                                        &
      CALL model_abort ( my_cart_id, 11702, yerrmsg, yroutine, izmplcode )

    DEALLOCATE (rsendbuf , STAT=nzerr)
    DEALLOCATE (isendbuf , STAT=nzerr)

! Store the received single-level reports in temporary arrays
! and make a list of these locally hosted reports
! -----------------------------------------------------------

    ALLOCATE ( zsghed (ndsizet , mxshed) , STAT=nzerr )
    ALLOCATE ( mzsghd (ndsizet , mxshdf) , STAT=nzerr )
    ALLOCATE ( zsgbdy (ndsizet , mxsbdy) , STAT=nzerr )
    ALLOCATE ( mzsgbd (ndsizet , mxsbdf) , STAT=nzerr )
    ALLOCATE ( mzslon (ndsizet         ) , STAT=nzerr )

    ntottsl = 0
    DO inode = 1, num_compute
      ntottsl = ntottsl + ircount(inode)
    ENDDO
    DO irep = 1 , ntottsl
      DO mvar = 1 , mxshed
        zsghed (irep,mvar) = rrecvbuf (mvar,irep)
      ENDDO
      DO mvar = 1 , mxshdf
        mzsghd (irep,mvar) = irecvbuf (mvar,irep)
      ENDDO
      DO mvar = 1 , mxsbdy
        zsgbdy (irep,mvar) = rrecvbuf (mxshed+mvar,irep)
      ENDDO
      DO mvar = 1 , mxsbdf
        mzsgbd (irep,mvar) = irecvbuf (mxshdf+mvar,irep)
      ENDDO
      mzslon (irep)     = irecvbuf (mxshdf+mxsbdf+1,irep)
      iairls (irep,1,1) = irep
      iairls (irep,2,1) = iarlls (mzslon(irep),2)
      yairls (irep  ,1) = yarlls (mzslon(irep)  )
    ENDDO

    DEALLOCATE (rrecvbuf , STAT=nzerr)
    DEALLOCATE (irecvbuf , STAT=nzerr)

  ELSE    !  num_compute = 1

    ndsizet = maxaidp
    ALLOCATE ( zsghed (ndsizet , mxshed) , STAT=nzerr )
    ALLOCATE ( mzsghd (ndsizet , mxshdf) , STAT=nzerr )
    ALLOCATE ( zsgbdy (ndsizet , mxsbdy) , STAT=nzerr )
    ALLOCATE ( mzsgbd (ndsizet , mxsbdf) , STAT=nzerr )
    ALLOCATE ( mzslon (ndsizet         ) , STAT=nzerr )

    istoffs = 0

  ENDIF

!-------------------------------------------------------------------------------
!  Section 3  : Production of multi-level reports
!-------------------------------------------------------------------------------

  ALLOCATE ( zmlhed (maxmll          , mxrhed) , STAT=nzerr )
  ALLOCATE ( mzmlhd (maxmll          , mxrhdf) , STAT=nzerr )
  ALLOCATE ( zmlbdy (maxmll , maxarl , mxrbdy) , STAT=nzerr )
  ALLOCATE ( mzmlbd (maxmll , maxarl , mxrbdf) , STAT=nzerr )
  ALLOCATE ( mzmlon (maxmll , maxarl         ) , STAT=nzerr )

  ntotaml = 0

  loop_over_station_ids: DO ista = 1 , nstaid

! select station id's for each node
! ---------------------------------

    IF (MOD( ista + istoffs , num_compute ) == my_cart_id) THEN
      IF (num_compute > 1) THEN
        ncandls = 1
        ncande  = ntottsl
      ELSE    !  num_compute = 1
        ncandls = 8
        ncande  = ntotlon
      ENDIF
      ntotssl = 0

! make list of reports with present station id
! --------------------------------------------

      CALL obs_air_list ( ncandls, 1, ncande, 2, ista, ista, 3 , ntotssl )
!     =================

! non-MPP runs: store only the reports with present sta. id in temporary array
! ----------------------------------------------------------------------------

      IF (num_compute == 1) THEN
        DO irep = 1 , ntotssl
          nlon  = iairls (irep,1,3)
          nsga  = iarlls (nlon,1)
          DO mvar = 1 , mxshed
            zsghed (irep,mvar) = osghed (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxshdf
            mzsghd (irep,mvar) = mosghd (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxsbdy
            zsgbdy (irep,mvar) = osgbdy (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxsbdf
            mzsgbd (irep,mvar) = mosgbd (nsga,mvar)
          ENDDO
          mzslon (irep) = nlon
! adjust 'iairls(:,1,3)' so that it contains the same as for MPP runs:
! the report index in the temporary arrays 'zsghed', 'mzsghd', 'zsgbdy', 'mzsgbd
          iairls (irep,1,3) = irep
        ENDDO
      ENDIF

! produce multi-level reports for one station id
! ----------------------------------------------

!     IF (ntotssl > 0)                                                         &
!       PRINT '("COLLECTED:",I4,4X,I4,1X,A ,I7   )'                            &
!             , my_cart_id, ista, yarlls(mzslon(iairls(1,1,3))), ntotssl

      IF (ntotssl > 0)                                                         &

        CALL obs_air_make_mult ( ntotssl, thairh , ntotaml, nexceair(2) )
!       ======================

    ENDIF
  ENDDO loop_over_station_ids


!-------------------------------------------------------------------------------
!  Section 4  : If the number of multi-level reports exceeded the size of the
!               local and global arrays, then convert the aircraft multi-level
!               reports with the smallest number of vertical levels back to
!               single-level reports
!-------------------------------------------------------------------------------

! determine: 'ntotoml' : local numbers of multi-level reports not produced here
!            'nglobaml': global number of multi-level reports produced in
!                        'obs_air_make_mult', used for array allocation
! -----------------------------------------------------------------------------

  DO inode = 1, num_compute
    ntotoml (inode) = 0
  ENDDO
  ntotoml (my_cart_id+1) = ntotml

  IF (num_compute > 1) THEN

    nglobaml = ntotaml

    CALL global_values ( nglobaml, 1, 'SUM', imp_integers                      &
                       , icomm_cart, -1, yerrmsg, nzerr )
!   ==================

    CALL global_values ( ntotoml, num_compute, 'MAX', imp_integers             &
                       , icomm_cart, -1, yerrmsg, nzerr )
!   ==================

  ELSE    !  num_compute = 1

    nglobaml = ntotaml

  ENDIF

! make a list over the (locally) produced multi-level reports containing the
! source node, the target node, the number of levels, a marker 'active', and
! the local report index
! (source node: the node which produced the multi-level report
!  target node: the node with the domain which the horizontal location of the
!               produced multi-level reports lies in)
! -----------------------------------------------------------------------------

  ALLOCATE ( imlls (MAX( nglobaml,1 ) , 5) , STAT=nzerr )

  IF (ntotaml > 0) THEN
    ALLOCATE ( irprcs (ntotaml+1) , STAT=nzerr )
    ALLOCATE ( irnode (ntotaml+1) , STAT=nzerr )
    ALLOCATE ( irsort (ntotaml+1) , STAT=nzerr )
    DO irps = 1 , ntotaml
      irprcs (irps) = irps
    ENDDO

    CALL obs_assign_sort_node ( maxmll, ntotaml, irprcs, mzmlhd(:,nhitot)      &
                              , mzmlhd(:,nhjtot), num_compute, i_subpos        &
                              , nboundlines, my_cart_id , irnode, irsort )
!   =========================

    DO irps = 1 , ntotaml
      irep  =  irsort(irps)
      imlls (irep,1) = my_cart_id
      imlls (irep,2) = irnode(irps)
      imlls (irep,3) = mzmlhd(irep,nhnlev)
      imlls (irep,4) = 1
      imlls (irep,5) = irep
    ENDDO
    DEALLOCATE ( irprcs , STAT=nzerr )
    DEALLOCATE ( irnode , STAT=nzerr )
    DEALLOCATE ( irsort , STAT=nzerr )
  ENDIF

! DO irep = 1 , ntotaml
!   io_tot = mzmlhd(irep,nhitot)
!   jo_tot = mzmlhd(irep,nhjtot)
!   Nodes: DO  inode = 0, num_compute-1
!     IF (     i_subpos(inode,1) <= io_tot .AND. i_subpos(inode,3) >= io_tot     &
!        .AND. i_subpos(inode,2) <= jo_tot .AND. i_subpos(inode,4) >= jo_tot) THEN
!       irnode = inode
!       EXIT Nodes
!     ENDIF
!   ENDDO Nodes
!   imlls (irep,1) = my_cart_id
!   imlls (irep,2) = irnode
!   imlls (irep,3) = mzmlhd(irep,nhnlev)
!   imlls (irep,4) = 1
!   imlls (irep,5) = irep
! ENDDO

! gather this list, and make a list with the local numbers of the old and/or
!                                            non-aircraft multi-level reports
! ---------------------------------------------------------------------------

  IF (num_compute > 1) THEN

    mxobs   = MAX( nglobaml , 1 )
    mxibls  = 5
    mxrbls  = 0
    mxybls  = 0
    ichlen  = ilstid

    IF (nzerr == 0) ALLOCATE (ibufls (    mxibls   , mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (rbufls (MAX(mxrbls,1), mxobs), STAT=nzerr)
    IF (nzerr == 0) ALLOCATE (ybufls (MAX(mxybls,1), mxobs), STAT=nzerr)
    jerr = ABS( nzerr )
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsg,'("* Allocation of xbufls failed *",2I5)' ) jerr, ierr
      CALL model_abort (my_cart_id, 11705, yerrmsg, yroutine)
    ENDIF

! write 'local' information to buffers
    DO   irep = 1 , ntotaml
      DO ilsq = 1 , 5
        ibufls (ilsq,irep) = imlls (irep,ilsq)
      ENDDO
    ENDDO
    kglobaml = ntotaml

! gather information globally

    CALL gather_all ( mxibls, mxrbls, mxybls, mxobs   , ichlen                 &
                    , ibufls, rbufls, ybufls, kglobaml, nall , yerrmsg, nzerr)
!   ===============

    IF ((kglobaml /= nglobaml) .OR. (nall > mxobs)) THEN
      WRITE( yerrmsg,'("* Error in gathering lists *",4I5)' ) kglobaml         &
                                                            , nglobaml, nall, mxobs
      CALL model_abort (my_cart_id, 11706, yerrmsg, yroutine)
    ENDIF

! read global information from buffers
    DO   irep = 1 , nglobaml
      DO ilsq = 1 , 5
        imlls (irep,ilsq) = ibufls (ilsq,irep)
      ENDDO
    ENDDO

    DEALLOCATE (ibufls, STAT=nzerr)
    DEALLOCATE (rbufls, STAT=nzerr)
    DEALLOCATE (ybufls, STAT=nzerr)
  ENDIF

! determine the total local number of multi-level reports for each target node
! separately; and if it exceeds the array size 'maxmll', set the (locally) ex-
! cessive multi-level reports with the smallest number of vertical levels to
! 'passive' (this is done by each node for each (sub-)domain)
! ----------------------------------------------------------------------------

  ntargaml (:) = 0

  DO   irep = 1 , nglobaml
    ntargaml (imlls(irep,2)+1) = ntargaml(imlls(irep,2)+1) + 1
  ENDDO

  loop_over_nodes: DO inode = 0, num_compute-1
    ntarget  =  ntotoml(inode+1) + ntargaml(inode+1)
    IF (inode == my_cart_id)  nexceair (1) = MAX( ntarget - maxmll , 0 )
    DO WHILE (ntarget > maxmll)
      minlev = maxarl + 1
      DO irep = 1 , nglobaml
        IF (imlls(irep,2) == inode) THEN
          IF ((imlls(irep,4) == 1) .AND. (imlls(irep,3) < minlev)) THEN
            irepass = irep
            minlev  = imlls(irep,3)
          ENDIF
        ENDIF
      ENDDO
      imlls (irepass,4) = 0
      ntarget  = ntarget - 1
    ENDDO
  ENDDO loop_over_nodes

! if my node is the source node, cancel any excessive report
! (and set the sourcing single-level reports on the list back to active)
! ----------------------------------------------------------------------

  DO   irep = 1 , nglobaml
    IF ((imlls(irep,1) == my_cart_id) .AND. (imlls(irep,4) == 0))              &
      mzmlhd (imlls(irep,5),nhaexi) = 0
  ENDDO
  ndiff  = 0
  DO irep = 1 , ntotaml
    IF (mzmlhd(irep,nhaexi) == 0) THEN
      ndiff = ndiff + 1
      kcdtyp = mzmlhd (irep,nhcode)
      kobtyp = mzmlhd (irep,nhobtp)
      icma  = i_cma ( kobtyp , kcdtyp )
      neventr (nenoml,icma) = neventr (nenoml,icma) + 1
! set single-level reports back to active on 'long' list
      nlev = mzmlhd(irep,nhnlev)
      DO ilev = 1 , nlev
        iarlls (mzmlon(irep,ilev),2) = kcdtyp
      ENDDO
    ELSEIF (ndiff > 0) THEN
      zmlbdy (irep-ndiff,1:maxarl,1:mxrbdy) = zmlbdy (irep,1:maxarl,1:mxrbdy)
      mzmlbd (irep-ndiff,1:maxarl,1:mxrbdf) = mzmlbd (irep,1:maxarl,1:mxrbdf)
      zmlhed (irep-ndiff,1:mxrhed)          = zmlhed (irep,1:mxrhed)
      mzmlhd (irep-ndiff,1:mxrhdf)          = mzmlhd (irep,1:mxrhdf)
      mzmlon (irep-ndiff,1:maxarl)          = mzmlon (irep,1:maxarl)
    ENDIF
  ENDDO
  ntotaml = ntotaml - ndiff

  DEALLOCATE ( imlls , STAT=nzerr )

  IF (MAX(ndiff,nexceair(1)) > 0)                                              &
    PRINT *,'CAUTION for maxmlo: nexceair ', nexceair(1), ndiff


!-------------------------------------------------------------------------------
!  Section 5  : Adaptation of the vertical correlation scales
!-------------------------------------------------------------------------------

  ALLOCATE ( ismls  (ntotaml) , STAT=nzerr )

! set passive bits in list 'tsl' for reports that are passive in the 'long' list

  IF (num_compute > 1) THEN
    DO irep = 1 , ntottsl
      iairls (irep,2,1)  =  iarlls (mzslon(irep),2)
    ENDDO
  ENDIF

  loop_over_station_ids2: DO ista = 1 , nstaid

! select station id's for each node
! ---------------------------------

    IF (MOD( ista + istoffs , num_compute ) == my_cart_id) THEN
      IF (num_compute > 1) THEN
        ncandls = -1
        ncande  = ntottsl
      ELSE
        ncandls = -8
        ncande  = ntotlon
      ENDIF
      ntotssl = 0

! make lists of single-level resp. multi-level reports with present station id
! ----------------------------------------------------------------------------

      CALL obs_air_list ( ncandls, 1, ncande, 2, ista, ista, 3 , ntotssl )
!     =================

      ntotsml = 0
      DO irep = 1 , ntotaml
        IF (yarlls(mzmlon(irep,1)) == yairls(ista,2)) THEN
          ntotsml = ntotsml + 1
          ismls (ntotsml) = irep
        ENDIF
      ENDDO

! non-MPP runs: store only the reports with present sta. id in temporary array
! ----------------------------------------------------------------------------

      IF (num_compute == 1) THEN
        DO irep = 1 , ntotssl
          nlon  = iairls (irep,1,3)
          nsga  = iarlls (nlon,1)
          DO mvar = 1 , mxshed
            zsghed (irep,mvar) = osghed (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxshdf
            mzsghd (irep,mvar) = mosghd (nsga,mvar)
          ENDDO
          DO mvar = 1 , mxsbdy
            zsgbdy (irep,mvar) = osgbdy (nsga,mvar)
          ENDDO
          mzslon (irep) = nlon
! adjust 'iairls(:,1,3)' so that it contains the same as for MPP runs:
! the report index in the temporary arrays 'zsghed', 'mzsghd', 'zsgbdy', 'mzsgbd
          iairls (irep,1,3) = irep
        ENDDO
      ENDIF

! adjust the vertical correlation scales
! --------------------------------------

      IF (lsvcorl) CALL obs_air_correl ( ntotsml , ntotssl                     &
                                       , vscale, tscale, hscale_c, hscale_v)
!                  ===================

! non-MPP runs: store scale corrections for the reports with present station id
! -----------------------------------------------------------------------------

      IF (num_compute == 1) THEN
        DO irep = 1 , ntotssl
          nlon  = mzslon (irep)
          nsga  = iarlls (nlon,1)
! the list 'ssl' now contains only active reports
! ==> reports set passive on the 'long' list are not on list 'ssl' anyway
!         IF (iarlls(nlon,2) > 0) THEN
            osghed (nsga,nhvcbu) = zsghed (irep,nhvcbu)
            osghed (nsga,nhvcbt) = zsghed (irep,nhvcbt)
            osghed (nsga,nhvcbq) = zsghed (irep,nhvcbq)
            osghed (nsga,nhvctu) = zsghed (irep,nhvctu)
            osghed (nsga,nhvctt) = zsghed (irep,nhvctt)
            osghed (nsga,nhvctq) = zsghed (irep,nhvctq)
!         ENDIF
        ENDDO
      ENDIF

    ENDIF
  ENDDO loop_over_station_ids2

  DEALLOCATE ( ismls   , STAT=nzerr )
  DEALLOCATE ( iairls  , STAT=nzerr )
  DEALLOCATE ( yairls  , STAT=nzerr )


!-------------------------------------------------------------------------------
!  Section 6  : Updating of the ODR for the single-level aircraft reports,
!               and of the counters for statistics
!-------------------------------------------------------------------------------
!  Section 6a : MPP-runs: Scattering and gathering of the passive bit and
!                         vertical correlation scale adjustment factors and
!                         and updating the (single-level) ODR by these factors;
!               (for non-MPP runs, this has already be done in section 5);
!               initiialise 'ntotlsr'
!-------------------------------------------------------------------------------

  IF (num_compute > 1) THEN

! Store the scales and passive bits (including the 'long' list index) in a
! sending buffer, sorted according to the target nodes of the reports.
! Note, that - the local reports 'zsghed' are already sorted (since the target
!              nodes are identical to the source nodes when 'zsghed' has been
!              received at the current node
!            - hence, the number of local reports to be sent to the different
!              target nodes is already known (see below)
! ------------------------------------------------------------------------------

    ALLOCATE (rsendbuf (6 , ndsizet) , STAT=nzerr)
    ALLOCATE (isendbuf (2 , ndsizet) , STAT=nzerr)

    DO irep = 1 , ntottsl
      nlon = mzslon (irep)
      rsendbuf (1,irep) = zsghed (irep,nhvcbu)
      rsendbuf (2,irep) = zsghed (irep,nhvcbt)
      rsendbuf (3,irep) = zsghed (irep,nhvcbq)
      rsendbuf (4,irep) = zsghed (irep,nhvctu)
      rsendbuf (5,irep) = zsghed (irep,nhvctt)
      rsendbuf (6,irep) = zsghed (irep,nhvctq)
      isendbuf (1,irep) = nlon
      isendbuf (2,irep) = iarlls (nlon,2)
    ENDDO

    DEALLOCATE ( zsghed , STAT=nzerr )
    DEALLOCATE ( mzsghd , STAT=nzerr )
    DEALLOCATE ( zsgbdy , STAT=nzerr )
    DEALLOCATE ( mzsgbd , STAT=nzerr )
    DEALLOCATE ( mzslon , STAT=nzerr )

! Get the numbers of local reports to be sent to the different target nodes,
! and the numbers of reports to be received from the different source nodes
! --------------------------------------------------------------------------

    DO inode = 1, num_compute
      ndiff             =  ircount (inode)
      ircount  (inode)  =  iscount (inode)
      iscount  (inode)  =  ndiff
    ENDDO

! Scatter and gather the single-level reports
! -------------------------------------------

    ALLOCATE (rrecvbuf (6 , ndsizel) , STAT=nzerr)
    ALLOCATE (irecvbuf (2 , ndsizel) , STAT=nzerr)

    CALL reports_all2all ( iscount , ircount , num_compute                     &
                         , ndsizet , isendbuf, rsendbuf , 2      , 6           &
                         , ndsizel , irecvbuf, rrecvbuf , yerrmsg, izmplcode )
  ! ====================

    IF (izmplcode /= 0)                                                        &
      CALL model_abort ( my_cart_id, 11707, yerrmsg, yroutine, izmplcode )

    DEALLOCATE (rsendbuf , STAT=nzerr)
    DEALLOCATE (isendbuf , STAT=nzerr)

! store the received vertical correlation scale adjustment factors in ODR
! -----------------------------------------------------------------------

    ntotlsr = 0
    DO inode = 1, num_compute
      ntotlsr = ntotlsr + ircount(inode)
    ENDDO
    IF ((ntotlsr /= ntotlsl) .AND. (lwonl))                                    &
      WRITE( nupr,'(''CAUTION: node'',I4,'' : # received reports'',I5,'' /=''  &
                  &,I4,'' (# reports sent originally'')' )                     &
             my_cart_id, ntotlsr, ntotlsl
    IF  (ntotlsr /= ntotlsl)                                                   &
      PRINT       '(''CAUTION: node'',I4,'' : # received reports'',I5,'' /=''  &
                  &,I4,'' (# reports sent originally'')' ,                     &
             my_cart_id, ntotlsr, ntotlsl

    DO irep = 1 , ntotlsr
      nlon  = irecvbuf(1,irep)
      nsga  = iarlls (nlon,1)
      IF ((iarlls(nlon,3) /= my_cart_id) .AND. (lwonl))                        &
        WRITE( nupr,'(''CAUTION: rep. no'',I5,'' id '',A ,'' of node'',I4      &
                    &,'' arrived at node'',I4)' )                              &
               nlon, yarlls(nlon), iarlls(nlon,3), my_cart_id
      IF  (iarlls(nlon,3) /= my_cart_id)                                       &
        PRINT       '(''CAUTION: rep. no'',I5,'' id '',A ,'' of node'',I4      &
                    &,'' arrived at node'',I4)' ,                              &
               nlon, yarlls(nlon), iarlls(nlon,3), my_cart_id
      IF (irecvbuf(2,irep) >= 1) THEN
        osghed (nsga,nhvcbu) = rrecvbuf (1,irep)
        osghed (nsga,nhvcbt) = rrecvbuf (2,irep)
        osghed (nsga,nhvcbq) = rrecvbuf (3,irep)
        osghed (nsga,nhvctu) = rrecvbuf (4,irep)
        osghed (nsga,nhvctt) = rrecvbuf (5,irep)
        osghed (nsga,nhvctq) = rrecvbuf (6,irep)
      ENDIF
    ENDDO

  ELSE    !  num_compute = 1

    DEALLOCATE ( zsghed , STAT=nzerr )
    DEALLOCATE ( mzsghd , STAT=nzerr )
    DEALLOCATE ( zsgbdy , STAT=nzerr )
    DEALLOCATE ( mzsgbd , STAT=nzerr )
    DEALLOCATE ( mzslon , STAT=nzerr )

    ntotlsr = ntotlon
  ENDIF

!-------------------------------------------------------------------------------
!  Section 6b : Updating of report status and flags, and counters for statistics
!-------------------------------------------------------------------------------

  DO irep = 1 , ntotlsr

    IF (num_compute > 1) THEN
      nlon  = irecvbuf(1,irep)
      nsga  = iarlls (nlon,1)
      npassiv = 1 - MIN( irecvbuf(2,irep) , 1 )
    ELSE    !  num_compute = 1
      nlon  = irep
      nsga  = iarlls (nlon,1)
      npassiv = 1 - MIN( iarlls(nlon,2) , 1 )
    ENDIF
! passive reports: 1 : single-level report used for multi-level report
!   'npassiv'  =   2 : passive (e.g. outside specified obs. area, etc.)
!                  3 : flagged by the flight track check
!                  4 : flagged by thinning
!                  0 : active report

    ! update report status counters and report event counters
    IF (      (mosghd(nsga,nhpass) /= MIN( npassiv,2 ))                        &
        .AND. ((lfirst) .OR. (osghed(nsga,nhtime) >= acthr-epsy))) THEN
      kcdtyp = mosghd (nsga,nhcode)
      kobtyp = mosghd (nsga,nhobtp)
      icma  = i_cma ( kobtyp , kcdtyp )
      IF ((mosghd(nsga,nhpass) <= 1) .AND. (npassiv >= 2)) THEN
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      ELSEIF ((mosghd(nsga,nhpass) == 2) .AND. (npassiv <= 1)) THEN
        cma(icma)%cnt_ac = cma(icma)%cnt_ac + 1
        cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      ENDIF
      ! report event counters
      IF ((mosghd(nsga,nhpass) /= 1) .AND. (npassiv == 1))                     &
        neventr (neslps,icma) = neventr (neslps,icma) + 1
      IF ((mosghd(nsga,nhpass) == 1) .AND. (npassiv /= 1))                     &
        neventr (neslps,icma) = neventr (neslps,icma) - 1
      IF (      (BTEST( mosghd(nsga,nhflag), FL_FLIGHTTRACK ))                 &
          .AND. (npassiv <= 1))                                                &
        neventr (netrac,icma) = neventr (netrac,icma) - 1
      IF ((mosghd(nsga,nhpass) <= 1) .AND. (npassiv == 3))                     &
        neventr (netrac,icma) = neventr (netrac,icma) + 1
      IF ((mosghd(nsga,nhpass) <= 1) .AND. (npassiv == 4))                     &
        neventr (nethin,icma) = neventr (nethin,icma) + 1
    ENDIF

    ! update report status and report flags
    IF (mosghd(nsga,nhpass) /= MIN( npassiv,2 )) THEN
      IF (      (BTEST( mosghd(nsga,nhflag), FL_FLIGHTTRACK ))                 &
          .AND. (npassiv <= 1)) THEN
        ! if the flight track flag has to be cancelled
        mosghd (nsga,nhschr) = IBCLR ( mosghd(nsga,nhschr), nvhtbp )
        mosghd (nsga,nhschr) = IBCLR ( mosghd(nsga,nhflag), FL_FLIGHTTRACK )
      ENDIF
      IF (npassiv == 0)  mosghd (nsga,nhpass) = 0
      IF (npassiv == 1)  mosghd (nsga,nhpass) = 1
      IF (npassiv == 3) THEN
        mosghd (nsga,nhpass) = 2
        mosghd (nsga,nhschr) = IBSET ( mosghd(nsga,nhschr), nvhtbp )
        mosghd (nsga,nhflag) = IBSET ( mosghd(nsga,nhflag), FL_FLIGHTTRACK )
      ENDIF
      IF (npassiv == 4) THEN
        mosghd (nsga,nhpass) = 2
        mosghd (nsga,nhschr) = IBSET ( mosghd(nsga,nhschr), nvhhbp )
        mosghd (nsga,nhflag) = IBSET ( mosghd(nsga,nhflag), FL_THIN )
      ENDIF
    ENDIF
  ENDDO

  IF (num_compute > 1)  DEALLOCATE (rrecvbuf , STAT=nzerr)
  IF (num_compute > 1)  DEALLOCATE (irecvbuf , STAT=nzerr)

  DEALLOCATE ( iarlls  , STAT=nzerr )
  DEALLOCATE ( rarlls  , STAT=nzerr )


!-------------------------------------------------------------------------------
!  Section 7  : Distribution and storage in the ODR
!               of the aircraft multi-level reports
!-------------------------------------------------------------------------------

! For MMP-runs
! ------------

  IF (num_compute > 1) THEN

! Store the local multi-level reports in a sending buffer, sorted according to
! the target nodes of the reports.
! Determine the number of local reports to be sent to the different target nodes
! ------------------------------------------------------------------------------

    nelemr = mxrhed + maxarl *mxrbdy
    nelemi = mxrhdf + maxarl *mxrbdf + 1
    ALLOCATE ( rsendbuf (nelemr , maxmll)            , STAT=nzerr )
    ALLOCATE ( isendbuf (nelemi , maxmll)            , STAT=nzerr )
    ALLOCATE ( alltoall (num_compute * num_compute)  , STAT=nzerr )

    rsendbuf (:,:) = rmdi
    isendbuf (:,:) = imdi
    alltoall   (:) = 0
    nrep           = 0
!   DO inode = 0, num_compute-1
!     alltoall (my_cart_id*num_compute+inode+1)  =  0
!   ENDDO

    IF (ntotaml > 0) THEN
      ALLOCATE ( irprcs (ntotaml+1) , STAT=nzerr )
      ALLOCATE ( irnode (ntotaml+1) , STAT=nzerr )
      ALLOCATE ( irsort (ntotaml+1) , STAT=nzerr )
      DO irps = 1 , ntotaml
        irprcs (irps) = irps
      ENDDO

      CALL obs_assign_sort_node ( maxmll, ntotaml, irprcs, mzmlhd(:,nhitot)    &
                                , mzmlhd(:,nhjtot), num_compute, i_subpos       &
                                , nboundlines, my_cart_id , irnode, irsort )
!     =========================

    ENDIF
    DO nrep = 1 , ntotaml
      irep = irsort(nrep)
      DO mvar = 1 , mxrhed
        rsendbuf (mvar,nrep) = zmlhed (irep,mvar)
      ENDDO
      DO mvar = 1 , mxrhdf
        isendbuf (mvar,nrep) = mzmlhd (irep,mvar)
      ENDDO
      DO ilev = 1 , mzmlhd(irep,nhnlev)
        iante = mxrhed + (ilev-1)*mxrbdy
        DO mvar = 1 , mxrbdy
          rsendbuf (iante+mvar,nrep) = zmlbdy (irep,ilev,mvar)
        ENDDO
        iante = mxrhdf + (ilev-1)*mxrbdf
        DO mvar = 1 , mxrbdf
          isendbuf (iante+mvar,nrep) = mzmlbd (irep,ilev,mvar)
        ENDDO
      ENDDO
      isendbuf (mxrhdf+maxarl*mxrbdf+1,nrep)  =  mzmlon (irep,1)
      alltoall (my_cart_id*num_compute+irnode(nrep)+1)  =                      &
      alltoall (my_cart_id*num_compute+irnode(nrep)+1)  +  1
    ENDDO
    IF (ntotaml > 0) THEN
      DEALLOCATE ( irprcs , STAT=nzerr )
      DEALLOCATE ( irnode , STAT=nzerr )
      DEALLOCATE ( irsort , STAT=nzerr )
    ENDIF

    DEALLOCATE ( zmlhed , STAT=nzerr )
    DEALLOCATE ( mzmlhd , STAT=nzerr )
    DEALLOCATE ( zmlbdy , STAT=nzerr )
    DEALLOCATE ( mzmlbd , STAT=nzerr )
    DEALLOCATE ( mzmlon , STAT=nzerr )

! Get the numbers of local reports to be sent to the different target nodes,
! and the numbers of reports to be received from the different source nodes
! --------------------------------------------------------------------------

    CALL global_values ( alltoall, num_compute*num_compute, 'MAX'              &
                       , imp_integers, icomm_cart, -1, yerrmsg, nzerr )
!   ==================

    DO inode = 1, num_compute
      iscount (inode) =  alltoall (my_cart_id*num_compute +inode       )
      ircount (inode) =  alltoall ((inode-1) *num_compute +my_cart_id+1)
    ENDDO
    DEALLOCATE (alltoall , STAT=nzerr)

! Scatter and gather the multi-level reports
! ------------------------------------------
    ALLOCATE ( rrecvbuf (nelemr , maxmll)            , STAT=nzerr )
    ALLOCATE ( irecvbuf (nelemi , maxmll)            , STAT=nzerr )

    CALL reports_all2all ( iscount , ircount , num_compute                     &
                         , maxmll  , isendbuf, rsendbuf , nelemi , nelemr      &
                         , maxmll  , irecvbuf, rrecvbuf , yerrmsg, izmplcode )
  ! ====================

    IF (izmplcode /= 0)                                                        &
      CALL model_abort ( my_cart_id, 11712, yerrmsg, yroutine, izmplcode )

! Store the received multi-level reports in the local ODR
! -------------------------------------------------------

    ntotaml = 0
    DO inode = 1, num_compute
      ntotaml = ntotaml + ircount(inode)
    ENDDO
    DO irep = 1 , ntotaml
      ntotml = ntotml + 1
      DO mvar = 1 , mxrhed
        omlhed (ntotml,mvar) = rrecvbuf (mvar,irep)
      ENDDO
      DO mvar = 1 , mxrhdf
        momlhd (ntotml,mvar) = irecvbuf (mvar,irep)
      ENDDO
!     DO ilev = 1 , maxarl
      DO ilev = 1 , momlhd(ntotml,nhnlev)
        iante = mxrhed + (ilev-1)*mxrbdy
        DO mvar = 1 , mxrbdy
          omlbdy (ntotml,ilev,mvar) = rrecvbuf (iante+mvar,irep)
        ENDDO
        iante = mxrhdf + (ilev-1)*mxrbdf
        DO mvar = 1 , mxrbdf
          momlbd (ntotml,ilev,mvar) = irecvbuf (iante+mvar,irep)
        ENDDO
!       omlbdy (ntotml,ilev,nbtzio) =  omlbdy(ntotml,ilev,nbtzio)              &
!                                   - i_subpos(my_cart_id,1) + nboundlines + 1
!       omlbdy (ntotml,ilev,nbtzjo) =  omlbdy(ntotml,ilev,nbtzjo)              &
!                                   - i_subpos(my_cart_id,2) + nboundlines + 1
      ENDDO
!     momlhd (ntotml,nhio) = momlhd(ntotml,nhitot) - i_subpos(my_cart_id,1)     &
!                                                + nboundlines + 1
!     momlhd (ntotml,nhjo) = momlhd(ntotml,nhjtot) - i_subpos(my_cart_id,2)     &
!                                                + nboundlines + 1
      yomlhd (ntotml)  =  yarlls (irecvbuf(mxrhdf+maxarl*mxrbdf+1,irep))
      kcdtyp = momlhd (ntotml,nhcode)
      kobtyp = momlhd (ntotml,nhobtp)
      icma  = i_cma ( kobtyp , kcdtyp )
      neventr (neslml,icma) = neventr (neslml,icma) + 1
    ENDDO

    DEALLOCATE (rsendbuf , STAT=nzerr)
    DEALLOCATE (isendbuf , STAT=nzerr)
    DEALLOCATE (rrecvbuf , STAT=nzerr)
    DEALLOCATE (irecvbuf , STAT=nzerr)

    DEALLOCATE ( yarlls , STAT=nzerr )

! For non-MPP runs: simply store the multi-level reports in the ODR
! -----------------------------------------------------------------

  ELSE    !  num_compute = 1

    DO irep = 1 , ntotaml
      ntotml = ntotml + 1
      DO mvar = 1 , mxrhed
        omlhed (ntotml,mvar) = zmlhed (irep,mvar)
      ENDDO
      DO mvar = 1 , mxrhdf
        momlhd (ntotml,mvar) = mzmlhd (irep,mvar)
      ENDDO
      DO ilev = 1 , maxarl
        DO mvar = 1 , mxrbdy
          omlbdy (ntotml,ilev,mvar) = zmlbdy (irep,ilev,mvar)
        ENDDO
        DO mvar = 1 , mxrbdf
          momlbd (ntotml,ilev,mvar) = mzmlbd (irep,ilev,mvar)
        ENDDO
      ENDDO
      yomlhd (ntotml)  =  yarlls (mzmlon(irep,1))
      kcdtyp = momlhd (ntotml,nhcode)
      kobtyp = momlhd (ntotml,nhobtp)
      icma  = i_cma ( kobtyp , kcdtyp )
      neventr (neslml,icma) = neventr (neslml,icma) + 1
    ENDDO

    DEALLOCATE ( zmlhed , STAT=nzerr )
    DEALLOCATE ( mzmlhd , STAT=nzerr )
    DEALLOCATE ( zmlbdy , STAT=nzerr )
    DEALLOCATE ( mzmlbd , STAT=nzerr )
    DEALLOCATE ( mzmlon , STAT=nzerr )

    DEALLOCATE ( yarlls , STAT=nzerr )

  ENDIF

ENDIF   !  (nstaid > 0)

!-------------------------------------------------------------------------------
! End of module procedure obs_air_org_mult
!-------------------------------------------------------------------------------

END SUBROUTINE obs_air_org_mult


!-------------------------------------------------------------------------------
!+ Module procedure for production of multi-level aircraft reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_air_make_mult ( ntotssl , thairh , ntotaml , nexcloc )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure produces multi-level aircraft reports from
!   single-level reports with one and the same flight number.
!   In addition, the flight phase (ascent, descent, or level flight) for each
!   report is determined, and checks for exaggerated colocation of the reports
!   as well as a general flight track check are performed. Finally, sequences of
!   frequent single-level reports at approximately the same pressure level are
!   thinned.
!   For Mode-S data, it is found that data with the same 'station id' are not
!   usually produced by one aircraft or flight. Therefore, everything is skipped
!   except for the check for exaggerated colocation and the thinning.
!
! Method:
!   1. The reports are sorted in a unique order according to the flight track.
!      This also results in independency from the domain decomposition.
!   2. The flight phases are determined by a forward and a backwards sweep
!      through the list of reports. Criteria are as given in section 2.
!   3. Checks for exaggerated horizontal or vertical colocation are made
!      independently.
!   4. The horizontal flight track check follows the ECMWF pre-processing
!      (Met. Bull. M1.4/3) roughly except for extension of the reference
!      distance used. It is complemented by a iterative check for missing
!      sign at the reported longitude, and by a vertical check.
!   5. The criteria for the production of multi-level reports are as follows:
!      - Only strictly 'active' single-level reports are used.
!      - No single-level report may be used twice.
!      - Each multi-level aircraft report is to start at as low a level as
!        possible and is to contain as many single-level reports as possible
!        within the limits as by parameters 'thairt', 'thairv' and namelist
!        variable 'thairt'.
!      - Each multi-level aircraft report must contain at least 'minslml' levels
!        levels (i.e. single-level reports).
!      - The time and horizontal grid point assignment of the multi-level report
!        is according to the single-level report from which the lowest level is
!        derived (similar as surface level (station) for TEMP / PILOT reports).
!   6. Thinning of vertically quasi-colocated (single-level) reports.
!
!   Note that unlike in the threshold quality control, reports set passive by
!   the flight track check, the colocation check, or the thinning are never
!   given the chance to become active again.
!   Here, in list 'iarlls', the code type of these reports is replaced by '-2'
!   for the track checks resp. '-3' for thinning, and the code type of reports
!   used for multi-level reports is replaced by 0.
!   27.01.98
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  Second part of flags replaced by threshold quality control flags.
! 1.27       1999/03/29 Christoph Schraff
!  Adjustments due to the modified ODR flag formats.
! 1.31       1999/07/01 Christoph Schraff
!  Removal of 'lreproduce'.
! 2.5        2001/06/01 Christoph Schraff
!  Introduction of the flight track check(s), the checks for exaggerated
!  colocation, the thinning, and of the determination of the flight phases.
! 2.6        2001/06/12 Christoph Schraff
!  Tuning of flight track check thresholds. Message on thinning aircraft reports.
! 2.12       2001/11/07 Christoph Schraff
!  Flight track check done also without production of multi-level aircraft rep.
! 2.13       2002/01/18 Christoph Schraff
!  Bug correction: savety limits for array indices.
! 3.6        2003/12/11 Christoph Schraff
!  Improved warning messages on horiz./vertic. collocation of reports.
!  Ensuring error message at model_abort.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments:
! --------------------

  INTEGER,        INTENT (IN)         ::       &
    ntotssl          ! number of single-level reports from the flight being
                     ! processed currently (i.e. with the same flight number)

  REAL (KIND=wp), INTENT (IN)         ::       &
    thairh           ! maximum horizontal distance for combining single level
                     ! aircraft reports to a multi-level report

  INTEGER,        INTENT (INOUT)      ::       &
    ntotaml       ,& ! number of multi-level aircraft reports produced
                     ! at the current timestep and node
    nexcloc          ! number of reports in excess of local array size

! Local parameters: None
! ----------------

  REAL    (KIND=wp   ),     PARAMETER           ::       &
    clon  = 10000.0_wp ,& ! latitudinal distance scalor
    c5000 =  5000.0_wp ,&
    c2500 =  2500.0_wp ,&
    c100  =   100.0_wp ,&
    c90   =    90.0_wp ,&
    c85   =    85.0_wp ,&
    c80   =    80.0_wp ,&
    c65   =    65.0_wp ,&
    c60   =    60.0_wp ,&
    c50   =    50.0_wp ,&
    c40   =    40.0_wp ,&
    c30   =    30.0_wp ,&
    c12p5 =    12.5_wp ,&
    c5    =     5.0_wp ,&
    c3    =     3.0_wp ,&
    c10r  =     0.1_wp

! Local scalars:
! -------------

  INTEGER         ::  &
    irep             ,& ! loop index over reports
    ivss             ,& ! loop index over group of reports
    nvss             ,& ! number of (temporal) group of reports
    kobtyp , kcdtyp  ,& ! observation type , observation code type
    ifb    , ife     ,& ! boolean factors used to compare vertical distances
    ireprb , irepre  ,& ! indices of reports to which the distance is measured
    idif1  , idif2   ,& ! non-zero distance to neighboring reports
    nob1   , nob2    ,& ! ODR indices of rep. for which distances are compared
    nobe   , nobb    ,& ! ODR indices of rep. to  which distances are compared
    izpmin , izpmax  ,& ! pressure range for current group of reports
    isame  , isamec  ,& ! number of partially colocated reports (12 min apart)
    itimee           ,& ! time of last colocated report
    iqcl             ,& ! (forward - backward) loop index for flight track check
    irepcm1, irepcm2 ,& ! confidence reports used to estimate current position
    irepcm3          ,& ! confidence report  used to estimate current position
    nest             ,& ! max. number of position estimates
    nob              ,& ! ODR index of report being flight track checked
    nobm1  , nobm2   ,& ! ODR indices of reports used to estimate the position
    irepm1 , irepm2  ,& ! indices of reports used to estimate the position
    nactssl          ,& ! number of available active reports (after the checks)
    nusedssl         ,& ! number of reports already processed
    nssl             ,& ! index in list 'ssl' (current station id)
    itsl             ,& ! index in temporary single-level report array
    nlon             ,& ! index in 'long' list
    ienext , jenext  ,& ! coord. of the next above the lowest active candidates
    nsnext           ,& ! number of the next above the lowest active candidates
    nlev             ,& ! number of levels in multi-level report
    klev             ,& ! vertical loop index
    itmp             ,& ! temporary station index in sorted list
    icma             ,& ! indices for statistics ('cma', events)
    ilen   , icl     ,& ! length of output record
    mprout           ,& ! mode of printing rejection message in
                        !   'obs_cdf_print_reject'
    nzerr               ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    thairh2          ,& ! (square of) max. horizontal distance [km2] between any
                        ! level and the base level of a multi-level air. report
    rdegkm2          ,& ! (square of) conversion factor degrees to 'km'
    zlat             ,& ! latitudinal distance scalor
    fidexy , fidez   ,& ! confidence  (horizontally / vertically)
    fidetw , fidew   ,& ! weight for (current / total sum of) confidence
    fidex2 , fidet2  ,& ! confidence / weight, with reversed sign for longitude
    fexpol           ,& ! degree of extrapolation
    xestim , yestim  ,& ! estimate of horizontal position
    zestim           ,& ! estimate of vertical position
    estdxy , estdz   ,& ! distance    current report  -  estimate (horiz / vert)
    estdx2           ,& ! distance to estimate (longitude: with reversed sign)
    refdxy , refdz   ,& ! distance confidence report  -  estimate (horiz / vert)
    exrefd           ,& ! confidence factor, used to extend the reference dist.
    fidxpr , fidzpr  ,& ! current (unweighted) confidence
    fidxp2           ,& ! current (unweighted) confidence (longit.: rever. sign)
    pslmax           ,& ! pressure of lowest active candidate
    timob1           ,& ! time     of lowest active candidate
    psnext           ,& ! pressure of the next above the lowest active candidate
    tmnext           ,& ! time     of the next above the lowest active candidate
    tmdist           ,& ! time  distance betw. lowest and next active candidates
    r2dist           ,& ! horiz distance betw. lowest and next active candidates
    obtimeb          ,& ! observation time at base of multi-level report
    obilonb, objlatb ,& ! longitude / latitude at base of multi-level report
    cosolat          ,& ! cos (latitude of report)
    rhzob            ,& ! (square of) horizontal distance between reports
    obbaltb          ,& ! pressure at base of multi-level report
    obbaltl          ,& ! pressure at top  of multi-level report so far
    pdiffmn          ,& ! pressure distance betw. top report and current report
    obtimec          ,& ! observation time of current report
    obilonc, objlatc ,& ! longitude / latitude of current report
    obbaltc          ,& ! pressure of current report
    rvtlmt              ! colocation threshold for vertical distance

  LOGICAL                  ::  &
    lmodes           ,& ! the current aircraft number relates to Mode-S data
    lchange          ,& ! the order of at least one pair of stations is changed
    lundef           ,& ! flight phase still undefined
    lchktrc          ,& ! check flight track
    lrejhor, lrejver ,& ! set report to passive (horiz / vertic. colocation)
    lrejrep          ,& ! set report to passive (if flight track check fails)
    lupdate, lupdat2 ,& ! update the total confidence by the new estimate
    lsamealt         ,& ! lowest active and present candidate have same altitude
    lobnext          ,& ! current report is the next above the lowest candidate
    lastfind         ,& ! another report was added to the multi-level report
    lactive             ! present report has not yet been processed

  CHARACTER (LEN=ilstidp)  ::  &
    ystid               ! station identity
  CHARACTER (LEN=3)        ::  &
    ystid3              ! first 3 letters of station identity
  CHARACTER (LEN=75) yerrmsg    ! error message
  CHARACTER (LEN=25) yroutine   ! name of this subroutine

! Local arrays:
! ------------

  INTEGER  ::  &
    iasmls  (maxarl,3)  ! lists of single-level reports used for
                        ! one multi-level report

  INTEGER, ALLOCATABLE :: &
    isortsl    (:)   ,& ! sorted list of reports with same station id
    irepb      (:)   ,& ! index of first report within a group
    irepe      (:)   ,& ! index of last  report within a group
    irepbb     (:)   ,& ! index of last  report before a group
    irepee     (:)   ,& ! index of first report before a group
    mphase     (:)   ,& ! phase of flight assigned to report
    izpp       (:)   ,& ! pressure of report
    iflftc     (:)   ,& ! flight track check flag
    mphaso     (:)      ! phase of flight as read from AOF

  REAL    (KIND=wp   )      , ALLOCATABLE :: &
    fideh      (:,:) ,& ! horizontal confidence
    fidev      (:,:) ,& ! vertical confidence
    zflat      (:)      ! cos (latitude of report)
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_air_make_mult
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1  : Sort reports for the quality control (flight tracking)
!               and for independency from the domain decomposition
!-------------------------------------------------------------------------------

! preliminaries
! -------------

  yroutine = 'obs_air_make_mult'
  thairh2  = thairh **2
  rdegkm2  = 111.1_wp **2

  ystid  = yarlls (mzslon(iairls(1,1,3))) (1:ilstidp)
  ystid3 = ystid (1:3)

  IF (ntotssl > 0) THEN
    nzerr = 0
    IF (nzerr == 0)  ALLOCATE ( isortsl (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( irepb   (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( irepe   (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( irepbb  (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( irepee  (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( mphase  (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( izpp    (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( iflftc  (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( mphaso  (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( zflat   (ntotssl) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( fideh   (ntotssl,4) , STAT = nzerr )
    IF (nzerr == 0)  ALLOCATE ( fidev   (ntotssl,2) , STAT = nzerr )
    IF (nzerr /= 0) THEN
      yerrmsg = ' *** Allocation failed ***'
      WRITE( yerrmsg,'(" *** Allocation failed ***",I5)' ) nzerr
      IF (my_cart_id /= 0)  PRINT '(A,2X,A)', yroutine, yerrmsg
      CALL model_abort (my_cart_id, 11888, yerrmsg, yroutine)
!     ================
    ENDIF
  ENDIF

  lmodes = .FALSE.
  DO irep = 1 , ntotssl
    IF (mzsghd(iairls(irep,1,3),nhcode) == nmodes)  lmodes = .TRUE.
  ENDDO

  DO irep = 1 , ntotssl
    isortsl (irep) = irep
  ENDDO

! sort temporally
! ---------------

  IF (ntotssl > 1) THEN
    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO irep = 1 , ntotssl-1
        IF (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )          &
             > NINT( zsghed(iairls(isortsl(irep+1),1,3),nhtime)*c60 )) THEN
          itmp             = isortsl (irep+1)
          isortsl (irep+1) = isortsl (irep)
          isortsl (irep)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO

! sub-sort horizontally (for independency from the domain decomposition only)
! ---------------------------------------------------------------------------

! sort longitudinally (stations are pre-sorted due to their origin nodes)
    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO irep = 1 , ntotssl-1
        IF (      (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )   &
                   == NINT( zsghed(iairls(isortsl(irep+1),1,3),nhtime)*c60 ))  &
            .AND. (   mzsghd(iairls(isortsl(irep  ),1,3),nhitot)               &
                    > mzsghd(iairls(isortsl(irep+1),1,3),nhitot))) THEN
          itmp             = isortsl (irep+1)
          isortsl (irep+1) = isortsl (irep)
          isortsl (irep)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO
! sort latitudinally
    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO irep = 1 , ntotssl-1
        IF (      (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )   &
                   == NINT( zsghed(iairls(isortsl(irep+1),1,3),nhtime)*c60 ))  &
            .AND. (   mzsghd(iairls(isortsl(irep  ),1,3),nhitot)               &
                   == mzsghd(iairls(isortsl(irep+1),1,3),nhitot))              &
            .AND. (   mzsghd(iairls(isortsl(irep  ),1,3),nhjtot)               &
                    > mzsghd(iairls(isortsl(irep+1),1,3),nhjtot))) THEN
          itmp             = isortsl (irep+1)
          isortsl (irep+1) = isortsl (irep)
          isortsl (irep)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO
  ENDIF

! re-sub-sort vertically, and sub-sub-sort horizontally
! -----------------------------------------------------

  IF (ntotssl >= 3) THEN
! group reports with same observation time
    nvss  = 1
    irepb = 0
    DO irep = 1 , ntotssl-1
      IF (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )            &
          == NINT( zsghed(iairls(isortsl(irep+1),1,3),nhtime)*c60 )) THEN
        IF (irepb(nvss) == 0) THEN
          irepb  (nvss) = irep
          irepbb (nvss) = 0
          IF (irep > 1) THEN
            IF (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )      &
                <= NINT( zsghed(iairls(isortsl(irep-1),1,3),nhtime)*c60 )+15)  &
              irepbb (nvss) = irep - 1
          ENDIF
        ENDIF
        irepe  (nvss) = irep + 1
        irepee (nvss) = irep + 2
        IF (irep == ntotssl-1)  irepee (nvss) = 0
        IF (irep == ntotssl-1)  nvss = nvss + 1
      ELSEIF (irepb(nvss) /= 0) THEN
        IF (   NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )          &
            <  NINT( zsghed(iairls(isortsl(irep+1),1,3),nhtime)*c60 )-15)      &
          irepee (nvss) = 0
        nvss  = nvss + 1
      ENDIF
    ENDDO
    nvss  = nvss - 1
! within each group, prepare decision to exchange reports for vertical sorting
    DO ivss = 1 , nvss
      ifb    = MIN( 1 , irepbb(ivss) )
      ife    = MIN( 1 , irepee(ivss) )
      ireprb = ifb *(irepbb(ivss) - irepb(ivss)) + irepb(ivss)
      irepre = ife *(irepee(ivss) - irepe(ivss)) + irepe(ivss)
      lchange = .TRUE.
      DO WHILE (lchange)
        lchange = .FALSE.
        DO irep = irepb(ivss) , irepe(ivss)-1
          idif1 = NINT( REAL(ife,wp)* ABS( zsgbdy(iairls(isortsl(irep  ),1,3),nbsp)     &
                                          -zsgbdy(iairls(isortsl(irepre),1,3),nbsp) )   &
                       -REAL(ifb,wp)* ABS( zsgbdy(iairls(isortsl(irep  ),1,3),nbsp)     &
                                          -zsgbdy(iairls(isortsl(ireprb),1,3),nbsp) ) )
          idif2 = NINT( REAL(ife,wp)* ABS( zsgbdy(iairls(isortsl(irep+1),1,3),nbsp)     &
                                          -zsgbdy(iairls(isortsl(irepre),1,3),nbsp) )   &
                       -REAL(ifb,wp)* ABS( zsgbdy(iairls(isortsl(irep+1),1,3),nbsp)     &
                                          -zsgbdy(iairls(isortsl(ireprb),1,3),nbsp) ) )
! for vertical colocation, prepare decision to exchange for horizontal sorting
          IF (idif1 == idif2) THEN
            nob1  = iairls(isortsl(irep  ),1,3)
            nob2  = iairls(isortsl(irep+1),1,3)
            nobe  = iairls(isortsl(irepre),1,3)
            nobb  = iairls(isortsl(ireprb),1,3)
            zlat  = c100 * COS( (  ABS( zsghed(nob1,nhjlat) )                  &
                                 + ABS( zsghed(nob2,nhjlat) )) * c05 *r_degrad )
            zlat  = zlat * zlat
            idif1 =NINT(ife*( zlat*(zsghed(nob1,nhilon) -zsghed(nobe,nhilon))  &
                                  *(zsghed(nob1,nhilon) -zsghed(nobe,nhilon))  &
                             +clon*(zsghed(nob1,nhjlat) -zsghed(nobe,nhjlat))  &
                                  *(zsghed(nob1,nhjlat) -zsghed(nobe,nhjlat))) &
                       -ifb*( zlat*(zsghed(nob1,nhilon) -zsghed(nobb,nhilon))  &
                                  *(zsghed(nob1,nhilon) -zsghed(nobb,nhilon))  &
                             +clon*(zsghed(nob1,nhjlat) -zsghed(nobb,nhjlat))  &
                                  *(zsghed(nob1,nhjlat) -zsghed(nobb,nhjlat))) )
            idif2 =NINT(ife*( zlat*(zsghed(nob2,nhilon) -zsghed(nobe,nhilon))  &
                                  *(zsghed(nob2,nhilon) -zsghed(nobe,nhilon))  &
                             +clon*(zsghed(nob2,nhjlat) -zsghed(nobe,nhjlat))  &
                                  *(zsghed(nob2,nhjlat) -zsghed(nobe,nhjlat))) &
                       -ifb*( zlat*(zsghed(nob2,nhilon) -zsghed(nobb,nhilon))  &
                                  *(zsghed(nob2,nhilon) -zsghed(nobb,nhilon))  &
                             +clon*(zsghed(nob2,nhjlat) -zsghed(nobb,nhjlat))  &
                                  *(zsghed(nob2,nhjlat) -zsghed(nobb,nhjlat))) )
          ENDIF
! exchange (order of) reports in sorted list, if needed
          IF (idif1 < idif2) THEN
            itmp             = isortsl (irep+1)
            isortsl (irep+1) = isortsl (irep)
            isortsl (irep)   = itmp
            lchange = .TRUE.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 2  : Determine flight phases: ascent, descent, or level flight
!-------------------------------------------------------------------------------

  IF (ntotssl > 0) THEN
! IF ((ntotssl >= 3) .AND. (ystid3 /= "XXX") .AND. (ystid3 /= "???")           &
!                    .AND. (ystid3 /= "///") .AND. (ystid3 /= "BBX")) THEN
    DO irep = 1 , ntotssl
      mphase (irep) = 0
      izpp   (irep) = NINT( zsgbdy(iairls(isortsl(irep),1,3),nbsp) )
    ENDDO
  ENDIF

  IF ((ntotssl > 0) .AND. (.NOT. lmodes)) THEN

! temporal split between flights
! ------------------------------

    nvss = 1
    irepb (nvss) = 1
    irepe (nvss) = 1
    DO irep = 2 , ntotssl
      IF ( NINT( zsghed(iairls(isortsl(irep  ),1,3),nhtime)*c60 )              &
          -NINT( zsghed(iairls(isortsl(irep-1),1,3),nhtime)*c60 ) > 30) THEN
        nvss         = nvss + 1
        irepb (nvss) = irep
        irepe (nvss) = irep
      ELSE
        irepe (nvss) = irep
      ENDIF
    ENDDO

! criteria to determine phases
! ----------------------------
! - within an  ascending phase, descending phases must be less than 10 hPa
!   ==> descents of at least 10 hPa are considered descents
!       (mphase == 6)
! - within an descending phase,  ascending phases must be less than 50 hPa
!   ==>  ascents of at least 50 hPa are considered  ascents
!       (mphase == 5)
! - level flight phases are whenever there are pairs of subsequent reports
!   above 500 hPa which are at most 1 hPa apart from each other, and
!   they include other reports above 500 hPa between such pairs of reports
!       (mphase == 3)

    DO ivss = 1 , nvss
      izpmin = izpp(irepb(ivss))
      izpmax = izpp(irepb(ivss))

! loop: once 'mphase /= 0' is set for one report then 'mphase /= 0' is set for
!       all subsequent reports
      DO irep = irepb(ivss) + 1 , irepe(ivss)
! level flight phases
        IF (      (izpp(irep) < 50000)                                         &
            .AND. (ABS( izpp(irep)-izpp(irep-1) ) <= 100)) THEN
          mphase (irep) = 3
          izpmin = MIN( izpp(irep-1) , izpp(irep) )
! within an undefined phase
        ELSEIF (mphase(irep-1) == 0) THEN
          IF (izpp(irep) >= izpmin+1000) THEN
            mphase (irep) = 6
            izpmax = izpp(irep)
          ELSEIF (izpp(irep) <= izpmax-5000) THEN
            mphase (irep) = 5
            izpmin = izpp(irep)
          ELSE
            izpmin = MIN( izpmin , izpp(irep) )
            izpmax = MAX( izpmax , izpp(irep) )
          ENDIF
! within a level flight phase  ('izpmin' is required)
        ELSEIF (ABS( mphase(irep-1) ) == 3) THEN
          if (izpp(irep) <= izpmin+100) then
!   --> further level flight, or further ascent also labelled as level flight
            mphase (irep) = 3
            izpmin = MIN( izpmin , izpp(irep) )
          ELSEIF (izpp(irep) < izpmin+1000) THEN
!   --> level flight or descent ?
            mphase (irep) = -3
          ELSE
!   --> turn to descent
            mphase (irep) = 6
            izpmax = izpp(irep)
          ENDIF
! within an ascending phase  ('izpmin' is required)
        ELSEIF (ABS( mphase(irep-1) ) == 5) THEN
          IF (izpp(irep) < izpmin) THEN
!   --> further ascent
            mphase (irep) = 5
            izpmin = izpp(irep)
          ELSEIF (izpp(irep) < izpmin+1000) THEN
!   --> 'waiting loop' before further ascent, or descent ?
            mphase (irep) = -5
          ELSE
!   --> turn to descent
            mphase (irep) = 6
            izpmax = izpp(irep)
          ENDIF
! within a descending phase  ('izpmax' is required)
        ELSEIF (ABS( mphase(irep-1) ) == 6) THEN
          IF (izpp(irep) > izpmax) THEN
!   --> further descent
            mphase (irep) = 6
            izpmax = izpp(irep)
          ELSEIF (izpp(irep) > izpmax-5000) THEN
!   --> waiting loop before landing, or re-ascent ?
            mphase (irep) = -6
          ELSE
!   --> turn to re-ascent
            mphase (irep) = 5
            izpmin = izpp(irep)
          ENDIF
        ENDIF
      ENDDO

! once 'mphase /= 0' for 1 report then 'mphase /= 0' for all subsequent reports
      lundef  = (mphase(irepe(ivss)) == 0)
! if no criterion for any phase has been met anywhere,
! then the phase for all reports is set to level flight or descent
      IF ((lundef) .AND. (izpmax <= 50000)) THEN
        DO irep = irepe(ivss) , irepb(ivss) , -1
          mphase (irep) = 3
        ENDDO
      ELSEIF (lundef) THEN
        DO irep = irepe(ivss) , irepb(ivss) , -1
          mphase (irep) = 6
        ENDDO
      ELSE
! if 'mphase < 0' for the last report then no criteria for a phase change is met
        mphase (irepe(ivss))  =  ABS( mphase(irepe(ivss)) )
! if 'mphase <= 0' for any other report (==> ambiguity) then the phase is set
! equal to the phase of the next subsequent report for which the phase could be
! determined (i.e. for which 'mphase > 0')
! (e.g. if report A has 'mphase == -6' then a subsequent report B will probably
!       have either 'mphase == 5' , in which case report A was the beginning of
!       the re-ascent, or 'mphase == 6', in which case report A was part of a
!       'waiting loop' within a descent.)
        DO irep = irepe(ivss) - 1 , irepb(ivss) , -1
          IF (mphase(irep) <= 0)  mphase (irep) = mphase(irep+1)
        ENDDO
      ENDIF

    ENDDO

! prepare writing of flight phase back to the ODR, if it is undefined there
! -------------------------------------------------------------------------

    DO irep = 1 , ntotssl
      mphaso (irep) = IBITS( mzsghd(iairls(isortsl(irep),1,3),nhschr)          &
                           , nvapbp, nvapoc )
! set flight phase to 7=missing , 3,4=flight level , 5=ascending , 6=descending
      IF ((mphaso(irep) <= 2) .OR. (mphaso(irep) >= 7))  mphaso (irep) = 7
      IF (mphaso(irep) == 7) THEN
        mphaso (irep) = mphase(irep) + 8
        CALL MVBITS ( mphaso(irep), 0, nvapoc                                  &
                    , mzsghd(iairls(isortsl(irep),1,3),nhschr), nvapbp )
      ENDIF
!     IF (((MOD(irep,2) == 0) .OR. (irep == ntotssl)) .AND. (acthr <= c05)) THEN
!       itmp = MAX( irep - 1 , 1 )
!       PRINT '(''air flt. phase: '',A ,F6.2,2(3I4,F6.2,F8.0,I4))',ystid,acthr &
!             , itmp, isortsl(itmp), iairls(isortsl(itmp),1,3)                 &
!             , zsghed(iairls(isortsl(itmp),1,3),nhtime)                       &
!             , zsgbdy(iairls(isortsl(itmp),1,3),nbsp)  , mphase(itmp)         &
!             , irep, isortsl(irep), iairls(isortsl(irep),1,3)                 &
!             , zsghed(iairls(isortsl(irep),1,3),nhtime)                       &
!             , zsgbdy(iairls(isortsl(irep),1,3),nbsp)  , mphase(irep)
!     ENDIF
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
!  Section 3  : Check for exaggerated colocation of reports
!-------------------------------------------------------------------------------

  lchktrc = ((ntotssl >= 3) .AND. (ystid3 /= "XXX") .AND. (ystid3 /= "???")    &
                            .AND. (ystid3 /= "///") .AND. (ystid3 /= "BBX"))
  lrejrep = .FALSE.
  IF (lchktrc) THEN

! horizontally: reject all, if at least 3 colocated reports are more than 12 min
! ------------  apart and if all colocated reports comprise at least 50 %
!               of all reports
!               (the time condition is included to account for colocation at
!                ascents or at waiting loops at descents.)

    lrejhor = .FALSE.
    DO irep = 1 , ntotssl-2
      IF (.NOT. lrejhor) THEN
        isame  = 1
        isamec = 1
        itimee = NINT( zsghed(iairls(isortsl(irep),1,3),nhtime)*c60 )
        DO itmp = irep + 1 , ntotssl
          IF (ABS( zsghed(iairls(isortsl(irep),1,3),nhilon)                    &
                  -zsghed(iairls(isortsl(itmp),1,3),nhilon))                   &
             +ABS( zsghed(iairls(isortsl(irep),1,3),nhjlat)                    &
                  -zsghed(iairls(isortsl(itmp),1,3),nhjlat)) < .005_wp) THEN
            isame = isame + 1
            IF ( NINT( zsghed(iairls(isortsl(itmp),1,3),nhtime)*c60 )          &
                -itimee > 12) THEN
              isamec = isamec + 1
              itimee = NINT( zsghed(iairls(isortsl(itmp),1,3),nhtime)*c60 )
            ENDIF
          ENDIF
        ENDDO
        IF ((isamec >= 3) .AND. (2*isame > ntotssl))  lrejhor = .TRUE.
      ENDIF
    ENDDO

    IF (lrejhor) PRINT '("aircraft ",A ," rejected:",I3," out of",I4," report" &
                       &,"s horizontally collocated")' , ystid, isame, ntotssl
    IF (lrejhor) PRINT '("aircraft ",A ,"         :",I3," collocated reports"  &
                       &," with dt > 12 min")' ,         ystid, isamec

! vertically, only below 500 hPa: reject all, if colocated reports (3 or more)
! ----------       comprise at least 50 % (below 700 hPa) or 75 % of all reports

    lrejver = .FALSE.
    DO irep = 1 , ntotssl-2
      IF ((izpp(irep) > 50000) .AND. (.NOT. lrejver)) THEN
        isamec = 1
        DO itmp = irep + 1 , ntotssl
          IF (ABS( izpp(irep)-izpp(itmp) ) == 0)  isamec = isamec + 1
        ENDDO
        IF (      (isamec >= 3)                                                &
            .AND. (     ((2*isamec >= ntotssl) .AND. (izpp(irep) > 70000))     &
                   .OR. ( 4*isamec >= 3*ntotssl)))   lrejver = .TRUE.
      ENDIF
    ENDDO

    IF (lrejver) PRINT '("aircraft ",A ," rejected:",I3," out of",I4," report" &
                       &,"s vertically collocated")' , ystid, isame, ntotssl
    IF (lrejver) PRINT '("aircraft ",A ,"         :",I3," collocated reports"  &
                       &," with dt > 12 min")' ,       ystid, isamec

! set the passive bit for the used single-level report in list 'asl'
! so that it can be set in the ODR later on
! ------------------------------------------------------------------

    lrejrep  =  (lrejhor) .OR. (lrejver)
    IF (lrejrep) THEN
      DO irep = 1 , ntotssl
        iarlls (mzslon(iairls(isortsl(irep),1,3)),2)  =  - 2
        iairls (isortsl(irep),2,3)  =  0
      ENDDO
      ilen = 6 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt18
        outbuf(nacout+3) = ntotssl
        outbuf(nacout+4) = 2
        IF (lrejhor) outbuf(nacout+4) = 1
        outbuf(nacout+5) = NINT( zsghed(iairls(isortsl(1),1,3),nhtime) * c100 )
        outbuf(nacout+6) = NINT( zsghed(iairls(isortsl(ntotssl),1,3),nhtime)   &
                               * c100 )
        DO icl = 1 , istrej
          outbuf(nacout+6+icl) = ICHAR( yarlls(mzslon(iairls(1,1,3))) (icl:icl))
        ENDDO
        nacout  = nacout + ilen
      ENDIF
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
!  Section 4  : Perform the flight track check
!-------------------------------------------------------------------------------

! initialisation
  IF ((lchktrc) .AND. (.NOT. lrejrep) .AND. (.NOT. lmodes)) THEN
    DO irep = 1 , ntotssl
      zflat  (irep) = COS( zsghed(iairls(isortsl(irep),1,3),nhjlat) *r_degrad)
      iflftc (irep) = iarlls(mzslon(iairls(isortsl(irep),1,3)),2)
    ENDDO
  ENDIF
  lchange = (lchktrc) .AND. (.NOT. lrejrep) .AND. (.NOT. lmodes)
  DO WHILE (lchange)
    lchange = .FALSE.
    DO irep = 1 , ntotssl
      fideh (irep,1) = c100 + c1
      fideh (irep,2) = fideh(irep,1)
      fideh (irep,3) = fideh(irep,1)
      fideh (irep,4) = fideh(irep,1)
      fidev (irep,1) = fideh(irep,1)
      fidev (irep,2) = fideh(irep,1)
    ENDDO

! loop: prepare to do the flight track check in forward and backward direction
    DO iqcl = 1 , 2
      irepcm1 = 0
      irepcm2 = 0
      irepcm3 = 0
      IF (iqcl == 1) THEN
        IF (iflftc(2) >= -1)  irepcm1 = 2
        IF (iflftc(1) >= -1)  irepcm2 = 1
        ireprb  = 3
        irepre  = ntotssl
      ELSE
        IF (iflftc(ntotssl-1) >= -1)  irepcm1 = ntotssl - 1
        IF (iflftc(ntotssl  ) >= -1)  irepcm2 = ntotssl
        ireprb  = ntotssl - 2
        irepre  = 1
      ENDIF
      IF (irepcm1 == 0)  irepcm1 = irepcm2
      IF (irepcm2 == irepcm1)  irepcm2 = 0

! loop over the reports (in foreward or backward direction)
      DO irep = ireprb , irepre , 3 - 2*iqcl
        nob   = iairls(isortsl(irep),1,3)
        zlat  = zflat(irep) * zflat(irep)
        fidexy = c0
        fidez  = c0
        fidetw = c0
        fidex2 = c0
        fidet2 = c0
        nest   = 4
!       IF (iarlls(mzslon(nob),2) == -2)  nest = 0

! loop over pairs of previous reports used to produce estimates
        DO itmp = 1 , nest

! determine the pair of previous reports
          nobm1 = 0
          nobm2 = 0
          IF (itmp == 1) THEN
            IF (irepcm1 > 0)  nobm1 = iairls(isortsl(MAX(irepcm1,1)),1,3)
            IF (irepcm2 > 0)  nobm2 = iairls(isortsl(MAX(irepcm2,1)),1,3)
          ELSEIF (itmp == 2) THEN
            IF (irepcm1 > 0)  nobm1 = iairls(isortsl(MAX(irepcm1,1)),1,3)
            IF (irepcm3 > 0)  nobm2 = iairls(isortsl(MAX(irepcm3,1)),1,3)
          ELSEIF (itmp == 3) THEN
            IF (irepcm2 > 0)  nobm1 = iairls(isortsl(MAX(irepcm2,1)),1,3)
            IF (irepcm3 > 0)  nobm2 = iairls(isortsl(MAX(irepcm3,1)),1,3)
          ELSEIF (itmp == 4) THEN
            irepm1 = irep + SIGN( 1 , 2*iqcl-3 )
            irepm2 = irep + SIGN( 2 , 2*iqcl-3 )
            nobm1 = iairls(isortsl(irepm1),1,3)
            nobm2 = iairls(isortsl(irepm2),1,3)
            IF (iflftc(irepm1) == -9)  nobm1 = 0
            IF (iflftc(irepm2) == -9)  nobm2 = 0
          ENDIF

! compute the horizontal and vertical confidences
! -----------------------------------------------

          IF ((nobm1 > 0) .AND. (nobm2 > 0)) THEN
! degree of the extrapolation
            fexpol =  ABS( zsghed(nobm1,nhtime) - zsghed(nobm2,nhtime) )
            IF (fexpol > epsy) THEN
!   regular case
              fexpol = ABS( zsghed(nob,nhtime) -zsghed(nobm2,nhtime) ) / fexpol
            ELSE
!   if predictors have same time, then add 1 minute to numerator and denominator
!   (--> if predictors and predicand have same time, the estimated position
!        of the predicand is set equal to the position of the 2nd predictor)
              fexpol = ABS( zsghed(nob,nhtime) -zsghed(nobm2,nhtime) ) *c60 + c1
            ENDIF
! compute estimates of the (horizontal / vertical) position for current report
            xestim = (c1-fexpol)* zsghed(nobm2,nhilon)                         &
                    +    fexpol * zsghed(nobm1,nhilon)
            yestim = (c1-fexpol)* zsghed(nobm2,nhjlat)                         &
                    +    fexpol * zsghed(nobm1,nhjlat)
!!    --> set 'fexpol = 1 + scal*(fexpol-1)' and replace 'scal=1' by 'scal=0.5'
!           zestim =  c05 *(c1-fexpol)* zsgbdy(nobm2,nbsp)                     &
!                   + c05 *(c1+fexpol)* zsgbdy(nobm1,nbsp)
            zestim = (c1-fexpol)* zsgbdy(nobm2,nbsp)                           &
                    +    fexpol * zsgbdy(nobm1,nbsp)
            zestim = MAX( 10000._wp , MIN( zestim , 105000._wp ) )
! compute the distance to the estimate
            estdxy = SQRT( (zlat* (zsghed(nob  ,nhilon) - xestim)              &
                                * (zsghed(nob  ,nhilon) - xestim)              &
                             +    (zsghed(nob  ,nhjlat) - yestim)              &
                                * (zsghed(nob  ,nhjlat) - yestim)) *rdegkm2 )
            estdz  = ABS( zsgbdy(nob  ,nbsp) - zestim )
!   test the observation location also with negative sign for the longitude
            estdx2 = SQRT( (zlat* (-ABS( zsghed(nob,nhilon) ) - xestim)        &
                                * (-ABS( zsghed(nob,nhilon) ) - xestim)        &
                             +    (zsghed(nob  ,nhjlat) - yestim)              &
                                * (zsghed(nob  ,nhjlat) - yestim)) *rdegkm2 )
! prepare an extension of the reference distance,
! if the extrapolation is not well conditioned
            exrefd = REAL (NINT(ABS( zsghed(nobm1,nhtime) -zsghed(nobm2,nhtime)) *c60), wp)
            exrefd = MAX( fexpol - c1 , c3 / (exrefd+c1) , c0 )
! compute the reference distance
            refdxy = SQRT( (zlat* (zsghed(nobm1,nhilon) - xestim)              &
                                * (zsghed(nobm1,nhilon) - xestim)              &
                             +    (zsghed(nobm1,nhjlat) - yestim)              &
                                * (zsghed(nobm1,nhjlat) - yestim)) *rdegkm2 )
!     (The last term allows increasing distance to the estimated location at
!      a speed of 750 km/h.)
            refdxy =   refdxy   +   c30   +      c5* exrefd   +   c12p5 *      &
                     (ABS( zsghed(nob,nhtime) -zsghed(nobm1,nhtime) ) *c60 + c1)
!     (The distance to the pressure estimate is not included in 'refdz' since
!      the vertical estimate error may be as large for an aircraft going from
!      level flight to descent as e.g. for the opposite case.
!      The last term allows (changes of) descent (or ascent) rates of 50 hPa
!      per minute , e.g. after the end of the level flight phase.)
!           refdz  = ABS( zsgbdy(nobm1,nbsp) - zestim )
            refdz  =              c5000   +   c2500* exrefd   +   c5000 *      &
                     (ABS( zsghed(nob,nhtime) -zsghed(nobm1,nhtime) ) *c60 + c1)
!        (at flight ascents, allow a 100 hPa decrease within the first minute)
            IF (mphase(irep) == 5)  refdz = refdz + c5000
! compute the confidence
!           fidexy = MAX( fidexy , c90 - estdxy /refdxy *c40 , c1 )
!           fidez  = MAX( fidez  , c90 - estdz  /refdz  *c30 , c1 )
            fidew  = c05 *(  ABS( zsghed(nob,nhtime)-zsghed(nobm1,nhtime) )    &
                           + ABS( zsghed(nob,nhtime)-zsghed(nobm2,nhtime) ))
!   the weight for the confidence from the current estimate is reduced, if
!   - the average of the predictors is > 0.5 hrs older than the predictand, or
!   - the extrapolation is not well conditioned,
            fidew  = MIN( c1 , c05 / MAX( fidew , epsy )                       &
                             , c1 / (c1 + c05*(exrefd-c1)) )
            lupdate = .TRUE.
            lupdat2 = .TRUE.
            IF ((itmp == 4) .AND. (fidetw > epsy)) THEN
!   update the total confidence by the current confidence deduced from the 2
!   immediately previous (and possibly erroneous) reports only if this increases
!   the minimum of the 2 (i.e. horizontal and vertical) confidences
              fidxpr  = MAX( c90 - estdxy /refdxy *c40 ,c0 )
              fidzpr  = MAX( c90 - estdz  /refdz  *c30 ,c0 )
              lupdate = (  MIN( fidexy /fidetw , fidez /fidetw )               &
                         < MIN( (fidexy +fidew*fidxpr) /(fidetw +fidew)        &
                              , (fidez  +fidew*fidzpr) /(fidetw +fidew) ) )
              fidxp2  = MAX( c90 - estdx2 /refdxy *c40 ,c0 )
              lupdat2 = (fidex2/fidet2 < (fidex2+fidew*fidxp2)/(fidet2+fidew))
            ENDIF
            IF (lupdate) THEN
              fidexy = fidexy + fidew* MAX( c90 - estdxy/refdxy*c40 ,c0 )
              fidez  = fidez  + fidew* MAX( c90 - estdz /refdz *c30 ,c0 )
              fidetw = fidetw + fidew
            ENDIF
            IF (lupdat2) THEN
              fidex2 = fidex2 + fidew* MAX( c90 - estdx2/refdxy*c40 ,c0 )
              fidet2 = fidet2 + fidew
            ENDIF
!           IF (      (     (estdxy/refdxy > 0.6_wp)                           &
!                      .OR. (estdz/refdz  > 0.75_wp))                          &
!               .AND. (acthr <= 0.01_wp)) THEN
!             fidxpr = MAX( c90 - estdxy /refdxy *c40 ,c0 )
!             fidzpr = MAX( c90 - estdz  /refdz  *c30 ,c0 )
!             fidxwp = fidexy / fidetw
!             fidzwp = fidez  / fidetw
!             PRINT '(''air flt. estm.: '',A ,2F6.2,F8.0,2F6.1,2F8.0,F5.2      &
!                   &,5F5.1)'    , ystid, acthr, zsghed(nob,nhtime)            &
!                   , zsgbdy(nob,nbsp), estdxy, refdxy, zestim, refdz, exrefd  &
!                   , fidxwp, fidzwp, fidetw, fidxpr, fidzpr
!           ENDIF
          ENDIF
        ENDDO
        IF (fidetw > epsy)  fideh (irep,iqcl)   = fidexy / fidetw
        IF (fidetw > epsy)  fidev (irep,iqcl)   = fidez  / fidetw
        IF (fidet2 > epsy)  fideh (irep,iqcl+2) = fidex2 / fidet2
!!      IF (      (lfirst) .AND. (fidetw > epsy)                               &
!       IF (      (acthr <= 0.01_wp) .AND. (fidetw > epsy)                     &
!           .AND. (MIN( fideh(irep,iqcl),fidev(irep,iqcl) ) < c60))            &
!         PRINT '(''air flt. trac : '',A ,2F6.2,F8.0,I2,3F5.1)'                &
!               , ystid, acthr, zsghed(nob,nhtime), zsgbdy(nob,nbsp)           &
!               , iqcl, fideh(irep,iqcl), fidev(irep,iqcl), fideh(irep,iqcl+2)
        IF ((     (MIN( fideh(irep,iqcl),fidev(irep,iqcl) ) >= c60)            &
             .OR. (fidetw <= epsy)) .AND. (nest > 0)) THEN
          irepcm3 = irepcm2
          irepcm2 = irepcm1
          IF (iflftc(irep) >= -1)  irepcm1 = irep
        ENDIF
      ENDDO
    ENDDO

! 1.: check iteratively for missing sign at reported longitude:
!     set reports with confidence < 65 % to passive
!     if it were > 80 % with negative sign for longitude
! -------------------------------------------------------------

    DO irep = 1 , ntotssl
      IF (      (     ((fideh(irep,1) < c50) .AND. (fideh(irep,3) > c80))      &
                 .OR. ((fideh(irep,2) < c50) .AND. (fideh(irep,4) > c80)))     &
          .AND. (iflftc(irep) >= -1)) THEN
        lchange = .TRUE.
        iflftc (irep)  =  - 9
!!      IF (lfirst) THEN
!       IF (acthr <= 0.01_wp) THEN
!         nob   = iairls(isortsl(irep),1,3)
!         PRINT '(''air flt. tracs: '',A ,2F6.2,F8.0,4F5.1)' , ystid, acthr    &
!             , zsghed(nob,nhtime), zsgbdy(nob,nbsp), (fideh(irep,icl), icl=1,4)
!       ENDIF
      ENDIF
    ENDDO
    IF (.NOT. lchange) THEN
      DO irep = 1 , ntotssl
        IF (      (     ((fideh(irep,1) < c65) .AND. (fideh(irep,3) > c80))    &
                   .OR. ((fideh(irep,2) < c65) .AND. (fideh(irep,4) > c80)))   &
            .AND. (iflftc(irep) >= -1)) THEN
          lchange = .TRUE.
          iflftc (irep)  =  - 9
        ENDIF
      ENDDO
    ENDIF
    IF (.NOT. lchange) THEN
      DO irep = 1 , ntotssl
        mprout = 0
        nob   = iairls(isortsl(irep),1,3)
        IF (      (iflftc(irep) == -9)                                         &
            .AND. (     (fideh(irep,1) < fideh(irep,3))                        &
                   .OR. (fideh(irep,2) < fideh(irep,4)))) THEN
          IF (iarlls(mzslon(nob),2) >= 1)  mprout = 1
          iarlls (mzslon(nob),2)  =  - 2

! 2.: check (non-iteratively) for other location errors: Set reports to passive,
!     - if the foreward or backward confidence is < 60 % and the other is < 85 %
!       (or the other
!     - or if both the foreward and backward confidences are < 65 %
! ------------------------------------------------------------------------------

        ELSEIF ((fideh(irep,1) <= c100) .AND. (fideh(irep,2) <= c100)) THEN
          IF (     (      (fideh(irep,1) < c60)                                &
                    .AND. (fideh(irep,2) < c85+(c60-fideh(irep,1))*c10r))      &
              .OR. (      (fideh(irep,2) < c60)                                &
                    .AND. (fideh(irep,1) < c85+(c60-fideh(irep,2))*c10r))      &
              .OR. (      (fidev(irep,1) < c60)                                &
                    .AND. (fidev(irep,2) < c85+(c60-fidev(irep,1))*c10r))      &
              .OR. (      (fidev(irep,2) < c60)                                &
                    .AND. (fidev(irep,1) < c85+(c60-fidev(irep,2))*c10r))      &
              .OR. ((fideh(irep,1) < c65) .AND. (fideh(irep,2) < c65))         &
              .OR. ((fidev(irep,1) < c65) .AND. (fidev(irep,2) < c65))) THEN
            IF (iarlls(mzslon(nob),2) >= 1)  mprout = 2
            iarlls (mzslon(nob),2)  =  - 2
!!          IF (lfirst)                                                        &
!           IF (acthr <= 0.01_wp)                                              &
!             PRINT '(''air flt. tracs: '',A ,2F6.2,F8.0,4F5.1)'               &
!                   , ystid, acthr, zsghed(nob,nhtime), zsgbdy(nob,nbsp)       &
!                   , (fideh(irep,icl),icl=1,2), (fidev(irep,icl),icl=1,2)
          ENDIF
        ELSEIF ((fideh(irep,1) < c60) .OR. (fidev(irep,1) < c60)) THEN
          IF (iarlls(mzslon(nob),2) >= 1)  mprout = 3
          iarlls (mzslon(nob),2)  =  - 2
        ELSEIF ((fideh(irep,2) < c60) .OR. (fidev(irep,2) < c60)) THEN
          IF (iarlls(mzslon(nob),2) >= 1)  mprout = 4
          iarlls (mzslon(nob),2)  =  - 2
        ENDIF
        IF (iarlls(mzslon(nob),2) <= -2)  iairls (isortsl(irep),2,3)  =  0
        IF (mprout > 0) THEN
          ilen = 7 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt19
            outbuf(nacout+3) = NINT( zsgbdy(nob,nbsp) )
            outbuf(nacout+4) = NINT( zsghed(nob,nhtime) * c100 )
            IF ((mprout == 1) .AND. (fideh(irep,1) < fideh(irep,2))) THEN
              outbuf(nacout+5) = 1
              outbuf(nacout+6) = NINT( fideh(irep,1) * c100 )
              outbuf(nacout+7) = NINT( fideh(irep,3) * c100 )
            ELSEIF (mprout == 1) THEN
              outbuf(nacout+5) = 2
              outbuf(nacout+6) = NINT( fideh(irep,2) * c100 )
              outbuf(nacout+7) = NINT( fideh(irep,4) * c100 )
            ELSEIF (  MIN( fideh(irep,1),fideh(irep,2) )                       &
                    < MIN( fidev(irep,1),fidev(irep,2) )) THEN
              outbuf(nacout+5) = 3
              outbuf(nacout+6) = NINT( fideh(irep,1) * c100 )
              outbuf(nacout+7) = NINT( fideh(irep,2) * c100 )
            ELSE
              outbuf(nacout+5) = 4
              outbuf(nacout+6) = NINT( fidev(irep,1) * c100 )
              outbuf(nacout+7) = NINT( fidev(irep,2) * c100 )
            ENDIF
            DO icl = 1 , istrej
              outbuf(nacout+7+icl) = ICHAR( yarlls(mzslon(nob)) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO


!-------------------------------------------------------------------------------
!  Section 5  : Prepare the production of multi-level reports by finding
!               the lowest 'active candidate' (report with highest pressure),
!               and put it on a list
!-------------------------------------------------------------------------------

! preset counters of used reports, and of available active reports
! ----------------------------------------------------------------

  nactssl = 0
  DO irep = 1 , ntotssl
    IF (iairls(isortsl(irep),2,3) >= 1)  nactssl = nactssl + 1
  ENDDO

  nusedssl = 0

! loop while not all active reports with the present station id are used
! ----------------------------------------------------------------------

DO WHILE ((nusedssl < nactssl) .AND. (thairh > epsy) .AND. (.NOT. lmodes))

! find a lowest 'active candidate'
! --------------------------------

  pslmax   = 2000.0_wp
  timob1   = -100.0_wp
  lsamealt = .FALSE.
  DO irep = 1 , ntotssl
    nssl   = isortsl (irep)
    itsl   = iairls (nssl,1,3)
    nlon   = mzslon (itsl)
    IF ((iairls(nssl,2,3) >= 1) .AND. (zsgbdy(itsl,nbsp) >= pslmax-epsy)) THEN
      lsamealt =       (zsgbdy(itsl,nbsp) <= pslmax+epsy)                      &
                 .AND. (ABS(timob1-zsghed(itsl,nhtime)) <= thairt)
      iasmls (1,1) = itsl
      iasmls (1,2) = nlon
      iasmls (1,3) = nssl
      pslmax   = zsgbdy(itsl,nbsp)
      timob1   = zsghed(itsl,nhtime)
    ENDIF
  ENDDO

! if there are several neighbouring lowest 'active candidates':
! find the most appropriate one: temporally or horizontally the
! closest one to the next report further above
! -------------------------------------------------------------
! --> a: find the next active report further above
!        (if the are several, take the one according to the
!         following criteria: 1: the earliest
!                             2: the westernmost
!                             3: the southernmost
!         this selection must be unique (if reproducible), and
!         in fact is unique, considering that reports assigned
!         to the same grid point are always ordered identically)

  IF (lsamealt) THEN
    psnext = 2000.0_wp
    tmnext = 1000.0_wp
    ienext = ABS( imdi )
    jenext = ABS( imdi )
! if no 'next' report is found, pre-set 'nsnext' so that the lowest report
! found finally below is the one for which no 'next' report exists
    nsnext = iasmls(1,1)
    DO irep = 1 , ntotssl
      nssl   = isortsl (irep)
      itsl   = iairls (nssl,1,3)
      IF (      (iairls(nssl,2,3) >= 1)                                        &
          .AND. (zsgbdy(itsl,nbsp) >= psnext-epsy)                             &
          .AND. (zsgbdy(itsl,nbsp) <  pslmax-epsy)                             &
          .AND. (ABS( zsghed(itsl,nhtime)-timob1 ) <= thairt)) THEN
        lobnext = .FALSE.
        IF (zsgbdy(itsl,nbsp) > psnext+epsy) THEN
          lobnext = .TRUE.
        ELSEIF (zsghed(itsl,nhtime) <  tmnext-epsy) THEN
          lobnext = .TRUE.
        ELSEIF (zsghed(itsl,nhtime) <= tmnext+epsy) THEN
          IF (mzsghd(itsl,nhitot) < ienext) THEN
            lobnext = .TRUE.
          ELSEIF (      (mzsghd(itsl,nhitot) == ienext)                        &
                  .AND. (mzsghd(itsl,nhjtot) <  jenext)) THEN
            lobnext = .TRUE.
          ENDIF
        ENDIF
        IF (lobnext) THEN
          psnext = zsgbdy (itsl,nbsp)
          tmnext = zsghed (itsl,nhtime)
          ienext = mzsghd (itsl,nhitot)
          jenext = mzsghd (itsl,nhjtot)
          nsnext = itsl
        ENDIF
      ENDIF
    ENDDO
    IF (lwonl)                                                                 &
      WRITE( nupr,'(''AMBIGUOUS base report; next report further above at''    &
                  &, F8.0, F6.2, 2I4, 2X, A )' )                               &
             psnext, tmnext, ienext, jenext, ystid
!   PRINT         '(''AMBIGUOUS base report; next report further above at''    &
!                 &, F8.0, F6.2, 2I4, 2X, A )' ,                               &
!            psnext, tmnext, ienext, jenext, ystid

! --> b: get the closest report (unique selection)

    tmdist = 1000.0_wp
    r2dist = 20000.0_wp **2
    ienext = ABS( imdi )
    jenext = ABS( imdi )
    obtimeb = zsghed (nsnext,nhtime)
    obilonb = zsghed (nsnext,nhilon)
    objlatb = zsghed (nsnext,nhjlat)
    cosolat = COS( objlatb * r_degrad )
    DO irep = 1 , ntotssl
      nssl   = isortsl (irep)
      itsl   = iairls (nssl,1,3)
      nlon   = mzslon (itsl)
      IF ((iairls(nssl,2,3) >= 1) .AND. (zsgbdy(itsl,nbsp) > psnext+epsy)) THEN
        lobnext = .FALSE.
        rhzob   = (  (cosolat *(zsghed(itsl,nhilon) - obilonb))**2             &
                   +           (zsghed(itsl,nhjlat) - objlatb) **2 ) * rdegkm2
        IF (ABS( zsghed(itsl,nhtime)-obtimeb ) < tmdist-epsy) THEN
          lobnext = .TRUE.
        ELSEIF (ABS( zsghed(itsl,nhtime)-obtimeb ) <= tmdist+epsy) THEN
          IF (rhzob < r2dist-epsy) THEN
            lobnext = .TRUE.
          ELSEIF (rhzob <= r2dist+epsy) THEN
            IF (mzsghd(itsl,nhitot) < ienext) THEN
              lobnext = .TRUE.
            ELSEIF (      (mzsghd(itsl,nhitot) == ienext)                      &
                    .AND. (mzsghd(itsl,nhjtot) <  jenext)) THEN
              lobnext = .TRUE.
            ENDIF
          ENDIF
        ENDIF
        IF (lobnext) THEN
          iasmls (1,1) = itsl
          iasmls (1,2) = nlon
          iasmls (1,3) = nssl
          tmdist = ABS( zsghed(itsl,nhtime) - obtimeb )
          r2dist = rhzob
          ienext = mzsghd (itsl,nhitot)
          jenext = mzsghd (itsl,nhjtot)
        ENDIF
      ENDIF
    ENDDO
    IF (lwonl)                                                                 &
      WRITE( nupr,'(''AMBIGUOUS lowest remaining single-lv. air report at''    &
                  &, F8.0, F6.2, 2I4, 2X, A )' )                               &
             pslmax, zsghed(iasmls(1,1),nhtime), ienext, jenext, ystid
  ENDIF

! get 4-d location, and cos(latitude) for distances
! -------------------------------------------------

  itsl  = iasmls (1,1)
  obtimeb = zsghed (itsl,nhtime)
  obilonb = zsghed (itsl,nhilon)
  objlatb = zsghed (itsl,nhjlat)
  obbaltb = zsgbdy (itsl,nbsp  )
  cosolat = COS( objlatb * r_degrad )
! PRINT '(''bottom aircraft report'', 2I5,I4, F8.0, F7.3, 2F7.2)'              &
!       , iasmls(1,1), iasmls(1,3), nusedssl, obbaltb, obtimeb, obilonb, objlatb

!-------------------------------------------------------------------------------
!  Section 6  : Find all single-level reports which will be used to produce a
!               multi-level report with the report found in section 2 as basis
!-------------------------------------------------------------------------------

  nlev     = 1
  lastfind = .TRUE.
  lsamealt = .FALSE.
  obbaltl  = obbaltb

! loop while further reports may be added to list from which the multi-level
! report is produced
! criteria: 1. the last trial for adding one more report has been successful
!           2. there is no other report in the horizontal and temporal
!              neighbourhood at the same altitude as the latest report added
!           3. the number of reports does not exceed the number of levels allowd
! ------------------------------------------------------------------------------

  DO WHILE ((lastfind) .AND. (.NOT. lsamealt) .AND. (nlev < maxarl))

    pdiffmn  = thairv
    lastfind = .FALSE.
    tmdist = 1000.0_wp
    r2dist = 20000.0_wp **2
    ienext = ABS( imdi )
    jenext = ABS( imdi )

! loop to find lowest 'active candidate' which
! can be added on top of the last report found
! --------------------------------------------

    DO irep = 1 , ntotssl

! get all required information on current report
      nssl    = isortsl (irep)
      itsl    = iairls (nssl,1,3)
      nlon    = mzslon (itsl)
      lactive = iairls (nssl,2,3)   >=   1
      obtimec = zsghed (itsl,nhtime)
      obilonc = zsghed (itsl,nhilon)
      objlatc = zsghed (itsl,nhjlat)
      obbaltc = zsgbdy (itsl,nbsp  )
!     IF ((nlev == 1) .AND. (nusedssl == 0))                                   &
!       PRINT '(''cand '', I5, 2F8.1)' , nlon, obbaltl, obbaltc

! check if the current report meets the requirements to be added on the list
! better than all preceedingly checked reports:   the report must
! - be active
! - be above the last accepted report but below all preceedingly checked reports
! - be within 'thairt' hrs of the lowermost report
! - be within 'thairh' km  of the lowermost report
!   (not required: be within 'rhzlair' km of the last accepted report)
! ------------------------------------------------------------------------------

      IF (lactive) THEN
        IF (      (obbaltc < obbaltl-epsy)                                     &
            .AND. (obbaltc > obbaltl-pdiffmn-epsy)) THEN
          IF (ABS(obtimec-obtimeb) <= thairt) THEN
            rhzob = (  (cosolat *(obilonc - obilonb))**2                       &
                     +           (objlatc - objlatb) **2 ) * rdegkm2
            IF (rhzob <= thairh2) THEN

! if all checks passed: list current report (or replace previously found report)
              lsamealt = (obbaltc  <  obbaltl-pdiffmn+epsy) .AND. (lastfind)
              IF (.NOT. lastfind) nlev = nlev + 1
              IF (.NOT. lsamealt) THEN
                lobnext = .TRUE.
              ELSEIF (ABS( obtimec-obtimeb ) < tmdist-epsy) THEN
                lobnext = .TRUE.
              ELSEIF (ABS( obtimec-obtimeb ) <= tmdist+epsy) THEN
                IF (rhzob < r2dist-epsy) THEN
                  lobnext = .TRUE.
                ELSEIF (rhzob <= r2dist+epsy) THEN
                  IF (mzsghd(itsl,nhitot) < ienext) THEN
                    lobnext = .TRUE.
                  ELSEIF (      (mzsghd(itsl,nhitot) == ienext)                &
                          .AND. (mzsghd(itsl,nhjtot) <  jenext)) THEN
                    lobnext = .TRUE.
                  ENDIF
                ENDIF
              ENDIF
              IF (lobnext) THEN
                iasmls (nlev,1) = itsl
                iasmls (nlev,2) = nlon
                iasmls (nlev,3) = nssl
                tmdist = ABS( zsghed(itsl,nhtime) - obtimeb )
                r2dist = rhzob
                ienext = mzsghd (itsl,nhitot)
                jenext = mzsghd (itsl,nhjtot)
                pdiffmn  = obbaltl - obbaltc
                lastfind = .TRUE.
              ENDIF
!             PRINT '(''cand all checks ok '',I4,I5,I4, F8.1,F7.3,2F7.2)'      &
!                   , nssl, nlon, nlev, obbaltc, obtimec, obilonc, objlatc
            ENDIF
          ENDIF
        ENDIF
      ENDIF

! close loop to find next 'lowest active candidate',
! if search has been successful, store its pressure and update counter for
!                                probably ascending or descending flight'
! and close inner while loop (for finding elements for one multi-level report)
! ----------------------------------------------------------------------------

    ENDDO
    IF (lastfind) obbaltl = zsgbdy (iasmls(nlev,1),nbsp)
  ENDDO

!-------------------------------------------------------------------------------
!  Section 7  : Produce the multi-level report (if the list contains at least
!               'minslml' reports),
!               and set used single-level reports to passive on list 'asl'
!-------------------------------------------------------------------------------

  IF ((nlev >= minslml) .AND. (ntotaml < maxmll)) THEN

    itsl    = iasmls (1,1)
    ntotaml = ntotaml + 1

! copy header of bottom single-level report onto header of multi-level report
! and complement it
! ---------------------------------------------------------------------------

    zmlhed (ntotaml,nhilon) = zsghed (itsl,nhilon)
    zmlhed (ntotaml,nhjlat) = zsghed (itsl,nhjlat)
    zmlhed (ntotaml,nhalt ) = zsghed (itsl,nhalt )
    zmlhed (ntotaml,nhtime) = zsghed (itsl,nhtime)
    zmlhed (ntotaml,nhsurf) = zsghed (itsl,nhsurf)
    zmlhed (ntotaml,nhzio ) = zsghed (itsl,nhzio )
    zmlhed (ntotaml,nhzjo ) = zsghed (itsl,nhzjo )
    zmlhed (ntotaml,nhsynt) = zsghed (itsl,nhsynt)
    zmlhed (ntotaml,nhtddb) = zsghed (itsl,nhtddb)
    zmlhed (ntotaml,nhsolz) = zsghed (itsl,nhsolz)
    zmlhed (ntotaml,nhssos) = zsghed (itsl,nhssos)
    mzmlhd (ntotaml,nhio  ) = mzsghd (itsl,nhio  )
    mzmlhd (ntotaml,nhjo  ) = mzsghd (itsl,nhjo  )
    mzmlhd (ntotaml,nhitot) = mzsghd (itsl,nhitot)
    mzmlhd (ntotaml,nhjtot) = mzsghd (itsl,nhjtot)
    mzmlhd (ntotaml,nhobtp) = mzsghd (itsl,nhobtp)
    mzmlhd (ntotaml,nhcode) = mzsghd (itsl,nhcode)
    mzmlhd (ntotaml,nhschr) = mzsghd (itsl,nhschr)
    mzmlhd (ntotaml,nhpass) = 0
!   mzmlhd (ntotaml,nhqcfw) = mzsghd (itsl,nhqcfw)
    mzmlhd (ntotaml,nhqcfw) = -1
    mzmlhd (ntotaml,nhflag) = mzsghd (itsl,nhflag)
    mzmlhd (ntotaml,nhflag) = IBCLR( mzmlhd(ntotaml,nhflag), FL_MERGE )
    mzmlhd (ntotaml,nhcorr) = mzsghd (itsl,nhcorr)
    mzmlhd (ntotaml,nhcat ) = mzsghd (itsl,nhcat )
    mzmlhd (ntotaml,nhcats) = mzsghd (itsl,nhcats)
    mzmlhd (ntotaml,nhkz  ) = mzsghd (itsl,nhkz  )
    mzmlhd (ntotaml,nhcent) = mzsghd (itsl,nhcent)
    mzmlhd (ntotaml,nhstid) = mzsghd (itsl,nhstid)
    mzmlhd (ntotaml,nhdate) = mzsghd (itsl,nhdate)
    mzmlhd (ntotaml,nhhrmn) = mzsghd (itsl,nhhrmn)
    mzmlhd (ntotaml,nhsyhr) = mzsghd (itsl,nhsyhr)
    mzmlhd (ntotaml,nhrtyp) = imdi
    mzmlhd (ntotaml,nhnlev) = nlev
    mzmlhd (ntotaml,nhaexi) = 0
    mzmlhd (ntotaml,nhuexi) = 0
    mzmlhd (ntotaml,nhtexi) = 0
    mzmlhd (ntotaml,nhqexi) = 0
    mzmlhd (ntotaml,nhtrac) = imdi
    mzmlhd (ntotaml,nhrad ) = imdi
    mzmlhd (ntotaml,nhna4 ) = imdi
    mzmlhd (ntotaml,nhwce ) = imdi
    mzmlhd (ntotaml,nhdt  ) = imdi
    zmlhed (ntotaml,nhvcbu) = c1
    zmlhed (ntotaml,nhvcbt) = c1
    zmlhed (ntotaml,nhvcbq) = c1
    zmlhed (ntotaml,nhvctu) = c1
    zmlhed (ntotaml,nhvctt) = c1
    zmlhed (ntotaml,nhvctq) = c1
    zmlhed (ntotaml,nhtvip) = rmdi
    zmlhed (ntotaml,nhtviz) = rmdi

! copy bodies of single-level reports onto body of multi-level report,
! complement this body,
! and list the list-'asl' indices of the single-level reports used
! --------------------------------------------------------------------

    DO klev = 1 , nlev
      itsl  = iasmls (klev,1)
      zmlbdy (ntotaml,klev,nbtu  ) = zsgbdy (itsl,nbsu  )
      zmlbdy (ntotaml,klev,nbtv  ) = zsgbdy (itsl,nbsv  )
      zmlbdy (ntotaml,klev,nbtt  ) = zsgbdy (itsl,nbst  )
      zmlbdy (ntotaml,klev,nbtrh ) = zsgbdy (itsl,nbsrh )
      zmlbdy (ntotaml,klev,nbtp  ) = zsgbdy (itsl,nbsp  )
      zmlbdy (ntotaml,klev,nbtz  ) = zsgbdy (itsl,nbsz  )
      zmlbdy (ntotaml,klev,nbtuer) = zsgbdy (itsl,nbsuer)
      zmlbdy (ntotaml,klev,nbtter) = zsgbdy (itsl,nbster)
      zmlbdy (ntotaml,klev,nbtqer) = zsgbdy (itsl,nbsqer)
      zmlbdy (ntotaml,klev,nbtzer) = zsgbdy (itsl,nbszer)
      zmlbdy (ntotaml,klev,nbtzio) = zsghed (itsl,nhilon)
      zmlbdy (ntotaml,klev,nbtzjo) = zsghed (itsl,nhjlat)
!     zmlbdy (ntotaml,klev,nbtzio) = zsghed (itsl,nhzio ) - mzsghd(itsl,nhio)  &
!                                                         + mzsghd(itsl,nhitot)
!     zmlbdy (ntotaml,klev,nbtzjo) = zsghed (itsl,nhzjo ) - mzsghd(itsl,nhjo)  &
!                                                         + mzsghd(itsl,nhjtot)
      zmlbdy (ntotaml,klev,nbttim) = zsghed (itsl,nhtime) - zmlhed(ntotaml     &
                                                                       ,nhtime)
      zmlbdy (ntotaml,klev,nbtlop) = LOG( zsgbdy(itsl,nbsp) )
      zmlbdy (ntotaml,klev,nbtplv) = zsgbdy (itsl,nbsplv)
      zmlbdy (ntotaml,klev,nbtdrh) = zsgbdy (itsl,nbsdrh)
      zmlbdy (ntotaml,klev,nbtw  ) = rmdi
      zmlbdy (ntotaml,klev,nbtsnr) = rmdi
      zmlbdy (ntotaml,klev,nbtuac) = rmdi
      mzmlbd (ntotaml,klev,nbtflg) = mzsgbd (itsl,nbsflg)
      mzmlbd (ntotaml,klev,nbterr) = mzsgbd (itsl,nbserr)
      mzmlbd (ntotaml,klev,nbtqcf) = mzsgbd (itsl,nbsqcf)
      mzmlbd (ntotaml,klev,nbtlid) = mzsgbd (itsl,nbslid)
      mzmlbd (ntotaml,klev,nbtlsg) = IBSET ( 0, LS_SIGN )
      IF (BTEST( mzsgbd(itsl,nbserr), nvru )) THEN
        mzmlhd (ntotaml,nhuexi) = 1
        mzmlhd (ntotaml,nhaexi) = 1
      ENDIF
      IF (BTEST( mzsgbd(itsl,nbserr), nvrt )) THEN
        mzmlhd (ntotaml,nhtexi) = 1
        mzmlhd (ntotaml,nhaexi) = 1
      ENDIF
      IF (BTEST( mzsgbd(itsl,nbserr), nvrq )) THEN
        mzmlhd (ntotaml,nhqexi) = 1
        mzmlhd (ntotaml,nhaexi) = 1
      ENDIF

! set the passive bit for the used single-level report in list 'asl'
! so that it can be set in the ODR later on
! ------------------------------------------------------------------

      iarlls (iasmls(klev,2),2)  = 0

! store 'asl' list index of the reports used for the current multi-level report,
! so that - the passive bits in list 'asl' can be unset if the multi-level
!           report has to be cancelled later on due to insufficient array sizes
!         - the station id of the multi-level report can be determined later on
! ------------------------------------------------------------------------------

      mzmlon (ntotaml,klev)      = iasmls(klev,2)
    ENDDO

! Warning and report event, if ODR is too small
! ---------------------------------------------

  ELSEIF (nlev >= minslml) THEN
    kcdtyp = mzsghd (iairls(1,1,3),nhcode)
    kobtyp = mzsghd (iairls(1,1,3),nhobtp)
    icma   = i_cma ( kobtyp , kcdtyp )
    neventr (nenoml,icma) = neventr (nenoml,icma) + 1
    nexcloc = nexcloc + 1
    IF (lwonl)                                                                 &
      WRITE( nupr,'(''CAUTION !!!!! flight number '',A,'' : multi-level ODR''  &
                  &,'' size'',I5,'' too small'')' )   ystid, maxmll
    PRINT         '(''CAUTION !!!!! flight number '',A,'' : multi-level ODR''  &
                  &,'' size'',I5,'' too small'')' ,   ystid, maxmll
  ENDIF

! print message
! -------------

    IF (nlev < minslml) THEN
      IF (lwonl)                                                               &
        WRITE( nupr,'(I2,'' coherent AIREPs:'',I2,'' < '',I1,'' ==> no multi'' &
                    &,''-level report at '',F6.1,F5.1,F8.0,F6.2)' )            &
               nlev, nlev, minslml, obilonb, objlatb, obbaltb, obtimeb
    ELSE
      IF (lwonl)                                                               &
        WRITE( nupr,'(I2,'' coherent AIREPs ==> make multi-level rep.''        &
                    &,I3,'', base at'', F5.1,F5.1, F8.0, F6.2)' )              &
               nlev, ntotaml, obilonb, objlatb, obbaltb, obtimeb
    ENDIF

!-------------------------------------------------------------------------------
!  Section 8  : Update lists and counters, and close loops
!-------------------------------------------------------------------------------

! if (nlev < minslml), account only for report klev=1
  IF (nlev < minslml) nlev = 1

! set (nlev) used single-level reports to passive on list 'ssl'

  DO klev = 1 , nlev
    nssl = iasmls (klev,3)
    iairls (nssl,2,3) = 0
  ENDDO

! update total number of 'used' single-level reports with present station id
  nusedssl = nusedssl + nlev

! close outer while loop
! ----------------------

ENDDO  ! WHILE

!-------------------------------------------------------------------------------
!  Section 9  : Thinning of sequence of single-level reports
!-------------------------------------------------------------------------------

  IF ((lchktrc) .AND. (.NOT. lrejrep) .AND. (.NOT. lmodes)) THEN
    rvtlmt  = rvtlair * c100
!   rvtlmt  = rvtlair * c100 * 5.1_wp
    irepcm1 = 0
    IF (iarlls(mzslon(iairls(isortsl(1),1,3)),2) >= 1)  irepcm1 = 1
    DO irep = 2 , ntotssl
      nob   = iairls(isortsl(irep),1,3)
      IF (iarlls(mzslon(nob),2) >= 1) THEN
        IF (irepcm1 > 0) THEN
          nobm1 = iairls(isortsl(irepcm1),1,3)
          IF (      (NINT( (zsghed(nob,nhtime)-zsghed(nobm1,nhtime))*c60) < 5) &
              .AND. (ABS( zsgbdy(nob,nbsp)-zsgbdy(nobm1,nbsp) ) < rvtlmt)) THEN
            iarlls (mzslon(nob),2)  =  - 3
            iairls (isortsl(irep),2,3)  =  0
            ilen = 5 + istrej
            IF (nacout+ilen <= nmxoln) THEN
              outbuf(nacout+1) = ilen
              outbuf(nacout+2) = nfmt17
              outbuf(nacout+3) = NINT( zsgbdy(nob,nbsp) )
              outbuf(nacout+4) = NINT( zsghed(nob,nhtime) * c100 )
              outbuf(nacout+5) = NINT( zsghed(nobm1,nhtime) * c100 )
              DO icl = 1 , istrej
                outbuf(nacout+5+icl) = ICHAR( yarlls(mzslon(nob)) (icl:icl) )
              ENDDO
              nacout  = nacout + ilen
            ENDIF
          ENDIF
        ENDIF
        IF (iarlls(mzslon(nob),2) >= 1)  irepcm1 = irep
      ENDIF
    ENDDO
  ENDIF

  IF (ntotssl > 0) THEN
    DEALLOCATE ( isortsl , STAT = nzerr )
    DEALLOCATE ( irepb   , STAT = nzerr )
    DEALLOCATE ( irepe   , STAT = nzerr )
    DEALLOCATE ( irepbb  , STAT = nzerr )
    DEALLOCATE ( irepee  , STAT = nzerr )
    DEALLOCATE ( mphase  , STAT = nzerr )
    DEALLOCATE ( izpp    , STAT = nzerr )
    DEALLOCATE ( iflftc  , STAT = nzerr )
    DEALLOCATE ( mphaso  , STAT = nzerr )
    DEALLOCATE ( zflat   , STAT = nzerr )
    DEALLOCATE ( fideh   , STAT = nzerr )
    DEALLOCATE ( fidev   , STAT = nzerr )
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_air_make_mult
!-------------------------------------------------------------------------------

END SUBROUTINE obs_air_make_mult


!-------------------------------------------------------------------------------
!+ Module procedure to shorten vertical correlation scale for dense reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_air_correl (ntotsml, ntotssl, vscale, tscale, hscale_c, hscale_v)

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure computes correction factors to shorten the vertical
!   correlation scales for aircraft reports, which are close to each other.
!
! Method:
!   Determination of correction factors to shorten the vertical correlation
!   scales at the top and/or bottom of aircraft reports if several of these
!   reports are from the same flight and are at close enough 4-D distance.
!   Such a factor is computed for each pair of reports:
!     f = 1 - w + w* min(cs/cn,1)
!     w : a gaussian weight in the temporal and horiz. distance betw. 2 reports
!     cn: correlation scale as given by the namelist
!     cs: a short correlation scale given by
!         cs = (max(dpp,dpm) /(2*p)) / (ln2)**.5
!         dpp: vertical distance betw. the 2 reports, in [Pa]
!         dpm: model layer thickness                , in [Pa]
!         p  : mean average between the 2 reports   , in [Pa]
!   With the use of cs, the vertical correlation takes the value 0.5 at half of
!   - the vert. distance between the 2 reports, or
!   - the model layer thickness
!
! Written by        :  Christoph Schraff, DWD  (original version: 25.06.98)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments:
! --------------------

  INTEGER, INTENT (IN)         ::  &
    ntotsml          ,& ! number of multi -level reports from the current flight
    ntotssl             ! number of single-level reports from the current flight

  REAL    (KIND=wp   )     , INTENT (IN)         ::  &
    vscale       (4) ,& ! scale (in [ln(p)] units) of the static Gaussian
                        ! vertical correlation function
                        ! --> usual and reasonable values: 0.577,0.577, 0.2, 0.2
    tscale           ,& ! scale [hrs] of the temporal Gaussian weight
    hscale_c     (4) ,& ! vertically invariant part [km]  \  of the scale of
    hscale_v     (4)    ! fraction of vertically variable  > the autoregressive
                        ! part given by tables 'rh?sond'  /  horizontal weight
                        ! ! where the total Gaussian weight (product of temporal
                        ! ! and horizontal weight) is given to the adapted
                        ! ! short correlation scale (which adapts to the
                        ! ! vertical distance between aircraft reports)
                        ! ! (1 minus the Gaussian weight is given to the large
                        ! !  static correlation scale)
                        ! --> usual and reasonable values for:
                        !     tscale_in   = 0.25
                        !     hscale_c_in = 0., 70., 0.  , 0.
                        !     hscale_v_in = 1.,  0., 0.83, 0.83

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER         ::  &
    ivar             ,& ! index of observed quantity: 1=(u,v); 3=T; 4=RH
    ivrs             ,& ! index of observed quantity: 1=(u,v); 2=T; 3=RH
    ncol   , ncom    ,& ! loop indices over reports
    nvrx             ,& ! bit position of (active) status for present variable
!   nbtxer , nbsxer  ,& ! ODR indices for mean obs. errors of present variable
    nhvcbx , nhvctx  ,& ! ODR indices for the correlation scale correction
    ntmb   , ntsb    ,& ! index of temporary multi- / single-level report b
    ntmc   , ntsc    ,& ! index of temporary multi- / single-level report c
    nlon             ,& ! report index in 'long' list
    nlev             ,& ! number of levels in multi-level report
    klev             ,& ! vertical loop index
    ncorlbb, ncorltb ,& ! base / top , report b:  \  indicates if correl. scale
    ncorlbc, ncorltc ,& ! base / top , report c:  /  need to be corrected
    nmean            ,& ! factor used to compute the relevant pressure level
    ilva   , ilvb    ,& ! vertical level indices
    itmp             ,& ! temporary station index in sorted list
    nzerr               ! status of memory (de-)allocation

  REAL    (KIND=wp   )     ::  &
    rinfl            ,& ! horiz. (correl.) scale relevant for the Gauss. weights
    rdegkm2          ,& ! (square of) conversion factor degrees to 'km'
    sq2, aln2, sqln2 ,& ! numerical constants
    halfwid          ,& ! half width in correlation scale units with the
                        ! horizontal correlation functions for the mass field
    obtimeb          ,& ! observation time of report b
    obtimec          ,& ! observation time of report c
    obbaltb, obtaltb ,& ! pressure at base / top of report b
    obbaltc, obtaltc ,& ! pressure at base / top of report c
    obilonb, objlatb ,& ! longitude / latitude of report b
    obilonc, objlatc ,& ! longitude / latitude of report c
    rdsprbb, rdsprtb ,& ! correlation scale correction at base / top of report b
    rdsprbc, rdsprtc ,& ! correlation scale correction at base / top of report c
    dppbb  , dpptb   ,& ! model layer thickness at base / top of report b
    dppbc  , dpptc   ,& ! model layer thickness at base / top of report c
    cosolat          ,& ! cos (latitude of report)
    tdiff            ,& ! time distance between reports (b and c)
    rhzob            ,& ! (square of) horizontal distance between reports
    zppobs           ,& ! pressure level where the correl. scale is adjusted
    dppobs           ,& ! layer thickness where the correl. scale is adjusted
    zlopco           ,& ! log( pressure ) where the correl. scale is adjusted
    rvinfl           ,& ! vertically dependent part of horiz. correlation funct.
    zf               ,& ! fraction for computing the vertically dependent part
    rhalfwid         ,& ! half width of horizontal correlation function [km]
    rhalfw2          ,& ! square of 'rhalfwid'
    fshorts          ,& ! weight given to the short correlation scale
    corshor          ,& ! the short vertical correlation scale
    frdpspr             ! correction factor to the vertical correlation scale
                        ! as from reports b and c

  LOGICAL                  ::  &
    lchange             ! the order of at least one pair of stations is changed

  CHARACTER (LEN=ilstidp)  ::  &
    ystid               ! station identity

! Local arrays:
! ------------

  INTEGER, ALLOCATABLE :: &
    isortsl    (:)      ! sorted list of reports with same station id
!
!------------ End of header ----------------------------------------------------


!-------------------------------------------------------------------------------
! Begin Subroutine obs_air_correl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1 : Preliminaries
!-------------------------------------------------------------------------------

! preset some variables
! ---------------------

  rdegkm2  = 111.1_wp **2
  sq2      = SQRT( c2 )
  aln2     = LOG ( c2 )
  sqln2    = SQRT( aln2 )
  halfwid  = 1.68_wp

! get station id
! --------------

  IF (ntotsml > 0) THEN
    ntmb  =  ismls  (1)
    nlon  =  mzmlon (ntmb,1)
    ystid =  yarlls (nlon)
  ELSEIF (ntotssl > 0) THEN
    ntsb  =  iairls (1,1,3)
    nlon  =  mzslon (ntsb)
    ystid =  yarlls (nlon)
  ENDIF


! sort reports so that their order is independent from the domain decomposition
! -----------------------------------------------------------------------------

  ALLOCATE ( isortsl (ntotsml + ntotssl) , STAT = nzerr )

! (Note: multi-level reports are already ordered independently
!        from the domain decomposition)
  DO ncom = 1 , ntotsml + ntotssl
    isortsl (ncom) = ncom
  ENDDO

  IF (ntotssl > 1) THEN
! sort longitudinally (stations are pre-sorted due to their origin nodes)
    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO ncom = ntotsml+1 , ntotsml+ntotssl-1
        IF (   mzsghd(iairls(isortsl(ncom  )-ntotsml,1,3),nhitot)              &
             > mzsghd(iairls(isortsl(ncom+1)-ntotsml,1,3),nhitot)) THEN
          itmp             = isortsl (ncom+1)
          isortsl (ncom+1) = isortsl (ncom)
          isortsl (ncom)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO
! sort latitudinally
    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO ncom = ntotsml+1 , ntotsml+ntotssl-1
        IF (      (   mzsghd(iairls(isortsl(ncom  )-ntotsml,1,3),nhitot)       &
                   == mzsghd(iairls(isortsl(ncom+1)-ntotsml,1,3),nhitot))      &
            .AND. (   mzsghd(iairls(isortsl(ncom  )-ntotsml,1,3),nhjtot)       &
                    > mzsghd(iairls(isortsl(ncom+1)-ntotsml,1,3),nhjtot))) THEN
          itmp             = isortsl (ncom+1)
          isortsl (ncom+1) = isortsl (ncom)
          isortsl (ncom)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
!  Section 2 : Provision of the required information for each pair of reports
!-------------------------------------------------------------------------------

! loop over reports (with present flight number)
! ----------------------------------------------

outer_loop_over_reports: DO ncol = 1 , ntotsml + ntotssl - 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! loop over all variables
! -----------------------

  DO ivrs = 1 , 3

    ivar = ivrs + MIN( ivrs-1 , 1 )
    IF (ivrs == 1) THEN
!     nbtxer = nbtuer
!     nbsxer = nbsuer
      nvrx   = nvru
      nhvcbx = nhvcbu
      nhvctx = nhvctu
    ELSEIF (ivrs == 2) THEN
!     nbtxer = nbtter
!     nbsxer = nbster
      nvrx   = nvrt
      nhvcbx = nhvcbt
      nhvctx = nhvctt
    ELSEIF (ivrs == 3) THEN
!     nbtxer = nbtqer
!     nbsxer = nbsqer
      nvrx   = nvrq
      nhvcbx = nhvcbq
      nhvctx = nhvctq
    ENDIF

! get all the required information on the current report
! ------------------------------------------------------

    obbaltb = rmdi
    IF (ncol <= ntotsml) THEN
      ntmb    = ismls  (ncol)
      obtimeb = zmlhed (ntmb,nhtime)
      obilonb = zmlhed (ntmb,nhilon)
      objlatb = zmlhed (ntmb,nhjlat)
      nlev    = mzmlhd (ntmb,nhnlev)
      DO klev = 1 , nlev
        IF (BTEST( mzmlbd(ntmb,klev,nbterr), nvrx )) THEN
          nlon      = mzmlon(ntmb,klev)
          IF (obbaltb <= rmdich) THEN
            obbaltb = zmlbdy(ntmb,klev,nbtp)
            dppbb   = rarlls(nlon)
          ENDIF
          obtaltb   = zmlbdy(ntmb,klev,nbtp)
          dpptb     = rarlls(nlon)
        ENDIF
      ENDDO
      rdsprbb = zmlhed (ntmb,nhvcbx)
      rdsprtb = zmlhed (ntmb,nhvctx)
    ELSE
      ntsb    = iairls (isortsl(ncol)-ntotsml,1,3)
      obtimeb = zsghed (ntsb,nhtime)
      obilonb = zsghed (ntsb,nhilon)
      objlatb = zsghed (ntsb,nhjlat)
      IF (BTEST( mzsgbd(ntsb,nbserr), nvrx )) obbaltb = zsgbdy (ntsb,nbsp)
      IF (BTEST( mzsgbd(ntsb,nbserr), nvrx )) obtaltb = zsgbdy (ntsb,nbsp)
      rdsprbb = zsghed (ntsb,nhvcbx)
      rdsprtb = zsghed (ntsb,nhvctx)
      nlon    = mzslon (ntsb)
      dppbb   = rarlls (nlon)
      dpptb   = rarlls (nlon)
    ENDIF
    IF (obbaltb > rmdi) THEN
      cosolat = cos( objlatb * r_degrad )

! inner loop: loop over all reports which have not yet
!             been checked against the current report
! ----------------------------------------------------

      inner_loop_over_reports: DO ncom = ncol + 1 , ntotsml + ntotssl
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! get all required information on second report
! ---------------------------------------------

        obbaltc = rmdi
        IF (ncom <= ntotsml) THEN
          ntmc    = ismls  (ncom)
          obtimec = zmlhed (ntmc,nhtime)
          obilonc = zmlhed (ntmc,nhilon)
          objlatc = zmlhed (ntmc,nhjlat)
          nlev    = mzmlhd (ntmc,nhnlev)
          DO klev = 1 , nlev
            IF (BTEST( mzmlbd(ntmc,klev,nbterr), nvrx )) THEN
              nlon      = mzmlon(ntmc,klev)
              IF (obbaltc <= rmdich) THEN
                obbaltc = zmlbdy(ntmc,klev,nbtp)
                dppbc   = rarlls(nlon)
              ENDIF
              obtaltc   = zmlbdy(ntmc,klev,nbtp)
              dpptc     = rarlls(nlon)
            ENDIF
          ENDDO
          rdsprbc = zmlhed (ntmc,nhvcbx)
          rdsprtc = zmlhed (ntmc,nhvctx)
        ELSE
          ntsc    = iairls (isortsl(ncom)-ntotsml,1,3)
          obtimec = zsghed (ntsc,nhtime)
          obilonc = zsghed (ntsc,nhilon)
          objlatc = zsghed (ntsc,nhjlat)
          IF (BTEST( mzsgbd(ntsc,nbserr), nvrx ))  obbaltc = zsgbdy (ntsc,nbsp)
          IF (BTEST( mzsgbd(ntsc,nbserr), nvrx ))  obtaltc = zsgbdy (ntsc,nbsp)
          rdsprbc = zsghed (ntsc,nhvcbx)
          rdsprtc = zsghed (ntsc,nhvctx)
          nlon    = mzslon (ntsc)
          dppbc   = rarlls (nlon)
          dpptc   = rarlls (nlon)
        ENDIF


!-------------------------------------------------------------------------------
!  Section 3 : Check if the reports of the current pair are sufficiently close
!              to each other, and preparation of the vertical correlation scales
!-------------------------------------------------------------------------------

! determine and check for temporal distance
! -----------------------------------------

        tdiff = c0
        IF (obbaltc > rmdi) tdiff = ABS( obtimec - obtimeb )
        IF ((obbaltc > rmdi) .AND. (tdiff < tscale*sq2)) THEN

! check for vertical distance, including determination of
! ---------------------------  - which vertical correlation scale (top / bottom)
!                                of which profile needs to be adjusted
!                              - the pressure level where this adjustment occurs

! ncorl?? = 0  ==>  correl. scale remains as specified by namelist
          ncorlbb = 0
          ncorltb = 0
          ncorlbc = 0
          ncorltc = 0
! the two profiles do not overlap
          IF (obbaltb < obtaltc-epsy) THEN
!   ncorl?? = 1  ==>  correl. scale depends on distance betw. profiles
            ncorlbb = 1
            ncorltc = 1
            zppobs  = c05 * (obbaltb + obtaltc)
            dppobs  = c05 * (dppbb   + dpptc  )
          ELSEIF (obtaltb > obbaltc+epsy) THEN
            ncorltb = 1
            ncorlbc = 1
            zppobs  = c05 * (obtaltb + obbaltc)
            dppobs  = c05 * (dpptb   + dppbc  )
          ELSE
! the two profiles overlap (but are not congruent)
            zppobs  = c0
            dppobs  = c0
!   ncorl?? = 2  ==>  correl. scale is set to a minimum value
            IF (obbaltb > obbaltc+epsy) ncorlbc = 2
            IF (obbaltb < obbaltc-epsy) ncorlbb = 2
            IF (obtaltb > obtaltc+epsy) ncorltb = 2
            IF (obtaltb < obtaltc-epsy) ncorltc = 2
            IF (ncorlbb == 2) THEN
              zppobs = obbaltb
              dppobs = dppbb
            ELSEIF (ncorlbc == 2) THEN
              zppobs = obbaltc
              dppobs = dppbc
            ENDIF
! assuming that 2 multi-level reports from the same flight do not usually
! overlap (substantially), only one (the mean) pressure level and model
! thickness is used (this is exact if one report is single-level)
            nmean = MAX( ncorlbb + ncorlbc , 1 )
            IF (ncorltb == 2) THEN
              zppobs = (obtaltb + zppobs) / nmean
              dppobs = (dpptb   + dppobs) / nmean
            ELSEIF (ncorltc == 2) THEN
              zppobs = (obtaltc + zppobs) / nmean
              dppobs = (dpptc   + dppobs) / nmean
            ENDIF
          ENDIF
          IF (MAX( ncorlbb, ncorlbc, ncorltb, ncorltc ) >= 1) THEN

! determine and check for horizontal distance
! -------------------------------------------

! horizontal correlation scale at pressure level 'zppobs'
            zlopco = LOG( zppobs )
            ilva = 1
            IF ((zlopco > tabcolp(1)) .OR. (zlopco <= tabcolp(ncolev))) THEN
              IF (zlopco <= tabcolp(ncolev)) ilva = ncolev
              IF (ivrs == 1) rvinfl = rhvsond(ilva)
              IF (ivrs == 2) rvinfl = rhtsond(ilva)
              IF (ivrs == 3) rvinfl = rhqsond(ilva)
            ELSE
              DO WHILE (zlopco <= tabcolp(ilva))
                ilva = ilva + 1
              ENDDO
              ilvb = ilva - 1
              zf   = (tabcolp(ilvb) - zlopco) / (tabcolp(ilvb) - tabcolp(ilva))
              IF (ivrs == 1) rvinfl = (c1-zf) *rhvsond(ilvb) + zf *rhvsond(ilva)
              IF (ivrs == 2) rvinfl = (c1-zf) *rhtsond(ilvb) + zf *rhtsond(ilva)
              IF (ivrs == 3) rvinfl = (c1-zf) *rhqsond(ilvb) + zf *rhqsond(ilva)
            ENDIF
            rinfl = MAX( hscale_c(ivar) + hscale_v(ivar) *rvinfl , epsy )
! check if the horizontal distance is small enough
            IF (ivar == 1) rhalfwid = rinfl * aln2
            IF (ivar >= 3) rhalfwid = rinfl * halfwid
            rhalfw2 = rhalfwid **2
            rhzob = (  (cosolat *(obilonc - obilonb))**2                       &
                     +           (objlatc - objlatb) **2 ) * rdegkm2
            IF (rhzob < c2*rhalfw2) THEN


!-------------------------------------------------------------------------------
!  Section 4 : Computation of the reduction factors to the correlation scales
!-------------------------------------------------------------------------------

! determine the weight to the short correlation scale
! ---------------------------------------------------

              fshorts =  EXP( - rhzob / rhalfw2 )                              &
                       * EXP( -(tdiff / tscale)**2 )

! determine the short correlation scale (which would result for zero temporal
! -------------------------------------  and horizontal distance betw. 2 rep.)

! determine the relevant vertical distance [Pa] (minimum: model layer thickness)
              IF (ncorlbb == 1) THEN
                dppobs  =  MAX( obtaltc - obbaltb , dppobs )
              ELSEIF (ncorlbc == 1) THEN
                dppobs  =  MAX( obtaltb - obbaltc , dppobs )
              ENDIF
! compute the short correlation scale
              corshor  =  c05 * dppobs / zppobs  /  sqln2

! determine the correction factor to the vertical correlation scale
! -----------------------------------------------------------------

              IF (corshor < vscale(ivar)-epsy) THEN
                frdpspr = (c1-fshorts) + fshorts *corshor /vscale(ivar)

! store the new factor where appropriate (i.e if it is smaller than the old one)
! --------------------------------------

                IF (ncorlbb >= 1) rdsprbb = MIN( rdsprbb , frdpspr )
                IF (ncorltb >= 1) rdsprtb = MIN( rdsprtb , frdpspr )
                IF ((ncorlbc >= 1) .AND. (frdpspr < rdsprbc)) THEN
                  IF (ncom <= ntotsml) THEN
                    zmlhed (ntmc,nhvcbx) = frdpspr
                  ELSE
                    zsghed (ntsc,nhvcbx) = frdpspr
                  ENDIF
                ENDIF
                IF ((ncorltc >= 1) .AND. (frdpspr < rdsprtc)) THEN
                  IF (ncom <= ntotsml) THEN
                    zmlhed (ntmc,nhvctx) = frdpspr
                  ELSE
                    zsghed (ntsc,nhvctx) = frdpspr
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF

! close inner loop over 'second reports'
! --------------------------------------

      ENDDO inner_loop_over_reports
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! store the correction factor in the 'current report'
! ---------------------------------------------------

      IF (ncol <= ntotsml) THEN
        zmlhed (ntmb,nhvcbx) = rdsprbb
        zmlhed (ntmb,nhvctx) = rdsprtb
      ELSE
        zsghed (ntsb,nhvcbx) = rdsprbb
        zsghed (ntsb,nhvctx) = rdsprtb
      ENDIF
      IF (lwonl)                                                               &
        WRITE( nupr,'(A,'': vert. correl.scale fac.'',I2,'', base/top:'',2F5.2 &
                    &,'', at'', F5.1,F5.1, F7.0, F7.0)' )                      &
               ystid, ivar, rdsprbb, rdsprtb, obilonb, objlatb, obbaltb, obtaltb
    ENDIF

! close outer loops over variables and 'current reports'
! ------------------------------------------------------

  ENDDO

ENDDO outer_loop_over_reports
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF (lwonl) THEN
    DO ncol = 1 , ntotsml + ntotssl
      IF (ncol <= ntotsml) THEN
        ntmb    = ismls  (ncol)
        WRITE( nupr,'(''crm'',I3,1X,A ,I2, 2I4,F7.0,F5.2, 4F8.5)' )            &
               my_cart_id, ystid  , mzmlhd(ntmb,nhpass)                        &
             , mzmlhd(ntmb,nhitot), mzmlhd(ntmb,nhjtot)                        &
             , zmlbdy(ntmb,1,nbtp), zmlhed(ntmb,nhtime)                        &
             , zmlhed(ntmb,nhvcbu), zmlhed(ntmb,nhvcbt)                        &
             , zmlhed(ntmb,nhvctu), zmlhed(ntmb,nhvctt)
      ELSE
        ntsb    = iairls (isortsl(ncol)-ntotsml,1,3)
        WRITE( nupr,'(''crs'',I3,1X,A ,I2, 2I4,F7.0,F5.2, 4F8.5)' )            &
               my_cart_id, ystid  , mzsghd(ntsb,nhpass)                        &
             , mzsghd(ntsb,nhitot), mzsghd(ntsb,nhjtot)                        &
             , zsgbdy(ntsb,  nbsp), zsghed(ntsb,nhtime)                        &
             , zsghed(ntsb,nhvcbu), zsghed(ntsb,nhvcbt)                        &
             , zsghed(ntsb,nhvctu), zsghed(ntsb,nhvctt)
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE ( isortsl , STAT = nzerr )

!-------------------------------------------------------------------------------
! End of module procedure obs_air_correl
!-------------------------------------------------------------------------------

END SUBROUTINE obs_air_correl


!-------------------------------------------------------------------------------
!+ Module procedure for lists of aircraft reports with specified requirements
!-------------------------------------------------------------------------------

SUBROUTINE obs_air_list ( ixcand , ncanda , ncande                             &
                        , ixrest , nresta , nreste , ixrequ , nrequ )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure produces lists of aircraft reports with specified
!   requirements.
!
! Method:
!   The elements of the lists point to aircraft reports which must meet the
!   requirements as specified by the input parameters.
!
! Written by        :  Christoph Schraff, DWD  (original version: 24.06.98)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Subroutine arguments:
! --------------------

  INTEGER, INTENT (IN)         ::       &
    ixcand           ,& ! specifies the 'candidate list' (CANDLS) from which
                        ! elements (reports) may be added to the 'required list'
                        ! if (ABS(ixcand) <= 7): CANDLS is 'iairls(,,ixcand)'
                        ! if (ABS(ixcand) == 8): CANDLS is 'iarlls' (long list)
                        ! if (ABS(ixcand) == 9): CANDLS is ODR
                        ! if (ixcand < 0) then only active reports from the
                        !   'candidate list' may be added to the 'required list'
                        ! (redundant reports are never used)
    ncanda , ncande  ,& ! specify beginning and end of the 'candidate list'
                        ! data at current timestep
    ixrest           ,& ! index of the 'restriction list'.
                        ! If (ixrest >= 1) then a 'candidate' may only be added
                        ! to the 'required list', if its station id. is in the
                        ! interval [nresta,nreste] of the 'restriction list'
    nresta , nreste  ,& ! active interval of the 'restriction list', cf ixrest
    ixrequ              ! index of the 'required list'

  INTEGER, INTENT (INOUT)      ::       &
    nrequ               ! number of elements in the 'required list'


! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER         ::  &
    mxcand           ,& ! ABS( ixcand )
    nrestev          ,& ! variable start of active interval in 'restriction list
    ixrestv          ,& ! variable index of the 'restriction list'
    ntoomany         ,& ! number of reports in excess of the list size
    ncand            ,& ! loop index over candidates' list
    nrest            ,& ! loop index over restriction list
    ncdodr           ,& ! ODR or list index of candidate
    kcdtyp              ! observation code type

  LOGICAL                  ::  &
    ldeflt           ,& ! empty 'restriction list'
                        ! (arbitrary candidate's station id)
    lsameid             ! candidate's station id's is in restriction list

  CHARACTER (LEN=ilstid)   ::  ycdstid     ! station id
  CHARACTER (LEN=50)       ::  ycomment    ! comment, 1st part
  CHARACTER (LEN=29)       ::  ycommen2    ! comment, 2nd part

! Local (automatic) arrays: None
! -------------------------
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_air_list
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1  : Preliminaries
!-------------------------------------------------------------------------------

  mxcand  = ABS( ixcand )

! ensure consistency of parameters defining the check list of station id's

  nrestev = nreste
  ixrestv = ixrest
  IF ((ixrestv < 1) .OR. (ixrestv > 3)) nrestev = nresta - 1
  IF (nrestev < nresta) ixrestv = 0

! for empty 'restriction list', the candidate's station id is arbitrary

  ldeflt  =  (ixrestv == 0)

! preset further variables

  ntoomany = 0

!-------------------------------------------------------------------------------
!  Section 2  : Production of the required list
!-------------------------------------------------------------------------------
!  Section 2.1: Search with inner loop: if the 'active' part of the
!               'restriction list' contains more than one element
!-------------------------------------------------------------------------------

  IF (nrestev > nresta) THEN

    loop_over_candidates: DO ncand = ncanda, ncande

! determine the 'candidate list', and get (ODR) index, code type,
! and station index of the candidate
! ---------------------------------------------------------------

      IF (mxcand == 9) THEN
        ncdodr   =  ncand
        kcdtyp   =  mosghd (ncand,nhcode)
        ycdstid  =  yosghd (ncand)
        IF (mosghd(ncand,nhpass) == 2)                      kcdtyp = -1
        IF (BTEST( mosghd(ncand,nhflag), FL_FLIGHTTRACK ))  kcdtyp = -2
        IF (BTEST( mosghd(ncand,nhflag), FL_THIN        ))  kcdtyp = -3
        IF (BTEST( mosghd(ncand,nhflag), FL_REDUNDANT   ))  kcdtyp = -4
        IF (BTEST( mosghd(ncand,nhflag), FL_HEIGHT      ))  kcdtyp = -4
        IF (mosghd(ncand,nhobtp) /= nairep)                 kcdtyp = -5
      ELSEIF (mxcand == 8) THEN
        ncdodr   =  iarlls (ncand,1)
        kcdtyp   =  iarlls (ncand,2)
        ycdstid  =  yarlls (ncand  )
      ELSE
        ncdodr   =  iairls (ncand,1,mxcand)
        kcdtyp   =  iairls (ncand,2,mxcand)
        ycdstid  =  yairls (ncand  ,mxcand)
      ENDIF

! check if the candidate is an 'active' aircraft report on the 'candidate list',
! and if the station id is admitted to the 'required list', i.e. if it is on the
! 'restriction list'
! ------------------------------------------------------------------------------

      IF (     (kcdtyp == naircd) .OR. (kcdtyp == ncodar)                      &
          .OR. (kcdtyp == namdar) .OR. (kcdtyp == nacar)                       &
          .OR. (kcdtyp == nmodes)                                              &
! reports which were set passive due to observation (code) type selection or
! flight track check are also put on the list (for further flight track checks)
! if list may not only include active reports
          .OR. ((kcdtyp <= 0) .AND. (kcdtyp >= -2) .AND. (ixcand > 0))) THEN
        lsameid  = ldeflt
        DO nrest = nresta , nrestev
          IF (ycdstid == yairls(nrest,ixrestv)) lsameid = .TRUE.
        ENDDO
        IF ((lsameid) .AND. (nrequ <= nairls)) THEN

! add the candidate to the 'required list',
! if required, set the successful candidate to 'passive' on the 'candidate list'
! and indicate if there is not enough space for the 'required list'

          nrequ = nrequ + 1
          iairls (nrequ,1,ixrequ) = ncdodr
          iairls (nrequ,2,ixrequ) = kcdtyp
          yairls (nrequ  ,ixrequ) = ycdstid
        ELSEIF (lsameid) THEN
          ntoomany = ntoomany + 1
        ENDIF
      ENDIF

    ENDDO loop_over_candidates

  ELSE

!-------------------------------------------------------------------------------
!  Section 2.2: Search without inner loop: if the elements of the 'required list
!                  - have arbitrary station id's (no active 'restriction list')
!               or - must have one specific station id. (1 active element in the
!                                                        'restriction list')
!-------------------------------------------------------------------------------

    lsameid = ldeflt
    DO ncand = ncanda, ncande
      IF (mxcand == 9) THEN
        ncdodr   =  ncand
        kcdtyp   =  mosghd (ncand,nhcode)
        ycdstid  =  yosghd (ncand)
        IF (mosghd(ncand,nhpass) == 2)                      kcdtyp = -1
        IF (BTEST( mosghd(ncand,nhflag), FL_FLIGHTTRACK ))  kcdtyp = -2
        IF (BTEST( mosghd(ncand,nhflag), FL_THIN        ))  kcdtyp = -3
        IF (BTEST( mosghd(ncand,nhflag), FL_REDUNDANT   ))  kcdtyp = -4
        IF (BTEST( mosghd(ncand,nhflag), FL_HEIGHT      ))  kcdtyp = -4
        IF (mosghd(ncand,nhobtp) /= nairep)                 kcdtyp = -5
      ELSEIF (mxcand == 8) THEN
        ncdodr   =  ncand
        kcdtyp   =  iarlls (ncand,2)
        ycdstid  =  yarlls (ncand  )
      ELSE
        ncdodr   =  iairls (ncand,1,mxcand)
        kcdtyp   =  iairls (ncand,2,mxcand)
        ycdstid  =  yairls (ncand  ,mxcand)
      ENDIF
      IF (     (kcdtyp == naircd) .OR. (kcdtyp == ncodar)                      &
          .OR. (kcdtyp == namdar) .OR. (kcdtyp == nacar)                       &
          .OR. (kcdtyp == nmodes)                                              &
! reports which were set passive due to observation (code) type selection or
! flight track check are also put on the list (for further flight track checks)
! if list may not only include active reports
          .OR. ((kcdtyp <= 0) .AND. (kcdtyp >= -2) .AND. (ixcand > 0))) THEN
        IF (.NOT. ldeflt) lsameid  =  ycdstid == yairls(nresta,ixrestv)
        IF ((lsameid) .AND. (nrequ <= nairls)) THEN
          nrequ = nrequ + 1
          iairls (nrequ,1,ixrequ) = ncdodr
          iairls (nrequ,2,ixrequ) = kcdtyp
          yairls (nrequ  ,ixrequ) = ycdstid
        ELSEIF (lsameid) THEN
          ntoomany = ntoomany + 1
        ENDIF
      ENDIF
    ENDDO

  ENDIF

!-------------------------------------------------------------------------------
!  Section 3  : Print message if there was not enough space in the required list
!-------------------------------------------------------------------------------

  IF (ntoomany >= 1) THEN
    IF (mxcand == 9) THEN
      WRITE(ycomment,'(I4,'' too many active single-level aircraft reports'')' &
                      )   ntoomany
    ELSEIF (mxcand == 8) THEN
      WRITE(ycomment,'(I4,'' too many aircraft reports in the "long" list '')' &
                      )   ntoomany
    ELSE
      WRITE(ycomment,'(I4,'' too many active airc. reports found in list''     &
                     &,I2)' ) ntoomany , ixcand
    ENDIF
    IF (nrestev < nresta) THEN
      WRITE( ycommen2, '('' with arbitrary station ID'')' )
    ELSEIF (nrestev == nresta) THEN
      WRITE( ycommen2, '('' with station ID :   '',A)' ) yairls(nresta,ixrestv)
    ELSE
      WRITE( ycommen2, '('' with same ID as in list '',I1)' ) ixrestv
    ENDIF
    PRINT '(''CAUTION with preparation of multi-level aircraft ''              &
          &,''reports (obs_air_list):'')'
    PRINT '(A50,A26)', ycomment, ycommen2
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_air_list
!-------------------------------------------------------------------------------

END SUBROUTINE obs_air_list

!===============================================================================
!+ Module procedure in "src_obs_proc_air.f90" for interface to model environmt.
!-------------------------------------------------------------------------------

SUBROUTINE obs_air_interface ( zacthr , imaxmll, imaxsgl, imaxmlv, llwonl      &
                             , zdegrad, nlbc_lines, npe, my_id, icomm          &
                             , mpi_reals, mpi_integers )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_proc_air.f90" makes available the values
!   of those input variables (and parameters) from module 'data_obs_lib_cosmo'
!   which are set in (and depend on) other modules from the model environment.
!   This step makes this module independent from the model environment.
!   Input variables are not only those variables used explicitly in this module
!   but also variables used in routines (i.e. 'obs_assign_sort_node') called
!   in this module.
!
! Method:
!   This subroutine has to be called once prior to calling 'obs_air_org_mult'.
!   It does not have to be called at all, only if routine 'obs_cdf_interface'
!   from module 'src_obs_cdfin_org.f90' has already been called previously.
!   has to be called once prior to calling 'obs_air_org_mult'.
!   In the call, the argument list consists of variables from the 'model'
!   environment. Inside this routine, the values in the argument list are
!   assigned to the variables from module 'data_obs_lib_cosmo'.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER, INTENT (IN)   ::  &
    imaxmlv        ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    imaxmll        ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl        ,& ! size (report dimension) of the single-level (s-l)  ODR
    ! variables related to parallelisation / domain decomposition
    nlbc_lines     ,& ! number of overlapping boundary lines of the subdomains
    npe            ,& ! number of compute PEs
    my_id          ,& ! rank of this subdomain in the cartesian communicator
    icomm          ,& ! communicator for the virtual cartesian topology
    mpi_reals      ,& ! REAL      type used for MPI
    mpi_integers      ! INTEGER   type used for MPI

  REAL    (KIND=wp)        , INTENT (IN)   ::  &
    zdegrad        ,& ! factor for transforming degree to rad
    zacthr            ! actual model time [hours] with respect to 'yydate_ref'

  LOGICAL                  , INTENT (IN)   ::  &
    llwonl            ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open

! Local parameters: None
! ----------------

! Local scalars: None
! -------------

! Local arrays: None
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_air_interface
!-------------------------------------------------------------------------------

! provide values to input variables before 'obs_air_org_mult' is called

  num_compute  =  npe
  nboundlines  =  nlbc_lines
  my_cart_id   =  my_id
  icomm_cart   =  icomm
  imp_reals    =  mpi_reals
  imp_integers =  mpi_integers
  r_degrad     =  zdegrad
  maxmlv       =  imaxmlv
  maxmll       =  imaxmll
  maxsgl       =  imaxsgl
  acthr        =  zacthr
  lwonl        =  llwonl

!-------------------------------------------------------------------------------
! End Subroutine obs_air_interface
!-------------------------------------------------------------------------------

END SUBROUTINE obs_air_interface

!-------------------------------------------------------------------------------

END MODULE src_obs_proc_air


