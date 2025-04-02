!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfin_gps

!-------------------------------------------------------------------------------
! Description:
!   This module performs the observation pre-processing of gps reports
!   read from NetCDF observation input files.
!   Special tasks are:
!    - reading reports from NetCDF files, selecting, assigning them to model
!      grid points and distributing them to the nodes (precessors) according
!      to the sub-domain which they lie on
!    - pre-processing, including some gross error checking, blacklist checking,
!      assigning observation errors, static bias correction, etc.
!    - storing the reports in the internal ODR arrays
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_read_gps   : (called by obs_cdf_read_org)
!    - obs_cdf_store_gps  :  called by obs_cdf_read_gps
!
!   This module also contains elemental functions, formerly statement functions:
!   - insert       : inserts bits set in one integer into another integer word
!                    at given bit position
!
!   It uses from:
!    - src_obs_cdfin_comhead: - obs_cdf_read_comhead
!                             - obs_cdf_buffer_comhead
!                             - obs_cdf_store_comhead
!    - src_obs_cdfin_blk:     - obs_cdf_whitelist_local
!                             - obs_cdf_blacklist_local
!    - src_obs_cdfin_util:    - obs_assign_sort_node
!                             - obs_cdf_distrib_reports
!                             - obs_td2rh , obs_qx2rh , obs_rhw2rh
!                             - obs_find_level
!                             - f_z2p
!    - parallel_utilities:    - distribute_values
!    - environment:           - model_abort
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin
!    - data_obs_record
!    - mo_fdbk_tables
!
! Current Code Owner (for COSMO and for DACE):
!  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012/01/31 Karolin Eichler + Christoph Schraff
!  Initial release by KE, based on src_obs_cdf_sing, modified by CS.
! V4_28        2013/07/12 Christoph Schraff
!  Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Bug fix: call of obs_rhw2rh made conditioned to avoid referencing
!           unallocated arrays ztt, zrhw.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  Variables related to the AOF interface removed.
! V5_4d        2016-12-12 Christoph Schraff
!  Discard reports completely for which both ZPD and IWV are missing.
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Editorial changes only.
! @VERSION@    @DATE@     Christoph Schraff
!  - Set subrout. (obs_cdf_read_gps) used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
!  18.04.23: 'rhtsat' moved from data_obs_cdfin to data_obs_lib_cosmo.
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

USE netcdf       ! NetCDF F90 interface

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Scalar variables originally defined in other data modules,
!    obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------

    ! horizontal and vertical sizes of the model fields
    ke             ,& ! number of grid pts in vertical direction (--> 'f_z2p')

    ! variables related to parallelisation / domain decomposition
    num_compute    ,& ! number of compute PEs
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers   ,& ! INTEGER   type used for MPI

    ! report dimensions of ODR (Observation Data Record, for observation storage
    ! on local sub-domains), and other variables related to namelist parameters
    madj_hum          ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable

USE data_obs_lib_cosmo, ONLY :  &

    ! constants for the horizontal rotated grid and related variables
    r_pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    r_pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
    r_polgam       ,& ! angle between the north poles of the systems
    r_dlon         ,& ! grid point distance in zonal direction (in degrees)
    r_dlat         ,& ! grid point distance in meridional direction (in degrees)
    r_startlat_tot ,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    r_startlon_tot ,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    r_degrad       ,& ! factor for transforming degree to rad

    ! variables related to namelist parameters, and other variables
    doromx         ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    altopsu        ,& ! SYNOP obs. above height 'altopsu' are set passive
    rhtsat         ,& ! relative humidity threshold above which obs is set =100%

    ! physical constants
!   r_g            ,& ! acceleration due to gravity
    tmelt          ,& ! melting temperature of ice
    r_d            ,& ! gas constant for dry air
!   rdv            ,& ! r_d / r_v
    b1             ,& ! variables for computing the saturation vapour pressure
    b2w            ,& ! over water (w) and ice (i)
    b2i            ,& !               -- " --
    b3             ,& !               -- " --
    b4w            ,& !               -- " --
    b4i            ,& !               -- " --

    ! switches related to namelist parameters and other
    lverpas        ,& ! write also passive reports on feedback files

! 2b. Pointers for arrays originally defined in other data modules
! ----------------------------------------------------------------

    ! array related to parallelisation / domain decomposition
    i_subpos       ,& ! positions of the subdomains in the total domain
                      ! (i-, j-indices of the lower left + upper right grid pt.
                      !  in the order: i_ll, j_ll, i_ur, j_ur; only the domain
                      !  interior is considered, not the boundary lines)

    ! model fields
    r_p            ,& ! pressure (at main levels)
    r_hhl          ,& ! height   (at half levels)
    r_t_ll         ,& ! temperature at lowest model (main) level
    r_ps              ! surface pressure

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    ngps       ,& ! GPS reports

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

    cma        ,& ! array of meta data on CMA observation and code types

! 7. Functions
! ------------

    i_cma         ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

! Section 1 : NetCDF Observation Input File formats
!             (this section is used ONLY IF obs data are read from NetCDF files)
!-------------------------------------------------------------------------------

!         1.1   internal attributes of the different NetCDF input files
!               -------------------------------------------------------

    icdfinlen      ,& ! maximum length of NetCDF observation input file name
    iannexlen      ,& ! maximum length of annex of NetCDF obs input file name
    ncdf_gps_zenith,& ! indicator for proc. NetCDF GPS (ZPD / IWV)    input
    ycdfin         ,& ! file names of NetCDF observation input files
    icdfin         ,& ! obs file type of NetCDF observation input files
    ncinid         ,& ! unit numbers of NetCDF observation input files
    yncannex       ,& ! annex of NetCDF observation input file names
    dimids            ! dimension IDs in NetCDF files

USE data_obs_cdfin, ONLY :  &

!         1.2   variables used to read from the data section of the NetCDF files
!               ----------------------------------------------------------------

!         1.2.1 'common' NetCDF header entries and derived variables
!               ----------------------------------------------------

    ilstidn        ,& ! character length of station identity from NetCDF files
! derived common header variables stored to ODR
    iobstot        ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot        ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc        ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc        ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc            ! indices of reports to be processed now

!         1.2.2 other NetCDF header entries
!               ---------------------------

!    nc_nix            ! NIX   , BUFR Table B 002001 : station type (man,auto,..)

USE data_obs_cdfin, ONLY :  &

!         1.3   NetCDF body entries
!               -------------------

!         1.3.1  frequent entries
!                ----------------
    rc_p           ,& ! MPN   , MPPP : pressure
    rc_t           ,& ! MTN   : temperature / dry-bulb temperature
    rc_lat0        ,& ! MLAH  : latitude
    rc_lon0        ,& ! ML0H  : longitude
    nc_rh          ,& ! MUUU  : relative humidity                           [%]
    nc_vdt         ,& ! NGGTP : time period of wind measurement           [min]
    rc_el          ,& ! MDE   : elevation                                   [°]
    rc_iwv         ,& ! NWLN  : integrated water vapour                 [kg/m²]
    rc_zwd         ,& ! NCZWV : zenith wet delay                            [m]
    rc_zpder       ,& ! NEERR : estimated error in zpd                      [m]
    rc_zpd         ,& ! NADES : zenith path delay                           [m]
    nc_flg         ,& ! NQFGD : quality flag


!         1.5   auxilliary buffer arrays and variables
!               --------------------------------------

    imiss          ,& ! missing value for integers in current NetCDF input file
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk          ! value smaller than 'rmiss', to check for missing value

USE data_obs_cdfin, ONLY :  &

! Section 2 : Blacklist and Whitelist
!-------------------------------------------------------------------------------

    ilstid_blk     ,& ! assume 8-character station-IDs in Black-/Whitelist
    blk_loc           ! blacklists for local reports

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
!-------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
    nezdif     ,& ! distance 'model orography - station altitude' too large
    nenoda     ,& ! no accepted data in report
    nenops     ,& ! pressure too small (< 20hPa), or missing in aircraft report

!         3.1.2   Data event counter array format
!                 -------------------------------
    mxdeve     ,& ! length of data event counter array
    nelodr     ,& ! level rejected: number of levels exceeding ODR size
    nelmis     ,& ! level rejected: pressure (PILOT: pressure + height) missing
    nelflg     ,& ! level rejected: pressure (PILOT: height) flagged
    nelsfc     ,& ! level rejected: too many surface levels
    nelnop     ,& ! level rejected: PILOT height level not in range of model lev
    nelext     ,& ! level rejected: pressure < 9hPa or level below station alt.
    nelsig     ,& ! level rejected: significant level above a specified limit
    nelrdn     ,& ! level rejected: redundant level in report  (not active yet)
    nepmis     ,& ! pressure (TEMP: height): missing
    nepflg     ,& ! pressure (TEMP: height): flagged
    neprac     ,& ! pressure: bad reporting practice
    nepalt     ,& ! pressure: sta height, or height dist. to orography too large
    nepsdt     ,& ! pressure tendency: flagged, or absolute value too large
    netmis     ,& ! temperature missing
    netflg     ,& ! temperature flagged
    netext     ,& ! temperature too low or too high
    netalt     ,& ! height (diff.) too large for 2m-temp.
    netlps     ,& ! lapse rate of multi-level temperature too large
    neqmis     ,& ! humidity missing
    neqflg     ,& ! humidity flagged
    neqlow     ,& ! humidity too low
    neq300     ,& ! humidity above 300 hpa
    neqbig     ,& ! humidity over allowed value (120%)
    neqsap     ,& ! humidity forced to be saturated (t>0)
    neqsam     ,& ! humidity forced to be saturated (t<0)
    neqclp     ,& ! humidity forced to be <=100% (t>0)
    neqclm     ,& ! humidity forced to be <=100% (t<0)
    neqalt     ,& ! height (diff.) too large for 2m-humid
    nevalt     ,& ! height (diff.) too large for 10m-wind / IWV
    negmis        ! ZPD missing or too small

USE data_obs_cdfin, ONLY :  &

!         3.2    Event counter arrays
!                --------------------
    neventr    ,& ! counter array of report events
    neventd    ,& ! counter array of data events

!         3.4    Variables' expectation table
!                ----------------------------
    nzex       ,& ! expect geopotential
    nuex       ,& ! expect horiz. wind
    ntex       ,& ! expect temperature
    ntdex         ! expect humidity

USE data_obs_cdfin, ONLY :  &

!         4.1  Observation error levels
!              ------------------------
    nerlev     ,& ! number of standard error levels
    rolnlv     ,& ! ln(rlevel(15))

!
!         4.2  Observation error constants
!              ---------------------------
    oezgps     ,& ! (root of) GPS height error variance (land)
    rherr1     ,& ! (root of) fixed    / normal conditions
    rherr2     ,& ! relative humidity <  if temperature below 233K
    rherr3     ,& ! error variances    \ if rel. humidity below 20%
    oezpd         ! GPS Zenith Path Delay default error variance

USE data_obs_cdfin, ONLY :  &

!         5.2    Temperature / humidity / pressure / height / fog limits
!                -------------------------------------------------------
    rttlim     ,& ! temperature limit below which rel. hum. error is increased
    rrhlim     ,& ! rel. hum. limit below which rel. hum. error is increased
    rtshlm     ,& ! gross error upper limit for relative humidity
    rerrpp     ,& ! msl pressure above which observed pressure is
                  ! reckoned to be erroneous
    pminsigt   ,& ! significant-level TEMP / PILOT data are neglected unless
    pminsigv   ,& !   p-obs >= pminsigt (10000 pa) and temperature obs exists or
                  !   p-obs >= pminsigv (20000 pa) and wind obs exists
    pqmin      ,& ! pressure [pa] of level above which moisture obs are not used
    rpplim     ,& ! pressure level obove which observations are not used

!         7.0    For data rejection messages: Output buffer, size and formats
!                ------------------------------------------------------------
    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt1      ,& ! no pressure
    nfmt3      ,& ! no accepted data
    nfmt6         ! excess of pressure tendency

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'ogphed'
!              ----------------------------------------------------
!   mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
!   nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
!   nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mogphd'
!              ----------------------------------------------------
!   mxghdf     ,& ! header length of GPS reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
!   nhitot     ,& ! global x-coord. of grid pt. assigned to obs
!   nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
!   nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
!   nhqcfw     ,& ! status of QC and of writing to feedobs files (and p-QC flag)
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
!   nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
!   nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre
!   nhstid     ,& ! station identity number
!   nhdate     ,& ! absolute exact observation date [yyyymmdd]
!   nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
!   nhstyp     ,& ! surface obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)


!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid        ! character length of the station identity
!   ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.3  Body format of ODR of surface reports: 'ogpbdy'
!              -----------------------------------------------
!   mxgbdy     ,& ! body length of gps reports
    nbgt       ,& ! temperature                                        [K]
    nbgrh      ,& ! relative humidity                                  [/]
    nbgp       ,& ! pressure                                           [Pa]
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay                                   [mm]
    nbgiwv     ,& ! integrated water vapour                    [kg/m2 ~ mm]
    nbgbia     ,& ! bias correction to integrated water vapour         [mm]
    nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour  [mm]
    nbgz       ,& ! height
    nbgzer     ,& ! error of observed height [m]
    nbgqer     ,& ! error of observed relative humidity
    nbgdrh     ,& ! bias correction for relative humidity
    nbgter     ,& ! error of observed temperature
    nbgviz     ,&

!       1.3.4  Body format of ODR of surface report flags: 'mogpbd'
!              ----------------------------------------------------
!   mxgbdf     ,& ! body length of GPS reports
    nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid     ,& ! level identity          (bit pattern, see below: 'nb?lid')
    nbgerr        ! pre-processing data flag


USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
!   nvru       ,& ! bit pos. for status/QC flags for horiz. wind   nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature         "
    nvrq       ,& ! bit pos. for status/QC flags for humidity            "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height     "
!   nvrw       ,& ! bit pos. for status/QC flags for vertical wind       "
    nvriwv     ,& ! bit pos. for status/QC flags for iwv                 "
    nvrzpd     ,& ! bit pos. for status/QC flags for zenith path delay   "
!   nvfubp     ,& ! bit pos. for main flag on wind                     nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature                "
    nvfqbp     ,& ! bit pos. for main flag on humidity                   "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfgbp     ,& ! bit pos. for main flag on IWV / ZPD                  "
!   nvfaoc     ,& ! no. of bits occ. by each main flag                   "
    nvfbps     ,& ! bit pattern for main flags:                          "
!   nvfboc     ,& ! no. of bits occ. for main flags                      "
!   nvflbp     ,& ! bit pos. for level flag: level below surface         "
    nvlidp        ! level id. bit pattern                              nb?lid
!   nvlido        ! no. bits occ. by each indicator in level id.         "

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotgp     ,& ! tot. number of stored single-level reports
    fdoro      ,& ! scaling factors to obs height minus model orography

!       2.     Observation data records (ODR) and associated arrays
!       ----------------------------------------------------------------------

!       2.1    Formats of observation data records
!              -----------------------------------
    ogpbdy     ,& ! body   of gps ODR
    ogphed     ,& ! header of gps ODR
    mogpbd     ,& ! body   of gps ODR
    mogphd     ,& ! header of gps ODR
    yogphd        ! header of gps ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_NO_OBS     ! no (active) observations in report

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
!   global_values,   & ! computes global values by operating on local arrays
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_sort_node   ,& ! assign node to reports and sort them accordingly
    obs_cdf_distrib_reports,& ! distribute reports to appropriate sub-domains
!   obs_td2rh              ,& ! convert dewpoint temperature to relat. humidity
!   obs_qx2rh              ,& ! convert mixing ratio to model-compatible rel hum
    obs_rhw2rh             ,& ! make (observed) relat. humidity model compatible
    obs_find_level         ,& ! get interpolation levels and factors
    f_z2p                     ! get (model) pressure for a specified height

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_comhead,       ONLY :  &
    obs_cdf_read_comhead   ,& ! reads and evaluates common header information,
                              !   including grid point assignment of reports
    obs_cdf_buffer_comhead ,& ! puts header info in buffer arrays for distribut.
    obs_cdf_store_comhead     ! stores common header info in ODR

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_blk,           ONLY :  &
    obs_cdf_whitelist_local,& ! indic. for local reports if missing on whitelist
    obs_cdf_blacklist_local   ! produces a list of blacklisted vertical
                              !   intervals for each of a set of local reports
!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_cdf_read_gps

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

!===============================================================================
!+ Module procedure in "src_obs_cdfin_gps" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_gps  ( min_sta , min_end , ilcdf , ngpnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_gps" organizes the reading,
!   pre-selection, distribution and storage in ODR of gps reports
!   from a given observation time period.
!
! Method:
!   The reports are read by 1 node (processor) and assigned each to its most
!   appropriate grid point of the total model domain. The reports then have
!   to be distributed to those sub-domains which contain the grid point they
!   are assigned to. To do this, the reports are temporarily stored in three
!   long buffer arrays (one each for integer, real, and character elements)
!   in the order according to these sub-domains, and then the appropriate
!   sections of these arrays are distributed to the nodes.
!   Pre-processing steps related to the observation body elements and storing
!   the data in long-term arrays (ODR) are then performed locally, i.e. in
!   parallel mode.
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

  INTEGER        , INTENT (IN)  ::  &
    min_sta       ,& ! start  \  of time interval to be read now
    min_end       ,& ! end    /  (in [minutes] of forecast time)
    ilcdf            ! index (number) of NetCDF observation input file

  INTEGER        , INTENT (INOUT)  ::  &
    ngpnew        ,& ! number of new GPS reports
    nexceed          ! number of reports in excess of array size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    nodrepn (num_compute) ,& ! number of            reports  \   to be
    nodleni (num_compute) ,& ! number of integer   elements   \  distributed
    nodlenr (num_compute) ,& ! number of real      elements   /  to the
    nodleny (num_compute) ,& ! number of character elements  /   different nodes
    nrepl                 ,& ! number of            reports  \  ( to be
    nlenli                ,& ! number of integer   elements   \   received )
    nlenlr                ,& ! number of real      elements   /  at the
    nlenly                ,& ! number of character elements  /   local node
    ibuflen  (3)          ,& ! length of (i/r/c) buffer arrays for distrib. rep.
    iloclen  (3)             ! length of (i/r/c) buffer arrays received locally

  INTEGER        ::  &
    kcdftyp        ,& ! type of NetCDF input file format (observation type):
                      !  = ncdf_gps_zenith     : GPS
    ncid           ,& ! file unit number of current NetCDF obs input file
    mrepsta        ,& ! index (in record dimension) of 1st report to be read
    nrep           ,& ! index interval (in record dimension) which contains
                      !   those reports that should be read now
    nreproc        ,& ! number of reports to be read and processed now
    noffscal (3,2) ,& ! number int/real/char - elements in common header part
                      !                      - scalar elements in total report
    iioffs         ,& ! offset of current report in long integer buffer
    iroffs         ,& ! offset of current report in long real    buffer
    niscal         ,& ! number of scalar integer   elements
    nrscal         ,& ! number of scalar real      elements
    nyscal         ,& ! number of scalar character elements
    inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    irps           ,& ! loop index over reports to be processed
    irep           ,& ! record index over reports to be processed
    status         ,& ! NetCDF status variable
!   ktimsig        ,& ! significance of report time
    istat  , ierr     ! error indicators

  INTEGER        ::  &          ! variable ID's in NetCDF file:
    varid_t                 ,& ! dry bulb temperature
    varid_rh                ,& ! relative humidity
    varid_p                 ,& ! station pressure
    varid_vdt               ,& !
    varid_el                ,& ! elevation
    varid_zpd               ,& ! zenith path delay
    varid_zpder             ,& ! error of zenith path delay
    varid_iwv               ,& ! integrated water vapour
    varid_zwd               ,& ! zenith wet delay
    varid_flg                  ! quality flag

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
! CHARACTER (LEN=30)       :: &
!   yerr              ! error message
  CHARACTER (LEN=70)       :: &
    yerrmsl           ! error message
! CHARACTER (LEN=80)       :: &
!   ymsg              ! control message
  CHARACTER (LEN= icdfinlen + iannexlen)  ::  &
    yfn               ! file name of NetCDF observation input file, with annex

! Local arrays:
! ------------

  INTEGER        , ALLOCATABLE :: &
    irnode     (:) ,& ! nodes to which the reports will be distributed
    irsort     (:) ,& ! report indices sorted according to 'irnode'
    iirlen     (:) ,& ! number of integer   elements  \   in each
    irrlen     (:) ,& ! number of real      elements   >  individual
    iyrlen     (:)    ! number of character elements  /   report

  INTEGER        , ALLOCATABLE :: &
    ibufsrt    (:) ,& ! integer buffer array with sorted reports read from file
    ibufloc    (:)    ! integer buffer array with local reports (at sub-domain)

  REAL (KIND=wp) , ALLOCATABLE :: &
    rbufsrt    (:) ,& ! real    buffer array with sorted reports read from file
    rbufloc    (:)    ! real    buffer array with local reports (at sub-domain)

  CHARACTER (LEN=ilstidn) , ALLOCATABLE :: &
    ybufsrt    (:) ,& ! character array with sorted reports read from file
    ybufloc    (:)    ! character array with local reports (at sub-domain)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_gps
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_gps'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

  DO ipe = 1 , num_compute
    nodrepn (ipe) = 0
    nodleni (ipe) = 0
    nodlenr (ipe) = 0
    nodleny (ipe) = 0
  ENDDO

! define number of int/real/char elements from common header part
! to be put into 'Xbufsrt'  (equal for all NetCDF observation input files)
! ------------------------------------------------------------------------
  noffscal (1,1)  =  3 + 18           ! 18 integer  elements + 3 length counters
  noffscal (2,1)  =  9                !  9 real     elements
  noffscal (3,1)  =  1                !  1 character element

! define total number of scalar int/real/char elements
! (i.e. elements without replication or with fixed replication)
! to be put into 'Xbufsrt'  (specific for the different file types)
! -----------------------------------------------------------------
  noffscal (1,2)  =  noffscal(1,1) + 3
  noffscal (2,2)  =  noffscal(2,1) + 6
  noffscal (3,2)  =  noffscal(3,1) + 0

!-------------------------------------------------------------------------------
! Section 1: Read and pre-process header information common to all obs types:
!            time / lat / lon / station altitude / obs type / station ID /
!            center of origin ...;
!            assign observations to grid points and to nodes (sub-domains)
!-------------------------------------------------------------------------------

  nrep = 0

  IF (my_cart_id == 0) THEN

    CALL obs_cdf_read_comhead ( min_sta , min_end , ilcdf                      &
                              , mrepsta , nrep , nreproc )
!   =========================
  ENDIF

! IF (num_compute > 1)  CALL global_values ( nrep, 1,'MAX',imp_integers        &
!                                                , icomm_cart, -1, yerr,ierr )
!                       ------------------
  IF (num_compute > 1) THEN
    CALL distribute_values (nrep    ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (imiss   ,1,0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (rmisschk,1,0,imp_reals   ,icomm_cart,ierr)
  ENDIF

  IF ((my_cart_id == 0) .AND. (nrep >= 1)) THEN

! assign the observations to the nodes (sub-domains), determine
! the local grid point indices, and produce a list of indices
! pointing to the reports sorted according to the nodes
! -------------------------------------------------------------

    ALLOCATE ( irsort (nreproc+1) , STAT=istat )
    ALLOCATE ( irnode (nreproc+1) , STAT=istat )

    CALL obs_assign_sort_node ( nrep, nreproc, irproc, iobstot, jobstot        &
                              , num_compute, i_subpos, nboundlines, my_cart_id &
                              , irnode, irsort, iobsloc, jobsloc )
!   =========================

! 'irnode' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get entries
!-------------------------------------------------------------------------------

! allocate arrays for elements without replication or with fixed replication
! --------------------------------------------------------------------------

    ALLOCATE ( rc_t   (1,nrep) , STAT=istat )
    ALLOCATE ( rc_p   (1,nrep) , STAT=istat )
    ALLOCATE ( nc_rh    (nrep) , STAT=istat )
    ALLOCATE ( nc_vdt   (nrep) , STAT=istat )
    ALLOCATE ( rc_el    (nrep) , STAT=istat )
    ALLOCATE ( rc_iwv (1,nrep) , STAT=istat )
    ALLOCATE ( rc_zwd (1,nrep) , STAT=istat )
    ALLOCATE ( rc_zpder (nrep) , STAT=istat )
    ALLOCATE ( rc_zpd   (nrep) , STAT=istat )
    ALLOCATE ( nc_flg   (nrep) , STAT=istat )


! get variable ID's in NetCDF file for elements with no or fixed replication
! --------------------------------------------------------------------------

                             status = nf90_inq_varid(ncid,'MPPP  ',varid_p   )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'MUUU  ',varid_rh  )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NGGTP ',varid_vdt )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'MDE   ',varid_el  )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NWLN  ',varid_iwv )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NCZWV ',varid_zwd )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NEERR ',varid_zpder)
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NADES ',varid_zpd )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'NQFGD ',varid_flg )
    IF(status == nf90_noerr) status = nf90_inq_varid(ncid,'MTN   ',varid_t   )

    IF (status /= nf90_noerr) THEN
      yerrmsl = 'STANDARD SURFACE HEADER INFO DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF

! get data for elements with no or fixed replication
! --------------------------------------------------

    status = nf90_get_var (ncid, varid_p    , rc_p    , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_rh   , nc_rh   , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_vdt  , nc_vdt  , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_el   , rc_el   , start=(/1,mrepsta/)     &
                                                      , count=(/1,nrep/))
    status = nf90_get_var (ncid, varid_iwv  , rc_iwv  , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_zwd  , rc_zwd  , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_zpder, rc_zpder, start=(/1,mrepsta/)     &
                                                      , count=(/1,nrep/))
    status = nf90_get_var (ncid, varid_zpd  , rc_zpd  , start=(/1,mrepsta/)     &
                                                      , count=(/1,nrep/))
    status = nf90_get_var (ncid, varid_flg  , nc_flg  , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_t    , rc_t    , (/mrepsta/), (/nrep/))


! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

    DO irps = 1 , nreproc
      irep  =  irsort(irps)
!     IF (ABS( rc_el(irep) ) < rmisschk) THEN
        IF (ABS( rc_el(irep)-90.0_wp ) > epsy) THEN
          rc_iwv (1,irep)  = rmiss
          rc_zwd (1,irep)  = rmiss
          rc_zpder (irep)  = rmiss
          rc_zpd   (irep)  = rmiss
        ENDIF
!     ENDIF
    ENDDO
    DEALLOCATE ( rc_el , STAT=istat )

!-------------------------------------------------------------------------------
! Section 3: Produce long buffer arrays 'Xbufsrt' in which the read reports are
!            sorted according to the nodes to which they will be distributed
!-------------------------------------------------------------------------------

    ALLOCATE ( iirlen (nreproc+1) , STAT=istat )
    ALLOCATE ( irrlen (nreproc+1) , STAT=istat )
    ALLOCATE ( iyrlen (nreproc+1) , STAT=istat )

! compute length (number of elements) for each individual TEMP report
! -------------------------------------------------------------------

! total number of scalar elements
    niscal = noffscal(1,2)
    nrscal = noffscal(2,2)
    nyscal = noffscal(3,2)

    DO irps = 1 , nreproc
      iirlen (irps)  =  niscal
      irrlen (irps)  =  nrscal
      iyrlen (irps)  =  nyscal
    ENDDO

! compute length of buffer arrays (i.e length of all reports together + 1)
! -------------------------------

    ibuflen (1)  =  1
    ibuflen (2)  =  1
    ibuflen (3)  =  1
    DO irps = 1 , nreproc
      ibuflen (1)  =  ibuflen(1) + iirlen (irps)
      ibuflen (2)  =  ibuflen(2) + irrlen (irps)
      ibuflen (3)  =  ibuflen(3) + iyrlen (irps)
    ENDDO

! allocate buffer arrays (only for my_cart_id == 0 here)
! ----------------------

    ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )

! fill the header part which is common to all obs types
! into the long buffer arrays 'Xbufsrt'
! -----------------------------------------------------

    CALL obs_cdf_buffer_comhead ( nreproc, irsort, iirlen, irrlen, iyrlen      &
                                , ibuflen, ibufsrt, rbufsrt, ybufsrt )
!   ===========================

! fill the file-specific scalar elements into the long buffer arrays 'Xbufsrt'
! ----------------------------------------------------------------------------

    iioffs = 0
    iroffs = 0

    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_rh   (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_vdt  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_flg  (irep)
      rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_t    (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_p    (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_iwv  (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 4)  =  rc_zwd  (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 5)  =  rc_zpder(irep)
      rbufsrt (iroffs+noffscal(2,1)+ 6)  =  rc_zpd  (irep)

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

! de-allocate reading arrays
    DEALLOCATE ( nc_rh    , STAT=istat )
    DEALLOCATE ( rc_t     , STAT=istat )
    DEALLOCATE ( rc_p     , STAT=istat )
    DEALLOCATE ( nc_vdt   , STAT=istat )
    DEALLOCATE ( rc_iwv   , STAT=istat )
    DEALLOCATE ( rc_zwd   , STAT=istat )
    DEALLOCATE ( rc_zpder , STAT=istat )
    DEALLOCATE ( rc_zpd   , STAT=istat )
    DEALLOCATE ( nc_flg   , STAT=istat )


!   PRINT '("buf0a ",8F10.2)' , (rbufsrt(istat),istat=1,8)
!   PRINT '("buf0b ",8F10.2)' , (rbufsrt(istat),istat=9,16)

!-------------------------------------------------------------------------------
! Section 4: Distribute the reports
!-------------------------------------------------------------------------------

! determine the number of int / real / char elements
! to be distributed to the different nodes

    DO irps = 1 , nreproc
      inode =  irnode(irps)
      nodrepn (inode+1) = nodrepn(inode+1) + 1
      nodleni (inode+1) = nodleni(inode+1) + iirlen(irps)
      nodlenr (inode+1) = nodlenr(inode+1) + irrlen(irps)
      nodleny (inode+1) = nodleny(inode+1) + iyrlen(irps)
    ENDDO
!   PRINT *,'distr3a ', nodrepn
!   PRINT *,'distr3b ', nodleni
!   PRINT *,'distr3c ', nodlenr
!   PRINT *,'distr3d ', nodleny

    DEALLOCATE ( iirlen  , STAT=istat )
    DEALLOCATE ( irrlen  , STAT=istat )
    DEALLOCATE ( iyrlen  , STAT=istat )
    DEALLOCATE ( irsort  , STAT=istat )
    DEALLOCATE ( irnode  , STAT=istat )

  ENDIF    ! (my_cart_id == 0)

  IF ((nrep >= 1) .AND. (num_compute > 1)) THEN

    iloclen(1)  =  MAXVAL( nodleni ) + 1
    iloclen(2)  =  MAXVAL( nodlenr ) + 1
    iloclen(3)  =  MAXVAL( nodleny ) + 1

    CALL distribute_values (ibuflen, 3, 0,imp_integers,icomm_cart,ierr)
    CALL distribute_values (iloclen, 3, 0,imp_integers,icomm_cart,ierr)

    IF (my_cart_id /= 0)  ALLOCATE ( ibufsrt (ibuflen(1)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( rbufsrt (ibuflen(2)) , STAT=istat )
    IF (my_cart_id /= 0)  ALLOCATE ( ybufsrt (ibuflen(3)) , STAT=istat )
                          ALLOCATE ( ibufloc (iloclen(1)) , STAT=istat )
                          ALLOCATE ( rbufloc (iloclen(2)) , STAT=istat )
                          ALLOCATE ( ybufloc (iloclen(3)) , STAT=istat )

    CALL obs_cdf_distrib_reports ( ibuflen, ibufsrt, rbufsrt, ybufsrt, ilstidn &
                                 , nodrepn, nodleni, nodlenr, nodleny          &
                                 , nrepl  , nlenli , nlenlr , nlenly           &
                                 , iloclen, ibufloc, rbufloc, ybufloc          &
                                 , num_compute , my_cart_id , icomm_cart       &
                                 , imp_integers, imp_reals )
!   ============================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

!-------------------------------------------------------------------------------
! Section 5: Store the local reports in the ODR
!-------------------------------------------------------------------------------

!   PRINT '("buf3a ",8F10.2)' , (rbufloc(istat),istat=1,8)
!   PRINT '("buf3b ",8F10.2)' , (rbufloc(istat),istat=9,16)
    CALL obs_cdf_store_gps    ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly   &
                              , noffscal,          ibufloc, rbufloc, ybufloc  &
                              , ngpnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_gps    ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly   &
                              , noffscal,          ibufsrt, rbufsrt, ybufsrt  &
                              , ngpnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_gps
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_gps


!===============================================================================
!+ Module procedure in "src_obs_cdfin_gps" for storing gps rep. in ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_store_gps ( kcdftyp, nrepl, nlenli, nlenlr, nlenly     &
                             , noffscal      , ibuf  , rbuf  , ybuf       &
                             , ngpnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_gps" stores the gps
!   reports in the internal ODR arrays.
!
! Method:
!   First, the header part common to all observation types is stored by calling
!   'obs_cdf_store_comhead', and then the rest of the header is added which
!   depends on the observation type.
!   Next, the observations are copied from the observation type
!   dependent buffers into standardised arrays. Subsequently, the quality flag
!   patterns are built (e.g. from blacklisting or gross error checking), and
!   the elements are stored in the ODR (internal arrays).
!   To evaluate the blacklist, first, a list of blacklisted vertical intervals
!   for each of the current reports is produced. Then, within the loop over
!   current reports, a blacklist flag 'kblk' is prepared for each variable
!   (according to the vertical level of the observation).
!   Finally, reports without observations are deleted in the ODR.
!
! Initial release: Christoph Schraff, 20.12.07
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)  ::  &
    kcdftyp          ,& ! type of NetCDF input file format (observation type):
                        !  = ncdf_gps_zenith : GPS
    nrepl            ,& ! number of            reports  \   at
    nlenli           ,& ! number of integer   elements   \  the
    nlenlr           ,& ! number of real      elements   /  local
    nlenly           ,& ! number of character elements  /   node
    noffscal (3,2)      ! number int/real/char - elements in common header part
                        !                      - scalar elements in total rep.

  INTEGER        , INTENT (IN)  ::  &
    ibuf (nlenli+1)     ! local buffer array containing all integer elements

  REAL (KIND=wp) , INTENT (IN)  ::  &
    rbuf (nlenlr+1)     ! local buffer array containing all real elements

  CHARACTER (LEN=10)  , INTENT (IN)  ::  &
    ybuf (nlenly+1)     ! local buffer array containing all character elements

  INTEGER        , INTENT (INOUT)  ::  &
    ngpnew           ,& ! number of new gps reports
    nexceed             ! number of reports in excess of array size

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER ::  &
    c1000 = 1000._wp ,& ! 1000.0
    c100r =  0.01_wp    !    0.01

! Local scalars:
! -------------

  LOGICAL        ::  &
    lblkw  (nrepl+1)    ! report blacklisted because sta-ID missing on whitelist
!   lrhtd            ,& ! report template contains: rel. hum. + dewpoint temp.
!   lrhqx            ,& ! report template contains: rel. hum. + mixing ratio
!   lrhtdqx          ,& ! report template contains: RH + dewpoint + mixing ratio
!   ltd              ,& ! report template contains: dewpoint temperature only
!   lqx              ,& ! report template contains: mixing ratio only
!   lexitd           ,& ! report template includes dewpoint temperature
!   lexirh           ,& ! report template includes relative humidity
!   lexiqx           ,& ! report template includes mixing ratio
!   lrej             ,& ! reject upper-air report (no valid pressure level)
!   liwv

  INTEGER        ::  &
    iioffs (nrepl+1) ,& ! offset of the integer   elements  \  for each
    iroffs (nrepl+1) ,& ! offset of the real      elements   > local
    iyoffs (nrepl+1) ,& ! offset of the character elements  /  report
    jobtyp (nrepl+1) ,& ! observation type (input array for local blacklist)
    jcdtyp (nrepl+1) ,& ! obs  code   type (input array for local whitelist)
    kblk  (nrepl+1,4),& ! local blacklist (0: ok , 1: blacklisted, for each pro-
                        !                 cessed report and variable separately)
    nflgx  (nrepl+1) ,& ! total flag for an observation of a variable
    icmaa  (nrepl+1) ,& ! pointer index for 'cma' / diagnostic arrays
    ioberr (nrepl+1) ,& ! active status
    icma             ,& ! pointer index for 'cma' / diagnostic arrays
    ngpob            ,& ! target ODR single-level report index
    nicomh , nrcomh  ,& ! number of common header integer / real  elements
    niscal , nrscal  ,& ! number of scalar integer / real  elements
    nrepgp           ,& ! number of reports which can be stored in the ODR now
    nexcegp          ,& ! no. of local sing-lev. reports in excess of array size
    irpl             ,& ! loop index over local reports
    intv             ,& ! loop index over blacklisted vertical intervals
    kobtyp , kcdtyp  ,& ! observation type / observation code type
    istalt           ,& ! station altitude
!   npfi             ,& ! pressure level [hPa] of reported height
    ilverr              ! nearest standard error level below observation level

  INTEGER        ::  &
    iob    , job     ,& ! indices of grid point assigned to observation
    iactr            ,& ! stepping for event counters
    ilen             ,& ! length of control message
    iactx (5)        ,& ! indicates presence of active  data (for variables)
    ipasx (5)        ,& ! indicates presence of passive data (for variables)
    nzaexi           ,& ! -1: only passive data ; 0: no data ; 1: active data
    ilstid_chk       ,& ! length of (part of) station ID that can be checked
    istat            ,& ! error indicators
    kcentsub (nrepl+1)

  REAL (KIND=wp) ::  &
    zobhr (nrepl+1)  ,& ! obs  time (input array for local whitelist)
    zpref (nrepl+1)  ,& ! 'reference pressure' for evaluating the blacklist
    zpuse (nrepl+1)  ,& ! reported pressure
    zzuse (nrepl+1)  ,& ! reported height
    zoberr(nrepl+1)  ,& ! observation error
    fisd  (nrepl+1)  ,& ! height diff. betw. model orography and sta. altitude
!   zrlat (nrepl+1,1),& ! latitude  of (processed) observation level
!   zrlon (nrepl+1,1),& ! longitude of (processed) observation level
!   ztmean           ,& ! mean temperature (betw. station height and msl)
!   zpred            ,& ! reported station pressure reduced to mean sea level
!   zplim            ,& ! threshold for reported station pressure check
!   zsurf            ,& ! height of model orography
!   zporo            ,& ! (approx.) model surface pressure
    zadder           ,& ! additional pressure obs error due to extrapolation
    fisdtt           ,& ! scaled extrapolation distance for 2-m temperature val
    fisdrh           ,& ! scaled extrapolation distance for 2-m humidity values
    fisdzz           ,& ! scaled extrapolation distance for surface press. val.
    fisdnu           ,& ! scaled extrapolation distance for full values and
                        !   increments of surface pressure (for nudging only)
    fisdpd           ,& ! scaled extrapolation distance for zenith path delay
    zzkeml           ,& ! height of the lowest main model level
    zlop             ,& ! LOG( pressure )
    fiperr           ,& ! interpolation weight factor for error level 'ilverr'
    zbias               ! IWV bias correction value

  CHARACTER (LEN=25)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid_blk)  :: &
    ystidl  (nrepl+1)   ! observation type (input array for local blacklist)

! Local arrays:
! ------------

  INTEGER        , ALLOCATABLE ::  &
!   ivdt         (:) ,& ! time period or displacement
    nflgg        (:)    ! quality flag for GNSS data

  REAL (KIND=wp) , ALLOCATABLE ::  &
    zrhc         (:) ,& ! relative humidity (model compatible, bias corrected)
    zrhw         (:) ,& ! relative humidity over water
!   ztd          (:) ,& ! dew point temperature
!   zqx          (:) ,& ! mixing ratio
    zpp          (:) ,& ! pressure
    ztt          (:) ,& ! temperature
!   zzz          (:) ,& ! height
!   zttq         (:) ,& ! precision of temperature observation
!   zzt          (:) ,& ! height of temperature sensor above ground
!   zpmsl        (:) ,& ! mean sea level pressure
!   zdpdt        (:) ,& ! 3-hour pressure change
!   zpfi         (:) ,& ! pressure level of reported height
!   zrhw1        (:) ,& ! relative humidity over water (converted from dewpoint)
!   zrhw2        (:) ,& ! relative humidity over water (conv. from mixing ratio)
!   zrh1         (:) ,& ! model compatible rel humid.  (converted from dewpoint)
!   zrh2         (:) ,& ! model compatible rel humid.  (conv. from mixing ratio)
!   zqvw         (:) ,& ! specific humidity over water (conv. from mixing ratio)
!   zqv          (:) ,& ! model compatible spec. humid (conv. from mixing ratio)
    iwv          (:) ,& ! integrated water vapour
    zwd          (:) ,& ! zenith wet delay
    zpder        (:) ,& ! error in zenith path delay
    zpd          (:)    ! zenith path delay

!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_store_gps
!-------------------------------------------------------------------------------

 yroutine = 'obs_cdf_store_gps'

! determine the offsets of the different reports in the long (local) array
! ------------------------------------------------------------------------

  iioffs (1) = 0
  iroffs (1) = 0
  iyoffs (1) = 0
  DO irpl = 2 , nrepl
    iioffs (irpl)  =  iioffs(irpl-1) + ibuf(iioffs(irpl-1)+1)
    iroffs (irpl)  =  iroffs(irpl-1) + ibuf(iioffs(irpl-1)+2)
    iyoffs (irpl)  =  iyoffs(irpl-1) + ibuf(iioffs(irpl-1)+3)
!   PRINT *,'bufoff ', irpl, iioffs(irpl), iroffs(irpl), iyoffs(irpl)
  ENDDO
! PRINT '("buf4a ",8F10.2)' , (rbuf(istat),istat=1,8)
! PRINT '("buf4b ",8F10.2)' , (rbuf(istat),istat=9,16)

! total number of common header elements
! --------------------------------------
  nicomh = noffscal(1,1)
  nrcomh = noffscal(2,1)

! total number of scalar elements
! -------------------------------
  niscal = noffscal(1,2)
  nrscal = noffscal(2,2)

!-------------------------------------------------------------------------------
! Section 1: Store report header
!-------------------------------------------------------------------------------

! store header part, which is common to all observation types, in ODR
! -------------------------------------------------------------------

  CALL obs_cdf_store_comhead ( 3, ntotgp, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, ngpnew, nexcegp )
! ==========================

! update number of reports which can be stored in the ODR now
! (the following line is equivalent to:  nrepgp = ngpnew)
  nrepgp = nrepl - nexcegp

! update counter for insufficient ODR size, for caution messages
  nexceed = nexceed + nexcegp

! store header part, which is specific to NetCDF
! observation input file type, in single-level ODR:
! ------------------------------------------------

!    iofstyp = 1
!    DO irpl = 1 , nrepgp
!      ngpob = ntotgp + irpl
!      mogphd (ngpob,nhstyp) = imdi
!      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ iofstyp) /= imiss)                  &
!        mogphd (ngpob,nhstyp) = ibuf (iioffs(irpl)+noffscal(1,1)+ iofstyp)
!    ENDDO

! get diagnostic array position (--> i_cma)
! -----------------------------------------

  DO irpl = 1 , nrepgp
    ngpob  = ntotgp + irpl
    kobtyp = mogphd(ngpob,nhobtp)
    kcdtyp = mogphd(ngpob,nhcode)
    icmaa (irpl)  =  i_cma ( kobtyp , kcdtyp )
!                    =====
  ENDDO

! get list of blacklisted vertical intervals for each report
! ----------------------------------------------------------

  ilstid_chk  =  MIN( ilstid, ilstid_blk )
  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    jobtyp (irpl)  =  mogphd(ngpob,nhobtp)
    jcdtyp (irpl)  =  mogphd(ngpob,nhcode)
    zobhr  (irpl)  =  ogphed(ngpob,nhtime)
    ystidl (irpl)  =  ' '
    ystidl (irpl)  =  yogphd(ngpob) (1:ilstid_chk)
    kcentsub (irpl)= mogphd(ngpob,nhcent)
  ENDDO

  IF (nrepgp >= 1) THEN
    ALLOCATE ( blk_loc (nrepgp) , STAT=istat )

    CALL obs_cdf_blacklist_local ( nrepgp , jobtyp , ilstid_chk , ystidl )
!   ============================

! determine which reports are missing on the whitelists
! -----------------------------------------------------

    CALL obs_cdf_whitelist_local ( nrepgp , jobtyp , jcdtyp , zobhr            &
                                 , ilstid_chk , ystidl , lblkw )
!   ============================

  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Report body: Put the observations from the
!            buffers into standardised arrays
!-------------------------------------------------------------------------------

! allocate standardised arrays
! ----------------------------

  IF (nrepgp>= 1) THEN
    ALLOCATE ( zrhw    (nrepgp) , STAT=istat )
    ALLOCATE ( zrhc    (nrepgp) , STAT=istat )
    ALLOCATE ( zpp     (nrepgp) , STAT=istat )
    ALLOCATE ( ztt     (nrepgp) , STAT=istat )
!   ALLOCATE ( ivdt    (nrepgp) , STAT=istat )
    ALLOCATE ( iwv     (nrepgp) , STAT=istat )
    ALLOCATE ( zwd     (nrepgp) , STAT=istat )
    ALLOCATE ( zpder   (nrepgp) , STAT=istat )
    ALLOCATE ( zpd     (nrepgp) , STAT=istat )
    ALLOCATE ( nflgg   (nrepgp) , STAT=istat )
  ENDIF

! initialise standardised arrays
! ------------------------------

  DO irpl = 1 , nrepgp
    zrhw   (irpl) = rmdi
    zrhc   (irpl) = rmdi

    zpp    (irpl) = rmdi
    ztt    (irpl) = rmdi
!   ivdt   (irpl) = imdi
    iwv    (irpl) = rmdi
    zwd    (irpl) = rmdi
    zpder  (irpl) = rmdi
    zpd    (irpl) = rmdi
    nflgg  (irpl) = imdi
  ENDDO

! fill standardised arrays
! ------------------------

  IF (nrepgp >= 1) THEN
    DO irpl = 1 , nrepgp
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        ztt   (irpl) =      rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (     ibuf(iioffs(irpl)+nicomh+ 1)  /= imiss   ) THEN
        zrhw  (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 1) , wp ) *c100r
      ENDIF
!     IF (     ibuf(iioffs(irpl)+nicomh+ 2)  /= imiss   )                      &
!         ivdt  (irpl) =    ibuf (iioffs(irpl)+nicomh+ 2)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        zpp  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 3) ) < rmisschk)                      &
        iwv  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 3)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 4) ) < rmisschk)                      &
        zwd  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 4)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 5) ) < rmisschk)                      &
        zpder(irpl) =       rbuf (iroffs(irpl)+nrcomh+ 5)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 6) ) < rmisschk)                      &
        zpd  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 6)
      IF (   ( ibuf(iioffs(irpl)+nicomh+ 3) ) /= imiss  )                      &
        nflgg(irpl) =       ibuf (iioffs(irpl)+nicomh+ 3)
    ENDDO
  ENDIF


!-------------------------------------------------------------------------------
! Section 3: Process pressure and geopotential, derive local blacklist, and
!            fill pressure and height into ODR arrays
!-------------------------------------------------------------------------------

! ---------------
! surface reports: get surface pressure, reference height level, and pressure code
! ---------------  and station height

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
! note: station height is always well-defined for GPS reports
!       (see routine 'obs_cdf_read_comhead')
    istalt        = NINT( ogphed(ngpob,nhalt) )
! use station pressure normally
    zpuse  (irpl) = zpp(irpl)
    zzuse  (irpl) = ogphed(ngpob,nhalt)

! get 'station height minus model orography' + land/sea obs indicator, used
! further below to decide whether different variables are active or passive
!   fisd   (irpl) = ogphed(ngpob,nhsurf) - ogphed(ngpob,nhalt)
    fisd   (irpl) = ogphed(ngpob,nhalt)  - ogphed(ngpob,nhsurf)
! get 'reference pressure' for evaluating the blacklist
! ( = zpuse , or ps_mod if zpuse is an extrapolated value or undefined)
    zpref  (irpl) = zpuse(irpl)
    IF (zpref(irpl) < rmdich) THEN
      iob    = mogphd(ngpob,nhio)
      job    = mogphd(ngpob,nhjo)
      zpref(irpl) = r_ps(iob,job)
    ENDIF
  ENDDO

! auxilliary: produce 2-dim obs blacklist field (report, variable v,p,T,q)
! ------------------------------------------------------------------------

  DO irpl = 1 , nrepgp
    kblk (irpl,1) = 0
    kblk (irpl,2) = 0
    kblk (irpl,3) = 0
    kblk (irpl,4) = 0
    DO intv = 1 , blk_loc(irpl)% ndim
! if pressure not defined: variables with a non-zero entry in the original
! global blacklist are always flagged here
      IF (zpuse(irpl) <= rmdich) THEN
        kblk (irpl,blk_loc(irpl)%kvar(intv)) = 1
      ELSE
! if pressure defined: variables are flagged only if pressure is in the
! interval given in the original blacklist
        IF (      (zpref(irpl) <= blk_loc(irpl)%plow(intv) +epsy)              &
            .AND. (zpref(irpl) >= blk_loc(irpl)%pup (intv) -epsy))             &
          kblk (irpl,blk_loc(irpl)%kvar(intv)) = 1
      ENDIF
    ENDDO
  ENDDO
  DO irpl = 1 , nrepgp
! if missing on whitelist: all variables are flagged
    IF (lblkw(irpl)) THEN
      kblk (irpl,1) = 1
      kblk (irpl,2) = 1
      kblk (irpl,3) = 1
      kblk (irpl,4) = 1
    ENDIF
  ENDDO

  IF (nrepgp >= 1) THEN
    DEALLOCATE ( blk_loc , STAT=istat )
  ENDIF

! ---------------------------------------
! pressure observation error and flagging
! ---------------------------------------

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    kobtyp = mogphd(ngpob,nhobtp)
    icma   = icmaa (irpl)
    iactr  = 0
    IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1
    nflgx  (irpl)  =   0
    zoberr (irpl)  =  c0
    ioberr (irpl)  =   1
    IF (mogphd(ngpob,nhpass) == -1)  THEN
      ioberr (irpl)  =  0
      fisdnu         =  rmdi
 !if pressure missing (no flag in 'nflgx' set)
    ELSEIF (zpuse(irpl) < rmdich) THEN
      ioberr (irpl)  =  0
      fisdnu         =  rmdi
      IF (nzex(kobtyp) >= 1)  neventd (nepmis,icma) =                          &
                              neventd (nepmis,icma) + iactr
    ELSE
! if pressure blacklisted
      IF (kblk(irpl,2) == 1) THEN
        IF (ioberr(irpl) == 1)  neventd (nepflg,icma) =                        &
                                neventd (nepflg,icma) + iactr
        nflgx  (irpl)  =  IBSET( nflgx(irpl), nvfbps(2) )
        ioberr (irpl)  =  0
      ENDIF

! rejection of observed pressure depending on assumed extrapolation errors
! (fisdzz shall account only for errors from extrapolation of full (model or
!  observed) values, not for interpolation errors)
! --> If a station reports station pressure and lies above the model orography
!     (or more precisely the height of the lowest main model level 'zzkeml')
!     then the observation increment will be computed at the reported pressure
!     level 'zzuse' by interpolation of model pressure values; and only in the
!     case of nudging, this increment has to be extrapolated further to the
!     lowest model level.
! --> If (zzuse - zzkeml) > doromx(2) then the increment is not used any more
!     in the nudging, but still used in the LETKF. Then, the status bit in word
!     'nbgerr' is set to 1 (active), but at the same time, the invalid height
!     range flag is set in 'nbgflg'.
! (for surface reports, 'zzuse' and station height are always well-defined)
      iob    = mogphd(ngpob,nhio)
      job    = mogphd(ngpob,nhjo)
      zzkeml = c05* (r_hhl(iob,job,ke+1)+r_hhl(iob,job,ke))
      IF (zzuse(irpl) >= zzkeml-epsy) THEN
        ! --> only extrapolate obs-value from station height to reported p-level
        !     (model-values are interpolated to reported p-level)
        fisdzz = ABS( ogphed(ngpob,nhalt) - zzuse(irpl) )
        fisdnu = fisdzz + fdoro(2) *(zzuse(irpl) - zzkeml)
      ELSEIF (zzuse(irpl) >= ogphed(ngpob,nhalt)-epsy) THEN
        ! i.e. zzkeml > zzuse >= sta_alt
        ! --> extrapolate mod-val from model orography down to reported p-level
        !            plus obs-val from station height  up   to reported p-level
        ! --> ABS(zzkeml-zzuse)+ABS(zzuse-sta_alt) = zzkeml-sta_alt
        fisdzz = zzkeml - ogphed(ngpob,nhalt)
        fisdnu = fisdzz
      ELSE
        ! i.e. zzkeml and sta_alt > zzuse
        ! --> extrapolate mod-val from model orography down to reported p-level
        !            plus obs-val from station height  down to reported p-level
        ! ==> errors from extrapolation of model and observed value in the same
        !     direction may (or may not) be similar (and then cancel each other
        !     at the computation of the observation increment)
        ! --> assume that the total extrapolation error (to compute increment)
        !     is close to the maximum of the 2 individual extrapolation errors
        fisdzz = MAX( ogphed(ngpob,nhalt) - zzuse(irpl)                        &
                    , zzkeml              - zzuse(irpl) )
        fisdnu = fisdzz
      ENDIF
! additional observation error assigned due to extrapolation
      zadder = 0.04_wp * fisdzz
! if station height or its diff. to model orography too large for surf. pressure
      IF ((ogphed(ngpob,nhalt) > altopsu(2))      &
          .OR. (fisdzz > doromx(2))) THEN
        IF (ioberr(irpl) == 1)  neventd (nepalt,icma) =                        &
                                neventd (nepalt,icma) + iactr
        nflgx  (irpl)  =  IBSET( nflgx(irpl), nvfbps(4) )
        ioberr (irpl)  =  0
      ELSEIF (fisdnu > doromx(2)) THEN
! use this obs for LETKF, but reject it in the nudging
! due to additional error from extrapolation of (obs) increment
        nflgx  (irpl)  =  IBSET( nflgx(irpl), nvfbps(4) )
      ENDIF

! observation error enhancement for surface pressure
!     zoberr (irpl) = obs_error ( nvrpoe(3) ) + zadder
! convey height error to pressure error (rough approximation)
      zoberr (irpl) = zadder * zpuse(irpl) / (r_d * tmelt)
    ENDIF

! fill ODR (pressure, height)
! --------
!   no input available for override, data base quality, and human monitoring flag
!     (override flag is set in obs header, not here)
!   nflgx(irpl)  =  0
!   CALL MVBITS ( nppfwr+izzfwr, nf1bps(4), nvfboc(1), nflgx(irpl), nvfbps(1) )
!   CALL MVBITS ( nppfwr+izzfwr, nf1bps(5), nvfboc(2), nflgx(irpl), nvfbps(2) )
!   CALL MVBITS ( nppfwr+izzfwr, nf1bps(3), nvfboc(3), nflgx(irpl), nvfbps(3) )
!     (for upper-air reports, 'zpuse' is always well-defined (i.e. not missing value))
    ogpbdy (ngpob,nbgp  )  =  zpuse (irpl)
!     (for surface reports, 'zzuse' and station height are always well-defined)
    ogpbdy (ngpob,nbgz  )  =  zzuse (irpl)
    ogpbdy (ngpob,nbgzer)  =  zoberr(irpl)
    ogpbdy (ngpob,nbgviz)  =  fisdnu
    mogpbd (ngpob,nbgerr)  =  insert( 0, ioberr(irpl), nvrz )
    mogpbd (ngpob,nbgflg)  =  insert( 0, nflgx(irpl), nvfzbp )
    mogpbd (ngpob,nbglid)  =  IBSET ( 0, nvlidp(7) )

    mogpbd (ngpob,nbgqcf)  =  0

  ENDDO

!-------------------------------------------------------------------------------
! Section 4: Determine and store the main elements in the single-level ODR
!            (temperature, humidity, wind)
!-------------------------------------------------------------------------------

! -----------
! temperature
! -----------

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    IF ((mogphd(ngpob,nhpass) /= -1)) THEN
      kobtyp = mogphd(ngpob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1
      nflgx  (irpl) =  0
      zoberr (irpl) = c0
      ioberr (irpl) =  1

! if temperature missing or below absolute minimum (no flag in 'nflgx' set)
      IF ((ztt(irpl) < rmdich) .OR. (ztt(irpl) < c1)) THEN
        ztt    (irpl) = rmdi
        ioberr (irpl) = 0
        IF (ntex(kobtyp) >= 1)  neventd (netmis,icma) =                        &
                                neventd (netmis,icma) + iactr
      ELSE
! if height of sensor above ground is not close to 2 m (ok: 1.5m <= zzt <= 2.5m)
!       IF ((zzt(irpl) > rmdich) .AND. (ABS( zzt(irpl)-c2 ) > c05+epsy)) THEN
!         IF (ioberr(irpl) == 1)  neventd (netflg,icma) =                      &
!                                 neventd (netflg,icma) + iactr
!         nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(5) )
!         ioberr (irpl) = 0
!       ENDIF
! if temperature blacklisted
        IF (kblk(irpl,3) == 1) THEN
          IF (ioberr(irpl) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(2) )
          ioberr (irpl) = 0
        ENDIF
! if station height or its diff to model orography too large for 2-m temperature
        fisdtt =  ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd(irpl) ))   &
                  * fisd(irpl)
        IF (     (ogphed(ngpob,nhalt) > altopsu(3))  &
            .OR. (fisdtt > doromx(3))) THEN
          IF (ioberr(irpl) == 1)  neventd (netalt,icma) =                      &
                                  neventd (netalt,icma) + iactr
          nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(4) )
          ioberr (irpl) = 0
        ENDIF
! if temperature gross error (note: zpref is always defined)
        IF (     (ztt(irpl) > tmelt+60._wp)                                    &
            .OR.((ztt(irpl) > tmelt+20._wp) .AND. (zpref(irpl) < 70000._wp))   &
            .OR.((ztt(irpl) > tmelt+ 5._wp) .AND. (zpref(irpl) < 50000._wp))   &
            .OR.((ztt(irpl) > tmelt- 5._wp) .AND. (zpref(irpl) < 40000._wp))   &
            .OR. (ztt(irpl) < tmelt-90._wp)) THEN
          IF (ioberr(irpl) == 1)  neventd (netext,icma) =                      &
                                  neventd (netext,icma) + iactr
          nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(3) )
          ioberr (irpl) = 0
        ENDIF
      ENDIF

! observation error for 2-m temperature
      zoberr (irpl) = c1 + 0.003_wp *ABS( fisd(irpl) )

! fill ODR (temperature)
      ogpbdy (ngpob,nbgt  ) = ztt   (irpl)
      ogpbdy (ngpob,nbgter) = zoberr(irpl)
      mogpbd (ngpob,nbgerr) = insert( mogpbd(ngpob,nbgerr), ioberr(irpl), nvrt )
      mogpbd (ngpob,nbgflg) = insert( mogpbd(ngpob,nbgflg), nflgx(irpl), nvftbp)
    ENDIF
  ENDDO

! --------
! humidity
! --------

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
!   PRINT *, 'zrh0a ',irpl, mosghd(nsgob,nhpass), zrhw(irpl)
    IF (mogphd(ngpob,nhpass) /= -1) THEN
!     PRINT *, 'zqx0  ',irpl, ystidl(irpl), zrhw(irpl), zqx(irpl), mosghd(nsgob,nhhrmn)
      kobtyp = mogphd(ngpob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1
      nflgx  (irpl) =  0
      zoberr (irpl) = c0
      ioberr (irpl) =  1

! reported humidity is fully discarded if negative or exactly zero
!   (sometimes reported zero values denote undefined values)
      IF (zrhw(irpl) <= 1.E-10_wp) THEN
        zrhw   (irpl)  =  rmdi
        ioberr (irpl)  =  0   ! is set further below if humidity is missing
        IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =       & ! also set
                                neventd (neqflg,icma) + iactr   ! below
      ENDIF
!     IF (ztd(irpl) <= 1.E-10_wp)  ztd (irpl)  =  rmdi
!     IF (zqx(irpl) <= 1.E-10_wp)  zqx (irpl)  =  rmdi

! as long as reported humidity of any kind is converted into relative humidity
! in order to be processed further (as is currently the case for nudging), then:
!   mixing ratio without temperature or pressure, and
!   dewpoint temperature without temperature      are completely discarded
!   (since these quantities are required for conversion into relative humidity)
!     IF ((zqx(irpl) > rmdich) .AND. (     (ztt  (irpl) < rmdich)              &
!                                     .OR. (zpuse(irpl) < rmdich))) THEN
!       zqx (irpl)  =  rmdi
!       IF ((zrhw(irpl) < rmdich) .AND. (ztd(irpl) < rmdich)) THEN
!         IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
!                                 neventd (neqflg,icma) + iactr
!         nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!         ioberr (irpl) = 0
!       ENDIF
!     ENDIF
!     IF ((ztd(irpl) > rmdich) .AND.       (ztt  (irpl) < rmdich) ) THEN
!       ztd (irpl)  =  rmdi
!       IF ((zrhw(irpl) < rmdich) .AND. (zqx(irpl) < rmdich)) THEN
!         IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
!                                 neventd (neqflg,icma) + iactr
!         nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!         ioberr (irpl) = 0
!       ENDIF
!     ENDIF

! if (relative or dewpoint) humidity missing (no flag in 'nflgx' set)
!     IF (    ((lrhtdqx) .AND. (MAX(zrhw(irpl),ztd(irpl),zqx(irpl)) < rmdich)) &
!        .OR. ((lrhtd  ) .AND. (MAX(zrhw(irpl),ztd(irpl)          ) < rmdich)) &
!        .OR. ((lrhqx  ) .AND. (MAX(zrhw(irpl)          ,zqx(irpl)) < rmdich)) &
!        .OR. ((ltd    ) .AND. (    ztd (irpl) < rmdich   ))) THEN
      IF (zrhw(irpl) < rmdich) THEN
        IF ((ioberr(irpl) == 1) .AND. (ntdex(kobtyp) >= 1))                    &
          neventd (neqmis,icma) = neventd(neqmis,icma) + iactr
!       zrhw   (irpl) = rmdi
!       ztd    (irpl) = rmdi
!       zqx    (irpl) = rmdi
        zrhc   (irpl) = rmdi
        ioberr (irpl) = 0
      ELSE
! if humidity blacklisted
        IF (kblk(irpl,4) == 1) THEN
          IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(2) )
          ioberr (irpl) = 0
        ENDIF
! if station height or its differ. to model orography too large for 2-m humidity
        fisdrh =  ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd(irpl) ))   &
                  * fisd(irpl)
        IF ((ogphed(ngpob,nhalt) > altopsu(4)) .OR. (fisdrh > doromx(4))) THEN
          IF (ioberr(irpl) == 1)  neventd (neqalt,icma) =                      &
                                  neventd (neqalt,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(4) )
          ioberr (irpl) = 0
        ENDIF
! derived rel. humidity is set passive if temperature passive and 'zrhw' missing
!       IF (      (.NOT. BTEST( mogpbd(ngpob,nbserr), nvrt ))                  &
!           .AND. (zrhw(irpl) < rmdich)                                        &
!           .AND. ((ztd(irpl) > rmdich) .OR. (zqx(irpl) > rmdich))) THEN
!         IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
!                                 neventd (neqflg,icma) + iactr
!         nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!         ioberr (irpl) = 0
!       ENDIF
! if dewpoint temperature gross error, or greater than temperature+2K
!      IF (lexitd) THEN    <-- not required (if 'ztd' is always initialised)
!       IF (      (ztd(irpl) > rmdich)                                         &
!           .AND. (      (ztd(irpl) < tmelt -150._wp)                          &
!                  .OR. ((ztd(irpl) < tmelt - 90._wp)                          &
!                  .OR.  (ztd(irpl) > tmelt + 40._wp)                          &
!                  .OR.  (ztd(irpl) > ztt(irpl)+c2)))) THEN
!         IF (MAX( zrhw(irpl),zqx(irpl) ) < rmdich) THEN
!           IF (ioberr(irpl) == 1)  neventd (neqlow,icma) =                    &
!                                   neventd (neqlow,icma) + iactr
!           nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
!           ioberr (irpl) = 0
! (minimum ztd must be > max( b4w,b4i ) to avoid division by zero)
!           ztd  (irpl)  =  MAX( ztd(irpl), b4w+c1 )
!           ztd  (irpl)  =  MIN( MIN( ztd(irpl), ztt(irpl)+c2 )                &
!                              , tmelt + 40._wp )
!         ELSE
! (if 'ztd' has gross error, use 'zrhw' or 'zqx' unless both are missing)
!           ztd  (irpl)  =  rmdi
!         ENDIF
!       ENDIF
! if mixing ratio gross error  , i.e. not in (1E-10 < zqx(w) <= 0.03)
!       IF (      (zqx(irpl) > rmdich)                                         &
!           .AND. (     (zqx(irpl) > 0.03_wp)                                  &
!                  .OR. (zqx(irpl) <= 1.E-10_wp))) THEN
!         IF (MAX( zrhw(irpl),ztd(irpl) ) < rmdich) THEN
!           IF (ioberr(irpl) == 1) THEN
!             IF (zqx(irpl) < 1.E-10_wp) neventd (neqlow,icma) =               &
!                                        neventd (neqlow,icma) + iactr
!             IF (zqx(irpl) > 0.03_wp)   neventd (neqbig,icma) =               &
!                                        neventd (neqbig,icma) + iactr
!           ENDIF
!           nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
!           ioberr (irpl) = 0
!           zqx (irpl)  =  MAX( MIN( zqx(irpl) , 0.03_wp ) , 1.E-10_wp )
!         ELSE
!    (if 'zqx' has gross error, use 'zrhw' or 'ztd' unless both are missing)
!           zqx (irpl)  =  rmdi
!         ENDIF
!       ENDIF
! if relative humidity gross error  , i.e. not in (epsy < zrhw <= 102%)
        IF (      (zrhw(irpl) > rmdich)                                        &
            .AND. (     (zrhw(irpl) > rtshlm)                                  &
                   .OR. (zrhw(irpl) <= epsy))) THEN
!         IF (MAX( ztd(irpl),zqx(irpl) ) < rmdich) THEN
            IF (ioberr(irpl) == 1) THEN
              IF (zrhw(irpl) <= epsy)   neventd (neqlow,icma) =                &
                                        neventd (neqlow,icma) + iactr
              IF (zrhw(irpl) > rtshlm)  neventd (neqbig,icma) =                &
                                        neventd (neqbig,icma) + iactr
            ENDIF
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
            ioberr (irpl) = 0
!           zrhw  (irpl)  =  MAX( MIN( zrhw(irpl) , c1 ) , c05*epsy )
            zrhw  (irpl)  =  MAX( MIN( zrhw(irpl) , c1 ) , c0 )
!         ELSE
!    (if 'zrh' has gross error, use 'ztd' or 'zqx' unless both are missing)
!           zrhw  (irpl)  =  rmdi
!         ENDIF
        ENDIF
      ENDIF
!     PRINT *, 'zqx2  ',irpl, ystidl(irpl), zqx(irpl), zrhw(irpl)
    ENDIF
    zrhc  (irpl)  =  zrhw (irpl)
!   PRINT *, 'zrh1 ',lexirh,lexitd,lexiqx,irpl, zrhw(irpl)
  ENDDO

! compute model compatible relative humidity
! ------------------------------------------
! at this point, 'zrhw', 'ztd', 'zqx' and 'zrhc' should have either reasonable
! or missing values, which should prevent lethal crashes in the following calls

  IF (nrepgp >= 1)  CALL obs_rhw2rh ( madj_hum, nrepgp, ztt, zrhw, rmdi        &
                                    , b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhc )
                  ! ===============

! IF ((lexitd) .AND. (nrepgp > 0)) THEN
!   ALLOCATE ( zrhw1 (nrepgp) , STAT=istat )
!   ALLOCATE ( zrh1  (nrepgp) , STAT=istat )
!
!   CALL obs_td2rh  ( madj_hum, nrepgp, ztt , ztd                              &
!                   , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhw1, zrh1 )
!   ==============
! ENDIF

! IF ((lexiqx) .AND. (nrepgp > 0)) THEN
!   ALLOCATE ( zrhw2 (nrepgp) , STAT=istat )
!   ALLOCATE ( zrh2  (nrepgp) , STAT=istat )
!   ALLOCATE ( zqvw  (nrepgp) , STAT=istat )
!   ALLOCATE ( zqv   (nrepgp) , STAT=istat )
!   zqvw = rmdi

!   CALL obs_qx2rh ( madj_hum, nrepgp, ztt, zpuse, zqx , zqvw                  &
!                  , rmdi, b1,b2w,b2i, b3,b4w,b4i, rdv, tmelt, zqv, zrhw2, zrh2)
!   ==============
! ENDIF

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    IF (mogphd(ngpob,nhpass) /= -1) THEN
!     PRINT *, 'zqx7b  ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1

! check consistency between reported and derived relative humidity values
! -----------------------------------------------------------------------
!     IF ((lexitd) .AND. (lexirh)) THEN
! if dewpoint-derived rel. humidity and reported rel. humidity differ by > 4%
!       IF ((zrh1(irpl) > rmdich) .AND. (zrhc(irpl) > rmdich)) THEN
!         IF (ABS( zrh1(irpl) - zrhc(irpl) ) > 0.04_wp) THEN
!           IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
!                                   neventd (neqflg,icma) + iactr
!           nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!           ioberr (irpl) = 0
!         ENDIF
!       ENDIF
!     ENDIF
!     IF ((lexiqx) .AND. (lexirh)) THEN
! if mixing-ratio-derived rel. humidity and reported rel. humid. differ by >4%
!       IF ((zrh2(irpl) > rmdich) .AND. (zrhc(irpl) > rmdich)) THEN
!         IF (ABS( zrh2(irpl) - zrhc(irpl) ) > 0.04_wp) THEN
!           IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
!                                   neventd (neqflg,icma) + iactr
!           nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!           ioberr (irpl) = 0
!         ENDIF
!       ENDIF
!     ENDIF
!     IF ((lexiqx) .AND. (lexitd)) THEN
! if mixing-ratio-derived and dewpoint-derived rel. humidity differ by >4%
!       IF ((zrh2(irpl) > rmdich) .AND. (zrh1(irpl) > rmdich)) THEN
!         IF (ABS( zrh2(irpl) - zrh1(irpl) ) > 0.04_wp) THEN
!           IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
!                                   neventd (neqflg,icma) + iactr
!           nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
!           ioberr (irpl) = 0
!         ENDIF
!       ENDIF
!     ENDIF
!     PRINT *, 'zqx8   ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)

! if relative humidity is not reported ('zrhc' = missing value)
! then use the dewpoint- or mixing-ratio-derived value 'zrh1' resp. 'zrh2'
! ------------------------------------------------------------------------
!     IF (lexitd) THEN
!      PRINT *, 'zrhtd ',zrh1(irpl)
!       IF ((zrh1(irpl) > rmdich) .AND. (zrhc(irpl) < rmdich)) THEN
!         zrhc  (irpl)  =  zrh1 (irpl)
!         zrhw  (irpl)  =  zrhw1(irpl)
!    (if 'zrhw1' has gross error, use 'zrh2' if well-defined)
!         IF ((zrhw1(irpl) > rtshlm) .AND. (zqx(irpl) > rmdich))               &
!           zrhc  (irpl)  =  rmdi
!       ENDIF
!     ENDIF
!     IF (lexiqx) THEN
!      PRINT *, 'zrhqx ',zrh2(irpl)
!       IF ((zrh2(irpl) > rmdich) .AND. (zrhc(irpl) < rmdich)) THEN
!         zrhc  (irpl)  =  zrh2 (irpl)
!         zrhw  (irpl)  =  zrhw2(irpl)
!       ENDIF
!     ENDIF
!     PRINT *, 'zrh3 ',irpl, zrhc(irpl)

! gross error if original (i.e. not model compatible) relative humidity > 102 %
! (previous version: if > 120 %)
      IF ((zrhw(irpl) > rtshlm) .AND. (zrhw(irpl) > rmdich)) THEN
        IF (ioberr(irpl) == 1) neventd (neqbig,icma) =                         &
                               neventd (neqbig,icma) + iactr
        nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
        ioberr (irpl) = 0
        zrhc  (irpl)  =  MAX( MIN( zrhc(irpl) , c1 ) , epsy )
      ENDIF

! bias correction of relative humidity (near saturation only)
! ------------------------------------
      IF ((zrhc(irpl) >= rhtsat(2)) .AND. (zrhc(irpl) < c1-epsy)) THEN
        IF ((ztt(irpl) <  tmelt) .AND. (ioberr(irpl) == 1))                    &
          neventd (neqsam,icma) = neventd(neqsam,icma) + iactr
        IF ((ztt(irpl) >= tmelt) .AND. (ioberr(irpl) == 1))                    &
          neventd (neqsap,icma) = neventd(neqsap,icma) + iactr
!       nflgx (irpl)  =  IBSET( nflgx(irpl), nvfbps(5) )
        nflgx (irpl)  =  IBSET( nflgx(irpl), nvfbps(6) )
      ELSEIF (zrhc(irpl) > c1+epsy) THEN
        IF ((ztt(irpl) <  tmelt) .AND. (ioberr(irpl) == 1))                    &
          neventd (neqclm,icma) = neventd(neqclm,icma) + iactr
        IF ((ztt(irpl) >= tmelt) .AND. (ioberr(irpl) == 1))                    &
          neventd (neqclp,icma) = neventd(neqclp,icma) + iactr
!       nflgx (irpl)  =  IBSET( nflgx(irpl), nvfbps(6) )
      ENDIF
      IF (zrhc(irpl) >= rhtsat(2)) THEN
        zrhc   (irpl)  =  c1
      ENDIF
!     PRINT *, 'zqx9b  ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)

! observation error for 2-m humidity
      zoberr (irpl) =  0.01_wp *rherr1  + 0.0004_wp *ABS( fisd(irpl) )
!     PRINT *, 'zrh4 ',irpl, zrhc(irpl)
!     PRINT *, 'zrh4  ',irpl, ystidl(irpl), zrhc(irpl)

! fill ODR (humidity)
! -------------------
      ogpbdy (ngpob,nbgrh ) = zrhc  (irpl)
      ogpbdy (ngpob,nbgqer) = zoberr(irpl)
      IF ((zrhc(irpl) > rmdich) .AND. (zrhw(irpl) > rmdich))                   &
      ogpbdy (ngpob,nbgdrh) = zrhc (irpl)  -  zrhw(irpl)
      mogpbd (ngpob,nbgerr) = insert( mogpbd(ngpob,nbgerr), ioberr(irpl), nvrq )
      mogpbd (ngpob,nbgflg) = insert( mogpbd(ngpob,nbgflg), nflgx(irpl), nvfqbp)
    ENDIF
  ENDDO

!    DEALLOCATE ( zrhw1   , STAT=istat )
!    DEALLOCATE ( zrh1    , STAT=istat )
!    DEALLOCATE ( zrhw2   , STAT=istat )
!    DEALLOCATE ( zrh2    , STAT=istat )
!    DEALLOCATE ( zqvw    , STAT=istat )
!    DEALLOCATE ( zqv     , STAT=istat )

! ---
! IWV
! ---

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    IF ((mogphd(ngpob,nhpass) /= -1)) THEN
      kobtyp = mogphd(ngpob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1
      nflgx  (irpl) =  0
      zoberr (irpl) = c0
      ioberr (irpl) =  1

! if ZPD < 0: bad reporting practice
      IF ((zpd(irpl) > rmdich) .AND. (zpd(irpl) < 0.0_wp)) THEN
        IF (ioberr(irpl) == 1)  neventd (negmis,icma) =                        &
                                neventd (negmis,icma) + iactr
        iwv   (irpl) = rmdi
        zpd   (irpl) = rmdi
        nflgx (irpl) = IBSET( nflgx(irpl), nvfbps(5) )
        ioberr(irpl) = 0
      ENDIF
!!! discard whole report if both ZPD and IWV are missing !!!
!     (so that application of observation operator to missing value is excluded)
      IF ((zpd(irpl) < rmdich) .AND. (iwv(irpl) < rmdich)) THEN
        mogphd(ngpob,nhpass) = -1
        IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
        IF (ioberr(irpl) == 1)  neventr (nenoda,icma) =                        &
                                neventr (nenoda,icma) + iactr
      ENDIF
    ENDIF
    IF ((mogphd(ngpob,nhpass) /= -1)) THEN
! if station height or its diff to model orography too large for ZPD
      fisdpd =  ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd(irpl) ))     &
                * fisd(irpl)
      IF ((ogphed(ngpob,nhalt) > altopsu(2)) .OR. (fisdpd > doromx(2))) THEN
        IF (ioberr(irpl) == 1)  neventd (nevalt,icma) =                        &
                                neventd (nevalt,icma) + iactr
        nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(4) )
        ioberr (irpl) = 0
      ENDIF
    ENDIF
    zbias = rmdi
    IF (iwv(irpl) > rmdich)  zbias = c0
    ! pwd, zpd are given in [m], iwv is given in [kg/m2]
    ogpbdy (ngpob,nbgiwv) = iwv   (irpl)
    IF (zwd(irpl) > rmdich)  ogpbdy (ngpob,nbgzwd) = zwd(irpl) *c1000
    IF (zpd(irpl) > rmdich)  ogpbdy (ngpob,nbgzpd) = zpd(irpl) *c1000
!   ogpbdy (ngpob,nbgtze) = zpder (irpl)
    ogpbdy (ngpob,nbgbia) = zbias
    mogpbd (ngpob,nbgflg) = insert( mogpbd(ngpob,nbgflg), nflgx(irpl), nvfgbp)
    mogpbd (ngpob,nbgerr) = insert( mogpbd(ngpob,nbgerr), ioberr(irpl), nvriwv )
    mogpbd (ngpob,nbgerr) = insert( mogpbd(ngpob,nbgerr), ioberr(irpl), nvrzpd )
  ENDDO

!-------------------------------------------------------------------------------
! Section 6: Get observation errors for ZPD, height
!            (T-2m, RH-2m obs errors are already assigned)
!-------------------------------------------------------------------------------

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
! height / pressure obs error
    IF ((mogphd(ngpob,nhpass) /= -1) .AND. (ogpbdy(ngpob,nbgp) > rmdich)) THEN
      zlop = LOG( ogpbdy(ngpob,nbgp) )

      CALL obs_find_level ( nerlev, rolnlv, zlop , ilverr, fiperr )
!     ===================

      ogpbdy(ngpob,nbgzer) = ogpbdy(ngpob,nbgzer) + fiperr * oezgps(ilverr)    &
                                              + (c1-fiperr)* oezgps(ilverr+1)
    ELSE
      ogpbdy(ngpob,nbgzer) = rmdi
      mogpbd(ngpob,nbgerr) = IBCLR ( mogpbd(ngpob,nbgerr), nvrz )
    ENDIF
  ENDDO

! zenith path delay error
! note: the obs error for IWV as required for the feedobs files is set
!       in routine 'obs_print_ncfeedobs' of module 'src_obs_pr_ncfeedobs.f90' !
  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl
    ! if no obs error is given in the dataset, assign default value
    IF (zpder(irpl) <= rmdich)  zpder (irpl) = oezpd
    ! increase ZPD (IWV) obs error by 0.5* default variance (oezpd)
    ! per 150 m of scaled difference between station height and model orography
    fisdpd =  ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd(irpl) ))       &
              * fisd(irpl)
    ogpbdy (ngpob,nbgtze) = (zpder(irpl) + c05* oezpd *fisdpd /doromx(2)) *c1000
  ENDDO

!-------------------------------------------------------------------------------
! Section 7: Check existence of data in report
!-------------------------------------------------------------------------------

  DO irpl = 1 , nrepgp
    ngpob = ntotgp + irpl

    IF (mogphd(ngpob,nhpass) /= -1) THEN
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mogphd(ngpob,nhpass) == 0)  iactr  = 1

      iactx = -1
      ipasx = -1
      IF (BTEST( mogpbd(ngpob,nbgerr), nvrz   )) iactx (2) = 1
      IF (BTEST( mogpbd(ngpob,nbgerr), nvrt   )) iactx (3) = 1
      IF (BTEST( mogpbd(ngpob,nbgerr), nvrq   )) iactx (4) = 1
      IF (BTEST( mogpbd(ngpob,nbgerr), nvriwv )) iactx (1) = 1
      IF (ogpbdy(ngpob,nbgz  ) > rmdich) ipasx (2) = 1
      IF (ogpbdy(ngpob,nbgt  ) > rmdich) ipasx (3) = 1
      IF (ogpbdy(ngpob,nbgrh ) > rmdich) ipasx (4) = 1
      IF (ogpbdy(ngpob,nbgiwv) > rmdich) ipasx (1) = 1
      iactx (5)             = MAX( iactx(1), iactx(2), iactx(3), iactx(4))
      ipasx (5)             = MAX( ipasx(1), iactx(2), ipasx(3), ipasx(4))
      nzaexi                = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
 !(nzaexi ==  1): active data present
 !(nzaexi == -1): only passive data present
 !(nzaexi ==  0): no data at all in total report
      IF (nzaexi <= 0)                                                         &
        neventr (nenoda,icma) = neventr(nenoda,icma) + iactr
      IF ((nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)) THEN
        mogphd (ngpob,nhpass) =  2
        mogphd (ngpob,nhflag) =  IBSET ( mogphd(ngpob,nhflag), FL_NO_OBS )
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1

! flag reports which are to be discarded completely
      ELSEIF ((nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))) THEN
        mogphd (ngpob,nhpass) = -1
        IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ENDIF
!     PRINT *,'noctrj_s6 ', icma, irpl, cma(icma)%cnt_ps, cma(icma)%cnt_rj     &
!            , iactr, nzaexi, mosghd(nsgob,nhpass)
!     PRINT *,'ztt31 ', yosghd(nsgob), nsgob, iactx, ipasx, iactr              &
!                     , mosghd(nsgob,nhpass)
! control message if no data in report
      IF (nzaexi == 0) THEN
        ilen = 2 + istrej
      ENDIF
    ENDIF
  ENDDO

! this is done after obs_cdf_redundancy
! TODO : delete reports with (moxxhd(nxxob,nhpass) == -1)
! CALL obs_del_old_reports
! ========================

! clean up: de-allocate standardised arrays
! -----------------------------------------

  IF (nrepgp >= 1) THEN
    DEALLOCATE ( zrhw   , STAT=istat )
    DEALLOCATE ( zrhc   , STAT=istat )
    DEALLOCATE ( zpp    , STAT=istat )
    DEALLOCATE ( ztt    , STAT=istat )
    DEALLOCATE ( iwv    , STAT=istat )
    DEALLOCATE ( zwd    , STAT=istat )
    DEALLOCATE ( zpder    , STAT=istat )
    DEALLOCATE ( zpd    , STAT=istat )
    DEALLOCATE ( nflgg  , STAT=istat )
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_store_gps
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_store_gps


!===============================================================================

ELEMENTAL INTEGER FUNCTION insert  ( invar, inval, ibp )
  !-----------------------------------------------------
  INTEGER         , INTENT (IN)  ::  invar, inval, ibp
  !----------------------------------------------------------------------------
  ! inserts bits set in 'inval' into integer word 'invar' at bit position 'ibp'
  !----------------------------------------------------------------------------
  !
  insert = IOR( invar , ISHFT( inval, ibp ) )
  !
END FUNCTION insert

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_gps
