!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_cdfin_comhead

!-------------------------------------------------------------------------------
! Description:
!   This module reads from a NetCDF observation input file and stores in the ODR
!   those values of those variables in the report headers, which are commonly
!   available for all observation types
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_read_comhead   : called by obs_cdf_read_temp_pilot etc.,
!                               reads common header information into temporary
!                               arrays and exploits it, including horizontal
!                               grid point assignment to obs reports
!    - obs_cdf_buffer_comhead : called by obs_cdf_read_temp_pilot etc.
!                               copies info into special buffer arrays used
!                               for distribution to sub-domains (nodes)
!    - obs_cdf_store_comhead  : called by obs_cdf_store_multilev etc.,
!                               reads info from buffer arrays and stores it in
!                               the ODR
!
!   This module also contains elemental functions, formerly statement functions:
!   - ibit1        : returns 1 bit at given bit position of given integer word
!
!   It uses from:
!    - src_obs_cdfin_util:    - obs_assign_gridpt
!    - utilities:             - diff_minutes
!                             - get_utc_date
!    - environment:           - model_abort
!
!   Data modules used:
!    - kind_parameters
!    - data_obs_lib_cosmo
!    - data_obs_cdfin
!    - data_obs_record
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
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted (e.g.
!  modified routine interfaces and diagnostic arrays, 'obs_pointrs --> i_cma').
!  Call of external routine 'difmin' replaced by new routine 'diff_minutes'
!  from utilities (type of routine arguments: iintegers instead of intgribf).
!  Karolin Eichler: each original GPS processing centre gets its own new CMA
!  observation code type.
! V4_23        2012/05/10 Ulrich Schaettler
!  Adapted interface to SR diff_minutes (added itype_calendar=0 hard coded)
! V4_27        2013/03/19 Christoph Schraff
!  (Non-zero) minutes of model reference time 'ydate_ref' considered.
! V4_28        2013/07/12 Christoph Schraff
!  Relaxed conditions on variable for station ID in NetCDF obs input files.
!  Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Extensions for reading / processing Mode-S aircraft observations, from NetCDF
!  files with 2 different templates: (i) converted from BUFR as received from
!  Mode-S processing centre KNMI, (ii) converted into DWD-ACARS template. (CS)
!  Also check for unrepresentable characters in the station ID and convert them
!  into blank.
!  Replaced ireals by wp (working precision) (OF)
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' reduced from 3 to 2.
!  Variables related to the AOF interface removed.
! V5_4f        2017-09-01 Christoph Schraff
!  Minor fix in redundancy check to prefer TEMP parts A, B over parts C, D.
!  Update from DACE V1_51: Adapted call of 'obs_assign_gridpt' (computation of
!  obs location in grid point coordinates now moved into 'obs_assign_gridpt').
! V5_4g        2017-11-13 Christoph Schraff
!  Fix to read reports with section1_data_category=6 and 
!    section1_int_data_sub_category=10 (which is not officially approved value).
! V5_4h        2017-12-15 Christoph Schraff
!  Adaptations (data category, station ID) to read AMV reports from cdfin_satob.
! V5_6         2019-02-27 Christoph Schraff
!  WMO marine obs platform id 'NMOPI' introduced for new buoy template.
!  Comments on buoy station type (MTODB instead of MQOBL) corrected.
! V5_6b        2019-10-16 Christoph Schraff
!  - Inclusion of processing of tower profile reports (with special treatment
!    for grid point assignment)  and descending TEMP reports.
!  - Extensions for descending TEMP reports with temporal and horizontal drift.
!  - Adaptations and new code types for high-resolution radiosonde reports.
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! V5_7a        2020-05-11 Christoph Schraff
!  - Processing of ground-based wind lidar reports added.
!  - Processing of test Synop surface reports added.
!  - Corrected search for the lowest level in a descending sonde report now
!    accounts for the pressure (and does not depend any more on the order
!    (ascending or descending) of the levels in the report).
!  - Int. data sub category added for descending sonde BUFR 3 reports.
! @VERSION@  @DATE@     Christoph Schraff
!  - Correction of lethal bug and more generalized formulation for tagging a
!    particular level as 'surface level' in the descending radiosonde report.
!    Now, time and horizontal shift values are allowed to be missing in any
!    combination. 
!  - Correction of a (minor) bug at horizonal shifting of descending TEMP.
!  - Adaptation to avoid a lethal crash and allow for the mandatory multi-level
!    variables being missing if all the reports in the file have zero observ.
!    levels; in this case, a Caution message is issued to file YUCAUTN.
!  - Some code clean-up.
!  - 16.05.2022: 'nhssos' added.
!  - 31.05.2023: Shift obs time of wind profiler by averaging time if required.
!                Assign code type ICOS tower if surface pressure reported.
!  - 07.10.2024: Bug fix: always use ABS in the context of rmisschk, assuming
!                that it's not known whether rmiss is huge positive or negative.
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

USE netcdf            ! NetCDF F90 interface
USE mo_netcdf_param   ! parameters in file 'netcdf.inc' of the NetCDF package

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
!   c2         ,& ! standard real constant 2.0
!   c05        ,& ! standard real constant 0.5
!   c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    ! horizontal and vertical sizes of the model fields
    ie_tot         ,& ! number of grid pts in zonal direction (in total domain)
    je_tot         ,& ! number of grid pts in meridional dir. (in total domain)

    ! variables related to parallelisation / domain decomposition
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator

    ! report dimensions of ODR (Observation Data Record, for observation storage
    ! on local sub-domains), and other variables related to namelist parameters
    maxmll         ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl         ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl         ,& ! size (report dimension) of the (ground-based) GPS  ODR
                      !
    nolbc          ,& ! number of grid rows at lateral boundaries
                      !   where obs are neglected
    itim_wp           ! mode of correction of obs time for wind profiler

USE data_obs_lib_cosmo, ONLY :  &

    ! variables related to namelist parameters, and other variables
    doromx         ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    acthr          ,& ! actual model time [hours] with respect to 'ydate_ref'

    ! switches related to namelist parameters and other
    lverpas        ,& ! write also passive reports on feedback files
    lwonl          ,& ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open
!   lgpsbias          ! .t. ==> bias correction to GPS IWV applied
    ydate_ref         ! reference date (e.g. start of the forecast)
                      ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

USE data_obs_lib_cosmo, ONLY :  &

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yucautn    ,& ! caution messages if too many obs for ODR size
    yurejct    ,& ! direct reporting of rejected obs. reports
    nucautn    ,& ! caution messages if too many obs for current ODR size
    nurej      ,& ! direct reporting of rejected obs. reports
    nupr       ,& ! all the remaining information
    lopen_rej     ! .true. if yurejct is open

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports)
    nsatob     ,& ! SATOB reports
    ndribu     ,& ! DRIBU reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    nsatem     ,& ! SATEM reports
    nsattv     ,& ! SATEM reports
    ngps       ,& ! GPS   reports
    nscatt     ,& ! SCATT reports (from NetCDF only)
    nsrscd     ,& !   synop surface report
    natscd     ,& !   automatic synop surface report
    nshscd     ,& !   ship synop report
!   nabscd     ,& !   ship synop abbreviated report
!   nshred     ,& !   shred report
    natshs     ,& !   automatic ship synop report
    nmetar     ,& !   Metar
    nsytst     ,& !   test synop
    naircd     ,& !   aircraft report
    ncodar     ,& !   codar report
!   ncolba     ,& !   colba report
    namdar     ,& !   amdar report
    nacar      ,& !   acar  report
    nmodes     ,& !   mode-s report
!   nstbcd     ,& !   satob report
!   nhrvis     ,& !   high-res. VIS wind report
    namv       ,& !   AMV   report
!   nsst       ,& !   sst report
    ndrbcd        !   dribu report
!   nbathy        !   bathy report
!   ntesac        !   tesac report

USE data_obs_lib_cosmo, ONLY :  &
    nldtcd     ,& !   temp land   report
    nshtcd     ,& !   temp ship   report
    nmotcd     ,& !   temp mobile report
    ntdrop     ,& !   temp drop   report
    nbtemp     ,& !   temp land   report (high-res. BUFR)
    nbship     ,& !   temp ship   report (high-res. BUFR)
    nbdrop     ,& !   temp drop   report (high-res. BUFR)
    ntdesc     ,& !   temp high-res. descending report (HA 02.04.19: as EC)
!   nrocob     ,& !   rocob      report
!   nrocsh     ,& !   rocob ship report
    nldpcd     ,& !   pilot land   report
    nshpcd     ,& !   pilot ship   report
    nmopcd     ,& !   pilot mobile report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    ntower     ,& !   tower profile  report
    ntowic     ,& !   icos tower profile report
!   npr_us     ,& !   US Wind Profiler/RASS report
    nwlidr     ,& !   ground-based wind lidar report
    nstmcd     ,& !   satem report
    nstovs     ,& !   high resolution ATOVS satellite data
    nascat     ,& !   ASCAT scatterometer report
    nqscat        !   QuickScat scatterometer report
!   ngpgfz     ,& !   GPS report processed by GFZ
!   ngpasi     ,& !   GPS report processed by ASI
!   ngpbkg     ,& !   GPS report processed by BKG
!   ngpgop     ,& !   GPS report processed by GOP
!   ngpieec    ,& !   GPS report processed by IEEC
!   ngpsgn     ,& !   GPS report processed by SGN
!   ngplpt     ,& !   GPS report processed by LPT
!   ngpmet     ,& !   GPS report processed by MET
!   ngprob     ,& !   GPS report processed by ROB
!   ngpige     ,& !   GPS report processed by IGE
!   ngpknmi    ,& !   GPS report processed by KNMI
!   ngpnga        !   GPS report processed by NGA

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &
!   n_cma      ,& ! number of CMA obs and code types
    cma        ,& ! array of meta data on CMA observation and code types
    n_gpc      ,& ! number of GPS processing centres
    gpc        ,& ! array of meta data on GPS processing centres
!   n_gpc_offs ,& ! index offset of GPS code type elements in array 'cma'

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
    ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_temphirs  ,& ! indicator for proc. NetCDF TEMP high-res BUFR input
    ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_tempdesc  ,& ! indicator for proc. NetCDF descending TEMP    input
    ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
    ncdf_amdar     ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_acars     ,& ! indicator for proc. NetCDF ACARS single-level input
    ncdf_modes     ,& ! indicator for proc. NetCDF MODE-S KNMI format input
    ncdf_modes_acr ,& ! indicator for proc. NetCDF MODE-S ACARS fmt.  input
    ncdf_wprof     ,& ! indicator for proc. NetCDF wind profiler      input
    ncdf_rass      ,& ! indicator for proc. NetCDF RASS profiler      input
    ncdf_radar_vad ,& ! indicator for proc. NetCDF radar wind prof.   input
    ncdf_synop     ,& ! indicator for proc. NetCDF SYNOP              input
    ncdf_synop_mob ,& ! indicator for proc. NetCDF SYNOP mobile       input
    ncdf_ship      ,& ! indicator for proc. NetCDF SHIP               input
    ncdf_buoy      ,& ! indicator for proc. NetCDF BUOY               input
    ncdf_metar     ,& ! indicator for proc. NetCDF METAR sfc aviation input
    ncdf_gps_zenith,& ! indicator for proc. NetCDF GPS (ZPD / IWV)    input
    ncdf_ascat     ,& ! indicator for proc. NetCDF ASCAT scatterom.   input
    ncdf_qscat     ,& ! indicator for proc. NetCDF QuickScat scatter. input
    ncdf_satob     ,& ! indicator for proc. NetCDF SATOB wind         input
    ncdf_acars_uk  ,& ! indicator for proc. NetCDF ACARS UK + Canada  input
    ncdf_acars_us  ,& ! indicator for proc. NetCDF ACARS US w. humid. input
    ncdf_tower     ,& ! indicator for proc. NetCDF tower profile      input
    ncdf_wlidar_wp ,& ! indicator for proc. NetCDF ground-b wind lidar input
    ncdf_synop_tst ,& ! indicator for proc. NetCDF test SYNOP         input
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
!   nc_dim_len     ,& ! length of a dimension          (only temporary variable)
! BUFR section 1 and section 2 entries  (common to all obs types, stored in ODR
!                                             except if only temporary variable)
    kcat           ,& ! data category     (WMO Common Code Table C13)
    kcatsub        ,& ! data sub-category (WMO Common Code Table C13)
    kcentre        ,& ! data centre       (WMO Common Code Table C11)
    kcensub        ,& ! data sub-centre   (WMO C12)    (only temporary variable)
    kupdate        ,& ! update sequence number (indicates station correction)
    kz_dwd         ,& ! DWD-internal classification number (Kennzahl KZ)
    nsynhr         ,& ! nominal (synoptic) hour [yymmddhh] (from s1_date,_time)
! common header entries in NetCDF file
    nc_mjjj        ,& ! MJJJ   : year                  (only temporary variable)
    nc_mmm         ,& ! MMM    : month                 (only temporary variable)
    nc_myy         ,& ! MYY    : day                   (only temporary variable)
    nc_mgg         ,& ! MGG    : hour                  (only temporary variable)
    nc_ngg         ,& ! NGG    : minute                (only temporary variable)
!   nc_msec        ,& ! MSEC   : second                (only temporary variable)
!   nc_tisi        ,& ! MTISI  : time significance (BUFR Table 008021)    (dito)
    nc_mii         ,& ! MII    : WMO block number      (only temporary variable)
    nc_niii        ,& ! NIII   : WMO station number    (only temporary variable)
    nc_alt         ,& ! MHP, MH: station altitude      (only temporary variable)
    nc_altq        ,& ! MSEQM  : altitude quality (mobil, BUFR 033024)    (dito)
    nc_locq        ,& ! MQOBL  : location quality (buoy , BUFR 033023)    (dito)
! derived common header variables stored to ODR
    iobtyp         ,& ! observation type
    icdtyp         ,& ! observation code type
    istidn         ,& ! station number
    nr_date        ,& ! absolute date [YYYYMMDD]
    nr_time        ,& ! absolute time [HHMM]
    isurfob        ,& ! indicat. if orographic surface report conditions are met
                      !  (for descend. TEMP: minus index of lowest usable level)
    mrepflg        ,& ! report flags
    kflag          ,& ! processing flag
    iobstot        ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot        ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc        ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc        ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc            ! indices of reports to be processed now

USE data_obs_cdfin, ONLY :  &

! common header entries in NetCDF file     (stored to ODR or temporary variable)
    rc_lat         ,& ! MLAH, MLALA : latitude
    rc_lon         ,& ! MLOH, MLOLO : longitude
    rc_altp        ,& ! MHOBNN      : barometer altitude    (temporary variable)
    rc_alt         ,& ! MHOSNN      : station   altitude    (temporary variable)
! derived common header variables stored to ODR
    zr_hour        ,& ! observation time relative to model initial time [hour]
    zsynhr         ,& ! nominal report time rel.  to model initial time [hour]
    zaltob         ,& ! obs station altitude
    zaltmo         ,& ! model orography (at the grid point assigned to the obs)
    rio_tot        ,& ! longitude in model grid point units
    rjo_tot        ,& ! latitude  in model grid point units
    ztdecdb        ,& ! data base decoding time rel to model initial time [hour]
!   yc_dim_name    ,& ! name of dimension in NetCDF file
    ync_stid       ,& ! YWIGLI, YDDDD, YSSOSN, YAIRN, YCCC8, YXXNN, YSOSN,YCCCC:
                      !   station id
    ystidn         ,& ! station identity

!         1.2.2 other NetCDF header entries
!               ---------------------------

    nc_nix         ,& ! NIX   , BUFR Table B 002001 : station type (man,auto,..)

!         1.5   auxilliary buffer arrays and variables
!               --------------------------------------

    imiss          ,& ! missing value for integers in current NetCDF input file
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk          ! value smaller than 'rmiss', to check for missing value

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
!-------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
!   mxreve     ,& ! length of report event counter array
    nedbfl     ,& ! data base flag on loc/tim/alt high
!   netime     ,& ! time out of range (too old)
    nenoal     ,& ! no station altitude
    neloca     ,& ! station location out of range
    nezdif     ,& ! distance 'model orography - station altitude' too large
    neblak     ,& ! blacklisted ship
    neobct     ,& ! observation or code type excluded on area with sta. location
    nesodr     ,& ! report number exceeding size of ODR (==> adjust Namelist)
    nenoda     ,& ! no accepted data in report
!   nenops     ,& ! pressure too small (< 20hPa), or missing in aircraft report
!   neredn     ,& ! redundancy between 2 multi-level, or 2 single-level reports
!   neredx     ,& ! redundancy between 1 multi- and 1 single-level report
!   neslml     ,& ! one multi-level report made from  single-level reports
!   neslps     ,& ! single-level report put in multi-level report and set passiv
!   nenoml     ,& ! multi-levl report not made due to ODR array size
!   netrac     ,& ! (flight) track suspicious
!   nethin     ,& ! thinning of aircraft (flight) track

!         3.2    Event counter arrays
!                --------------------
    neventr    ,& ! counter array of report events

! Section 8 : Temporary global model fields
!-------------------------------------------------------------------------------

    hsurf_tot  ,& ! total array of model surface height
    fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
!   mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
    nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
    nhzjo      ,& ! latitude  of obs. station in grid pt. units
    nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
    nhtddb     ,& ! data base decoding time in forecast hours
    nhsolz     ,& ! solar zenith angle [deg]
    nhssos     ,& ! standard devation of sub-grid scale model orogrophy [m]
    nhvcbu     ,& ! correction factor to vertical correlation scale for wind
                  ! at base of report
    nhvcbt     ,& ! as 'nhvcbu', but for temperature
    nhvcbq     ,& ! as 'nhvcbu', but for humidity
    nhvctu     ,& ! correction factor to vertical correlation scale for wind
                  ! at top of report
    nhvctt     ,& ! as 'nhvctu', but for temperature
    nhvctq        ! as 'nhvctu', but for humidity

USE data_obs_record, ONLY :   &

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
!   mxghdf     ,& ! header length of GPS reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! threshold quality control flags for pressure, status of
                  ! verification output
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
    nhcorr     ,& ! update sequence number (station correction indicator)
    nhcat      ,& ! data     category (from BUFR Section 1)
    nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre
    nhstid     ,& ! station identity number / satellite ID (WMO Table C-5)
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
    nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
!   nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
    nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
    nhuexi     ,& ! flag for existence of wind data        in multi-level report
    nhtexi     ,& ! flag for existence of temperature data in multi-level report
    nhqexi     ,& ! flag for existence of humidity data    in multi-level report
!   nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
!   nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
!   nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
!   nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
!   nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
!   nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]
!   nhstyp     ,& ! surface obs: station type (buoy: MTODB, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)


!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid        ! character length of the station identity
!   ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
!   nvpabp     ,& ! bit pos. for report set passive since it is     nvschr
                  !              used in a multi-level pseudo report
    nvpsbp     ,& ! bit pos. for report set passive since at least    "
                  !              1 of the following 5 flags or
                  !              the flight track flag applies
    nvobbp     ,& ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
    nvalbp     ,& ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
    nvbkbp     ,& ! bit pos. for flag: 'blacklisted station (ship)'   "
    nvexbp     ,& ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
!   nvrdbp     ,& ! bit pos. for flag: 'redundant report'             "
    nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvscbp     ,& ! bit pos. for station correction indicator         "
    nvssbp        ! bit pos. for station suspicion indicator          "
!   nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
!   nvapoc     ,& ! no. of bits occ. by phase of flight               "
!   nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
!   nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    fdoro      ,& ! scaling factors to obs height minus model orography

!       2.     Observation data records (ODR) and associated arrays
!       ----------------------------------------------------------------------

!       2.1    Formats of observation data records
!              -----------------------------------
    omlhed     ,& ! header of multi-level ODR
    osghed     ,& ! header of single-level ODR
    ogphed     ,& ! header of GPS ODR
    momlhd     ,& ! header of multi-level ODR
    mosghd     ,& ! header of single-level ODR
    mogphd     ,& ! header of GPS ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd        ! header of GPS ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_OBSTYPE    ,& ! passive report type (at obs. location)
    FL_BLACKLIST  ,& ! blacklist (or not whitelist)  (in particular:
                     !   bad (hard blacklisted) aircraft station identity)
    FL_SUSP_LOCT  ,& ! suspicious location (3D)
    FL_AREA       ,& ! location not in valid area (i.e. outside user-def. area)
    FL_HEIGHT     ,& ! location not in valid height range
    FL_SURF          ! incorrect surface (land,ice,etc), e.g. sea obs over land

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    diff_minutes    ,& ! compute difference in minutes between 2 dates/times
    get_utc_date       ! get absolute date / time given a relative time

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_gridpt  ! horiz. assignment of an obs report to a grid pt.

!-------------------------------------------------------------------------------

! USE mo_satid,                 ONLY :  &
!   satname            ! derive satellite name (8char mnemonic) from satid

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE
PUBLIC :: obs_cdf_read_comhead,  obs_cdf_buffer_comhead,  obs_cdf_store_comhead

! LOGICAL  ,  PARAMETER    ::  &
!   lax  =  .TRUE.    ! (re-)lax(ed) conditions on required variables in
!                     !   NetCDF feedobs files

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_comhead" for reading/evaluat. header info
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_comhead ( min_sta , min_end , ilcdf                    &
                                , mrepsta , nrep , nreproc )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_comhead" reads the values
!   of those variables in the report header, which are commonly available
!   for all observation types, from a NetCDF observation input file.
!   This includes BUFR section 1 and 2 variables (e.g. data category),
!   station identity (report identity), observation time and location,
!   as well as station altitude.
!   Only for those reports for which the observation time is in the wanted
!   time range, the read pieces of information are stored in temporary arrays,
!   and some information is already evaluated, e.g. in order to assign each
!   report to a model grid point.
!
! Method:
!   First, the record dimension is determined, and the observation time is read
!   for all reports. The interval in the record dimension is determined which
!   contains all reports within the wanted time range, and for subsequent
!   variables, data are read only from this interval.
!   ODR observation type and code type is derived from BUFR / NetCDF data (sub)
!   category. The station identity is obtained by checking a sequence of
!   possible station identity NetCDF variables for availability and reading
!   the first one which is actually available. Barometer altitude is preferred
!   over standard station altitude, etc.
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

  INTEGER        , INTENT (OUT) ::  &
    mrepsta       ,& ! index (in record dimension) of first report to be read
    nrep          ,& ! index interval (in record dimension) which contains
                     !   those reports that should be read now
    nreproc          ! number of reports to be read and processed now

! Local parameters:
! ----------------

  INTEGER        , PARAMETER  :: &
    itype_calendar = 0    ! for specifying the calendar used in 'get_utc_date'

! Local scalars:
! -------------

  LOGICAL                  ::  &
    lvarid           ,& ! variable for station ID exists in NetCDF file
    lsfcob           ,& ! orographic surface report conditions are met
    loricen          ,& ! entry 'MMIOGC' for original data center exists
    lupair           ,& ! upper-air report
    lseaobs          ,& ! sea observation type
    lseaob2          ,& ! sea observation type
!   ldex1   ,ldex    ,& ! horizonal shift exists for this level / for any level
    lvertfit         ,& ! .true, if then minimum vertical distance 'dz' between
                        !           station height and model orography has given
                        !           absolute priority even if an 'appropriate
                        !           (|dz| < 40m, zdo < 0.25) grid point nearby'
                        !           is available
    lindenb          ,& ! true if station Lindenberg: special grid pt assignment
    lwrej               ! true if messages are written to file unit nurej

  INTEGER        ::  &
    ianyy   ,imoyy   ,& ! year     / of the date and time
    ianmm   ,imomm   ,& ! month   /  of the current
    iandd   ,imodd   ,& ! day    <   observation
    ianhh   ,imohh   ,& ! hour    \  resp. of the initial
    ianmin  ,imomin  ,& ! minute   \ model date and time
    imindif          ,& ! time difference [min] between these two times
    ierrf               ! error status for time difference routine

  INTEGER        ::  &
    kcdftyp          ,& ! type of NetCDF input file format (observation type):
                        !  = ncdf_temp      : TEMP
                        !  = ncdf_tempship  : TEMPSHIP
                        !  = ncdf_temphirs  : high-resolution BUFR TEMP
                        !  = ncdf_tempdesc  : descending TEMP
                        !  = ncdf_pilot     : PILOT (z-levels)
                        !  = ncdf_pilot_p   : PILOT (p-levels)
                        !  = ncdf_wprof     : WIND PROFILER
                        !  = ncdf_synop     : SYNOP
                        !  etc.
    ncid             ,& ! file unit number of current NetCDF obs input file
    dimid_reps       ,& ! dim. ID of report (record) dimension in NetCDF file
    dimid_mxlv       ,& ! dimension ID for vertical levels in NetCDF file
    ndims            ,& ! number of dimensions (for 'edition_number' in NetCDF)
    natts            ,& ! number of attributes (for 'edition_number' in NetCDF)
    status           ,& ! NetCDF status variable
    ntotrep          ,& ! length of record dimension in NetCDF file
    ireptot          ,& ! NetCDF file record index over reports
    mrepend          ,& ! index (in record dimension) of last report to be read
    mrepstao         ,& ! index (in record dimension) of first report to be read
    irep   , irepn   ,& ! record index over reports
    irpa             ,& ! loop   index over reports to be processed further
    ily              ,& ! loop   index over single characters
    icma             ,& ! pointer index for 'cma' / diagnostic arrays
    kobtyp , kcdtyp  ,& ! observation type , observation code type
    nerrobt          ,& ! number of reports with unknown data (sub-)category
    ireperr          ,& ! index of last report with unknown data (sub-)category
    dimid_ystid      ,& ! dimension ID for character station ID in NetCDF file
    ilen             ,& ! character length of station ID entries in NetCDF file
    nsynmin          ,& ! synoptic time [hhmm]
!   mbaro            ,& ! = 1 if barometer altitude is available
    nelevt           ,& ! station altitude
    iob_tot, job_tot ,& ! grid point assigned to obs (total model area)
    nzday            ,& ! Julian day
    ic               ,& ! loop index (over GPS processing centres)
    maxlev           ,& ! max number of vertical levels in all processd reports
    iext             ,& ! for descend. TEMP: indicator on time and position info
    ndsec            ,& ! for descend. TEMP: time shift of surface level [sec]
    nsta2 (2)        ,& ! start indices       \  of the reports to be read
    ncnt2 (2)        ,& ! number of elements  /  in the external NetCDF
    ilev             ,& ! loop index over vertical levels
    istat               ! error indicator

  INTEGER        ::  &  !   variable ID's in NetCDF file for:
    varid_edno                 ,& ! (BUFR) edition number
                                  ! - BUFR section 1 and 2 variables:
    varid_s1date, varid_s1time ,& ! nominal (synoptic) date , time
    varid_s1cat , varid_s1cats ,& ! data category , data sub category
    varid_s1cent, varid_s1cents,& ! data centre   , data sub centre
    varid_s2ikz , varid_s1updat,& ! DWD-internal classifier, update sequence no.
                                  ! - BUFR section 3 for data variables:
    varid_s2dbdd, varid_s2dbdt ,& ! data base decoding date , time
    varid_mjjj  , varid_mmm    ,& ! year  , month   \  of exact
    varid_myy   , varid_mgg    ,& ! day   , hour     > observation
    varid_ngg                  ,& ! minute, second  /  time
    varid_dts                  ,& ! obs period [s]  (for wind profiler)
    varid_mii   , varid_niii   ,& ! WMO block number , WMO station number
    varid_mabnn , varid_nmopi  ,& ! Buoy/platform id., marine obs platform id.
    varid_mi1i2                ,& ! WMO satellite identifier
    varid_ystid                ,& ! any type of character station ID
    varid_lat   , varid_lon    ,& ! latitude , longitude
    varid_altp  , varid_alt    ,& ! barometer altitude , (real) station altitude
    varid_altq  , varid_alti   ,& ! altitude quality mark , int station altitude
    varid_locq                 ,& ! location quality mark
    varid_oricen               ,& ! data centre of original ACARS report
    varid_dlat  , varid_dlon   ,& ! latitude shift, longitude shift
    varid_dt    , varid_p      ,& ! time shift    , pressure
    varid_t     , varid_ff     ,& ! temperature   , wind speed
    varid_nix                  ,& ! stati type (0:auto,1:manned,2:hybrid,3:miss)
    varid_mand                    ! a mandatory variable in multi-level reports

  REAL (KIND=wp) ::  &
    zsurf            ,& ! height of model grid pt. to which obs. is assigned
    zfdoro(5)        ,& ! scaling factor to obs height minus model orography
    zio_tot, zjo_tot ,& ! observation location in grid point units (total area)
    roblat , roblon  ,& ! geographical latitude and longitude
    rlat   , rlon    ,& ! latitude and longitude in the rotated system
    zhr              ,& ! hour
    pmax                ! maximum pressure in descent report
! REAL (KIND=wp) ::  &
!   rmod             ,& ! name of statement function
!   ztti   , ziv        ! dummy arguments of statement function

  CHARACTER (LEN=nf90_max_name)    :: &
    ydim_reps           ! name of report (record) dimension in NetCDF file
  CHARACTER (LEN=12)       :: &
    yobdat              ! date and time of current observation
  CHARACTER (LEN= 7)       :: &
    ymnem               ! mnemonics for mandatory multi-level variable
  CHARACTER (LEN= 6)       :: &
    ymnem_id            ! mnemonics for character station ID
  CHARACTER (LEN=13)       :: &
    ymnem_strl          ! mnemonics for character station ID string length
  CHARACTER (LEN= 1)       :: &
    yfillchar           ! character fill value of NetCDF file
  CHARACTER (LEN=25)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=90)       :: &
    yerrmsl             ! error message
  CHARACTER (LEN=30)       :: &
    yerr                ! error message
  CHARACTER (LEN= icdfinlen + iannexlen)  ::  &
    yfn                 ! file name of NetCDF observation input file, with annex
  CHARACTER (LEN=14)       :: &
    yakdat              ! date (format: yyyymmddhhmmss: year,month,day,hr,min,s)
  CHARACTER (LEN=28)       :: &
    yakdat2             ! date in the form: wd   dd.mm.yy  hh mm ss UTC

! Local arrays:
! ------------

  LOGICAL        , ALLOCATABLE :: &
    lwrsfc       (:)    ! writing message into file 'nurej'

  INTEGER        , ALLOCATABLE :: &
    kproc        (:) ,& ! flags which indicate no further processing
    nt_date      (:) ,& ! absolute date [YYYYMMDD]
    nt_time      (:) ,& ! absolute time [HHMM]
    nsyndat      (:) ,& ! absolute nominal (synoptic) obs date [YYYYMMDD]
    nsyntim      (:) ,& ! absolute nominal (synoptic) obs time [HHMM]
    ndbddat      (:) ,& ! data base report decoding date [YYYYMMDD]
    ndbdtim      (:) ,& ! data base report decoding time [HHMM]
    min_ob       (:) ,& ! obs time relative to model initial time [min]
    min_dtl      (:) ,& ! obs time minus time when buoy location is determined
    nedno        (:) ,& ! (BUFR) edition number
    mpassiv      (:) ,& ! passive report flag word
    icmaa        (:) ,& ! pointer index for 'cma' / diagnostic arrays
    koricen      (:) ,& ! data sub-centre   (WMO C12)  (only temporary variable)
    isurflv      (:) ,& ! for descend. TEMP: minus index of lowest usable level
    iexl         (:) ,& ! for descend. TEMP: indicator on time and position info
    nc_dts       (:) ,& ! obs period [s] (for wind profiler)
    nc_dt      (:,:)    ! time displacement of obs relative to report

  REAL (KIND=wp) , ALLOCATABLE :: &
    rc_p       (:,:) ,& ! obs pressure
    rc_t       (:,:) ,& ! obs temperature
    rc_ff      (:,:) ,& ! obs wind speed
    rc_dlat    (:,:) ,& ! latitude  displacement of obs relative to report
    rc_dlon    (:,:) ,& ! longitude displacement of obs relative to report
    rc_dlatr     (:) ,& ! latitude  displacement of obs relative to report
    rc_dlonr     (:)    ! longitude displacement of obs relative to report
!
!============ End of header ====================================================

!CSC:   statement function, to remove (see below)
!   rmod (ztti,ziv) = ztti - ziv *INT( ztti/ziv + epsy ) + epsy
!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_comhead
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_comhead'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

!-------------------------------------------------------------------------------
! Section 1: Get - number of dimensions, variables, and attributes
!                - dimension names and lengths used in the NetCDF file
!                - number of reports and corresponding dimension ID in the file
!                - value of '_FillValue' ( = missing value) for integers
!-------------------------------------------------------------------------------

! inquire about number of dimensions, variables, and attributes, and
!         about umlimited dimension
!-------------------------------------------------------------------

! status = nf90_Inquire (ncid, ncdims, ncvars, ncatts, nc_unlim_dimid)
! IF (status /= nf90_noerr) CALL handle_err (status)

! get dimension names and lengths used in the NetCDF file
!--------------------------------------------------------

! ALLOCATE (yc_dim_name (ncdims) , STAT=istat )
! ALLOCATE (nc_dim_len  (ncdims) , STAT=istat )

! DO idim = 1 , ncdims
!   status = nf90_Inquire_Dimension (ncid, idim, name=yname, len=ilen)
!   yc_dim_name (idim)  =  yname
!   nc_dim_len  (idim)  =  ilen
! ENDDO

! assume that the variable 'edition_number' always exists and
! has 1 dimension equal to the number of reports in the file;
! the dimension ID of this dimension can thus be determined
! -----------------------------------------------------------

  ndims  = 0
  status = nf90_inq_varid (ncid, 'edition_number', varid_edno)
  IF  (status == nf90_noerr)                                                   &
    status = nf90_Inquire_Variable ( ncid, varid_edno, ndims=ndims             &
                                   , dimids=dimids   , natts=natts )
  IF ((status /= nf90_noerr) .OR. (ndims /= 1)) THEN
    yerrmsl = 'VARIBLE edition_number DOES NOT EXIST OR HAVE DIMENSION 1 IN '  &
             // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF
  dimid_reps = dimids(1)
  status = nf90_Inquire_Dimension ( ncid, dimid_reps, name=ydim_reps           &
                                                    , len =ntotrep )

! get value of '_FillValue' for integers
! --------------------------------------
! ('_Fillvalue' denotes missing value in the NetCDF file; assume here that
!  the same value for _FillValue is used for all integers in the file,
!  and that 'edition_number' has an scalar attribute '_FillValue')

  status = nf90_get_att (ncid, varid_edno, '_FillValue', imiss)

! get value of '_FillValue' for floats
! ------------------------------------
! ('_Fillvalue' denotes missing value in the NetCDF file; assume here that
!  the same value for _FillValue is used for all floats in the file)

                            status = nf90_inq_varid (ncid,'MLAH'  ,varid_lat  )
  IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MLALA' ,varid_lat  )
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'STANDARD LATITUDE VARIABLE DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF
  status = nf90_get_att (ncid, varid_lat, '_FillValue', rmiss)

  rmisschk = ABS( rmiss ) * 0.9_wp

! if number of reports in the NetCDF file is zero or undefined
! then set (nrep=0) to indicate to leave this and calling routines
! ----------------------------------------------------------------
  IF ((ntotrep <= 0) .OR. (ntotrep == imiss)) THEN
    nrep  =  0
    PRINT '(" CAUTION: number of reports in file ",A," is zero or undefined")' &
           , yfn (1:LEN_TRIM(yfn))
    RETURN
!   ~~~~~~
  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Get observation time and
!            the interval of the record dimension which has to be read now
!-------------------------------------------------------------------------------

!   (these variables are de-allocated in this section again)
  ALLOCATE ( nc_mjjj   (ntotrep) , STAT=istat )
  ALLOCATE ( nc_mmm    (ntotrep) , STAT=istat )
  ALLOCATE ( nc_myy    (ntotrep) , STAT=istat )
  ALLOCATE ( nc_mgg    (ntotrep) , STAT=istat )
  ALLOCATE ( nc_ngg    (ntotrep) , STAT=istat )
! ALLOCATE ( nc_msec   (ntotrep) , STAT=istat )

! get NetCDF ID's of date / time header variables
                            status = nf90_inq_varid (ncid, 'MJJJ',  varid_mjjj)
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MMM',   varid_mmm)
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MYY',   varid_myy)
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MGG',   varid_mgg)
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'NGG',   varid_ngg)
! IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MSEC',  varid_msec)
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'STANDARD DATE / TIME VARIABLES DO NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF

! get date / time values from NetCDF file
  status = nf90_get_var (ncid, varid_mjjj, nc_mjjj)
  status = nf90_get_var (ncid, varid_mmm , nc_mmm )
  status = nf90_get_var (ncid, varid_myy , nc_myy )
  status = nf90_get_var (ncid, varid_mgg , nc_mgg )
  status = nf90_get_var (ncid, varid_ngg , nc_ngg )
! status = nf90_get_var (ncid, varid_msec, nc_msec)

! get obs period for wind profiler
  ALLOCATE ( nc_dts    (ntotrep) , STAT=istat )
  nc_dts = 0
  status = nf90_inq_varid (ncid, 'MSETP', varid_dts)
  IF (status == nf90_noerr) status = nf90_get_var (ncid, varid_dts , nc_dts )

! Initial date and time of the forecast run
  READ (ydate_ref,'(I4,4I2)') imoyy, imomm, imodd, imohh, imomin

!   (these variables are de-allocated in this section again)
  ALLOCATE ( nt_date   (ntotrep) , STAT=istat )
  ALLOCATE ( nt_time   (ntotrep) , STAT=istat )
  ALLOCATE ( min_ob    (ntotrep) , STAT=istat )

  mrepsta = 0
  mrepend = -1

  loop_rep_time:  DO irep = 1 , ntotrep

! observation time
! ----------------

    IF (     (nc_mjjj(irep) < 0) .OR. (nc_mjjj(irep) > 2099)                   &
        .OR. (nc_mmm (irep) < 1) .OR. (nc_mmm (irep) > 12)                     &
        .OR. (nc_myy (irep) < 1) .OR. (nc_myy (irep) > 31)                     &
        .OR. (nc_mgg (irep) < 0) .OR. (nc_mgg (irep) > 24)                     &
        .OR. (nc_ngg (irep) < 0) .OR. (nc_ngg (irep) > 60)) THEN
      min_ob  (irep) = imdi
      CYCLE loop_rep_time
    ENDIF

    IF (nc_mjjj(irep) < 100) THEN
      nc_mjjj (irep) = 2000 + nc_mjjj(irep)
      IF (nc_mjjj(irep) > 2059)  nc_mjjj (irep) = nc_mjjj(irep) - 100
    ENDIF

    !    wind profiler: shift of obs time in minutes
    IF ((nc_dts(irep) == imiss) .OR. (itim_wp == 0))   nc_dts (irep) = 0
    IF (nc_dts(irep) /= 0) THEN
      IF (itim_wp == 1)  nc_dts (irep) =          nc_dts(irep) /120
      IF (itim_wp == 2)  nc_dts (irep) = MAX( 0 , nc_dts(irep) / 60 )
      IF (itim_wp == 3)  nc_dts (irep) = MIN( 0 , nc_dts(irep) / 60 )
    ENDIF

!CS
!   IF ((kcdftyp == ncdf_modes     ) .OR. (kcdftyp == ncdf_modes_acr )) THEN
!   nc_mjjj (irep) = 2009
!   nc_mmm  (irep) = 2
!   nc_myy  (irep) = 24
!!  nc_mgg  (irep) = 9
!   nc_mgg  (irep) = nc_mgg(irep) + 4
!   ENDIF
!   IF ((kcdftyp == ncdf_buoy      )                                  ) THEN
!   IF ((kcdftyp == ncdf_tempdesc  )                                  ) THEN
!     nc_mjjj (irep) = 2020
!     nc_mmm  (irep) = 08
!     nc_myy  (irep) = 03
!     nc_mgg  (irep) = 06
!     IF (nc_ngg(irep) >= 23)  nc_mgg  (irep) = 16
!!    nc_mgg  (irep) = nc_mgg(irep) + 4
!   ENDIF
!   IF ((kcdftyp == ncdf_tempdesc  )                                  ) THEN
!     nc_mjjj (irep) = 2022
!     nc_mmm  (irep) = 04
!     nc_myy  (irep) = 04
!     nc_mgg  (irep) = 11
!   ENDIF
!   IF ((kcdftyp == ncdf_temphirs  )                                  ) THEN
!     nc_mjjj (irep) = 2017
!     nc_mmm  (irep) = 10
!     nc_myy  (irep) = 20
!     nc_mgg  (irep) = 15
!!    nc_mgg  (irep) = nc_mgg(irep) + 4
!   ENDIF
!   IF ((kcdftyp == ncdf_satob     ) .OR. (kcdftyp == ncdf_pilot_p   )) THEN
!!    start date of my run: 20171020 1500 ; date of obs: 20171219 0430 - 0600
!!    nc_mjjj (irep) = 2017
!     nc_mmm  (irep) = nc_mmm(irep) -  2
!     nc_myy  (irep) = nc_myy(irep) +  1
!     nc_mgg  (irep) = nc_mgg(irep) + 11
!   ENDIF
!   IF (kcdftyp == ncdf_tower     ) THEN
!!    start date of my run: 20171020 1500 ; date of obs: 20171113 1500 - 1800
!     nc_mjjj (irep) = 2017
!     nc_mmm  (irep) = 10
!     nc_myy  (irep) = 20
!!    nc_mmm  (irep) = nc_mmm(irep) -  1
!!    nc_myy  (irep) = nc_myy(irep) +  7
!!    nc_mgg  (irep) = nc_mgg(irep) + 11
!   ENDIF
!CS

    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin                     &
                      , nc_mjjj(irep), nc_mmm(irep), nc_myy(irep)              &
                      , nc_mgg (irep), nc_ngg(irep), 0, imindif, ierrf )
!   =================
!   IF (irep < 10) &
!     PRINT   *, 'ZZ3 ', kcdftyp, irep, imindif, nc_mjjj(irep), nc_mmm(irep)   &
!                               , nc_myy(irep) , nc_mgg(irep) , nc_ngg(irep)

!   observation time in forecast minute units
    min_ob  (irep) = imindif + nc_dts(irep)

!   determine the interval of the record dimension which has to be read now
!   (for descending radiosonde reports, the obs time thresholds are extended
!    because their obs time is modified below)
    IF (    ((min_ob(irep) >= min_sta) .AND. (min_ob(irep) <= min_end))        &
        .OR.((kcdftyp == ncdf_tempdesc).AND. (min_ob(irep) >= min_sta-180)     &
                                       .AND. (min_ob(irep) <= min_end+30))) THEN
      IF (mrepsta < 1) mrepsta = irep
      mrepend = irep
    ENDIF
!  (else reject report if too old, and update report event counter)
!   WRITE( yobdat,'(I4,4I2.2)' )  nc_mjjj(irep), nc_mmm(irep), nc_myy(irep)    &
!                                              , nc_mgg(irep), nc_ngg(irep)
!   READ ( yobdat,'(I8,I4)' )  nt_date(irep), nt_time(irep)
    nt_date (irep)  =  10000* nc_mjjj(irep) + 100* nc_mmm(irep) + nc_myy(irep)
    nt_time (irep)  =                         100* nc_mgg(irep) + nc_ngg(irep)

    IF (nc_dts(irep) /= 0) THEN

      CALL get_utc_date ( min_ob(irep), ydate_ref, 60._wp, itype_calendar      &
                        , yakdat, yakdat2, nzday, zhr )
!     =================
      yobdat = yakdat(1:12)
      READ( yobdat,'(I8,I4)' )  nt_date(irep), nt_time(irep)
    ENDIF

  ENDDO loop_rep_time

  DEALLOCATE ( nc_mjjj  , STAT=istat )
  DEALLOCATE ( nc_mmm   , STAT=istat )
  DEALLOCATE ( nc_myy   , STAT=istat )
  DEALLOCATE ( nc_mgg   , STAT=istat )
  DEALLOCATE ( nc_ngg   , STAT=istat )
! DEALLOCATE ( nc_msec  , STAT=istat )
  DEALLOCATE ( nc_dts   , STAT=istat )

! get number of reports that should be read now
! ---------------------------------------------
  nrep = mrepend - mrepsta + 1

!-------------------------------------------------------------------------------
! Section 2c: For multi-level reports, mandatory multi-level variables are
!             allowed to be missing, if all the reports in the file have
!             zero observation levels. 
!             If a particular mandatory multi-level variable is missing,
!             the code therefore refrains from aborting. Instead, it writes
!             a Caution message and sets the number of reports 'nrep' to be
!             read and processed to zero.
!             The particular mandatory multi-level variable is chosen to be
!             equal to the one used to check the length of the vertical
!             dimension in 'src_obs_cdfin_mult.f90'.
!-------------------------------------------------------------------------------

  IF (nrep > 0) THEN
    status = nf90_noerr
    IF   (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)      &
          .OR. (kcdftyp == ncdf_pilot   ) .OR. (kcdftyp == ncdf_pilot_p )) THEN
      ymnem = 'NFNFN  '
      status = nf90_inq_varid (ncid, ymnem, varid_mand )
    ELSEIF (     (kcdftyp == ncdf_wprof ) .OR. (kcdftyp == ncdf_radar_vad)     &
            .OR. (kcdftyp == ncdf_rass  ) .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
      ymnem = 'MH     '
      status = nf90_inq_varid (ncid, ymnem, varid_mand )
    ELSEIF (     (kcdftyp == ncdf_tower ) ) THEN
      ymnem = 'MHOSEN '
      status = nf90_inq_varid (ncid, ymnem, varid_mand )
    ENDIF
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'MANDATORY VARIABLE "' // TRIM( ymnem ) //                     &
                '" DOES NOT EXIST IN ' // yfn
!     CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=istat)
      IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! ",A )' )  yerrmsl
      CLOSE (nucautn)
      !   set 'nrep' to zero as an indicator to skip the rest of the routine
      nrep = 0
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 2c: For descending TEMP,
!             replace the report time by the obs time of the lowest obs level,
!             re-determine the number 'nrep' of reports that should be read now,
!             and store the level index and dislocation of lowest obs level
!-------------------------------------------------------------------------------

  ALLOCATE ( isurfob  (MAX( nrep, 1 ) ) , STAT=istat )
  isurfob (:) = 0

  IF ((kcdftyp == ncdf_tempdesc) .AND. (nrep > 0)) THEN
    ! find level index of lowest obs level
                              status = nf90_inq_varid (ncid,'MLADH',varid_dlat)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MLODH',varid_dlon)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NLTPD',varid_dt  )
!   IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NHHHN',varid_z   )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MPN'  ,varid_p   )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MTDBT',varid_t   )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NFNFN',varid_ff  )
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'STANDARD dLAT / dLON VARIABLES DO NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
!     OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
!                                , POSITION='APPEND', IOSTAT=istat)
!     IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
!     IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
!     WRITE( nucautn,'("CAUTION !!!!! ",A )' )  yerrmsl
!     CLOSE (nucautn)
!     nrep = 0
    ENDIF
! ENDIF
! IF ((kcdftyp == ncdf_tempdesc) .AND. (nrep > 0)) THEN
    ! get length of vertical dimension in NetCDF file
    status = nf90_Inquire_Variable (ncid, varid_ff, ndims=ndims, dimids=dimids)
    dimid_mxlv = dimids(1)
    status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=maxlev)
    ALLOCATE ( rc_dlat  (maxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_dlon  (maxlev,nrep) , STAT=istat )
    ALLOCATE ( nc_dt    (maxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_p     (maxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_ff    (maxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_t     (maxlev,nrep) , STAT=istat )
    ALLOCATE ( isurflv         (nrep) , STAT=istat )
    isurflv (:) = 0
    nsta2 (1) = 1
    ncnt2 (1) = maxlev
    nsta2 (2) = mrepsta
    ncnt2 (2) = nrep
    status = nf90_get_var (ncid, varid_dlat , rc_dlat , nsta2, ncnt2)
    status = nf90_get_var (ncid, varid_dlon , rc_dlon , nsta2, ncnt2)
    status = nf90_get_var (ncid, varid_dt   , nc_dt   , nsta2, ncnt2)
    status = nf90_get_var (ncid, varid_p    , rc_p    , nsta2, ncnt2)
    status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
    status = nf90_get_var (ncid, varid_t    , rc_t    , nsta2, ncnt2)
    ! find level index 'isurflv' of lowest obs level, which is used to re-assign
    ! the time and position of the whole report;
    ! for 'isurflv', use a level with maximum information on time and position,
    ! i.e. if there are levels with time and position use such a level,
    ! else if there are levels with time only         use such a level
    ! else if there are levels with position only     use such a level
    ! else use a level without time and position;
    ! if there is no time and/or position at the selected 'isurflv' then the
    ! respective info related to the whole report will not be corrected;
    ! note that time is mandatory at 'isurflv' but lat/lon is not mandatory(!);
    ! if no such level is found (i.e. isurflv == isurfob == 0, because no level
    ! is found with temperature or wind at a given pressure) then the report
    ! will be discarded further below (kproc == 64) 
    ALLOCATE ( iexl (maxlev) , STAT=istat )
    DO irep = 1 , nrep
      iexl (:) = 0
      !   get info on time and position for each pressure level with T or wind:
      !     iexl(ilev) == 4 : time + position exist
      !                   3 : time exists
      !                   2 : position exists
      !                   1 : neither time nor position exist
      !                   0 : not even pressure and (T or wind) exist
      DO ilev = 1, maxlev
        IF (      (ABS( rc_p(ilev,irep) )                   < rmisschk)        &
            .AND. (     (ABS( rc_t (ilev,irep) ) < rmisschk)                   &
                   .OR. (ABS( rc_ff(ilev,irep) ) < rmisschk))) THEN
          IF (      (ABS( rc_dlat(ilev,irep) ) < rmisschk)                     &
              .AND. (ABS( rc_dlon(ilev,irep) ) < rmisschk)                     &
              .AND. (nc_dt(ilev,irep) /= imiss)) THEN
            iexl (ilev) = 4
          ELSEIF (nc_dt(ilev,irep) /= imiss) THEN
            iexl (ilev) = 3
          ELSEIF (      (ABS( rc_dlat(ilev,irep) ) < rmisschk)                 &
                  .AND. (ABS( rc_dlon(ilev,irep) ) < rmisschk)) THEN
            iexl (ilev) = 2
          ELSE
            iexl (ilev) = 1
          ENDIF
        ENDIF
      ENDDO
      !   determine maximum information on time and position
      iext = MAXVAL( iexl )
      pmax = -1._wp
      !   if there are levels with T or wind:
      IF (iext > 0) THEN
        DO ilev = 1, maxlev
          !   from all levels with max. info on time and position,
          !   get the level with max. pressure
          IF ((iexl(ilev) == iext) .AND. (rc_p(ilev,irep) > pmax)) THEN
            pmax           = rc_p(ilev,irep)
            isurflv (irep) = -ilev
          ENDIF
        ENDDO
      ENDIF
!     WRITE(0,*) my_cart_id, 'ZZq2 ',irep, isurflv(irep), iext                 &
!                                 , nc_dt(-isurflv(irep), irep)
    ENDDO
    DEALLOCATE ( iexl     , STAT=istat )
    DEALLOCATE ( rc_p     , STAT=istat )
    DEALLOCATE ( rc_ff    , STAT=istat )
    DEALLOCATE ( rc_t     , STAT=istat )
    ALLOCATE ( rc_dlatr  (nrep) , STAT=istat )
    ALLOCATE ( rc_dlonr  (nrep) , STAT=istat )
    rc_dlatr = 0.0_wp
    rc_dlonr = 0.0_wp

    ! re-determine 'min_ob', 'nt_date', 'nt_time[1:ntotrep]'
    mrepstao = mrepsta
    nrep     = mrepend - mrepstao + 1
    mrepsta  = 0
    mrepend  = -1
    DO irep = 1 , nrep
      ireptot = irep + mrepstao - 1
      ndsec   = nc_dt(MAX(-isurflv(irep),1),irep)
      IF ((isurflv(irep) /= 0) .AND. (ndsec /= imiss))                         &
        min_ob (ireptot) =  min_ob(ireptot) + NINT( REAL( ndsec, wp) /60._wp)

      CALL get_utc_date ( min_ob(ireptot), ydate_ref, 60._wp, itype_calendar   &
                        , yakdat, yakdat2, nzday, zhr )
!     =================
      yobdat = yakdat(1:12)
      READ( yobdat,'(I8,I4)' )  nt_date(ireptot), nt_time(ireptot)
      ! re-determine 'mrepsta', 'mrepend', and further below: 'nrep'
      IF (      (min_ob(ireptot) >= min_sta)                                   &
          .AND. (min_ob(ireptot) <= min_end)) THEN
        IF (mrepsta < 1) mrepsta = ireptot
        mrepend = ireptot
        irepn   = irep + mrepstao - mrepsta  !  = ireptot - mrepsta + 1
        ! re-map 'rc_dlat, rc_dlon , isurflv' (size: old value of nrep)
        !  onto 'rc_dlatr, rc_dlonr, isurfob' (size: new value of nrep as below)
        IF (isurflv(irep) /= 0) THEN
          ilev = - isurflv(irep)
          IF (      (ABS( rc_dlat(ilev,irep) ) < rmisschk)                     &
              .AND. (ABS( rc_dlon(ilev,irep) ) < rmisschk)) THEN
            rc_dlatr (irepn)  =  rc_dlat(ilev,irep)
            rc_dlonr (irepn)  =  rc_dlon(ilev,irep)
          ENDIF
          isurfob  (irepn)  =  -ilev
      ! else if (isurflv == 0) kproc = kproc + 64, and report discarded later on
        ENDIF
      ENDIF
    ENDDO
    DEALLOCATE ( nc_dt    , STAT=istat )
    DEALLOCATE ( rc_dlat  , STAT=istat )
    DEALLOCATE ( rc_dlon  , STAT=istat )
    DEALLOCATE ( isurflv  , STAT=istat )

    nrep = mrepend - mrepsta + 1
  ENDIF

!-------------------------------------------------------------------------------
! Section 2d: Get the indices (related to the NetCDF file record dimension)
!             of the reports to be read now
!-------------------------------------------------------------------------------

! if number of reports to be read now is zero then clean up and return
! --------------------------------------------------------------------
! PRINT *, 'obs file type , number of reports ', yroutine, '  ', kcdftyp, nrep

  IF (nrep <= 0) THEN
    DEALLOCATE ( nt_date  , STAT=istat )
    DEALLOCATE ( nt_time  , STAT=istat )
    DEALLOCATE ( min_ob   , STAT=istat )
    DEALLOCATE ( isurfob  , STAT=istat )
    IF (kcdftyp == ncdf_tempdesc) THEN
      DEALLOCATE ( rc_dlatr , STAT=istat )
      DEALLOCATE ( rc_dlonr , STAT=istat )
    ENDIF
    nrep  =  0
    PRINT '(" processing",I6," reports from",I5," -",I5," [min] from file ",A)'&
          ,nrep, min_sta, min_end, yfn(1:LEN_TRIM(yfn))
    RETURN
!   ~~~~~~
  ENDIF

! get the chunk of the record dimension
! which contains those reports that should be read now
! ----------------------------------------------------

!   (these variables are de-allocated in 'obs_cdf_buffer_comhead'
!    (if (nrep>0).and.(my_cart_id==0)))
  ALLOCATE ( zr_hour  (nrep) , STAT=istat )
  ALLOCATE ( nr_date  (nrep) , STAT=istat )
  ALLOCATE ( nr_time  (nrep) , STAT=istat )

  DO irep = 1 , nrep
    zr_hour (irep)  =  min_ob (irep+mrepsta-1) / 60.0_wp
    nr_date (irep)  =  nt_date(irep+mrepsta-1)
    nr_time (irep)  =  nt_time(irep+mrepsta-1)
  ENDDO

! get the vector of indices (in the interval 1...nrep) of reports
! for which the observation time lies within the wanted interval
! and which therefore have to be processed now
! ---------------------------------------------------------------

!   ('irproc' is de-allocated in 'obs_cdf_read_*' (if nrep>0, my_cart_id=0))
  ALLOCATE ( irproc    (nrep) , STAT=istat )
!   (this variable is de-allocated at the end of this subroutine again)
  ALLOCATE ( kproc     (nrep) , STAT=istat )
  irproc  (:) = 0
  nreproc     = 0

  DO irep = 1 , nrep
    IF (      (min_ob(irep+mrepsta-1) >= min_sta)                              &
        .AND. (min_ob(irep+mrepsta-1) <= min_end)) THEN
      nreproc = nreproc + 1
      irproc (nreproc)  =  irep
! (kproc(irep) == 0) means that the report has to be processed further
      kproc (irep)      =  0
    ELSE
      kproc (irep)      =  1
    ENDIF
  ENDDO
! PRINT *, 'comnrep0 ', nrep, nreproc

  DEALLOCATE ( nt_date  , STAT=istat )
  DEALLOCATE ( nt_time  , STAT=istat )
  DEALLOCATE ( min_ob   , STAT=istat )

! allocate header variables which are read / derived in this routine
! ------------------------------------------------------------------
!   (these variables are de-allocated in 'obs_cdf_buffer_comhead')
!    (if (nrep>0).and.(my_cart_id==0)))

  ALLOCATE ( kflag    (nrep) , STAT=istat )
  ALLOCATE ( nsynhr   (nrep) , STAT=istat )
  ALLOCATE ( kcat     (nrep) , STAT=istat )
  ALLOCATE ( kcatsub  (nrep) , STAT=istat )
  ALLOCATE ( iobtyp   (nrep) , STAT=istat )
  ALLOCATE ( icdtyp   (nrep) , STAT=istat )
  ALLOCATE ( istidn   (nrep) , STAT=istat )
  ALLOCATE ( ystidn   (nrep) , STAT=istat )
  ALLOCATE ( kcentre  (nrep) , STAT=istat )
  ALLOCATE ( kupdate  (nrep) , STAT=istat )
  ALLOCATE ( kz_dwd   (nrep) , STAT=istat )
  ALLOCATE ( iobstot  (nrep) , STAT=istat )
  ALLOCATE ( jobstot  (nrep) , STAT=istat )
  ALLOCATE ( mrepflg  (nrep) , STAT=istat )
  ALLOCATE ( iobsloc  (nrep) , STAT=istat )
  ALLOCATE ( jobsloc  (nrep) , STAT=istat )
! ALLOCATE ( isurfob  (nrep) , STAT=istat )

  ALLOCATE ( zsynhr   (nrep) , STAT=istat )
  ALLOCATE ( rc_lat   (nrep) , STAT=istat )
  ALLOCATE ( rc_lon   (nrep) , STAT=istat )
  ALLOCATE ( zaltob   (nrep) , STAT=istat )
  ALLOCATE ( zaltmo   (nrep) , STAT=istat )
  ALLOCATE ( rio_tot  (nrep) , STAT=istat )
  ALLOCATE ( rjo_tot  (nrep) , STAT=istat )
  ALLOCATE ( ztdecdb  (nrep) , STAT=istat )

  DO irep = 1 , nrep
    kflag   (irep) =  0
    kcat    (irep) = imdi
    kcatsub (irep) = imdi
    iobtyp  (irep) = -9
    icdtyp  (irep) = -9
    istidn  (irep) = imdi
    ystidn  (irep) = ' '
    kcentre (irep) = imdi
    kupdate (irep) = imdi
    kz_dwd  (irep) = imdi
    zaltob  (irep) = rmdi
    zaltmo  (irep) = rmdi
    iobstot (irep) = imdi
    jobstot (irep) = imdi
!   isurfob (irep) = 0
    mrepflg (irep) = 0
  ENDDO

!-------------------------------------------------------------------------------
! Section 2e: Get nominal (synoptic) observation time
!-------------------------------------------------------------------------------

!   (these variables are de-allocated in this section again)
  ALLOCATE ( nsyndat  (nrep) , STAT=istat )
  ALLOCATE ( nsyntim  (nrep) , STAT=istat )

  status   = nf90_inq_varid (ncid,'section1_date' ,varid_s1date )
  IF (status == nf90_noerr)                                                    &
    status = nf90_inq_varid (ncid,'section1_time' ,varid_s1time )
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'section1_date /-_time DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF

  status = nf90_get_var (ncid, varid_s1date, nsyndat, (/mrepsta/), (/nrep/))
  status = nf90_get_var (ncid, varid_s1time, nsyntim, (/mrepsta/), (/nrep/))

  DO irep = 1 , nrep
    IF ((nsyndat(irep) /= imiss) .AND. (nsyntim(irep) /= imiss)) THEN
      nsynmin = nsyntim(irep) / 100

!CS
!     IF (     (kcdftyp == ncdf_tempdesc) .OR. (kcdftyp == ncdf_temphirs)) THEN
!       nsyndat (irep) = nr_date(irep)
!       nsyntim (irep) = nr_time(irep)
!     ENDIF
!CS

      WRITE( yobdat,'(I8,I4.4)' )  nsyndat(irep), nsynmin
      READ ( yobdat,'(I4,4I2)' )  ianyy, ianmm, iandd, ianhh, ianmin

!     Time difference in minutes between observation time and initial model time
      CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                  &
                          ianyy, ianmm, iandd, ianhh, ianmin, 0, imindif, ierrf)
!     =================
!     nominal (synoptic) observation time in [yymmddhh]
      nsynhr  (irep)  =  MOD( nsyndat(irep), 1000000 ) *100                    &
                        +     nsyntim(irep)/ 10000
!     nominal (synoptic) observation time in forecast hour units
      zsynhr  (irep)  =  imindif / 60.0_wp
    ELSE
      nsynhr  (irep)  =  imdi
      zsynhr  (irep)  =  rmdi
    ENDIF
  ENDDO
  DEALLOCATE ( nsyndat  , STAT=istat )
  DEALLOCATE ( nsyntim  , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2f: Prepare check of time of last known position for buoys (section 6)
!-------------------------------------------------------------------------------

!   (this variable is de-allocated in section 6 again)
  ALLOCATE ( min_dtl   (nrep) , STAT=istat )
  min_dtl = 0

  IF (kcdftyp == ncdf_buoy) THEN
                              status = nf90_inq_varid (ncid, 'MJJJ0',varid_mjjj)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MMM0', varid_mmm)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MYY0', varid_myy)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'MGG0', varid_mgg)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid, 'NGG0', varid_ngg)
!   IF (status /= nf90_noerr) THEN
!     yerrmsl = 'POSITION DATE / TIME VARIABLES DO NOT EXIST IN ' // yfn
!     CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
!   ENDIF
    IF (status == nf90_noerr) THEN

      !   (these variables are de-allocated in this section again)
      ALLOCATE ( nc_mjjj   (nrep) , STAT=istat )
      ALLOCATE ( nc_mmm    (nrep) , STAT=istat )
      ALLOCATE ( nc_myy    (nrep) , STAT=istat )
      ALLOCATE ( nc_mgg    (nrep) , STAT=istat )
      ALLOCATE ( nc_ngg    (nrep) , STAT=istat )
      status = nf90_get_var (ncid, varid_mjjj, nc_mjjj, (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_mmm , nc_mmm , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_myy , nc_myy , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_mgg , nc_mgg , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_ngg , nc_ngg , (/mrepsta/), (/nrep/))

      DO irpa = 1 , nreproc
        irep = irproc(irpa)

        IF (      (nc_mjjj(irep) >= 0) .AND. (nc_mjjj(irep) <= 2099)           &
            .AND. (nc_mmm (irep) >= 1) .AND. (nc_mmm (irep) <= 12)             &
            .AND. (nc_myy (irep) >= 1) .AND. (nc_myy (irep) <= 31)             &
            .AND. (nc_mgg (irep) >= 0) .AND. (nc_mgg (irep) <= 24)             &
            .AND. (nc_ngg (irep) >= 0) .AND. (nc_ngg (irep) <= 60)) THEN
          IF (nc_mjjj(irep) < 100) THEN
            nc_mjjj (irep) = 2000 + nc_mjjj(irep)
            IF (nc_mjjj(irep) > 2059)  nc_mjjj (irep) = nc_mjjj(irep) - 100
          ENDIF

          !   difference in minutes between location time and initial model time
          CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin               &
                            , nc_mjjj(irep), nc_mmm(irep), nc_myy(irep)        &
                            , nc_mgg (irep), nc_ngg(irep), 0, imindif, ierrf )
!         =================

          !   observation time minus location time
          min_dtl (irep)  =  NINT( zr_hour(irep) * 60.0_wp )  -  imindif

        ENDIF
      ENDDO
      DEALLOCATE ( nc_mjjj  , STAT=istat )
      DEALLOCATE ( nc_mmm   , STAT=istat )
      DEALLOCATE ( nc_myy   , STAT=istat )
      DEALLOCATE ( nc_mgg   , STAT=istat )
      DEALLOCATE ( nc_ngg   , STAT=istat )
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Get observation and code type
!            from data category and sub category (BUFR 4, Section 1)
!-------------------------------------------------------------------------------
!  data category:                       data sub category
!  0: land surface:                     0: SYNOP (hourly fixed land)
!                                       1: SYNOP (intermediate fixed land)
!                                       2: SYNOP (main synoptic fixed land)
!                                       3: SYNOP MOBIL (hourly)
!                                       4: SYNOP MOBIL (intermediate)
!                                       5: SYNOP MOBIL (main synoptic)
!                                       6: 1-hour obs from automated stations
!                                       7: n-minute obs from AWS stations
!                                      10: METAR (routine aeronautical obs)
!                                      11: SPECI (special aeronautical obs)
!  1: sea surface:                      0: SHIP (synoptic observations)
!                                       6: 1-hour obs from automated stations
!                                       7: n-minute obs from AWS stations
!                                      25: BUOY
!  2: vertical soundings (no sat):      1: PILOT (fixed land)
!                                       2: PILOT SHIP
!                                       3: PILOT MOBIL
!                                       4: TEMP (fixed land)
!                                       5: TEMP SHIP
!                                       6: TEMP MOBIL
!                                       7: TEMP DROP
!                                      10: Wind Profiler
!                                      11: RASS temperature profiles
!                                      20: AMDAR / ACARS profiles
!                       not yet defined  : Tower (Lindenberg)
!  3: satellite soundings:              0: temperature (SATEM)
!                                       1: ATOVS
!  4: single-level upper-air (no sat):  0: AMDAR / ACARS
!                                       1: Manual (AIREP, PIREP)
!  5: single-level upper-air satellite: 0: cloud wind data (AMV)
!  6: radar:                            0: reflectivity
!                                       1: Doppler wind profiles
!                                       2: derived products
!                                       3: ground radar weather (RADOB)
! 12: satellite surface data:           4: SSM/I radiometer
!                                       5: Quickscat
!                                       6: surface temperat./radiation (SATOB)
!-------------------------------------------------------------------------------

  ALLOCATE ( nedno    (nrep) , STAT=istat )

  status   = nf90_inq_varid (ncid,'section1_data_category'        ,varid_s1cat )
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'section1_data_category DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF
    status = nf90_inq_varid (ncid,'section1_int_data_sub_category',varid_s1cats)
  IF (status /= nf90_noerr)                                                    &
    status = nf90_inq_varid (ncid,'section1_data_sub_category'    ,varid_s1cats)
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'section1_data_(sub)_category DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF

  status = nf90_get_var (ncid, varid_edno  , nedno  , (/mrepsta/), (/nrep/))
  status = nf90_get_var (ncid, varid_s1cat , kcat   , (/mrepsta/), (/nrep/))
  status = nf90_get_var (ncid, varid_s1cats, kcatsub, (/mrepsta/), (/nrep/))
! manned (nc_nix = 1) or automated (nc_nix = 0) station ?
!   (this variable is de-allocated in this section again)
  ALLOCATE ( nc_nix   (nrep) , STAT=istat )
  ALLOCATE ( rc_p   (1,nrep) , STAT=istat )
  DO irep = 1 , nrep
    nc_nix  (irep) =  imdi
    rc_p  (1,irep) =  rmiss
  ENDDO
                               status = nf90_inq_varid (ncid,'NIX  ',varid_nix )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_nix   , nc_nix , (/mrepsta/), (/nrep/))
                               status = nf90_inq_varid (ncid,'MPPP ',varid_p   )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_p     , rc_p   , (/mrepsta/), (/nrep/))

! determine observation and code types
  DO irep = 1 , nrep
    kcdtyp  =  icdtyp(irep)

    ! create value for 'int_data_sub_category' for BUFR edition 3 data
    ! because this entry must be set only for BUFR edition 4 ff
    ! (some sub-categories (e.g. TEMP MOBILE) cannot be identified
    !  and are assigned to the most common sub-category (e.g. TEMP))
    IF ((nedno(irep) <= 3) .OR. (nedno(irep) > 9)) THEN
      IF (kcdftyp == ncdf_amdar     )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_acars     )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_acars_us  )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_acars_uk  )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_modes     )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_modes_acr )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_amdar_ml  )  kcatsub (irep) = 20
      IF (kcdftyp == ncdf_amdar_vp  )  kcatsub (irep) = 20
      IF (kcdftyp == ncdf_buoy      )  kcatsub (irep) = 25
      IF (kcdftyp == ncdf_ship      )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_synop     )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_synop_mob )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_synop_tst )  kcatsub (irep) = 0
      IF (kcdftyp == ncdf_temp      )  kcatsub (irep) = 4
      IF (kcdftyp == ncdf_temphirs  )  kcatsub (irep) = 4
      IF (kcdftyp == ncdf_tempship  )  kcatsub (irep) = 5
      IF (kcdftyp == ncdf_tempdesc  )  kcatsub (irep) = 14
      IF (kcdftyp == ncdf_tempdrop  )  kcatsub (irep) = 7
      IF (kcdftyp == ncdf_pilot     )  kcatsub (irep) = 1
      IF (kcdftyp == ncdf_pilot_p   )  kcatsub (irep) = 1
      IF (kcdftyp == ncdf_wprof     )  kcatsub (irep) = 10
      IF (kcdftyp == ncdf_rass      )  kcatsub (irep) = 11
      IF (kcdftyp == ncdf_radar_vad )  kcatsub (irep) = 1
      IF (kcdftyp == ncdf_wlidar_wp )  kcatsub (irep) = 10
      IF (kcdftyp == ncdf_metar     )  kcatsub (irep) = 10
      ! GPS ZPD (ZTD) data are obtained from UKMO as BUFR edition 3 to date
      IF (kcdftyp == ncdf_gps_zenith)  kcatsub (irep) = 14
      IF (kcdftyp == ncdf_ascat     )  kcatsub (irep) = 7
      IF (kcdftyp == ncdf_qscat     )  kcatsub (irep) = 5

      ! before 1 August 2008, AMDAR data is BUFR edition 3 with erroneous
      ! 'section1_data_category == 2' instead of == 4 in the DWD data base
      ! --> correct this to avoid problems in re-analysis
      IF ((kcdftyp == ncdf_amdar) .AND. (kcat(irep) == 2))  kcat (irep) = 4

!CS: begin DWD data base bug fixes
    ! between 1 - 24 July 2008 there are some BUFR edition 4 data with erroneous
    ! 'section1_int_data_sub_category == 255' in the DWD data base
    ! --> correct this to avoid problems in re-analysis
    ! ELSEIF (kcatsub(irep) == 255) THEN
    !   IF (kcdftyp == ncdf_acars     )  kcatsub(irep) = 0
    !   IF (kcdftyp == ncdf_wprof     )  kcatsub(irep) = 10
    !   IF (kcdftyp == ncdf_rass      )  kcatsub(irep) = 11
    !   IF (kcdftyp == ncdf_radar_vad )  kcatsub(irep) = 1
    ! some VADs had wrong 'section1_int_data_sub_category' (0 instead of 1)
    ! in the DWD data base (e.g. in 2009)
    ! --> correct this to avoid problems in re-analysis
    ELSEIF (kcdftyp == ncdf_radar_vad) THEN
      kcatsub(irep) = 1
    ! RASS can have wrong 'section1_int_data_sub_category' (10 instead of 11)
    ! if bufrxform chain is applied to combine wind profiler - RASS reports
    ! --> correct this to avoid problems in re-analysis
    ELSEIF (kcdftyp == ncdf_rass) THEN
      IF (kcatsub(irep) == 10)  kcatsub(irep) = 11
    ELSEIF (kcdftyp == ncdf_tower) THEN
      kcatsub(irep) = 75
    ELSEIF (kcdftyp == ncdf_tempdesc) THEN
      kcatsub(irep) = 14
!CS: end DWD data base bug fixes
!CS: begin COMET data base bug fixes
    ELSEIF (kcdftyp == ncdf_satob) THEN
      IF (kcatsub(irep) == imiss) kcatsub(irep) = 0
!CS: end COMET data base bug fixes
    ENDIF
!   IF (kcdftyp == ncdf_ascat     )  kcatsub (irep) = 7
!   IF (kcdftyp == ncdf_qscat     )  kcatsub (irep) = 5

    IF (kcat(irep) == 0) THEN
      iobtyp (irep) = nsynop
      IF ((kcatsub(irep) >= 0) .AND. (kcatsub(irep) <= 5))   kcdtyp = nsrscd
      IF ((kcatsub(irep) == 6)  .OR. (kcatsub(irep) == 7))   kcdtyp = natscd
!     IF  (kcatsub(irep) == 10)                              kcdtyp = nmetar
      IF ((kcatsub(irep) == 10) .OR. (kcatsub(irep) == 193)) kcdtyp = nmetar
      IF ((nc_nix(irep) == 0) .AND. (kcdtyp == nsrscd))      kcdtyp = natscd
      IF ((nc_nix(irep) == 1) .AND. (kcdtyp == natscd))      kcdtyp = nsrscd
      IF (kcdftyp == ncdf_synop_tst)                         kcdtyp = nsytst
      IF (kcatsub(irep) == 14) THEN
        iobtyp (irep) = ngps
        ! code type for yet unknown processing centre
        kcdtyp        = gpc(n_gpc) %cdtyp
      ENDIF
    ELSEIF (kcat(irep) == 1) THEN
      IF (kcatsub(irep) == 25) THEN
        iobtyp (irep) = ndribu
        kcdtyp        = ndrbcd
      ELSE
        iobtyp (irep) = nsynop
!       IF  (kcatsub(irep) == 0)                             kcdtyp = nshscd
        IF ((kcatsub(irep) == 0) .OR. (kcatsub(irep) == 1))  kcdtyp = nshscd
        IF ((kcatsub(irep) == 6) .OR. (kcatsub(irep) == 7))  kcdtyp = natshs
        IF ((nc_nix(irep) == 0) .AND. (kcdtyp == nshscd))    kcdtyp = natshs
        IF ((nc_nix(irep) == 1) .AND. (kcdtyp == natshs))    kcdtyp = nshscd
      ENDIF
    ELSEIF (kcat(irep) == 2) THEN
      IF (     ((kcatsub(irep) >= 4) .AND. (kcatsub(irep) <= 7))               &
          .OR. (kcatsub(irep) == 14) .OR.  (kcatsub(irep) == 213)) THEN
        iobtyp (irep) = ntemp
        IF (kcatsub(irep) == 4)  kcdtyp = nldtcd
        IF (kcatsub(irep) == 5)  kcdtyp = nshtcd
        IF (kcatsub(irep) == 6)  kcdtyp = nmotcd
        IF (kcatsub(irep) == 7)  kcdtyp = ntdrop
        IF (kcatsub(irep) == 14) kcdtyp = ntdesc
        IF (kcatsub(irep) == 213) kcdtyp = nldtcd ! (where does this come from?)
        ! high-res TEMP can only be distinguished from TEMP by file name (or KZ)
        IF (kcdftyp == ncdf_temphirs) THEN
          IF (kcatsub(irep) == 4)  kcdtyp = nbtemp
          IF (kcatsub(irep) == 5)  kcdtyp = nbship
          IF (kcatsub(irep) == 7)  kcdtyp = nbdrop
        ENDIF
      ELSEIF ((kcatsub(irep) >= 1) .AND. (kcatsub(irep) <= 3)) THEN
        iobtyp (irep) = npilot
        IF (kcatsub(irep) == 1)  kcdtyp = nldpcd  ! wind profiler in BUFR 3 !
        IF (kcatsub(irep) == 2)  kcdtyp = nshpcd
        IF (kcatsub(irep) == 3)  kcdtyp = nmopcd
        IF (kcatsub(irep) == 211)  kcdtyp = nldpcd
      ELSEIF ((kcatsub(irep) == 10) .OR. (kcatsub(irep) == 11)) THEN
        iobtyp (irep) = npilot
        IF (kcatsub(irep) == 10) kcdtyp = nwp_eu
        IF (kcatsub(irep) == 11) kcdtyp = nra_eu
        IF (kcdftyp == ncdf_wlidar_wp) kcdtyp = nwlidr
      ELSEIF (kcatsub(irep) == 20) THEN
        iobtyp (irep) = nairep
!CSC: to be modified (multi-level ACARS denoted as CODAR)
        kcdtyp        = ncodar
      ELSEIF (kcatsub(irep) == 75) THEN
        iobtyp (irep) = npilot
        kcdtyp        = ntower
        !   if surface pressure obs present then assign ICOS tower
        IF (ABS( rc_p(1,irep) ) < rmisschk)  kcdtyp = ntowic
      ENDIF
    ELSEIF (kcat(irep) == 3) THEN
      IF (kcatsub(irep) == 0) THEN
        iobtyp (irep) = nsatem
        kcdtyp        = nstmcd
      ELSEIF (kcatsub(irep) == 1) THEN
        iobtyp (irep) = nsattv
        kcdtyp        = nstovs
      ENDIF
    ELSEIF (kcat(irep) == 4) THEN
      iobtyp (irep) = nairep
      IF (kcatsub(irep) == 0)        kcdtyp = namdar
      IF (kcdftyp == ncdf_acars)     kcdtyp = nacar
      IF (kcdftyp == ncdf_acars_uk)  kcdtyp = nacar
      IF (kcdftyp == ncdf_acars_us)  kcdtyp = nacar
      IF (kcdftyp == ncdf_modes)     kcdtyp = nmodes
      IF (kcdftyp == ncdf_modes_acr) kcdtyp = nmodes
      IF (kcatsub(irep) == 1)        kcdtyp = naircd
!CSC: to be modified as soon as WMO defines ACARS data_sub_categories
      IF (kcatsub(irep) == imiss) kcdtyp = nacar
    ELSEIF (kcat(irep) == 5) THEN
      iobtyp (irep) = nsatob
!     IF (kcatsub(irep) == 0)  kcdtyp = nstbcd
      IF (kcatsub(irep) == 0)  kcdtyp = namv
    ELSEIF (kcat(irep) == 6) THEN
      iobtyp (irep) = npilot
      IF (kcatsub(irep) ==  1)  kcdtyp = nravad
      IF (kcatsub(irep) == 10)  kcdtyp = nwp_eu
    ELSEIF (kcat(irep) == 12) THEN
      iobtyp (irep) = nscatt
      IF (kcatsub(irep) == 5)  kcdtyp = nqscat
      IF (kcatsub(irep) == 7)  kcdtyp = nascat
    ENDIF
    icdtyp (irep) = kcdtyp
    IF ((iobtyp(irep) == -9) .OR. (kcdtyp == -9)) kproc(irep) = kproc(irep) + 2
  ENDDO
  DEALLOCATE ( nedno  , STAT=istat )
  DEALLOCATE ( nc_nix , STAT=istat )
  DEALLOCATE ( rc_p   , STAT=istat )

! write alerting CAUTION messages onto file YUCAUTN if unknown obs types occur
! ----------------------------------------------------------------------------

  nerrobt = 0
  DO irep = 1 , nrep
    IF ((MOD(kproc(irep),4) >= 2) .AND. (MOD(kproc(irep), 2) < 1)) THEN
      nerrobt = nerrobt + 1
      ireperr = irep
    ENDIF
  ENDDO
  IF (nerrobt > 0) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
    WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," TIMES WRONG DATA (SUB-)" &
                   &," CATEGORY IN FILE:",A)' )  acthr, nerrobt, yfn
    WRITE( nucautn,'("CAUTION !!!!!   ",6X,"  E.G. DATA CATEGORY:",I11         &
                   &,", SUB CATEGORY:",I11)' )  kcat(ireperr), kcatsub(ireperr)
    WRITE( nucautn,'("        ==>  CHECK   NetCDF - FILE ",A)' ) yfn
    PRINT          '("CAUTION !!!!! t=",F6.3,":",I5," TIMES WRONG DATA (SUB-)" &
                   &," CATEGORY IN FILE:",A)' ,  acthr, nerrobt, yfn
    PRINT          '("CAUTION !!!!!   ",6X,"  E.G. DATA CATEGORY:",I11         &
                   &,", SUB CATEGORY:",I11)' ,  kcat(ireperr), kcatsub(ireperr)
    PRINT          '("        ==>  CHECK   NetCDF - FILE ",A)' , yfn
    CLOSE (nucautn)
  ENDIF

! re-determine the vector of indices of reports which have to be processed now
! (observation time within wanted interval, known data (sub-)category)
! ----------------------------------------------------------------------------

  nreproc = 0
  DO irep = 1 , nrep
    IF (kproc(irep) == 0) THEN
      nreproc = nreproc + 1
      irproc (nreproc)  =  irep
    ENDIF
  ENDDO
! PRINT *, 'comnrep2 ', nrep, nreproc

!-------------------------------------------------------------------------------
! Section 4: Get originating center, sub center, and update sequence number
!            (correction indicator) from (BUFR) Section 1,
!            and DWD data base 'Kennzahl' (internal classification number)
!            and data base report decoding time from (BUFR) Section 2.
!-------------------------------------------------------------------------------
!  originating center                  sub center
!  74: UKMO, Exeter (RSMC)             21: Agenzia Spaziale Italiana       (I)
!                                      22: CNRS                            (F)
!                                      23: GFZ Potsdam                     (D)
!                                      24: Geodetic Observatory Pecny      (CZ)
!                                      25: Inst. Estud. Espaci. Catalunya  (E)
!                                      26: Swiss Topo                      (CH)
!                                      27: Nordic Commission of Geodesy    (N)
!                                      28: Nordic Commission of Geodesy    (S)
!                                      29: Inst. Geodesie National         (F)
!                                      30: B.amt Kartographie u. Geodaesie (D)
!                                      31: Inst. Engin. Sat. Survey Geodesy(UK)
!                                      33:
!                                      34:
!  76: Moscow (RSMC)
!  78: Offenbach (RSMC)
!  80: Rome (RSMC)
!  82: Norrkoping
!  84, 85: Toulouse (RSMC)
!  86: Helsinki                       some other centers:
!  87: Belgrade                       210: Frascati (ESA / ESRIN)
!  88: Oslo                           211: Lannion
!  89: Prague                         214: Madrid
!  90: Episkopi                       215: Zurich
!  91: Ankara                         216: Service ARGOS Toulouse
!  92: Frankfurt /Main                220: Warsaw
!  93: London (WAFC)                  242: Romania NMC
!  94: Copenhagen                      10: Cairo (RSMC)       20: Las Palmas
!  95: Rota                            16: Casablanca (RSMC)  21: Algiers (RSMC)
!  96: Athens                          34: Tokyo (RSMC), JMA  38: Beijing (RSMC)
!  97: European Space Agency (ESA)      7: NCEP
!  98: ECMWF, RSMC                      8: US NWS Telecomm. Gateway (NWSTG)
!  99: De Bilt                          9: US NWS (other)
!  160: US NOAA / NESDIS               56: ARINC Centre
!  161: US NOAA Oce. Atmos. Research   57: US Air Force
!  173: US NASA                        59: NOAA NFL, Boulder
!  254: EUMETSAT Operation Centre      60: NCAR
!  65535: missing
!-------------------------------------------------------------------------------

!   (this variable is de-allocated in this section again)
  ALLOCATE ( kcensub  (nrep) , STAT=istat )
  ALLOCATE ( koricen  (nrep) , STAT=istat )

  kcensub = imdi
  koricen = imdi

  status   = nf90_inq_varid (ncid,'section1_centre'    ,varid_s1cent )
  IF (status == nf90_noerr)                                                    &
    status = nf90_inq_varid (ncid,'section1_subcentre' ,varid_s1cents)
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'section1_(sub)centre DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF

  status = nf90_get_var (ncid, varid_s1cent , kcentre, (/mrepsta/), (/nrep/))
  status = nf90_get_var (ncid, varid_s1cents, kcensub, (/mrepsta/), (/nrep/))

  loricen = (nf90_inq_varid (ncid,'MMIOGC',varid_oricen) == nf90_noerr)
  IF ((.NOT. loricen) .AND. (kcdftyp == ncdf_modes_acr))                       &
    loricen = (nf90_inq_varid (ncid,'MOC11',varid_oricen) == nf90_noerr)
  IF (loricen)                                                                 &
    status = nf90_get_var (ncid, varid_oricen, koricen, (/mrepsta/), (/nrep/))

! if 'sub-centre' is given (for ground-based GPS, this should indicate the
! processing centre), multiply by 1000 and store it also in the 'centre' word
  DO irpa = 1 , nreproc
    irep = irproc(irpa)
    IF (      (iobtyp(irep) == ngps)                                           &
        .AND. (kcentre(irep) /= imiss) .AND. (kcentre(irep) < 655)             &
        .AND. (kcensub(irep) /= imiss) .AND. (kcensub(irep) < 100)) THEN
      kcentre (irep)  =  kcentre(irep)  +  1000* kcensub(irep)
!     IF (kcensub(irep) == 21) kcdtyp = ngpasi
!     IF (kcensub(irep) == 23) kcdtyp = ngpgfz
!     IF (kcensub(irep) == 24) kcdtyp = ngpgop
!     IF (kcensub(irep) == 25) kcdtyp = ngpieec
!     IF (kcensub(irep) == 26) kcdtyp = ngplpt
!     IF (kcensub(irep) == 29) kcdtyp = ngpsgn
!     IF (kcensub(irep) == 30) kcdtyp = ngpbkg
!     IF((kcensub(irep) == 32) .OR. (kcensub(irep) == 37)) kcdtyp = ngprob
!     IF (kcensub(irep) == 33) kcdtyp = ngpknmi
!     IF (kcensub(irep) == 34) kcdtyp = ngpnga
!     IF (kcensub(irep) == 35) kcdtyp = ngpige
!     IF (kcensub(irep) ==  0) kcdtyp = ngpmet
!     icdtyp (irep)  =  kcdtyp

! if 'original centre' is given (for standarised ACARS)
! multiply by 1000 and store it also in the 'centre' word
    ELSEIF (      (loricen)                                                    &
            .AND. (kcentre(irep) /= imiss) .AND. (kcentre(irep) < 655)         &
            .AND. (koricen(irep) /= imiss) .AND. (koricen(irep) < 655)) THEN
      kcentre (irep)  =  kcentre(irep)  +  1000* koricen(irep)
    ENDIF
  ENDDO

  DEALLOCATE ( koricen , STAT=istat )

! update sequence number (used as correction indicator)
! ----------------------
  status = nf90_inq_varid (ncid,'section1_update_sequence_nr',varid_s1updat)
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'section1_update_sequence_nr DOES NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF
  IF (status == nf90_noerr) THEN
    status = nf90_get_var (ncid, varid_s1updat, kupdate, (/mrepsta/), (/nrep/))
  ENDIF

! DWD-Kennzahl (internal observation classification number from DWD data base)
! ------------
  status = nf90_inq_varid (ncid,'section2_ikz'        ,varid_s2ikz  )
  IF (status == nf90_noerr) THEN
    status = nf90_get_var (ncid, varid_s2ikz  , kz_dwd , (/mrepsta/), (/nrep/))
  ENDIF

! data base decoding time (i.e. time when a report has been decoded and written
! -----------------------       to the data base; this is required to mimic
!                               (operational) data cut-off times)

!   (these variables are de-allocated in this section again)
  ALLOCATE ( ndbddat  (nrep) , STAT=istat )
  ALLOCATE ( ndbdtim  (nrep) , STAT=istat )
  DO irep = 1 , nrep
    ndbddat (irep) = imiss
    ndbdtim (irep) = imiss
  ENDDO

  status = nf90_inq_varid (ncid,'section2_decoding_date' ,varid_s2dbdd )
  IF (status == nf90_noerr) THEN
    status = nf90_get_var (ncid, varid_s2dbdd , ndbddat , (/mrepsta/), (/nrep/))
  ENDIF
  status = nf90_inq_varid (ncid,'section2_decoding_time' ,varid_s2dbdt )
  IF (status == nf90_noerr) THEN
    status = nf90_get_var (ncid, varid_s2dbdt , ndbdtim , (/mrepsta/), (/nrep/))
  ENDIF

  DO irep = 1 , nrep
    IF ((ndbddat(irep) /= imiss) .AND. (ndbdtim(irep) /= imiss)) THEN
      ndbdtim (irep) = ndbdtim(irep) / 100
!CS
!     IF (     (kcdftyp == ncdf_tempdesc) .OR. (kcdftyp == ncdf_temphirs)) THEN
!       ndbddat (irep) = nr_date(irep)
!       ndbdtim (irep) = nr_time(irep)
!     ENDIF
!CS
      WRITE( yobdat,'(I8,I4.4)' )  ndbddat(irep), ndbdtim(irep)
      READ ( yobdat,'(I4,4I2)' )  ianyy, ianmm, iandd, ianhh, ianmin

      ! time difference in minutes between obs time and initial model time
      CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                  &
                          ianyy, ianmm, iandd, ianhh, ianmin, 0, imindif, ierrf)
!     =================
      ztdecdb (irep)  =  imindif / 60.0_wp
    ELSE
      ztdecdb (irep)  =  rmdi
    ENDIF
  ENDDO

  DEALLOCATE ( ndbddat  , STAT=istat )
  DEALLOCATE ( ndbdtim  , STAT=istat )

!-------------------------------------------------------------------------------
! Section 5: Get station ID
!-------------------------------------------------------------------------------

  lvarid     = .FALSE.
  yfillchar  =  ACHAR( NF_FILL_CHAR )

! get integer station ID's
! ------------------------

! get integer WMO block and station number if present
                            status = nf90_inq_varid (ncid,'MII'  ,varid_mii  )
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NIII' ,varid_niii )
  IF (status == nf90_noerr) THEN
    lvarid  = .TRUE.
!     (these variables are de-allocated in this section again)
    ALLOCATE ( nc_mii   (nrep) , STAT=istat )
    ALLOCATE ( nc_niii  (nrep) , STAT=istat )
    status = nf90_get_var (ncid, varid_mii  , nc_mii , (/mrepsta/), (/nrep/))
    status = nf90_get_var (ncid, varid_niii , nc_niii, (/mrepsta/), (/nrep/))
    DO irpa = 1 , nreproc
      irep = irproc(irpa)
      IF ((nc_mii(irep) /= imiss) .AND. (nc_niii(irep) /= imiss)) THEN
        istidn (irep)  =  nc_mii(irep) *1000  +  nc_niii(irep)
        WRITE( ystidn(irep),'(I5.5)' )  istidn(irep)
      ENDIF
    ENDDO
    DEALLOCATE ( nc_mii , nc_niii   , STAT=istat )
  ENDIF
! else get integer WMO BUOY identifier if present
  IF (      ((status /= nf90_noerr) .OR. (kcdftyp == ncdf_buoy))               &
      .AND. (nf90_inq_varid (ncid,'MABNN',varid_mabnn ) == nf90_noerr)) THEN
    lvarid  = .TRUE.
    status = nf90_get_var (ncid, varid_mabnn, istidn, (/mrepsta/), (/nrep/))
    DO irpa = 1 , nreproc
      irep = irproc(irpa)
      IF ((istidn(irep) /= imiss) .AND. (istidn(irep) < 99999)) THEN
        WRITE( ystidn(irep),'(I5.5)' )  istidn(irep)
      ENDIF
    ENDDO
! else get integer WMO marine obs platform extended identifier if present
  ELSEIF (      (status /= nf90_noerr)                                         &
          .AND. (nf90_inq_varid (ncid,'NMOPI',varid_nmopi ) == nf90_noerr)) THEN
    lvarid  = .TRUE.
    status = nf90_get_var (ncid, varid_nmopi, istidn, (/mrepsta/), (/nrep/))
    DO irpa = 1 , nreproc
      irep = irproc(irpa)
      IF ((istidn(irep) /= imiss) .AND. (istidn(irep) < 9999999))              &
        WRITE( ystidn(irep),'(I7.7)' )  istidn(irep)
    ENDDO
! else get integer WMO satellite identifier if present
  ELSEIF (      (status /= nf90_noerr)                                         &
          .AND. (nf90_inq_varid (ncid,'MI1I2',varid_mi1i2 ) == nf90_noerr)) THEN
    lvarid  = .TRUE.
    status = nf90_get_var (ncid, varid_mi1i2, istidn, (/mrepsta/), (/nrep/))
    DO irpa = 1 , nreproc
      irep = irproc(irpa)
      IF ((istidn(irep) /= imiss) .AND. (istidn(irep) < 999))                  &
        WRITE( ystidn(irep),'(I3.3)' )  istidn(irep)
    ENDDO
  ENDIF

! get character station ID's
! --------------------------

! IF (lax) THEN

  !   get NetCDF variable ID for station ID's
  !   and NetCDF dimension ID for string length of station ID's
  DO ic = 1 , 8
    IF (ic == 1)  ymnem_id = 'YWIGLI'       ! WIGOS (e.g. descend. TEMP)
    IF (ic == 2)  ymnem_id = 'YDDDD '
    IF (ic == 3)  ymnem_id = 'YSSOSN'
    IF (ic == 4)  ymnem_id = 'YAIRN '       ! for AMDAR, ACARS
    IF (ic == 5)  ymnem_id = 'YCCC8 '       ! for METAR
    IF (ic == 6)  ymnem_id = 'YXXNN '       ! for TEMPDROP
    IF (ic == 7)  ymnem_id = 'YSOSN '
    IF (ic == 8)  ymnem_id = 'YCCCC '
    ymnem_strl  =  '             '
    ymnem_strl  =  TRIM( ymnem_id ) // '_strlen'
    dimid_ystid = 0
    status = imdi
    IF (nf90_inq_varid (ncid, ymnem_id ,varid_ystid ) == nf90_noerr)           &
      status = nf90_inq_dimid (ncid, ymnem_strl ,dimid_ystid)
    IF ((status == nf90_noerr) .AND. (dimid_ystid >= 1)) THEN
      lvarid  = .TRUE.
      !   get NetCDF string length of station ID's
      status = nf90_Inquire_Dimension (ncid, dimid_ystid, len=ilen)
      !     (this variable is de-allocated in this section again)
      ALLOCATE ( ync_stid (ilen,nrep) , STAT=istat )
      !   get character station ID's
      ync_stid = ' '
      status = nf90_get_var (ncid, varid_ystid, ync_stid , start=(/1,mrepsta/) &
                                                         , count=(/ilen,nrep/))
      DO irpa = 1 , nreproc
        irep = irproc(irpa)
        IF (      (ystidn(irep) (1:1) == ' ')                                  &
            .AND. (ync_stid(1,irep) (1:1) /= yfillchar)) THEN
          DO ily = 1 , MIN( ilen, ilstidn )
            ystidn (irep) (ily:ily)  =  ync_stid(ily,irep)
            !   check for unrepresentable characters
            IF (ICHAR( ync_stid(ily,irep) ) <= 31) ystidn (irep) (ily:ily) = ' '
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE ( ync_stid , STAT=istat )
    ENDIF
  ENDDO
! ENDIF

! write alerting CAUTION message to file YUCAUTN if station ID variable missing
! -----------------------------------------------------------------------------

  IF (.NOT. lvarid) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
    WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,": STATION IDENTITY VARIABLE "    &
                   &,"MISSING IN FILE:",A)' )  acthr, yfn
    WRITE( nucautn,'("        ==>  CHECK   NetCDF - FILE ",A)' ) yfn
    PRINT          '("CAUTION !!!!! t=",F6.3,": STATION IDENTITY VARIABLE "    &
                   &,"MISSING IN FILE:",A)' ,  acthr, yfn
    PRINT          '("        ==>  CHECK   NetCDF - FILE ",A)' , yfn
    CLOSE (nucautn)
! ELSE
!TODO     WRITE : NO STATION ID'S AVAILABLE !
  ENDIF
! for bad station ID's (of AMDARs), see below

! for GPS reports: assign CMA code type according to sub-centre and station ID
! ----------------------------------------------------------------------------

  IF (kcdftyp == ncdf_gps_zenith) THEN
    DO ic = 1 , n_gpc
      DO irpa = 1 , nreproc
        irep = irproc(irpa)
!       IF (kcensub(irep)     == gpc(ic) %gpscen) icdtyp (irep) = gpc(ic) %cdtyp
        IF (ystidn(irep)(6:9) == gpc(ic) %name)   icdtyp (irep) = gpc(ic) %cdtyp
      ENDDO
    ENDDO
  ENDIF
  DEALLOCATE ( kcensub , STAT=istat )

! for AMV: assign character station ID according to integer satellite ID
! ----------------------------------------------------------------------
  IF (kcdftyp == ncdf_satob) THEN
    DO irpa = 1 , nreproc
      irep = irproc(irpa)
!CS:  the following line would have been the more elegant and general solution
!     but would have required the use of 'mo_satid.f90' in 'libradiance'
!     which is available only for users that need to apply RTTOV
!     ystidn (irep) (1:8) = satname( istidn(irep) )
      IF (istidn(irep) ==  54)  ystidn (irep) = 'METEOS_7  '
      IF (istidn(irep) ==  55)  ystidn (irep) = 'METEOS_8  '
      IF (istidn(irep) ==  56)  ystidn (irep) = 'METEOS_9  '
      IF (istidn(irep) ==  57)  ystidn (irep) = 'METEOS10  '
      IF (istidn(irep) ==  70)  ystidn (irep) = 'METEOS11  '
      IF (istidn(irep) ==  71)  ystidn (irep) = 'MSG-1     '
      IF (istidn(irep) ==  72)  ystidn (irep) = 'MSG-2     '
      IF (istidn(irep) ==  73)  ystidn (irep) = 'MSG-3     '
      IF (istidn(irep) == 257)  ystidn (irep) = 'GOES_13   '
      IF (istidn(irep) == 258)  ystidn (irep) = 'GOES_14   '
      IF (istidn(irep) == 259)  ystidn (irep) = 'GOES_15   '
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 6: Get latitude / longitude / station height
!-------------------------------------------------------------------------------

                            status = nf90_inq_varid (ncid,'MLAH'  ,varid_lat  )
  IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MLALA' ,varid_lat  )
  IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MLOH'  ,varid_lon  )
  IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MLOLO' ,varid_lon  )
  IF (status /= nf90_noerr) THEN
    yerrmsl = 'STANDARD LAT / LON VARIABLES DO NOT EXIST IN ' // yfn
    CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
  ENDIF

!   (these variables are de-allocated in this section again)
  ALLOCATE ( rc_altp  (nrep) , STAT=istat )
  ALLOCATE ( rc_alt   (nrep) , STAT=istat )
  ALLOCATE ( nc_alt   (nrep) , STAT=istat )
!   (these variables are de-allocated in section 7 again)
  ALLOCATE ( nc_altq  (nrep) , STAT=istat )
  ALLOCATE ( nc_locq  (nrep) , STAT=istat )
  rc_lat  = rmiss
  rc_lon  = rmiss
  rc_altp = rmiss
  rc_alt  = rmiss
  nc_alt  = imiss
  nc_altq = imiss
  nc_locq = imiss

! read lat / lon
! --------------

  status = nf90_get_var (ncid, varid_lat , rc_lat , (/mrepsta/), (/nrep/))
  status = nf90_get_var (ncid, varid_lon , rc_lon , (/mrepsta/), (/nrep/))

! quality of location
! -------------------

                            status = nf90_inq_varid (ncid,'MQOBL ',varid_locq )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_locq, nc_locq, (/mrepsta/), (/nrep/))

! read barometer or station altitude
! ----------------------------------

                            status = nf90_inq_varid (ncid,'MHOBNN',varid_altp )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_altp, rc_altp, (/mrepsta/), (/nrep/))
                            status = nf90_inq_varid (ncid,'MHOSNN',varid_alt  )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_alt , rc_alt , (/mrepsta/), (/nrep/))
                            status = nf90_inq_varid (ncid,'MHP   ',varid_alti )
  IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MH    ',varid_alti )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_alti, nc_alt , (/mrepsta/), (/nrep/))

                            status = nf90_inq_varid (ncid,'MSEQM ',varid_altq )
  IF (status == nf90_noerr)                                                    &
    status = nf90_get_var (ncid, varid_altq, nc_altq, (/mrepsta/), (/nrep/))

!   (these variables are de-allocated at the end of this subroutine again)
  ALLOCATE ( mpassiv  (nrep) , STAT=istat )
  ALLOCATE ( icmaa    (nrep) , STAT=istat )

! ----------------------------------
! loops over reports to be processed (in wanted time interval + valid obs type)
! ----------------------------------

  DO irpa = 1 , nreproc
    irep = irproc(irpa)
    kobtyp = iobtyp(irep)
    kcdtyp = icdtyp(irep)

!   get diagnostic array position (--> i_cma)
!   -----------------------------------------
    icmaa (irep)  =  i_cma ( kobtyp , kcdtyp )
!                    =====
    icma   = icmaa(irep)

!   pre-set passive flag to zero
    mpassiv (irep) = 0
    mrepflg (irep) = 0

!   lupair = (kobtyp == ntemp) .OR. (kobtyp == npilot) .OR. (kobtyp == nairep) &
!                              .OR. (kobtyp == nsatem) .OR. (kobtyp == nsatob)
    lupair =       (kobtyp /= nsynop) .AND. (kobtyp /= ngps)                   &
             .AND. (kobtyp /= ndribu) .AND. (kobtyp /= nscatt)

!   (for AOF data base flag on lat / lon / date / time, MLAHQ, MLOHQ, MGGQ and
!    NGGQ have been used in makeaof; currently, these Q-bits are not used here)
!   use MSEQM for altitude (WMO code table 033024):
!             1,5: excellent (within 3m); 2,6: good (<10m); 3,7: fair ( <20m);
!             4,8: poor ( >20m); 15: missing

!   station altitude  (station altitude fits better to reported surface
!   ----------------   pressure obs in TEMPs than barometer altitude (!)
!                      but for surface observations (Synops), barometer
!                      altitude has to be taken for reported surface pressure)
    zaltob (irep) = rmdi
!   mbaro  = 0
    IF (      (ABS( rc_alt (irep) ) < rmisschk)                                &
        .AND. (     (kcdftyp == ncdf_temp) .OR. (kcdftyp == ncdf_tempship)     &
               .OR. (kcdftyp == ncdf_tower))) THEN
      zaltob (irep) = rc_alt (irep)
    ELSEIF (ABS( rc_altp(irep) ) < rmisschk) THEN
!     mbaro  = 1
      zaltob (irep) = rc_altp(irep)
    ELSEIF (ABS( rc_alt (irep) ) < rmisschk) THEN
      zaltob (irep) = rc_alt (irep)
    ELSEIF (     nc_alt(irep)   /= imiss) THEN
      zaltob (irep) = REAL( nc_alt (irep) , wp )
    ENDIF
    IF (kcdftyp == ncdf_tempdesc)  zaltob (irep) = rmdi

!   set missing station altitude of a sea surface obs (ship, buoy, scat) to zero
!   ----------------------------------------------------------------------------
    IF ((zaltob(irep) <= rmdich) .AND. (     (kcat(irep) == 1)                 &
                                        .OR. (kcat(irep) == 12))) THEN
      zaltob (irep) = c0
    ENDIF

!   reject (surface reports) or flag report if station altitude is missing
!   ----------------------------------------------------------------------
    IF (zaltob(irep) <= rmdich) THEN
      IF (     ((.NOT. lupair) .OR. (kobtyp == ntemp) .OR. (kobtyp == npilot)) &
         .AND. (mpassiv(irep) == 0))                                           &
        neventr (nenoal,icma) = neventr(nenoal,icma) + 1
      IF (.NOT. lupair) THEN
        mpassiv (irep) = 2
        kproc (irep)  =  kproc(irep) + 32
      ENDIF
    ENDIF
    nelevt = imdi
    IF (zaltob(irep) > rmdich)  nelevt = NINT( zaltob(irep) )

!   observation location  (allow longitude values between -180 and +360)
!   --------------------
    IF (      (ABS( rc_lat(irep) ) < rmisschk)                                 &
        .AND. (ABS( rc_lon(irep) ) < rmisschk)                                 &
        .AND. (ABS( rc_lat(irep) ) <=   90._wp)                                &
        .AND. (     rc_lon(irep)   <=  360._wp +epsy)                          &
        .AND. (     rc_lon(irep)   >= -180._wp -epsy)) THEN

      IF ((kcdftyp == ncdf_tempdesc) .AND. (isurfob(irep) < 0)) THEN
        rc_lat (irep)  =  rc_lat(irep) + rc_dlatr(irep)
        rc_lon (irep)  =  rc_lon(irep) + rc_dlonr(irep)
        IF (rc_lon(irep) > 360._wp +epsy)  rc_lon(irep) = rc_lon(irep) - 360._wp
        IF (rc_lon(irep) <-180._wp -epsy)  rc_lon(irep) = rc_lon(irep) + 360._wp
      ELSEIF ((kcdftyp == ncdf_tempdesc) .AND. (isurfob(irep) == 0)) THEN
        IF (mpassiv(irep) == 0)                                                &
          neventr (nenoda,icma) = neventr (nenoda,icma) + 1
        mpassiv (irep) = 2
        kproc (irep)  =  kproc(irep) + 64
      ENDIF

!CSC: to remove  (manipulates obs from COSMO-EU domain into COSMO-DE domain)
!     IF (kcdftyp == ncdf_buoy) THEN
!       rc_lon(irep)  =  rmod( rc_lon(irep)      , 6._wp ) +  2._wp
!       rc_lat(irep)  =  rmod( rc_lat(irep)      , 2._wp ) + 54._wp
!!      rc_lon(irep)  =  - rc_lon(irep)
!!      rc_lat(irep)  =    rc_lat(irep) - c1
!!    ELSE
!!      rc_lon(irep)  =  rmod( rc_lon(irep)      , 8._wp ) + 10._wp
!!      rc_lat(irep)  =  rmod( rc_lat(irep)      , 5._wp ) + 50._wp
!!    ELSE
!!      rc_lon(irep)  =  rmod( rc_lon(irep)+78.0 , 16._wp ) +  2._wp
!!      rc_lat(irep)  =  rmod( rc_lat(irep)+5.0  , 10._wp ) + 45._wp
!     ELSEIF (kcdftyp == ncdf_pilot_p) THEN
!       rc_lon(irep)  =  rmod( rc_lon(irep)      , 1._wp ) +  9._wp
!       rc_lat(irep)  =  rmod( rc_lat(irep)      , 5._wp ) + 50._wp
!     ELSEIF (kcdftyp == ncdf_satob) THEN
!       rc_lon(irep)  =  rmod( ABS(rc_lon(irep)) , 30._wp ) +  0._wp
!       rc_lat(irep)  =  rmod( rc_lat(irep)      , 30._wp ) + 30._wp
!     ELSEIF (kcdftyp == ncdf_tower) THEN
!       zaltob(irep)  =  zaltob(irep) + 40._wp
!     ENDIF
!CSC: to remove

      roblat     =  rc_lat(irep)
      roblon     =  rc_lon(irep)

      lindenb = .FALSE.
      IF ((ystidn(irep)(1:5) == '10393').OR.(ystidn(irep)(1:5) == '10394')) THEN
!       Lindenberg Met Observatory reports must be assigned to a specific grid
!       point, that is used by Lindenberg staff for diagnostic purposes and has
!       appropriate hand-made values in the external soil / surface parameters.
!       This is done by assigning the report to pre-specified geographical
!       coordinates and then assigning this to the horizontally nearest model
!       grid point (by setting nelevt=imdi for obs_assign_gridpt)
        IF (      (roblat > 52.1_wp) .AND. (roblat < 52.35_wp)                &
            .AND. (roblon > 14.1_wp) .AND. (roblon < 14.25_wp)) THEN
          lindenb = .TRUE.
          roblat  =  52.220_wp
          roblon  =  14.135_wp
          nelevt  =  imdi
        ENDIF
      ENDIF

!     Assign observation to grid point
!     --------------------------------
      lseaobs   =      (kobtyp == ndribu) .OR. (kobtyp == nscatt)              &
                  .OR. (kcdtyp == nshscd) .OR. (kcdtyp == natshs)              &
                  .OR. (kcdtyp == nshtcd) .OR. (kcdtyp == nshpcd)              &
                  .OR. (kcdtyp == nbship)
      lseaob2   = lseaobs
      lvertfit  = .FALSE.
      zfdoro(:) = fdoro(:)
      ! tower obs: find grid point with smallest absolute height diff.
      IF (kcdtyp == ntower)  zfdoro (5) = 1
      IF (kcdtyp == ntower)  lvertfit   = .TRUE.

      CALL obs_assign_gridpt ( roblat, roblon, nelevt, lupair, .FALSE.         &
                             , lvertfit, ie_tot, je_tot, hsurf_tot, fland_tot  &
                             , doromx, zfdoro, imdi, nboundlines , nolbc       &
                             , lseaobs , iob_tot, job_tot, zio_tot, zjo_tot    &
                             , zsurf, lsfcob )
!     ======================

!                       isurfob (irep) = 0
      IF (lsfcob)       isurfob (irep) = 1
      IF (lseaobs)      mrepflg (irep) = IBSET( mrepflg(irep) , nvsebp )
      IF (iob_tot /= 0) zaltmo  (irep) = zsurf

!     special Lindenberg treatment:
!     report is never rejected here even if "model orography - station altitude"
!     is too large (which could occur if only if doromx(.) is chosen very small)
      IF (lindenb) THEN
        iob_tot = ABS( iob_tot )
        job_tot = ABS( job_tot )
        IF (kcdftyp /= ncdf_tempdesc)  isurfob (irep) = 1
      ENDIF

!     for bogus observations in semi-idealised experiments only
!     --> select it through ob?lat, ob?lon instead of ionl(2)
!     IF (     (ABS(fbogus) > epsy)                                            &
!         .AND.((iob_tot /= ionl ) .OR. (job_tot /= jonl ))                    &
!         .AND.((iob_tot /= ionl2) .OR. (job_tot /= jonl2)))  iob_tot = 0

!     reject report if observation location outside model domain
!     ----------------------------------------------------------
      IF   (iob_tot == 0) THEN
        IF (mpassiv(irep) == 0)                                                &
          neventr (neloca,icma) = neventr (neloca,icma) + 1
        mpassiv (irep) = 2
        kproc (irep)  =  kproc(irep) + 16
      ENDIF

!     set report to passive if (buoy) observation location doubtful
!     (flagged, or time of last known position outside [-6h,1h] of the obs time)
!     --------------------------------------------------------------------------
      IF ((nc_locq(irep) == 2) .OR. (min_dtl(irep) > 360)                      &
                               .OR. (min_dtl(irep) < -60)) THEN
        IF (mpassiv(irep) == 0)                                                &
          neventr (nedbfl,icma) = neventr (nedbfl,icma) + 1
        mpassiv (irep) = 2
        mrepflg (irep) = IBSET( mrepflg(irep) , nvssbp )
        kflag   (irep) = IBSET( kflag  (irep) , FL_SUSP_LOCT )
      ENDIF

!     set report to passive if sea obs is over model land
!     ---------------------------------------------------
      IF ((lseaob2) .AND. (.NOT. lseaobs)) THEN
        IF (mpassiv(irep) == 0)                                                &
          neventr (nezdif,icma) = neventr (nezdif,icma) + 1
        mpassiv (irep) = 2
!       mrepflg (irep) = IBSET( mrepflg(irep) , nvsebp )
        kflag   (irep) = IBSET( kflag  (irep) , FL_SURF )
      ENDIF

      rio_tot (irep) = zio_tot
      rjo_tot (irep) = zjo_tot
      iobstot (irep) = iob_tot
      jobstot (irep) = job_tot

! write alerting CAUTION messages onto file YUCAUTN
! if obs location undefined or invalid
! -------------------------------------------------

    ELSE
      IF (mpassiv(irep) == 0)                                                  &
        neventr (neloca,icma) = neventr (neloca,icma) + 1
      mpassiv (irep) = 2
      kproc (irep)  =  kproc(irep) + 8
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=istat)
      IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,", FILE:",A,": OBS LOCATION UN" &
                     &,"DEFINED, STA:",A,F6.1," HRS")' )                       &
             acthr, yfn, ystidn(irep), zr_hour(irep)
      PRINT          '("CAUTION !!!!! t=",F6.3,", FILE:",A,": OBS LOCATION UN" &
                     &,"DEFINED, STA:",A,F6.1," HRS")' ,                       &
             acthr, yfn, ystidn(irep), zr_hour(irep)
      CLOSE (nucautn)
    ENDIF
  ENDDO

  DEALLOCATE ( rc_altp  , STAT=istat )
  DEALLOCATE ( rc_alt   , STAT=istat )
  DEALLOCATE ( nc_alt   , STAT=istat )
  DEALLOCATE ( min_dtl  , STAT=istat )
  IF (kcdftyp == ncdf_tempdesc) THEN
    DEALLOCATE ( rc_dlatr , STAT=istat )
    DEALLOCATE ( rc_dlonr , STAT=istat )
  ENDIF
  ALLOCATE ( lwrsfc  (nrep) , STAT=istat )
  lwrsfc (:) = .FALSE.

! re-determine the vector of indices of reports which have to be processed now
! (obs time ok, known data category, obs location defined and in model domain)
! ----------------------------------------------------------------------------

  nreproc = 0
  DO irep = 1 , nrep
    IF (kproc(irep) == 0) THEN
      nreproc = nreproc + 1
      irproc (nreproc)  =  irep
    ENDIF
  ENDDO
! PRINT *, 'comnrep3 ', nrep, nreproc

! second, vectorisable loop over reports to be processed

  DO irpa = 1 , nreproc
    irep = irproc(irpa)

    kobtyp = iobtyp(irep)
    kcdtyp = icdtyp(irep)
    roblat = rc_lat(irep)
    roblon = rc_lon(irep)
    icma   = icmaa (irep)

!   observation location outside user-defined area
!   ----------------------------------------------
    IF (     (roblat > cma(icma)%obnlat) .OR. (roblat < cma(icma)%obslat)      &
        .OR. (      (cma(icma)%obwlon >= cma(icma)%obelon)                     &
              .AND. (roblon < cma(icma)%obwlon)                                &
              .AND. (roblon > cma(icma)%obelon))                               &
        .OR. (      (cma(icma)%obwlon <  cma(icma)%obelon)                     &
              .AND. (     (roblon < cma(icma)%obwlon)                          &
                     .OR. (roblon > cma(icma)%obelon)))) THEN
      IF (mpassiv(irep) == 0)                                                  &
        neventr (neloca,icma) = neventr (neloca,icma) + 1
      mpassiv (irep) = 2
      mrepflg (irep) = IBSET( mrepflg(irep) , nvobbp )
      kflag   (irep) = IBSET( kflag  (irep) , FL_AREA )
    ENDIF

!   distance "model orography - station altitude" too large
!   -------------------------------------------------------
!   SYNOP / DRIBU
    IF (iobstot(irep) <  0) THEN
      IF (mpassiv(irep) == 0)                                                  &
        neventr (nezdif,icma) = neventr (nezdif,icma) + 1
      mpassiv (irep) = 2
      mrepflg (irep) = IBSET( mrepflg(irep) , nvalbp )
      kflag   (irep) = IBSET( kflag  (irep) , FL_HEIGHT )
      iobstot (irep) = ABS( iobstot(irep) )
      jobstot (irep) = ABS( jobstot(irep) )
    ENDIF
!   TEMP / PILOT
    IF (      ((kobtyp == ntemp) .OR. (kobtyp == npilot))                      &
        .AND. (isurfob(irep) == 0) .AND. (kcdtyp /= ntdesc)) THEN
      IF ((zaltob(irep) > rmdich) .AND. (mpassiv(irep) == 0)) THEN
        neventr (nezdif,icma) = neventr (nezdif,icma) + 1
        mrepflg (irep) = IBSET( mrepflg(irep) , nvalbp )
      ENDIF
      lwrsfc  (irep) = .TRUE.
!     kflag   (irep) = IBSET( kflag  (irep) , FL_HEIGHT )
    ELSEIF ((nc_altq(irep) == 4) .OR. (nc_altq(irep) == 8)) THEN
!     large inaccuracy ( > 20 m) of station altitude: no surface report
      IF (mpassiv(irep) == 0) THEN
        neventr (nedbfl,icma) = neventr (nedbfl,icma) + 1
        mrepflg (irep) = IBSET( mrepflg(irep) , nvalbp )
        IF ((kobtyp == nsynop) .OR. (kobtyp == ngps)) THEN
          mpassiv (irep) = 2
          kflag   (irep) = IBSET( kflag(irep) , FL_SUSP_LOCT )
        ELSE
          lwrsfc  (irep) = .TRUE.
        ENDIF
      ENDIF
      IF ((kobtyp == ntemp) .OR. (kobtyp == npilot))  lwrsfc (irep) = .TRUE.
!   ELSEIF (nc_altq(irep) == 2,3,6,7) ... inaccuracy of station altitude > 3m,
!     then no surface pressure obs should be derived, but it is not done anyway
    ENDIF

!-------------------------------------------------------------------------------
! Section 7: Check for station ID, observation type, code type
!-------------------------------------------------------------------------------

!   Check for bad aircraft station id
!   ---------------------------------
    IF (      (kobtyp == nairep)                                               &
        .AND. (     (ystidn(irep)(1:3) == "XXX")                               &
               .OR. (ystidn(irep)(1:3) == "???")                               &
               .OR. (ystidn(irep)(1:3) == "///")                               &
               .OR. (ystidn(irep)(1:3) == "***")                               &
               .OR. (ystidn(irep)(1:3) == "   "))) THEN
      IF (mpassiv(irep) == 0)                                                  &
        neventr (neblak,icma) = neventr(neblak,icma) + 1
      mpassiv (irep) = 2
      mrepflg (irep) = IBSET( mrepflg(irep) , nvbkbp )
      kflag   (irep) = IBSET( kflag  (irep) , FL_BLACKLIST )
    ENDIF

!   Check whether in exclusion area for current observation or code type
!   --------------------------------------------------------------------
    IF (      (roblat <= cma(icma)%exnlat+epsy)                                &
        .AND. (roblat >= cma(icma)%exslat-epsy)                                &
        .AND. (roblon >= cma(icma)%exwlon-epsy)                                &
        .AND. (roblon <= cma(icma)%exelon+epsy))  THEN
      IF (mpassiv(irep) == 0)                                                  &
        neventr (neobct,icma) = neventr(neobct,icma) + 1
      mpassiv (irep) = 2
      mrepflg (irep) = IBSET( mrepflg(irep) , nvexbp )
      kflag   (irep) = IBSET( kflag  (irep) , FL_OBSTYPE )
    ENDIF

!CSC: to be removed when multi-level AMDAR are split up into single-level
!     reports or when the redundancy checking is extended
    IF ((kcdftyp == ncdf_amdar_ml) .OR. (kcdftyp == ncdf_amdar_vp)) THEN
      IF (mpassiv(irep) == 0)                                                  &
        neventr (neobct,icma) = neventr(neobct,icma) + 1
      mpassiv (irep) = 2
      mrepflg (irep) = IBSET( mrepflg(irep) , nvexbp )
      kflag   (irep) = IBSET( kflag  (irep) , FL_OBSTYPE )
    ENDIF

!   passive report     indicator put into report flag word
    IF (mpassiv(irep) == 2)  mrepflg(irep) = IBSET( mrepflg(irep) , nvpsbp )
!   station correction indicator put into report flag word
    IF ((kupdate(irep) >= 1) .AND. (kupdate(irep) /= imiss))                   &
                             mrepflg(irep) = IBSET( mrepflg(irep) , nvscbp )
!   PRINT *, 'comflag ', irep, nreproc, ystidn(irep), kproc(irep), kflag(irep) &
!                      , iobstot(irep), jobstot(irep)
  ENDDO

  DEALLOCATE ( nc_altq  , STAT=istat )
  DEALLOCATE ( nc_locq  , STAT=istat )

!-------------------------------------------------------------------------------
! Section 8: Write messages or rejecting of flagging reports,
!            update counters for processed, active, passive, discarded reports,
!            and create new list of indices of those reports
!            which have to be processed further now
!-------------------------------------------------------------------------------

! (re-)determine the vector of indices of reports for which statistics has to
! be updated and possibly error messages written (i.e. report with obs time ok)

  nreproc = 0
  DO irep = 1 , nrep
    IF (MOD(kproc(irep), 2) < 1) THEN
      nreproc = nreproc + 1
      irproc (nreproc)  =  irep
    ENDIF
  ENDDO

! open file unit 'nurej' if requied (if messages are written and file not open)
  lwrej = .FALSE.
  DO irpa = 1 , nreproc
    irep = irproc(irpa)
    IF (     (kproc(irep)/2 >= 1) .OR. (lwrsfc(irep))                          &
        .OR. (kflag(irep) /= ISHFT( ibit1( kflag(irep),FL_OBSTYPE )            &
                                  , FL_OBSTYPE )))                             &
      lwrej = .TRUE.
  ENDDO
  IF ((lwrej) .AND. (.NOT. lopen_rej)) THEN
    OPEN (nurej ,FILE=yurejct ,FORM='FORMATTED',STATUS='UNKNOWN'               &
                              ,POSITION='APPEND',IOSTAT=istat)
    lopen_rej = .TRUE.
    IF (istat /= 0)  yerr = 'OPENING OF FILE yurejct FAILED'
    IF (istat /= 0)  CALL model_abort (my_cart_id, 7004, yerr, yroutine)
  ENDIF

! write messages
  DO irpa = 1 , nreproc
    irep = irproc(irpa)
    IF (kproc(irep) > 0) THEN
      IF     (MOD(kproc(irep),  4) >=  2) THEN
        WRITE( nurej,'(" STATION ",A ," : UKNOWN OBSERVATION / CODE TYPES",    &
                      &2I11,",",F6.1," HRS")' )                                &
               ystidn(irep), iobtyp(irep), icdtyp(irep), zr_hour(irep)
      ELSEIF (MOD(kproc(irep), 64) >= 32) THEN
        WRITE( nurej,'(" STATION ",A ," : SYNOP STA. HEIGHT MISSING",          &
                      &F6.1," HRS")' ) ystidn(irep), zr_hour(irep)
      ELSEIF (MOD(kproc(irep), 32) >= 16) THEN
        WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF DOMAIN ",2F7.1  &
                     &," , ",F6.1," HRS")' )                                   &
               ystidn(irep), rc_lon(irep), rc_lat(irep), zr_hour(irep)
      ELSEIF (MOD(kproc(irep),128) >= 64) THEN
        WRITE( nurej,'(" STATION ",A ," : DESCENDING TEMP WITHOUT OBS",        &
                      &F6.1," HRS")' ) ystidn(irep), zr_hour(irep)
      ENDIF
    ENDIF
    IF (kflag(irep) > 0) THEN
      IF     (ibit1( kflag(irep),FL_SURF      ) ==  1) THEN
        WRITE( nurej,'(" STATION ",A ," : SEA OBS. IS LOCATED OVER (MODEL) LA" &
                     &,"ND AT ",2F7.1," , ",F6.1," HRS")' )                    &
               ystidn(irep), rc_lon(irep), rc_lat(irep), zr_hour(irep)
      ELSEIF (ibit1( kflag(irep),FL_AREA      ) ==  1) THEN
        WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF USER-SPECIFIED" &
                     &," AREA ",2F7.1," , ",F6.1," HRS")' )                    &
               ystidn(irep), rc_lon(irep), rc_lat(irep), zr_hour(irep)
      ELSEIF (ibit1( kflag(irep),FL_HEIGHT    ) ==  1) THEN
        WRITE( nurej,'(" STATION ",A ," : HEIGHT ",F5.0," DIFF. TO MODEL "     &
                     &,F5.0," TOO LARGE, ",F7.1," HRS")')                      &
               ystidn(irep), zaltob(irep), zaltmo(irep), zr_hour(irep)
      ELSEIF (ibit1( kflag(irep),FL_SUSP_LOCT ) ==  1) THEN
        WRITE( nurej,'(" STATION ",A ," : STA HEIGHT FLAGGED, ",F7.1," HRS")') &
               ystidn(irep), zr_hour(irep)
      ELSEIF (ibit1( kflag(irep),FL_BLACKLIST ) ==  1) THEN
        WRITE( nurej,'(" STATION ",A ," : BAD AIRCRAFT ID")')  ystidn(irep)
        WRITE( nurej ,'(A20,": CAUTION: bad aircraft station ID: ",A ," ==> "  &
                      &,"problem in data base ???")' )  yroutine, ystidn(irep)
        PRINT         '(A20,": CAUTION: bad aircraft station ID: ",A ," ==> "  &
                      &,"problem in data base ???")' ,  yroutine, ystidn(irep)
      ENDIF
    ENDIF
    IF (lwrsfc(irep)) WRITE( nurej,'("STATION ",A," : NO SURFACE-LEVEL REPORT "&
                                   &,"DERIVED FROM TEMP / PILOT")') ystidn(irep)
  ENDDO
  DEALLOCATE ( lwrsfc  , STAT=istat )

! update counters for processed / active / passive / discarded reports
! --------------------------------------------------------------------

  DO irpa = 1 , nreproc
    irep = irproc(irpa)
    icma   = icmaa (irep)
    cma(icma)%cnt_pr = cma(icma)%cnt_pr + 1
    IF (kproc(irep) /= 0) THEN
!     means: report is outside model domain, or a synop without station height
!     ==> report is always thrown away (not only set passive)
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ELSEIF ((kflag(irep) /= 0) .AND. (.NOT. lverpas)) THEN
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ELSEIF  (kflag(irep) /= 0)  THEN
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
    ELSE
      cma(icma)%cnt_ac = cma(icma)%cnt_ac + 1
    ENDIF
  ENDDO

! create new list of indices for those reports which need to be processed now
! ---------------------------------------------------------------------------

  irproc  (:) = 0
  nreproc     = 0

  DO irep = 1 , nrep
    IF (       (kproc(irep) == 0)                                              &
        .AND. ((kflag(irep) == 0) .OR. (lverpas))) THEN
      nreproc = nreproc + 1
      irproc (nreproc)  =  irep
    ENDIF
  ENDDO

  PRINT '(" processing",I6," reports from",I5," - ",I4," [min] from file ",A)' &
        ,nreproc, min_sta, min_end, yfn(1:LEN_TRIM(yfn))

  DEALLOCATE ( mpassiv  , STAT=istat )
  DEALLOCATE ( icmaa    , STAT=istat )
  DEALLOCATE ( kproc    , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_comhead
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_comhead


!===============================================================================
!+ Module procedure in "src_obs_cdfin_comhead" for buffering common header info
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_buffer_comhead ( nreproc, irsort, iirlen, irrlen, iyrlen    &
                                  , ibuflen, ibufsrt, rbufsrt, ybufsrt )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_comhead" fills those elements
!   of the report headers which are common to all observation types, into
!   long buffer arrays ('Xbufsrt' which are used later on to distribute the
!   reports to the different nodes resp. sub-domains).
!
! Method:
!   As input, the numbers of elements for each individual report is used to
!   determine the correct offsets in the long buffer arrays.
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
    nreproc             ,& ! number of reports
    ibuflen (3)            ! length of int/real/char input buffer arrays

  INTEGER        , INTENT (IN)  ::  &
    irsort  (nreproc+1) ,& ! report indices sorted according to 'irnode'
    iirlen  (nreproc+1) ,& ! number of integer   elements  \   in each
    irrlen  (nreproc+1) ,& ! number of real      elements   >  individual
    iyrlen  (nreproc+1)    ! number of character elements  /   report

  INTEGER        , INTENT (OUT) ::  &
    ibufsrt (ibuflen(1))   ! integer   input buffer array

  REAL (KIND=wp) , INTENT (OUT) ::  &
    rbufsrt (ibuflen(2))   ! real      input buffer array

  CHARACTER (LEN=*)  , INTENT (OUT) ::  &
    ybufsrt (ibuflen(3))   ! character input buffer array

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    irps   , irep      ,& ! loop indices
    iioffs , iroffs    ,& ! offset for integer / real / character elements
             iyoffs    ,& ! of current report in the long target buffer arrays
    istat                 ! error status variable

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_buffer_comhead
!-------------------------------------------------------------------------------

  ibufsrt =  0
  rbufsrt = c0
  ybufsrt = ' '

  iioffs = 0
  iroffs = 0
  iyoffs = 0

  DO irps = 1 , nreproc
    irep  =  irsort(irps)
    ybufsrt (iyoffs+ 1)  =  ystidn  (irep)
    rbufsrt (iroffs+ 1)  =  zr_hour (irep)
    rbufsrt (iroffs+ 2)  =  rc_lat  (irep)
    rbufsrt (iroffs+ 3)  =  rc_lon  (irep)
    rbufsrt (iroffs+ 4)  =  zaltob  (irep)
    rbufsrt (iroffs+ 5)  =  zaltmo  (irep)
    rbufsrt (iroffs+ 6)  =  rio_tot (irep)
    rbufsrt (iroffs+ 7)  =  rjo_tot (irep)
    rbufsrt (iroffs+ 8)  =  zsynhr  (irep)
    rbufsrt (iroffs+ 9)  =  ztdecdb (irep)
    ibufsrt (iioffs+ 1)  =  iirlen  (irps)
    ibufsrt (iioffs+ 2)  =  irrlen  (irps)
    ibufsrt (iioffs+ 3)  =  iyrlen  (irps)
    ibufsrt (iioffs+ 4)  =  nr_date (irep)
    ibufsrt (iioffs+ 5)  =  nr_time (irep)
    ibufsrt (iioffs+ 6)  =  kflag   (irep)
    ibufsrt (iioffs+ 7)  =  kcat    (irep)
    ibufsrt (iioffs+ 8)  =  kcatsub (irep)
    ibufsrt (iioffs+ 9)  =  iobtyp  (irep)
    ibufsrt (iioffs+10)  =  icdtyp  (irep)
    ibufsrt (iioffs+11)  =  istidn  (irep)
    ibufsrt (iioffs+12)  =  kcentre (irep)
    ibufsrt (iioffs+13)  =  kupdate (irep)
    ibufsrt (iioffs+14)  =  kz_dwd  (irep)
    ibufsrt (iioffs+15)  =  iobstot (irep)
    ibufsrt (iioffs+16)  =  jobstot (irep)
    ibufsrt (iioffs+17)  =  isurfob (irep)
    ibufsrt (iioffs+18)  =  mrepflg (irep)
    ibufsrt (iioffs+19)  =  iobsloc (irep)
    ibufsrt (iioffs+20)  =  jobsloc (irep)
    ibufsrt (iioffs+21)  =  nsynhr  (irep)

    iioffs  =  iioffs + iirlen(irps)
    iroffs  =  iroffs + irrlen(irps)
    iyoffs  =  iyoffs + iyrlen(irps)
  ENDDO

  DEALLOCATE ( zr_hour  , STAT=istat )
  DEALLOCATE ( nr_date  , STAT=istat )
  DEALLOCATE ( nr_time  , STAT=istat )

  DEALLOCATE ( kflag    , STAT=istat )
  DEALLOCATE ( nsynhr   , STAT=istat )
  DEALLOCATE ( kcat     , STAT=istat )
  DEALLOCATE ( kcatsub  , STAT=istat )
  DEALLOCATE ( iobtyp   , STAT=istat )
  DEALLOCATE ( icdtyp   , STAT=istat )
  DEALLOCATE ( istidn   , STAT=istat )
  DEALLOCATE ( ystidn   , STAT=istat )
  DEALLOCATE ( kcentre  , STAT=istat )
  DEALLOCATE ( kupdate  , STAT=istat )
  DEALLOCATE ( kz_dwd   , STAT=istat )
  DEALLOCATE ( iobstot  , STAT=istat )
  DEALLOCATE ( jobstot  , STAT=istat )
  DEALLOCATE ( mrepflg  , STAT=istat )
  DEALLOCATE ( iobsloc  , STAT=istat )
  DEALLOCATE ( jobsloc  , STAT=istat )
  DEALLOCATE ( isurfob  , STAT=istat )

  DEALLOCATE ( zsynhr   , STAT=istat )
  DEALLOCATE ( rc_lat   , STAT=istat )
  DEALLOCATE ( rc_lon   , STAT=istat )
  DEALLOCATE ( zaltob   , STAT=istat )
  DEALLOCATE ( zaltmo   , STAT=istat )
  DEALLOCATE ( rio_tot  , STAT=istat )
  DEALLOCATE ( rjo_tot  , STAT=istat )
  DEALLOCATE ( ztdecdb  , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_buffer_comhead
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_buffer_comhead


!===============================================================================
!+ Module procedure in "src_obs_cdfin_comhead" for storing common header info
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_store_comhead (kodrtyp, noba, nrepl, nlenli, nlenlr, nlenly &
                                 ,ibuf, rbuf, ybuf, nodrnew, nexcess, ksurfob )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_comhead" fills those elements
!   of the report headers which are common to all observation types, into
!   the ODR long-term storage arrays. In addition, some other ODR header
!   elements are given default values.
!
! Method:
!   As input, use buffer arrays 'ibuf', 'rbuf', 'ybuf' which have been filled
!   in routine 'obs_cdf_buffer_comhead'.
!   Note that header elements specific to observation types are filled outside
!   of this routine.
!   Also determine numbers 'nodrnew' of new reports stored in the ODR and
!   'nexcess' of reports to be discarded due to insufficient ODR size.
!
! Initial release: Christoph Schraff, 02.01.08
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
    kodrtyp           ,& ! type of ODR arrays: 1: conventional multi-level
                         !                     2: conventional single-level
                         !                     3: ground-based GPS
                         !                     4: surface-level from TEMP
    noba              ,& ! number of previous reports in ODR
    nrepl             ,& ! number of            reports  \   at
    nlenli            ,& ! number of integer   elements   \  the
    nlenlr            ,& ! number of real      elements   /  local
    nlenly               ! number of character elements  /   node

  INTEGER        , INTENT (IN)  ::  &
    ibuf (nlenli+1)      ! local buffer array containing all integer elements

  REAL (KIND=wp) , INTENT (IN)  ::  &
    rbuf (nlenlr+1)      ! local buffer array containing all real elements

  CHARACTER (LEN=10)  , INTENT (IN)  ::  &
    ybuf (nlenly+1)      ! local buffer array containing all character elements

  INTEGER        , INTENT (OUT) ::  &
    nodrnew           ,& ! number of new reports to be put in ODR array
    nexcess              ! number of reports in excess of array size

  INTEGER        , INTENT (INOUT) , OPTIONAL ::  &
    ksurfob (nrepl+1)    ! surface report indicator and report index
                         !   (used if kodrtyp = 1 or 4)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    irpl              ,& ! loop index used for buffer arrays
    nob               ,& ! report index in ODR arrays
    mpassiv           ,& ! passive report flag
    kobtyp  , kcdtyp  ,& ! observation type , observation code type
    icma                 ! pointer index for 'cma' / diagnostic arrays

  CHARACTER (LEN=80)       :: &
    ymsg                 ! control message
  CHARACTER (LEN=13)       :: &
    yaux                 ! auxilliary text

! Local arrays:
! ------------

  INTEGER        ::  &
    iioffs (nrepl+1)  ,& ! offset of the integer   elements  \  for each
    iroffs (nrepl+1)  ,& ! offset of the real      elements   > local
    iyoffs (nrepl+1)     ! offset of the character elements  /  report
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_store_comhead
!-------------------------------------------------------------------------------

! determine the offsets of the different reports in the long (local) array
! ------------------------------------------------------------------------
  iioffs (1) = 0
  iroffs (1) = 0
  iyoffs (1) = 0
  DO irpl = 2 , nrepl
    iioffs (irpl)  =  iioffs(irpl-1) + ibuf(iioffs(irpl-1)+1)
    iroffs (irpl)  =  iroffs(irpl-1) + ibuf(iioffs(irpl-1)+2)
    iyoffs (irpl)  =  iyoffs(irpl-1) + ibuf(iioffs(irpl-1)+3)
  ENDDO

! conventional multi-level reports
! --------------------------------
  IF (kodrtyp == 1) THEN
    nexcess  =  MAX( 0 , noba + nrepl - maxmll )
    nodrnew  =  nrepl - nexcess
!   PRINT *,'storeco0 ', noba, nrepl, maxmll, nexcess, nodrnew
    DO irpl = 1 , nodrnew
      nob = noba + irpl
      yomlhd (nob)        = ybuf (iyoffs(irpl)+ 1) (1:MIN(ilstidn,ilstid))
      omlhed (nob,nhtime) = rbuf (iroffs(irpl)+ 1)
      omlhed (nob,nhjlat) = rbuf (iroffs(irpl)+ 2)
      omlhed (nob,nhilon) = rbuf (iroffs(irpl)+ 3)
      omlhed (nob,nhalt ) = rbuf (iroffs(irpl)+ 4)
      omlhed (nob,nhsurf) = rbuf (iroffs(irpl)+ 5)
      omlhed (nob,nhzio ) = rbuf (iroffs(irpl)+ 6)
      omlhed (nob,nhzjo ) = rbuf (iroffs(irpl)+ 7)
      omlhed (nob,nhsynt) = rbuf (iroffs(irpl)+ 8)
      omlhed (nob,nhtddb) = rbuf (iroffs(irpl)+ 9)
      omlhed (nob,nhsolz) = rmdi
      omlhed (nob,nhssos) = rmdi
      momlhd (nob,nhdate) = ibuf (iioffs(irpl)+ 4)
      momlhd (nob,nhhrmn) = ibuf (iioffs(irpl)+ 5)
      momlhd (nob,nhflag) = ibuf (iioffs(irpl)+ 6)
      momlhd (nob,nhpass) = MIN( ibuf(iioffs(irpl)+ 6) , 1 ) *2
      momlhd (nob,nhcat ) = ibuf (iioffs(irpl)+ 7)
      momlhd (nob,nhcats) = ibuf (iioffs(irpl)+ 8)
      momlhd (nob,nhobtp) = ibuf (iioffs(irpl)+ 9)
      momlhd (nob,nhcode) = ibuf (iioffs(irpl)+10)
      momlhd (nob,nhstid) = ibuf (iioffs(irpl)+11)
      momlhd (nob,nhcent) = ibuf (iioffs(irpl)+12)
      momlhd (nob,nhcorr) = ibuf (iioffs(irpl)+13)
      momlhd (nob,nhkz  ) = ibuf (iioffs(irpl)+14)
      momlhd (nob,nhitot) = ibuf (iioffs(irpl)+15)
      momlhd (nob,nhjtot) = ibuf (iioffs(irpl)+16)
      ksurfob (irpl)      = ibuf (iioffs(irpl)+17)
      momlhd (nob,nhschr) = ibuf (iioffs(irpl)+18)
      momlhd (nob,nhio  ) = ibuf (iioffs(irpl)+19)
      momlhd (nob,nhjo  ) = ibuf (iioffs(irpl)+20)
      momlhd (nob,nhsyhr) = ibuf (iioffs(irpl)+21)
! default values in some header entries of conventional multi-level reports
      omlhed (nob,nhvcbu) = c1
      omlhed (nob,nhvcbt) = c1
      omlhed (nob,nhvcbq) = c1
      omlhed (nob,nhvctu) = c1
      omlhed (nob,nhvctt) = c1
      omlhed (nob,nhvctq) = c1
      momlhd (nob,nhqcfw) = -1
      momlhd (nob,nhaexi) = 0
      momlhd (nob,nhuexi) = 0
      momlhd (nob,nhtexi) = 0
      momlhd (nob,nhqexi) = 0
!     PRINT *,'storecom ', yomlhd(nob), irpl, nob, omlhed(nob,nhtime)          &
!                        , momlhd(nob,nhhrmn), momlhd(nob,nhjtot)
    ENDDO
! if multi-level report cannot be produced due to insufficient ODR size
! then set ksurfob = 0 in order to inhibit creation of surface report
    DO irpl = nrepl - nexcess + 1 , nrepl
      ksurfob (irpl)      = 0
    ENDDO

! conventional single-level reports
! ---------------------------------
  ELSEIF (kodrtyp == 2) THEN
    nexcess  =  MAX( 0 , noba + nrepl - maxsgl )
    nodrnew  =  nrepl - nexcess
    DO irpl = 1 , nodrnew
      nob = noba + irpl
      yosghd (nob)        = ybuf (iyoffs(irpl)+ 1) (1:MIN(ilstidn,ilstid))
      osghed (nob,nhtime) = rbuf (iroffs(irpl)+ 1)
      osghed (nob,nhjlat) = rbuf (iroffs(irpl)+ 2)
      osghed (nob,nhilon) = rbuf (iroffs(irpl)+ 3)
      osghed (nob,nhalt ) = rbuf (iroffs(irpl)+ 4)
      osghed (nob,nhsurf) = rbuf (iroffs(irpl)+ 5)
      osghed (nob,nhzio ) = rbuf (iroffs(irpl)+ 6)
      osghed (nob,nhzjo ) = rbuf (iroffs(irpl)+ 7)
      osghed (nob,nhsynt) = rbuf (iroffs(irpl)+ 8)
      osghed (nob,nhtddb) = rbuf (iroffs(irpl)+ 9)
      osghed (nob,nhsolz) = rmdi
      osghed (nob,nhssos) = rmdi
      mosghd (nob,nhdate) = ibuf (iioffs(irpl)+ 4)
      mosghd (nob,nhhrmn) = ibuf (iioffs(irpl)+ 5)
      mosghd (nob,nhflag) = ibuf (iioffs(irpl)+ 6)
      mosghd (nob,nhpass) = MIN( ibuf(iioffs(irpl)+ 6) , 1 ) *2
      mosghd (nob,nhcat ) = ibuf (iioffs(irpl)+ 7)
      mosghd (nob,nhcats) = ibuf (iioffs(irpl)+ 8)
      mosghd (nob,nhobtp) = ibuf (iioffs(irpl)+ 9)
      mosghd (nob,nhcode) = ibuf (iioffs(irpl)+10)
      mosghd (nob,nhstid) = ibuf (iioffs(irpl)+11)
      mosghd (nob,nhcent) = ibuf (iioffs(irpl)+12)
      mosghd (nob,nhcorr) = ibuf (iioffs(irpl)+13)
      mosghd (nob,nhkz  ) = ibuf (iioffs(irpl)+14)
      mosghd (nob,nhitot) = ibuf (iioffs(irpl)+15)
      mosghd (nob,nhjtot) = ibuf (iioffs(irpl)+16)
!!    ksurfob (irpl)      = ibuf (iioffs(irpl)+17)
      mosghd (nob,nhschr) = ibuf (iioffs(irpl)+18)
      mosghd (nob,nhio  ) = ibuf (iioffs(irpl)+19)
      mosghd (nob,nhjo  ) = ibuf (iioffs(irpl)+20)
      mosghd (nob,nhsyhr) = ibuf (iioffs(irpl)+21)
! default values in some header entries of conventional single-level reports
      osghed (nob,nhvcbu) = c1
      osghed (nob,nhvcbt) = c1
      osghed (nob,nhvcbq) = c1
      osghed (nob,nhvctu) = c1
      osghed (nob,nhvctt) = c1
      osghed (nob,nhvctq) = c1
      mosghd (nob,nhqcfw) = -1
    ENDDO

! single-level surface reports from TEMP
! --------------------------------------
  ELSEIF (kodrtyp == 4) THEN
    nodrnew = 0
    nexcess = 0
    nob = noba
    DO irpl = 1 , nrepl
      IF ((ksurfob(irpl) >= 1) .AND. (nob >= maxsgl)) THEN
        nexcess  =  nexcess + 1
        ksurfob (irpl)      = 0
      ELSEIF (ksurfob(irpl) >= 1) THEN
        nodrnew  =  nodrnew + 1
        nob      =  noba    + nodrnew
        yosghd (nob)        = ybuf (iyoffs(irpl)+ 1) (1:MIN(ilstidn,ilstid))
        osghed (nob,nhtime) = rbuf (iroffs(irpl)+ 1)
        osghed (nob,nhjlat) = rbuf (iroffs(irpl)+ 2)
        osghed (nob,nhilon) = rbuf (iroffs(irpl)+ 3)
        osghed (nob,nhalt ) = rbuf (iroffs(irpl)+ 4)
        osghed (nob,nhsurf) = rbuf (iroffs(irpl)+ 5)
        osghed (nob,nhzio ) = rbuf (iroffs(irpl)+ 6)
        osghed (nob,nhzjo ) = rbuf (iroffs(irpl)+ 7)
        osghed (nob,nhsynt) = rbuf (iroffs(irpl)+ 8)
        osghed (nob,nhtddb) = rbuf (iroffs(irpl)+ 9)
        osghed (nob,nhsolz) = rmdi
        osghed (nob,nhssos) = rmdi
        mosghd (nob,nhdate) = ibuf (iioffs(irpl)+ 4)
        mosghd (nob,nhhrmn) = ibuf (iioffs(irpl)+ 5)
        mosghd (nob,nhflag) = ibuf (iioffs(irpl)+ 6)
        mosghd (nob,nhpass) = MIN( ibuf(iioffs(irpl)+ 6) , 1 ) *2
        mosghd (nob,nhcat ) = ibuf (iioffs(irpl)+ 7)
        mosghd (nob,nhcats) = ibuf (iioffs(irpl)+ 8)
        mosghd (nob,nhobtp) = ibuf (iioffs(irpl)+ 9)
        mosghd (nob,nhcode) = ibuf (iioffs(irpl)+10)
        mosghd (nob,nhstid) = ibuf (iioffs(irpl)+11)
        mosghd (nob,nhcent) = ibuf (iioffs(irpl)+12)
        mosghd (nob,nhcorr) = ibuf (iioffs(irpl)+13)
        mosghd (nob,nhkz  ) = ibuf (iioffs(irpl)+14)
        mosghd (nob,nhitot) = ibuf (iioffs(irpl)+15)
        mosghd (nob,nhjtot) = ibuf (iioffs(irpl)+16)
!       ksurfob (irpl)      = ibuf (iioffs(irpl)+17)
        ksurfob (irpl)      = nob
        mosghd (nob,nhschr) = ibuf (iioffs(irpl)+18)
        mosghd (nob,nhio  ) = ibuf (iioffs(irpl)+19)
        mosghd (nob,nhjo  ) = ibuf (iioffs(irpl)+20)
        mosghd (nob,nhsyhr) = ibuf (iioffs(irpl)+21)
! default values in some header entries of conventional single-level reports
        osghed (nob,nhvcbu) = c1
        osghed (nob,nhvcbt) = c1
        osghed (nob,nhvcbq) = c1
        osghed (nob,nhvctu) = c1
        osghed (nob,nhvctt) = c1
        osghed (nob,nhvctq) = c1
        mosghd (nob,nhqcfw) = -1
      ENDIF
    ENDDO

! ground-based GPS reports
! ------------------------
  ELSEIF (kodrtyp == 3) THEN
    nexcess  =  MAX( 0 , noba + nrepl - maxgpl )
    nodrnew  =  nrepl - nexcess
    DO irpl = 1 , nodrnew
      nob = noba + irpl
      yogphd (nob)        = ybuf (iyoffs(irpl)+ 1) (1:MIN(ilstidn,ilstid))
!     yogphd (nob)        = ybuf (iyoffs(irpl)+ 1) (1:4)
      ogphed (nob,nhtime) = rbuf (iroffs(irpl)+ 1)
      ogphed (nob,nhjlat) = rbuf (iroffs(irpl)+ 2)
      ogphed (nob,nhilon) = rbuf (iroffs(irpl)+ 3)
      ogphed (nob,nhalt ) = rbuf (iroffs(irpl)+ 4)
      ogphed (nob,nhsurf) = rbuf (iroffs(irpl)+ 5)
      ogphed (nob,nhzio ) = rbuf (iroffs(irpl)+ 6)
      ogphed (nob,nhzjo ) = rbuf (iroffs(irpl)+ 7)
      ogphed (nob,nhsynt) = rbuf (iroffs(irpl)+ 8)
      ogphed (nob,nhtddb) = rbuf (iroffs(irpl)+ 9)
      ogphed (nob,nhsolz) = rmdi
      ogphed (nob,nhssos) = rmdi
      mogphd (nob,nhdate) = ibuf (iioffs(irpl)+ 4)
      mogphd (nob,nhhrmn) = ibuf (iioffs(irpl)+ 5)
      mogphd (nob,nhflag) = ibuf (iioffs(irpl)+ 6)
      mogphd (nob,nhpass) = MIN( ibuf(iioffs(irpl)+ 6) , 1 ) *2
      mogphd (nob,nhcat ) = ibuf (iioffs(irpl)+ 7)
      mogphd (nob,nhcats) = ibuf (iioffs(irpl)+ 8)
      mogphd (nob,nhobtp) = ibuf (iioffs(irpl)+ 9)
      mogphd (nob,nhcode) = ibuf (iioffs(irpl)+10)
      mogphd (nob,nhstid) = ibuf (iioffs(irpl)+11)
      mogphd (nob,nhcent) = ibuf (iioffs(irpl)+12)
      mogphd (nob,nhcorr) = ibuf (iioffs(irpl)+13)
      mogphd (nob,nhkz  ) = ibuf (iioffs(irpl)+14)
      mogphd (nob,nhitot) = ibuf (iioffs(irpl)+15)
      mogphd (nob,nhjtot) = ibuf (iioffs(irpl)+16)
!     ksurfob (irpl)      = ibuf (iioffs(irpl)+17)
      mogphd (nob,nhschr) = ibuf (iioffs(irpl)+18)
      mogphd (nob,nhio  ) = ibuf (iioffs(irpl)+19)
      mogphd (nob,nhjo  ) = ibuf (iioffs(irpl)+20)
      mogphd (nob,nhsyhr) = ibuf (iioffs(irpl)+21)
! default values in some header entries of ground-based GPS reports
      mogphd (nob,nhqcfw) = -1
    ENDDO
  ENDIF

! if ODR size is too small: update statstics counters, write control messages
! ---------------------------------------------------------------------------
  DO irpl = nrepl - nexcess + 1 , nrepl
    mpassiv = MIN( ibuf(iioffs(irpl)+ 6) , 1 ) *2
    kobtyp  =      ibuf(iioffs(irpl)+ 9)
    kcdtyp  =      ibuf(iioffs(irpl)+10)
    icma    = i_cma ( kobtyp , kcdtyp )
!             =====
!   PRINT *,'noctrj_m3 ', icma, irpl, cma(icma)%cnt_ps, cma(icma)%cnt_rj       &
!                       , nexcess, nrepl, kodrtyp,kobtyp,kcdtyp
    IF (kodrtyp <= 3) THEN
! report event 'insufficient ODR size' is updated even for passive reports !
      neventr (nesodr,icma) = neventr(nesodr,icma) + 1
      IF (mpassiv == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF
    IF ((kodrtyp == 1) .AND. (lwonl))                                          &
      WRITE( nupr,'(" CAUTION !!!!! MULTI-LEVEL REPORT ",A ,": EXCEEDS ODR"    &
                  &," SIZE MAXMLL",I6)' ) ybuf(iyoffs(irpl)+1), maxmll
    IF  (kodrtyp == 1)                                                         &
      PRINT       '(" CAUTION !!!!! MULTI-LEVEL REPORT ",A ,": EXCEEDS ODR"    &
                  &," SIZE MAXMLL",I6)' , ybuf(iyoffs(irpl)+1), maxmll
    IF (irpl == nrepl) THEN
      IF (kodrtyp == 1) WRITE( yaux,'("MAXMLL=",I6)' ) maxmll
      IF (kodrtyp == 2) WRITE( yaux,'("MAXSGL=",I6)' ) maxsgl
      IF (kodrtyp == 3) WRITE( yaux,'("MAXGPL=",I6)' ) maxgpl
      IF (kodrtyp == 4) WRITE( yaux,'("MAXSGL=",I6)' ) maxsgl
      WRITE( ymsg,'(" CAUTION !!!!! ODR SIZE ",A13," EXCEEDED BY",I6," FOR"    &
                  &," OBS/CODE TYPE",I3,I4)' ) yaux, nexcess, kobtyp, kcdtyp
      IF (lwonl) WRITE( nupr,'(A )' )  ymsg
                 PRINT       '(A )' ,  ymsg
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_store_comhead
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_store_comhead

!===============================================================================

ELEMENTAL INTEGER FUNCTION ibit1  ( invar, ibp )
  !---------------------------------------------
  INTEGER         , INTENT (IN)  ::  invar, ibp
  !---------------------------------------------------------------------------
  ! returns 1 bit at bit position 'ibp' of integer word 'invar'
  !---------------------------------------------------------------------------
  !
  ibit1 = IAND( ISHFT( invar,-ibp ), 1 )
  !
END FUNCTION ibit1

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_comhead
