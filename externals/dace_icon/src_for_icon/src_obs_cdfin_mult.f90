!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------
#if __NEC_VERSION__ >= 30005    /* options directive supported since 3.0.5 */
!NEC$ options "-O2 -fno-move-loop-invariants"
#endif

MODULE src_obs_cdfin_mult

!-------------------------------------------------------------------------------
! Description:
!   This module performs the observation pre-processing of multi-level reports
!   read from NetCDF observation input files. The current types of multi-level
!   reports comprise of radiosonde (TEMP and PILOT) and ground-based remote-
!   sensing profiler (Wind Profiler, RASS, Radar VAD) reports.
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
!    - obs_cdf_read_temp_pilot: (called by obs_cdf_read_org)
!    - obs_cdf_read_profiler  : (called by obs_cdf_read_org)
!    - obs_cdf_store_multilev :  called by obs_cdf_read_temp_pilot
!                                      and obs_cdf_read_profiler
!    - check_data_in_report   :  called by obs_cdf_store_multilev
!
!   This module also contains elemental functions, formerly statement functions:
!   - insert       : inserts bits set in one integer into another integer word
!                    at given bit position
!   - ibit1        : returns 1 bit at given bit position of given integer word
!
!   It uses from:
!    - src_obs_cdfin_comhead: - obs_cdf_read_comhead
!                             - obs_cdf_buffer_comhead
!                             - obs_cdf_store_comhead
!    - src_obs_cdfin_blk:     - obs_cdf_whitelist_local
!                             - obs_cdf_blacklist_local
!    - src_obs_cdfin_util:    - obs_assign_sort_node
!                             - obs_cdf_distrib_reports
!                             - obs_average_layers
!                             - obs_rm_duplicate_levels
!                             - std_atmosphere
!                             - std_atm_p2z_fast
!                             - obs_td2rh
!                             - obs_qx2rh
!                             - obs_rhw2rh
!                             - obs_find_level
!                             - f_z2p
!    - utilities:           \ - phi2phirot
!                           \ - rla2rlarot
!                             - uv2uvrot_vec
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
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted (e.g.
!    modified routine interfaces and diagnostic arrays, 'obs_pointrs'-->'i_cma',
!    modified ODR (cloud group words, observation status word)).
!  - Call of external routine 'atmhp' (libmisc) replaced by 'std_atmosphere'.
!  - To make observation input more flexible, some variables in the NetCDF
!    files are changed from mandatory to optional.
!  - Observed 'temperature variable' (e.g. virtual temperature) not converted
!    into temperature any more.
!  - Bug correction at selection of surface observations (usage of 'fdoro').
!  - Bug correction allowing for active use of pressure obs from TEMP surface
!    level (in multi-level report only).
!  - Pressure derived from height flagged in main flag word.
!  - Screen-level obs within multi-level reports set to passive (and flagged).
! V4_28        2013/07/12 Christoph Schraff
!  All height obs flagged except lowest active TEMP z-obs (usually surf. level).
!  Surface pressure obs error converted from height to pressure units.
!  Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Oliver Fuhrer
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Christoph Schraff
!  More flexible variable names for replication factor and height in wind
!  profilers.
! V5_4         2016-03-10 Christoph Schraff
!  Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  Variables related to AOF interface removed.
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Transformat. of observed wind into rotated coordinate
!  system and a posteriori shift of location executed only if run in COSMO.
! V5_4h        2017-12-15 Christoph Schraff
!  For COMET: Reading of 'cdfin_pilot' and 'cdfin_pilot_p' is modified such that
!  both types of files are free to contain (at least) either pressure (MPN),
!  geopotential height (NHHH) or geopotential (NHNHN), and either 'MEDRE' or
!  'MDREP' or both, one of which is used as delayed replication factor for the
!  number of vertical levels.
! V5_6b        2019-10-16 Christoph Schraff
!  - Introduction of tower profile reports, high-resolution BUFR radiosonde
!    profile reports and descending TEMP reports.
!  - Superobbing of high-resolution profiles introduced.
!  - New subroutine 'check_data_in_report' introduced which checks the existence
!    of active or passive data in a multi-level report.
!  - Processing of relative humidity and mixing ratio is added to dewpoint
!    temperature as observed humidity variable for multi-level reports.
!  - Discarding any obs levels with height above 'model top' (defined here as
!    the mean of top main and half model level) or with pressure smaller than
!    certain limits (e.g. p <= 30 hPa with model top <= 22 km).
!  - For everything (e.g. blacklisting, gross error checking, superobbing)
!    affecting the state (active, passive...) of a report or observation,
!    'observed pressure' (nbtp) is replaced by 'pressure level independent
!    from model state' (nbtplv) in order to render the number of reports and
!    obs independent from the model state and thus of the member of an ensemble
!    (this is required for 'online LETKF').
!    Set additional report flag 'nvlzbp': pressure derived from height by
!    using model state.
!  - Bug fix: intialising 'rc_tv' with 'rmiss' instead of 'imiss'.
!  - Bug fix: 'zww', 'zsinor', 'zstdff', 'iqci' were allocated multiple times
!             but never de-allocated.
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - Other levels with the same pressure as the surface level are discarded.
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! V5_7a        2020-05-11 Christoph Schraff
!  - Fix of a lethal bug for 'ltempb' if DWD KZ is missing and lowest level
!    in report is not below 300 hPa.
!  - Fix to allow for TEMP input files with only zero level reports and hence
!    missing NFNFN wind variable.
!  - Correction: omit search for surface level for descending and drop sondes.
!  - Introduction of ground-based wind lidar reports from files with wind
!    profiler template.
! @VERSION@    @DATE@     Christoph Schraff
!  - Meaning of entry 'nbtzio' changed from grid pt. to geographic coordinates.
!  - Variables related to structured grid of COSMO such as 'r_dlat' removed
!    except in code lines conditioned to ifdef COSMO.
!  - Additional safety if 'surface level' was not defined for a descending TEMP.
!    Bug fixes for the case that a surface level is missing (ilevsfc < 1).
!  - Some code clean-up.
!  01.11.2021, CS: Bug fixes: If a superobs is missing value, the corresponding
!                  flag (nbterr) for active obs has to be unset explicitly after
!                  calling obs_average_layers.
!                  Also produce a (single) superobs for very short obs profiles.
!  16.05.2022, CS: Std. dev. of sub-grid scale model orography added as
!                  criterion for 10-m wind
!  18.04.2023, CS: 'rhtsat' is namelist parameter instead of fixed value and
!                  moved from data_obs_cdfin to data_obs_lib_cosmo.
!  24.04.2023, CS: Introduction of ICOS tower profile reports.
!  16.05.2023, CS: Derivation of single-level report from (ICOS) tower surface
!                  level (T2M, RH2M) + surface pressure + global radiation
!                  (compiling 1-hrly sums in new routine 'obs_radsum_1h').
!  05.06.2023, CS: Bug fixes to avoid array bound violations.
!  22.02.2024: CS: Set entry 'nbsvip' to model surface pressure in derived
!                  single-level TEMP reports for later use as 'plevel'
!                  (for localisation in LETKF)
!  14.06.2024: CS: re-set 'nbslid' to allow for identifying single-level reports
!                  from towers
!  04.10.2024: CS: for radiosondes, always take pressure as level info
!                  (indlev=2); p-z gross error check added.
!  07.10.2024: CS: obs type, code type, and status of radiosonde-derived
!                  surface report made namelist-dependent (icdt_rss)
!  07.10.2024: CS: Bug fix: always use ABS in the context of rmisschk, assuming
!                  it is not known whether rmiss is huge positive or negative.
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
    luse_mlz   ,& ! if false then use multi-level T, not z, and set z-obs to
                  !          passive except for the lowest z-p-obs

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

    ! other variables related to namelist parameters
    maxmlv         ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    maxmll         ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    icdt_tws       ,& !   0 : mode for obs/code type of tower surface-level rep.
                      !       = 0 : no correction (remains PILOT/ ICOS tower)
                      !       = +/- 1 : automatic Synop
                      !       = +/- 2 : test      Synop
                      !       > 0: obstype status active if Synop to be active
                      !       < 0: active only if both Synop + tower active
    icdt_rss       ,& !   0 : mode for radiosonde-derived surface-level report
                      !       =-1 : no single-level surface report (SR) created
                      !       = 0 : SR: type + status like parent TEMP
                      !       = 1 : SR: obs type SYNOP, code 835, active only if
                      !             parent TEMP code type active and lcd835 true
                      !       = 2 : SR: obs type SYNOP, code 835, status lcd835
    madj_hum          ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable

USE data_obs_lib_cosmo, ONLY :  &

    ! constants for the horizontal rotated grid and related variables
    r_pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    r_pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
!   r_polgam       ,& ! angle between the north poles of the systems
!   r_dlon         ,& ! grid point distance in zonal direction (in degrees)
!   r_dlat         ,& ! grid point distance in meridional direction (in degrees)
!   r_startlat_tot ,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
!   r_startlon_tot ,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    r_degrad       ,& ! factor for transforming degree to rad

    ! variables related to namelist parameters, and other variables
    doromx         ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    altopsu        ,& ! SYNOP obs. above height 'altopsu' are set passive
    rhtsat         ,& ! relative humidity threshold above which obs is set =100%
    av_reso        ,& ! apply superobbing if the averaged resolution of the obs
                      !   profile exceeds 'av_reso' times the model resolution

    ! physical constants
    r_g            ,& ! acceleration due to gravity
    tmelt          ,& ! melting temperature of ice
    r_d            ,& ! gas constant for dry air
    rdv            ,& ! r_d / r_v
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
    r_ps           ,& ! surface pressure
    r_sso_sd       ,& ! std. dev. of sub-grid scale model orography [m]

! 2d: Quantities related to above pointers
!-----------------------------------------

    pref           ,& ! reference interpolation pressure levels for
                      ! superobbing of high-resolution radiosonde reports
    nref              ! actual number of interpolation reference levels

USE data_obs_lib_cosmo, ONLY :  &

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------
    nurej      ,& ! direct reporting of rejected obs. reports

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
!   nairep     ,& ! AIREP reports (all aircraft reports)
    ntemp      ,& ! TEMP  reports
!   npilot     ,& ! PILOT reports
    natscd     ,& !   automatic synop surface report
    nsytmp     ,& !   surface report from TEMP
    nsytow     ,& !   surface report from tower
    nldtcd     ,& !   temp land   report
!   nmotcd     ,& !   temp mobile report
!   ntdrop     ,& !   temp drop   report
    nbtemp     ,& !   temp land   report (high-res. BUFR)
    ntdesc     ,& !   temp descending report
!   nmopcd     ,& !   pilot mobile report
!   nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
!   npr_us     ,& !   US Wind Profiler/RASS report
    ntower     ,& !   tower profile  report
    ntowic     ,& !   icos tower profile  report

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

!   n_cma      ,& ! number of CMA obs and code types
!   t_cmatyp   ,& ! data type for information on CMA observation and code types
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
    ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_temphirs  ,& ! indicator for proc. NetCDF TEMP high-res BUFR input
    ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_tempdesc  ,& ! indicator for proc. NetCDF descending TEMP    input
    ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
!   ncdf_amdar     ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_wprof     ,& ! indicator for proc. NetCDF wind profiler      input
    ncdf_rass      ,& ! indicator for proc. NetCDF RASS profiler      input
    ncdf_radar_vad ,& ! indicator for proc. NetCDF radar wind prof.   input
    ncdf_tower     ,& ! indicator for proc. NetCDF tower profile      input
    ncdf_tower_icos,& ! indicator for proc. NetCDF tower profile      input
    ncdf_wlidar_wp ,& ! indicator for proc. NetCDF ground-b wind lidar input
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
! common header entries in NetCDF file
    nc_tisi        ,& ! MTISI  : time significance (BUFR Table 008021)    (dito)
! derived common header variables stored to ODR
    iobstot        ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot        ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc        ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc        ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc         ,& ! indices of reports to be processed now
    rc_altp        ,& ! MHOBNN      : barometer altitude    (temporary variable)

!         1.2.2 other NetCDF header entries
!               ---------------------------

    nc_rstyp       ,& ! NRARA , WMO Common Table C2 : radiosonde type/system
    nc_rad         ,& ! NSR   , BUFR Table B 002013 : solar + IR radiation corr.
    nc_track       ,& ! NSASA , WMO Common Table C7 : tracking technique, status
    nc_na4         ,& ! NA4   , BUFR Table B 002003 : type of measur. equipment
!   nc_nix         ,& ! NIX   , BUFR Table B 002001 : station type (man,auto,..)

!         1.3   NetCDF body entries
!               -------------------

!         1.3.1  frequent entries
!                ----------------
    nc_nlev        ,& ! MEDRE / MDREP: delayed descriptor replication factor
                      !                (e.g. number of vertical levels)
    nc_dt          ,& ! NLTPD : time [sec] since launch time
    nc_lvtyp       ,& ! MEVSS , BUFR Tab 008042 : extended vertical sounding
                      !                           significance (level identity)
    nc_z           ,& ! NHHHN, NHHH : geopotential height [gpm]    (upper-air)
    nc_pp          ,& ! MPN   : pressure                    (int or real)
    nc_dd          ,& ! NDNDN : wind direction      [degree]
    rc_p           ,& ! MPN   , MPPP : pressure             (int or real)
    rc_zz          ,& ! NHHHN : geopotential height [gpm]   (int or real)
    rc_dlat        ,& ! MLADH : latitude  displacement since launch site
    rc_dlon        ,& ! MLODH : longitude displacement since launch site
    rc_t           ,& ! MTDBT : temperature / dry-bulb temperature
    rc_td          ,& ! MTDNH : dew-point temperature
    rc_ff             ! NFNFN : wind speed

USE data_obs_cdfin, ONLY :  &

!         1.3.3  additional profiler or PILOT elements
!                ----------------------------

    nc_sinor       ,& ! MSINOR, BUFR Tab 002064 : signal to noise ratio
    nc_qci         ,& ! MQINZ , BUFR Tab 033002 : quality information
    nc_qci2        ,& ! NWPQ  , BUFR Tab 025034 : NOAA QC results   (temporary)
    nc_wce         ,& ! MWCE  , BUFR Tab 025021 : wind computation enhancement
    nc_dto         ,& ! MSETP , time period of measurement [sec]
    nc_dto2        ,& ! NGGTP , time period of measurement [min]    (temporary)
    nc_uuu2        ,& ! MUUU  : relative humidity            (tower)
    nc_tisi2       ,& ! MTISI : time signific. (Table 008021)(tower, temporary)
    nc_peri        ,& ! NGGTP : time period [min]            (tower, temporary)
    nc_qct         ,& ! MADDF*: assoc. signif. (Table 031021)(tower, temporary)
    nc_qc          ,& ! *Q    : DWD quality bits             (tower, temporary)
    nc_dd3         ,& ! NDNDN : wind direction [deg]         (tower, temporary)
    rc_fi          ,& ! NHNHN : PILOT   : geopotential                  [m2/s2]
    rc_w           ,& ! MWMPS : vertical velocity (w-component of wind)   [m/s]
    rc_tv          ,& ! MTVIR : virtual temperature                         [K]
    rc_stdff       ,& ! NSTDFF: standard deviation wind speed             [m/s]
    rc_dz          ,& ! MHOSEN: sensor height a. ground [m]  (tower, temporary)
    rc_azi         ,& ! MDA   : azimuth                      (tower)      [deg]
    rc_ff3         ,& ! NFNFN : wind speed                   (tower, temporary)
    rc_azi3        ,& ! MDA   : bearing or azimuth [deg]     (tower, temporary)

!         1.3.4  additional synoptic elements
!                ----------------------------
    nc_clsig       ,& ! MVTSU , BUFR Tab 008002 : vertical significance
    nc_clclm       ,& ! MNH   , BUFR Tab 020011 :(low or mid-level) cloud amount
    nc_ccl         ,& ! MCC   , BUFR Tab 020012 : cloud type (low clouds)
    nc_ccm         ,& ! MCC0  , BUFR Tab 020012 : cloud type (middle clouds)
    nc_cch         ,& ! MCC1  , BUFR Tab 020012 : cloud type (high clouds)
    rc_cbase       ,& ! NH    : height of base of cloud
    rc_radgl       ,& ! MGLSR : global  solar radiation                  [J/m2]

!         1.4.1  Bit positions for level significance (MEVSS, BUFR Tab 008042)
!                ------------------------------------

    ilv_sfc        ,& ! surface                 level bit
    ilv_std        ,& ! standard                level bit
    ilv_tropo      ,& ! tropopause              level bit
    ilv_max        ,& ! maximum wind            level bit
    ilv_sigt       ,& ! significant temperature level bit
    ilv_sigq       ,& ! significant humidity    level bit
    ilv_sigv       ,& ! significant wind        level bit
    ilv_miss       ,& ! missing value indicator       bit

!         1.5   auxilliary buffer arrays and variables
!               --------------------------------------

    imiss          ,& ! missing value for integers in current NetCDF input file
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk       ,& ! value smaller than 'rmiss', to check for missing value

! Section 2 : Blacklist and Whitelist
!-------------------------------------------------------------------------------

    ilstid_blk     ,& ! assume 8-character station-IDs in Black-/Whitelist
!   maxintv        ,& ! max. number of vertical blacklist intervals per 1 stat.
    blk_loc           ! blacklists for local reports

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
!-------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
    neobct     ,& ! observation or code type excluded on area with sta. location
    nenoda     ,& ! no accepted data in report

!         3.1.2   Data event counter array format
!                 -------------------------------
!   mxdeve     ,& ! length of data event counter array
    nelodr     ,& ! level rejected: number of levels exceeding ODR size
    nelmis     ,& ! level rejected: pressure (PILOT: pressure + height) missing
    nelflg     ,& ! level rejected: pressure (PILOT: height) flagged
    nelsfc     ,& ! level rejected: too many surface levels
    nelnop     ,& ! level rejected: PILOT height level not in range of model lev
    nelext     ,& ! level rejected: pressure < 9hPa or level below station alt.
    nelsig     ,& ! level rejected: significant level above a specified limit
!   nelrdn     ,& ! level rejected: redundant level in report  (not active yet)
    nepmis     ,& ! pressure (TEMP: height): missing
    nepflg     ,& ! pressure (TEMP: height): flagged
!   neprac     ,& ! pressure: bad reporting practice
    nepalt     ,& ! pressure: sta height, or height dist. to orography too large
!   nepsdt     ,& ! pressure tendency: flagged, or absolute value too large
    netmis     ,& ! temperature missing
    netflg     ,& ! temperature flagged
    netext     ,& ! temperature too low or too high
    netalt     ,& ! height (diff.) too large for 2m-temp.
!   netlps     ,& ! lapse rate of multi-level temperature too large
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
    nedmis     ,& ! wind direction missing
    nefmis     ,& ! wind speed missing
    nedflg     ,& ! wind direction flagged
    nefflg     ,& ! wind speed flagged
    nefneg     ,& ! wind speed too small  ( < 0 ; DRIBU <= 0 ; VAD < 3m/s )
    nevalt        ! height (diff.) too large for 10m-wind
!   nefshr     ,& ! wind speed shear too large
!   nedshr     ,& ! directional wind shear too large
!   nerlim        ! prec.amount exceeds threshold limit

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
    oevsond    ,& ! (root of) radiosonde (TEMP, PILOT) wind error variance
    oezsond    ,& ! (root of) radiosonde (TEMP, PILOT) height error variance
    oetsond    ,& ! (root of) radiosonde temperature error variance
    oeairep    ,& ! (root of) AIREP wind error variance
    oetairp    ,& ! (root of) AIREP temperature error variance
!   oevsynp    ,& ! (root of) SYNOP wind error variance
    oezsynp    ,& ! (root of) SYNOP height error variance (land)
!   oesatob    ,& ! (root of) SATOB wind error variance
!   oevdrib    ,& ! (root of) DRIBU wind   height error variance
!   oezdrib    ,& ! (root of) DRIBU height height error variance
!   oezship    ,& ! (root of) SHIP (sea SYNOP) height error variance
    rherr1     ,& ! (root of) fixed    / normal conditions
    rherr2     ,& ! relative humidity <  if temperature below 233K
    rherr3        ! error variances    \ if rel. humidity below 20%

USE data_obs_cdfin, ONLY :  &

!         5.2    Temperature / humidity / pressure / height / fog limits
!                -------------------------------------------------------
    rttlim     ,& ! temperature limit below which rel. hum. error is increased
    rrhlim     ,& ! rel. hum. limit below which rel. hum. error is increased
    rtshlm     ,& ! gross error upper limit for relative humidity
!   rerrpp     ,& ! msl pressure above which observed pressure is
!                 ! reckoned to be erroneous
    pminsigt   ,& ! significant-level TEMP / PILOT data are neglected unless
    pminsigv   ,& !   p-obs >= pminsigt (10000 pa) and temperature obs exists or
                  !   p-obs >= pminsigv (20000 pa) and wind obs exists
    pqmin      ,& ! pressure [pa] of level above which moisture obs are not used
    rpplim     ,& ! pressure level obove which observations are not used
!   vfoglim    ,& ! visibility threshold [m] below which the existence of low
!                 !      cloud (fog) is assumed in the presence of precipitation
    fflim_vad  ,& ! lower limit for accepting VAD wind speed
    picoshift  ,& ! shift [pa] applied to colocated ICOS tower wind obs

!         7.0    For data rejection messages: Output buffer, size and formats
!                ------------------------------------------------------------
    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
!   nfmt3      ,& ! no accepted data
    nfmt4      ,& ! excess of levels
    nfmt5         ! several surface levels

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.     Header formats of ODR reports
!       ------------------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
!   mxrhed     ,& ! header length of multi-level reports
!   mxshed     ,& ! header length of single-level reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
    nhsurf     ,& ! height of model grid pt. to which obs. is assigned
!   nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
!   nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
    nhssos     ,& ! standard devation of sub-grid scale model orogrophy [m]
!   nhvcbu     ,& ! correction factor to vertical correlation scale for wind
                  ! at base of report
!   nhvcbt     ,& ! as 'nhvcbu', but for temperature
!   nhvcbq     ,& ! as 'nhvcbu', but for humidity
!   nhvctu     ,& ! correction factor to vertical correlation scale for wind
                  ! at top of report
!   nhvctt     ,& ! as 'nhvctu', but for temperature
!   nhvctq     ,& ! as 'nhvctu', but for humidity

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
!   nhitot     ,& ! global x-coord. of grid pt. assigned to obs
!   nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
!   nhqcfw     ,& ! status of QC and of writing to feedobs files (and p-QC flag)
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
    nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
!   nhcent     ,& ! originating centre
!   nhstid     ,& ! station identity number
!   nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
    nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
    nhuexi     ,& ! flag for existence of wind data        in multi-level report
    nhtexi     ,& ! flag for existence of temperature data in multi-level report
    nhqexi     ,& ! flag for existence of humidity data    in multi-level report
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
    nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
    nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
    nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
    nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]
!   nhstyp     ,& ! surface obs: station type (buoy: MQOBL, BUFR Table 002149,
                  !                            else: NIX  , BUFR Table 002001)

!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid     ,& ! character length of the station identity
!   ilstidp    ,& ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
    nvexbp     ,& ! bit pos. for flag: 'obs or code type excluded   nhschr
                  !                     at station location'
    nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvlzbp     ,& ! bit pos. for vertical levels reported as height   "
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
    nvapoc     ,& ! no. of bits occ. by phase of flight               "
    nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
    nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
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
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control
    nbtzio     ,& ! longitude in geographical coordinates
    nbtzjo     ,& ! latitude  in geographical coordinates
    nbttim     ,& ! observation time relative to report (header) time
    nbtlop     ,& ! LOG( pressure )
    nbtplv     ,& ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbtdrh     ,& ! bias correction for relative humidity [/]
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! profiler: signal to noise ratio
                  ! (ICOS) tower: 0.5* deviation (deg.) of obs device from lee
    nbtuac     ,& ! profiler: wind accuracy (std dev from data provider) [m/s]
                  ! (ICOS) tower: azimuth of wind sensor


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
    mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
!   mxsbdy     ,& ! body length of single-level reports
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
!   nbspst     ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    nbscbs     ,& ! (lowest) cloud base height CBH above surface       [m]
    nbscl      ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm      ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
!   nbsch      ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
!   nbsct      ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
!   nbsvis     ,& ! (horizontal) visibility                            [m]
!   nbsrr1     ,& ! precipitation amount over 1 hour                   [mm]
!   nbsrr6     ,& ! precipitation amount over 6 hours                  [mm]
!   nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
!   nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
!   nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
!   nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
!   nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
!   nbstn      ,& ! minimum temperature (at 2m during past 12 hrs)     [K]
!   nbstx      ,& ! maximum temperature (at 2m during past 12 hrs)     [K]
    nbsrad     ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
!   nbshsw     ,& ! total snow depth                                   [m]
    nbsplv     ,& ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbsdrh     ,& ! bias correction for relative humidity [/]
    nbsvip     ,& ! model surface pressure (for LETKF localisation)    [Pa]

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
    nbswwe     ,& ! NetCDF read, SYNOP: weather and ground group word  (below)
!   nbstur     ,& ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)
    nbsclg        ! general           cloud       group (code)
!   nbscl1     ,& ! first  individual cloud layer group (code)
!   nbscl2     ,& ! second individual cloud layer group (code)
!   nbscl3     ,& ! third  individual cloud layer group (code)
!   nbscl4        ! forth  individual cloud layer group (code)

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status/QC flags for horiz. wind   nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature         "
    nvrq       ,& ! bit pos. for status/QC flags for humidity            "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height     "
!   nvrw       ,& ! bit pos. for status/QC flags for vertical wind       "
    nvfubp     ,& ! bit pos. for main flag on wind                     nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature                "
    nvfqbp     ,& ! bit pos. for main flag on humidity                   "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfaoc     ,& ! no. of bits occ. by each main flag                   "
    nvfbps     ,& ! bit pattern for main flags:                          "
    nvfboc     ,& ! no. of bits occ. for main flags                      "
    nvflbp     ,& ! bit pos. for level flag: level below surface         "
    nvlidp        ! level id. bit pattern                              nb?lid
!   nvlido        ! no. bits occ. by each indicator in level id.         "

USE data_obs_record, ONLY :   &

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------
                  ! combined cloud and weather group                   nbscwg
    nvchbp     ,& ! bit position for ch  (type of high cloud)
    nvcmbp     ,& !         "        cm  (type of middle cloud)
    nvclbp     ,& !         "        cl  (type of low cloud)
    nvnhbp     ,& !         "        nh  (cover of low, else of middle cloud)
    nvhbp      ,& !         "        h   (cloud base height)
!   nvnbp      ,& !         "        n   (total cloud cover)
!   nvwwbp     ,& !         "        ww  (present weather)
                  !                      (see VUB WMO Code tables:)
    nvchoc     ,& ! no. of bits occupied by ch    [Code table 0509]
    nvcmoc     ,& !           "             cm    [Code table 0515]
    nvcloc     ,& !           "             cl    [Code table 0513]
    nvnhoc     ,& !           "             nh    [Code table 2700]
    nvhoc      ,& !           "             h     [Code table 1600]
    nvnoc      ,& !           "             n     [Code table 2700]
!   nvwwoc     ,& !           "             ww    [Code table 4677]
                  ! --> weather and ground group word                  nbswwe
    nvcqbp     ,& ! bit position for refined quality flag on ccl
    nvcqoc     ,& ! no. bits occupied by refined quality flag on ccl
                  !
                  ! --> general    cloud       group word              nbsclg
    nctlbp     ,& ! bit position for low    cloud type code (WMO Table 020012)
    nctmbp     ,& ! bit position for middle cloud type code (WMO Table 020012)
    ncthbp     ,& ! bit position for high   cloud type code (WMO Table 020012)
                  !
                  ! --> individual cloud layer group words             nbscl?
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxclbp     ,& ! bit position for cloud amount      code (WMO Table 020011)
    nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
    nxcloc     ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc        ! no. bits occupied for cloud type   code (WMO Table 020012)

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    fdoro      ,& ! scaling factors to obs height minus model orography

!       2.     Observation data records (ODR)
!       -------------------------------------

    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    momlbd     ,& ! body   of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    mosgbd     ,& ! body   of single-level ODR
    mosghd     ,& ! header of single-level ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR

!       3.     Masking constants
!       ------------------------

    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_OBSTYPE    ,& ! passive report type (at obs. location)
    FL_MERGE      ,& ! merged reports (e.g. TEMP ABCD)
    FL_NO_OBS     ,& ! no (active) observations in report
    LS_SURFACE    ,& ! surface
    LS_STANDARD   ,& ! standard level
    LS_TROPO      ,& ! tropopause level
    LS_MAX        ,& ! maximum wind level
    LS_SIGN       ,& ! significant level
    LS_SUPEROBS      ! superobservation (layer average)

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
!   global_values,   & ! computes global values by operating on local arrays
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
!   phi2phirot,      & ! converts phi from the real to the rotated system
!   rla2rlarot,      & ! converts lambda from the real to the rotated system
    uv2uvrot_vec,    & ! converts wind components from normal to rotated system
    uv2df              ! converts the wind components to wind direction and speed

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_sort_node   ,& ! assign node to reports and sort them accordingly
    obs_cdf_distrib_reports,& ! distribute reports to appropriate sub-domains
    obs_average_layers     ,& ! compute superobbed profile from high-res profile
    obs_rm_duplicate_levels,& ! remove duplicate levels
    std_atmosphere         ,& ! convert variables accord. to standard atmosphere
    std_atm_p2z_fast       ,& ! convert p to z    accord. to standard atmosphere
    obs_td2rh              ,& ! convert dewpoint temperature to relat. humidity
    obs_qx2rh              ,& ! convert mixing ratio to model-compatible rel hum
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
PUBLIC :: obs_cdf_read_temp_pilot, obs_cdf_read_profiler, obs_cdf_store_multilev

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_temp_pilot ( min_sta , min_end , ilcdf                 &
                                   , nodrnew , nexceed)

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" organizes the reading,
!   pre-selection, distribution and storage in ODR of radiosonde (TEMP and
!   PILOT) reports from a given observation time period.
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
    nodrnew (2)   ,& ! number of new reports
    nexceed (2)      ! number of reports in excess of array size

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
                      !  = ncdf_temp      : TEMP
                      !  = ncdf_tempship  : TEMPSHIP
                      !  = ncdf_temphirs  : high-resolution BUFR TEMP
                      !  = ncdf_tempdesc  : descending TEMP
                      !  = ncdf_pilot     : PILOT  (height levels)
                      !  = ncdf_pilot_p   : PILOT  (pressure levels)
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
    nrealml        ,& ! number of multi-level real elements
    inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    irps           ,& ! loop index over reports to be processed
    irep           ,& ! record index over reports to be processed
    ilev           ,& ! loop index over vertical levels
    nlev           ,& ! number of vertical levels
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    natts          ,& ! number of attributes
    iztype         ,& ! variable type (int or float) of element 'NHHHN'
    iptype         ,& ! variable type (int or float) of element 'MPN'
    jj             ,& ! other loop index
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    status , nfsta2,& ! NetCDF status variable
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
    ktimsig        ,& ! significance of report time
    istat  , ierr     ! error indicators

  INTEGER        ::  &          ! variable ID's in NetCDF file for:
    varid_nlev, varid_lvtyp, varid_z   ,& ! number of levels, level type, height
    varid_dlat, varid_dlon , varid_dt  ,& ! latitude / longitude / time shifts
    varid_p   , varid_t    , varid_td  ,& ! pressure, temperature, dew point
    varid_dd  , varid_ff   , varid_fi  ,& ! wind direction/ speed / geopotential
                varid_clsig            ,& ! vertical significance for cloud
    varid_mnh , varid_nh               ,& ! low cloud amount, cloud base height
    varid_mccl, varid_mccm , varid_mcch,& ! type of low / mid-level / high cloud
    varid_na4 , varid_rstyp            ,& ! measurem. equipment, radiosonde type
    varid_rad , varid_track            ,& ! radiation correction, tracking tech.
                varid_mtisi            ,& ! time significance
                varid_ma_ml               ! a(ny) mandatory multi-level variable

  CHARACTER (LEN=180)      :: &
    yattname          ! name of attribute of a particular NetCDF variable
  CHARACTER (LEN=9  )      :: &
    ydrep             ! name of variable for replication factor for profile obs
  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=150)      :: &
    yerrmsl           ! error message
  CHARACTER (LEN=80)       :: &
    ymsg              ! control message
  CHARACTER (LEN=80)       :: &
    ytxtob            ! observation type in text messages
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
! Begin Subroutine obs_cdf_read_temp_pilot
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_temp_pilot'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

! PRINT *, yroutine, kcdftyp, min_sta, min_end

  IF (kcdftyp == ncdf_temp)      ytxtob = 'TEMP '
  IF (kcdftyp == ncdf_tempship)  ytxtob = 'TEMPSHIP'
  IF (kcdftyp == ncdf_temphirs)  ytxtob = 'TEMPHIRS'
  IF (kcdftyp == ncdf_tempdrop)  ytxtob = 'TEMPDROP'
  IF (kcdftyp == ncdf_tempdesc)  ytxtob = 'TEMPDESC'
  IF (kcdftyp == ncdf_pilot)     ytxtob = 'PILOT'
  IF (kcdftyp == ncdf_pilot_p)   ytxtob = 'PILOT (p-levels)'

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
! to be put into 'Xbufsrt'  (specific for TEMPs here)
! ----------------------------------------------------
  IF (    (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)           &
     .OR. (kcdftyp == ncdf_temphirs)) THEN
    noffscal (1,2)  =  noffscal(1,1) + 11
    noffscal (2,2)  =  noffscal(2,1) + 1
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF ((kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
    noffscal (1,2)  =  noffscal(1,1) + 6
    noffscal (2,2)  =  noffscal(2,1) + 0
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF ((kcdftyp == ncdf_pilot) .OR. (kcdftyp == ncdf_pilot_p)) THEN
    noffscal (1,2)  =  noffscal(1,1) + 5
    noffscal (2,2)  =  noffscal(2,1) + 0
    noffscal (3,2)  =  noffscal(3,1) + 0
  ENDIF

!-------------------------------------------------------------------------------
! Section 1: Read and pre-process header information common to all obs types:
!            time / lat / lon / station altitude / obs type / station ID /
!            center of origin ...;
!            assign observations to grid points and to nodes (sub-domains)
!-------------------------------------------------------------------------------

  nrep    = 0

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

! 'irproc' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get entries specific to observation type (TEMP)
!-------------------------------------------------------------------------------

    ALLOCATE ( nc_nlev  (nrep) , STAT=istat )
    ALLOCATE ( nc_tisi  (nrep) , STAT=istat )
    ALLOCATE ( nc_rstyp (nrep) , STAT=istat )
    ALLOCATE ( nc_track (nrep) , STAT=istat )
    ALLOCATE ( nc_na4   (nrep) , STAT=istat )
    DO irep = 1 , nrep
      nc_rstyp(irep)  =  imiss
      nc_track(irep)  =  imiss
      nc_na4  (irep)  =  imiss
      nc_tisi (irep)  =  imiss
!     nc_nlev (irep)  =  imiss
    ENDDO
    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)) THEN
      ALLOCATE ( nc_rad   (nrep) , STAT=istat )
      ALLOCATE ( nc_clsig (nrep) , STAT=istat )
      ALLOCATE ( nc_clclm (nrep) , STAT=istat )
      ALLOCATE ( rc_cbase (nrep) , STAT=istat )
      ALLOCATE ( nc_ccl   (nrep) , STAT=istat )
      ALLOCATE ( nc_ccm   (nrep) , STAT=istat )
      ALLOCATE ( nc_cch   (nrep) , STAT=istat )
      DO irep = 1 , nrep
        nc_rad   (irep)  =  imiss
        nc_clsig (irep)  =  imiss
        nc_clclm (irep)  =  imiss
        rc_cbase (irep)  =  rmiss
        nc_ccl   (irep)  =  imiss
        nc_ccm   (irep)  =  imiss
        nc_cch   (irep)  =  imiss
      ENDDO
    ELSEIF ((kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
      ALLOCATE ( nc_rad   (nrep) , STAT=istat )
    ENDIF

! get header info on type of radiosonde / tracking / radiation correction ...
! ---------------------------------------------------------------------------
! Time significance: 2: time averaged; 3: accumulated; 18: radiosonde launch t.;
!   (MTISI, 008021)  23: monitoring period; 25: nominal reporting time;
!                    26: time of last known position; 31: missing value
!                    (for other values, see WMO BUFR code Table 008021)
! Radiosonde type  : 17: Graw G. (D); 26: Meteolabor Basora (CH);
!   (NRARA)          61: Vaisala RS80 Loran, Digicora I,II or Marwin;
!                    71: Vaisala RS90 Digicora I,II or Marwin; 34: Vinohrady(CZ)
!                    79: Vaisala RS90 Digicora I,II or Marwin;
!                    80: Vaisala RS92 Digicora III; 81: Vaisala RS92 Autosonde
!                    255: missing value
!                    (for other values, see WMO common code Table C2)
! Tracking technique:0: no windfinding; 3: automatic with auxilliary ranging;
!   (NSASA)          2: automatic with auxilliary radio detection finding;
! (Common Table C7)  6: automatic cross chain Loran-C; 8: automatic satellite
!                       navigation; 19: tracking technique not specified;
!                    127: missing value;
!                    (for other values, see WMO common code Table C7)
! Type of measuring: 0: pressure instr. associat. with wind measuring equipment;
!    equipment used: 1: optical theodolite; 2: radio theodolite; 3: radar;
!   (NA4, 002003)    4: VLF-Omega; 5: Loran-C; 6: wind profiler;
!                    7: satellite navigation; 8: RASS; 9: Sodar; 15: missing
! Solar & infrared : 0: no correction; 1: CIMO solar + CIMO infrared corrected;
! radiation correct: 2: CIMO solar + infrared corr.; 3: CIMO solar corr. only;
!   (NSR, 002013)    4: solar + infrared corr. automatically by radiosonde syst;
!                    5: solar corr. autom. by radios.; 6: solar + infrared corr.
!                       as specified by country; 7: solar corr. by country;
!                    15: missing

! treat radiosonde type as mandatory variable
                            status = nf90_inq_varid (ncid,'NRARA',varid_rstyp)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'RADIOSONDE TYPE NRARA DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_rstyp, nc_rstyp, (/mrepsta/), (/nrep/))

! treat other variables as optional
    status = nf90_inq_varid (ncid,'NSASA',varid_track)
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_track, nc_track, (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'NA4'  ,varid_na4  )
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_na4  , nc_na4  , (/mrepsta/), (/nrep/))
    status = nf90_inq_varid (ncid,'MTISI',varid_mtisi)
    IF (status == nf90_noerr)                                                  &
      status = nf90_get_var (ncid, varid_mtisi, nc_tisi , (/mrepsta/), (/nrep/))
    IF (     (kcdftyp == ncdf_temp )    .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)                                        &
        .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
      status = nf90_inq_varid (ncid,'NSR'  ,varid_rad  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_rad, nc_rad  , (/mrepsta/), (/nrep/))
    ENDIF

! control message if some / all TEMPs report nominal time instead of launch time
    ktimsig = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      IF (nc_tisi(irep) == 18) ktimsig = 2* (ktimsig /2) + 1
      IF (nc_tisi(irep) == 25) ktimsig = MOD(ktimsig, 2) + 2
      IF (ktimsig == 2) THEN
        WRITE( ymsg,'("NOTE: all ",A ,"s report NOMINAL TIME only (forecast"   &
                    &," minutes:",I5," - ",I5,")")' )                          &
               ytxtob(1:LEN_TRIM(ytxtob)), min_sta, min_end
      ELSEIF (ktimsig == 3) THEN
        WRITE( ymsg,'("NOTE: some ",A ,"s report NOMINAL TIME only (forecast"  &
                    &," minutes:",I5," - ",I5,")")' )                          &
               ytxtob(1:LEN_TRIM(ytxtob)), min_sta, min_end
      ENDIF
      IF (ktimsig >= 2)  PRINT         '(A)' , ymsg
      IF (ktimsig >= 2)  WRITE( nurej ,'(A)' ) ymsg
    ENDDO

! get cloud information from TEMP: treat as optional
! -------------------------------

    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)) THEN
      status = nf90_inq_varid (ncid,'MVTSU',varid_clsig)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_clsig, nc_clsig, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MNH'  ,varid_mnh  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mnh  , nc_clclm, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NH'   ,varid_nh   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nh   , rc_cbase, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC'  ,varid_mccl )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mccl , nc_ccl  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC0' ,varid_mccm )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mccm , nc_ccm  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC1' ,varid_mcch )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mcch , nc_cch  , (/mrepsta/),(/nrep/))
    ENDIF

! get number of vertical levels
! -----------------------------
    ! name of variable for replication factor used for profile obs
    !   (usually, this should be 'MEDRE' for TEMP/PILOT, but may also be 'MDREP')
    ydrep = 'MEDRE    '
    !   get variable ID of a mandatory (!) mulit-level field
                             nfsta2 = nf90_inq_varid (ncid,'NFNFN',varid_ma_ml )
    IF (nfsta2 == nf90_noerr) THEN
      status = nf90_Inquire_Variable ( ncid, varid_ma_ml, ndims=ndims          &
                                     , dimids=dimids, natts=natts )
      IF (status == nf90_noerr) THEN
        DO jj = 1, natts
          status = nf90_inq_attname ( ncid, varid_ma_ml, jj, yattname )
          IF (TRIM( yattname ) == 'dim1_length') &
            status = nf90_get_att ( ncid, varid_ma_ml, yattname, ydrep )
        ENDDO
        IF (TRIM( ydrep ) /= 'MDREP')  ydrep = 'MEDRE    '
      ENDIF
    ENDIF
!   PRINT *,'ydrep: ', ydrep
    ! get replication factor, i.e. number of vertical levels
    status = nf90_inq_varid (ncid, ydrep, varid_nlev)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'NUMBER OF VERTICAL LEVELS DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_nlev, nc_nlev, (/mrepsta/), (/nrep/))
    ! get maximum number of vertical levels (according to entry 'MEDRE')
    maxlev = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      maxlev = MAX( nc_nlev(irep) , maxlev )
    ENDDO
    IF ((maxlev >= 1) .AND. (nfsta2 /= nf90_noerr)) THEN
      yerrmsl = 'mandatory variable "NFNFN" is NOT correct IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF

! check maximum number of vertical levels
! ---------------------------------------

    IF (maxlev >= 1) THEN
!                              status =nf90_inq_varid (ncid,'MEVSS',varid_lvtyp)
!     IF (status /= nf90_noerr) THEN
!       yerrmsl = 'VERT. SOUND. SIGNIF. MEVSS DOES NOT EXIST IN ' // yfn
!       CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
!     ENDIF
! get length of vertical dimension in NetCDF file
!     status = nf90_Inquire_Variable ( ncid, varid_lvtyp, ndims=ndims          &
!                                    , dimids=dimids)
      dimid_mxlv = dimids(1)
      status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
      IF ((maxlev > mxlev) .OR. (ndims /= 2)) THEN
        PRINT *, 'ERROR with dimensions in file ', ndims, mxlev, maxlev
        yerrmsl = 'NUMBER OF LEVELS EXCEEDS DIMENSION IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
    ENDIF

! get profile observations
! ------------------------

!   PRINT *,'varid3 ', ncid, maxlev, ncdf_temp, kcdftyp, nrep
    IF (maxlev >= 1) THEN
      ALLOCATE ( nc_lvtyp (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_dt    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_dd    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_z     (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_pp    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_dlat  (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_dlon  (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_ff    (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_p     (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_zz    (maxlev,nrep) , STAT=istat )
!     nc_z    (:,:) = imiss            !  need not exist for (ncdf_pilot_p)
      rc_p  (:,:)  =  rmiss
      rc_zz (:,:)  =  rmiss
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        ALLOCATE ( rc_t     (maxlev,nrep) , STAT=istat )
        ALLOCATE ( rc_td    (maxlev,nrep) , STAT=istat )
      ELSEIF ((kcdftyp == ncdf_pilot ) .OR. (kcdftyp == ncdf_pilot_p)) THEN
        ALLOCATE ( rc_fi    (maxlev,nrep) , STAT=istat )
        rc_fi (:,:)  =  rmiss
!       rc_p  (:,:)  =  rmiss
        nc_z  (:,:)  =  imiss
      ENDIF

! 1. variables treated as mandatory
      iztype    = -999
      iptype    = -999
                               status =nf90_inq_varid (ncid,'NDNDN',varid_dd   )
      IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NFNFN',varid_ff   )
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MEVSS',varid_lvtyp)
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'NHHHN',varid_z  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MPN'  ,varid_p  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MTDBT',varid_t  )
        IF(status == nf90_noerr) status =nf90_inq_varid (ncid,'MTDNH',varid_td )
        IF(status == nf90_noerr)                                               &
          status = nf90_Inquire_Variable ( ncid, varid_z, xtype=iztype )
        IF(status == nf90_noerr)                                               &
          status = nf90_Inquire_Variable ( ncid, varid_p, xtype=iptype )
      ENDIF
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'STANDARD '// ytxtob // ' PROFILE DATA DO NOT EXIST IN '// yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF

      nsta2 (1) = 1
      ncnt2 (1) = maxlev
      nsta2 (2) = mrepsta
      ncnt2 (2) = nrep
      status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        IF (iztype == nf90_int) THEN
          status = nf90_get_var (ncid, varid_z  , nc_z    , nsta2, ncnt2)
        ELSEIF (iztype == nf90_float) THEN
          status = nf90_get_var (ncid, varid_z  , rc_zz   , nsta2, ncnt2)
        ENDIF
        IF (iptype == nf90_int) THEN
          status = nf90_get_var (ncid, varid_p  , nc_pp   , nsta2, ncnt2)
        ELSEIF (iptype == nf90_float) THEN
          status = nf90_get_var (ncid, varid_p  , rc_p    , nsta2, ncnt2)
        ENDIF
!       status = nf90_get_var (ncid, varid_z    , nc_z    , nsta2, ncnt2)
!       status = nf90_get_var (ncid, varid_p    , rc_p    , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_lvtyp, nc_lvtyp, nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_t    , rc_t    , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_td   , rc_td   , nsta2, ncnt2)
!       IF (status /= nf90_noerr) PRINT *,'ppm3 ', TRIM(nf90_strerror(status))
      ENDIF
!     PRINT *,'press0 ', varid_p, nsta2, ncnt2, maxlev, nrep
!     PRINT '("press1 ",10F8.0)' , (rc_p(ilev,1),ilev=1,5)                     &
!                                , (rc_p(1,irep),irep=1,5)

! 2. group of variables at least one of which must exist
      IF ((kcdftyp == ncdf_pilot ) .OR. (kcdftyp == ncdf_pilot_p)) THEN
                                status = nf90_inq_varid (ncid,'NHHH ',varid_z  )
        IF (status /= nf90_noerr) varid_z      = imiss
                                status = nf90_inq_varid (ncid,'NHNHN',varid_fi )
        IF (status /= nf90_noerr) varid_fi     = imiss
                                status = nf90_inq_varid (ncid,'MPN  ',varid_p  )
        IF (status /= nf90_noerr) varid_p      = imiss
        IF (varid_z     /= imiss)                                              &
          status = nf90_get_var (ncid, varid_z   , nc_z   , nsta2, ncnt2)
        IF (varid_fi    /= imiss)                                              &
          status = nf90_get_var (ncid, varid_fi  , rc_fi  , nsta2, ncnt2)
        IF (varid_p     /= imiss)                                              &
          status = nf90_get_var (ncid, varid_p   , rc_p   , nsta2, ncnt2)
        IF (      (varid_z == imiss) .AND. (varid_fi == imiss)                 &
            .AND. (varid_p == imiss)) THEN
          yerrmsl = 'NO VERTICAL COORDINATE EXISTS IN ' // yfn
          CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
        ENDIF
      ENDIF

! 3. variables treated as optional
      nc_dt   (:,:) = imiss
      rc_dlat (:,:) = rmiss
      rc_dlon (:,:) = rmiss
      status = nf90_inq_varid (ncid,'NLTPD',varid_dt   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dt   , nc_dt   , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MLADH',varid_dlat )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dlat , rc_dlat , nsta2, ncnt2)
      status = nf90_inq_varid (ncid,'MLODH',varid_dlon )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dlon , rc_dlon , nsta2, ncnt2)
      IF ((kcdftyp == ncdf_pilot ) .OR. (kcdftyp == ncdf_pilot_p)) THEN
        nc_lvtyp(:,:) = imiss
        status = nf90_inq_varid (ncid,'MEVSS',varid_lvtyp)
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_lvtyp, nc_lvtyp, nsta2, ncnt2)
      ENDIF

! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          IF (iztype == nf90_int) THEN
!NEC$ ivdep
!NEC$ nomove
            DO ilev = 1 , nc_nlev(irep)
              IF (nc_z(ilev,irep) /= imiss)                                    &
                rc_zz (ilev,irep)  =  REAL( nc_z(ilev,irep) , wp )
            ENDDO
          ENDIF
          IF (iptype == nf90_int) THEN
!NEC$ ivdep
!NEC$ nomove
            DO ilev = 1 , nc_nlev(irep)
              IF (nc_pp(ilev,irep) /= imiss)                                   &
                rc_p (ilev,irep)  =  REAL( nc_pp(ilev,irep) , wp )
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      IF ((kcdftyp == ncdf_pilot ) .OR. (kcdftyp == ncdf_pilot_p)) THEN
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
!NEC$ ivdep
!NEC$ nomove
          DO ilev = 1 , nc_nlev(irep)
            IF (nc_z(ilev,irep) /= imiss) THEN
              rc_zz (ilev,irep)  =  REAL( nc_z(ilev,irep) , wp )
            ELSEIF (ABS(rc_fi(ilev,irep)) < rmisschk) THEN
              rc_zz (ilev,irep)  =        rc_fi(ilev,irep) / r_g
            ENDIF
!           IF (      (    nc_z (ilev,irep) == imiss)                          &
!               .AND. (ABS(rc_fi(ilev,irep)) < rmisschk))                      &
!             nc_z (ilev,irep)  =  NINT( rc_fi(ilev,irep) / r_g )
          ENDDO
        ENDDO
      ENDIF
    ENDIF

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

    IF (   (kcdftyp == ncdf_temp )   .OR.(kcdftyp == ncdf_tempship)            &
       .OR.(kcdftyp == ncdf_temphirs)                                          &
       .OR.(kcdftyp == ncdf_tempdrop).OR.(kcdftyp == ncdf_tempdesc)) nrealml = 7
    IF    ((kcdftyp == ncdf_pilot)   .OR.(kcdftyp == ncdf_pilot_p )) nrealml = 5
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
! determine length of integer / real buffer for each report
      iirlen (irps)  =  niscal +      3 *nc_nlev(irep)
      irrlen (irps)  =  nrscal + nrealml*nc_nlev(irep)
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

! fill the remaining scalar elements into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------------

    iioffs = 0
    iroffs = 0

    DO irps = 1 , nreproc
! the following part is common to all observation types
      irep  =  irsort(irps)
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nlev (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_rstyp(irep)
      ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_track(irep)
      ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_na4  (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_tisi (irep)
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_rad  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 7)  =  nc_clsig(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 8)  =  nc_clclm(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 9)  =  nc_ccl  (irep)
        ibufsrt (iioffs+noffscal(1,1)+10)  =  nc_ccm  (irep)
        ibufsrt (iioffs+noffscal(1,1)+11)  =  nc_cch  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_cbase(irep)
      ELSEIF ((kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_rad  (irep)
      ENDIF

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

    DEALLOCATE ( nc_tisi  , STAT=istat )
    DEALLOCATE ( nc_rstyp , STAT=istat )
    DEALLOCATE ( nc_track , STAT=istat )
    DEALLOCATE ( nc_na4   , STAT=istat )
    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)) THEN
      DEALLOCATE ( nc_rad   , STAT=istat )
      DEALLOCATE ( nc_clsig , STAT=istat )
      DEALLOCATE ( nc_clclm , STAT=istat )
      DEALLOCATE ( rc_cbase , STAT=istat )
      DEALLOCATE ( nc_ccl   , STAT=istat )
      DEALLOCATE ( nc_ccm   , STAT=istat )
      DEALLOCATE ( nc_cch   , STAT=istat )
    ELSEIF ((kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
      DEALLOCATE ( nc_rad   , STAT=istat )
    ENDIF

! fill in the multi-level info into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------

    IF (maxlev >= 1) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev  = nc_nlev(irep)
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_lvtyp(ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_dd   (ilev,irep)
          ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_dt   (ilev,irep)
          rbufsrt (iroffs+nrscal       +ilev)  =  rc_dlat (ilev,irep)
          rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_dlon (ilev,irep)
          rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_ff   (ilev,irep)
          rbufsrt (iroffs+nrscal+3*nlev+ilev)  =  rc_p    (ilev,irep)
          rbufsrt (iroffs+nrscal+4*nlev+ilev)  =  rc_zz   (ilev,irep)
          IF (   (kcdftyp == ncdf_temp    ).OR. (kcdftyp == ncdf_tempship)     &
             .OR.(kcdftyp == ncdf_temphirs)                                    &
             .OR.(kcdftyp == ncdf_tempdrop).OR. (kcdftyp == ncdf_tempdesc)) THEN
            rbufsrt (iroffs+nrscal+5*nlev+ilev)  =  rc_t    (ilev,irep)
            rbufsrt (iroffs+nrscal+6*nlev+ilev)  =  rc_td   (ilev,irep)
          ENDIF
        ENDDO
!       PRINT '("bufa0 ",4I7   )' , irps, irep, nlev, iroffs
!       PRINT '("bufaa ",8F10.2)' , (rbufsrt(iroffs+istat),istat=1,8)
!       PRINT '("bufab ",8F10.2)' , (rbufsrt(iroffs+istat),istat=9,16)
!       PRINT '("bufac ",8F10.2)' , (rbufsrt(iroffs+nrscal+3*nlev+istat),istat=9,16)
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO

      DEALLOCATE ( nc_lvtyp , STAT=istat )
      DEALLOCATE ( nc_dt    , STAT=istat )
      DEALLOCATE ( nc_dd    , STAT=istat )
      DEALLOCATE ( nc_z     , STAT=istat )
      DEALLOCATE ( nc_pp    , STAT=istat )
      DEALLOCATE ( rc_dlat  , STAT=istat )
      DEALLOCATE ( rc_dlon  , STAT=istat )
      DEALLOCATE ( rc_ff    , STAT=istat )
      DEALLOCATE ( rc_p     , STAT=istat )
      DEALLOCATE ( rc_zz    , STAT=istat )
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        DEALLOCATE ( rc_t     , STAT=istat )
        DEALLOCATE ( rc_td    , STAT=istat )
      ELSEIF ((kcdftyp == ncdf_pilot ) .OR. (kcdftyp == ncdf_pilot_p)) THEN
        DEALLOCATE ( rc_fi    , STAT=istat )
      ENDIF
    ENDIF
    DEALLOCATE ( nc_nlev  , STAT=istat )
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
!     PRINT *,'distr2 ', irps, irnode(irps), iirlen(irps), irrlen(irps), iyrlen(irps)
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

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufloc, rbufloc, ybufloc &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufsrt, rbufsrt, ybufsrt &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_temp_pilot
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_temp_pilot


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_profiler ( min_sta , min_end , ilcdf                   &
                                 , nodrnew , nexceed)

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" organizes the reading,
!   pre-selection, distribution and storage in ODR of remote-sensing profiler
!   (wind profiler, RASS radio acoustic sounding system virtual temperature
!   profiler, radar VAD wind profiler) and tower reports from a given
!   observation time period.
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
    nodrnew (2)   ,& ! number of new reports
    nexceed (2)      ! number of reports in excess of array size

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
                      !  = ncdf_wprof     : Wind Profiler
                      !  = ncdf_radar_vad : Radar VAD
                      !  = ncdf_rass      : RASS
                      !  = ncdf_tower     : tower profile
                      !  = ncdf_tower_icos: ICOS tower profile
                      !  = ncdf_wlidar_wp : ground-b wind lidar (wprof template)
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
    nrealml        ,& ! number of multi-level real elements
    nintml         ,& ! number of multi-level int  elements
    inode          ,& ! node to which the current report is assigned
    ipe            ,& ! loop index over nodes
    irps           ,& ! loop index over reports to be processed
    irep           ,& ! record index over reports to be processed
    ivar           ,& ! loop index over (observed) variables
    ilev           ,& ! loop index over vertical levels
    ilvx           ,& ! loop index over generalized levels ((ICOS) tower)
    iazi           ,& ! loop index over azimuths           ((ICOS) tower)
    nlev           ,& ! number of vertical levels
    nlev_n         ,& ! number of generalized levels ((ICOS) tower)
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    mxazi          ,& ! length of azimuth dimension in tower NetCDF file
    maxklv         ,& ! maxlev * mxazi
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    dimid_mxaz     ,& ! dimension ID for azimuths in (ICOS) tower NetCDF file
    status         ,& ! NetCDF status variable
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
    nsta3 (3)      ,& ! start indices       \  of the reports to be read
    ncnt3 (3)      ,& ! number of elements  /  in the external NetCDF
    istat  , ierr     ! error indicators

  INTEGER        ::  & ! variable ID's in NetCDF file for:
    varid_nlev , varid_z      ,& ! number of levels, height (MH)
    varid_dd   , varid_ff     ,& ! wind direction / speed
    varid_w    , varid_tv     ,& ! vertical velocity / virtual temperature
    varid_t    , varid_td     ,& ! temperature / dew point
    varid_rh   , varid_dz     ,& ! relative humidity / height above ground
    varid_sinor               ,& ! signal to noise ratio (NSINOR)
    varid_stdff               ,& ! standard deviation wind speed (NSTDFF)
    varid_qci  , varid_qci2   ,& ! quality indices (MQINZ / NWPQ)
    varid_na4                 ,& ! measurement. equipment
    varid_wce                 ,& ! wind computation enhancement (MWCE)
    varid_mtisi               ,& ! time significance
    varid_dto  , varid_dto2   ,& ! time period of obs. [sec] / [min]
    varid_qct  , varid_qc     ,& ! assoc. field significance / Q-bits
    varid_azi                 ,& ! azimuth of obs device in tower report
    varid_altp , varid_p      ,& ! barometer height / pressure
    varid_radgl                  ! global radiation

  LOGICAL                  ::  &
    lprof                     ,& ! profiler (wind profiler, RASS, radar VAD,
                                 !           or wind lidar)
    ltower                    ,& ! tower profile (tower, ICOS tower)
    ldim3                     ,& ! tower wind variables have 3rd dimension
                                 !   for azimuth
    lexi_tp1                     ! NGGTP1 exists (in tower report)

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
  CHARACTER (LEN=70)       :: &
    yerrmsl           ! error message
  CHARACTER (LEN= icdfinlen + iannexlen)  ::  &
    yfn               ! file name of NetCDF observation input file, with annex

! Local arrays:
! ------------
  INTEGER        , ALLOCATABLE :: &
    klev       (:) ,& ! indices of generalized levels ((ICOS) tower)
    kazi       (:)    ! indices of azimuths           ((ICOS) tower)

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
! Begin Subroutine obs_cdf_read_profiler
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_profiler'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))
  lprof    =       (kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad )   &
              .OR. (kcdftyp == ncdf_rass ) .OR. (kcdftyp == ncdf_wlidar_wp )
  ltower   =       (kcdftyp == ncdf_tower) .OR. (kcdftyp == ncdf_tower_icos)

! PRINT *, yroutine, kcdftyp, min_sta, min_end

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
! to be put into 'Xbufsrt'  (specific for Profilers here)
! -------------------------------------------------------
  IF (lprof) THEN
    noffscal (1,2)  =  noffscal(1,1) + 5
    noffscal (2,2)  =  noffscal(2,1) + 0
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (ltower) THEN
    noffscal (1,2)  =  noffscal(1,1) + 2
    noffscal (2,2)  =  noffscal(2,1) + 3
    noffscal (3,2)  =  noffscal(3,1) + 0
  ENDIF

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

! 'irproc' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get header entries specific to remote-sensing profilers
!-------------------------------------------------------------------------------

    ALLOCATE ( nc_dto   (nrep) , STAT=istat )
    DO irep = 1 , nrep
      nc_dto  (irep)  =  imiss
    ENDDO

    IF (lprof) THEN    !  i.e. .not.ltower
      ALLOCATE ( nc_na4   (nrep) , STAT=istat )
      ALLOCATE ( nc_tisi  (nrep) , STAT=istat )
      ALLOCATE ( nc_wce   (nrep) , STAT=istat )
      ALLOCATE ( nc_dto2  (nrep) , STAT=istat )
      DO irep = 1 , nrep
        nc_na4  (irep)  =  imiss
        nc_tisi (irep)  =  imiss
        nc_wce  (irep)  =  imiss
        nc_dto2 (irep)  =  imiss
      ENDDO

! get header info on type of radiosonde / tracking / radiation correction ...
! ---------------------------------------------------------------------------
! Time significance: 2: time averaged; 31: missing value
!   (MTISI, 008021)  (for other values, see WMO BUFR code Table 008021)
! Type of measuring: 0: pressure instr. associat. with wind measuring equipment;
!    equipment used: 1: optical theodolite; 2: radio theodolite; 3: radar;
!   (NA4, 002003)    4: VLF-Omega; 5: Loran-C; 6: wind profiler;
!                    7: satellite navigation; 8: RASS; 9: Sodar; 15: missing

      status = nf90_inq_varid (ncid,'NA4'  ,varid_na4  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_na4  , nc_na4 , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTISI',varid_mtisi)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_mtisi, nc_tisi, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MSETP',varid_dto  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dto  , nc_dto , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP',varid_dto2 )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dto2 , nc_dto2, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MWCE' ,varid_wce  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_wce  , nc_wce , (/mrepsta/),(/nrep/))

! if time significance is not 'time averaged' then unset 'time period of obs'
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF ((nc_dto(irep) == imiss) .AND. (nc_dto2(irep) /= imiss))            &
          nc_dto (irep) = nc_dto2(irep) * 60
        IF (      (nc_tisi(irep) /= imiss) .AND. (nc_tisi(irep) /= 31)         &
            .AND. (nc_tisi(irep) /= 2))   nc_dto (irep) = imiss
      ENDDO

! get non-profile info from ICOS tower reports: pressure, global solar radiation
! ------------------------------------------------------------------------------
    ELSEIF (ltower) THEN
      ALLOCATE ( rc_altp  (nrep) , STAT=istat )
      ALLOCATE ( rc_p   (1,nrep) , STAT=istat )
      ALLOCATE ( rc_radgl (nrep) , STAT=istat )
      DO irep = 1 , nrep
        rc_altp  (irep)  =  rmiss
        rc_p   (1,irep)  =  rmiss
        rc_radgl (irep)  =  rmiss
      ENDDO
      status = nf90_inq_varid (ncid,'MHOBNN',varid_altp )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_altp , rc_altp, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MPPP'  ,varid_p    )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_p    , rc_p   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP' ,varid_dto2 )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dto2 , nc_dto , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MGLSR  ',varid_radgl )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_radgl, rc_radgl,(/mrepsta/),(/nrep/))
    ENDIF

! get number of vertical levels ('MDREP' or 'MEDRE' is assumed to exist always)
! -----------------------------

    ALLOCATE ( nc_nlev  (nrep) , STAT=istat )
                              status = nf90_inq_varid (ncid,'MDREP' ,varid_nlev)
    IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MEDRE' ,varid_nlev)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'NUMBER OF VERTICAL LEVELS DOES NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_nlev, nc_nlev, (/mrepsta/), (/nrep/))
! get maximum number of vertical levels (according to entry 'MDREP')
    maxlev = 0
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
      maxlev = MAX( nc_nlev(irep) , maxlev )
    ENDDO

! check maximum number of vertical levels  ('MH', 'MHVEH' is assumed to exist
! ---------------------------------------   (for towers: 'MHOSEN'))
! and check for a 3rd dimension for azimuth (towers)
! -----------------------------------------

    ldim3  = .FALSE.
    IF (maxlev >= 1) THEN
                                status = nf90_inq_varid (ncid,'MH'    ,varid_z)
      IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MHVEH' ,varid_z)
      IF (ltower)               status = nf90_inq_varid (ncid,'MHOSEN',varid_z)
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'LEVEL HEIGHT DOES NOT EXIST IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
! get length of vertical dimension in NetCDF file
      status = nf90_Inquire_Variable ( ncid, varid_z    , ndims=ndims          &
                                     , dimids=dimids)
      dimid_mxlv = dimids(1)
      status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
      IF (maxlev > mxlev) THEN
        PRINT *, 'ERROR with dimensions in file ', ndims, mxlev, maxlev
        yerrmsl = 'NUMBER OF LEVELS EXCEEDS DIMENSION IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
      ! (ICOS) towers: get length of azimuth obs device dimension in NetCDF file
      mxazi = 1
      IF (ltower) THEN
                                status = nf90_inq_varid (ncid,'NFNFN',varid_ff )
        IF (status /= nf90_noerr)  varid_ff  = imiss
                                status = nf90_inq_varid (ncid,'MDA'  ,varid_azi)
        IF (status /= nf90_noerr)  varid_azi = imiss
        ndims = 0
        IF (varid_azi /= imiss) THEN
          status = nf90_Inquire_Variable ( ncid, varid_azi  , ndims=ndims      &
                                         , dimids=dimids)
        ELSEIF (varid_ff /= imiss) THEN
          status = nf90_Inquire_Variable ( ncid, varid_ff   , ndims=ndims      &
                                         , dimids=dimids)
        ENDIF
        ldim3  = (ndims == 3)
        IF (ldim3) THEN
          dimid_mxaz = dimids(1)
          status = nf90_Inquire_Dimension (ncid, dimid_mxaz, len=mxazi)
          ! if   the azimuth dimension 'mxazi' > 1,
          ! then the azimuth variable MDA must exist
          IF ((mxazi > 1) .AND. (varid_azi == imiss)) THEN
            yerrmsl = 'AZIMUTH VARIABLE MDA DOES NOT EXIST IN ' // yfn
            CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

! get profile observations
! ------------------------

!   PRINT *,'varid3 ', ncid, maxlev, kcdftyp, nrep
    IF (maxlev >= 1) THEN
      IF (lprof) THEN
        ALLOCATE ( nc_z     (maxlev,nrep) , STAT=istat )
        ALLOCATE ( rc_w     (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_qci   (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_sinor (maxlev,nrep) , STAT=istat )
        IF (     (kcdftyp == ncdf_wprof    ) .OR. (kcdftyp == ncdf_radar_vad)  &
            .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
          ALLOCATE ( nc_dd    (maxlev,nrep) , STAT=istat )
          ALLOCATE ( rc_ff    (maxlev,nrep) , STAT=istat )
          ALLOCATE ( rc_stdff (maxlev,nrep) , STAT=istat )
        ELSEIF (kcdftyp == ncdf_rass ) THEN
          ALLOCATE ( rc_tv    (maxlev,nrep) , STAT=istat )
        ENDIF
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          DO ilev = 1 , maxlev
            nc_z    (ilev,irep) = imiss
            rc_w    (ilev,irep) = rmiss
            nc_qci  (ilev,irep) = imiss
            nc_sinor(ilev,irep) = imiss
            IF (     (kcdftyp == ncdf_wprof ) .OR. (kcdftyp == ncdf_radar_vad) &
                .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
              nc_dd   (ilev,irep) = imiss
              rc_ff   (ilev,irep) = rmiss
              rc_stdff(ilev,irep) = rmiss
            ELSEIF (kcdftyp == ncdf_rass) THEN
              rc_tv   (ilev,irep) = rmiss
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF (ltower) THEN
        maxklv = maxlev *mxazi
        ALLOCATE ( nc_z     (maxklv,nrep) , STAT=istat )
        ALLOCATE ( nc_dd    (maxklv,nrep) , STAT=istat )
        ALLOCATE ( rc_ff    (maxklv,nrep) , STAT=istat )
        ALLOCATE ( rc_azi   (maxklv,nrep) , STAT=istat )
        ALLOCATE ( rc_t     (maxklv,nrep) , STAT=istat )
        ALLOCATE ( rc_td    (maxklv,nrep) , STAT=istat )
        ALLOCATE ( nc_uuu2  (maxklv,nrep) , STAT=istat )
        ALLOCATE ( nc_qci   (maxklv,nrep) , STAT=istat )
        ALLOCATE ( rc_dz    (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_tisi2 (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_peri  (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_qct   (maxlev,nrep) , STAT=istat )
        ALLOCATE ( nc_qc    (maxlev,nrep) , STAT=istat )
        IF (ldim3) THEN
          ALLOCATE ( nc_dd3 (mxazi,maxlev,nrep) , STAT=istat )
          ALLOCATE ( rc_ff3 (mxazi,maxlev,nrep) , STAT=istat )
          ALLOCATE ( rc_azi3(mxazi,maxlev,nrep) , STAT=istat )
        ENDIF
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          DO ilev = 1 , maxklv
            nc_z    (ilev,irep) = imiss
            nc_dd   (ilev,irep) = imiss
            rc_ff   (ilev,irep) = rmiss
            rc_azi  (ilev,irep) = rmiss
            rc_t    (ilev,irep) = rmiss
            rc_td   (ilev,irep) = rmiss
            nc_uuu2 (ilev,irep) = imiss
            nc_qci  (ilev,irep) = 0
          ENDDO
          DO ilev = 1 , maxlev
            rc_dz   (ilev,irep) = rmiss
            nc_tisi2(ilev,irep) = imiss
            nc_peri (ilev,irep) = imiss
            nc_qct  (ilev,irep) = imiss
            nc_qc   (ilev,irep) = imiss
          ENDDO
        ENDDO
        IF (ldim3) THEN
          DO irps = 1 , nreproc
            irep  =  irsort(irps)
            DO ilev = 1 , maxlev
              DO iazi = 1 , mxazi
                nc_dd3  (iazi,ilev,irep) = imiss
                rc_ff3  (iazi,ilev,irep) = rmiss
                rc_azi3 (iazi,ilev,irep) = rmiss
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      nsta2 (1) = 1
      ncnt2 (1) = maxlev
      nsta2 (2) = mrepsta
      ncnt2 (2) = nrep
      IF (lprof) THEN
        !   for profiler data, height is integer
        status = nf90_get_var   (ncid, varid_z    , nc_z    , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'MWMPS' ,varid_w    )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_w    , rc_w    , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'MQINZ' ,varid_qci  )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_qci  , nc_qci  , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'NSINOR',varid_sinor)
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_sinor, nc_sinor, nsta2, ncnt2)
      ENDIF
      IF (     (kcdftyp == ncdf_wprof    ) .OR. (kcdftyp == ncdf_radar_vad)    &
          .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
        status = nf90_inq_varid (ncid,'NDNDN' ,varid_dd   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'NFNFN' ,varid_ff   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
        status = nf90_inq_varid (ncid,'NSTDFF',varid_stdff)
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_stdff, rc_stdff, nsta2, ncnt2)
      ELSEIF (kcdftyp == ncdf_rass) THEN
        status = nf90_inq_varid (ncid,'MTVIR' ,varid_tv   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_tv   , rc_tv   , nsta2, ncnt2)
      ENDIF

      IF (ltower) THEN
        IF (ldim3) THEN
          nsta3 (1) = 1
          ncnt3 (1) = mxazi
          nsta3 (2) = 1
          ncnt3 (2) = maxlev
          nsta3 (3) = mrepsta
          ncnt3 (3) = nrep
          status = nf90_inq_varid (ncid,'NDNDN' ,varid_dd   )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_dd   , nc_dd3  , nsta3, ncnt3)
          status = nf90_inq_varid (ncid,'NFNFN' ,varid_ff   )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_ff   , rc_ff3  , nsta3, ncnt3)
          status = nf90_inq_varid (ncid,'MDA'   ,varid_azi  )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_azi  , rc_azi3 , nsta3, ncnt3)
        ELSE
          status = nf90_inq_varid (ncid,'NDNDN' ,varid_dd   )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
          status = nf90_inq_varid (ncid,'NFNFN' ,varid_ff   )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
          status = nf90_inq_varid (ncid,'MDA'   ,varid_azi  )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_azi  , rc_azi  , nsta2, ncnt2)
        ENDIF
        !   for tower data, height (above sensor) is float
        status = nf90_inq_varid (ncid,'MHOSEN',varid_dz   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_dz, rc_dz            , nsta2,ncnt2)
        status = nf90_inq_varid (ncid,'MTDBT ',varid_t    )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_t , rc_t (1:maxlev,:), nsta2,ncnt2)
        status = nf90_inq_varid (ncid,'MTDNH ',varid_td   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_td, rc_td(1:maxlev,:), nsta2,ncnt2)
        status = nf90_inq_varid (ncid,'MUUU  ',varid_rh   )
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_rh,nc_uuu2(1:maxlev,:),nsta2,ncnt2)
        status = nf90_inq_varid (ncid,'MTISI ',varid_mtisi)
        IF (status == nf90_noerr)                                              &
          status = nf90_get_var (ncid, varid_mtisi, nc_tisi2, nsta2, ncnt2)
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
        ENDDO
        ! map onto 1-dim 'nc_tisi':
        ! all integer:
        ! NGGTP MADDF MTDBTQ MTDNHQ MUUUQ MUUU MTISI NGGTP0 MADDF0 NDNDNQ NFNFNQ
        !   1. check if NGGPT1 exists: then it refers to wind, NGGTP0 to T,q
        !                                    (and NGGTP to global radiation),
        !                     otherwise NGGTP0 refers to wind, NGGTP  to T,q
                            status = nf90_inq_varid (ncid,'NGGTP1',varid_dto2 )
        lexi_tp1  =  (status == nf90_noerr)
        DO ivar = 1 , 5
          IF (lexi_tp1) THEN
            IF (ivar <= 2)  status = nf90_inq_varid (ncid,'NGGTP1',varid_dto2 )
            IF (ivar >= 3)  status = nf90_inq_varid (ncid,'NGGTP0',varid_dto2 )
          ELSE
            IF (ivar <= 2)  status = nf90_inq_varid (ncid,'NGGTP0',varid_dto2 )
            IF (ivar >= 3)  status = nf90_inq_varid (ncid,'NGGTP ',varid_dto2 )
          ENDIF
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_dto2 , nc_peri , nsta2, ncnt2)
          IF (ivar <= 2)  status = nf90_inq_varid (ncid,'MADDF0',varid_qct  )
          IF (ivar >= 3)  status = nf90_inq_varid (ncid,'MADDF ',varid_qct  )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_qct  , nc_qct  , nsta2, ncnt2)
          IF (ivar == 1)  status = nf90_inq_varid (ncid,'NDNDNQ',varid_qc   )
          IF (ivar == 2)  status = nf90_inq_varid (ncid,'NFNFNQ',varid_qc   )
          IF (ivar == 3)  status = nf90_inq_varid (ncid,'MTDBTQ',varid_qc   )
          IF (ivar == 4)  status = nf90_inq_varid (ncid,'MUUUQ ',varid_qc   )
          IF (ivar == 5)  status = nf90_inq_varid (ncid,'MTDNHQ',varid_qc   )
          IF (status == nf90_noerr)                                            &
            status = nf90_get_var (ncid, varid_qc   , nc_qc   , nsta2, ncnt2)

! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

          DO irps = 1 , nreproc
            irep  =  irsort(irps)
            DO ilev = 1 , nc_nlev(irep)
              !   set bad reporting practice flag if time averaging > 60 min
              IF (   (      (ABS( nc_peri (ilev,irep) ) > 60)                  &
                      .AND. (     nc_peri (ilev,irep)  /= imiss))              &
                 .OR.(      (     nc_tisi2(ilev,irep)  /= 2)                   &
                      .AND. (     nc_tisi2(ilev,irep)  >= 1)                   &
                      .AND. (     nc_tisi2(ilev,irep)  <= 30)))                &
                nc_qci (ilev,irep) = IBSET( nc_qci(ilev,irep) , ivar )
              !   set dataset flag if flagged bad according to WMO table 031021
              IF (   (     (nc_qct(ilev,irep)== 6)                             &
                     .AND.((nc_qc(ilev,irep) == 3).OR.(nc_qc(ilev,irep) == 4)))&
                 .OR.(     (nc_qct(ilev,irep)== 2)                             &
                     .AND.((nc_qc(ilev,irep) == 2).OR.(nc_qc(ilev,irep) == 3)))&
                 .OR.(    ((nc_qct(ilev,irep)== 1).OR.(nc_qct(ilev,irep)== 8)) &
                     .AND. (nc_qc(ilev,irep) == 1)))                           &
                nc_qci (ilev,irep) = IBSET( nc_qci(ilev,irep) , ivar + 10 )
            ENDDO
          ENDDO
          IF (ivar == 1) THEN
            DO irps = 1 , nreproc
              irep  =  irsort(irps)
              DO ilev = 1 , nc_nlev(irep)
                IF ((nc_dto(irep) ==imiss) .AND. (nc_peri(ilev,irep) /=imiss)) &
                  nc_dto (irep)  =  nc_peri(ilev,irep)
                IF (ABS( rc_dz(ilev,irep) ) < rmisschk)                        &
                  nc_z (ilev,irep)  =  NINT( rc_dz(ilev,irep) )
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE ( rc_dz     , STAT=istat )
        DEALLOCATE ( nc_tisi2  , STAT=istat )
        DEALLOCATE ( nc_peri   , STAT=istat )
        DEALLOCATE ( nc_qct    , STAT=istat )
        DEALLOCATE ( nc_qc     , STAT=istat )
      ENDIF  !  (ltower)

! exploit 'NOAA WIND PROFILER, QUAL. CONTR. RESULTS' bits (BUFR Table 025034)
      status = nf90_inq_varid (ncid,'NWPQ'  ,varid_qci2 )
      IF ((status == nf90_noerr) .AND. (lprof)) THEN
        ALLOCATE ( nc_qci2  (maxlev,nrep) , STAT=istat )
        status = nf90_get_var (ncid, varid_qci2 , nc_qci2 , nsta2, ncnt2)
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          DO ilev = 1 , nc_nlev(irep)
            IF ((nc_qci2(ilev,irep) >= 1) .AND. (nc_qci2(ilev,irep) <= 3))     &
              nc_qci (ilev,irep) = 1
          ENDDO
        ENDDO
        DEALLOCATE ( nc_qci2  , STAT=istat )
      ENDIF

!     PRINT *,'press0 ', varid_p, nsta2, ncnt2, maxlev, nrep
!     PRINT '("press1 ",10F8.0)' , (rc_p(ilev,1),ilev=1,5)                     &
!                                , (rc_p(1,irep),irep=1,5)

! ICOS tower profiles: convert the 2 dimensions of level and azimuth
!                      into 1 dimension of generalized level with index 'klev';
!                      for several (up to 4) 'klev', the height can be the same
!                      but then the azimuth must be different
      IF ((ltower) .AND. (ldim3)) THEN
        ALLOCATE ( klev     (maxklv) , STAT=istat )
        ALLOCATE ( kazi     (maxklv) , STAT=istat )
        DO irps = 1 , nreproc
          irep  =  irsort(irps)
          ilvx  =  0
          DO ilev = 1 , maxlev
!           n_azi (ilev,irep) = 0
            IF (ilev <= nc_nlev(irep)) THEN
              ! at each obs level, there must be at least 1 (wind or other) obs
              ilvx = ilvx + 1 
              klev (ilvx) = ilev
              kazi (ilvx) = 1
!             n_azi (ilev,irep) = 1
              DO iazi = 2 , mxazi
                IF (ABS( rc_azi3(iazi,ilev,irep) ) < rmisschk) THEN
                  ilvx = ilvx + 1
                  klev (ilvx) = ilev
                  kazi (ilvx) = iazi
                ENDIF
!               IF (rc_azi3(iazi,ilev,irep) < rmisschk)  n_azi (ilev,irep) = iazi
              ENDDO
            ENDIF
          ENDDO
          nlev_n = ilvx

          DO ilvx = 1 , nlev_n
            ilev = klev(ilvx)
            iazi = kazi(ilvx)
            ! conversion of 3-D into 2-D arrays
            nc_dd  (ilvx,irep)  =  nc_dd3 (iazi,ilev,irep)
            rc_ff  (ilvx,irep)  =  rc_ff3 (iazi,ilev,irep)
            rc_azi (ilvx,irep)  =  rc_azi3(iazi,ilev,irep)
          ENDDO
          IF (nlev_n > nc_nlev(irep)) THEN
            ! this loop must run top down because ilvx >= ilev
            DO ilvx = nlev_n , 1 , -1
              ! for quantities with already 2-D arrays, the values
              ! are added (duplicated) to all 'klev' with the same height
              ilev = klev(ilvx)
              nc_z   (ilvx,irep)  =  nc_z   (ilev,irep) 
              rc_t   (ilvx,irep)  =  rc_t   (ilev,irep) 
              rc_td  (ilvx,irep)  =  rc_td  (ilev,irep) 
              nc_uuu2(ilvx,irep)  =  nc_uuu2(ilev,irep) 
              nc_qci (ilvx,irep)  =  nc_qci (ilev,irep) 
            ENDDO
          ENDIF
          ! update nc_nlev(irep) !
          nc_nlev (irep)  =  nlev_n
        ENDDO
        DEALLOCATE ( klev     , STAT=istat )
        DEALLOCATE ( kazi     , STAT=istat )
        DEALLOCATE ( nc_dd3   , STAT=istat )
        DEALLOCATE ( rc_ff3   , STAT=istat )
        DEALLOCATE ( rc_azi3  , STAT=istat )
      ENDIF  !  ICOS tower
    ENDIF  !  (maxlev >= 1)

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

    IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad )              &
                                .OR. (kcdftyp == ncdf_wlidar_wp )) nintml  = 4
    IF ((kcdftyp == ncdf_wprof) .OR. (kcdftyp == ncdf_radar_vad )              &
                                .OR. (kcdftyp == ncdf_wlidar_wp )) nrealml = 3
    IF  (kcdftyp == ncdf_rass )  nintml  = 3
    IF  (kcdftyp == ncdf_rass )  nrealml = 2
    IF  (ltower)                 nintml  = 4
    IF  (ltower)                 nrealml = 4
    DO irps = 1 , nreproc
      irep  =  irsort(irps)
! determine length of integer / real buffer for each report
      iirlen (irps)  =  niscal + nintml *nc_nlev(irep)
      irrlen (irps)  =  nrscal + nrealml*nc_nlev(irep)
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

! fill the remaining scalar elements into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------------

    iioffs = 0
    iroffs = 0

    DO irps = 1 , nreproc
! the following part is common to all observation types
      irep  =  irsort(irps)
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nlev (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_dto  (irep)
      IF (lprof) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_wce  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_na4  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_tisi (irep)
      ENDIF
      IF (ltower) THEN
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_altp (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_p  (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_radgl(irep)
      ENDIF

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

    DEALLOCATE ( nc_dto   , STAT=istat )
    IF (lprof) THEN
      DEALLOCATE ( nc_na4   , STAT=istat )
      DEALLOCATE ( nc_tisi  , STAT=istat )
      DEALLOCATE ( nc_dto2  , STAT=istat )
      DEALLOCATE ( nc_wce   , STAT=istat )
    ELSEIF (ltower) THEN
      DEALLOCATE ( rc_altp  , STAT=istat )
      DEALLOCATE ( rc_p     , STAT=istat )
      DEALLOCATE ( rc_radgl , STAT=istat )
    ENDIF

! fill in the multi-level info into the long buffer arrays 'Xbufsrt'
! ------------------------------------------------------------------

    IF (maxlev >= 1) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev = nc_nlev(irep)
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_z    (ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_qci  (ilev,irep)
          IF (lprof) THEN
            ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_sinor(ilev,irep)
            rbufsrt (iroffs+nrscal       +ilev)  =  rc_w    (ilev,irep)
          ENDIF
          IF (     (kcdftyp == ncdf_wprof    ).OR. (kcdftyp == ncdf_radar_vad) &
              .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
            ibufsrt (iioffs+niscal+3*nlev+ilev)  =  nc_dd   (ilev,irep)
            rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_ff   (ilev,irep)
            rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_stdff(ilev,irep)
          ELSEIF (kcdftyp == ncdf_rass ) THEN
            rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_tv   (ilev,irep)
          ENDIF
          IF (ltower) THEN
            ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_uuu2 (ilev,irep)
            ibufsrt (iioffs+niscal+3*nlev+ilev)  =  nc_dd   (ilev,irep)
            rbufsrt (iroffs+nrscal       +ilev)  =  rc_t    (ilev,irep)
            rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_ff   (ilev,irep)
            rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_td   (ilev,irep)
            rbufsrt (iroffs+nrscal+3*nlev+ilev)  =  rc_azi  (ilev,irep)
          ENDIF
        ENDDO
!       PRINT '("bufa0 ",4I7   )' , irps, irep, nlev, iroffs
!       PRINT '("bufaa ",8F10.2)' , (rbufsrt(iroffs+istat),istat=1,8)
!       PRINT '("bufab ",8F10.2)' , (rbufsrt(iroffs+istat),istat=9,16)
!       PRINT '("bufac ",8F10.2)' , (rbufsrt(iroffs+nrscal+3*nlev+istat),istat=9,16)
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO

      IF (lprof) THEN
        DEALLOCATE ( nc_z      , STAT=istat )
        DEALLOCATE ( rc_w      , STAT=istat )
        DEALLOCATE ( nc_qci    , STAT=istat )
        DEALLOCATE ( nc_sinor  , STAT=istat )
        IF (     (kcdftyp == ncdf_wprof    ) .OR. (kcdftyp == ncdf_radar_vad)  &
            .OR. (kcdftyp == ncdf_wlidar_wp)) THEN
          DEALLOCATE ( nc_dd     , STAT=istat )
          DEALLOCATE ( rc_ff     , STAT=istat )
          DEALLOCATE ( rc_stdff  , STAT=istat )
        ELSEIF (kcdftyp == ncdf_rass ) THEN
          DEALLOCATE ( rc_tv     , STAT=istat )
        ENDIF
      ENDIF
      IF (ltower) THEN
        DEALLOCATE ( nc_z      , STAT=istat )
        DEALLOCATE ( nc_qci    , STAT=istat )
        DEALLOCATE ( nc_dd     , STAT=istat )
        DEALLOCATE ( rc_ff     , STAT=istat )
        DEALLOCATE ( rc_azi    , STAT=istat )
        DEALLOCATE ( rc_t      , STAT=istat )
        DEALLOCATE ( rc_td     , STAT=istat )
        DEALLOCATE ( nc_uuu2   , STAT=istat )
      ENDIF
    ENDIF
    DEALLOCATE ( nc_nlev  , STAT=istat )
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
!     PRINT *,'distr2 ', irps, irnode(irps), iirlen(irps), irrlen(irps), iyrlen(irps)
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

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufloc, rbufloc, ybufloc &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_multilev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly  &
                                , noffscal,          ibufsrt, rbufsrt, ybufsrt &
                                , nodrnew , nexceed )
!   ===========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_profiler
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_profiler


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for storing multi-level repo. in ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_store_multilev ( kcdftyp, nrepl, nlenli, nlenlr, nlenly     &
                                  , noffscal      , ibuf  , rbuf  , ybuf       &
                                  , nodrnew, nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" stores the multi-level
!   reports in the internal ODR arrays.
!
! Method:
!   First, the header part common to all observation types is stored by calling
!   'obs_cdf_store_comhead', and then the rest of the header is added which
!   depends on the observation type.
!   Next, after storing also the cloud information from TEMPs in a surface
!   report, the multi-level observations are copied from the observation type
!   dependent buffers into standardised arrays. Using these, it is decided
!   which vertical levels are discarded depending on pressure (or height),
!   and a vertically sorted list of those levels is compiled which are processed
!   further. Subsequently, the flag patterns are built, and the elements are
!   stored in the multi-level ODR.
!   To evaluate the blacklist, first, a list of blacklisted vertical intervals
!   for each of the current reports is produced. Then, within the loop over
!   current reports, a blacklist flag 'kblk' is prepared for each level and
!   variable.
!   In addition, a separate surface report is created to accommodate the
!   surface level. (In the multi-level report, 10-m wind, 2-m temperature and
!   humidity are set passive.)
!   Finally, reports without observations are deleted in the ODR.
!
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
                        !  = ncdf_temp      : TEMP
                        !  = ncdf_tempship  : TEMPSHIP
                        !  = ncdf_temphirs  : high-resolution BUFR TEMP
                        !  = ncdf_tempdrop  : TEMP Dropsonde
                        !  = ncdf_tempdesc  : descending TEMP
                        !  = ncdf_pilot     : PILOT (z-levels)
                        !  = ncdf_pilot_p   : PILOT (p-levels)
                        !  = ncdf_wprof     : WIND PROFILER
                        !  = ncdf_radar_vad : RADAR VAD (wind profiles)
                        !  = ncdf_rass      : RASS (virt. temperature profiler)
                        !  = ncdf_tower     : tower profile
                        !  = ncdf_tower_icos: ICOS tower profile
                        !  = ncdf_wlidar_wp : ground-b. wind lidar (wprof templ)
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
    nodrnew (2)      ,& ! number of new reports
    nexceed (2)         ! number of reports in excess of array size

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    dp_surf = 200._wp,& ! dist. of lowest superobbing layer from surface (2 hPa)
    fbogus  =   0._wp   ! if /= 0 then bogus data are produced by applying
                        ! factor 'fbogus' to pre-defined obs 'corrections'

  LOGICAL        , PARAMETER ::  &
    lmult_shft =.FALSE. ! upper-air observations are not horizontally shifted

  REAL (KIND=wp) , PARAMETER ::  &
    c100r = 0.01_wp     ! maximum number of elements to be received

! Local scalars:
! -------------

  LOGICAL        ::  &
    lblkw  (nrepl+1) ,& ! report blacklisted because sta-ID missing on whitelist
    lflg0  (mxrbdf)  ,& ! .true. then use 'IOR' for superobbing of flags (bits)
    lchange          ,& ! indicates changes in latest loop of bubblesorting
    lseaobs          ,& ! sea report (ship, buoy, not land report)
    lseaonly         ,& ! use actively only observations from sea platforms
    ltempasc         ,& ! ascending TEMP data are processed
    lamdarml         ,& ! multi-level AMDAR data are processed
    lprofw           ,& ! remote-sensing wind        profiler are processed
    lproft           ,& ! remote-sensing temperature profiler are processed
    ltower           ,& ! tower profile reports are processed
    lexitd           ,& ! dewpoint temperature obs exist
    lexirh           ,& ! relative humidity    obs exist
    lexiqx           ,& ! mixing ratio         obs exist
    lbogus           ,& ! .true. if bogus data are produced
    lcheck_zp        ,& ! apply height gross error check for pressure levels
    ltempb           ,& ! probably TEMP part B
    ll1    , ll2        ! local logical switches

  INTEGER        ::  &
    nsgob  , nmlob   ,& ! target ODR report indices (single-/ multi-level reps.)
    iioffs (nrepl+1) ,& ! offset of the integer   elements  \  for each
    iroffs (nrepl+1) ,& ! offset of the real      elements   > local
    iyoffs (nrepl+1) ,& ! offset of the character elements  /  report
    jobtyp (nrepl+1) ,& ! observation type (input array for local blacklist)
    jcdtyp (nrepl+1) ,& ! obs  code   type (input array for local whitelist)
    icma   , icma2   ,& ! pointer index for 'cma' / diagnostic arrays
    niscal , nrscal  ,& ! number of scalar integer / real  elements
    nrepml           ,& ! number of reports which can be stored in the ODR now
    nexceml, nexcesg ,& ! no. of local mult./sfc reports in excess of array size
    nrepadd          ,& ! number of 'merged' reports after superobbing into ODR
    nmladd           ,& ! ODR report index for 'merged' report
    irpl             ,& ! loop index over local reports
    ilev             ,& ! loop index over vertical levels (in input reports)
    klev             ,& ! loop index over vertical levels (in ODR)
    intv             ,& ! loop index over blacklisted vertical intervals
    icl              ,& ! loop index over single characters
    kk               ,& ! loop index over model levels
    istep            ,& ! loop index over directions (ICOS)
    maxlev           ,& ! max. number of vertical levels in all local reports
    nlev             ,& ! number of vertical levels in current report
    nlvp             ,& ! number of vertical levels to be processed (for ODR)
    nlvp2            ,& ! number of vertical levels >= 200 hPa to be processed
    nrejlv           ,& ! number of empty    levels to be rejected
    nlvz             ,& ! number of vertical levels with reported height
    indlev           ,& ! indication of type of level (pressure and/or height)
    kobtyp , kcdtyp  ,& ! observation type / observation code type
    mphase           ,& ! phase of flight                     (WMO Table 008004)
    nccl             ,& ! (low or mid-level) cloud cover    (WMO Table 020011)
    nhkey            ,& ! class of cloud base height (VUB WMO Code table 1600)
    ncsig            ,& ! vertical significance of cloud    (WMO Table 008002)
    nclcl  , nclcm   ,& ! derived low / middle cloud amount            (octas)
    nclqf            ,& ! quality flag for derived low cloud amount
    ncc              ,& ! low/middle/high cloud type (VUB T. 0513, 0515, 0509)
    nclwgn , nclwgr  ,& ! combined cloud and weather group (new/old, for verif)
    iob    , job     ,& ! indices of grid point assigned to observation
    iactr            ,& ! stepping for event counters
    ixrhs               ! index for 'rhtsat'

  INTEGER        ::  &
    ilevsfc          ,& ! level index of surface level (in input array)
    klevsfc          ,& ! level index of surface level (in ODR)
    nlidvof          ,& ! level identity (format of VOF)
    nlidcdf          ,& ! level identity (format of feedback NetCDF)
    nlidin           ,& ! level identity (format of input NetCDF)
    idsurf , idstd   ,& ! indicators for: surface level    / standard level
    idtropo, iduvmax ,& ! indicators for: tropopause level / maximum wind level
    idsigv , idsigq  ,& ! indicators for significant level of: wind / humidity
    idsigt           ,& ! indicators for significant level of: temperature
    iqcbit , iqcbit2 ,& ! quality control bit
    nflag            ,& ! total observation flag in ODR
    klvz             ,& ! level with active height obs
    nzpp             ,& ! pressure [pa]
    ilen             ,& ! length of control message
    ilverr           ,& ! nearest standard error level below observation level
    nhsl             ,& ! (ICOS) tower: number of heights with several levels
    ihsl             ,& ! (ICOS) tower: index over heights with several levels
    ilact            ,& ! (ICOS) tower: index of active level
    kbot   , ktop    ,& ! indices of bottom / top interpolation levels used
    nlimsup          ,& ! lower limit for number of levels for doing superobbing
    nlint            ,& ! number of interpolation levels used for current report
    nlvsup           ,& ! number of superobbed levels
    ictot            ,& ! number of levels below 300 hPa
    icmiss(3)        ,& ! number of levels with missing obs
    necnt            ,& ! index of event counter
    iactx (5)        ,& ! indicates presence of active  data (for variables)
    ipasx (5)        ,& ! indicates presence of passive data (for variables)
    nzaexi           ,& ! -1: only passive data ; 0: no data ; 1: active data
    ilstid_chk       ,& ! length of (part of) station ID that can be checked
    istat  , irm     ,& ! error indicators
    itmp                ! auxilliary buffer

  REAL (KIND=wp) ::  &
    zobhr  (nrepl+1) ,& ! obs  time (input array for local whitelist)
    htop             ,& ! 'top' of model atmosphere [m] in the sense that obs
                        !   levels further above can / should be neglected
                        !   (here: htop = mean of top main and half model level)
    zclz             ,& ! cloud base height
    zstalt , zzaltsy ,& ! station altitude
    roblat , roblon  ,& ! latitude and longitude in the geograhical system
!   rlat   , rlon    ,& ! latitude and longitude in the rotated system
    roberr           ,& ! observation error
    fiperr           ,& ! interpolation weight factor for error level 'ilverr'
    fisd             ,& ! height diff. betw. model orography and sta. altitude
    fisduv           ,& ! modified height difference for 10-m wind data
    fisdtt           ,& ! modified height difference for 2-m temperature data
    fisdrh           ,& ! modified height difference for 2-m humidity data
    fisdzz           ,& ! scaled extrapolation distance for surface press. val.
    zzkeml           ,& ! height of the lowest main model level
    p_bot  , p_top   ,& ! pressure of bottom and of top level in report
    zddmean          ,& ! (ICOS) tower: upwind direction averaged over azimuths
    zdazi            ,& ! (ICOS) tower: difference of azimuth and upwind direct.
    zdazimn          ,& ! (ICOS) tower: minimum 'zdazi'
    zp200            ,& ! 20000. Pa 
    zzmis            ,& ! missing height value
    zvalue           ,& ! any value
    c3600r              ! 1 / 3600.0_wp

! CHARACTER (LEN=25)       :: &
!   yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid_blk)  :: &
    ystidl  (nrepl+1)   ! observation type (input array for local blacklist)

! Local allocatable arrays: (they include the dimension length 'maxlev+1')
! ------------------------

  LOGICAL        , ALLOCATABLE ::  &
    lneedt       (:) ,& ! indic. temperature must be available at current level
    lneedq       (:) ,& ! indicates humidity must be available at current level
    lneedv       (:) ,& ! indicates   wind   must be available at current level
    llevsfc      (:)    ! indicates this is a surface level (in (ICOS) tower)

  INTEGER        , ALLOCATABLE ::  &
    mzmlbd     (:,:) ,& ! observation body (flags) of report as read / processed
    msupbd     (:,:) ,& ! observation body (flags) of superobbed report
    ksurfob      (:) ,& ! surface report indicator and report index
    kproclv      (:) ,& ! processing of level: -9: discarded; -1,-2: below surf;
                        !   1: level given by height; 2: normal pressure level
                        !  -7: p-z gross error; -8: p-z gross error < 200hPa
    isrtlvp      (:) ,& ! sorted list of indices of good vertical levels
    izdt         (:) ,& ! time since launch
    izlv         (:) ,& ! level identity (vertical sounding significance)
    iroll        (:) ,& ! roll angle quality flag  (aircraft, BUFR Table 002064)
    iqci         (:) ,& ! quality index            (profiler, BUFR Table 033002,
                        !                                    with use of 025034)
    isls         (:) ,& ! (ICOS) tower: first level with same height
    isle         (:) ,& ! (ICOS) tower: last  level with same height
    kblk       (:,:) ,& ! local blacklist (0: ok , 1: blacklisted, for each pro-
                        !                  cessed level and variable separately)
    nflgx        (:) ,& ! total flag for an observation of a variable
    ioberr       (:)    ! active status                  (processed levels only)

  REAL (KIND=wp) , ALLOCATABLE ::  &
    zmlbdy     (:,:) ,& ! observation body of report as read / processed
    supbdy     (:,:) ,& ! observation body of superobbed report
    zpp          (:) ,& ! pressure
    zppl         (:) ,& ! pressure as vertical level indep. from model state
                        !   (if level by height, then zppl from std. atmosphere)
    zdlat        (:) ,& ! latitude  displacement
    zdlon        (:) ,& ! longitude displacement
    ztt          (:) ,& ! temperature
    ztd          (:) ,& ! dew-point temperature
    zrh          (:) ,& ! relative humidity
    zqx          (:) ,& ! mixing ratio
    zff          (:) ,& ! wind speed
    zzz          (:) ,& ! height
    zdd          (:) ,& ! wind direction
    zww          (:) ,& ! vertical velocity
    zsinor       (:) ,& ! signal to noise ratio
    zstdff       (:) ,& ! standard deviation wind speed
    zazi         (:) ,& ! azimuth (ICOS tower)
!   ztv          (:) ,& ! virtual temperature            (processed levels only)
    zttl         (:) ,& ! temperature                    (processed levels only)
    ztdl         (:) ,& ! dew-point temperature          (processed levels only)
    zqxl         (:) ,& ! mixing ratio                   (processed levels only)
    zrhw         (:) ,& ! relative humidity over water   (processed levels only)
    zrhc         (:) ,& ! model compatible rel. humidity (processed levels only)
    zrhw2        (:) ,& ! rel. humidity over water from dewpoint (proc. levels)
    zrhc2        (:) ,& ! compatible rel. humidity from dewpoint (proc. levels)
    zrhw3        (:) ,& ! rel. hum. over water from mixing ratio (proc. levels)
    zrhc3        (:) ,& ! compatible rel. hum. from mixing ratio (proc. levels)
    zqvw         (:) ,& ! specific humidity over water (from 'zqxl', proc. lev.)
    zqv          (:) ,& ! model compatible spec. humid (from 'zqxl', proc. lev.)
    zuu        (:,:) ,& ! zonal      wind                (processed levels only)
    zvv        (:,:) ,& ! meridional wind                (processed levels only)
    zrlat      (:,:) ,& ! latitude  of (processed) observation level
    zrlon      (:,:) ,& ! longitude of (processed) observation level
    zoberr       (:) ,& ! observation error              (processed levels only)
    zstd         (:) ,& ! height from pressure with standard atmosphere
    p_int        (:)    ! pressure levels for superobbing of current report
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_store_multilev
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_store_multilev'
  c3600r = c1 / c3600
  lbogus = (ABS( fbogus ) > epsy)

! lrs_temp  =     (kcdftyp == ncdf_temp )    .OR. (kcdftyp == ncdf_tempship)  &
!            .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)
! lrs_pilot =     (kcdftyp == ncdf_pilot)    .OR. (kcdftyp == ncdf_pilot_p )
  ltempasc  =     (kcdftyp == ncdf_temp)     .OR. (kcdftyp == ncdf_tempship)  &
             .OR. (kcdftyp == ncdf_temphirs )
  lamdarml  =     (kcdftyp == ncdf_amdar_ml ) .OR. (kcdftyp == ncdf_amdar_vp)
  lprofw    =     (kcdftyp == ncdf_wprof    ) .OR. (kcdftyp == ncdf_radar_vad) &
             .OR. (kcdftyp == ncdf_wlidar_wp)
  lproft    =     (kcdftyp == ncdf_rass     )
  ltower    =     (kcdftyp == ncdf_tower    ) .OR. (kcdftyp == ncdf_tower_icos)
  ixrhs      = 4
  IF ((ltempasc) .OR. (kcdftyp == ncdf_tempdesc)                               &
                 .OR. (kcdftyp == ncdf_tempdrop))  ixrhs = 3
  IF (lamdarml)  ixrhs = 5

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

! total number of scalar elements
! -------------------------------
  niscal = noffscal(1,2)
  nrscal = noffscal(2,2)

!-------------------------------------------------------------------------------
! Section 1: Store report header
!-------------------------------------------------------------------------------

  ALLOCATE ( ksurfob (nrepl+1) , STAT=istat )
  ksurfob (:)  =  0

! store header part, which is common to all observation types, in ODR
! -------------------------------------------------------------------

  CALL obs_cdf_store_comhead ( 1, ntotml, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, nodrnew(1), nexceml, ksurfob )
! ==========================

! store the same header part for single-level surface reports derived from TEMP
!   (here, a surface report header is created for all obs types,
!    but is prepared for cancellation in Section 5 if (ilevsfc == 0); this
!    always applies if obs/code type is not TEMP, PILOT balloon, or ICOS tower)
! -----------------------------------------------------------------------------

  CALL obs_cdf_store_comhead ( 4, ntotsg, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, nodrnew(2), nexcesg, ksurfob )
! ==========================

! update number of reports which can be stored in the ODR now
! (the following line is equivalent to:  nrepml = nodrnew(1))
  nrepml = nrepl - nexceml

! update counters for insufficient ODR size, for caution messages
  nexceed (1) = nexceed(1) + nexceml
  nexceed (2) = nexceed(2) + nexcesg

! store model background info which is available (only) on local sub-domain
! -------------------------------------------------------------------------
! - standard deviation of sub-grid scale model orography (for 10-m wind)
  DO irpl = 1 , nrepml
    nmlob = ntotml + irpl
    omlhed (nmlob,nhssos) = r_sso_sd(momlhd(nmlob,nhio),momlhd(nmlob,nhjo))
  ENDDO

! store header part, which is specific to NetCDF
! observation input file type, in multi-level ODR
! -----------------------------------------------

! radiosondes (TEMPs and PILOTs)
  IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)          &
      .OR. (kcdftyp == ncdf_temphirs)                                          &
      .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)          &
      .OR. (kcdftyp == ncdf_pilot)    .OR. (kcdftyp == ncdf_pilot_p )) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                        &
        momlhd (nmlob,nhrtyp) = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 3) /= imiss)                        &
        momlhd (nmlob,nhtrac) = ibuf (iioffs(irpl)+noffscal(1,1)+ 3)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 4) /= imiss)                        &
        momlhd (nmlob,nhna4 ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 4)
!   time significance is not used / stored here currently
!     IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 5) /= imiss)                        &
!       mtimsig (irpl)        = ibuf (iioffs(irpl)+noffscal(1,1)+ 5)
    ENDDO
    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)                                        &
        .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
      DO irpl = 1 , nrepml
        nmlob = ntotml + irpl
        IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 6) /= imiss)                      &
          momlhd (nmlob,nhrad ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 6)
      ENDDO
    ENDIF
  ENDIF
! aircrafts (phase of flight)
  IF (lamdarml) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      mphase = nibits(nvapoc)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                        &
        mphase                = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
      CALL MVBITS( mphase , 0 , nvapoc , momlhd(nmlob,nhschr) , nvapbp )
    ENDDO
  ENDIF
! ground-based remote-sensing profilers (wind profiler, radar VAD, RASS)
  IF (     (lprofw) .OR. (lproft)                                              &
      .OR. (ltower)) THEN
    DO irpl = 1 , nrepml
      nmlob = ntotml + irpl
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        momlhd (nmlob,nhnlev) = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      IF ((lprofw) .OR. (lproft)) THEN
        IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                      &
          momlhd (nmlob,nhdt  ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
        IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 3) /= imiss)                      &
          momlhd (nmlob,nhwce ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 3)
        IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 4) /= imiss)                      &
          momlhd (nmlob,nhna4 ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 4)
!   time significance is not used / stored here currently
!       IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 5) /= imiss)                      &
!         mtimsig (irpl)        = ibuf (iioffs(irpl)+noffscal(1,1)+ 5)
      ENDIF
    ENDDO
  ENDIF

! get list of blacklisted vertical intervals for each report
! ----------------------------------------------------------

  ilstid_chk  =  MIN( ilstid, ilstid_blk )
  DO irpl = 1 , nrepml
    nmlob = ntotml + irpl
    jobtyp (irpl)  =  momlhd(nmlob,nhobtp)
    jcdtyp (irpl)  =  momlhd(nmlob,nhcode)
    zobhr  (irpl)  =  omlhed(nmlob,nhtime)
    ystidl (irpl)  =  ' '
    ystidl (irpl)  =  yomlhd(nmlob) (1:ilstid_chk)
  ENDDO
  IF (nrepml >= 1) THEN
    ALLOCATE ( blk_loc (nrepml) , STAT=istat )

    CALL obs_cdf_blacklist_local ( nrepml , jobtyp , ilstid_chk , ystidl )
!   ============================

! determine which reports are missing on the whitelists
! -----------------------------------------------------

    CALL obs_cdf_whitelist_local ( nrepml , jobtyp , jcdtyp , zobhr            &
                                 , ilstid_chk , ystidl , lblkw )
!   ============================

  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Report body: Put the observations from the observation type
!            dependent buffers into standardised arrays
!-------------------------------------------------------------------------------

! store cloud part, which is specific to TEMP, in single-level ODR
! ----------------------------------------------------------------

! IF (     (kcdftyp == ncdf_temp ) .OR. (kcdftyp == ncdf_tempship)             &
!     .OR. (kcdftyp == ncdf_temphirs)) THEN

!   DO irpl = 1 , nrepml
!     IF (ksurfob(irpl) >= 1) THEN
!       nsgob = ksurfob(irpl)
! restrict to standard observing rules for base of lowest cloud and cloud types
! (yet to be done if required: low / middle / high cloud type)
!       IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 7) == 0) THEN
!         zclz = rbuf(iroffs(irpl)+noffscal(2,1)+ 1)
!         nccl = ibuf(iioffs(irpl)+noffscal(1,1)+ 8)
!         IF ((nccl /= imiss) .AND. (nccl >= 0) .AND. (nccl <= 9)) THEN
!           nhkey = nibits(nvhoc)
!           IF (ABS(zclz) < rmisschk) nhkey = 9
!           IF (zclz < 2500._wp-epsy) nhkey = 8
!           IF (zclz < 2000._wp-epsy) nhkey = 7
!           IF (zclz < 1500._wp-epsy) nhkey = 6
!           IF (zclz < 1000._wp-epsy) nhkey = 5
!           IF (zclz <  600._wp-epsy) nhkey = 4
!           IF (zclz <  300._wp-epsy) nhkey = 3
!           IF (zclz <  200._wp-epsy) nhkey = 2
!           IF (zclz <  100._wp-epsy) nhkey = 1
!           IF (zclz <   50._wp-epsy) nhkey = 0
!           IF (nhkey <= 7) THEN
!             osgbdy (nsgob,nbscl ) = REAL( MIN( nccl, 8 ) , wp )
!           ELSEIF (nhkey <= 9) THEN
!             osgbdy (nsgob,nbscl ) = c0
!           ENDIF
!           nclwgr = imdi
!           CALL MVBITS( nhkey , 0 , nvhoc  , nclwgr , nvhbp  )
!           CALL MVBITS( nccl  , 0 , nvnhoc , nclwgr , nvnhbp )
! write combined cloud and weather group to ODR
!           mosgbd (nsgob,nbscwg) = nclwgr
!         ENDIF
!       ENDIF
!     ENDIF
!   ENDDO
! ENDIF

! store cloud part, which is specific to TEMP, in single-level ODR
! ----------------------------------------------------------------

  IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)          &
      .OR. (kcdftyp == ncdf_temphirs)) THEN

    DO irpl = 1 , nrepml
      IF (ksurfob(irpl) >= 1) THEN
        nsgob = ksurfob(irpl)
        ncsig  = ibuf(iioffs(irpl)+noffscal(1,1)+ 7)
        nccl   = ibuf(iioffs(irpl)+noffscal(1,1)+ 8)
        zclz   = rbuf(iroffs(irpl)+noffscal(2,1)+ 1)
! restrict to vertical significance = standard observing rule or clear sky or
!                                     low, middle, or high cloud
!         and cloud cover between 0 and 9
        IF (      (nccl /= imiss) .AND. (nccl >= 0) .AND. (nccl <= 9)          &
            .AND. (     (ncsig == 0) .OR. (ncsig == 7) .OR. (ncsig == 8)       &
                   .OR. (ncsig == 9) .OR. (ncsig == 62))) THEN
          nclwgn = imdi
          nclwgr = imdi
          CALL MVBITS( ncsig , 0 , nxsgoc , nclwgn , nxsgbp )
! cloud base height: convert 'm' into key code
!   (lower limits of bins: 0, 50, 100, 200, 300, 600, 1000, 1500, 2000, 2500 m)
          IF (     (ABS(zclz) > rmisschk)                                      &
              .OR. (zclz < c0) .OR. (zclz > 16380._wp)) THEN
            zclz = rmdi
            nhkey = nibits(nvhoc)
          ELSE
                                      nhkey = 9
            IF (zclz < 2500._wp-epsy) nhkey = 8
            IF (zclz < 2000._wp-epsy) nhkey = 7
            IF (zclz < 1500._wp-epsy) nhkey = 6
            IF (zclz < 1000._wp-epsy) nhkey = 5
            IF (zclz <  600._wp-epsy) nhkey = 4
            IF (zclz <  300._wp-epsy) nhkey = 3
            IF (zclz <  200._wp-epsy) nhkey = 2
            IF (zclz <  100._wp-epsy) nhkey = 1
            IF (zclz <   50._wp-epsy) nhkey = 0
            osgbdy (nsgob,nbscbs) = zclz
          ENDIF
          CALL MVBITS( nhkey , 0 , nvhoc , nclwgr , nvhbp )
! low / middle cloud cover
          IF ((nccl == imdi) .OR. (nccl > nibits(nxcloc)) .OR. (nccl < 0))     &
            nccl  =  nibits(nxcloc)
          CALL MVBITS( nccl , 0 , nxcloc , nclwgn , nxclbp )
          CALL MVBITS( nccl , 0 , nvnhoc , nclwgr , nvnhbp )
! derive low and middle cloud cover
          nclcl = nccl
          nclcm = nccl
          nclqf = 0
!     if clear sky (vert. signif. = value not applicable)
          IF ((ncsig == 62) .AND. (nccl == nibits(nvnoc))) THEN
            nclcl = 0
            nclcm = 0
!     if 'MNH'=9 (sky invisible)
          ELSEIF (nclcl == 9) THEN
            nclcl = 8
            nclqf = 1
            nclcm = nibits(nxcloc)
          ELSE
            IF (nhkey <= 7) THEN           ! cloud base height below 2000 m
              nclcm = nibits(nxcloc)
            ELSEIF (nhkey <= 9) THEN       ! cloud base height above 2000 m
              nclcl = 0
            ENDIF
                                           ! vertical significance :
            IF (ncsig == 7) THEN           !   low    cloud  (WMO Table 008002)
              nclcm = nibits(nxcloc)
            ELSEIF (ncsig == 8) THEN       !   middle cloud
              nclcl = 0
            ELSEIF (ncsig == 9) THEN       !   high   cloud
              nclcl = 0
!     (inconsistency between vertical significance and 'MNH')
              IF (nclcm >= 1) nclcm = nibits(nxcloc)
            ENDIF
          ENDIF
!   if cloud base height and vertical significance not consistent:
!   put missing value for low and middle cloud cover
          IF (     ((ncsig == 9) .AND. (nhkey <= 8))                           &
              .OR. ((ncsig == 8) .AND. (nhkey <= 7))                           &
              .OR. ((ncsig <= 7) .AND. (nhkey == 9))) THEN
            nclcl = nibits(nxcloc)
            nclcm = nibits(nxcloc)
          ENDIF
!   broken, scattered, or few cloud
          IF ((nclcl >= 11) .AND. (nclcl <= 13))  nclqf = 1
          IF (nclcl == 12)  nclcl = 6
          IF (nclcl == 11)  nclcl = 2
          IF (nclcl == 13)  nclcl = 1
          IF (nclcm == 12)  nclcm = 6
          IF (nclcm == 11)  nclcm = 2
          IF (nclcm == 13)  nclcm = 1
!   obscured or undefined
          IF (nclcl >   8)  nclcl = nibits(nxcloc)
          IF (nclcm >   8)  nclcm = nibits(nxcloc)
!     (at this point 0 <= nclcl, nclcm <= 8 , or == nibits(nxcloc))
!   write low and middle cloud cover to ODR,
!   and include low cloud flag in combined flag word
          IF (nclcl < nibits(nxcloc)) osgbdy (nsgob,nbscl) = REAL(nclcl, wp)
          IF (nclcm < nibits(nxcloc)) osgbdy (nsgob,nbscm) = REAL(nclcm, wp)
          CALL MVBITS( nclqf , 0 , nvcqoc , mosgbd(nsgob,nbswwe) , nvcqbp )
! low cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+ 9)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , nctlbp )
          IF ((ncc >= 30) .AND. (ncc <= 39)) THEN
            ncc = ncc - 30
          ELSEIF ((ncc == 59) .OR. (ncc == 62)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvcloc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvcloc , nclwgr , nvclbp )
! middle cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+10)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , nctmbp )
          IF ((ncc >= 20) .AND. (ncc <= 29)) THEN
            ncc = ncc - 20
          ELSEIF ((ncc == 59) .OR. (ncc == 61)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvcmoc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvcmoc , nclwgr , nvcmbp )
! high cloud type
          ncc = ibuf(iioffs(irpl)+noffscal(1,1)+11)
          CALL MVBITS( ncc , 0 , nxctoc , nclwgn , ncthbp )
          IF ((ncc >= 10) .AND. (ncc <= 19)) THEN
            ncc = ncc - 10
          ELSEIF ((ncc == 59) .OR. (ncc == 60)) THEN
            ncc = 10
          ELSE
            ncc = nibits(nvchoc)
          ENDIF
          CALL MVBITS( ncc , 0 , nvchoc , nclwgr , nvchbp )
! write combined cloud and weather group to ODR
          mosgbd (nsgob,nbsclg) = nclwgn
          mosgbd (nsgob,nbscwg) = nclwgr
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! store specific surface-level part for ICOS tower reports
! --------------------------------------------------------

  IF (ltower) THEN
    DO irpl = 1 , nrepml
      IF (ksurfob(irpl) >= 1) THEN
        nsgob = ksurfob(irpl)
        nmlob = ntotml + irpl
        !   pressure at certain height
        IF (ABS( rbuf(iroffs(irpl)+noffscal(2,1)+ 1) ) < rmisschk)             & 
          osgbdy (nsgob,nbsz  ) = rbuf(iroffs(irpl)+noffscal(2,1)+ 1)
        IF (ABS( rbuf(iroffs(irpl)+noffscal(2,1)+ 2) ) < rmisschk) THEN
          osgbdy (nsgob,nbsp  ) = rbuf(iroffs(irpl)+noffscal(2,1)+ 2)
          osgbdy (nsgob,nbsplv) = rbuf(iroffs(irpl)+noffscal(2,1)+ 2)
        ENDIF
        !   global radiation and its accumulation period
        IF (ABS( rbuf(iroffs(irpl)+noffscal(2,1)+ 3) ) < rmisschk) THEN
          osgbdy (nsgob,nbsrad) = rbuf(iroffs(irpl)+noffscal(2,1)+ 3)
          IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 2) /= imiss)                    &
            mosghd (nsgob,nhdt  ) = ibuf (iioffs(irpl)+noffscal(1,1)+ 2)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! allocate standardised arrays
! ----------------------------

! get maximum number of vertical levels
  maxlev = 0
  DO irpl = 1 , nrepml
    nmlob = ntotml + irpl
    maxlev = MAX( maxlev, momlhd(nmlob,nhnlev) )
  ENDDO

  ALLOCATE ( zpp     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zppl    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdlat   (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdlon   (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztt     (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztd     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zqx     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zff     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zzz     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zdd     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zww     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zsinor  (maxlev+1)   , STAT=istat )
  ALLOCATE ( zstdff  (maxlev+1)   , STAT=istat )
  ALLOCATE ( zazi    (maxlev+1)   , STAT=istat )
  ALLOCATE ( izdt    (maxlev+1)   , STAT=istat )
  ALLOCATE ( izlv    (maxlev+1)   , STAT=istat )
  ALLOCATE ( iroll   (maxlev+1)   , STAT=istat )
  ALLOCATE ( iqci    (maxlev+1)   , STAT=istat )
  ALLOCATE ( isls    (maxlev+1)   , STAT=istat )
  ALLOCATE ( isle    (maxlev+1)   , STAT=istat )
  ALLOCATE ( llevsfc (maxlev+1)   , STAT=istat )
  ALLOCATE ( kproclv (maxlev+1)   , STAT=istat )
  ALLOCATE ( isrtlvp (maxlev+1)   , STAT=istat )

! the following variables are used in section 4 only, and could alternatively be
! re-allocated there for each report with length 'nlvp' (within the report loop)
  ALLOCATE ( kblk    (maxlev+1,4) , STAT=istat )
  ALLOCATE ( nflgx   (maxlev+1)   , STAT=istat )
  ALLOCATE ( ioberr  (maxlev+1)   , STAT=istat )
! ALLOCATE ( ztv     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zttl    (maxlev+1)   , STAT=istat )
  ALLOCATE ( ztdl    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zqxl    (maxlev+1)   , STAT=istat )
! ALLOCATE ( ztqx    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zrhw    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zrhc    (maxlev+1)   , STAT=istat )
  ALLOCATE ( zrh     (maxlev+1)   , STAT=istat )
  ALLOCATE ( zuu     (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zvv     (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zrlat   (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zrlon   (maxlev+1,1) , STAT=istat )
  ALLOCATE ( zoberr  (maxlev+1)   , STAT=istat )
  ALLOCATE ( zstd    (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedt  (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedq  (maxlev+1)   , STAT=istat )
  ALLOCATE ( lneedv  (maxlev+1)   , STAT=istat )
! IF (my_cart_id == 0) PRINT *,'rmisschk', rmisschk, imiss, nrepml

  ! first fill temporary array 'zmlbdy', 'mzmlbd', which can have thousands of
  ! levels, then apply superobbing for high-resolution profiles (with more than
  ! 'maxrtv' levels in arrays 'zmlbdy', 'mzmlbd'), and finally copy these
  ! profiles into 'omlbdy', 'momlbd' (which have <= maxrtv levels)
  ALLOCATE ( zmlbdy  (maxlev+1,mxrbdy)   , STAT=istat )
  ALLOCATE ( mzmlbd  (maxlev+1,mxrbdf)   , STAT=istat )

  nrepadd = 0

! ----------------------
! huge loop over reports
! ----------------------

  DO irpl = 1 , nrepml
! ~~~~~~~~~~~~~~~~~~~~

    nmlob = ntotml + irpl
    nlev   = momlhd(nmlob,nhnlev)
    kobtyp = momlhd(nmlob,nhobtp)
    kcdtyp = momlhd(nmlob,nhcode)
    iob    = momlhd(nmlob,nhio)
    job    = momlhd(nmlob,nhjo)
    icma   = i_cma ( kobtyp , kcdtyp )
!            =====
    iactr  = 0
    IF (momlhd(nmlob,nhpass) == 0)  iactr  = 1
    zstalt = omlhed(nmlob,nhalt)
    zmlbdy (:,:) = rmdi
    ! additional safety even though all elements are set explicitly below
    mzmlbd (:,:) = 0

! store NetCDF file type dependent observations in standarised temporary arrays
! -----------------------------------------------------------------------------
    DO ilev = 1 , nlev
      zpp   (ilev) = rmdi
      zppl  (ilev) = rmdi
      zdlat (ilev) = rmdi
      zdlon (ilev) = rmdi
      ztt   (ilev) = rmdi
      ztd   (ilev) = rmdi
      zqx   (ilev) = rmdi
      zrh   (ilev) = rmdi
      zrhw  (ilev) = rmdi
      zrhc  (ilev) = rmdi
      zff   (ilev) = rmdi
      zzz   (ilev) = rmdi
      zdd   (ilev) = rmdi
      zww   (ilev) = rmdi
      zsinor(ilev) = rmdi
      zstdff(ilev) = rmdi
      zazi  (ilev) = rmdi
      izdt  (ilev) = 0
      izlv  (ilev) = imdi
      iroll (ilev) = imdi
      iqci  (ilev) = imdi
      isls  (ilev) = imdi
      isle  (ilev) = imdi
      llevsfc(ilev) = .FALSE.
    ENDDO
    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)                                        &
        .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)        &
        .OR. (kcdftyp == ncdf_pilot)    .OR. (kcdftyp == ncdf_pilot_p )) THEN
      DO ilev = 1 , nlev
        IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
          zdlat (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
          zdlon (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
          zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+3*nlev+ilev) ) < rmisschk) THEN
          zpp   (ilev) =       rbuf (iroffs(irpl)+nrscal+3*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+4*nlev+ilev) ) < rmisschk) THEN
          zzz   (ilev) =       rbuf (iroffs(irpl)+nrscal+4*nlev+ilev)
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+  nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
          izdt  (ilev) =       ibuf (iioffs(irpl)+niscal+2*nlev+ilev)
        ENDIF
        IF (                  ( ibuf(iioffs(irpl)+niscal+ilev) /= imiss  )     &
            .AND. (.NOT. BTEST( ibuf(iioffs(irpl)+niscal+ilev),ilv_miss ))) THEN
          izlv  (ilev) =       ibuf (iioffs(irpl)+niscal       +ilev)
        ENDIF
      ENDDO
      IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)      &
          .OR. (kcdftyp == ncdf_temphirs)                                      &
          .OR. (kcdftyp == ncdf_tempdrop) .OR. (kcdftyp == ncdf_tempdesc)) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+5*nlev+ilev) ) < rmisschk) THEN
            ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+5*nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+6*nlev+ilev) ) < rmisschk) THEN
            ztd   (ilev) =       rbuf (iroffs(irpl)+nrscal+6*nlev+ilev)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    IF (lamdarml) THEN
      DO ilev = 1 , nlev
        IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
          zdlat (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
          zdlon (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
          zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+3*nlev+ilev) ) < rmisschk) THEN
          ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+3*nlev+ilev)
        ENDIF
        IF (ABS( rbuf(iroffs(irpl)+nrscal+4*nlev+ilev) ) < rmisschk) THEN
          ztd   (ilev) =       rbuf (iroffs(irpl)+nrscal+4*nlev+ilev)
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+  nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
          zzz   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+2*nlev+ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal       +ilev)  /= imiss   ) THEN
          iroll (ilev) =       ibuf (iioffs(irpl)+niscal       +ilev)
        ENDIF
      ENDDO
    ENDIF
    IF ((lprofw) .OR. (lproft) .OR. (ltower)) THEN
      DO ilev = 1 , nlev
        IF (     ibuf(iioffs(irpl)+niscal       +ilev)  /= imiss   ) THEN
          zzz   (ilev) = REAL( ibuf (iioffs(irpl)+niscal       +ilev) , wp )
        ENDIF
        IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   ) THEN
          iqci  (ilev) =       ibuf (iioffs(irpl)+niscal+  nlev+ilev)
        ENDIF
      ENDDO
      IF ((lprofw) .OR. (lproft)) THEN
        DO ilev = 1 , nlev
          IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
            zsinor(ilev) = REAL( ibuf (iioffs(irpl)+niscal+2*nlev+ilev) , wp )
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
            zww   (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
          ENDIF
        ENDDO
      ENDIF
      IF (lprofw) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
            zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
            zstdff(ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
          ENDIF
          IF (     ibuf(iioffs(irpl)+niscal+3*nlev+ilev)  /= imiss   ) THEN
            zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+3*nlev+ilev) ,wp)
          ENDIF
        ENDDO
      ELSEIF (lproft) THEN
        DO ilev = 1 , nlev
          IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
            ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
          ENDIF
        ENDDO
      ELSEIF (ltower) THEN
        DO ilev = 1 , nlev
          IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   ) THEN
            zrh   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+2*nlev+ilev) ,wp)   &
                           * 0.01_wp
          ENDIF
          IF (     ibuf(iioffs(irpl)+niscal+3*nlev+ilev)  /= imiss   ) THEN
            zdd   (ilev) = REAL( ibuf (iioffs(irpl)+niscal+3*nlev+ilev) ,wp)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk) THEN
            ztt   (ilev) =       rbuf (iroffs(irpl)+nrscal       +ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+  nlev+ilev) ) < rmisschk) THEN
            zff   (ilev) =       rbuf (iroffs(irpl)+nrscal+  nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+2*nlev+ilev) ) < rmisschk) THEN
            ztd   (ilev) =       rbuf (iroffs(irpl)+nrscal+2*nlev+ilev)
          ENDIF
          IF (ABS( rbuf(iroffs(irpl)+nrscal+3*nlev+ilev) ) < rmisschk) THEN
            zazi  (ilev) =       rbuf (iroffs(irpl)+nrscal+3*nlev+ilev)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Decide which vertical levels shall be discarded, and
!            compile a sorted list of levels which shall be processed further
!-------------------------------------------------------------------------------

! decide whether levels are reported in pressure or height or mixed
! -----------------------------------------------------------------
!   --> indlev = -1 : all levels with neither pressure nor height reported
!                 0 : only levels with both pressure and height reported
!                       (after the first assignment, this is set to '2')
!                 1 : height reported at all levels    ;  pressure \ at some or
!                 2 : pressure reported at all levels  ;    height / no levels
!                 3 : mixed sequence of height and pressure levels: this is
!                       explicitly excluded in a second step by removing levels
!               (0-3: levels with neither pressure not height may also occur)

    indlev = -1
    DO ilev = 1 , nlev
      IF     ((zzz(ilev) <= rmdich) .AND. (zpp(ilev) > rmdich)) THEN
        IF (indlev <= 0)  indlev = 2
        IF (indlev == 1)  indlev = 3
      ELSEIF ((zpp(ilev) <= rmdich) .AND. (zzz(ilev) > rmdich)) THEN
        IF (indlev <= 0)  indlev = 1
        IF (indlev == 2)  indlev = 3
      ELSEIF ((zpp(ilev) >  rmdich) .AND. (zzz(ilev) > rmdich)) THEN
        IF (indlev == -1) indlev = 0
      ENDIF
    ENDDO
    ! from this point, (indlev == 0) can be treated like (indlev == 2)
    IF (indlev == 0)      indlev = 2

    ! for TEMP, always require pressure levels (assume that if pressure levels
    !           are missing, then something is wrong with the levels (or the
    !           whole report); this can happen e.g. with descending TEMPs)
    IF (kobtyp == ntemp)  indlev = 2

    ! in principle, (indlev == 3) should never occur (at least if using obs
    ! from the DWD data base);
    ! however for safety, it should be excluded explicitly because ordering
    ! levels in the vertical (and some other processing steps) independently
    ! from the model state would imply the need to compare reported pressure
    ! for some levels with pressure derived from height by standard atmosphere
    ! for other levels, and this does not make sense
    IF (indlev == 3) THEN
      ! check if there are more pressure or height levels below 300 hPa
      ! (= 9164 m in std. atmos.), and remove the levels of the more rare type
      nlvp = 0
      nlvz = 0
      DO ilev = 1 , nlev
        IF ((zpp(ilev) > rmdich) .AND. (zpp(ilev) >= 30000._wp)) nlvp = nlvp + 1
        IF ((zzz(ilev) > rmdich) .AND. (zzz(ilev) <   9164._wp)) nlvz = nlvz + 1
      ENDDO
      indlev = 2
      IF (nlvz > nlvp) indlev = 1
    ENDIF
    ! from this point, we need to discriminate only between 'indlev' = 1, 2

    DO ilev = 1 , nlev
      IF ((indlev == 2) .AND. (zpp(ilev) <= rmdich))  zzz (ilev) = rmdi
      IF ((indlev == 1) .AND. (zzz(ilev) <= rmdich))  zpp (ilev) = rmdi
    ENDDO

! check vertical coordinate and level (pressure or height)
! --------------------------------------------------------

    ! set 'zppl' (pressure level independent from model state) to 'zpp'
    !   (for (indlev == 3), 'zppl' will remain undefined here for some levels)
    ! except for (indlev == 1) where 'zppl' must be derived from 'zzz' for all
    !   levels (even if some levels also have 'zpp' reported)
    IF (indlev /= 1) THEN
      DO ilev = 1 , nlev
        zppl (ilev) = zpp(ilev)
      ENDDO
    ENDIF

    DO ilev = 1 , nlev
      kproclv (ilev) =  2

      ! if height is given instead of pressure as a vertical coordinate:
!     IF (      (zpp(ilev) <= rmdich)                                          &
!         .AND. (zzz(ilev) >  rmdich) .AND. (kobtyp /= ntemp)) THEN
      IF (      ( zzz(ilev) >  rmdich)                                         &
          .AND. ((zpp(ilev) <= rmdich) .OR. (indlev == 1)                      &
                                       .OR. (indlev == 4))) THEN
        kproclv (ilev) = 1
        ! tower data: 'zzz' contains height of sensor above ground;
        ! ----------  --> allow only levels up to 999 m and
        !                 convert into altitude by adding station altitude;
        !             later on in the obs operator, there are 2 possibilities:
        !              a) apply the standard obs operator, such that the model
        !                 equivalents are computed and compared with the obs
        !                 at equal height / pressure level (but not at equal
        !                 'height of sensor above ground')
        !              b) first shift the model column vertically according to
        !                 the difference between station height and model
        !                 orography and apply a height correction to T, qv, and
        !                 then apply the standard obs operator, such that the
        !                 model equivalents are computed and compared with the
        !                 obs at equal 'height of sensor above ground'.
        !                 (An alternative implementation for option b) by adding
        !                  the model orography to 'zzz' here and later on simply
        !                  applying the standard obs operator has been removed
        !                  (due to issues in implementing a height correction to
        !                   the obs and inconsistencies in the feedback files).)
        !             obs levels with height of sensor height above ground
        !!!            -  < 1.5 m are discarded !!!
        !              -  between 1.5 m and 2.5 m are treated as 'surface level'
        IF (ltower) THEN
          IF ((zzz(ilev) < 0._wp) .OR. (zzz(ilev) > 999._wp)) kproclv(ilev) = -9
          IF                 (zzz(ilev) <= 1.49_wp)           kproclv(ilev) = -9
          llevsfc (ilev)  = ((zzz(ilev) >  1.49_wp) .AND. (zzz(ilev) < 2.51_wp))
          zzz (ilev)  =  zzz(ilev) + zstalt
!         zzz (ilev)  =  zzz(ilev) + omlhed(nmlob,nhsurf)
        ENDIF

        !   derive 'zppl' (pressure as vertical level independ.from model state)
        !   from standard atmosphere

        CALL std_atmosphere ( 'z', zzz(ilev), 'p', zppl(ilev), r_g, r_d, irm )
!       ===================
        IF (irm /= 0)  zppl (ilev) = rmdi

        !   AMDAR: zzz is flight level (which was derived from measured pressure
        !          using std. atmos.), therefore set obs pressure 'zpp' = zppl
        !          and set temporally: indlev == 4
!       IF ((lamdarml) .AND. (zpp(ilev) <= rmdich)) THEN
        IF (lamdarml) THEN
          zpp (ilev) = zppl(ilev)
          indlev = 4
        !   others: use model fields to determine pressure at given height level
        ELSEIF ((zppl(ilev) > rmdich) .AND. (zpp(ilev) <= rmdich)) THEN
          zpp (ilev) = f_z2p ( zzz(ilev), ke, r_hhl(iob,job,:), r_p(iob,job,:) &
                             , r_t_ll(iob,job), r_g, r_d, rmdi, ltower )
!                      =====
          ! note: this value of 'zpp' depends on the model state; however
          !       whether 'zpp' has missing value or not depends on model height
          !       rather than pressure and hence is not dependent on model state
          ! note: for towers, pressure is extrapolated beyond model orography
        ENDIF
        !   discard z-levels below model orography (except towers) or above top
        !   model level   (this if clause is not dependent on model state)
        IF ((zpp(ilev) <= rmdich) .OR. (zppl(ilev) <= rmdich)) THEN
          kproclv (ilev) = -9
          neventd (nelnop,icma) = neventd(nelnop,icma) + iactr
        ENDIF

      ELSEIF (zpp(ilev) <= rmdich) THEN
        !   vertical coordinate does not exist
        kproclv (ilev) = -9
        neventd (nelmis,icma) = neventd(nelmis,icma) + iactr
      ENDIF
      ! Note: from this point,
      !       - 'zppl' is independent from the model state
      !                         (i.e. from the ensemble member),
      !       - for all levels within this same report,
      !         'zppl' is - either the reported pressure   (indlev == 2)
      !                   - or derived from flight level   (indlev == 4)
      !                   - or derived from reported height (altitude)
      !                     using the standard atmosphere  (indlev == 1),
      !       - and (zppl > rmdich) only if (zpp > rmdich) for (kproclv > -9)

      ! discard any levels with reported (!) pressure above 'rpplim'
      !               or height above the mean of top main and half model levels
      !               or height more than 40 m below station height
      htop  =  0.75_wp* r_hhl(iob,job,1)  +  0.25_wp* r_hhl(iob,job,2)
      IF (    (     (zpp(ilev) > rmdich) .AND. (kproclv(ilev) == 2)            &
              .AND. (     (zpp(ilev) < rpplim)                                 &
                     .OR. ((zpp(ilev) < 3001._wp) .AND. (htop < 22001._wp))    &
                     .OR. ((zpp(ilev) < 2001._wp) .AND. (htop < 24001._wp))    &
                     .OR. ((zpp(ilev) < 1001._wp) .AND. (htop < 28001._wp))))  &
         .OR. (      (zzz(ilev) > rmdich)                                      &
              .AND. (     ((zzz(ilev)+40._wp < zstalt) .AND.(zstalt > rmdich)) &
                     .OR. (zzz(ilev) > htop)))) THEN
        kproclv (ilev) = -9
        neventd (nelext,icma) = neventd(nelext,icma) + iactr

      ! discard non-mandatory reported (!) pressure levels except for levels
      !   - below 'pminsigt' (i.e. p-obs >= pminsigt) which contain T obs, or
      !   - below 'pminsigv' (i.e. p-obs >= pminsigv) which contain wind obs
      !   (for high-resolution BUFR and for descending sondes, assume high
      !    resolution and apply superobbing without discarding any levels)
      ELSEIF (    (zpp(ilev) > rmdich) .AND. (     (zpp(ilev) < pminsigv-c05)  &
                                              .OR. (zpp(ilev) < pminsigt-c05)) &
             .AND.(kproclv(ilev) == 2)                                         &
             .AND.(   (kcdftyp == ncdf_temp    ).OR.(kcdftyp == ncdf_tempship) &
!                 .OR.(kcdftyp == ncdf_temphirs).OR.(kcdftyp == ncdf_tempdesc) &
                  .OR.(kcdftyp == ncdf_tempdrop)                               &
                  .OR.(kcdftyp == ncdf_pilot)   .OR.(kcdftyp == ncdf_pilot_p)) &
             ) THEN
        nzpp = NINT( zpp(ilev) )
        IF (      (nzpp-85000 /= 0) .AND. (nzpp-70000 /= 0)                    &
            .AND. (nzpp-50000 /= 0) .AND. (nzpp-40000 /= 0)                    &
            .AND. (nzpp-30000 /= 0) .AND. (nzpp-25000 /= 0)                    &
            .AND. (nzpp-20000 /= 0) .AND. (nzpp-15000 /= 0)                    &
            .AND. (nzpp-10000 /= 0) .AND. (nzpp- 7000 /= 0)                    &
            .AND. (nzpp- 5000 /= 0) .AND. (nzpp- 3000 /= 0)                    &
            .AND. (nzpp- 2500 /= 0) .AND. (nzpp- 2000 /= 0)                    &
            .AND. (nzpp- 1000 /= 0)) THEN
          IF (      ((ztt(ilev) < rmdich) .OR. (zpp(ilev) < pminsigt-c05))     &
              .AND. ((zff(ilev) < rmdich) .OR. (zpp(ilev) < pminsigv-c05))) THEN
            kproclv (ilev) = -9
            neventd (nelsig,icma) = neventd(nelsig,icma) + iactr
          ENDIF
        ENDIF
      ! do not (!): discard aircraft data below 'model surface pressure + 1hPa'
      ! because this would depend on model state and thus on ensemble member
!     ELSEIF ((lamdarml) .AND. (zpp(ilev) > rmdich)                            &
!                        .AND. (zpp(ilev) > r_ps(iob,job)+100._wp)) THEN
!       kproclv (ilev) = -9
!       neventd (nelsig,icma) = neventd(nelsig,icma) + iactr
      ENDIF
    ENDDO
    !   from now on, treat AMDAR with pressure derived from flight level
    !   as report with reported pressure (AMDAR originally measured pressure)
    IF (indlev == 4)  indlev = 2

! tower profiles: in case of several 'levels' at the same height
! --------------  choose the one with azimuth closest to upwind direction
    IF (ltower) THEN
      !   find indices of first and last 'level' with the same height
      !     (but different azimuth), possibly several times
      nhsl     = 0
      isle (1) = 0
      DO ilev = 2 , nlev
        IF (ABS( zzz(ilev)-zzz(ilev-1) ) < 0.4_wp) THEN
          IF ((nhsl == 0) .OR. (isle(MAX(nhsl,1)) < ilev-1)) THEN
            nhsl        = nhsl + 1
            isls (nhsl) = ilev - 1
          ENDIF
          isle (nhsl) = ilev
        ENDIF
      ENDDO
      !   loop over heights with several levels
      DO ihsl = 1 , nhsl
        !   compute upwind direction as a mean over individual wind directions
        !   (use only wind directions with wind speeds > 1 m/s)
        zuu (1,1) = 0
        zvv (1,1) = 0
        DO ilev = isls(ihsl) , isle(ihsl)
          IF ((zff(ilev) > rmdich) .AND. (zff(ilev) > 1._wp)) THEN
            zuu(1,1) = zuu(1,1) - SIN( zdd(ilev)*r_degrad )
            zvv(1,1) = zvv(1,1) - COS( zdd(ilev)*r_degrad )
          ENDIF
        ENDDO
        CALL uv2df ( zuu(1,1), zvv(1,1), zddmean, zvalue )
      ! ==========
        !   find level index 'ilact' with azimuth closest to upwind direction
        !     ('zazi' and 'zddmean' are (assumed to be) between 0 and 360
        !     (alternatively, 'zddmean' might be replaced by 'zdd(ilev)')
        !     (if no azimuths are reported take the first level as active one)
        zdazimn = 361._wp
        ilact = isls(ihsl)
        DO ilev = isls(ihsl) , isle(ihsl)
          IF (zazi(ilev) > rmdich) THEN
            zdazi  =  MIN( ABS( zazi(ilev) - zddmean           )               &
                         , ABS( zazi(ilev) - zddmean + 360._wp )               &
                         , ABS( zazi(ilev) - zddmean - 360._wp ) )
            IF (zdazi < zdazimn) THEN
              ilact   = ilev
              zdazimn = zdazi
            ENDIF
            zsinor (ilev) = 90.0_wp  -  0.5_wp* zdazi
          ENDIF
        ENDDO
        !   for level indices other than 'ilact'
        !    - decrease the level pressure / increase the height by small, but
        !        distiguishable amounts (dp: picoshift; dz=-0.08*dp)
        !        (for convenience, a fixed value of 0.08 is used
        !                          to convert from pa to m)
        !    - remove all obs except horizontal wind
        !    - flag level in order to set the wind to rejected later on
        DO ilev = isls(ihsl) , isle(ihsl)
          IF (ilev /= ilact) THEN
            istep = ilev - isls(ihsl) + 1
!           IF (ilev > ilact)  istep = istep - 1
            IF (zpp (ilev) > rmdich)  zpp (ilev) = zpp (ilev) - istep *picoshift
            IF (zppl(ilev) > rmdich)  zppl(ilev) = zppl(ilev) - istep *picoshift
            zzz    (ilev)  =  zzz(ilev) + istep *picoshift *0.08_wp
            ztt    (ilev)  =  rmdi
            ztd    (ilev)  =  rmdi
            zrh    (ilev)  =  rmdi
            iqci   (ilev)  =  IBSET( iqci(ilev) , 11 )
            llevsfc(ilev)  =  .FALSE.
          ENDIF
        ENDDO
      ENDDO
    ENDIF

! check for several surface levels    (assume no surface levels for
! --------------------------------     descending sondes and dropsondes)

    ilevsfc = 0
    IF (     (kcdftyp == ncdf_temp    ) .OR. (kcdftyp == ncdf_tempship)        &
        .OR. (kcdftyp == ncdf_temphirs)                                        &
!       .OR. (kcdftyp == ncdf_tempdesc) .OR. (kcdftyp == ncdf_tempdrop)        &
        .OR. (kcdftyp == ncdf_pilot)    .OR. (kcdftyp == ncdf_pilot_p )) THEN
      DO ilev = 1 , nlev
        !   if surface level
        IF (     (BTEST( izlv(ilev),ilv_sfc )) .AND. (izlv(ilev) /= imdi)      &
           .AND. (kproclv(ilev) /= -9)) THEN
          IF (ilevsfc == 0) THEN
            ilevsfc = ilev
          ELSE
            !   several surface levels: take the one
            !                           which is closer to the station altitude
            IF (     (zzz(ilev) <= rmdich) .OR. (zstalt <= rmdich)             &
                .OR. (ABS(zzz(ilev)-zstalt) > ABS(zzz(ilevsfc)-zstalt))) THEN
              kproclv (ilev)    = -9
            ELSE
              kproclv (ilevsfc) = -9
              ilevsfc           = ilev
            ENDIF
            !   for statistics and control message only
            neventd(nelsfc,icma) = neventd(nelsfc,icma) + iactr
            ilen = 2 + istrej
            IF (nacout+ilen <= nmxoln) THEN
              outbuf(nacout+1) = ilen
              outbuf(nacout+2) = nfmt5
              DO icl = 1 , istrej
                outbuf(nacout+2+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
              ENDDO
              nacout  = nacout + ilen
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      !   allow a level be the surface level only if within 100 m of sta height
      IF (ilevsfc >= 1) THEN
        IF (MIN( zzz(ilevsfc),zstalt ) > rmdich) THEN
          IF (ABS( zzz(ilevsfc) -zstalt ) > 99.9_wp)  ilevsfc = 0
        ENDIF
      ENDIF
    ENDIF

    IF (ilevsfc >= 1) THEN
      !   set missing height of surface level equal to station height (TEMP)
      IF ((ltempasc) .AND. (zzz(ilevsfc) < rmdich) .AND. (zstalt > rmdich))    &
        zzz (ilevsfc) = zstalt

    !   discard any (non-surface) levels below or at the reported surface level
    !     (note: to be independent from model state, 'zppl' must be used)
    !     (note: this would not work properly for the excluded case indlev==3)
      DO ilev = 1 , nlev
        IF (      (zppl(ilev) > zppl(ilevsfc)-50._wp) .AND. (ilev /= ilevsfc)  &
            .AND. (kproclv(ilev) /= -9)) THEN
          kproclv (ilev) = -9
          neventd (nelext,icma) = neventd(nelext,icma) + iactr
        ENDIF
      ENDDO
    ENDIF

! get surface level from tower reports
! ------------------------------------
    IF (ltower) THEN
      DO ilev = 1 , nlev
        IF (llevsfc(ilev))  ilevsfc = ilev
      ENDDO
    ENDIF

! compile a list of indices of those vertical levels which will be processed
! --------------------------------------------------------------------------

    nlvp  = 0
    nlvp2 = 0
    DO ilev = 1 , nlev
      IF (kproclv(ilev) > -9) THEN
        nlvp  =  nlvp + 1
        isrtlvp (nlvp) = ilev
        IF (zpp(ilev) > zp200) nlvp2 = nlvp2 + 1
      ENDIF
    ENDDO

! sort this list of indices in the vertical
! (iterative sorting with bubblesort algorithm)
! ---------------------------------------------
! (this replaces 'obs_sort_levels' in 'src_obs_proc_aof.f90')
! (note: either 'zppl' or 'zpp' can be used for the vertical ordering)

    lchange = .TRUE.
    DO WHILE (lchange)
      lchange = .FALSE.
      DO klev = 1 , nlvp-1
        IF (zppl(isrtlvp(klev)) < zppl(isrtlvp(klev+1))) THEN
          itmp             = isrtlvp (klev+1)
          isrtlvp (klev+1) = isrtlvp (klev)
          isrtlvp (klev)   = itmp
          lchange = .TRUE.
        ENDIF
      ENDDO
    ENDDO

! get surface level
    klevsfc = 0
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (ilev == ilevsfc)  klevsfc = klev
!     PRINT *,'press6 ', yomlhd(nmlob), nlvp, klev, ilev, zpp(ilev), ztt(ilev)
    ENDDO
!   PRINT *,'klevsfc ',yomlhd(nmlob),klevsfc,ilevsfc, zpp(MAX(ilevsfc,1))      &
!                                                   , ztt(MAX(ilevsfc,1))

!!!!!!!!!!!!!!!!!!!!!
!!!   10 hPa: 27050 - 32400 m: 31050m - 4000 +1350
!!!  100 hPa: 14400 - 16800 m: 16180m - 1780 + 620
!!!  200 hPa: 10300 - 12700 m: 11785m - 1485 + 915
!!!  250 hPa:  9200 - 11000 m: 10360m - 1160 + 640
!!!  300 hPa:  7800 -  9900 m:  9164m - 1364 + 734
!!!  500 hPa:  4600 -  6100 m:  5574m -  974 + 526
!!!  700 hPa:  2200 -  3400 m:  3012m -  812 + 388
!!!  850 hPa:   700 -  1800 m:  1457m -  757 + 343
!!! 1000 hPa:  -400 -   300 m:   110m -  510 + 190

! flag pressure levels if incompatible with reported height (gross error check)
! -----------------------------------------------------------------------------
!!! here are min/max height values found in some global ICON forecasts:
!!! pressure level    1000   850   700   500   300   250   200   100    10 hPa
!!! std atm: height    110  1457  3012  5574  9164 10360 11785 16180 31050 m
!!! min(z) at p-lev   -400   700  2200  4600  7800  9200 10300 14400 27050 m
!!! max(z) at p-lev    300  1800  3400  6100  9900 11000 12700 16800 32400 m
!!! max neg. distance -510  -757  -812  -974 -1364 -1160 -1485 -1780 -4000 m
!!! max pos. distance  190   343   388   526   734   640   915   620  1350 m
!!! gross error limits set to:
!!! high pressure: zstd +  500m + 0.04*zstd 
!!! low  pressure: zstd - 1000m - 0.05*zstd for p >= 200 hPa
!!!                zstd         - 0.15*zstd for p <  200 hPa
!!! multi-level check (note: active obs layer in KENDA is at >= 200hPa):
!!! - if >20% of levels >= 200hPa flagged, flag all levels;
!!! - if >20% of levels <  200hPa flagged, flag all levels at < 200hPa

    lcheck_zp = .FALSE.
    IF (indlev == 2) THEN
      DO klev = 1 , nlvp
        ilev = isrtlvp(klev)
        IF ((zzz(ilev) > rmdich) .AND. (ilev /= ilevsfc))  lcheck_zp = .TRUE.
      ENDDO
    ENDIF
    IF (lcheck_zp) THEN
      zstd  (ilev) = rmdi
      zzmis        = -9999._wp
      CALL std_atm_p2z_fast ( nlev, nlvp, isrtlvp, zpp, r_g, r_d, zzmis , zstd )
    ! =====================
      ! gross error check for individual levels as in the table above
      DO klev = 1 , nlvp
        ilev = isrtlvp(klev)
        IF ((zzz(ilev) > rmdich) .AND. (zstd(ilev) /= zzmis)) THEN
          IF (zpp(ilev) >= 19999._wp) THEN
            IF (     (zzz(ilev) < zstd(ilev)-1000._wp -0.05_wp*zstd(ilev))     &
                .OR. (zzz(ilev) > zstd(ilev)+ 500._wp +0.04_wp*zstd(ilev)))    &
              kproclv (ilev) = -7
          ELSE
            IF (     (zzz(ilev) < zstd(ilev)          -0.15_wp*zstd(ilev))     &
                .OR. (zzz(ilev) > zstd(ilev)+ 500._wp +0.04_wp*zstd(ilev)))    &
              kproclv (ilev) = -8
          ENDIF
        ENDIF
      ENDDO
      ! if > 20% of all obs >= 200 hPa are flagged, all obs are flagged
      IF (COUNT( kproclv(1:nlev) == -7 ) > 0.2_wp*nlvp2)                       &
        kproclv (1:nlev) = MIN( kproclv(1:nlev), -7 )
      ! if > 20% of all obs <  200 hPa are flagged, all obs < 200hPa are flagged
      IF (COUNT( kproclv(1:nlev) == -8 ) > 0.2_wp*(nlvp-nlvp2)) THEN
        WHERE (zpp(1:nlev) < zp200)                                            &
          kproclv (1:nlev) = MIN( kproclv(1:nlev), -8 )
      ENDIF
      neventd (nelflg,icma) = neventd(nelflg,icma) + COUNT( kproclv == -7 )
      neventd (nelflg,icma) = neventd(nelflg,icma) + COUNT( kproclv == -8 )
!     WHERE (kproclv(1:nlev) <= -7)  kproclv (1:nlev) = -9
    ENDIF

! blacklist
! ---------
! TYPE blacklist_loc
!   INTEGER :: ndim            ! number of blacklisted intervals in report
!   INTEGER :: kvar (maxintv)  ! variable type of blacklisted interval
!                              !   1: wind, 2: geopot, 3: temperat., 4: humidity
!   REAL    :: plow (maxintv)  ! lower boundary (pressure) of blacklist interval
!   REAL    :: pup  (maxintv)  ! upper boundary (pressure) of blacklist interval
! END TYPE blacklist_loc
! TYPE (blacklist_loc) , DIMENSION(:), ALLOCATABLE  :: blk_loc(nrepml)
! ALLOCATE ( blk_loc (nrepml) , STAT=istat )

!-------------------------------------------------------------------------------
! Section 4: Determine and store the elements of the multi-level ODR
!            (for the 'good' vertical levels only)
!-------------------------------------------------------------------------------

! ----------------------------------------------------
! auxilliary: observation blacklist for current report
! ----------------------------------------------------
! for each level and variable, determine whether it is in a blacklisted interval
! important: to be independent from model state, 'zppl' must be used here
!            (for wind profilers and VAD radars, it also makes more sense to use
!             height (converted into 'zppl' by standard atmosphere) than
!             pressure for blacklisting (as height is more closely related to
!             radar range, which should be the most relevant criterion here))

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      kblk (klev,1) = 0
      kblk (klev,2) = 0
      kblk (klev,3) = 0
      kblk (klev,4) = 0
      DO intv = 1 , blk_loc(irpl)% ndim
        IF (      (zppl(ilev) <= blk_loc(irpl)%plow(intv) +epsy)               &
            .AND. (zppl(ilev) >= blk_loc(irpl)%pup (intv) -epsy)) THEN
          kblk (klev,blk_loc(irpl)%kvar(intv)) = 1
        ENDIF
      ENDDO
    ENDDO
    IF (lblkw(irpl)) THEN
      DO klev = 1 , nlvp
        kblk (klev,1) = 1
        kblk (klev,2) = 1
        kblk (klev,3) = 1
        kblk (klev,4) = 1
      ENDDO
    ENDIF

! -----------------------------
! level ID / level significance
! -----------------------------
! convert level identity from BUFR Table 0 08 042 into:
! - level identity as defined in VOF
! - level identity as defined for NetCDF feedback file format
! -----------------------------------------------------------

    DO klev = 1 , nlvp
      nlidvof  =  0
      nlidcdf  =  0
      lneedt (klev)  =  .FALSE.
      lneedq (klev)  =  .FALSE.
      lneedv (klev)  =  .FALSE.
      IF (izlv(isrtlvp(klev)) /= imdi) THEN
        !   get    level significance
        nlidin  =  izlv(isrtlvp(klev))
        idsurf   =  ibit1 ( nlidin, ilv_sfc  )
        idstd    =  ibit1 ( nlidin, ilv_std  )
        idtropo  =  ibit1 ( nlidin, ilv_tropo)
        iduvmax  =  ibit1 ( nlidin, ilv_max  )
        idsigt   =  ibit1 ( nlidin, ilv_sigt )
        idsigq   =  ibit1 ( nlidin, ilv_sigq )
        idsigv   =  ibit1 ( nlidin, ilv_sigv )
        !   NetCDF level significance
        nlidcdf  =  insert( 0      , idsurf , LS_SURFACE )
        nlidcdf  =  insert( nlidcdf, idstd  , LS_STANDARD)
        nlidcdf  =  insert( nlidcdf, idtropo, LS_TROPO   )
        nlidcdf  =  insert( nlidcdf, iduvmax, LS_MAX     )
        nlidcdf  =  insert( nlidcdf, IOR( IOR(idsigt,idsigq), idsigv ), LS_SIGN)
        !   VOF    level significance
        nlidvof  =  insert( 0      , idsurf , nvlidp(7) )
        IF (zppl(isrtlvp(klev)) >= 10000.0_wp-epsy) THEN
          nlidvof = insert( nlidvof, idstd  , nvlidp(6) )
        ELSE
          nlidvof = insert( nlidvof, idstd  , nvlidp(4) )
        ENDIF
        nlidvof  =  insert( nlidvof, idtropo, nvlidp(2) )
        nlidvof  =  insert( nlidvof, iduvmax, nvlidp(1) )
        nlidvof  =  insert( nlidvof, idsigv , nvlidp(8) )
        nlidvof  =  insert( nlidvof, IOR( idsigt,idsigq ), nvlidp(9) )
        ! decide whether a certain variable must be available at a certain level
        !   (used only for event counters)
        IF (idsurf == 1)  idstd = 1
        lneedt (klev)  =  (idsigt == 1) .OR. (idtropo == 1) .OR. (idstd == 1)
        lneedq (klev)  =  (idsigq == 1)                     .OR. (idstd == 1)
        lneedv (klev)  =  (idsigv == 1) .OR. (iduvmax == 1) .OR. (idstd == 1)
      ENDIF
!CS:  !   temporal solution: set standard level for towers to indicate that
      !            simulated FG is computed at height of sensor above ground
      IF (ltower)           nlidvof  =  insert( nlidvof, 1, nvlidp(6) )
      IF (ltower)           nlidcdf  =  insert( nlidcdf, 1, LS_STANDARD)
      !   always set surface significance for surface level
      IF (klev == klevsfc)  nlidvof  =  insert( nlidvof, 1, nvlidp(7) )
      IF (klev == klevsfc)  nlidcdf  =  insert( nlidcdf, 1, LS_SURFACE)

      ! fill ODR : level ID and level significance
      ! --------
      mzmlbd (klev,nbtlid) = nlidvof
      mzmlbd (klev,nbtlsg) = nlidcdf

      ! fill ODR : initialize status and quality control flag words
      ! --------
      mzmlbd (klev,nbterr) =  0
      mzmlbd (klev,nbtqcf) =  0
    ENDDO
    ! fill ODR station characteristics: vertical levels reported as height
    ! --------
    IF (indlev == 1)  CALL MVBITS( 1 , 0 , 1 , momlhd(nmlob,nhschr), nvlzbp )

! --------------------------------------------
! latitude / longitude [grid pt. units] / time
! --------------------------------------------
! Note: Consideration of changes in the horizontal observation location
! ----                   as a function of the vertical level
!       - is always done here (i.e. for 'nbtzio', 'nbtzjo')
!       - is not done when pressure is computed from height levels (call'f_z2p')
!       - depends on 'lmult_shft' when computing the rotation of the wind
!         components (call 'uv2uvrot_vec', see below) !

    ll2 = .FALSE.
    IF ((kcdftyp == ncdf_tempdesc) .AND. (ksurfob(irpl) < 0)) THEN
      ll2 =       (zdlat(-ksurfob(irpl)) > rmdich)                             &
            .AND. (zdlon(-ksurfob(irpl)) > rmdich)
    ENDIF
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)

!     IF ((r_dlon > c0) .AND. (r_dlat > c0)) THEN
      !---------------------------------------------------------
      ! A posteriori shift of horizontal location
      ! This part is applicable only to the COSMO (regular grid)
      !---------------------------------------------------------
      roblat  =  omlhed(nmlob,nhjlat)
      roblon  =  omlhed(nmlob,nhilon)

      IF ((zdlat(ilev) > rmdich) .AND. (zdlon(ilev) > rmdich) .AND. (ll2)) THEN
        !   (AMDAR: zdlat, zdlon are full lat./long., i.e. not a displacement)
        IF (lamdarml)  roblat  =  c0
        IF (lamdarml)  roblon  =  c0
        IF (kcdftyp /= ncdf_tempdesc) THEN
          roblat  =  roblat + zdlat(ilev)
          roblon  =  roblon + zdlon(ilev)
        ELSEIF (ll2) THEN
          !  (for descending TEMP: entry 'nhjlat' is lat of 'surface' obs level;
          !   the index of that level is -ksurfob);
          !   the displacement has to be relative to this surface level but
          !   'zdlat(ilev)' is reported relative to another level (typically
          !   the uppermost one))
          roblat  =  roblat + zdlat(ilev) - zdlat(-ksurfob(irpl))
          roblon  =  roblon + zdlon(ilev) - zdlon(-ksurfob(irpl))
        ENDIF
      ENDIF

!       rlon = rla2rlarot (roblat, roblon, r_pollat, r_pollon, r_polgam)
!       rlat = phi2phirot (roblat, roblon, r_pollat, r_pollon)
!              ==========
      ! fill ODR     (observation position in geographical coordinates)
      ! -------- (previously: obs osition in grid point units (total area !))
!     IF ((r_dlon > c0) .AND. (r_dlat > c0)) THEN
!       zmlbdy (klev,nbtzio) = c1 + (rlon - r_startlon_tot) /r_dlon
!       zmlbdy (klev,nbtzjo) = c1 + (rlat - r_startlat_tot) /r_dlat
!     ENDIF
      zmlbdy (klev,nbtzio) = roblon
      zmlbdy (klev,nbtzjo) = roblat
      zmlbdy (klev,nbttim) = omlhed(nmlob,nhtime) + izdt(ilev) *c3600r
      !  (for descending TEMP: entry 'nhtime' relates to lowest obs level)
      IF (kcdftyp == ncdf_tempdesc)                                            &
        zmlbdy (klev,nbttim) =  omlhed(nmlob,nhtime)                           &
                              + (izdt(ilev) - izdt(-ksurfob(irpl))) *c3600r
    ENDDO

! -------------------
! height and pressure
! -------------------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflag   =  0
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

      IF (     (zzz(ilev) <= rmdich) .OR. (kproclv(ilev) <= -7)                &
          .OR. (zpp(ilev) <= rmdich)) THEN
!         .OR. (zpp(ilev) <= rmdich) .OR. (zstalt <= rmdich)) THEN
!       IF ((nzex(kobtyp) >= 1) .AND. (ioberr(klev) == 1))                     &
        IF ((nzex(kobtyp) >= 1) .AND. (kproclv(ilev) <= -9))                   &
          neventd (nepmis,icma) = neventd(nepmis,icma) + iactr
        ioberr (klev)  =  0
      ELSE
        !---   if level (reported height) is below surface (station height)
        IF (zstalt > rmdich) THEN
          IF (zzz(ilev) <= zstalt-c2) THEN
            IF (ilev /= ilevsfc) THEN
              !   level below surface, means that also T,q,uv shall
              !   not be used actively !  ( --> set kproclv=-kproclv )
              IF (ioberr(klev) == 1)  neventd (nelext,icma) =                  &
                                      neventd (nelext,icma) + iactr
              nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
              nflag          =  IBSET( nflag      , nvflbp )
              kproclv (ilev) =  - kproclv(ilev)
            ENDIF
            ioberr (klev)  =  0
          ENDIF
        ENDIF
        !---   if pressure is derived from height (or vice versa)
        IF (ABS( kproclv(ilev) ) == 1) THEN
          IF (lamdarml)                                                        &
            nflgx (klev) =  IBSET( nflgx(klev), nvfbps(5) )
          IF ((ABS( kproclv(ilev) ) == 1) .AND. (.NOT. lamdarml))              &
            nflgx (klev) =  IBSET( nflgx(klev), nvfbps(6) )
          ioberr (klev)  =  0
        ENDIF
        IF (ABS( kproclv(ilev) ) == 2) THEN
          !---   geopotential without temperature is set passive
          IF ((ztt(ilev) < rmdich) .OR. (ztt(ilev) < c1)) THEN
            IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                    &
                                    neventd (nepflg,icma) + iactr
            nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
            ioberr (klev)  =  0
          !---   if geopotential blacklisted
          ELSEIF (kblk(klev,2) == 1) THEN
            IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                    &
                                    neventd (nepflg,icma) + iactr
            nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
            ioberr (klev)  =  0
          ENDIF
        ENDIF
      ENDIF

      ! fill ODR (except for flags)
      ! --------
      zmlbdy (klev,nbtp  ) = zpp   (ilev)
      zmlbdy (klev,nbtz  ) = zzz   (ilev)
      zmlbdy (klev,nbtplv) = zppl  (ilev)
      zmlbdy (klev,nbtzer) = zoberr(klev)
      mzmlbd (klev,nbtflg) = nflag
      !   bogus observations for semi-idealised tests
      IF (lbogus)  zmlbdy (klev,nbtp)   = zmlbdy(klev,nbtp)   + 100._wp*fbogus
      IF (lbogus)  zmlbdy (klev,nbtplv) = zmlbdy(klev,nbtplv) + 100._wp*fbogus
      !   auxilliary ODR element
      zmlbdy (klev,nbtlop) = LOG ( zmlbdy(klev,nbtp) )

!     IF (zpp(ilev) > 90000.0_wp)                                              &
!       PRINT *,'zzz4 ', yomlhd(nmlob), ilev, zmlbdy(klev,nbtzer)              &
!                                           , zmlbdy(klev,nbtz  )
    ENDDO

    !   if (.not. luse_mlz) then
    !   set 'height'/passive flag for all obs except for TEMP surface level,
    !        or (if no surface level exists) for lowest active TEMP z-obs
    klvz  =  klevsfc
    IF ((.NOT. luse_mlz) .AND. (klevsfc == 0) .AND. (ltempasc)                 &
                         .AND. (nlvp > 0)) THEN
      !   if no surface obs level, get lowest obs level with active p-z obs
      klvz = 1
      DO WHILE ((klvz < nlvp) .AND.(.NOT.BTEST(mzmlbd(klvz,nbterr),nvrz)))
        klvz  =  klvz + 1
      ENDDO
      IF ((klvz == nlvp) .AND. (.NOT.BTEST( mzmlbd(klvz,nbterr),nvrz )))       &
        klvz  =  0
    ENDIF
    IF (.NOT. ltempasc)  klvz = 0
    DO klev = 1 , nlvp
      IF ((.NOT. luse_mlz) .AND. (klev /= klvz)) THEN
        ilev = isrtlvp(klev)
        IF ((zzz(ilev) > rmdich) .AND. (ABS( kproclv(ilev) ) == 2)) THEN
          IF (ioberr(klev) == 1)  neventd (nepflg,icma) =                      &
                                  neventd (nepflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF

      ! fill ODR: flags only
      ! --------
      mzmlbd (klev,nbterr) = insert( mzmlbd(klev,nbterr), ioberr(klev), nvrz   )
      mzmlbd (klev,nbtflg) = insert( mzmlbd(klev,nbtflg), nflgx (klev), nvfzbp )
    ENDDO

! -----------
! temperature
! -----------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

! if temperature missing or below absolute minimum (no flag in 'nflgx' set)
      IF ((ztt(ilev) < rmdich) .OR. (ztt(ilev) < c1)) THEN
        IF ((ntex(kobtyp) >= 1) .AND. (lneedt(klev)))                          &
          neventd (netmis,icma) = neventd(netmis,icma) + iactr
        ztt    (ilev)  =  rmdi
        ioberr (klev)  =  0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if temperature blacklisted
        IF (kblk(klev,3) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if temperature gross error
        IF (     (ztt(ilev) < tmelt-90._wp)                                    &
           .OR.  (ztt(ilev) > tmelt+60._wp)                                    &
           .OR. ((ztt(ilev) > tmelt+20._wp).AND. (zppl(ilev) < 70000._wp))     &
           .OR. ((ztt(ilev) > tmelt+ 5._wp).AND. (zppl(ilev) < 50000._wp))     &
           .OR. ((ztt(ilev) > tmelt- 5._wp).AND. (zppl(ilev) < 40000._wp))) THEN
          IF (ioberr(klev) == 1)  neventd (netext,icma) =                      &
                                  neventd (netext,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      IF (      (kcdtyp == nra_eu)                                             &
          .AND. ((ztt(ilev) > rmdich) .OR. (zww(ilev) > rmdich))) THEN
! RASS profiler quality information (also valid for vertical wind)
        IF ((iqci (ilev) < 0) .OR. (iqci(ilev) > nibits(nvfboc(1))))           &
          iqci (ilev) = nibits(nvfboc(1))
        CALL MVBITS( iqci(ilev) , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
        IF (iqci (ilev) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          ioberr (klev)  =  0
        ENDIF
! RASS profiler signal to noise ratio (also valid for vertical wind)
        IF ((zsinor(ilev) > rmdich) .AND. (zsinor(ilev) < -20.0_wp)) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          CALL MVBITS( 1 , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      IF ((ltower) .AND. (ztt(ilev) > rmdich)) THEN
! tower report quality information
        iqcbit   =  ibit1 ( iqci(ilev),  3 )
        iqcbit2  =  ibit1 ( iqci(ilev), 13 )
        CALL MVBITS( iqcbit  , 0 , nvfboc(5) , nflgx(klev) , nvfbps(5) )
        CALL MVBITS( iqcbit2 , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
        IF ((iqcbit == 1) .OR. (iqcbit2 == 1)) THEN
          IF (ioberr(klev) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          ioberr (klev)  =  0
        ENDIF
      ENDIF
    ENDDO

! Convert RASS / SODAR vitual temperatures into normal air temperatures: NO !!!
! ---------------------------------------------------------------------  --
!TODO ( --> 'tv2t' is not yet vectorised)
!   IF ((kcdtyp = nra_eu) THEN
!     DO klev = 1 , nlvp
!       ilev = isrtlvp(klev)
!       IF ((ztt(ilev) > rmdich) .AND. (ztt(ilev) > tmelt-90._wp)              &
!                                .AND. (ztt(ilev) < tmelt+60._wp)) THEN
!         ztv (ilev)  =  ztt(ilev)

!         ztt (ilev)  =  tv2t ( iob, job, ztv(ilev), zzz(ilev) )
!                        ====
!       ENDIF
!     ENDDO
!   ENDIF

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)

! bogus observations for semi-idealised tests
      IF (      (lbogus)             .AND. (zppl(ilev) >= 66999._wp)           &
          .AND. (ztt(ilev) > rmdich) .AND. (zppl(ilev) <= 86000._wp))          &
        ztt (ilev)  =  MAX( 100._wp, ztt(ilev) + 1._wp*fbogus )

      ! the following is obsolete since 'shift_profile'
      ! (in src_obs_operator_conv.f90) is applied !!
      ! height correction for tower temperature profile observations
      ! in order to adjust them to their re-assigned height (such that the
      ! height of sensor about ground is the same for model equiv. and obs)
      ! (basic climatological lapse rate is considered good enough, assuming
      !  distance of model orography to station height is small for tower sta.)
!     IF ((ltower) .AND. (ztt(ilev) > rmdich))                                 &
!       ztt (ilev)  =  ztt(ilev) - .0065_wp *(omlhed(nmlob,nhsurf) - zstalt)

! fill ODR
! --------
      zmlbdy (klev,nbtt  ) = ztt   (ilev)
      zmlbdy (klev,nbtter) = zoberr(klev)
      mzmlbd (klev,nbterr) = insert( mzmlbd(klev,nbterr), ioberr(klev), nvrt  )
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      mzmlbd (klev,nbtflg) = insert( mzmlbd(klev,nbtflg), nflgx(klev), nvftbp )
    ENDDO

! --------
! humidity     (dew point, relative humidity and/or mixing ratio as input,
! --------      relative humidity is written to ODR)
!              (if input is Td or qx, then temperature is also required )
!               note: processing below is done as if qx could be non-missing,
!                     but currently, it is not implemented to be read

    lexitd = .FALSE.
    lexiqx = .FALSE.
    lexirh = .FALSE.
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1

! reported humidity is fully discarded if negative or exactly zero
!   (sometimes reported zero values denote undefined values)
      IF (zrh(ilev) <= epsy) THEN
        zrh    (ilev)  =  rmdi
!       ioberr (ilev)  =  0   ! is set further below if humidity is missing
!       IF (ioberr(ilev) == 1)  neventd (neqflg,icma) =       & ! also set
!                               neventd (neqflg,icma) + iactr   ! below
      ENDIF
      IF (ztd(ilev) <= b4w+c1   )  ztd (ilev)  =  rmdi
      IF (zqx(ilev) <= 1.E-10_wp)  zqx (ilev)  =  rmdi

! as long as reported humidity of any kind is converted into relative humidity
! in order to be processed further (as is currently the case), then:
! dewpoint temp. (mixing ratio) w/o temperature (T,p) is completely discarded
      ll1  =  (ztd(ilev) > rmdich) .AND.      (ztt(ilev) < rmdich)
      ll2  =  (zqx(ilev) > rmdich) .AND. (MIN( ztt(ilev),zpp(ilev) ) < rmdich)
      IF ((ll1) .OR. (ll2)) THEN
        IF (ll1)  ztd (ilev)  =  rmdi
        IF (ll2)  zqx (ilev)  =  rmdi
        IF (MAX( ztd(ilev), zqx(ilev), zrh(ilev) ) < rmdich) THEN
          IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (klev) = IBSET ( nflgx(klev), nvfbps(5) )
          ioberr (klev) = 0
        ENDIF
      ENDIF

! if (relative or dewpoint) humidity missing (no flag in 'nflgx' set)
      IF (MAX( zrh(ilev), ztd(ilev), zqx(ilev) ) < rmdich) THEN
        IF (      (ioberr(klev) == 1) .AND. (ntdex(kobtyp) >= 1)               &
            .AND. (lneedq(klev)) .AND. (zppl(ilev) >= pqmin))                  &
          neventd (neqmis,icma) = neventd(neqmis,icma) + iactr
        zrh    (ilev) = rmdi
        ztd    (ilev) = rmdi
        zqx    (ilev) = rmdi
!       zrhc   (ilev) = rmdi
        ioberr (klev) = 0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
        IF (ltower) THEN
          IF ((ztd(ilev) > rmdich) .AND. (zrh(ilev) > rmdich)) THEN
            ! if RH bad and Td good, then discard RH, and vice versa
            IF (      (     (BTEST( iqci(ilev), 4 ))                           &
                       .OR. (BTEST( iqci(ilev),14 )))                          &
                .AND. (.NOT. BTEST( iqci(ilev), 5 ))                           &
                .AND. (.NOT. BTEST( iqci(ilev),15 )) ) THEN
              zrh (ilev) = rmdich
            ELSEIF (      (     (BTEST( iqci(ilev), 5 ))                       &
                           .OR. (BTEST( iqci(ilev),15 )))                      &
                    .AND. (.NOT. BTEST( iqci(ilev), 4 ))                       &
                    .AND. (.NOT. BTEST( iqci(ilev),14 )) ) THEN
              ztd (ilev) = rmdich
            ENDIF
          ENDIF
          IF (ztd(ilev) <= rmdich)  iqci(ilev) = IBCLR( iqci(ilev), 5 )
          IF (ztd(ilev) <= rmdich)  iqci(ilev) = IBCLR( iqci(ilev),15 )
          IF (zrh(ilev) <= rmdich)  iqci(ilev) = IBCLR( iqci(ilev), 4 )
          IF (zrh(ilev) <= rmdich)  iqci(ilev) = IBCLR( iqci(ilev),14 )
          ! at this point, if 1 QC bit is set, then neither RH nor TD is good
          IF (     (BTEST( iqci(ilev),4 )) .OR. (BTEST( iqci(ilev),14 ))       &
              .OR. (BTEST( iqci(ilev),5 )) .OR. (BTEST( iqci(ilev),15 ))) THEN
            IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                    &
                                    neventd (neqflg,icma) + iactr
            ioberr (klev)  =  0
            IF ((BTEST( iqci(ilev), 4 )) .OR. (BTEST( iqci(ilev), 5 )))        &
              CALL MVBITS( 1 , 0 , nvfboc(5) , nflgx(klev) , nvfbps(5) )
            IF ((BTEST( iqci(ilev),14 )) .OR. (BTEST( iqci(ilev),15 )))        &
              CALL MVBITS( 1 , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
          ENDIF
        ENDIF
! if humidity blacklisted
        IF (kblk(klev,4) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if humidity not in valid height range
        IF (zppl(ilev) < pqmin) THEN
          IF (ioberr(klev) == 1)  neventd (neq300,icma) =                      &
                                  neventd (neq300,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if dew point gross error
        IF (      (ztd(ilev) > rmdich)                                         &
            .AND. (      (ztd(ilev) < tmelt-150._wp)                           &
                   .OR.  (ztd(ilev) > tmelt+ 40._wp)                           &
                   .OR. ((ztd(ilev) < tmelt- 90._wp) .AND. (ilev == ilevsfc))  &
                   .OR.  (ztd(ilev) > ztt(ilev)+2._wp))) THEN
          IF (MAX( zrh(ilev),zqx(ilev) ) < rmdich) THEN
            IF (ioberr(klev) == 1)  neventd (neqlow,icma) =                    &
                                    neventd (neqlow,icma) + iactr
            nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
            ioberr (klev)  =  0
            ! (minimum ztd must be > max( b4w,b4i ) to avoid division by zero)
            ztd    (ilev)  =  MAX( ztd(ilev), b4w+c1 )
            ztd    (ilev)  =  MIN( MIN( ztd(ilev), ztt(ilev)+c2 ), tmelt+40._wp)
          ELSE
            ztd    (ilev)  =  rmdi
          ENDIF
        ENDIF
! if mixing ratio gross error  , i.e. not in (1E-10 < zqx(w) <= 0.03)
        IF ((zqx(ilev) > rmdich) .AND. (     (zqx(ilev) > 0.03_wp)             &
                                        .OR. (zqx(ilev) <= 1.E-10_wp))) THEN
          IF (MAX( zrh(ilev),ztd(ilev) ) < rmdich) THEN
                                      necnt = neqlow
            IF (zqx(ilev) > 0.03_wp)  necnt = neqbig
            IF (ioberr(klev) == 1)  neventd (necnt ,icma) =                    &
                                    neventd (necnt ,icma) + iactr
            nflgx  (klev)  =  IBSET ( nflgx(klev), nvfbps(3) )
            ioberr (klev)  =  0
            zqx    (ilev)  =  MAX( MIN( zqx(ilev) , 0.03_wp ) , 1.E-10_wp )
          ELSE
            zqx    (ilev)  =  rmdi
          ENDIF
        ENDIF
! if relative humidity gross error  , i.e. not in (epsy < zrh <= 102%)
        IF ((zrh(ilev) > rmdich) .AND. (     (zrh(ilev) > rtshlm)              &
                                        .OR. (zrh(ilev) <= epsy) )) THEN
          IF (MAX( ztd(ilev),zqx(ilev) ) < rmdich) THEN
                                      necnt = neqlow
            IF (zrh(ilev) >  rtshlm)  necnt = neqbig
            IF (ioberr(klev) == 1)  neventd (necnt ,icma) =                    &
                                    neventd (necnt ,icma) + iactr
            nflgx  (klev)  =  IBSET ( nflgx(klev), nvfbps(3) )
            ioberr (klev)  =  0
!           zrh    (ilev)  =  MAX( MIN( zrh(ilev) , c1 ) , c05*epsy )
            zrh    (ilev)  =  MAX( MIN( zrh(ilev) , c1 ) , c0 )
          ELSE
            zrh    (ilev)  =  rmdi
          ENDIF
        ENDIF
! set dew point / mixing ratio to passive if temperature is not active
        IF (      (.NOT. BTEST( mzmlbd(klev,nbterr),nvrt ))                    &
            .AND. (zrh(ilev) < rmdich) .AND. (     (ztd(ilev) > rmdich)        &
                                              .OR. (zqx(ilev) > rmdich))) THEN
          IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
          ioberr (klev)  =  0
        ENDIF
      ENDIF
    ENDDO

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      zttl (klev)  =  ztt (ilev)
      ztdl (klev)  =  ztd (ilev)
      zqxl (klev)  =  zqx (ilev)
      zrhw (klev)  =  zrh (ilev)
      zrhc (klev)  =  zrh (ilev)
      IF (ztdl(klev) > rmdich)  lexitd = .TRUE.
      IF (zqxl(klev) > rmdich)  lexiqx = .TRUE.
      IF (zrhw(klev) > rmdich)  lexirh = .TRUE.
    ENDDO
! at this point, 'zrhw', 'ztdl', 'zqxl' and 'zrhc' should have either reasonable
! or missing values, which should prevent lethal crashes in the following calls

! compute model compatible relative humidity
! ------------------------------------------

    IF (nlvp > 0) THEN
      ALLOCATE ( zrhc2 (nlvp) , STAT=istat )
      ALLOCATE ( zrhc3 (nlvp) , STAT=istat )
      zrhc2 (:) = rmdi
      zrhc3 (:) = rmdi
    ENDIF

    IF ((lexirh) .AND. (nlvp > 0)) THEN

      CALL obs_rhw2rh ( madj_hum, nlvp, zttl, zrhw                             &
                      , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhc )
!     ===============
    ENDIF

    IF ((lexitd) .AND. (nlvp > 0)) THEN
      ALLOCATE ( zrhw2 (nlvp) , STAT=istat )

      CALL obs_td2rh  ( madj_hum, nlvp, zttl , ztdl                            &
                      , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhw2, zrhc2 )
!     ==============
    ENDIF

    IF ((lexiqx) .AND. (nlvp > 0)) THEN
      ALLOCATE ( zrhw3 (nlvp) , STAT=istat )
      ALLOCATE ( zqvw  (nlvp) , STAT=istat )
      ALLOCATE ( zqv   (nlvp) , STAT=istat )
      zqvw (:) = rmdi
      !   here, 'zpp' instead of 'zppl' should (and can) be used

      CALL obs_qx2rh ( madj_hum, nlvp, zttl, zpp, zqxl , zqvw                  &
                     , rmdi,b1,b2w,b2i,b3,b4w,b4i,rdv, tmelt, zqv, zrhw3, zrhc3)
!     ==============
    ENDIF

!NEC$ nomove
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
! check consistency between reported and derived relative humidity values
! -----------------------------------------------------------------------
      ll1 = .FALSE.
!     IF ((lexitd) .AND. (lexirh)) THEN
!   if dewpoint-derived rel. humidity and reported rel. humidity differ by > 4%
        IF ((zrhc2(klev) > rmdich) .AND. (zrhc(klev) > rmdich)) THEN
          IF (ABS( zrhc2(klev) - zrhc(klev) ) > 0.04_wp)  ll1 = .TRUE.
        ENDIF
!     ENDIF
!     IF ((lexiqx) .AND. (lexirh)) THEN
!   if mixing-ratio-derived rel. humidity and reported rel. humid. differ by >4%
        IF ((zrhc3(klev) > rmdich) .AND. (zrhc(klev) > rmdich)) THEN
          IF (ABS( zrhc3(klev) - zrhc(klev) ) > 0.04_wp)  ll1 = .TRUE.
        ENDIF
!     ENDIF
!     IF ((lexiqx) .AND. (lexitd)) THEN
!   if mixing-ratio-derived and dewpoint-derived rel. humidity differ by >4%
        IF ((zrhc3(klev) > rmdich) .AND. (zrhc2(klev) > rmdich)) THEN
          IF (ABS( zrhc3(klev) - zrhc2(klev) ) > 0.04_wp)  ll1 = .TRUE.
        ENDIF
!     ENDIF
      IF (ll1) THEN
        IF (ioberr(klev) == 1)  neventd (neqflg,icma) =                        &
                                neventd (neqflg,icma) + iactr
        nflgx  (klev) = IBSET ( nflgx(klev), nvfbps(5) )
        ioberr (klev) = 0
      ENDIF

! if relative humidity is not reported ('zrhc' = missing value)
! then use the dewpoint- or mixing-ratio-derived value 'zrhc2' resp. 'zrhc3'
! --------------------------------------------------------------------------
!     IF (lexitd)  PRINT *, 'zrhtd ',zrhc2(klev)
!     IF (lexiqx)  PRINT *, 'zrhqx ',zrhc3(klev)
      IF ((lexitd) .AND. (zrhc2(klev) > rmdich).AND. (zrhc(klev) < rmdich)) THEN
        zrhc  (klev)  =  zrhc2(klev)
        zrhw  (klev)  =  zrhw2(klev)
        !   (if 'zrhw2' has gross error, use 'zrhc3' if well-defined)
        IF ((zrhw2(klev) > rtshlm) .AND. (zrhc3(klev) > rmdich))               &
          zrhc  (klev)  =  rmdi
      ENDIF
      IF ((lexiqx) .AND. (zrhc3(klev) > rmdich).AND. (zrhc(klev) < rmdich)) THEN
        zrhc  (klev)  =  zrhc3(klev)
        zrhw  (klev)  =  zrhw3(klev)
      ENDIF

! gross error if original (i.e. not model compatible) relative humidity > 102 %
!                         (previous version: if > 120 %)
      IF ((zrhw(klev) >  rtshlm) .AND. (zrhw(klev) > rmdich)) THEN
        IF (ioberr(klev) == 1)  neventd (neqbig,icma) =                        &
                                neventd (neqbig,icma) + iactr
        nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
        ioberr (klev)  =  0
        zrhc   (klev)  =  MAX( MIN( zrhc(klev) , c1 ) , epsy )

! bias correction of relative humidity (near saturation only)
! ------------------------------------
      ELSEIF ((zrhc(klev) >= rhtsat(ixrhs)) .AND. (zrhc(klev) < c1-epsy)) THEN
        IF ((zttl(klev) <  tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqsam,icma) = neventd(neqsam,icma) + iactr
        IF ((zttl(klev) >= tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqsap,icma) = neventd(neqsap,icma) + iactr
!       nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
        nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(6) )
      ELSEIF (zrhc(klev) > c1+epsy) THEN
        IF ((zttl(klev) <  tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqclm,icma) = neventd(neqclm,icma) + iactr
        IF ((zttl(klev) >= tmelt) .AND. (ioberr(klev) == 1))                   &
          neventd (neqclp,icma) = neventd(neqclp,icma) + iactr
!       nflgx (klev)  =  IBSET( nflgx(klev), nvfbps(6) )
      ENDIF
      IF (zrhc(klev) >= rhtsat(ixrhs)) THEN
        zrhc  (klev)  =  c1
!       ztdl  (klev)  =  zttl(klev)  ! ztdl is not model compatible / bias corr.
      ENDIF

! bogus observations for semi-idealised tests
      IF (      (lbogus)              .AND. (zppl(ilev) >= 66999._wp)          &
          .AND. (zrhc(klev) > rmdich) .AND. (zppl(ilev) <= 86000._wp))         &
        zrhc(klev) = MAX( 0.01_wp , MIN( c1, zrhc(klev) +0.1_wp*fbogus ) )

! fill ODR
! --------
      zmlbdy (klev,nbtrh ) = zrhc  (klev)
      zmlbdy (klev,nbtqer) = zoberr(klev)
      mzmlbd (klev,nbterr) = insert( mzmlbd(klev,nbterr), ioberr(klev), nvrq )
      IF ((zrhc(klev) > rmdich) .AND. (zrhw(klev) > rmdich)) THEN
        zmlbdy (klev,nbtdrh) = zrhc (klev)  -  zrhw(klev)
      ELSEIF (zrhc(klev) > rmdich) THEN
        zmlbdy (klev,nbtdrh) = c0
      ENDIF
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      mzmlbd (klev,nbtflg) = insert( mzmlbd(klev,nbtflg), nflgx(klev), nvfqbp )
    ENDDO

    IF  (nlvp > 0)                 DEALLOCATE ( zrhc2 , STAT=istat )
    IF  (nlvp > 0)                 DEALLOCATE ( zrhc3 , STAT=istat )
    IF ((nlvp > 0) .AND.(lexitd))  DEALLOCATE ( zrhw2 , STAT=istat )
    IF ((nlvp > 0) .AND.(lexiqx)) THEN
      DEALLOCATE ( zrhw3 , STAT=istat )
      DEALLOCATE ( zqvw  , STAT=istat )
      DEALLOCATE ( zqv   , STAT=istat )
    ENDIF

! ---------------
! horizontal wind
! ---------------

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      nflgx  (klev)  =   0
      zoberr (klev)  =  c0
      ioberr (klev)  =   1
      nflag          =   0

      ! vertical velocity
      IF ((lprofw) .AND. (zww(ilev) > rmdich)) THEN
        !   (bad reporting practice: zero wind probably means missing value)
        IF (ABS( zww(ilev) ) < epsy)  zww (ilev) = rmdi
      ENDIF

! if wind speed or direction missing (no flag in 'nflgx' set)
      IF ((zff(ilev) <= rmdich) .OR. (zdd(ilev) <= rmdich)) THEN
        IF ((nuex(kobtyp) >= 1) .AND. (lneedv(klev))) THEN
          IF (zff(ilev) <= rmdich) neventd (nefmis,icma) =                     &
                                   neventd (nefmis,icma) + iactr
          IF (zdd(ilev) <= rmdich) neventd (nedmis,icma) =                     &
                                   neventd (nedmis,icma) + iactr
        ENDIF
        zff    (ilev)  =  rmdi
        zdd    (ilev)  =  rmdi
        ioberr (klev)  =  0
      ELSE
! if level below surface
        IF ((kproclv(ilev) == -2) .OR. (kproclv(ilev) == -1)) THEN
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(4) )
          ioberr (klev)  =  0
        ENDIF
! if wind blacklisted
        IF (kblk(klev,1) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (nedflg,icma) =                      &
                                  neventd (nedflg,icma) + iactr
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(2) )
          ioberr (klev)  =  0
        ENDIF
! if wind speed / direction gross error
        IF (      (ABS( zdd(ilev) ) > 360._wp)                                 &
            .OR.  (zff(ilev) > 150._wp)                                        &
            .OR. ((zff(ilev) >  90._wp) .AND. (zppl(ilev) < 70000._wp))        &
            .OR.  (zff(ilev) < -epsy)) THEN
          IF (ioberr(klev) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(3) )
          ioberr (klev)  =  0
          zff    (ilev)  =  MAX( zff(ilev) , c0 )
        ENDIF
        IF (zff(ilev) < epsy)  zdd (ilev)  =  c0
! if VAD wind speed bad reporting practice
        IF ((zff(ilev) < fflim_vad) .AND. (kcdtyp == nravad)) THEN
          IF (ioberr(klev) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (klev)  =  IBSET( nflgx(klev), nvfbps(5) )
          ioberr (klev)  =  0
        ENDIF
! aircraft roll angle as revised quality flag for wind
        IF (lamdarml) THEN
          IF ((iroll(ilev) < 0) .OR. (iroll(ilev) > nibits(nvfboc(1))))        &
            iroll (ilev) = nibits(nvfboc(1))
!         CALL MVBITS( iroll(ilev) , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
          nflag            =  iroll(ilev)
          IF (iroll(ilev) == 1)  ioberr (klev)  =  0
          !  roll angle in station characteristics header word ('nvaaoc' bits)
          IF (iroll(ilev) == nibits(nvfboc(1)))  iroll (ilev) = nibits(nvaaoc)
          CALL MVBITS( iroll(ilev), 0 , nvaaoc , momlhd(nmlob,nhschr), nvaabp )
        ELSEIF (ltower) THEN
          iqcbit   =  MIN( 1 , ibit1( iqci(ilev), 1 ) + ibit1( iqci(ilev), 2 ) )
          iqcbit2  =  MIN( 1 , ibit1( iqci(ilev),11 ) + ibit1( iqci(ilev),12 ) )
          CALL MVBITS( iqcbit  , 0 , nvfboc(5) , nflgx(klev) , nvfbps(5) )
          IF (iqcbit2 == 1)  nflag  =  1
          IF ((iqcbit == 1) .OR. (iqcbit2 == 1)) THEN
            IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                    &
                                    neventd (nefflg,icma) + iactr
            ioberr (klev)  =  0
          ENDIF
        ELSEIF (zstdff(ilev) > 0.1_wp) THEN
! wind profiler / radar VAD standard deviation wind speed defined
          IF (zstdff(ilev) > 10._wp) THEN
            IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                    &
                                    neventd (nefflg,icma) + iactr
            nflag          =  1
            ioberr (klev)  =  0
          ENDIF
          !  adjust observation error
          zoberr (klev)  =  zstdff(ilev)
        ENDIF
      ENDIF
      IF ((lprofw) .AND. ((zff(ilev) > rmdich) .OR. (zww(ilev) > rmdich))) THEN
! wind profiler / radar VAD quality information (also valid for vertical wind)
        IF ((iqci (ilev) < 0) .OR. (iqci(ilev) > nibits(nvfboc(1))))           &
          iqci (ilev) = nibits(nvfboc(1))
        IF (nflag /= 1)  nflag  =  iqci(ilev)
        IF (iqci (ilev) == 1) THEN
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          ioberr (klev)  =  0
        ENDIF
! wind profiler (radar VAD) signal to noise ratio (also valid for vertical wind)
        IF ((zsinor(ilev) > rmdich) .AND. (zsinor(ilev) < -20.0_wp)) THEN
          IF (ioberr(klev) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          nflag          =  1
          ioberr (klev)  =  0
        ENDIF
      ENDIF
      CALL MVBITS( nflag , 0 , nvfboc(1) , nflgx(klev) , nvfbps(1) )
    ENDDO

! Transformation of wind speed and direction to wind components,
! within COSMO (not DACE) also convert to the rotated model coordinate system
! ---------------------------------------------------------------------------
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (zff(ilev) > rmdich) THEN
        zuu (klev,1)  =  zff(ilev) * (-SIN( zdd(ilev)*r_degrad ))
        zvv (klev,1)  =  zff(ilev) * (-COS( zdd(ilev)*r_degrad ))
      ELSE
        zuu (klev,1)  =  c0
        zvv (klev,1)  =  c0
      ENDIF
      zrlat (klev,1)  =  omlhed(nmlob,nhjlat)
      zrlon (klev,1)  =  omlhed(nmlob,nhilon)
      IF ((zdlat(ilev) > rmdich) .AND. (lmult_shft))                           &
        zrlat (klev,1)  =  zrlat(klev,1) + zdlat(ilev)
      IF ((zdlon(ilev) > rmdich) .AND. (lmult_shft))                           &
        zrlon (klev,1)  =  zrlon(klev,1) + zdlon(ilev)
    ENDDO

    IF (nlvp > 0) THEN
#ifdef __COSMO__
      CALL uv2uvrot_vec ( zuu, zvv, zrlat, zrlon, r_pollat, r_pollon, nlvp, 1 )
!     =================
#endif
    ENDIF

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)

! bogus observations for semi-idealised tests
      IF (lbogus) THEN
        IF ((zppl(ilev) > 55999._wp) .AND. (zppl(ilev) < 86001._wp)) THEN
          zuu (klev,1)  =  zuu (klev,1) + 1._wp*fbogus
          zvv (klev,1)  =  zvv (klev,1) + 1._wp*fbogus
        ELSEIF (zppl(ilev) > 25999._wp .AND. zppl(ilev) < 54001._wp) THEN
          zuu (klev,1)  =  zuu (klev,1) + 1._wp*fbogus
        ELSEIF (zppl(ilev) < 24001._wp) THEN
          zuu (klev,1)  =  zuu (klev,1) - 1._wp*fbogus
          zvv (klev,1)  =  zvv (klev,1) - 1._wp*fbogus
        ENDIF
      ENDIF

! fill ODR
! --------
      IF (zff(ilev) <= rmdich)  zuu (klev,1)  =  rmdi
      IF (zff(ilev) <= rmdich)  zvv (klev,1)  =  rmdi
      zmlbdy (klev,nbtu  ) = zuu   (klev,1)
      zmlbdy (klev,nbtv  ) = zvv   (klev,1)
      zmlbdy (klev,nbtuer) = zoberr(klev)
      mzmlbd (klev,nbterr) = insert( mzmlbd(klev,nbterr), ioberr(klev), nvru )
! note: no data base flag set, override flag is set in obs header, not here;
!       therefore only data base flag and level below surface flag are set
      mzmlbd (klev,nbtflg) = insert( mzmlbd(klev,nbtflg), nflgx(klev), nvfubp )
      zmlbdy (klev,nbtw  ) = zww   (ilev)
!     zmlbdy (klev,nbtsnr) = zsinor(ilev)
      IF (ltower) THEN
        zmlbdy (klev,nbtuac) = zazi  (ilev)
        zmlbdy (klev,nbtsnr) = zsinor(ilev)
      ELSEIF (zstdff(ilev) > 0.1_wp) THEN
        zmlbdy (klev,nbtuac) = zstdff(ilev)
        zmlbdy (klev,nbtsnr) = zsinor(ilev)
      ENDIF
    ENDDO

! ------------------
! observation errors
! ------------------
! get the observation errors only at the end for those observations
! for which the errors need to be specified
! (i.e. for which the errors have been set to zero previously
! -----------------------------------------------------------------

    DO klev = 1 , nlvp

      CALL obs_find_level ( nerlev, rolnlv, zmlbdy(klev,nbtlop), ilverr, fiperr)
!     ===================

      IF (zmlbdy(klev,nbtz ) > rmdich) THEN
        zmlbdy(klev,nbtzer) =       fiperr * oezsond(ilverr)                   &
                              + (c1-fiperr)* oezsond(ilverr+1)
        IF (lamdarml)                                                          &
          mzmlbd (klev,nbterr) = IBCLR( mzmlbd(klev,nbterr), nvrz )
      ENDIF
      IF (zmlbdy(klev,nbtu ) > rmdich) THEN
        IF (lamdarml) THEN
          zmlbdy(klev,nbtuer) =       fiperr * oeairep(ilverr)                 &
                                + (c1-fiperr)* oeairep(ilverr+1)
        ELSEIF (zmlbdy(klev,nbtuer) < +epsy) THEN
          zmlbdy(klev,nbtuer) =       fiperr * oevsond(ilverr)                 &
                                + (c1-fiperr)* oevsond(ilverr+1)
        ELSE
!    standard deviation wind speed + 0.5 * radiosonde wind error
          zmlbdy(klev,nbtuer) =   zmlbdy(klev,nbtuer)                          &
                                +     fiperr * oevsond(ilverr)   *c05          &
                                + (c1-fiperr)* oevsond(ilverr+1) *c05
        ENDIF
      ENDIF
      IF (zmlbdy(klev,nbtt ) > rmdich) THEN
        IF (lamdarml) THEN
          zmlbdy(klev,nbtter) =       fiperr * oetairp(ilverr)                 &
                                + (c1-fiperr)* oetairp(ilverr+1)
        ELSE
          zmlbdy(klev,nbtter) =       fiperr * oetsond(ilverr)                 &
                                + (c1-fiperr)* oetsond(ilverr+1)
        ENDIF
      ENDIF
      IF (zmlbdy(klev,nbtrh) > rmdich) THEN
                                          roberr  =  rherr1 * c100r
        IF (zmlbdy(klev,nbtrh) < rrhlim)  roberr  =  rherr2 * c100r
        IF (zmlbdy(klev,nbtt ) < rttlim)  roberr  =  rherr3 * c100r
        IF (klev == klevsfc)              roberr  =  rherr1 * c100r *4._wp
        zmlbdy(klev,nbtqer) = roberr
      ENDIF
!     PRINT '("nbtqer ",A,2I3,2I5,F8.0,2F8.2,I5)'                              &
!           , yomlhd(nmlob), ilevsfc, klevsfc, klev, mzmlbd(klev,nbtlid) &
!           , zmlbdy(klev,nbtp) , zmlbdy(klev,nbtqer)              &
!           , zmlbdy(klev,nbtrh), mzmlbd(klev,nbterr)
    ENDDO

!-------------------------------------------------------------------------------
! Section 5: Create surface report from surface-level TEMP / TOWER observations
!-------------------------------------------------------------------------------
!   PRINT *,'zqq1 ', yomlhd(nmlob), klevsfc, ilevsfc, ksurfob(irpl)            &
!                  , mzmlbd(klevsfc,nbterr)

    IF ((ksurfob(irpl) >= 1) .AND. ((ilevsfc > 0) .OR. (ltower))) THEN
      nsgob = ksurfob(irpl)

! copy surface-level data from multi-level report into single-level ODR
! ---------------------------------------------------------------------

      IF (klevsfc > 0) THEN
        osgbdy (nsgob,nbsu  ) = zmlbdy(klevsfc,nbtu  )
        osgbdy (nsgob,nbsv  ) = zmlbdy(klevsfc,nbtv  )
        osgbdy (nsgob,nbst  ) = zmlbdy(klevsfc,nbtt  )
        osgbdy (nsgob,nbsrh ) = zmlbdy(klevsfc,nbtrh )
        osgbdy (nsgob,nbsuer) = zmlbdy(klevsfc,nbtuer)
        osgbdy (nsgob,nbster) = zmlbdy(klevsfc,nbtter)
        osgbdy (nsgob,nbsqer) = zmlbdy(klevsfc,nbtqer)
        osgbdy (nsgob,nbsdrh) = zmlbdy(klevsfc,nbtdrh)
!       mosgbd (nsgob,nbslsg) = mzmlbd(klevsfc,nbtlsg)
        mosgbd (nsgob,nbslid) = mzmlbd(klevsfc,nbtlid)
        mosgbd (nsgob,nbserr) = mzmlbd(klevsfc,nbterr)
        mosgbd (nsgob,nbsflg) = mzmlbd(klevsfc,nbtflg)  ! ???
        mosgbd (nsgob,nbsqcf) = 0
        osgbdy (nsgob,nbsvip) = r_ps(iob,job) 
      ELSE
        mosgbd (nsgob,nbslid) = 0
        mosgbd (nsgob,nbserr) = 0
        mosgbd (nsgob,nbsflg) = 0
        mosgbd (nsgob,nbsqcf) = 0
      ENDIF
! (near-)surface pressure from TEMP (all code types except towers)
      IF (.NOT. ltower) THEN
        osgbdy (nsgob,nbsp  ) = zmlbdy(klevsfc,nbtp  )
        osgbdy (nsgob,nbsplv) = zmlbdy(klevsfc,nbtplv)
        osgbdy (nsgob,nbsz  ) = zmlbdy(klevsfc,nbtz  )
        osgbdy (nsgob,nbszer) = zmlbdy(klevsfc,nbtzer)
        !   cancel (TEMP) surface pressure obs if converted from reported height
        !                                   or if (barometric) height is missing
        IF (     (ABS( kproclv(ilevsfc) ) /= 2)                                &
            .OR. (osgbdy(nsgob,nbsz) < rmdich)) THEN
          osgbdy (nsgob,nbsp  ) = rmdi
          osgbdy (nsgob,nbsz  ) = rmdi
          osgbdy (nsgob,nbszer) = rmdi
          CALL MVBITS( 0 , 0 , nvfaoc , mosgbd(nsgob,nbsflg) , nvfzbp )
        ENDIF
        !   'ps' from a TEMP derived surface-level report is never used actively
        !   because a 'ps' increment is derived from the multi-level report
        mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrz )
        IF (osgbdy(nsgob,nbsp) > rmdich) THEN
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvfzbp+nvfbps(5))
          IF (osgbdy(nsgob,nbszer) < epsy)  osgbdy (nsgob,nbszer) = oezsynp(1)
          osgbdy (nsgob,nbszer) =   osgbdy(nsgob,nbszer)                       &
                                  * osgbdy(nsgob,nbsp) *r_g /(r_d * tmelt)
        ENDIF
      ENDIF
! surface pressure from (header part of ICOS) tower
!   (height + pressure have already been assigned)
      IF (ltower) THEN
        mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrz )
        CALL MVBITS( 0 , 0 , nvfaoc , mosgbd(nsgob,nbsflg) , nvfzbp )
        IF (      (osgbdy(nsgob,nbsp) >  50000._wp)                            &
            .AND. (osgbdy(nsgob,nbsp) < 110000._wp)                            &
            .AND. (osgbdy(nsgob,nbsz) <   5000._wp)                            &
            .AND. (osgbdy(nsgob,nbsz) >   -400._wp)) THEN
          zzkeml = c05* (r_hhl(iob,job,ke+1)+r_hhl(iob,job,ke))
          fisdzz = MAX( zzkeml - osgbdy(nsgob,nbsz) , 0._wp )
          IF ((osgbdy(nsgob,nbsz) > altopsu(2)) .OR. (fisdzz > doromx(2))) THEN
            neventd (nepalt,icma) = neventd (nepalt,icma) + iactr
            mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg)                &
                                         , nvfzbp+nvfbps(4))
          ELSE
            mosgbd (nsgob,nbserr) = IBSET( mosgbd(nsgob,nbserr) , nvrz )
!           mosgbd (nsgob,nbslid) = 1
          ENDIF
          !   ps obs error estimate as for Synop
          osgbdy (nsgob,nbszer) =  (oezsynp(1) + 0.04_wp* fisdzz)              &
                                  * osgbdy(nsgob,nbsp) *r_g /(r_d * tmelt)
        ELSE
          osgbdy (nsgob,nbsp) = rmdi
          osgbdy (nsgob,nbsz) = rmdi
        ENDIF
!     PRINT *,'klevsfc2 ', yomlhd(nmlob), klevsfc, zmlbdy(klevsfc,nbtp)
      ENDIF

! re-set obs type + status of TEMP or tower surface report to Synop, if required
! ------------------------------------------------------------------------------
      IF (     (ltower .AND. ((ABS(icdt_tws) == 1) .OR. (ABS(icdt_tws) == 2))) &
          .OR. ((kobtyp == ntemp) .AND. (icdt_rss >= 1))) THEN
        mosghd (nsgob,nhobtp) = nsynop
        IF (kobtyp == ntemp)                    mosghd (nsgob,nhcode) = nsytmp
        IF (ltower)                             mosghd (nsgob,nhcode) = natscd
        IF (ltower .AND. (ABS(icdt_tws) == 2))  mosghd (nsgob,nhcode) = nsytow
        IF (     ((kobtyp == ntemp) .AND. (icdt_rss >= 2))                     &
            .OR. ((ltower)          .AND. (icdt_tws >  0))) THEN
          !    remove obstype flag for TEMP resp. tower
          mosghd (nsgob,nhflag) = IBCLR( mosghd(nsgob,nhflag), FL_OBSTYPE )
          mosghd (nsgob,nhschr) = IBCLR( mosghd(nsgob,nhschr), nvexbp )
          IF (mosghd(nsgob,nhflag) == 0)  mosghd (nsgob,nhpass) = 0
        ENDIF
        !    set obstype flag for synop
        icma2   = i_cma ( mosghd(nsgob,nhobtp), mosghd(nsgob,nhcode) )
!                 =====
        IF (      (osghed(nsgob,nhjlat) <= cma(icma2)%exnlat+epsy)             &
            .AND. (osghed(nsgob,nhjlat) >= cma(icma2)%exslat-epsy)             &
            .AND. (osghed(nsgob,nhilon) >= cma(icma2)%exwlon-epsy)             &
            .AND. (osghed(nsgob,nhilon) <= cma(icma2)%exelon+epsy)) THEN
          IF ((mosghd(nsgob,nhpass) == 0) .AND. (mosgbd(nsgob,nbserr) > 0))    &
            neventr (neobct,icma2) = neventr(neobct,icma2) + 1
          mosghd (nsgob,nhpass) = 2
          mosghd (nsgob,nhschr) = IBSET( mosghd(nsgob,nhschr), nvexbp )
          mosghd (nsgob,nhflag) = IBSET( mosghd(nsgob,nhflag), FL_OBSTYPE )
        ENDIF
      ENDIF

! re-set 'nbslid' to allow for identifying tower-/TEMP-derived surface reports
! ----------------------------------------------------------------------------
      mosgbd (nsgob,nbslid) = insert( 0, 1, nvlidp(10) )

! set observation to passive to corresponding station selection criteria fail
! ---------------------------------------------------------------------------

! get observation height and difference to model orography
      zzaltsy = osghed(nsgob,nhalt)
      IF (ilevsfc > 0) THEN
        IF (zzz(ilevsfc) > rmdich)  zzaltsy = zzz(ilevsfc)
      ENDIF
!     fisd    = osghed(nsgob,nhsurf) - osghed(nsgob,nhalt)
      fisd    = osghed(nsgob,nhalt) - osghed(nsgob,nhsurf)
      fisdtt  = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd )) * fisd
      fisdrh  = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd )) * fisd
      fisduv  = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd )) * fisd
      lseaobs = BTEST( mosghd(nsgob,nhschr),nvsebp )

! 2-m temperature
      IF (osgbdy(nsgob,nbst) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbst) = MAX( 100._wp,  osgbdy(nsgob,nbst)     &
                                                      - 1._wp*fbogus )
        lseaonly = (ABS(altopsu(3)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(3)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisdtt > doromx(3))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvrt ))                              &
            neventd(netalt,icma) = neventd(netalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrt )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvftbp+nvfbps(4))
        ENDIF
      ENDIF

! 2-m relative humidity
      IF (osgbdy(nsgob,nbsrh) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbsrh) = MAX( MIN( c1 ,  osgbdy(nsgob,nbsrh)  &
                                                        + 0.1_wp*fbogus )      &
                                             , 0.01_wp )
        lseaonly = (ABS(altopsu(4)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(4)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisdrh > doromx(4))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvrq ))                              &
            neventd(neqalt,icma) = neventd(neqalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvrq )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvfqbp+nvfbps(4))
        ENDIF
      ENDIF

! 10-m horizontal wind
      IF (osgbdy(nsgob,nbsu) > rmdich) THEN
        IF (lbogus) osgbdy(nsgob,nbsu) = osgbdy(nsgob,nbsu) - 1._wp*fbogus
        IF (lbogus) osgbdy(nsgob,nbsv) = osgbdy(nsgob,nbsv) + 1._wp*fbogus
        lseaonly = (ABS(altopsu(1)) < epsy)
        IF (     ((.NOT. lseaonly) .AND. (zzaltsy > altopsu(1)))               &
            .OR. ((      lseaonly) .AND. (.NOT. lseaobs))                      &
            .OR. (fisduv > doromx(1))) THEN
          IF (BTEST( mosgbd(nsgob,nbserr),nvru ))                              &
            neventd(nevalt,icma) = neventd(nevalt,icma) + iactr
          mosgbd (nsgob,nbserr) = IBCLR( mosgbd(nsgob,nbserr) , nvru )
          mosgbd (nsgob,nbsflg) = IBSET( mosgbd(nsgob,nbsflg), nvfubp+nvfbps(4))
        ENDIF
      ENDIF

! set screen-level observations in multi-level report to passive
! (note: the obs operator for multi-level reports in COSMO
!        inter-/extra-polates between model levels and never
!        uses T2m etc.)
! --------------------------------------------------------------

      IF (klevsfc > 0) THEN
        IF (zmlbdy(klevsfc,nbtt) > rmdich) THEN
          mzmlbd (klevsfc,nbterr) = IBCLR( mzmlbd(klevsfc,nbterr), nvrt )
          mzmlbd (klevsfc,nbtflg) = IBSET( mzmlbd(klevsfc,nbtflg)              &
                                                            , nvftbp+nvfbps(4) )
!         CALL MVBITS( 1,0, nvfboc(4),mzmlbd(klevsfc,nbtflg), nvftbp+nvfbps(4) )
        ENDIF
        IF (zmlbdy(klevsfc,nbtrh) > rmdich) THEN
          mzmlbd (klevsfc,nbterr) = IBCLR( mzmlbd(klevsfc,nbterr), nvrq )
          mzmlbd (klevsfc,nbtflg) = IBSET( mzmlbd(klevsfc,nbtflg)              &
                                                            , nvfqbp+nvfbps(4) )
        ENDIF
        IF (zmlbdy(klevsfc,nbtu) > rmdich) THEN
          mzmlbd (klevsfc,nbterr) = IBCLR( mzmlbd(klevsfc,nbterr), nvru )
          mzmlbd (klevsfc,nbtflg) = IBSET( mzmlbd(klevsfc,nbtflg)              &
                                                            , nvfubp+nvfbps(4) )
        ENDIF
      ENDIF

! avoid active status flags for empty reports
! -------------------------------------------

    ELSEIF (ksurfob(irpl) >= 1) THEN
      nsgob = ksurfob(irpl)
      mosghd (nsgob,nhpass) = -1
      mosgbd (nsgob,nbserr) =  0
    ENDIF


!-------------------------------------------------------------------------------
! Section 6: Discard levels which do not contain any data
!            except for pressure derived from height
!-------------------------------------------------------------------------------
! note that 'isrtlvp' must not be used any more after the following loop !

    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (      (zmlbdy(klev,nbtu  ) <= rmdich)                                &
          .AND. (zmlbdy(klev,nbtt  ) <= rmdich)                                &
          .AND. (zmlbdy(klev,nbtrh ) <= rmdich)                                &
!         .AND. (zmlbdy(klev,nbtw  ) <= rmdich)                                &
          .AND. (     (zmlbdy(klev,nbtz  ) <= rmdich)                          &
                 .OR. (kproclv(ilev) <= 1)                                     &
                 .OR. (kcdtyp == ntdesc)))   kproclv (ilev) = -9
    ENDDO

    !  also discard levels below surface level
    !  (in fact, this has already been done further above)
    DO klev = 1 , klevsfc - 1
      ilev = isrtlvp(klev)
      kproclv (ilev) = -9
    ENDDO
    !   from now on, 'klevsfc' is either 0 (no surface level)
    !                                 or 1 (level index of surface level in ODR)
    klevsfc  =  MIN( klevsfc , 1 )

    nrejlv = 0
    DO klev = 1 , nlvp
      ilev = isrtlvp(klev)
      IF (kproclv(ilev) <= -7) THEN
        nrejlv = nrejlv + 1
      ELSEIF (nrejlv >= 1) THEN
        zmlbdy (klev-nrejlv,:) = zmlbdy(klev,:)
        mzmlbd (klev-nrejlv,:) = mzmlbd(klev,:)
      ENDIF
    ENDDO
    nlvp = nlvp - nrejlv
    momlhd (nmlob,nhnlev) = nlvp

!-------------------------------------------------------------------------------
! Section 7: Apply superobbing for high-resolution profiles
!-------------------------------------------------------------------------------
!            A weighted average of all reported levels encompassing the
!            superobbing layer is done;
!            if at least (only) one of the reported layers has missing value
!            then the whole superobs has missing value as well.
!            Note: In this whole section, use index 'nbtplv' instead of 'nbtp'
!                  so that the vertical level information (level type 'pressure'
!                  for reported pressure (indlev==2), 'height' for reported
!                  height (indlev == 1)) of the superobs will be independent
!                  from the model state and therefore identical for all
!                  members of an ensemble.

    ! 1st criterion for superobbing: profile contains more levels than a
    !               threshold which depends on namelist parameter 'av_reso'
    ! ---------------------------------------------------------------------

    nlimsup = NINT( av_reso* ke )
    IF (nlvp >= 2) THEN
      ! superobb. threshold: av_reso*ke levels scaled by range of profile + 2
      nlimsup  =  NINT( av_reso* ke *(zmlbdy(1,nbtplv) - zmlbdy(nlvp,nbtplv))  &
                                    *1E-5_wp )  +  2
      IF ((zmlbdy(1,nbtz) > rmdich) .AND. (zmlbdy(nlvp,nbtz) > rmdich)) THEN
        !   (kbot - ktop) is the (approx.) number of model layers
        !   covered by the observed profile
        ktop = 1
        kbot = 1
        DO kk = 1 , ke
          IF (r_hhl(iob,job,kk) > zmlbdy(nlvp,nbtz))  ktop = kk
          IF (r_hhl(iob,job,kk) > zmlbdy(1   ,nbtz))  kbot = kk
        ENDDO
        nlimsup  =  NINT( av_reso *(kbot - ktop) )  +  2
      ENDIF
      !   for practical reasons, 'nlimsup' may be reduced depending on ODR size
      nlimsup  =  MIN( nlimsup , 2*MAX( maxmlv , nref ) )
!     PRINT '(A4,A8,2I4,F6.2,2I5,2F8.0,2I5)', "SUP ", yomlhd(nmlob)            &
!           , kobtyp, kcdtyp, omlhed(nmlob,nhtime), nlvp, nlimsup              &
!           , zmlbdy(1,nbtplv), zmlbdy(nlvp,nbtplv), kbot, ktop
    ENDIF
    ! 2nd criterion for superobbing: 
    !     - report is not TEMP part B (which contains a lot of missing values
    !       (e.g. T at significant wind levels) and would result in a lot of
    !        missing superobbed data)
    !       (for other TEMP parts superobbing should not be done either), and
    !     - report is not TOWER
    ! --------------------------------------------------------------------------
    ltempb = .FALSE.
    IF (kobtyp == ntemp) THEN
      ! using Kennzahl of DWD data base
      IF ((momlhd(nmlob,nhkz) > 500) .AND. (momlhd(nmlob,nhkz) < 800)) THEN
        ltempb = .TRUE.
      ! otherwise: the fraction of missing data (below the 300 hPa level)
      !            should be < 20 % for each variabele
      ELSE
        ictot        = 0
        icmiss (1:3) = 0
        DO klev = 1 , nlvp
          IF (zmlbdy(klev,nbtplv) > 30000._wp) THEN
                                              ictot     = ictot     + 1
            IF (zmlbdy(klev,nbtu ) < rmdich)  icmiss(1) = icmiss(1) + 1
            IF (zmlbdy(klev,nbtt ) < rmdich)  icmiss(2) = icmiss(2) + 1
            IF (zmlbdy(klev,nbtrh) < rmdich)  icmiss(3) = icmiss(3) + 1
          ENDIF
        ENDDO
        IF (ictot > 0)  ltempb = ( REAL( MAXVAL( icmiss(1:3) ),wp )            &
                                  /REAL(         ictot        ,wp ) > 0.2_wp)
      ENDIF
    ENDIF
    nlint = 0
!   IF (nlvp > MAX( 100 , MIN( 300 , maxmlv+50 ) )) THEN
    IF ((nlvp > nlimsup) .AND. (.NOT. ltempb) .AND. (.NOT. ltower)) THEN
      !   re-set code type if superobbing is applied to land radiosonde profile
      IF ((kcdtyp == nldtcd) .AND. (nlvp > nlimsup)) THEN
        momlhd(nmlob,nhcode)  =  nbtemp
        kcdtyp                =  momlhd(nmlob,nhcode)
      ENDIF

      ! derive list of interpolation levels 'p_int' for current report
      ! --------------------------------------------------------------
      !   get pressure of bottom and of top level in report
      !   (note: levels without obs have already been discarded above)
      p_top = zmlbdy(nlvp,nbtplv)
      p_bot = zmlbdy(1   ,nbtplv)
      IF (klevsfc >= 1)  p_bot = MIN( zmlbdy(klevsfc+1,nbtplv)                 &
                                    , zmlbdy(klevsfc  ,nbtplv) - dp_surf )
      !   get interpolation levels to be used for current report
      ktop  =  MINLOC( ABS( pref - p_top ), dim=1 )
      kbot  =  MINLOC( ABS( pref - p_bot ), dim=1 )
      nlint =  ktop - kbot + 1
    ENDIF
    IF (nlint >= 1) THEN
      ALLOCATE ( p_int (MAX(nlint,2)) , STAT=istat )
      IF (nlint >= 2) THEN
        p_int (1:nlint)  =  pref(kbot:ktop)
        !   adjust interpolation top / bottom levels to avoid extrapolation
        !   (in this way, the top and lowermost superobs layer is at least
        !    as large as half of the regular superobs layer given by pref)
        p_int (1)        =  MIN( p_int(1)    , p_bot )
        p_int (nlint)    =  MAX( p_int(nlint), p_top )
      ELSE
        !   for very short obs profiles (which do not cover half of any of the
        !   reference superobs layers), produce 1 superobs over all obs levels
        nlint            =  2
        p_int (1)        =  p_bot
        p_int (nlint)    =  p_top
      ENDIF
!     PRINT '(A4,A8,2I4,F6.2,2I4,2F8.0,2I5)', "SUR ", yomlhd(nmlob)            &
!           , kobtyp, kcdtyp, omlhed(nmlob,nhtime), nlvp, nlint                &
!           , p_int(1), p_int(MAX(nlint,1)), kbot, ktop

      ! compute (nlint-1) superobs (as averages between interpolation levels)
      ! --------------------------
      ALLOCATE ( supbdy (nlint, mxrbdy) , STAT=istat )
      ALLOCATE ( msupbd (nlint, mxrbdf) , STAT=istat )
      supbdy (:,:) = rmdi
      !   set bit (flag) in superobs as soon as bit (flag) is set in one of the
      !   original levels --> use 'IOR'   (for 'nbtflg', 'nbtqcf')
      lflg0 (:)      = .TRUE.
      !   set superobs to passive bs as soon as one original obs is passive
      !                   --> use 'IAND'  (for 'nbterr', 'nbtlsg', 'nbtlid')
      lflg0 (nbterr) = .FALSE.
      lflg0 (nbtlsg) = .FALSE.
      lflg0 (nbtlid) = .FALSE.
      !   (at most) nlint-1 superobs are derived in obs_average_layers

      CALL obs_average_layers ( nlint, p_int, mxrbdy, mxrbdf, nlvp-klevsfc     &
                              , nbtp, nbtplv, lflg0, rmdi, rmdich              &
                              , mzmlbd(klevsfc+1:nlvp,:)                       &
                              , zmlbdy(klevsfc+1:nlvp,:) , msupbd, supbdy )
    ! =======================
      DEALLOCATE ( p_int , STAT=istat )

      ! discard levels which do not contain any data (as in previous section)
      nrejlv = 0
!     DO klev = 1 , nlint-1
      DO klev = 1 , nlint
        IF (      (supbdy(klev,nbtu  ) <= rmdich)                              &
            .AND. (supbdy(klev,nbtt  ) <= rmdich)                              &
            .AND. (supbdy(klev,nbtrh ) <= rmdich)                              &
            .AND. (supbdy(klev,nbtz  ) <= rmdich)                              &
            .AND. (supbdy(klev,nbtw  ) <= rmdich)) THEN
          nrejlv = nrejlv + 1
        ELSEIF (nrejlv >= 1) THEN
          supbdy (klev-nrejlv,:) = supbdy(klev,:)
          msupbd (klev-nrejlv,:) = msupbd(klev,:)
        ENDIF
      ENDDO
      nlvsup = nlint - 1 - nrejlv
      ! adjust flags for missing superobs
      DO klev = 1 , nlvsup
        IF (MIN( supbdy(klev,nbtu),supbdy(klev,nbtv) ) <= rmdich)              &
          msupbd (klev,nbterr) = IBCLR( msupbd(klev,nbterr), nvru )
        IF (supbdy(klev,nbtt ) <= rmdich)                                      &
          msupbd (klev,nbterr) = IBCLR( msupbd(klev,nbterr), nvrt )
        IF (supbdy(klev,nbtrh) <= rmdich)                                      &
          msupbd (klev,nbterr) = IBCLR( msupbd(klev,nbterr), nvrq )
        IF (supbdy(klev,nbtz ) <= rmdich)                                      &
          msupbd (klev,nbterr) = IBCLR( msupbd(klev,nbterr), nvrz )
      ENDDO
      !   set level significance to superobs
      DO klev = 1 , nlvsup
        msupbd (klev,nbtlsg) = IBSET( 0 , LS_SUPEROBS )
        msupbd (klev,nbtlid) = IBSET( 0 , nvlidp(3)   )
      ENDDO

      ! if the original levels fit into the ODR (level dimension)
      ! and there is space in the ODR (report dimension) then
      ! append (copy) the original report, flagged as 'merged',
      ! to the end of the ODR
      ! -------------------------------------------------------
      IF (nlvp <= maxmlv) THEN
        IF (ntotml + nrepml + nrepadd + 1 > maxmll) THEN
          nexceed (1) = nexceed(1) + 1
        ELSE
          nrepadd = nrepadd + 1
          nmladd  = ntotml + nrepml + nrepadd
          yomlhd (nmladd         ) = yomlhd(nmlob         )
          omlhed (nmladd       ,:) = omlhed(nmlob       ,:)
          momlhd (nmladd       ,:) = momlhd(nmlob       ,:)
          omlbdy (nmladd,1:nlvp,:) = zmlbdy(      1:nlvp,:)
          momlbd (nmladd,1:nlvp,:) = mzmlbd(      1:nlvp,:)
          ! flag report as 'merged' (set active state to merged)
          momlhd (nmladd,nhpass) = MAX  ( momlhd(nmladd,nhpass) , 1 )
          momlhd (nmladd,nhflag) = IBSET( momlhd(nmladd,nhflag) , FL_MERGE )

          CALL obs_rm_duplicate_levels ( nlvp, mxrbdy, mxrbdf, nbtplv, rmdich  &
                                       , omlbdy(nmladd,1:nlvp,:)               &
                                       , momlbd(nmladd,1:nlvp,:), nrejlv )
        ! ============================
          momlhd (nmladd,nhnlev) = nlvp - nrejlv

          ! adjust report state if no data in report, without updating counters

          CALL check_data_in_report ( nmladd, -1, iactr )
        ! =========================
        ENDIF
      ENDIF

      ! copy superobbed levels into 'zmlbdy', 'mzmlbd'
      ! ----------------------------------------------
      zmlbdy (klevsfc+1:klevsfc+nlvsup,:) = supbdy(1:nlvsup,:)
      mzmlbd (klevsfc+1:klevsfc+nlvsup,:) = msupbd(1:nlvsup,:)

      ! increase correction indicator to impose preference in redundancy check
      momlhd (nmlob,nhcorr) = momlhd(nmlob,nhcorr) + 10

      DEALLOCATE ( supbdy  , STAT=istat )
      DEALLOCATE ( msupbd  , STAT=istat )

      ! re-set number of processed levels (report with super-obs)
      nlvp  =  nlvsup + klevsfc
    ENDIF

!-------------------------------------------------------------------------------
! Section 8: Remove duplicate levels and
!            fill temporary observation body into ODR (obs data record)
!-------------------------------------------------------------------------------

    ! first remove duplicate levels
    ! -----------------------------
    !   use 'nbtplv' instead of 'nbtp' to get independence from model state

    CALL obs_rm_duplicate_levels ( nlvp, mxrbdy, mxrbdf, nbtplv, rmdich        &
                                 , zmlbdy(1:nlvp,:), mzmlbd(1:nlvp,:), nrejlv )
  ! ============================
!   momlhd (nmlob,nhnlev) = nlvp - nrejlv
    nlvp                  = nlvp - nrejlv

    ! check if number of levels exceeds ODR size, and adjust number of levels
    ! -----------------------------------------------------------------------
    IF (nlvp > maxmlv) THEN
      !   insufficient ODR size is the only report and data events, which is
      !   updated even for reports set to passive previously
      neventd (nelodr,icma) = neventd(nelodr,icma) + nlvp - maxmlv
      ilen = 3 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt4
        outbuf(nacout+3) = nlvp
        DO icl = 1 , istrej
          outbuf(nacout+3+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF
    ENDIF
    nlvp  =  MIN( nlvp , maxmlv )

    !   in order to accommodate TEMP (or PILOT) Parts A, B, C, D in 1 report
    !   (in obs_cdf_redundancy) limit each part (mainly part B !) to 'maxmlv-15'
    !   (note: this makes use of the DWD-specific element 'nhkz')
    IF ((momlhd(nmlob,nhkz) > 500) .AND. (momlhd(nmlob,nhkz) < 800))           &
      nlvp  =  MIN( nlvp , maxmlv - 15 )

    ! adjust the number of vertical levels in the ODR header
    ! to those levels which are filled in the ODR body
    ! ------------------------------------------------------
    momlhd (nmlob,nhnlev) = nlvp

    ! copy temporary array into ODR
    ! -----------------------------
    omlbdy (nmlob,1:nlvp,:) = zmlbdy(1:nlvp,:)
    momlbd (nmlob,1:nlvp,:) = mzmlbd(1:nlvp,:)

!-------------------------------------------------------------------------------
! Section 9: Check existence of data in multi-level report and surface report
!-------------------------------------------------------------------------------

! multi-level report
! ------------------

    CALL check_data_in_report ( nmlob, icma, iactr )
  ! =========================

!   PRINT *,'noctrj_m1 ', icma, irpl, cma(icma)%cnt_ps, cma(icma)%cnt_rj       &
!                       , iactr, nzaexi, momlhd(nmlob,nhpass), nlvp

! surface report
! --------------

    IF (ksurfob(irpl) >= 1) THEN
      nsgob = ksurfob(irpl)
      iactx = -1
      ipasx = -1
      IF (BTEST( mosgbd(nsgob,nbserr),nvru )) iactx (1) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrz )) iactx (2) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrt )) iactx (3) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrq )) iactx (4) = 1
      IF (osgbdy(nsgob,nbsu  ) > rmdich) ipasx (1) = 1
      IF (osgbdy(nsgob,nbsp  ) > rmdich) ipasx (2) = 1
      IF (osgbdy(nsgob,nbst  ) > rmdich) ipasx (3) = 1
      IF (osgbdy(nsgob,nbsrh ) > rmdich) ipasx (4) = 1
      iactx (5)             = MAX( iactx(1), iactx(2), iactx(3), iactx(4) )
      ipasx (5)             = MAX( ipasx(1), ipasx(2), ipasx(3), ipasx(4) )
      nzaexi                = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
      IF (      (nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)              &
          .AND. (ilevsfc > 0)) THEN
        mosghd (nsgob,nhpass)  =  2
        mosghd (nsgob,nhflag)  =  IBSET( mosghd(nsgob,nhflag) , FL_NO_OBS )
      ENDIF

      ! flag reports which are to be discarded completely
      IF (     (nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))       &
          .OR. ((ilevsfc == 0) .AND. (.NOT. ltower))                           &
          .OR. ((kobtyp == ntemp) .AND. (icdt_rss <= -1)))                     &
        mosghd (nsgob,nhpass) = -1

      ! update counters for obs type Synop (surface report from TEMP/tower)
      IF ((mosghd(nsgob,nhobtp) == nsynop).AND.(mosghd(nsgob,nhpass) >= 0)) THEN
        icma2   = i_cma ( mosghd(nsgob,nhobtp), mosghd(nsgob,nhcode) )
!                 =====
        cma(icma2)%cnt_pr = cma(icma2)%cnt_pr + 1
        IF (mosghd(nsgob,nhpass) == 0) cma(icma2)%cnt_ac = cma(icma2)%cnt_ac + 1
        IF (mosghd(nsgob,nhpass) == 2) cma(icma2)%cnt_ps = cma(icma2)%cnt_ps + 1
      ENDIF
    ENDIF

  ENDDO  ! huge loop over report
! ~~~~~

  nodrnew(1)  =  nodrnew(1) + nrepadd

! this is done after obs_cdf_redundancy
! TODO : delete reports with (moxxhd(nxxob,nhpass) == -1)
! CALL obs_del_old_reports
! ========================

! clean up
! --------
  DEALLOCATE ( zmlbdy  , STAT=istat )
  DEALLOCATE ( mzmlbd  , STAT=istat )

  IF (nrepml >= 1) DEALLOCATE ( blk_loc , STAT=istat )

  DEALLOCATE ( kblk    , STAT=istat )
  DEALLOCATE ( nflgx   , STAT=istat )
  DEALLOCATE ( ioberr  , STAT=istat )
! DEALLOCATE ( ztv     , STAT=istat )
  DEALLOCATE ( zttl    , STAT=istat )
  DEALLOCATE ( ztdl    , STAT=istat )
  DEALLOCATE ( zqxl    , STAT=istat )
! DEALLOCATE ( ztqx    , STAT=istat )
  DEALLOCATE ( zrhw    , STAT=istat )
  DEALLOCATE ( zrhc    , STAT=istat )
  DEALLOCATE ( zuu     , STAT=istat )
  DEALLOCATE ( zvv     , STAT=istat )
  DEALLOCATE ( zrlat   , STAT=istat )
  DEALLOCATE ( zrlon   , STAT=istat )
  DEALLOCATE ( zoberr  , STAT=istat )
  DEALLOCATE ( zstd    , STAT=istat )
  DEALLOCATE ( lneedt  , STAT=istat )
  DEALLOCATE ( lneedq  , STAT=istat )
  DEALLOCATE ( lneedv  , STAT=istat )

  DEALLOCATE ( zpp     , STAT=istat )
  DEALLOCATE ( zppl    , STAT=istat )
  DEALLOCATE ( zdlat   , STAT=istat )
  DEALLOCATE ( zdlon   , STAT=istat )
  DEALLOCATE ( ztt     , STAT=istat )
  DEALLOCATE ( ztd     , STAT=istat )
  DEALLOCATE ( zqx     , STAT=istat )
  DEALLOCATE ( zrh     , STAT=istat )
  DEALLOCATE ( zff     , STAT=istat )
  DEALLOCATE ( zzz     , STAT=istat )
  DEALLOCATE ( zdd     , STAT=istat )
  DEALLOCATE ( zww     , STAT=istat )
  DEALLOCATE ( zsinor  , STAT=istat )
  DEALLOCATE ( zstdff  , STAT=istat )
  DEALLOCATE ( zazi    , STAT=istat )
  DEALLOCATE ( izdt    , STAT=istat )
  DEALLOCATE ( izlv    , STAT=istat )
  DEALLOCATE ( iroll   , STAT=istat )
  DEALLOCATE ( iqci    , STAT=istat )
  DEALLOCATE ( isls    , STAT=istat )
  DEALLOCATE ( isle    , STAT=istat )
  DEALLOCATE ( llevsfc , STAT=istat )
  DEALLOCATE ( kproclv , STAT=istat )
  DEALLOCATE ( isrtlvp , STAT=istat )

! IF ((ltower) .AND. (ydate_ref(11:14) == '0000')) THEN
  IF ((ltower) .AND. (nrepml > 0))  CALL obs_radsum_1h ( nrepml , ksurfob )
                                  ! ==================

  DEALLOCATE ( ksurfob , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_store_multilev
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_store_multilev


!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for checking existence of data
!-------------------------------------------------------------------------------

SUBROUTINE obs_radsum_1h ( nrepml , ksurfob )

!-------------------------------------------------------------------------------
! Description:
!   This subroutine produces 1-hour sums of observed surface global radiation
!   (valid a full hour) from sub-hourly values. 
!   Remaining sums over periods other than 1 hour are removed.
!
! Method: Straightforward
! Written by: Christoph Schraff, 08.05.23
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    nrepml           ,& ! number of multi-level reports
    ksurfob (nrepml+1)  ! surface report indicator and report index

! Local scalars or fixed arrays:
! -----------------------------

  INTEGER                  ::  &
    ixrob (nrepml+1) ,& ! index list of obs with sub-hourly radiation sums
    ixro1h(nrepml+1) ,& ! index list of obs with radiation valid at full hours
    irpl             ,& ! loop index over multi-level reports
    nob              ,& ! number of reports with sub-hourly radiation sums
    nos              ,& ! number of reports with radiation valid at full hours
    iob              ,& ! loop index   over 'nob'
    ios    , io2     ,& ! loop indices over 'nos'
    nsgob  , nsgo2   ,& ! ODR report index for single-level reports
    minob  , mino2   ,& ! obs time in [minutes]
    mino2s           ,& ! start of period [min] of observed radiation
    minlist (60)        ! indicates which minutes are already covered to sum
                        !   up a 1-hour radiation sum

  LOGICAL                  ::  &
    lnewos              ! new report with radiation valid at full hours
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_radsum_1h
!-------------------------------------------------------------------------------

  !   make an (index) list of obs reports with sub-hourly radiation sums
  !     (and 60 min divisible by obs period)
  nob = 0
  DO irpl = 1 , nrepml
    IF (ksurfob(irpl) >= 1) THEN
      nsgob = ksurfob(irpl)
      IF (osgbdy(nsgob,nbsrad) > rmdich) THEN
        IF (     (ABS( mosghd(nsgob,nhdt) ) > 60)                              &
            .OR. (ABS( mosghd(nsgob,nhdt) ) == 0)                              &
            .OR. (MOD( 60 , MAX(1,ABS( mosghd(nsgob,nhdt) )) ) /= 0)) THEN
          osgbdy (nsgob,nbsrad) = rmdi
        ELSEIF (ABS( mosghd(nsgob,nhdt) ) < 60) THEN
          nob = nob + 1
          ixrob (nob) = nsgob
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  !   make list of obs valid at full hours and with unique (sta ID & obs time)
  !   (this is the list of obs which can be complemented to 1-hrly radiation)
  !     (use entry 'nhhrmn' here as ref time does not need to be at full hour)
  ixro1h = 1
  nos    = 0
  DO iob = 1 , nob
    nsgob = ixrob(iob)
    IF (MOD( mosghd(nsgob,nhhrmn), 100 ) == 0) THEN
      lnewos = .TRUE.
      DO ios = 1 , nos
        IF (      (yosghd(nsgob)        == yosghd(ixro1h(ios)))                &
            .AND. (mosghd(nsgob,nhhrmn) == mosghd(ixro1h(ios),nhhrmn)))        &
          !   if same station ID and obs time: obs redundant and not used
          lnewos = .FALSE.
      ENDDO
      IF (lnewos) THEN
        nos = nos + 1
        ixro1h (nos) = nsgob
      ENDIF
    ENDIF
  ENDDO
  DO ios = 1 , nos
    nsgob = ixro1h(ios)
    !   the end [min] of the full hour for the 1-hrly sum
    minob = NINT( osghed(nsgob,nhtime) *60._wp )
    minlist = 0
    DO iob = 1 , nob
      nsgo2  = ixrob(iob)
      mino2  = NINT( osghed(nsgo2,nhtime) *60._wp )
      mino2s = mino2 - ABS( mosghd(nsgo2,nhdt) ) + 1
      IF ((yosghd(nsgo2) == yosghd(nsgob)) .AND. (mino2  <= minob)             &
                                           .AND. (mino2s >= minob-59)) THEN
        !   period of obs 'nsgo2' lies within required 1-hr interval
        IF (MAXVAL( minlist(mino2s-minob+60:mino2-minob+60) ) == 0) THEN
          !   if this part of the 1-hr interval is not yet covered then add the
          !   radiation (unless nsgo2==nsgob for which the radiation is already
          !   'added') and indicate this part of the interval as being covered
          minlist (mino2s-minob+60:mino2-minob+60) = 1
          IF (nsgo2 /= nsgob)  osgbdy (nsgob,nbsrad) =   osgbdy(nsgob,nbsrad)  &
                                                       + osgbdy(nsgo2,nbsrad)
        ENDIF
      ENDIF
    ENDDO
    !   set the obs period to 60 min if the full hour is covered
    !     (double counting (double coverage) is excluded above)
    IF (MINVAL( minlist ) == 1)  mosghd (nsgob,nhdt) = -60
  ENDDO
  !   only keep the 1-hr sums of radiation obs (set other sums to missing value)
  DO iob = 1 , nob
    nsgob = ixrob(iob)
    IF (ABS( mosghd(nsgob,nhdt) ) /= 60)  osgbdy (nsgob,nbsrad) = rmdi
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_radsum_1h
!-------------------------------------------------------------------------------
END SUBROUTINE obs_radsum_1h

!===============================================================================
!+ Module procedure in "src_obs_cdfin_mult" for checking existence of data
!-------------------------------------------------------------------------------

SUBROUTINE check_data_in_report ( nmlrp, icma, iactr )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_mult" checks the existence
!   of active or passive data in multi-level report.
!
! Method:
!   Straightforward
!
! Initial release: Christoph Schraff, 26.04.19
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER                  , INTENT (IN)    ::  &
    nmlrp            ,& ! target ODR report index (multi-level report)
    icma             ,& ! >= 0 : pointer index for 'cma' / diagnostic arrays
                        ! == -1: do not update     'cma' / diagnostic arrays
    iactr               ! stepping for event counters

! Local scalars or fixed arrays:
! -----------------------------

  INTEGER                  ::  &
    klev             ,& ! loop index over vertical levels (in ODR)
    iactx (5)        ,& ! indicates presence of active  data (for variables)
    ipasx (5)        ,& ! indicates presence of passive data (for variables)
    nzaexi              ! -1: only passive data ; 0: no data ; 1: active data
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine check_data_in_report
!-------------------------------------------------------------------------------

  iactx = -1
  ipasx = -1
  DO klev = 1 , momlhd(nmlrp,nhnlev)
!   ilev = isrtlvp(klev)
    IF (BTEST( momlbd(nmlrp,klev,nbterr),nvru )) iactx (1) = 1
    IF (BTEST( momlbd(nmlrp,klev,nbterr),nvrz )) iactx (2) = 1
    IF (BTEST( momlbd(nmlrp,klev,nbterr),nvrt )) iactx (3) = 1
    IF (BTEST( momlbd(nmlrp,klev,nbterr),nvrq )) iactx (4) = 1
    IF (omlbdy(nmlrp,klev,nbtu  ) > rmdich) ipasx (1) = 1
    IF (omlbdy(nmlrp,klev,nbtz  ) > rmdich) ipasx (2) = 1
    IF (omlbdy(nmlrp,klev,nbtt  ) > rmdich) ipasx (3) = 1
    IF (omlbdy(nmlrp,klev,nbtrh ) > rmdich) ipasx (4) = 1
  ENDDO
  iactx (5)             = MAX( iactx(1), iactx(3), iactx(4) )
  ipasx (5)             = MAX( ipasx(1), ipasx(3), ipasx(4) )
  iactx (2)             = MAX( iactx(2), iactx(5) )
  ipasx (2)             = MAX( ipasx(2), ipasx(5) )
  momlhd (nmlrp,nhuexi) = MAX( MIN( 0 , -ipasx(1) ) , iactx(1) )
  momlhd (nmlrp,nhtexi) = MAX( MIN( 0 , -ipasx(3) ) , iactx(3) )
  momlhd (nmlrp,nhqexi) = MAX( MIN( 0 , -ipasx(4) ) , iactx(4) )
  momlhd (nmlrp,nhaexi) = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
  nzaexi                = MAX( MIN( 0 , -ipasx(2) ) , iactx(2) )

  IF ((nzaexi <= 0) .AND. (icma > -1))                                         &
    neventr (nenoda,icma) = neventr(nenoda,icma) + iactr
  IF ((nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)) THEN
    momlhd (nmlrp,nhpass)  =  2
    momlhd (nmlrp,nhflag)  =  IBSET( momlhd(nmlrp,nhflag) , FL_NO_OBS )
    IF (icma > -1)  cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
    IF (icma > -1)  cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1

! flag reports which are to be discarded completely
  ELSEIF ((nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))) THEN
    momlhd (nmlrp,nhpass)  = -1
    IF (icma > -1) THEN
      IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine check_data_in_report
!-------------------------------------------------------------------------------
END SUBROUTINE check_data_in_report

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

END MODULE src_obs_cdfin_mult
