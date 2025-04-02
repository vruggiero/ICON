!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------
#if __NEC_VERSION__ >= 30005    /* options directive supported since 3.0.5 */
!NEC$ options "-fno-move-loop-invariants"
#endif

MODULE src_obs_cdfin_sing

!-------------------------------------------------------------------------------
! Description:
!   This module performs the observation pre-processing of single-level reports
!   read from NetCDF observation input files. The current types of single-level
!   reports comprise of surface-level reports and upper-air aircraft reports,
!   however multi-level aircraft reports are also read and pre-processed here.
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
!   Known bug: In MEC, the option based on roughness length and/or Laplacian
!   =========  of orography to select 10-m wind stations is hard-set to false,
!              because in MEC, r_hhl is not a global field and r_z0 not read
!              so far!
!
! Method:
!   This module contains the following module procedures:
!    - obs_cdf_read_surface   : (called by obs_cdf_read_org)
!    - obs_cdf_read_aircraft  : (called by obs_cdf_read_org)
!    - obs_cdf_read_amv       : (called by obs_cdf_read_org)
!    - obs_cdf_store_singlev  :  called by obs_cdf_read_surface
!                                      and obs_cdf_read_aircraft
!
!   This module also contains elemental functions, formerly statement functions:
!   - ireplace     : replaces a bit pattern in an integer by another bit pattern
!   - insert       : inserts bits set in one integer into another integer word
!                    at given bit position
!
!   It uses from:
!    - src_obs_cdfin_mult:    - obs_cdf_store_multilev
!    - src_obs_cdfin_comhead: - obs_cdf_read_comhead
!                             - obs_cdf_buffer_comhead
!                             - obs_cdf_store_comhead
!    - src_obs_cdfin_blk:     - obs_cdf_whitelist_local
!                             - obs_cdf_blacklist_local
!    - src_obs_cdfin_util:    - obs_assign_sort_node
!                             - obs_cdf_distrib_reports
!                             - std_atmosphere
!                             - obs_td2rh , obs_qx2rh , obs_rhw2rh, obs_rh2td
!                             - obs_find_level
!                             - f_z2p
!    - utilities:             - uv2uvrot_vec
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
!  Initial release, extracted from module 'src_obs_proc_cdf' and adapted (e.g
!    modified routine interfaces and diagnostic arrays, 'obs_pointrs'-->'i_cma',
!    modified ODR (high cloud cover and observation status word introduced,
!    cloud group words modified)).
!  - Call of external routine 'atmhp' (libmisc) replaced by 'std_atmosphere'.
!  - To make observation input more flexible, some variables in the NetCDF
!    files are changed from mandatory to optional.
!  - Bug correction at selection of surface obs (using sign of 'fdoro').
!  - Modified criteria to reject surface pressure obs due to vertical
!    extrapolation errors.
!  - Bug correction: condition to reject AMDAR obs with reported height below
!    model orography is replaced by a pressure criterion because reported height
!    is ficticious (in fact, it is flight level, derived from pressure using
!    the ICAO standard atmosphere).
!  - Bug correction at vertical significance of individual cloud layers.
! V4_28        2013/07/12 Christoph Schraff
!  - Relaxed conditions on required variables in NetCDF obs input files ('lax':
!    e.g. variable names for main observables, variable for vertical coordinate
!    (e.g. pressure allowed also for AMDAR), roll angle not mandatory variable).
!  - Surface pressure obs error converted from height to pressure units.
!  - Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Extensions for reading / processing Mode-S aircraft observations, from NetCDF
!  files with 2 different templates: (i) converted from BUFR as received from
!  Mode-S processing centre KNMI, (ii) converted into DWD-ACARS template. (CS)
!  Setting rel. humidity to passive because of passive temperature obs is now
!  imposed only if the T obs has been set passive for other reasons than not
!  being in valid height range. (Required for KENDA) (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_3         2015-10-09 Lucio Torrisi, Christoph Schraff
!  - Adaptations to read and process AMDAR data according to WIGOS template and
!    scatterometer (METOP) data according to Eumetsat template.
!  - Additional quantities written to single-level ODR for verification
!    (model equivalent calculator) purposes: wind speed + direction, dewpoint
!    temperature, radiation, 3-hourly precip and wind gusts. Max. and min. T-2m
!    over arbitrary period instead of fixed 12 hours.
! V5_4         2016-03-10 Christoph Schraff
!  - Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  - Variables related to AOF interface removed, and code clean-up.
!  - For VOF: if height flag set, set passive bit in station characteristics.
! V5_4d        2016-12-12 Christoph Schraff
!  - Additional limits 'zlimv10' for the use of 10-m wind observations.
!  - Reading of relative humidity (erroneously) defined as Real variable in
!    some Mode-S reports.
!  - Surface pressure <1Pa is set to missing to avoid division by zero later on.
!  - Zero wind speed from all obs types (not only buoys and scatterometer)
!    flagged as bad reporting practice (for safety).
!  - Processing added for observation input files 'cdfin_nat_dwd' with national
!    surface reports from DWD containing 1-hourly precip and wind gusts.
! V5_4e        2017-03-23 Ulrich Schaettler
!  Initialized variable status to nf90_noerr
! V5_4f        2017-09-01 Christoph Schraff
!  Update from DACE V1_51: Transformation of observed wind into rotated
!  coordinate system executed only if run as part of COSMO code.
! V5_4h        2017-12-15 Christoph Schraff
!  - Adaptations (e.g. introduction of subr. 'obs_cdf_read_amv') to read,
!    process and store SATOB / AMV reports.
!  - Variable for dewpoint ('MTDNH' or 'MTDN') made optional for AMDAR (if lax).
! V5_6         2019-02-27 Christoph Schraff
!  - Satellite frequency digit added to wind computation method entry for AMV.
!  - 'NIX' changed to a non-mandatory variable for in-situ data.
!  - Comments on buoy station type (MTODB instead of MQOBL) corrected.
! V5_6b        2019-10-16 Christoph Schraff
!  - For processing steps (e.g. blacklisting, gross error checking, superobbing)
!    affecting the state (active, passive...) of a report or observation,
!    'observed pressure' (nbsp) is replaced by 'pressure level independent
!    from model state' (nbsplv) in order to render the number of reports and
!    obs independent from the model state and thus of the member of an ensemble
!    (this is required for 'online LETKF').
!    There is a single exception: Setting aircraft and AMV reports to passive
!    close to the model orography deploys the model state,
!    i.e. surface pressure.
!    Set additional report flag 'nvlzbp': pressure derived from height by
!    using model state.
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! @VERSION@    @DATE@     Christoph Schraff
!  Introduction of synop report from alternative networks for testing.
!  25.01.2022: Corrected unit conversion for relative humidity obs.
!  16.05.2022: Std. dev. of sub-grid scale model orography added as criterion
!              for 10-m wind
!  20.05.2022: Definition of 'octas' adapted to WMO Table 020011.
!  01.08.2022: Complete revision of cloud-related observation pre-processing,
!              in particular of derivation of low, mid-level, and high cloud
!              cover (amount). Derivation of cloud ceiling is added.
!              For ceiling and cloud base height, preference is given to the
!              info from the individual cloud layer group over that of the
!              general cloud group.
!  11.10.2022: Further revision of cloud-related observation pre-processing:
!              - more comprehensive use of present weather 'ww' for
!                fog (probability) indicator;
!              - improved consistency check in cloud-free case (nclct == 0);
!              - non-integer cloud cover 'clcd' (for few / scattered  cloud);
!  13.01.2023: Results from revised cloud-related diagnostics coincide with
!              those from (mainly 'cloudcover_synop_lam' in) mo_synop.f90.
!  18.04.2023: 'rhtsat' is namelist parameter instead of fixed value and
!              moved from data_obs_cdfin to data_obs_lib_cosmo.
!  24.04.2023: Adapted interface of 'f_z2p'.
!  20.02.2024: Set entry 'nbsvip' to model surface pressure for later use as
!              'plevel' (for localisation in LETKF)
!
! ==============================================================================
! DACE history:
! ------------
!  1.51       2017-02-24 Andreas Rhodin : update shared modules from COSMO 5.04d
!             2018-04-13 Christoph Schraff
!  update shared modules from COSMO 5.05, except for option based on roughness
!  length and/or Laplacian of orography to select 10-m wind stations!
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
    zlimv10        ,& ! additional limits for the use of 10-m wind observations:
                      !   pos./neg. Laplacian of orography, surface roughness
    rhtsat         ,& ! relative humidity threshold above which obs is set =100%

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
    r_z0              ! surface roughness length

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP report
    nairep     ,& ! AIREP report (all aircraft reports)
    nsatob     ,& ! SATOB report
    ndribu     ,& ! DRIBU report
    nscatt     ,& ! SCATT report
!   namdar     ,& !   amdar report
!   nacar      ,& !   acar  report

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
!   ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
!   ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
!   ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
!   ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
!   ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
    ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
    ncdf_amdar     ,& ! indicator for proc. NetCDF AMDAR single-level input
    ncdf_acars     ,& ! indicator for proc. NetCDF ACARS single-level input
    ncdf_modes     ,& ! indicator for proc. NetCDF MODE-S KNMI format input
    ncdf_modes_acr ,& ! indicator for proc. NetCDF MODE-S ACARS fmt.  input
!   ncdf_wprof     ,& ! indicator for proc. NetCDF wind profiler      input
!   ncdf_rass      ,& ! indicator for proc. NetCDF RASS profiler      input
!   ncdf_radar_vad ,& ! indicator for proc. NetCDF radar wind prof.   input
    ncdf_synop     ,& ! indicator for proc. NetCDF SYNOP              input
    ncdf_synop_mob ,& ! indicator for proc. NetCDF SYNOP mobile       input
    ncdf_ship      ,& ! indicator for proc. NetCDF SHIP               input
    ncdf_buoy      ,& ! indicator for proc. NetCDF BUOY               input
    ncdf_metar     ,& ! indicator for proc. NetCDF METAR sfc aviation input
!   ncdf_gps_zenith,& ! indicator for proc. NetCDF GPS (ZPD / IWV)    input
    ncdf_ascat     ,& ! indicator for proc. NetCDF ASCAT scatterom.   input
    ncdf_qscat     ,& ! indicator for proc. NetCDF QuickScat scatter. input
    ncdf_satob     ,& ! indicator for proc. NetCDF SATOB wind         input
    ncdf_acars_uk  ,& ! indicator for proc. NetCDF ACARS UK + Canada  input
    ncdf_acars_us  ,& ! indicator for proc. NetCDF ACARS US w. humid. input
    ncdf_nat_dwd   ,& ! indicator for proc. NetCDF DWD national surf. input
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
! derived common header variables stored to ODR
    iobstot        ,& ! longitudinal index of grid pt. to which obs is assigned
    jobstot        ,& ! latitudinal  index of grid pt. to which obs is assigned
    iobsloc        ,& ! longitudinal index of grid pt. in local sub-domain
    jobsloc        ,& ! latitudinal  index of grid pt. in local sub-domain
! auxilliary variable, only temporarily available in reader routine
    irproc         ,& ! indices of reports to be processed now

!         1.2.2 other NetCDF header entries
!               ---------------------------

    nc_nix         ,& ! NIX   , BUFR Table B 002001 : station type (man,auto,..)
    nc_tbuoy       ,& ! MTODB , BUFR Table B 002149 : type of data buoy
    nc_phase       ,& ! MPHAI , BUFR Table B 008004 : phase of aircraft flight
    nc_phasd          ! NDEPF , BUFR Table B 008009 : detailed phase of flight

USE data_obs_cdfin, ONLY :  &

!         1.3   NetCDF body entries
!               -------------------

!         1.3.1  frequent entries
!                ----------------
    nc_nlev        ,& ! MEDRE : delayed descriptor replication factor
                      !         (e.g. number of vertical levels)
    nc_dt          ,& ! NLTPD : time [sec] since launch time
    nc_lvtyp       ,& ! MEVSS , BUFR Tab 008042 : extended vertical sounding
                      !                           significance (level identity)
    nc_z           ,& ! NHHHN : geopotential height [gpm]    (upper-air)
    nc_dd          ,& ! NDNDN : wind direction      [degree]
    rc_p           ,& ! MPN   , MPPP : pressure
!   rc_dlat        ,& ! MLADH : latitude  displacement since launch site
!   rc_dlon        ,& ! MLODH : longitude displacement since launch site
    rc_t           ,& ! MTDBT : temperature / dry-bulb temperature
    rc_td          ,& ! MTDNH : dew-point temperature
    rc_ff          ,& ! NFNFN : wind speed  (or NFF for scatterometer)
    rc_dd          ,& ! NDD   : wind direction (for scatterometer only)
    rc_lat0        ,& ! MLAH0 : latitude
    rc_lon0        ,& ! ML0H0 : longitude

!         1.3.2  additional aircraft elements
!                ----------------------------
    nc_rolla       ,& ! MQARA , BUFR Tab 002064 : aircraft roll angle quality
    nc_turb        ,& ! MB    , BUFR Tab 011031 : degree of turbulence
    nc_qxq         ,& ! MMRQ  , BUFR Tab 033026 : mixing ratio quality
    nc_rh          ,& ! MUUU  : relative humidity                           [%]
    rc_rh          ,& ! MUUU  : relative humidity (for MODE-S)              [%]
    rc_qx          ,& ! MMIXR : mixing ratio                            [kg/kg]
    rc_tq          ,& ! MPOTO : precision of temperature observation        [K]
    rc_vgust       ,& ! NMDEWX: maximum derived equivalent vertical gust  [m/s]

!         1.3.2  additional aircraft elements
!                ----------------------------
    rc_z              ! MHHH  : height or altitude (vertical location)

USE data_obs_cdfin, ONLY :  &

!         1.3.4  additional synoptic elements
!                ----------------------------
    nc_uuu         ,& ! MUUU  : relative humidity
    nc_fi          ,& ! NHHHN : geopotential height [gpm] of standard level
    nc_vtisi       ,& ! MTISI*: BUFR Tab 008021 : time significance (wind obs)
    nc_vdt         ,& ! NGGTP : time period of wind measurement           [min]
    nc_gudt        ,& ! NGGTP0: time period of (1st) max. wind gust speed [min]
    nc_gudt2       ,& ! NGGTP0: time period of  2nd  max. wind gust speed [min]
    nc_tmxt1       ,& ! MGGTP*: time displacement: start of period of T-max [h]
    nc_tmxt2       ,& ! MGGTP*: time displacement: end   of period of T-max [h]
    nc_tmnt1       ,& ! MGGTP*: time displacement: start of period of T-min [h]
    nc_tmnt2       ,& ! MGGTP*: time displacement: end   of period of T-min [h]
    nc_rrdt        ,& ! MGGTP*: time period of (first) precipitation obs    [h]
    nc_rrdt2       ,& ! MGGTP*: time period of  second precipitation obs    [h]
                      !         (for DWD national reports, the unit is minutes)
    nc_clct        ,& ! MN                      : total cloud amount [%]
    nc_clsig       ,& ! MVTSU , BUFR Tab 008002 : vertical significance
    nc_clxsg       ,& ! MVTSU*, BUFR Tab 008002 :    dito    (additional levels)
    nc_clclm       ,& ! MNH   , BUFR Tab 020011 :(low or mid-level) cloud amount
    nc_clxcl       ,& ! MNH*  , BUFR Tab 020011 :    dito    (additional levels)
    nc_ccl         ,& ! MCC   , BUFR Tab 020012 : cloud type (low clouds)
    nc_clxct       ,& ! MCC*  , BUFR Tab 020012 :    dito    (additional levels)
    nc_ccm         ,& ! MCC0  , BUFR Tab 020012 : cloud type (middle clouds)
    nc_cch         ,& ! MCC1  , BUFR Tab 020012 : cloud type (high clouds)
    nc_clev        ,& ! MDREP*, replication factor (for additional cloud levels)
    nc_wdt         ,& ! MGGTP : time period [h] for past weather
    nc_ww          ,& ! WW    , BUFR Tab 020003 : present weather
    nc_w1          ,& ! W1    , BUFR Tab 020004 : past    weather (1)
    nc_w2          ,& ! W2    , BUFR Tab 020005 : past    weather (2)
    nc_e              ! ME    , BUFR Tab 020062 : state of ground (w/wo snow)

USE data_obs_cdfin, ONLY :  &

    rc_pfi         ,& ! MPN   : pressure of standard level  (surface report)
    rc_pmsl        ,& ! MPPPP : pressure reduced to mean sea level
    rc_dpdt        ,& ! NPPP  : 3-hour pressure change
    rc_gust        ,& ! NFXGU : maximum wind speed of gusts
    rc_gust2       ,& ! NFXGU : maximum wind speed of gusts, second period
    rc_tmax        ,& ! MTXTXH: maximum temperature \ height, period specified,
    rc_tmin        ,& ! MTNTNH: minimum temperature / processed if: 2-m , 12-h
    rc_zt          ,& ! MHOSEN*:height of T-sensor above local ground       [m]
    rc_vis         ,& ! MVV(VV):horizontal visibility
    rc_rr24        ,& ! MRR24 : total precipitation over past 24 hours
    rc_rr          ,& ! MRRR  : total precipitation (first period)
    rc_rr2         ,& ! MRRR  : total precipitation, second period
    rc_cbase       ,& ! NH    : height of base of cloud
    rc_clxbs       ,& ! NH*   : height of base of cloud (additional levels)
    rc_radgl       ,& ! MGLSR : global  solar radiation                  [J/m2]
    rc_raddf       ,& ! MDSRH : diffuse solar radiation                  [J/m2]
    rc_radlw       ,& ! MLWR  : long-wave     radiation                  [J/m2]
    rc_hsnow       ,& ! NSSS  : total snow depth                            [m]

!         1.3.5  additional SATOB/AMV elements
!                -----------------------------
    nc_rtyp        ,& ! MSDWCM, BUFR Tab 002023: sat.-derived wind comp. method
    nc_qual        ,& ! MPCCO : per cent confidence                         [%]
    rc_sccf        ,& ! MSCCF : satellite channel centre frequency    [dp: 1/s]

!         1.5   auxilliary buffer arrays and variables
!               --------------------------------------

    imiss          ,& ! missing value for integers in current NetCDF input file
    rmiss          ,& ! missing value for reals    in current NetCDF input file
    rmisschk          ! value smaller than 'rmiss', to check for missing value

USE data_obs_cdfin, ONLY :  &

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
    nedmis     ,& ! wind direction missing
    nefmis     ,& ! wind speed missing
    nedflg     ,& ! wind direction flagged
    nefflg     ,& ! wind speed flagged
    nefneg     ,& ! wind speed too small  ( < 0 ; DRIBU: <= 0)
    nevalt     ,& ! height (diff.) too large for 10m-wind
!   nefshr     ,& ! wind speed shear too large
!   nedshr     ,& ! directional wind shear too large
    nerlim        ! prec.amount exceeds threshold limit

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
    oeairep    ,& ! (root of) AIREP wind   error variance
    oetairp    ,& ! (root of) AIREP temperature error variance
    oevsynp    ,& ! (root of) SYNOP wind   error variance
    oezsynp    ,& ! (root of) SYNOP height error variance (land)
    oesatob    ,& ! (root of) SATOB wind   error variance
    oevscat    ,& ! (root of) SCATT wind   error variance
    oevdrib    ,& ! (root of) DRIBU wind   error variance
    oezdrib    ,& ! (root of) DRIBU height error variance
    oezship    ,& ! (root of) SHIP (sea SYNOP) height error variance
    rherr1     ,& ! (root of) fixed    / normal conditions
    rherr2     ,& ! relative humidity <  if temperature below 233K
    rherr3        ! error variances    \ if rel. humidity below 20%

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
    vfoglim    ,& ! visibility threshold [m] below which the existence of low
                  !      cloud (fog) is assumed in the presence of precipitation
!   ccl_lim    ,& ! height limit below which cloud is considered low  cloud [m]
!   cch_lim    ,& ! height limit above which cloud is considered high cloud [m]
    h_low      ,& ! height limit below which cloud is low cloud       [m]
    h_mid      ,& ! height limit below which cloud is mid-level cloud [m]
    cbh_clr    ,& ! cloud base height value being assigned if clear sky [m]

!         7.0    For data rejection messages: Output buffer, size and formats
!                ------------------------------------------------------------
    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt1      ,& ! no pressure
    nfmt2      ,& ! excess of precipitation
    nfmt3      ,& ! no accepted data
    nfmt6      ,& ! excess of pressure tendency
    nfmt20        ! message only: fog and precipitation

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

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
!   mxrhdf     ,& ! header length of multi-level reports
!   mxshdf     ,& ! header length of single-level reports
    nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw     ,& ! status of QC and of writing to feedobs files (and p-QC flag)
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
!   nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
!   nhcent     ,& ! originating centre
    nhstid     ,& ! station identity number
!   nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
!   nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
!   nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
!   nhuexi     ,& ! flag for existence of wind data        in multi-level report
!   nhtexi     ,& ! flag for existence of temperature data in multi-level report
!   nhqexi     ,& ! flag for existence of humidity data    in multi-level report
!   nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
!   nhtrac     ,& ! tracking technique (NSASA, see WMO common code Table C7)
!   nhrad      ,& ! solar and IR radiation correction (NSR, BUFR Table 002013)
!   nhna4      ,& ! instrument type                   (NA4, BUFR Table 002003)
!   nhwce      ,& ! wind comput. enhancement (w-prof, MWCE, BUFR Table 025021)
!   nhdt       ,& ! time period of measurement (e.g. w-prof)               [s]
    nhstyp     ,& ! surface obs: station type (buoy:  MTODB, BUFR Table 002149,
                  !                            satob: MSDWCM,BUFR Table 002023)
                  !                            else:  NIX  , BUFR Table 002001)


!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid     ,& ! character length of the station identity
    ilstidp    ,& ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
!   nvpabp     ,& ! bit pos. for report set passive since it is     nhschr
                  !              used in a multi-level pseudo report
    nvpsbp     ,& ! bit pos. for report set passive since 1 of next   "
                  !              5 flags or flight track flag applies
!   nvobbp     ,& ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
    nvalbp     ,& ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
!   nvbkbp     ,& ! bit pos. for flag: 'blacklisted station (ship)'   "
!   nvexbp     ,& ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
!   nvrdbp     ,& ! bit pos. for flag: 'redundant report'             "
    nvsebp     ,& ! bit pos. for report located at sea grid pt.       "
    nvlzbp     ,& ! bit pos. for vertical levels reported as height   "
    nvapbp     ,& ! bit pos. for phase of flight (aircraft)           "
    nvapoc     ,& ! no. of bits occ. by phase of flight               "
    nvaabp     ,& ! bit pos. for aircraft roll angle (code)           "
    nvaaoc        ! no. of bits occ. by aircraft roll angle           "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

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
    nbspst     ,& ! (3-hourly) pressure tendency                       [Pa/3h]
    nbsviz     ,& ! scaled extrapolat. distance for surf. pressure obs [m]
    nbsvip     ,& ! model surface pressure (for LETKF localisation)    [Pa]
    nbscbs     ,& ! (lowest) cloud base height CBH above surface       [m]
    nbscil     ,& ! ceiling (CBH of lowest cloud layer >= 5 octas)     [m]
    nbscl      ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm      ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch      ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct      ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis     ,& ! (horizontal) visibility                            [m]
    nbsff      ,& ! wind speed                                         [m/s]
    nbsdd      ,& ! wind direction                                     [deg]
    nbstd         ! dewpoint temperature                               [K]

USE data_obs_record, ONLY :   &

    nbsrr1     ,& ! precipitation amount over 1 hour                   [mm]
    nbsrr3     ,& ! precipitation amount over 3 hours                  [mm]
    nbsrr6     ,& ! precipitation amount over 6 hours                  [mm]
    nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
    nbsfgv     ,& ! max. derived equivalent vertical gust (aircraft)   [m/s]
    nbsfg1     ,& ! max. wind speed of gusts over 1 hour               [m/s]
    nbsfg3     ,& ! max. wind speed of gusts over 3 hour               [m/s]
    nbsfg6     ,& ! max. wind speed of gusts over 6 hours              [m/s]
    nbstn      ,& ! minimum temperature (at 2m, in period 'nbsttr')    [K]
    nbstx      ,& ! maximum temperature (at 2m, in period 'nbsttr')    [K]
    nbsrad     ,& ! global    solar    radiation, sum over 1 hour      [J/m2]
    nbsrdd     ,& ! diffuse   solar    radiation, sum over 1 hour      [J/m2]
    nbsrdl     ,& ! long-wave downward radiation, sum over 1 hour      [J/m2]
    nbshsw     ,& ! total snow depth                                   [m]
    nbsplv     ,& ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbsdrh        ! bias correction for relative humidity [/]

USE data_obs_record, ONLY :   &

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
!   mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbserr     ,& ! pre-processing status flags  (bit p., see below: 'nb?err')
    nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg     ,& ! combined cloud and weather group (set of classes, s below)
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

USE data_obs_record, ONLY :   &

!       1.3.9  Bit patterns for packed information in body 'mosgbd'
!              ----------------------------------------------------
                  ! --> weather and ground group word                  nbswwe
    nvw0bp     ,& ! bit position for ww (present wea.) code (WMO Table 020003)
    nvw1bp     ,& ! bit position for w  (past weather) code (WMO Table 020004)
    nvwtbp     ,& ! bit position for time period of w                  [h]
    nvcqbp     ,& ! bit position for refined quality flag on ccl
    nveebp     ,& ! bit position for state of ground        (WMO Table 020062)
    nrtrbp     ,& ! bit position for code of precipitation
                  !     measurement duration       [Code table 4019, keys 0-7]
    nvw0oc     ,& ! no. bits occupied for ww code           (WMO Table 020003)
    nvw1oc     ,& ! no. bits occupied for w  code           (WMO Table 020004)
    nvwtoc     ,& ! no. bits occupied for time period of w             [h]
    nvcqoc     ,& ! no. bits occupied by refined quality flag on ccl
    nveeoc     ,& ! no. bits occupied for state of ground   (WMO Table 020062)
    nrtroc     ,& ! no. bits occupied by precip obs. duration code
                  !
                  ! --> general    cloud       group word              nbsclg
!   nxclbp     ,& ! bit position for low/middle cloud amount(WMO Table 020011)
    nctlbp     ,& ! bit position for low    cloud type code (WMO Table 020012)
    nctmbp     ,& ! bit position for middle cloud type code (WMO Table 020012)
    ncthbp     ,& ! bit position for high   cloud type code (WMO Table 020012)
    nxsgbp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
!   nxsgoc     ,& ! no. bits occupied for vert. signf. code (WMO Table 008002)
!   nxcloc     ,& ! no. bits occupied for low/mid cld amount(WMO Table 020011)
!   nxctoc     ,& ! no. bits occupied for   cloud type code (WMO Table 020012)
                  !
                  ! --> individual cloud layer group words             nbscl?
    nxclbp     ,& ! bit position for cloud amount      code (WMO Table 020011)
    nxctbp     ,& ! bit position for cloud type        code (WMO Table 020012)
    nxbsbp     ,& ! bit position for cloud base height                 [m]
    nxsibp     ,& ! bit position for vertic. signific. code (WMO Table 008002)
    nxcloc     ,& ! no. bits occupied for cloud amount code (WMO Table 020011)
    nxctoc     ,& ! no. bits occupied for cloud type   code (WMO Table 020012)
    nxbsoc     ,& ! no. bits occupied for cloud base height            [m]
    nxsgoc        ! no. bits occupied for vert. signf. code (WMO Table 008002)

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status/QC flags for horiz. wind   nb?err/nb?qcf
    nvrt       ,& ! bit pos. for status/QC flags for temperature         "
    nvrq       ,& ! bit pos. for status/QC flags for humidity            "
    nvrz       ,& ! bit pos. for status/QC flags for pressure/height     "
    nvrw       ,& ! bit pos. for status/QC flags for vertical wind       "
    nvfubp     ,& ! bit pos. for main flag on wind                     nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature                "
    nvfqbp     ,& ! bit pos. for main flag on humidity                   "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.         "
    nvfaoc     ,& ! no. of bits occ. by each main flag                   "
    nvfbps     ,& ! bit pattern for main flags:                          "
    nvfboc     ,& ! no. of bits occ. for main flags                      "
!   nvflbp     ,& ! bit pos. for level flag: level below surface         "
    ntxbp      ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp      ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc      ,& ! no. of bits occ. by each time period              "
    nvlidp     ,& ! level id. bit pattern                              nb?lid
!   nvlido     ,& ! no. bits occ. by each indicator in level id.         "

!       1.4.3  Bit patterns for 'optional groups' in ODR body 'mosgbd' (and VOF)
!              -----------------------------------------------------------------
                  ! combined cloud and weather group                   nvbcwg
    nvchbp     ,& ! bit position for ch  (type of high cloud)
    nvcmbp     ,& !         "        cm  (type of middle cloud)
    nvclbp     ,& !         "        cl  (type of low cloud)
    nvnhbp     ,& !         "        nh  (cover of low, else of middle cloud)
    nvhbp      ,& !         "        h   (cloud base height)
    nvnbp      ,& !         "        n   (total cloud cover)
    nvwwbp     ,& !         "        ww  (present weather)
                  !                      (see VUB WMO Code tables:)
    nvchoc     ,& ! no. of bits occupied by ch    [Code table 0509]
    nvcmoc     ,& !           "             cm    [Code table 0515]
    nvcloc     ,& !           "             cl    [Code table 0513]
    nvnhoc     ,& !           "             nh    [Code table 2700]
    nvhoc      ,& !           "             h     [Code table 1600]
    nvnoc      ,& !           "             n     [Code table 2700]
    nvwwoc        !           "             ww    [Code table 4677]

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotsg     ,& ! tot. number of stored single-level reports
    fdoro      ,& ! scaling factors to obs height minus model orography

!       2.     Observation data records (ODR) and associated arrays
!       ----------------------------------------------------------------------

!       2.1    Formats of observation data records
!              -----------------------------------
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

    FL_HEIGHT     ,& ! location not in valid height range
    FL_NO_OBS        ! no (active) observations in report

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
    uv2uvrot_vec       ! converts wind components from normal to rotated system

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_assign_sort_node   ,& ! assign node to reports and sort them accordingly
    obs_cdf_distrib_reports,& ! distribute reports to appropriate sub-domains
    std_atmosphere         ,& ! convert variables accord. to standard atmosphere
    obs_td2rh              ,& ! convert dewpoint temperature to relat. humidity
    obs_qx2rh              ,& ! convert mixing ratio to model-compatible rel hum
    obs_rhw2rh             ,& ! make (observed) relat. humidity model compatible
    obs_rh2td              ,& ! convert relat. humidity to dewpoint temperature
    obs_find_level         ,& ! get interpolation levels and factors
    check_wind_method      ,& ! add sat frequency digit to wind comput. method
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

 USE src_obs_cdfin_mult,          ONLY :  &
    obs_cdf_store_multilev    ! storing multi-level reports in ODR

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
PUBLIC :: obs_cdf_read_surface,  obs_cdf_read_aircraft,  obs_cdf_read_amv

  LOGICAL  ,  PARAMETER    ::  &
    lax  =  .TRUE.    ! (re-)lax(ed) conditions on required variables in
                      !   NetCDF feedobs files

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_sing" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_surface ( min_sta , min_end , ilcdf                    &
                                , nsgnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_sing" organizes the reading,
!   pre-selection, distribution and storage in ODR of surface (SYNOP, SHIP,
!   BUOY, METAR) reports from a given observation time period.
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
    nsgnew        ,& ! number of new single-level reports
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
                      !  = ncdf_synop     : SYNOP
                      !  = ncdf_synop_mob : SYNOP mobile
                      !  = ncdf_ship      : SHIP
                      !  = ncdf_buoy      : BUOY
                      !  = ncdf_metar     : METAR (surface aviation)
                      !  = ncdf_synop_mob : test SYNOP
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
    ilev           ,& ! loop index over vertical levels
    nlev           ,& ! number of vertical levels
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    status         ,& ! NetCDF status variable
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
!   ktimsig        ,& ! significance of report time
    istat  , ierr     ! error indicators

  INTEGER        ::  &          ! variable ID's in NetCDF file for:
!                              !   all file types (synop, ship, buoy, metar):
    varid_nix               ,& ! station type  (0:auto,1:manned,2:hybrid,3:miss)
    varid_t    , varid_td   ,& ! dry bulb temperature , dew point temperature
    varid_dd   , varid_ff   ,& ! wind direction       , wind speed
    varid_niswv             ,& ! index of selected wind vector (scatterometer)
    varid_gust              ,& ! max. wind gust speed
    varid_zt                ,& ! T-sensor height above ground
!                              !   synop, ship, buoy (but not metar):
    varid_rh   , varid_gudt ,& ! relative humidity   , time period of max. gusts
    varid_vtisi, varid_vdt  ,& ! time significance   , time period of wind obs
    varid_rr   , varid_rrdt ,& ! total precipitation , period of precip obs
    varid_rr2  , varid_rrdt2,& ! total precipitation , period of precip obs
    varid_p    , varid_pmsl ,& ! station pressure    , mean sea level pressure
    varid_dpdt              ,& ! pressure tendency
!                              !   synop, ship, metar (but not buoy):
    varid_vis  , varid_clev ,& ! visibility,    , cloud layer replication factor
    varid_clxsg, varid_clxcl,& ! vertical significance, cloud amount
    varid_clxct, varid_clxbs,& ! cloud type           , cloud base height
!                              !   synop and ship:
    varid_rr24 , varid_clct ,& ! 24-hr precipitation  , total cloud cover [%]
    varid_clsig, varid_clclm,& ! vertical significance, low/middle cloud amount
    varid_cbase, varid_ctcl ,& ! cloud base height [m], low  cloud type [020012]
    varid_ctcm , varid_ctch ,& ! middle cloud type    , high cloud type
    varid_ww   , varid_wdt  ,& ! present weather  , period [hrs] of past weather
    varid_w1   , varid_w2   ,& ! past    weather (1+2)
    varid_tmxt1, varid_tmxt2,& ! time interval [hrs] for maximum temperature
    varid_tmnt1, varid_tmnt2,& ! time interval [hrs] for minimum temperature
    varid_tmax , varid_tmin ,& ! maximum temperature  , minimum temperature
!                              !   synop only:
    varid_radgl, varid_raddf,& ! global solar radiat. , diffuse solar radiation
    varid_radlw, varid_raddt,& ! long-wave radiation  , time period of radiation
    varid_e    , varid_hsnow,& ! state of ground [020062] , total snow depth [m]
    varid_pfi  , varid_fi   ,& ! pressure of , geopot. height at standard level
!                              !   buoy  only:
    varid_tbuoy                ! type of data buoy [002149]
!
!                              !   for fixed land synop, ship, buoy:
!   varid_sst  , varid_sstz ,& ! sea / water temperature, depth of SST obs
!   varid_zv                   ! wind sensor above ground or moving platform
!                              !   for ship and buoy:
!   varid_stadd, varid_staff,& ! direction, speed of moving platform
!   varid_ht   , varid_hv      ! T-, wind sensor heights above sea level
!                              !   for buoy:
!   varid_sstrp                ! replication factor for SST obs
!                              !   for fixed land synop only (not mobile synop):
!   varid_fmedt, varid_fme10,& ! time period [min], max. 10-min mean wind speed
!   varid_tgmin                ! minimum (12-hr) ground temperature

  LOGICAL        ::  &
    linsitu        ,& ! file with in-situ observations
    lscatto           ! file with scatterometer data

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
    ibufloc    (:) ,& ! integer buffer array with local reports (at sub-domain)
    nc_niswv   (:) ,& ! NISWV: index of selected wind vector (for scatterometer)
    nc_raddt (:,:)    ! MGGTP7: time period of radiation                    [h]

  REAL (KIND=wp) , ALLOCATABLE :: &
    rbufsrt    (:) ,& ! real    buffer array with sorted reports read from file
    rbufloc    (:) ,& ! real    buffer array with local reports (at sub-domain)
    rc_ff4   (:,:) ,& ! NFF  : wind speed (4 possible values, for scatterometer)
    rc_dd4   (:,:) ,& ! NDD  : wind direction (4 values, for scatterometer)
    rc_radg2   (:) ,& ! MGLSR : global  solar radiation, second period   [J/m2]
    rc_radd2   (:) ,& ! MDSRH : diffuse solar radiation, second period   [J/m2]
    rc_radl2   (:)    ! MLWR  : long-wave     radiation, second period   [J/m2]

  CHARACTER (LEN=ilstidn) , ALLOCATABLE :: &
    ybufsrt    (:) ,& ! character array with sorted reports read from file
    ybufloc    (:)    ! character array with local reports (at sub-domain)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_surface
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_surface'
  status   = nf90_noerr

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

! no differences between synop and synop mobile, except for station ID and
! station height quality - which is considered in 'obs_cdf_read_comhead'
  IF (kcdftyp == ncdf_synop_mob)  kcdftyp = ncdf_synop
  IF (kcdftyp == ncdf_synop_tst)  kcdftyp = ncdf_synop

  lscatto  =      (kcdftyp == ncdf_ascat) .OR. (kcdftyp == ncdf_qscat)
  linsitu  =      (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )         &
             .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_metar)

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
! (i.e. elements without replication or with fixed replication)
! to be put into 'Xbufsrt'  (specific for the different file types)
! -----------------------------------------------------------------
  IF     (kcdftyp == ncdf_synop) THEN
    noffscal (1,2)  =  noffscal(1,1) + 24
    noffscal (2,2)  =  noffscal(2,1) + 21
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_ship ) THEN
    noffscal (1,2)  =  noffscal(1,1) + 22
    noffscal (2,2)  =  noffscal(2,1) + 16
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_buoy ) THEN
    noffscal (1,2)  =  noffscal(1,1) + 8
    noffscal (2,2)  =  noffscal(2,1) + 9
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_metar) THEN
    noffscal (1,2)  =  noffscal(1,1) + 3
    noffscal (2,2)  =  noffscal(2,1) + 6
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_nat_dwd) THEN
    noffscal (1,2)  =  noffscal(1,1) + 8
    noffscal (2,2)  =  noffscal(2,1) + 5
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (lscatto) THEN
    noffscal (1,2)  =  noffscal(1,1) + 0
    noffscal (2,2)  =  noffscal(2,1) + 2
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

! 'irnode' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get entries specific to observation type
!-------------------------------------------------------------------------------

! allocate arrays for elements without replication or with fixed replication
! --------------------------------------------------------------------------

    IF (lscatto) THEN
      ALLOCATE ( rc_dd  (1,nrep) , STAT=istat )
      ALLOCATE ( rc_ff  (1,nrep) , STAT=istat )
      ALLOCATE ( rc_ff4 (4,nrep) , STAT=istat )
      ALLOCATE ( rc_dd4 (4,nrep) , STAT=istat )
      ALLOCATE ( nc_niswv (nrep) , STAT=istat )
    ENDIF

    IF (linsitu) THEN
      ALLOCATE ( nc_nix   (nrep) , STAT=istat )
      ALLOCATE ( rc_t   (1,nrep) , STAT=istat )
      ALLOCATE ( rc_td  (1,nrep) , STAT=istat )
      ALLOCATE ( nc_dd  (1,nrep) , STAT=istat )
      ALLOCATE ( rc_ff  (1,nrep) , STAT=istat )
      ALLOCATE ( rc_gust  (nrep) , STAT=istat )
      ALLOCATE ( rc_zt    (nrep) , STAT=istat )

      ALLOCATE ( nc_clev  (nrep) , STAT=istat )

      DO irep = 1 , nrep
        rc_gust  (irep)  =  rmiss
        rc_zt    (irep)  =  rmiss
        nc_clev  (irep)  =  imiss
      ENDDO
    ENDIF

    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_buoy )) THEN
      ALLOCATE ( rc_p   (1,nrep) , STAT=istat )
      ALLOCATE ( rc_pmsl  (nrep) , STAT=istat )
      ALLOCATE ( rc_dpdt  (nrep) , STAT=istat )
      ALLOCATE ( rc_rr    (nrep) , STAT=istat )
      ALLOCATE ( nc_uuu   (nrep) , STAT=istat )
      ALLOCATE ( nc_rrdt  (nrep) , STAT=istat )
      ALLOCATE ( nc_vdt   (nrep) , STAT=istat )
      ALLOCATE ( nc_gudt  (nrep) , STAT=istat )
      ALLOCATE ( nc_vtisi (nrep) , STAT=istat )
      DO irep = 1 , nrep
!       rc_p   (1,irep)  =  rmiss
        rc_pmsl  (irep)  =  rmiss
        rc_dpdt  (irep)  =  rmiss
        rc_rr    (irep)  =  rmiss
        nc_uuu   (irep)  =  imiss
        nc_rrdt  (irep)  =  imiss
        nc_vdt   (irep)  =  imiss
        nc_gudt  (irep)  =  imiss
        nc_vtisi (irep)  =  imiss
      ENDDO
    ENDIF

    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_metar)) THEN
      ALLOCATE ( rc_vis   (nrep) , STAT=istat )
      DO irep = 1 , nrep
        rc_vis   (irep)  =  rmiss
      ENDDO
    ENDIF

    IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
      ALLOCATE ( rc_gust2 (nrep) , STAT=istat )
      ALLOCATE ( nc_gudt2 (nrep) , STAT=istat )
      ALLOCATE ( rc_rr24  (nrep) , STAT=istat )
      ALLOCATE ( rc_rr2   (nrep) , STAT=istat )
      ALLOCATE ( nc_rrdt2 (nrep) , STAT=istat )
      ALLOCATE ( nc_clct  (nrep) , STAT=istat )
      ALLOCATE ( nc_clsig (nrep) , STAT=istat )
      ALLOCATE ( nc_clclm (nrep) , STAT=istat )
      ALLOCATE ( rc_cbase (nrep) , STAT=istat )
      ALLOCATE ( nc_ccl   (nrep) , STAT=istat )
      ALLOCATE ( nc_ccm   (nrep) , STAT=istat )
      ALLOCATE ( nc_cch   (nrep) , STAT=istat )
      ALLOCATE ( nc_ww    (nrep) , STAT=istat )
      ALLOCATE ( nc_wdt   (nrep) , STAT=istat )
      ALLOCATE ( nc_w1    (nrep) , STAT=istat )
      ALLOCATE ( nc_w2    (nrep) , STAT=istat )
      ALLOCATE ( rc_tmax  (nrep) , STAT=istat )
      ALLOCATE ( rc_tmin  (nrep) , STAT=istat )
      ALLOCATE ( nc_tmxt1 (nrep) , STAT=istat )
      ALLOCATE ( nc_tmxt2 (nrep) , STAT=istat )
      ALLOCATE ( nc_tmnt1 (nrep) , STAT=istat )
      ALLOCATE ( nc_tmnt2 (nrep) , STAT=istat )
      DO irep = 1 , nrep
        rc_gust2 (irep)  =  rmiss
        nc_gudt2 (irep)  =  imiss
        rc_rr24  (irep)  =  rmiss
        rc_rr2   (irep)  =  rmiss
        nc_rrdt2 (irep)  =  imiss
        nc_clct  (irep)  =  imiss
        nc_clsig (irep)  =  imiss
        nc_clclm (irep)  =  imiss
        rc_cbase (irep)  =  rmiss
        nc_ccl   (irep)  =  imiss
        nc_ccm   (irep)  =  imiss
        nc_cch   (irep)  =  imiss
        nc_ww    (irep)  =  imiss
        nc_wdt   (irep)  =  imiss
        nc_w1    (irep)  =  imiss
        nc_w2    (irep)  =  imiss
        rc_tmax  (irep)  =  rmiss
        rc_tmin  (irep)  =  rmiss
        nc_tmxt1 (irep)  =  imiss
        nc_tmxt2 (irep)  =  imiss
        nc_tmnt1 (irep)  =  imiss
        nc_tmnt2 (irep)  =  imiss
      ENDDO
    ENDIF

    IF     (kcdftyp == ncdf_synop) THEN
      ALLOCATE ( nc_e      (nrep) , STAT=istat )
      ALLOCATE ( nc_raddt(2,nrep) , STAT=istat )
      ALLOCATE ( rc_radgl  (nrep) , STAT=istat )
      ALLOCATE ( rc_raddf  (nrep) , STAT=istat )
      ALLOCATE ( rc_radlw  (nrep) , STAT=istat )
      ALLOCATE ( rc_hsnow  (nrep) , STAT=istat )
      ALLOCATE ( rc_pfi    (nrep) , STAT=istat )
      ALLOCATE ( nc_fi     (nrep) , STAT=istat )
      ALLOCATE ( rc_radg2  (nrep) , STAT=istat )
      ALLOCATE ( rc_radd2  (nrep) , STAT=istat )
      ALLOCATE ( rc_radl2  (nrep) , STAT=istat )
      DO irep = 1 , nrep
        nc_e     (irep)  =  imiss
        rc_radgl (irep)  =  rmiss
        rc_raddf (irep)  =  rmiss
        rc_radlw (irep)  =  rmiss
        rc_hsnow (irep)  =  rmiss
        rc_pfi   (irep)  =  rmiss
        nc_fi    (irep)  =  imiss
        rc_radg2 (irep)  =  rmiss
        rc_radd2 (irep)  =  rmiss
        rc_radl2 (irep)  =  rmiss
      ENDDO
      nc_raddt = imiss
    ELSEIF (kcdftyp == ncdf_buoy ) THEN
      ALLOCATE ( nc_tbuoy (nrep) , STAT=istat )
      DO irep = 1 , nrep
        nc_tbuoy (irep)  =  imiss
      ENDDO
    ELSEIF (kcdftyp == ncdf_nat_dwd ) THEN
      ALLOCATE ( rc_rr    (nrep) , STAT=istat )
      ALLOCATE ( rc_rr2   (nrep) , STAT=istat )
      ALLOCATE ( rc_gust  (nrep) , STAT=istat )
      ALLOCATE ( rc_gust2 (nrep) , STAT=istat )
      ALLOCATE ( rc_ff  (1,nrep) , STAT=istat )
      ALLOCATE ( nc_nix   (nrep) , STAT=istat )
      ALLOCATE ( nc_rrdt  (nrep) , STAT=istat )
      ALLOCATE ( nc_rrdt2 (nrep) , STAT=istat )
      ALLOCATE ( nc_gudt  (nrep) , STAT=istat )
      ALLOCATE ( nc_gudt2 (nrep) , STAT=istat )
      ALLOCATE ( nc_dd  (1,nrep) , STAT=istat )
      ALLOCATE ( nc_vdt   (nrep) , STAT=istat )
      ALLOCATE ( nc_vtisi (nrep) , STAT=istat )
      DO irep = 1 , nrep
        rc_rr    (irep)  =  rmiss
        rc_rr2   (irep)  =  rmiss
        rc_gust  (irep)  =  rmiss
        rc_gust2 (irep)  =  rmiss
        rc_ff  (1,irep)  =  rmiss
        nc_nix   (irep)  =  imiss
        nc_rrdt  (irep)  =  imiss
        nc_rrdt2 (irep)  =  imiss
        nc_gudt  (irep)  =  imiss
        nc_gudt2 (irep)  =  imiss
        nc_dd  (1,irep)  =  imiss
        nc_vdt   (irep)  =  imiss
        nc_vtisi (irep)  =  imiss
      ENDDO
    ENDIF

! get mandatory (!) variable ID's in NetCDF file and get data
! -----------------------------------------------------------

    IF (lscatto) THEN
                                status = nf90_inq_varid (ncid,'NDD  ',varid_dd )
      IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NFF  ',varid_ff )
    ENDIF
    IF ((lscatto) .AND. (status == nf90_noerr)) THEN
                              status = nf90_inq_varid (ncid,'NISWV',varid_niswv)
      IF (status == nf90_noerr) THEN
        !   (METOP) data read from file according to Eumetsat template
        status = nf90_get_var (ncid, varid_ff  , rc_ff4  , start=(/1,mrepsta/) &
                                                         , count=(/4,nrep/))
        status = nf90_get_var (ncid, varid_dd  , rc_dd4  , start=(/1,mrepsta/) &
                                                         , count=(/4,nrep/))
        status = nf90_get_var (ncid, varid_niswv,nc_niswv,(/mrepsta/), (/nrep/))
      ! PRINT *,'rmiss', rmiss, rc_dd4(:,1),rc_dd4(:,2), rc_ff4(:,1),rc_ff4(:,2)
        DO irep = 1, nrep
          IF ((nc_niswv(irep) >= 1) .AND. (nc_niswv(irep) <= 4)) THEN
            rc_dd (1,irep) = rc_dd4(nc_niswv(irep),irep)
            rc_ff (1,irep) = rc_ff4(nc_niswv(irep),irep)
          ELSE
            rc_dd (1,irep) = rmiss
            rc_ff (1,irep) = rmiss
          ENDIF
        ENDDO
      ELSE
        !   data read from file written by H.-W. Bitzer's preprocessing tool
        status = nf90_get_var (ncid, varid_dd  , rc_dd  , (/mrepsta/), (/nrep/))
        status = nf90_get_var (ncid, varid_ff  , rc_ff  , (/mrepsta/), (/nrep/))
      ENDIF
    ENDIF

    IF (linsitu) THEN
                                status = nf90_inq_varid (ncid,'MTDBT',varid_t  )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                                status = nf90_inq_varid (ncid,'MTTT ',varid_t  )
        !   for METAR
!       IF ((status /= nf90_noerr) .AND. (lax))                                &
!                               status = nf90_inq_varid (ncid,'MTN  ',varid_t  )
      IF (status == nf90_noerr) THEN
                                status = nf90_inq_varid (ncid,'MTDNH',varid_td )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                                status = nf90_inq_varid (ncid,'MTDTD',varid_td )
        !   for METAR
!       IF ((status /= nf90_noerr) .AND. (lax))                                &
!                               status = nf90_inq_varid (ncid,'MTDN ',varid_td )
      ENDIF
      IF (status == nf90_noerr) THEN
                                status = nf90_inq_varid (ncid,'NDNDN',varid_dd )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                                status = nf90_inq_varid (ncid,'NDD  ',varid_dd )
      ENDIF
      IF (status == nf90_noerr) THEN
                                status = nf90_inq_varid (ncid,'NFNFN',varid_ff )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                                status = nf90_inq_varid (ncid,'NFF  ',varid_ff )
      ENDIF
      IF ((status == nf90_noerr) .AND. (     (kcdftyp == ncdf_synop)           &
                                        .OR. (kcdftyp == ncdf_ship )           &
                                        .OR. (kcdftyp == ncdf_buoy )))         &
                                status = nf90_inq_varid(ncid,'MPPP  ',varid_p  )
    ENDIF
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'MANDATORY SURFACE BODY DATA DO NOT EXIST IN ' // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    IF (linsitu) THEN
      status = nf90_get_var (ncid, varid_t    , rc_t    , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_td   , rc_td   , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_dd   , nc_dd   , (/mrepsta/), (/nrep/))
      status = nf90_get_var (ncid, varid_ff   , rc_ff   , (/mrepsta/), (/nrep/))
    ENDIF
    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_buoy ))                                          &
      status = nf90_get_var (ncid, varid_p    , rc_p    , (/mrepsta/), (/nrep/))

! get optional (!) variable ID's in NetCDF file and
! data for elements without or with fixed replication
! ---------------------------------------------------

    IF (linsitu) THEN
      status = nf90_inq_varid (ncid,'NIX   ',varid_nix )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nix , nc_nix , (/mrepsta/), (/nrep/))
    ENDIF

    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_buoy )) THEN
      status = nf90_inq_varid (ncid,'MPPPP ',varid_pmsl)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_pmsl, rc_pmsl, (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'NPPP  ',varid_dpdt)
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dpdt, rc_dpdt, (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MUUU  ',varid_rh  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_rh  , nc_uuu , (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MHOSEN',varid_zt  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_zt  , rc_zt  , (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP ',varid_vdt )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_vdt , nc_vdt , (/mrepsta/), (/nrep/))
    ENDIF

    IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
      status = nf90_inq_varid (ncid,'MTISI ',varid_vtisi )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_vtisi, nc_vtisi, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MVV   ',varid_vis   )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_vis  , rc_vis  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MRR24 ',varid_rr24  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rr24 , rc_rr24 , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MN    ',varid_clct  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_clct , nc_clct , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MVTSU ',varid_clsig )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_clsig, nc_clsig, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MNH   ',varid_clclm )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_clclm, nc_clclm, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NH    ',varid_cbase )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_cbase, rc_cbase, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC   ',varid_ctcl  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_ctcl , nc_ccl  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC0  ',varid_ctcm  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_ctcm , nc_ccm  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MCC1  ',varid_ctch  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_ctch , nc_cch  , (/mrepsta/),(/nrep/))
      ! individual cloud layers are processed further below
      ! ('MDREP', 'MVTSU0', 'MNH0', 'MCC2', 'NH0': _clev, _clx..)
      status = nf90_inq_varid (ncid,'NWW   ',varid_ww    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_ww   , nc_ww   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MGGTP ',varid_wdt   )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_wdt  , nc_wdt  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MW1   ',varid_w1    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_w1   , nc_w1   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MW2   ',varid_w2    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_w2   , nc_w2   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTXTXH',varid_tmax  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmax , rc_tmax , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTNTNH',varid_tmin  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmin , rc_tmin , (/mrepsta/),(/nrep/))
      IF(kcdftyp == ncdf_synop) status= nf90_inq_varid(ncid,'MGGTP1',varid_rrdt)
      IF(kcdftyp == ncdf_ship ) status= nf90_inq_varid(ncid,'MGGTP0',varid_rrdt)
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_rrdt , nc_rrdt , start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_rrdt , nc_rrdt2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'MRRR  ',varid_rr    )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_rr, rc_rr , (/1,mrepsta/), (/1,nrep/))
        status= nf90_get_var (ncid, varid_rr, rc_rr2, (/2,mrepsta/), (/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'NFXGU ',varid_gust  )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_gust , rc_gust , start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_gust , rc_gust2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'NGGTP0',varid_gudt  )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_gudt , nc_gudt , start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_gudt , nc_gudt2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      IF (kcdftyp==ncdf_synop) status= nf90_inq_varid(ncid,'MGGTP2',varid_tmxt1)
      IF (kcdftyp==ncdf_ship ) status= nf90_inq_varid(ncid,'MGGTP1',varid_tmxt1)
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmxt1, nc_tmxt1, (/mrepsta/),(/nrep/))
      IF (kcdftyp==ncdf_synop) status= nf90_inq_varid(ncid,'MGGTP3',varid_tmxt2)
      IF (kcdftyp==ncdf_ship ) status= nf90_inq_varid(ncid,'MGGTP2',varid_tmxt2)
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmxt2, nc_tmxt2, (/mrepsta/),(/nrep/))
      IF (kcdftyp==ncdf_synop) status= nf90_inq_varid(ncid,'MGGTP4',varid_tmnt1)
      IF (kcdftyp==ncdf_ship ) status= nf90_inq_varid(ncid,'MGGTP3',varid_tmnt1)
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmnt1, nc_tmnt1, (/mrepsta/),(/nrep/))
      IF (kcdftyp==ncdf_synop) status= nf90_inq_varid(ncid,'MGGTP5',varid_tmnt2)
      IF (kcdftyp==ncdf_ship ) status= nf90_inq_varid(ncid,'MGGTP4',varid_tmnt2)
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tmnt2, nc_tmnt2, (/mrepsta/),(/nrep/))
    ENDIF

    IF     (kcdftyp == ncdf_synop) THEN
      status = nf90_inq_varid (ncid,'ME     ',varid_e     )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_e    , nc_e    , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MGGTP7 ',varid_raddt )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_raddt, nc_raddt, start=(/1,mrepsta/) &
                                                         , count=(/2,nrep/))
      status = nf90_inq_varid (ncid,'MGLSR  ',varid_radgl )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_radgl, rc_radgl, start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_radgl, rc_radg2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'MDSRH  ',varid_raddf )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_raddf, rc_raddf, start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_raddf, rc_radd2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'MLWR   ',varid_radlw )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_radlw, rc_radlw, start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_radlw, rc_radl2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'NSSS   ',varid_hsnow )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_hsnow, rc_hsnow, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MPN    ',varid_pfi   )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_pfi  , rc_pfi  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NHHHN  ',varid_fi    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_fi   , nc_fi   , (/mrepsta/),(/nrep/))
    ELSEIF (kcdftyp == ncdf_buoy ) THEN
      status = nf90_inq_varid (ncid,'MRRR   ',varid_rr    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rr   , rc_rr   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MGGTP  ',varid_rrdt  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rrdt , nc_rrdt , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NFXGU  ',varid_gust  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_gust , rc_gust , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP0 ',varid_gudt  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_gudt , nc_gudt , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTISI1 ',varid_vtisi )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_vtisi, nc_vtisi, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTODB  ',varid_tbuoy )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_tbuoy, nc_tbuoy, (/mrepsta/),(/nrep/))
    ELSEIF (kcdftyp == ncdf_nat_dwd) THEN
      status = nf90_inq_varid (ncid,'MRRR   ',varid_rr    )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rr   , rc_rr   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP8 ',varid_rrdt  )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rrdt , nc_rrdt , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MRRR0  ',varid_rr2   )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rr2  , rc_rr2  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTP9 ',varid_rrdt2 )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_rrdt2, nc_rrdt2, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NFXGU  ',varid_gust  )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_gust , rc_gust , start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_gust , rc_gust2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'NGGTPG ',varid_gudt  )
      IF (status == nf90_noerr) THEN
        status= nf90_get_var (ncid, varid_gudt , nc_gudt , start=(/1,mrepsta/) &
                                                         , count=(/1,nrep/))
        status= nf90_get_var (ncid, varid_gudt , nc_gudt2, start=(/2,mrepsta/) &
                                                         , count=(/1,nrep/))
      ENDIF
      status = nf90_inq_varid (ncid,'NDNDN',varid_dd )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_dd  , nc_dd   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NFNFN',varid_ff )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_ff  , rc_ff   , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MTISI ',varid_vtisi )
      IF (status == nf90_noerr)                                                &
        status= nf90_get_var (ncid, varid_vtisi, nc_vtisi, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NGGTPF',varid_vdt )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_vdt , nc_vdt  , (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'NIX   ',varid_nix )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nix , nc_nix  , (/mrepsta/),(/nrep/))
    ELSEIF (kcdftyp == ncdf_metar) THEN
      status = nf90_inq_varid (ncid,'NFXGU  ',varid_gust  )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_gust, rc_gust, (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MHOSEN0',varid_zt    )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_zt  , rc_zt  , (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MVVVV  ',varid_vis   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_vis , rc_vis , (/mrepsta/), (/nrep/))
    ENDIF

! variables not used:
!   ! sea / water (surface) temperature (not for mobile synop !!!)
!   IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
!     status = nf90_inq_varid (ncid,'NDBSS ',varid_sstz  )
!     status = nf90_inq_varid (ncid,'MTN00 ',varid_sst   )
!   ELSEIF (kcdftyp == ncdf_buoy ) THEN
!     status = nf90_inq_varid (ncid,'MDREP ',varid_sstrp )
!     status = nf90_inq_varid (ncid,'NZNZN ',varid_sstz  )
!     status = nf90_inq_varid (ncid,'MTN00 ',varid_sst   )
!   ENDIF
!   ! direction + speed of moving platforms (data in ship, missing in buoy reps)
!   IF ((kcdftyp == ncdf_ship ) .OR. (kcdftyp == ncdf_buoy ))                  &
!     status = nf90_inq_varid (ncid,'MDS   ',varid_stadd )
!   IF (kcdftyp == ncdf_ship) status = nf90_inq_varid (ncid,'NVS  ',varid_staff)
!   IF (kcdftyp == ncdf_buoy) status = nf90_inq_varid (ncid,'MDSDS',varid_staff)
!   ! height of temperature resp. wind sensor above water surface (missing)
!   IF ((kcdftyp == ncdf_ship ) .OR. (kcdftyp == ncdf_buoy ))                  &
!      status = nf90_inq_varid(ncid,'MHAWAS ',varid_ht)
!   IF (kcdftyp == ncdf_ship ) status = nf90_inq_varid (ncid,'MHAWAS3',varid_hv)
!   IF (kcdftyp == ncdf_buoy ) status = nf90_inq_varid (ncid,'MHAWAS1',varid_hv)
!   ! height of wind sensor above ground or platform deck (missing)
!   IF (kcdftyp == ncdf_synop) status = nf90_inq_varid (ncid,'MHOSEN5',varid_zv)
!   IF (kcdftyp == ncdf_ship ) status = nf90_inq_varid (ncid,'MHOSEN5',varid_zv)
!   IF (kcdftyp == ncdf_buoy ) status = nf90_inq_varid (ncid,'MHOSEN0',varid_zv)

! get additional cloud levels (variable replication, for synop, ship, metar)
! ---------------------------
! requirement: the replication factor for the number of cloud levels is
!              optional, but if it is present then the whole section
!              (cloud type, amount, base height, significance) must be present

    IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
      status = nf90_inq_varid (ncid, 'MDREP ', varid_clev )
    ELSEIF (kcdftyp == ncdf_metar) THEN
      status = nf90_inq_varid (ncid, 'MDREP2', varid_clev )
    ENDIF

    maxlev = 0
! get maximum number of cloud levels (according to entry 'MDREP')
    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_metar)) THEN
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_clev, nc_clev, (/mrepsta/), (/nrep/))
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF (nc_clev(irep) /= imiss)  maxlev = MAX( nc_clev(irep) , maxlev )
      ENDDO
    ENDIF

! if maximum number > 0: get NetCDF variable id's of additional cloud entries
    IF ((maxlev >= 1) .AND. (     (kcdftyp == ncdf_synop)                      &
                             .OR. (kcdftyp == ncdf_ship ))) THEN
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MVTSU0',varid_clxsg)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MNH0  ',varid_clxcl)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MCC2  ',varid_clxct)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'NH0   ',varid_clxbs)
    ELSEIF ((maxlev >= 1) .AND.   (kcdftyp == ncdf_metar))  THEN
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MVTSU ',varid_clxsg)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MNH   ',varid_clxcl)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'MCC   ',varid_clxct)
      IF(status == nf90_noerr) status= nf90_inq_varid(ncid,'NH    ',varid_clxbs)
    ENDIF
    IF ((maxlev >= 1) .AND. (status /= nf90_noerr)) THEN
      yerrmsl = 'ADDITIONAL CLOUD ELEMENTS DO NOT EXIST IN ' // yfn
      IF (.NOT. lax) THEN
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ELSE
        PRINT *, yerrmsl
        maxlev = 0
      ENDIF
    ENDIF
    IF (maxlev >= 1) THEN
! get dimension length in NetCDF file and check max. number of cloud levels
      status = nf90_Inquire_Variable ( ncid, varid_clxsg, ndims=ndims          &
                                     , dimids=dimids)
      dimid_mxlv = dimids(1)
      status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
      IF (maxlev > mxlev) THEN
        yerrmsl = 'NUMBER OF CLOUD LEVELS EXCEEDS DIMENSION IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
      maxlev = MIN( maxlev , mxlev )
! 'nc_clev' is adjusted below so as not to exceed the revised value of 'maxlev'
    ENDIF

! read elements of additional cloud levels
    IF (maxlev >= 1) THEN
      ALLOCATE ( nc_clxsg (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_clxcl (maxlev,nrep) , STAT=istat )
      ALLOCATE ( nc_clxct (maxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_clxbs (maxlev,nrep) , STAT=istat )
      nsta2 (1) = 1
      ncnt2 (1) = maxlev
      nsta2 (2) = mrepsta
      ncnt2 (2) = nrep
      status = nf90_get_var (ncid, varid_clxsg, nc_clxsg, nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_clxcl, nc_clxcl, nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_clxct, nc_clxct, nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_clxbs, rc_clxbs, nsta2, ncnt2)
    ENDIF

! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

! T-min and T-max are discarded, unless time period is in the range [-1,-24] hrs
!                                               and ends at reported obs time
    IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF (      (nc_tmxt1(irep) <  -24  ) .OR.  (nc_tmxt1(irep) >  -1)       &
            .OR. ((nc_tmxt2(irep) /= imiss) .AND. (nc_tmxt2(irep) /=  0))) THEN
          rc_tmax  (irep)  =  rmiss
          nc_tmxt1 (irep)  =  imiss
        ENDIF
        IF (      (nc_tmnt1(irep) <  -24  ) .OR.  (nc_tmnt1(irep) >  -1)       &
            .OR. ((nc_tmnt2(irep) /= imiss) .AND. (nc_tmnt2(irep) /=  0))) THEN
          rc_tmin  (irep)  =  rmiss
          nc_tmnt1 (irep)  =  imiss
        ENDIF
      ENDDO
      DEALLOCATE ( nc_tmxt2 , STAT=istat )
      DEALLOCATE ( nc_tmnt2 , STAT=istat )
    ENDIF
! radiation is discarded, unless time period is -1 hrs
    IF (kcdftyp == ncdf_synop) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF (ABS( nc_raddt(1,irep) ) /= 1) THEN
          rc_radgl (irep)  =  rmiss
          rc_raddf (irep)  =  rmiss
          rc_radlw (irep)  =  rmiss
        ENDIF
        IF (ABS( nc_raddt(2,irep) ) == 1) THEN
          rc_radgl (irep)  =  rc_radg2(irep)
          rc_raddf (irep)  =  rc_radd2(irep)
          rc_radlw (irep)  =  rc_radl2(irep)
        ENDIF
      ENDDO
      DEALLOCATE ( nc_raddt , STAT=istat )
      DEALLOCATE ( rc_radg2 , STAT=istat )
      DEALLOCATE ( rc_radd2 , STAT=istat )
      DEALLOCATE ( rc_radl2 , STAT=istat )
    ENDIF
! radiation is discarded, unless time period is -1 hrs
    IF (kcdftyp == ncdf_nat_dwd) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF (nc_rrdt (irep) /= imiss)  nc_rrdt (irep) = nc_rrdt (irep) / 60
        IF (nc_rrdt2(irep) /= imiss)  nc_rrdt2(irep) = nc_rrdt2(irep) / 60
      ENDDO
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

    IF (maxlev >= 1) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
! adjust 'nc_clev'
        nc_clev(irep)  =  MIN( nc_clev(irep) , maxlev )
! determine length of integer / real buffer for each report
        iirlen (irps)  =  niscal + 3 *nc_clev(irep)
        irrlen (irps)  =  nrscal +    nc_clev(irep)
        iyrlen (irps)  =  nyscal
      ENDDO
    ELSE
      DO irps = 1 , nreproc
        iirlen (irps)  =  niscal
        irrlen (irps)  =  nrscal
        iyrlen (irps)  =  nyscal
      ENDDO
    ENDIF

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
      IF (lscatto) THEN
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_dd (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_ff (1,irep)
      ENDIF
      IF (linsitu) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nix  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_dd (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_t  (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_td (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_ff (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 4)  =  rc_gust (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 5)  =  rc_zt   (irep)
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy )) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_uuu   (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_vdt   (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_gudt  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_vtisi (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 7)  =  nc_rrdt  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 6)  =  rc_p   (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 7)  =  rc_pmsl  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 8)  =  rc_dpdt  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 9)  =  rc_rr    (irep)
      ENDIF
      IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 8)  =  nc_gudt2(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 9)  =  nc_rrdt2(irep)
        ibufsrt (iioffs+noffscal(1,1)+10)  =  nc_clev (irep)
        ibufsrt (iioffs+noffscal(1,1)+11)  =  nc_clct (irep)
        ibufsrt (iioffs+noffscal(1,1)+12)  =  nc_clsig(irep)
        ibufsrt (iioffs+noffscal(1,1)+13)  =  nc_clclm(irep)
        ibufsrt (iioffs+noffscal(1,1)+14)  =  nc_ccl  (irep)
        ibufsrt (iioffs+noffscal(1,1)+15)  =  nc_ccm  (irep)
        ibufsrt (iioffs+noffscal(1,1)+16)  =  nc_cch  (irep)
        ibufsrt (iioffs+noffscal(1,1)+17)  =  nc_ww   (irep)
        ibufsrt (iioffs+noffscal(1,1)+18)  =  nc_wdt  (irep)
        ibufsrt (iioffs+noffscal(1,1)+19)  =  nc_w1   (irep)
        ibufsrt (iioffs+noffscal(1,1)+20)  =  nc_w2   (irep)
        ibufsrt (iioffs+noffscal(1,1)+21)  =  nc_tmxt1(irep)
        ibufsrt (iioffs+noffscal(1,1)+22)  =  nc_tmnt1(irep)
        rbufsrt (iroffs+noffscal(2,1)+10)  =  rc_gust2(irep)
        rbufsrt (iroffs+noffscal(2,1)+11)  =  rc_rr2  (irep)
        rbufsrt (iroffs+noffscal(2,1)+12)  =  rc_vis  (irep)
        rbufsrt (iroffs+noffscal(2,1)+13)  =  rc_rr24 (irep)
        rbufsrt (iroffs+noffscal(2,1)+14)  =  rc_cbase(irep)
        rbufsrt (iroffs+noffscal(2,1)+15)  =  rc_tmax (irep)
        rbufsrt (iroffs+noffscal(2,1)+16)  =  rc_tmin (irep)
      ENDIF
      IF     (kcdftyp == ncdf_synop) THEN
        ibufsrt (iioffs+noffscal(1,1)+23)  =  nc_e    (irep)
        ibufsrt (iioffs+noffscal(1,1)+24)  =  nc_fi   (irep)
        rbufsrt (iroffs+noffscal(2,1)+17)  =  rc_radgl(irep)
        rbufsrt (iroffs+noffscal(2,1)+18)  =  rc_raddf(irep)
        rbufsrt (iroffs+noffscal(2,1)+19)  =  rc_radlw(irep)
        rbufsrt (iroffs+noffscal(2,1)+20)  =  rc_hsnow(irep)
        rbufsrt (iroffs+noffscal(2,1)+21)  =  rc_pfi  (irep)
      ELSEIF (kcdftyp == ncdf_buoy ) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 8)  =  nc_tbuoy(irep)
      ELSEIF (kcdftyp == ncdf_metar) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_clev (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 6)  =  rc_vis  (irep)
      ELSEIF (kcdftyp == ncdf_nat_dwd) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nix  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_rrdt (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_rrdt2(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_gudt (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_gudt2(irep)
        ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_dd (1,irep)
        ibufsrt (iioffs+noffscal(1,1)+ 7)  =  nc_vdt  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 8)  =  nc_vtisi(irep)
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_rr   (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_rr2  (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_gust (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 4)  =  rc_gust2(irep)
        rbufsrt (iroffs+noffscal(2,1)+ 5)  =  rc_ff (1,irep)
      ENDIF

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

! de-allocate reading arrays
    IF (lscatto) THEN
      DEALLOCATE ( rc_dd    , STAT=istat )
      DEALLOCATE ( rc_ff    , STAT=istat )
      DEALLOCATE ( rc_ff4   , STAT=istat )
      DEALLOCATE ( rc_dd4   , STAT=istat )
      DEALLOCATE ( nc_niswv , STAT=istat )
    ENDIF
    IF (linsitu) THEN
      DEALLOCATE ( nc_nix   , STAT=istat )
      DEALLOCATE ( rc_t     , STAT=istat )
      DEALLOCATE ( rc_td    , STAT=istat )
      DEALLOCATE ( nc_dd    , STAT=istat )
      DEALLOCATE ( rc_ff    , STAT=istat )
      DEALLOCATE ( rc_gust  , STAT=istat )
      DEALLOCATE ( rc_zt    , STAT=istat )
    ENDIF
    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_buoy )) THEN
      DEALLOCATE ( rc_p     , STAT=istat )
      DEALLOCATE ( rc_pmsl  , STAT=istat )
      DEALLOCATE ( rc_dpdt  , STAT=istat )
      DEALLOCATE ( nc_uuu   , STAT=istat )
      DEALLOCATE ( rc_rr    , STAT=istat )
      DEALLOCATE ( nc_rrdt  , STAT=istat )
      DEALLOCATE ( nc_vdt   , STAT=istat )
      DEALLOCATE ( nc_gudt  , STAT=istat )
      DEALLOCATE ( nc_vtisi , STAT=istat )
    ENDIF
    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_metar)) THEN
      DEALLOCATE ( rc_vis   , STAT=istat )
!     DEALLOCATE ( nc_clev  , STAT=istat )
    ENDIF
    IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
      DEALLOCATE ( rc_gust2 , STAT=istat )
      DEALLOCATE ( nc_gudt2 , STAT=istat )
      DEALLOCATE ( rc_rr24  , STAT=istat )
      DEALLOCATE ( rc_rr2   , STAT=istat )
      DEALLOCATE ( nc_rrdt2 , STAT=istat )
      DEALLOCATE ( nc_clct  , STAT=istat )
      DEALLOCATE ( nc_clsig , STAT=istat )
      DEALLOCATE ( nc_clclm , STAT=istat )
      DEALLOCATE ( rc_cbase , STAT=istat )
      DEALLOCATE ( nc_ccl   , STAT=istat )
      DEALLOCATE ( nc_ccm   , STAT=istat )
      DEALLOCATE ( nc_cch   , STAT=istat )
      DEALLOCATE ( nc_ww    , STAT=istat )
      DEALLOCATE ( nc_wdt   , STAT=istat )
      DEALLOCATE ( nc_w1    , STAT=istat )
      DEALLOCATE ( nc_w2    , STAT=istat )
      DEALLOCATE ( rc_tmax  , STAT=istat )
      DEALLOCATE ( rc_tmin  , STAT=istat )
      DEALLOCATE ( nc_tmxt1 , STAT=istat )
      DEALLOCATE ( nc_tmnt1 , STAT=istat )
    ENDIF
    IF     (kcdftyp == ncdf_synop) THEN
      DEALLOCATE ( nc_e     , STAT=istat )
      DEALLOCATE ( rc_radgl , STAT=istat )
      DEALLOCATE ( rc_raddf , STAT=istat )
      DEALLOCATE ( rc_radlw , STAT=istat )
      DEALLOCATE ( rc_hsnow , STAT=istat )
      DEALLOCATE ( rc_pfi   , STAT=istat )
      DEALLOCATE ( nc_fi    , STAT=istat )
    ELSEIF (kcdftyp == ncdf_buoy ) THEN
      DEALLOCATE ( nc_tbuoy , STAT=istat )
    ELSEIF (kcdftyp == ncdf_nat_dwd ) THEN
      DEALLOCATE ( rc_rr    , STAT=istat )
      DEALLOCATE ( rc_rr2   , STAT=istat )
      DEALLOCATE ( rc_gust  , STAT=istat )
      DEALLOCATE ( rc_gust2 , STAT=istat )
      DEALLOCATE ( rc_ff    , STAT=istat )
      DEALLOCATE ( nc_nix   , STAT=istat )
      DEALLOCATE ( nc_rrdt  , STAT=istat )
      DEALLOCATE ( nc_rrdt2 , STAT=istat )
      DEALLOCATE ( nc_gudt  , STAT=istat )
      DEALLOCATE ( nc_gudt2 , STAT=istat )
      DEALLOCATE ( nc_dd    , STAT=istat )
      DEALLOCATE ( nc_vdt   , STAT=istat )
      DEALLOCATE ( nc_vtisi , STAT=istat )
    ENDIF

! fill in the variable-level cloud info into the long buffer arrays 'Xbufsrt'
! ---------------------------------------------------------------------------

    IF (maxlev >= 1) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev = nc_clev(irep)
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_clxsg(ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_clxcl(ilev,irep)
          ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_clxct(ilev,irep)
          rbufsrt (iroffs+nrscal       +ilev)  =  rc_clxbs(ilev,irep)
        ENDDO
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO

      DEALLOCATE ( nc_clxsg , STAT=istat )
      DEALLOCATE ( nc_clxcl , STAT=istat )
      DEALLOCATE ( nc_clxct , STAT=istat )
      DEALLOCATE ( rc_clxbs , STAT=istat )
    ENDIF
    DEALLOCATE ( nc_clev  , STAT=istat )
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

    CALL obs_cdf_store_singlev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly   &
                               , noffscal,          ibufloc, rbufloc, ybufloc  &
                               , nsgnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_singlev ( kcdftyp , nrepl  , nlenli , nlenlr , nlenly   &
                               , noffscal,          ibufsrt, rbufsrt, ybufsrt  &
                               , nsgnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_surface
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_surface


!===============================================================================
!+ Module procedure in "src_obs_cdfin_sing" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_aircraft ( min_sta , min_end , ilcdf                   &
                                 , nodrnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_sing" organizes the reading,
!   pre-selection, distribution and storage in ODR of aircraft (AMDAR (single-
!   level and multi-level), ACARS (USA ( ARINC Centre 56); UKMO and Canada)
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

  LOGICAL                  ::  &
    lamdarml       ,& ! multi-level AMDAR (with or w/o lat-lon profile)
    lamdar         ,& ! (multi- or single-level) AMDAR
    lacars            ! ACARS

  INTEGER        ::  &
    kcdftyp        ,& ! type of NetCDF input file format (observation type):
                      !  = ncdf_amdar     : AMDAR single-level
                      !  = ncdf_acars     : ACARS single-level
                      !  = ncdf_acars_uk  : ACARS UK + Canada
                      !  = ncdf_acars_us  : ACARS USA (with humidity)
                      !  = ncdf_amdar_ml  : AMDAR multi-level
                      !  = ncdf_amdar_vp  : AMDAR vertical profile
                      !  = ncdf_modes     : MODE-S single-level KNMI format
                      !  = ncdf_modes_acr : MODE-S single-level ACARS format
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
    ilev           ,& ! loop index over vertical levels
    nlev           ,& ! number of vertical levels
    maxlev         ,& ! max number of vertical levels in all processd reports
    mxlev          ,& ! length of vertical level dimension in NetCDF file
    ndims          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    dimid_mxlv     ,& ! dimension ID for vertical levels in NetCDF file
    status         ,& ! NetCDF status variable
    muuu_type      ,& ! type of variable MUUU in NetCDF file
    nsta2 (2)      ,& ! start indices       \  of the reports to be read
    ncnt2 (2)      ,& ! number of elements  /  in the external NetCDF
    istat  , ierr     ! error indicators

  INTEGER        ::  &          ! variable ID's in NetCDF file for:
!                              !   all file types (synop, ship, buoy, metar):
    varid_nix               ,& ! station type  (0:auto,1:manned,2:hybrid,3:miss)
    varid_nlev              ,& ! number of vertical levels
    varid_lat  , varid_lon  ,& ! latitude             , longitude
    varid_phase, varid_rolla,& ! flight phase [008004], roll angle [002064]
    varid_dd   , varid_ff   ,& ! wind direction       , wind speed
    varid_t    , varid_td   ,& ! dry bulb temperature , dew point temperature
    varid_rh   , varid_qx   ,& ! relative humidity    , mixing ratio
    varid_qxq  , varid_tq   ,& ! mixing ration quality, temperat. obs precision
    varid_z_i  , varid_z_r  ,& ! height (integer resp. real)
    varid_p    , varid_phasd,& ! pressure , detailed phase of flight [008009]
    varid_turb , varid_vgust   ! turbulence [011031]  , max. vertical gust

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
! CHARACTER (LEN=30)       :: &
!   yerr              ! error message
  CHARACTER (LEN=70)       :: &
    yerrmsl           ! error message
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
! Begin Subroutine obs_cdf_read_aircraft
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_aircraft'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

  lamdarml =      (kcdftyp == ncdf_amdar_ml) .OR. (kcdftyp == ncdf_amdar_vp)
  lamdar   =      (kcdftyp == ncdf_amdar   ) .OR. (lamdarml)
  lacars   =      (kcdftyp == ncdf_acars   ) .OR. (kcdftyp == ncdf_acars_uk)   &
                                             .OR. (kcdftyp == ncdf_acars_us)   &
                                             .OR. (kcdftyp == ncdf_modes_acr)

! PRINT *,yroutine, kcdftyp, min_sta,min_end, lamdarml, lamdar, lacars

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

  IF (kcdftyp == ncdf_amdar) THEN
    noffscal (1,2)  =  noffscal(1,1) + 3 + 2
    noffscal (2,2)  =  noffscal(2,1) + 2 + 5
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (lacars) THEN
    noffscal (1,2)  =  noffscal(1,1) + 5 + 2
    noffscal (2,2)  =  noffscal(2,1) + 3 + 5
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (kcdftyp == ncdf_modes) THEN
    noffscal (1,2)  =  noffscal(1,1) + 3 + 2
    noffscal (2,2)  =  noffscal(2,1) + 1 + 5
    noffscal (3,2)  =  noffscal(3,1) + 0
  ELSEIF (lamdarml) THEN
    noffscal (1,2)  =  noffscal(1,1) + 2
    noffscal (2,2)  =  noffscal(2,1) + 0
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

! 'irnode' has been allocated in 'obs_cdf_read_comhead' and is not used any more
    DEALLOCATE ( irproc , STAT=istat )

!-------------------------------------------------------------------------------
! Section 2: Get entries specific to observation type (AMDAR, ACARS, etc.)
!-------------------------------------------------------------------------------

! get dimension of fields (notably the vertical dimension for multi-level AMDAR)
! -----------------------

    ALLOCATE ( nc_nlev         (nrep) , STAT=istat )

    IF (lamdarml) THEN
                              status = nf90_inq_varid (ncid,'MDREP',varid_nlev )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nlev, nc_nlev, (/mrepsta/), (/nrep/))
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'NUMBER OF VERTICAL LEVELS DOES NOT EXIST IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
! get maximum number of vertical levels (according to entry 'MDREP')
      maxlev = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        maxlev = MAX( nc_nlev(irep) , maxlev )
      ENDDO
      IF (maxlev >= 1) THEN
! get and check (length of) vertical dimension in NetCDF file
                              status = nf90_inq_varid (ncid,'NDNDN',varid_dd   )
        IF (status == nf90_noerr)                                              &
          status = nf90_Inquire_Variable ( ncid, varid_dd, ndims=ndims         &
                                         , dimids=dimids)
        IF ((status /= nf90_noerr) .OR. (ndims /= 2)) THEN
          WRITE( yerrmsl,'("NO VALID DIMENSION FOR WIND (",I2,I6,") IN ",A)')  &
                 ndims, status,  yfn
          CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
        ENDIF
        dimid_mxlv = dimids(1)
        status = nf90_Inquire_Dimension (ncid, dimid_mxlv, len=mxlev)
        IF ((maxlev > mxlev) .OR. (status /= nf90_noerr)) THEN
          yerrmsl = 'NUMBER OF LEVELS EXCEEDS DIMENSION IN ' // yfn
          CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
        ENDIF
!  if multi-level report: the second dimension is the report dimension
        nsta2 (1) = 1
        ncnt2 (1) = maxlev
        nsta2 (2) = mrepsta
        ncnt2 (2) = nrep
      ENDIF

    ELSE  !  if single-level report: the first dimension is the report dimension
      nc_nlev (:) = 1
      maxlev      = 1
      nsta2 (1) = mrepsta
      ncnt2 (1) = nrep
      nsta2 (2) = 1
      ncnt2 (2) = 1
    ENDIF

! allocate arrays required for reading
! ------------------------------------
!   PRINT *,'varid3 ', ncid, maxlev, ncdf_temp, kcdftyp, nrep
    mxlev = MAX( maxlev , 1 )

    ! initialise mandatory elements
    ALLOCATE ( nc_dd    (mxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_ff    (mxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_t     (mxlev,nrep) , STAT=istat )
    ALLOCATE ( nc_phase       (nrep) , STAT=istat )
    ! (initialize, since it is not mandatory / present for MODE-S)
    nc_phase (:)   = imiss

    ! initialise elements which are mandatory only for AMDAR or if not (lax)
    ALLOCATE ( nc_rolla (mxlev,nrep) , STAT=istat )
    ALLOCATE ( nc_z     (mxlev,nrep) , STAT=istat )
    ALLOCATE ( rc_td    (mxlev,nrep) , STAT=istat )
    nc_rolla (:,:) = imiss
    nc_z     (:,:) = imiss
    rc_td    (:,:) = rmiss


    ! initialise optional elements
    IF     (lamdarml) THEN
      ALLOCATE ( rc_lat0  (mxlev,nrep) , STAT=istat )
      ALLOCATE ( rc_lon0  (mxlev,nrep) , STAT=istat )
      rc_lat0 (:,:) = rmiss
      rc_lon0 (:,:) = rmiss
    ENDIF
    IF ((kcdftyp == ncdf_amdar) .OR. (lacars) .OR. (kcdftyp == ncdf_modes)) THEN
      ALLOCATE ( nc_rh          (nrep) , STAT=istat )
      ALLOCATE ( nc_turb        (nrep) , STAT=istat )
      ALLOCATE ( nc_qxq         (nrep) , STAT=istat )
      ALLOCATE ( rc_qx          (nrep) , STAT=istat )
      ALLOCATE ( rc_p         (1,nrep) , STAT=istat )
      ALLOCATE ( rc_rh          (nrep) , STAT=istat )
      ALLOCATE ( rc_vgust       (nrep) , STAT=istat )
      ALLOCATE ( nc_phasd       (nrep) , STAT=istat )
      nc_rh    (:) = imiss
      nc_turb  (:) = imiss
      nc_qxq   (:) = imiss
      rc_qx    (:) = rmiss
      ! (rc_p is defined as a 2-dimensional array)
      rc_p   (1,:) = rmiss
      rc_rh    (:) = rmiss
      rc_vgust (:) = rmiss
      nc_phasd (:) = imiss
      IF (lacars) THEN
        ALLOCATE ( nc_nix         (nrep) , STAT=istat )
        ALLOCATE ( rc_tq          (nrep) , STAT=istat )
        nc_nix   (:) = imiss
        rc_tq    (:) = rmiss
      ENDIF
    ENDIF
    IF (     (lacars) .OR. (kcdftyp == ncdf_modes)                             &
        .OR. ((lax) .AND. (.NOT. lamdarml))) THEN
      ALLOCATE ( rc_z           (nrep) , STAT=istat )
      rc_z     (:) = rmiss
    ENDIF

! get mandatory (!) variable ID's in NetCDF file and get data
! -----------------------------------------------------------

    IF (maxlev >= 1) THEN
                              status = nf90_inq_varid (ncid,'NDNDN',varid_dd )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                              status = nf90_inq_varid (ncid,'NDD  ',varid_dd )
      IF (status == nf90_noerr) THEN
                              status = nf90_inq_varid (ncid,'NFNFN',varid_ff )
        IF ((status /= nf90_noerr) .AND. (lax))                                &
                              status = nf90_inq_varid (ncid,'NFF  ',varid_ff )
      ENDIF
      IF (status == nf90_noerr) THEN
                              status = nf90_inq_varid (ncid,'MTDBT',varid_t  )
        IF ((status /= nf90_noerr) .AND. ((lax) .OR. (lacars)))                &
                              status = nf90_inq_varid (ncid,'MTN  ',varid_t  )
      ENDIF
      IF ((.NOT. lax) .AND. (lamdar) .AND. (status == nf90_noerr))             &
                              status = nf90_inq_varid (ncid,'MTDNH',varid_td )
!       IF ((status /= nf90_noerr) .AND. (lax))                                &
!                             status = nf90_inq_varid (ncid,'MTDN ',varid_td )
!     ENDIF
      IF ((.NOT. lax) .AND. (lamdar .OR. lacars) .AND. (status == nf90_noerr)) &
                              status = nf90_inq_varid (ncid,'MPHAI',varid_phase) ! WMO 008004
      IF (      (((.NOT. lax) .AND. (kcdftyp == ncdf_amdar)) .OR. (lamdarml))  &
          .AND. (status == nf90_noerr))                                        &
                              status = nf90_inq_varid (ncid,'NFLEV',varid_z_i)
      IF ((.NOT. lax) .AND. (status == nf90_noerr))                            &
                              status = nf90_inq_varid (ncid,'MQARA',varid_rolla) ! WMO 002064
      IF (status /= nf90_noerr) THEN
        yerrmsl = 'STANDARD INFO DOES NOT EXIST IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
      status = nf90_get_var (ncid, varid_dd   , nc_dd   , nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_ff   , rc_ff   , nsta2, ncnt2)
      status = nf90_get_var (ncid, varid_t    , rc_t    , nsta2, ncnt2)
      IF ((.NOT. lax) .AND. (lamdar))                                          &
        status = nf90_get_var (ncid, varid_td , rc_td   , nsta2, ncnt2)
      IF ((.NOT. lax) .AND. (lamdar .OR. lacars))                              &
        status = nf90_get_var (ncid, varid_phase,nc_phase, (/mrepsta/),(/nrep/))
      IF (((.NOT. lax) .AND. (kcdftyp == ncdf_amdar)) .OR. (lamdarml))         &
        status = nf90_get_var (ncid, varid_z_i, nc_z    , nsta2, ncnt2)
      IF  (.NOT. lax)                                                          &
        status = nf90_get_var (ncid, varid_rolla, nc_rolla, nsta2, ncnt2)
    ENDIF

    !   for ACARS or MODE-S, or if (lax .and. single-level AMDAR),
    !   at least one of the vertical coordinates must be present
    !     (but 'MIAA' and 'NFLEV' must not be present both in the same file !)
!   IF ((maxlev >= 1) .AND. (.NOT. lamdarml) .AND. ((lacars) .OR. (lax))) THEN
!   IF ((maxlev >= 1) .AND. ((lacars) .OR. ((lax) .AND. (.NOT. lamdarml)))) THEN
    IF ((maxlev >= 1) .AND. (     (lacars) .OR. (kcdftyp == ncdf_modes)        &
                             .OR. ((lax) .AND. (.NOT. lamdarml)))) THEN
                                status = nf90_inq_varid (ncid,'MIAA ',varid_z_i)
      IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'NFLEV',varid_z_i)
      IF (status /= nf90_noerr) varid_z_i    = imiss
                                status = nf90_inq_varid (ncid,'MHHH ',varid_z_r)
      IF (status /= nf90_noerr) varid_z_r    = imiss
                                status = nf90_inq_varid (ncid,'MPN  ',varid_p  )
      IF (status /= nf90_noerr) varid_p      = imiss
      IF (      (varid_z_i == imiss) .AND. (varid_z_r == imiss)                &
          .AND. (varid_p   == imiss)) THEN
        yerrmsl = 'NO VERTICAL COORDINATE EXISTS IN ' // yfn
        CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
      ENDIF
      IF (varid_z_i   /= imiss)                                                &
        status = nf90_get_var (ncid, varid_z_i , nc_z , nsta2, ncnt2)
      IF (varid_z_r   /= imiss)                                                &
        status = nf90_get_var (ncid, varid_z_r , rc_z , (/mrepsta/), (/nrep/))
      IF (varid_p     /= imiss)                                                &
        status = nf90_get_var (ncid, varid_p   , rc_p , (/mrepsta/), (/nrep/))
    ENDIF

! get optional (!) variable ID's in NetCDF file and data
! ------------------------------------------------------
!     (currently for ncdf_acars_uk, 'MB', 'NMDEWX', and 'MHHH' exist,
!                for ncdf_acars_us, 'MUUU', 'NIX', 'MPN', 'MIAA', 'MMRQ' exist,
!            and for both types,    'MMIXR' and 'MPOTO' exist)

    IF     ((maxlev >= 1) .AND. (lax)) THEN
      status = nf90_inq_varid (ncid,'MQARA ',varid_rolla)           ! WMO 002064
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_rolla, nc_rolla, nsta2, ncnt2)
    ENDIF

    IF ((maxlev >= 1) .AND. ((lacars) .OR. ((lamdar) .AND. (lax)))) THEN
                                status = nf90_inq_varid (ncid,'MTDNH',varid_td )
      IF (status /= nf90_noerr) status = nf90_inq_varid (ncid,'MTDN ',varid_td )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_td  , rc_td  , nsta2, ncnt2)
    ENDIF

    IF ((maxlev >= 1) .AND. (lax) .AND. ((lamdar) .OR. (lacars))) THEN
      status = nf90_inq_varid (ncid,'MPHAI',varid_phase)            ! WMO 008004
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_phase,nc_phase, (/mrepsta/),(/nrep/))
    ENDIF

    ! all single-level (AMDAR + ACARS + MODE-S)
    IF     ((maxlev >= 1) .AND. (     (kcdftyp == ncdf_amdar) .OR. (lacars)    &
                                 .OR. (kcdftyp == ncdf_modes))) THEN
      status = nf90_inq_varid (ncid,'NDEPF ',varid_phasd)           ! WMO 008009
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_phasd,nc_phasd, (/mrepsta/),(/nrep/))
      status = nf90_inq_varid (ncid,'MMIXR ',varid_qx    )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_qx  , rc_qx  , (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MUUU  ',varid_rh    )
      IF (status == nf90_noerr) THEN
        ! for MODE-S, MUUU may be real instead of integer
        status = nf90_inquire_variable (ncid, varid_rh, xtype=muuu_type)
        IF ((status == nf90_noerr) .AND. (muuu_type == NF90_FLOAT)) THEN
          status = nf90_get_var (ncid, varid_rh , rc_rh , (/mrepsta/), (/nrep/))
        ELSEIF ((status == nf90_noerr) .AND. (muuu_type == NF90_INT)) THEN
          status = nf90_get_var (ncid, varid_rh , nc_rh , (/mrepsta/), (/nrep/))
        ENDIF
      ENDIF
    ENDIF

    ! single-level AMDAR or ACARS
    IF     ((maxlev >= 1) .AND. ((kcdftyp == ncdf_amdar) .OR. (lacars))) THEN
      status = nf90_inq_varid (ncid,'MB    ',varid_turb  )          ! WMO 011031
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_turb, nc_turb, (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'NMDEWX',varid_vgust )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_vgust,rc_vgust,(/mrepsta/), (/nrep/))
    ENDIF

    ! ACARS or MODE-S
    IF     ((maxlev >= 1) .AND. ((lacars) .OR. (kcdftyp == ncdf_modes))) THEN
      status = nf90_inq_varid (ncid,'MMRQ  ',varid_qxq   )          ! WMO 033026
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_qxq , nc_qxq , (/mrepsta/), (/nrep/))
    ENDIF

    ! ACARS
    IF     ((maxlev >= 1) .AND. (lacars)) THEN
      status = nf90_inq_varid (ncid,'NIX   ',varid_nix   )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_nix , nc_nix , (/mrepsta/), (/nrep/))
      status = nf90_inq_varid (ncid,'MPOTO ',varid_tq    )
      IF (status == nf90_noerr)                                                &
        status = nf90_get_var (ncid, varid_tq  , rc_tq  , (/mrepsta/), (/nrep/))
      ! Td: for AMDAR, it is mandatory and read further above
!     status = nf90_inq_varid (ncid,'MTDNH ',varid_td    )
!     IF (status == nf90_noerr)                                                &
!       status = nf90_get_var (ncid, varid_td  , rc_td  , nsta2, ncnt2)
    ENDIF

    ! multi-level AMDAR :  latitude / longitude
    IF     ((maxlev >= 1) .AND. (lamdarml)) THEN
                                status= nf90_inq_varid (ncid,'MLAH0',varid_lat )
      IF (status == nf90_noerr) status= nf90_inq_varid (ncid,'MLOH0',varid_lon )
      IF (status == nf90_noerr) THEN
        status = nf90_get_var (ncid, varid_lat  , rc_lat0 , nsta2, ncnt2)
        status = nf90_get_var (ncid, varid_lon  , rc_lon0 , nsta2, ncnt2)
      ENDIF
    ENDIF

! preliminary cheap evaluation steps (to reduce number of elements
! ----------------------------------  to be sent to other nodes)

! transfer from detailed code (WMO Table 008009) to ordinary code (Table 008004)
! for phase of flight, if ordinary code is not level flight, asc. or descending
    IF     ((maxlev >= 1) .AND. (     (kcdftyp == ncdf_amdar) .OR. (lacars)    &
                                 .OR. (kcdftyp == ncdf_modes))) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        ! do not change value of phase if already level flight, asc. or descend.
        IF ((nc_phase(irep) <= 2) .OR. (nc_phase(irep) >= 7)) THEN
          IF ((nc_phasd(irep) >=  2) .AND. (nc_phasd(irep) <=  6)) THEN
            !    (3,4 level flight, 5 ascending, 6 descending)
            nc_phase(irep) = nc_phasd(irep)
          ELSEIF ((nc_phasd(irep) ==  7) .OR. (nc_phasd(irep) ==  9)) THEN
            !    (ascending)
            nc_phase(irep) = 5
          ELSEIF ((nc_phasd(irep) == 11) .OR. (nc_phasd(irep) == 13)) THEN
            !    (descending)
            nc_phase(irep) = 6
          ELSEIF (nc_phase(irep) /= 2) THEN
            IF (     (nc_phasd(irep) ==  0) .OR. (nc_phasd(irep) ==  1)        &
                .OR. (nc_phasd(irep) ==  8) .OR. (nc_phasd(irep) == 10)        &
                .OR. (nc_phasd(irep) == 12) .OR. (nc_phasd(irep) == 14)) THEN
              !  (unsteady)
              nc_phase(irep) = 2
            ELSE
              !  (missing)
              nc_phase(irep) = 7
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!NEC$ ivdep
!NEC$ nomove
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF ((nc_rh(irep) /= imiss) .AND. (ABS(rc_rh(irep)) >= rmisschk))       &
          rc_rh (irep)  =  REAL( nc_rh(irep) , wp )
      ENDDO
    ELSEIF (maxlev >= 1) THEN
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF ((nc_phase(irep) <= 1) .OR. (nc_phase(irep) > 7))  nc_phase(irep) = 7
      ENDDO
    ENDIF

! homogenise vertical height to be integer
    IF ((maxlev >= 1) .AND. (     (lacars) .OR. (kcdftyp == ncdf_modes)        &
                             .OR. ((lax) .AND. (.NOT. lamdarml)))) THEN
      !  (if (.NOT. lamdarml) then maxlev = 1 currently)
      ilev = 1
!NEC$ ivdep
!NEC$ nomove
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        IF (      (nc_z(ilev,irep) == imiss)                                   &
            .AND. (ABS(rc_z(irep)) <  rmisschk)                                &
            .AND. (ABS(rc_z(irep)) <  REAL( imdi, wp )))                       &
          nc_z (ilev,irep)  =  NINT( rc_z(irep) )
      ENDDO
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

    IF ((maxlev >= 1) .AND. (lamdarml)) THEN
!     PRINT *,' AMDARML nlev2 ',maxlev, nc_nlev
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
! adjust 'nc_nlev'
        nc_nlev(irep)  =  MIN( nc_nlev(irep) , maxlev )
! determine length of integer / real buffer for each report
        iirlen (irps)  =  niscal + 3 *nc_nlev(irep)
        irrlen (irps)  =  nrscal + 5 *nc_nlev(irep)
        iyrlen (irps)  =  nyscal
      ENDDO
    ELSE
      DO irps = 1 , nreproc
        iirlen (irps)  =  niscal
        irrlen (irps)  =  nrscal
        iyrlen (irps)  =  nyscal
      ENDDO
    ENDIF

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
    ibufsrt (:) = imiss
    rbufsrt (:) = rmiss

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

!     PRINT *,'nc_rh2 ',irps, irep, kcdftyp, lamdar, lacars, nc_rh(irep)
      IF (     (kcdftyp == ncdf_amdar) .OR. (lacars)                           &
          .OR. (kcdftyp == ncdf_modes)) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_phase  (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_rolla(1,irep)
        ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_dd   (1,irep)
        ibufsrt (iioffs+noffscal(1,1)+ 4)  =  nc_z    (1,irep)
!       ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_rh     (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_ff   (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_t    (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_td   (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 4)  =  rc_qx     (irep)
        rbufsrt (iroffs+noffscal(2,1)+ 5)  =  rc_p    (1,irep)
        rbufsrt (iroffs+noffscal(2,1)+ 6)  =  rc_rh     (irep)
        IF ((kcdftyp == ncdf_amdar) .OR. (lacars)) THEN
          ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_turb   (irep)
          rbufsrt (iroffs+noffscal(2,1)+ 7)  =  rc_vgust  (irep)
        ENDIF
        IF (lacars) THEN
          ibufsrt (iioffs+noffscal(1,1)+ 6)  =  nc_qxq    (irep)
          rbufsrt (iroffs+noffscal(2,1)+ 8)  =  rc_tq     (irep)
          ibufsrt (iioffs+noffscal(1,1)+ 7)  =  nc_nix    (irep)
        ENDIF
        IF (kcdftyp == ncdf_modes) THEN
          ibufsrt (iioffs+noffscal(1,1)+ 5)  =  nc_qxq    (irep)
        ENDIF
      ELSEIF (lamdarml) THEN
        ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_nlev   (irep)
        ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_phase  (irep)
      ENDIF

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

! fill the multi-level info into the long buffer arrays 'Xbufsrt'
! ---------------------------------------------------------------

    IF ((maxlev >= 1) .AND. (lamdarml)) THEN
      iioffs = 0
      iroffs = 0
      DO irps = 1 , nreproc
        irep  =  irsort(irps)
        nlev = nc_nlev(irep)
!       PRINT *,'AMDARML nlev4 ',irps, nlev, irep
        DO ilev = 1 , nlev
          ibufsrt (iioffs+niscal       +ilev)  =  nc_rolla(ilev,irep)
          ibufsrt (iioffs+niscal+  nlev+ilev)  =  nc_dd   (ilev,irep)
          ibufsrt (iioffs+niscal+2*nlev+ilev)  =  nc_z    (ilev,irep)
          rbufsrt (iroffs+nrscal       +ilev)  =  rc_lat0 (ilev,irep)
          rbufsrt (iroffs+nrscal+  nlev+ilev)  =  rc_lon0 (ilev,irep)
          rbufsrt (iroffs+nrscal+2*nlev+ilev)  =  rc_ff   (ilev,irep)
          rbufsrt (iroffs+nrscal+3*nlev+ilev)  =  rc_t    (ilev,irep)
          rbufsrt (iroffs+nrscal+4*nlev+ilev)  =  rc_td   (ilev,irep)
        ENDDO
        iioffs  =  iioffs + iirlen(irps)
        iroffs  =  iroffs + irrlen(irps)
      ENDDO
    ENDIF
!   PRINT '("buf0a ",8F10.2)' , (rbufsrt(istat),istat=1,8)
!   PRINT '("buf0b ",8F10.2)' , (rbufsrt(istat),istat=9,16)

! de-allocate reading arrays
    DEALLOCATE ( nc_phase  , STAT=istat )
    DEALLOCATE ( nc_dd     , STAT=istat )
    DEALLOCATE ( nc_rolla  , STAT=istat )
    DEALLOCATE ( nc_z      , STAT=istat )
    DEALLOCATE ( rc_ff     , STAT=istat )
    DEALLOCATE ( rc_t      , STAT=istat )
    DEALLOCATE ( rc_td     , STAT=istat )
    IF (lamdarml) THEN
      DEALLOCATE ( rc_lat0   , STAT=istat )
      DEALLOCATE ( rc_lon0   , STAT=istat )
    ENDIF
    IF ((kcdftyp == ncdf_amdar) .OR. (lacars) .OR. (kcdftyp == ncdf_modes)) THEN
      DEALLOCATE ( nc_rh     , STAT=istat )
      DEALLOCATE ( nc_turb   , STAT=istat )
      DEALLOCATE ( nc_qxq    , STAT=istat )
      DEALLOCATE ( rc_qx     , STAT=istat )
      DEALLOCATE ( rc_p      , STAT=istat )
      DEALLOCATE ( rc_rh     , STAT=istat )
      DEALLOCATE ( rc_vgust  , STAT=istat )
      DEALLOCATE ( nc_phasd  , STAT=istat )
      IF (lacars) THEN
        DEALLOCATE ( nc_nix    , STAT=istat )
        DEALLOCATE ( rc_tq     , STAT=istat )
      ENDIF
    ENDIF
    IF (     (lacars) .OR. (kcdftyp == ncdf_modes)                             &
        .OR. ((lax) .AND. (.NOT. lamdarml)))                                   &
      DEALLOCATE ( rc_z      , STAT=istat )

    DEALLOCATE ( nc_nlev  , STAT=istat )

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

    IF ((kcdftyp == ncdf_amdar) .OR. (lacars) .OR. (kcdftyp == ncdf_modes)) THEN

      CALL obs_cdf_store_singlev ( kcdftyp , nrepl , nlenli , nlenlr , nlenly  &
                                 , noffscal,         ibufloc, rbufloc, ybufloc &
                                 , nodrnew(2), nexceed(2) )
!     ==========================
    ELSE

      CALL obs_cdf_store_multilev( kcdftyp , nrepl , nlenli , nlenlr , nlenly  &
                                 , noffscal,         ibufloc, rbufloc, ybufloc &
                                 , nodrnew , nexceed )
!     ===========================
    ENDIF


    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    IF ((kcdftyp == ncdf_amdar) .OR. (lacars) .OR. (kcdftyp == ncdf_modes)) THEN

      CALL obs_cdf_store_singlev ( kcdftyp , nrepl , nlenli , nlenlr , nlenly  &
                                 , noffscal,         ibufsrt, rbufsrt, ybufsrt &
                                 , nodrnew(2), nexceed(2) )
!     ==========================
    ELSE

      CALL obs_cdf_store_multilev( kcdftyp , nrepl , nlenli , nlenlr , nlenly  &
                                 , noffscal,         ibufsrt, rbufsrt, ybufsrt &
                                 , nodrnew , nexceed )
!     ===========================
    ENDIF

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF


!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_aircraft
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_aircraft


!===============================================================================
!+ Module procedure in "src_obs_cdfin_sing" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_amv ( min_sta , min_end , ilcdf , nsgnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_sing" organizes the reading,
!   pre-selection, distribution and storage in ODR of AMV (Atmospheric Motion
!   Vector derived from (sequences of) satellite VIS, IR, or WV radiance data)
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
    nsgnew        ,& ! number of new single-level reports
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
                      !  = ncdf_amdar     : AMDAR single-level
                      !  = ncdf_acars     : ACARS single-level
                      !  = ncdf_acars_uk  : ACARS UK + Canada
                      !  = ncdf_acars_us  : ACARS USA (with humidity)
                      !  = ncdf_amdar_ml  : AMDAR multi-level
                      !  = ncdf_amdar_vp  : AMDAR vertical profile
                      !  = ncdf_modes     : MODE-S single-level KNMI format
                      !  = ncdf_modes_acr : MODE-S single-level ACARS format
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
    ndimq          ,& ! number of dimensions (for 'edition_number' in NetCDF)
    status         ,& ! NetCDF status variable
    istat  , ierr     ! error indicators

  INTEGER        ::  &          ! variable ID's in NetCDF file for:
!                              !   all file types (synop, ship, buoy, metar):
    varid_rtyp     ,& ! satellite derived wind comput. method [002023]
    varid_qual     ,& ! per cent confidence                   [033007]
    varid_sccf     ,& ! satellite channel centre frequency    [1/s]
    varid_dd       ,& ! wind direction
    varid_ff       ,& ! wind speed
    varid_p           ! pressure

  CHARACTER (LEN=25)       :: &
    yroutine          ! name of this subroutine
! CHARACTER (LEN=30)       :: &
!   yerr              ! error message
  CHARACTER (LEN=70)       :: &
    yerrmsl           ! error message
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
! Begin Subroutine obs_cdf_read_amv
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_amv'

  ncid     =  ncinid (ilcdf)
  kcdftyp  =  icdfin (ilcdf)
  yfn      =  ycdfin (kcdftyp) (1:LEN_TRIM( ycdfin (kcdftyp) )) //             &
              yncannex (ilcdf) (1:LEN_TRIM( yncannex (ilcdf) ))

! lamdarml =      (kcdftyp == ncdf_amdar_ml) .OR. (kcdftyp == ncdf_amdar_vp)
! lamdar   =      (kcdftyp == ncdf_amdar   ) .OR. (lamdarml)
! lacars   =      (kcdftyp == ncdf_acars   ) .OR. (kcdftyp == ncdf_acars_uk)   &
!                                            .OR. (kcdftyp == ncdf_acars_us)   &
!                                            .OR. (kcdftyp == ncdf_modes_acr)

! PRINT *,yroutine, kcdftyp, min_sta,min_end, lamdarml, lamdar, lacars

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
  noffscal (2,2)  =  noffscal(2,1) + 3
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
! Section 2: Get entries specific to observation type (AMDAR, ACARS, etc.)
!-------------------------------------------------------------------------------

! get dimension of element MPCCO
    ndimq = 2
                          status = nf90_inq_varid (ncid,'MPCCO',varid_qual )
    IF (status == nf90_noerr) status = nf90_Inquire_Variable ( ncid, varid_qual &
                                                      , ndims=ndimq, dimids=dimids)
!   status = nf90_Inquire_Dimension (ncid, dimids(1), len=nqual)

! allocate arrays required for reading
! ------------------------------------

    ! initialise mandatory elements
    ALLOCATE ( nc_dd    (1,nrep) , STAT=istat )
    ALLOCATE ( rc_ff    (1,nrep) , STAT=istat )
    ALLOCATE ( rc_p     (1,nrep) , STAT=istat )
    ALLOCATE ( nc_rtyp  (  nrep) , STAT=istat )
    ALLOCATE ( rc_sccf  (  nrep) , STAT=istat )
    ALLOCATE ( nc_qual  (  nrep) , STAT=istat )
    nc_qual  (:) = imiss

! get mandatory (!) variable ID's in NetCDF file and get data
! -----------------------------------------------------------
                              status = nf90_inq_varid (ncid,'NDNDN' ,varid_dd  )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'NFNFN' ,varid_ff  )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MPN  ' ,varid_p   )
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MSDWCM',varid_rtyp)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MSCCF' ,varid_sccf)
    IF (status == nf90_noerr) status = nf90_inq_varid (ncid,'MPCCO' ,varid_qual)
    IF (status /= nf90_noerr) THEN
      yerrmsl = 'NDNDN, NFNFN, MPN, MSDWCM, MSCCF, or MPCCO do not exist in '  &
                // yfn
      CALL model_abort (my_cart_id, 11001, yerrmsl, yroutine)
    ENDIF
    status = nf90_get_var (ncid, varid_dd  , nc_dd  , (/  mrepsta/), (/  nrep/))
    status = nf90_get_var (ncid, varid_ff  , rc_ff  , (/  mrepsta/), (/  nrep/))
    status = nf90_get_var (ncid, varid_p   , rc_p   , (/  mrepsta/), (/  nrep/))
    status = nf90_get_var (ncid, varid_rtyp, nc_rtyp, (/  mrepsta/), (/  nrep/))
    status = nf90_get_var (ncid, varid_sccf, rc_sccf, (/  mrepsta/), (/  nrep/))
    IF (ndimq == 1) THEN
      status= nf90_get_var(ncid, varid_qual, nc_qual, (/  mrepsta/), (/  nrep/))
    ELSEIF (ndimq == 2) THEN
      !   from the 2-dim array 'MOGC0', always take the 1st element for each obs
      status= nf90_get_var(ncid, varid_qual, nc_qual, (/1,mrepsta/), (/1,nrep/))
!     status= nf90_get_var(ncid, varid_qual, nc_qual, start=(/1,mrepsta/)      &
!                                                   , count=(/1,nrep/))
    ENDIF

!-------------------------------------------------------------------------------
! Section 3: Produce long buffer arrays 'Xbufsrt' in which the read reports are
!            sorted according to the nodes to which they will be distributed
!-------------------------------------------------------------------------------

    ALLOCATE ( iirlen (nreproc+1) , STAT=istat )
    ALLOCATE ( irrlen (nreproc+1) , STAT=istat )
    ALLOCATE ( iyrlen (nreproc+1) , STAT=istat )

! compute length (number of elements) for each individual report
! --------------------------------------------------------------

    DO irps = 1 , nreproc
      iirlen (irps)  =  noffscal(1,2)
      irrlen (irps)  =  noffscal(2,2)
      iyrlen (irps)  =  noffscal(3,2)
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
    ibufsrt (:) = imiss
    rbufsrt (:) = rmiss

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
      ibufsrt (iioffs+noffscal(1,1)+ 1)  =  nc_rtyp (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 2)  =  nc_qual (irep)
      ibufsrt (iioffs+noffscal(1,1)+ 3)  =  nc_dd (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 1)  =  rc_ff (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 2)  =  rc_p  (1,irep)
      rbufsrt (iroffs+noffscal(2,1)+ 3)  =  rc_sccf (irep)

      iioffs  =  iioffs + iirlen(irps)
      iroffs  =  iroffs + irrlen(irps)
    ENDDO

! de-allocate reading arrays
    DEALLOCATE ( nc_dd     , STAT=istat )
    DEALLOCATE ( rc_ff     , STAT=istat )
    DEALLOCATE ( rc_p      , STAT=istat )
    DEALLOCATE ( nc_rtyp   , STAT=istat )
    DEALLOCATE ( rc_sccf   , STAT=istat )
    DEALLOCATE ( nc_qual   , STAT=istat )

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

    CALL obs_cdf_store_singlev ( kcdftyp , nrepl , nlenli , nlenlr , nlenly    &
                               , noffscal,         ibufloc, rbufloc, ybufloc   &
                               , nsgnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufloc , STAT=istat )
    DEALLOCATE ( rbufloc , STAT=istat )
    DEALLOCATE ( ybufloc , STAT=istat )

  ELSEIF (nrep >= 1) THEN
    nrepl   =  nodrepn(1)
    nlenli  =  nodleni(1)
    nlenlr  =  nodlenr(1)
    nlenly  =  nodleny(1)

    CALL obs_cdf_store_singlev ( kcdftyp , nrepl , nlenli , nlenlr , nlenly    &
                               , noffscal,         ibufsrt, rbufsrt, ybufsrt   &
                               , nsgnew  , nexceed )
!   ==========================

    DEALLOCATE ( ibufsrt , STAT=istat )
    DEALLOCATE ( rbufsrt , STAT=istat )
    DEALLOCATE ( ybufsrt , STAT=istat )

  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_amv
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_amv


!===============================================================================
!+ Module procedure in "src_obs_cdfin_sing" for storing single-level rep. in ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_store_singlev  ( kcdftyp, nrepl, nlenli, nlenlr, nlenly     &
                                  , noffscal      , ibuf  , rbuf  , ybuf       &
                                  , nsgnew , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_sing" stores the single-level
!   reports (except for GPS reports) in the internal ODR arrays.
!
! Method:
!   First, the header part common to all observation types is stored by calling
!   'obs_cdf_store_comhead', and then the rest of the header is added which
!   depends on the observation type.
!   Next, the single-level observations are copied from the observation type
!   dependent buffers into standardised arrays. Subsequently, the quality flag
!   patterns are built (e.g. from blacklisting or gross error checking), and
!   the elements are stored in the single-level ODR (internal arrays). Pressure,
!   if missing, is derived from height for upper-air reports. Low, mid-level,
!   and high cloud cover is derived from reported information on cloud, weather,
!   and visibility.
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
                        !  = ncdf_amdar     : AMDAR single-level
                        !  = ncdf_acars     : ACARS single-level
                        !  = ncdf_acars_uk  : ACARS UK + Canada
                        !  = ncdf_acars_us  : ACARS USA with humidity
                        !  = ncdf_modes     : MODE-S single-level KNMI format
                        !  = ncdf_modes_acr : MODE-S single-level ACARS format
                        !  = ncdf_synop     : SYNOP
                        !  = ncdf_synop_mob : SYNOP mobile
                        !  = ncdf_ship      : SHIP
                        !  = ncdf_buoy      : BUOY
                        !  = ncdf_metar     : METAR (surface aviation)
                        !  = ncdf_satob     : SATOB / AMV
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
    nsgnew           ,& ! number of new single-level reports
    nexceed             ! number of reports in excess of array size

! Local parameters:
! ----------------
! INTEGER        , PARAMETER ::  &
!   h_low = 1999     ,& ! height limit below which cloud is low cloud       [m]
!   h_mid = 6999        ! height limit below which cloud is mid-level cloud [m]

  REAL (KIND=wp) , PARAMETER ::  &
    c100r = 0.01_wp     ! maximum number of elements to be received

! Local scalars:
! -------------

  LOGICAL        , SAVE      ::  &
    lfirst = .TRUE.     ! routine called the first time (for extra crtl. output)

  LOGICAL        ::  &
    lblkw  (nrepl+1) ,& ! report blacklisted because sta-ID missing on whitelist
    lsurfob          ,& ! surface-level report (SYNOP, SHIP, BUOY, or METAR)
    lscatto          ,& ! scatterometer report (ASCAT, QuickScat)
    lacars           ,& ! ACARS
    lupair           ,& ! upper-air     report (AMDAR, ACARS, or SATOB)
    lrhtdqx          ,& ! report template contains: RH + dewpoint + mixing ratio
    lrhtd            ,& ! report template contains: rel. hum. + dewpoint temp.
    ltd              ,& ! report template contains: dewpoint temperature only
    lexitd           ,& ! report template includes dewpoint temperature
    lexirh           ,& ! report template includes relative humidity
    lexiqx           ,& ! report template includes mixing ratio
    lseaobs          ,& ! sea report (ship, buoy, not land report)
    lseaonly         ,& ! use actively only observations from sea platforms
    loro             ,& ! set passive: upper-air report too close to orography,
                        !   or 10-m wind in too large sub-grid scale orography
    lexcil           ,& ! ceiling is based on 'MNH'-based cloud cover >= 5 octa
    lrej                ! reject upper-air report (no valid pressure level)

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
    nsgob            ,& ! target ODR single-level report index
    nicomh , nrcomh  ,& ! number of common header integer / real  elements
    niscal , nrscal  ,& ! number of scalar integer / real  elements
    nrepsg           ,& ! number of reports which can be stored in the ODR now
    nexcesg          ,& ! no. of local sing-lev. reports in excess of array size
    irpl             ,& ! loop index over local reports
    ilev             ,& ! loop index over cloud layers
    intv             ,& ! loop index over blacklisted vertical intervals
    icl              ,& ! loop index over single characters
    maxlev           ,& ! max. number of cloud layers
    nlev   , nlv     ,& ! number of vertical levels in current report
    kobtyp , kcdtyp  ,& ! observation type / observation code type
    iofstyp          ,& ! index offset for station type element
    istalt           ,& ! station altitude
    npfi             ,& ! pressure level [hPa] of reported height
    mpassiv          ,& ! active / passive status of report
    ilverr           ,& ! nearest standard error level below observation level
    ixrhs               ! index for 'rhtsat'

  INTEGER        ::  &
    mphase           ,& ! phase of flight                     (WMO Table 008004)
    nrtr             ,& ! time period                  (VUB WMO Code Table 4019)
    nww    , nw2     ,& ! present or past weather / secondary past weather
    nclct            ,& ! total cloud cover            (octas, WMO Table 020011)
    nhkey            ,& ! class of cloud base height   (VUB WMO Code Table 1600)
    ncsig            ,& ! vertical significance rel. to cloud (WMO Table 008002)
!   nclcl            ,& ! low    cloud cover           (octas, WMO Table 020011)
!   nclcm            ,& ! middle cloud cover           (octas, WMO Table 020011)
!   nclch            ,& ! high   cloud cover           (octas, WMO Table 020011)
    nclqf            ,& ! quality flag for derived low cloud cover
    ncc              ,& ! (low, middle, or high) cloud type
    nicl             ,& ! individual cloud layer group (ICLG) (code)
    iclg_cb          ,& ! cloud base height      from ICLG [m]
    iclg_ci          ,& ! ceiling                from ICLG [m]
    ivis             ,& ! visibility
    iob    , job     ,& ! indices of grid point assigned to observation
    iactr            ,& ! stepping for event counters
    ilen             ,& ! length of control message
    iactx (5)        ,& ! indicates presence of active  data (for variables)
    ipasx (5)        ,& ! indicates presence of passive data (for variables)
    nzaexi           ,& ! -1: only passive data ; 0: no data ; 1: active data
    ilstid_chk       ,& ! length of (part of) station ID that can be checked
    istat  , irm        ! error indicators

  REAL (KIND=wp) ::  &
    zobhr (nrepl+1)  ,& ! obs  time (input array for local whitelist)
    zpref (nrepl+1)  ,& ! 'reference pressure' for evaluating the blacklist
                        !   (this must be independent from model state)
    zpuse (nrepl+1)  ,& ! reported pressure
    zzuse (nrepl+1)  ,& ! reported height
    zoberr(nrepl+1)  ,& ! observation error
    fisd  (nrepl+1)  ,& ! height diff. betw. model orography and sta. altitude
    zuu   (nrepl+1,1),& ! zonal      wind                (processed levels only)
    zvv   (nrepl+1,1),& ! meridional wind                (processed levels only)
    zrlat (nrepl+1,1),& ! latitude  of (processed) observation level
    zrlon (nrepl+1,1),& ! longitude of (processed) observation level
    ztmean           ,& ! mean temperature (betw. station height and msl)
    zpred            ,& ! reported station pressure reduced to mean sea level
    zplim            ,& ! threshold for reported station pressure check
    zporo            ,& ! (approx.) model surface pressure
    zadder           ,& ! additional pressure obs error due to extrapolation
    fisduv           ,& ! scaled extrapolation distance for 10-m wind values
    fisdtt           ,& ! scaled extrapolation distance for 2-m temperature val
    fisdrh           ,& ! scaled extrapolation distance for 2-m humidity values
    fisdzz           ,& ! scaled extrapolation distance for surface press. val.
    fisdnu           ,& ! scaled extrapolation distance for full values and
                        !   increments of surface pressure (for nudging only)
    zlaplac          ,& ! (proportional to) Laplacian of orography
    zzkeml           ,& ! height of the lowest main model level
    zz2p             ,& ! (approx.) conversion from height to pressure units
    zrr1   , zrr3    ,& ! precipitation sum over 1 hour resp. 3  hours
    zrr6             ,& ! precipitation sum over 6 hours
    zr12   , zr24    ,& ! precipitation sum over 12     resp. 24 hours
    clcd (4)         ,& ! low / mid / high / total cloud amount from GCG [octas]
    clg_cl           ,& ! low       cloud amount from ICLG [octas]
    clg_cm           ,& ! mid-level cloud amount from ICLG [octas]
    clg_ch           ,& ! high      cloud amount from ICLG [octas]
    clg_cml          ,& ! cloud amount of lowest mid-level layer in ICLG [octas]
    zlop             ,& ! LOG( pressure )
    zstalt           ,& ! station altitude
    fiperr              ! interpolation weight factor for error level 'ilverr'

! CHARACTER (LEN=25)       :: &
!   yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid_blk)  :: &
    ystidl  (nrepl+1)   ! observation type (input array for local blacklist)

! Local arrays:
! ------------

  LOGICAL        , ALLOCATABLE ::  &
    lpract       (:) ,& ! bad reporting practice for surface pressure
                        !                   (and for aircraft height)
    landsy       (:)    ! land synoptic report

  INTEGER        , ALLOCATABLE ::  &
    npcode       (:) ,& ! pressure code
    kfog         (:) ,& ! fog (probability) indicator:
                        ! = 3: yes  : reported flag   ; = 2: likely: low visib.
                        ! = 1: doubt: low vis + precip; = 0: no: no indication
                        ! =-1: dust, sand --> no fog even if vis is low
    icbh         (:) ,& ! cloud base height from general cloud group
    iroll        (:) ,& ! aircraft roll angle quality         (WMO Table 002064)
    iturb        (:) ,& ! aircraft degree of turbulence       (WMO Table 011031)
    iqual        (:) ,& ! per cent confidence                 (%)
    iqxq         (:) ,& ! mixing ratio quality                (WMO Table 033026)
    ivdt         (:) ,& ! time period of wind measurement
    ivtisi       (:) ,& ! time significance of wind obs       (WMO Table 008021)
    igudt        (:) ,& ! first  time period of max. wind gust speed
    igudt2       (:) ,& ! second time period of max. wind gust speed
    irrdt        (:) ,& ! first  time period of precipitation sum
    irrdt2       (:) ,& ! second time period of precipitation sum
    iclev        (:) ,& ! number of reported individual cloud layers
    iclct        (:) ,& ! total cloud amount
    iclsig       (:) ,& ! vertical significance rel. to cloud (WMO Table 008002)
    iclclm       (:) ,& ! lowest (low or middle) cloud amount (WMO Table 020011)
    iccl         (:) ,& ! low       cloud type                (WMO Table 020012)
    iccm         (:) ,& ! mid-level cloud type                (WMO Table 020012)
    icch         (:) ,& ! high      cloud type                (WMO Table 020012)
    iww          (:) ,& ! present weather                     (WMO Table 020003)
    iw1          (:) ,& ! past    weather, first  value       (WMO Table 020004)
    iw2          (:) ,& ! past    weather, second value       (WMO Table 020004)
    iwdt         (:) ,& ! time period related to past weather report
    itmxt        (:) ,& ! (negative) time period of maximum 2-m temperature
    itmnt        (:) ,& ! (negative) time period of minimum 2-m temperature
    iee          (:) ,& ! state of ground                     (WMO Table 020062)
    ifi          (:) ,& ! geopotential height (at specified pressure level)
    icbase       (:) ,& ! cloud base height of indiv. cloud layer
    iclxsg     (:,:) ,& ! vert. signif. of indiv. cloud layer (WMO Table 008002)
    iclxcl     (:,:) ,& ! cloud amount  of indiv. cloud layer (WMO Table 020011)
    iclxct     (:,:)    ! cloud type    of indiv. cloud layer (WMO Table 020012)

  REAL (KIND=wp) , ALLOCATABLE ::  &
    zrhc         (:) ,& ! relative humidity (model compatible, bias corrected)
    zrhw         (:) ,& ! relative humidity over water
    ztd          (:) ,& ! dew point temperature
    zqx          (:) ,& ! mixing ratio
    zpp          (:) ,& ! pressure
    zdd          (:) ,& ! wind direction
    zff          (:) ,& ! wind speed
    ztt          (:) ,& ! temperature
    zzz          (:) ,& ! height
    zttq         (:) ,& ! precision of temperature observation
    zvgust       (:) ,& ! maximum derived equivalent vertical gust
    zgust        (:) ,& ! maximum 10-m wind gust speed
    zgust2       (:) ,& ! maximum 10-m wind gust speed, second period
    zzt          (:) ,& ! height of temperature sensor above ground
    zpmsl        (:) ,& ! mean sea level pressure
    zdpdt        (:) ,& ! 3-hour pressure change
    zrr          (:) ,& ! precipitation sum
    zrr2         (:) ,& ! precipitation sum, second period
    zrr24        (:) ,& ! 24-hour precipitation sum
    zvis         (:) ,& ! (horizontal) visibility
    zcbase       (:) ,& ! height of base of lowest cloud
    ztmax        (:) ,& ! maximum 2-m temperature over time period 'itmxt'
    ztmin        (:) ,& ! minimum 2-m temperature over time period 'itmnt'
    zradgl       (:) ,& ! global  solar radiation, sum over 1 hour        [J/m2]
    zraddf       (:) ,& ! diffuse solar radiation, sum over 1 hour        [J/m2]
    zradlw       (:) ,& ! long-wave     radiation, sum over 1 hour        [J/m2]
    zhsnow       (:) ,& ! total snow depth
    zpfi         (:) ,& ! pressure level of reported height
    zclxbs     (:,:) ,& ! height of base of individual cloud layer
    zclxc        (:) ,& ! re-assigned cloud cover from indiv. cloud layer
    zrhw1        (:) ,& ! relative humidity over water (converted from dewpoint)
    zrhw2        (:) ,& ! relative humidity over water (conv. from mixing ratio)
    zrh1         (:) ,& ! model compatible rel humid.  (converted from dewpoint)
    zrh2         (:) ,& ! model compatible rel humid.  (conv. from mixing ratio)
    zqvw         (:) ,& ! specific humidity over water (conv. from mixing ratio)
    zqv          (:)    ! model compatible spec. humid (conv. from mixing ratio)
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_store_singlev
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_store_singlev'

  lsurfob  =      (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )         &
             .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_metar)         &
             .OR. (kcdftyp == ncdf_nat_dwd)

  lscatto  =      (kcdftyp == ncdf_ascat) .OR. (kcdftyp == ncdf_qscat)

  lacars   =      (kcdftyp == ncdf_acars) .OR. (kcdftyp == ncdf_acars_uk)      &
                                          .OR. (kcdftyp == ncdf_acars_us)      &
                                          .OR. (kcdftyp == ncdf_modes_acr)

  lupair   =      (kcdftyp == ncdf_amdar) .OR. (lacars)                        &
             .OR. (kcdftyp == ncdf_satob) .OR. (kcdftyp == ncdf_modes)
  ixrhs    = 2
  IF (lupair) ixrhs = 5

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

  CALL obs_cdf_store_comhead ( 2, ntotsg, nrepl, nlenli, nlenlr, nlenly        &
                             , ibuf, rbuf, ybuf, nsgnew, nexcesg )
! ==========================

! update number of reports which can be stored in the ODR now
! (the following line is equivalent to:  nrepsg = nsgnew)
  nrepsg = nrepl - nexcesg

! update counter for insufficient ODR size, for caution messages
  nexceed = nexceed + nexcesg

! store model background info which is available (only) on local sub-domain
! -------------------------------------------------------------------------
! - standard deviation of sub-grid scale model orography (for 10-m wind)
  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    osghed (nsgob,nhssos) = r_sso_sd(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))
  ENDDO

! store header part, which is specific to NetCDF
! observation input file type, in single-level ODR:
! ------------------------------------------------
! 'instrument specification': - NIX    (station type) for synop, ship, metar
!                             - MTODB  (type of data buoy,        Table 002149)
!                           ! - MQOBL  (quality of buoy location, Table 033023)
!                             - MSDWCM (wind comput. method) for satob (AMV)
                              iofstyp = 1
  IF (kcdftyp == ncdf_buoy )  iofstyp = 8
  IF (lacars)                 iofstyp = 8
  IF ((lsurfob) .OR. (lacars)) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      mosghd (nsgob,nhstyp) = imdi
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ iofstyp) /= imiss)                  &
        mosghd (nsgob,nhstyp) = ibuf (iioffs(irpl)+noffscal(1,1)+ iofstyp)
    ENDDO
  ENDIF
  IF (kcdftyp == ncdf_satob) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      mosghd (nsgob,nhstyp) = imdi
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ iofstyp) /= imiss) THEN
        mosghd (nsgob,nhstyp) = ibuf (iioffs(irpl)+noffscal(1,1)+ iofstyp)
        CALL check_wind_method ( mosghd(nsgob,nhstid)                          &
                               , rbuf (iroffs(irpl)+noffscal(2,1)+ 3)          &
                               , mosghd(nsgob,nhstyp) )
      ! ======================
      ENDIF
    ENDDO
  ENDIF
! aircraft: phase of flight
  IF ((kcdftyp == ncdf_amdar) .OR. (lacars) .OR. (kcdftyp == ncdf_modes)) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      mphase = nibits(nvapoc)
      IF (ibuf(iioffs(irpl)+noffscal(1,1)+ 1) /= imiss)                        &
        mphase                = ibuf (iioffs(irpl)+noffscal(1,1)+ 1)
      CALL MVBITS( mphase, 0, nvapoc, mosghd(nsgob,nhschr), nvapbp )
!     mosghd (nsgob,nhschr) = ireplace( mosghd(nsgob,nhschr)                   &
!                                     , nvapbp, nvapoc, mphase , 0 )
    ENDDO
  ENDIF

! get diagnostic array position (--> i_cma)
! -----------------------------------------

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    kobtyp = mosghd(nsgob,nhobtp)
    kcdtyp = mosghd(nsgob,nhcode)
    icmaa (irpl)  =  i_cma ( kobtyp , kcdtyp )
!                    =====
  ENDDO

! get list of blacklisted vertical intervals for each report
! ----------------------------------------------------------

  ilstid_chk  =  MIN( ilstid, ilstid_blk )
  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    jobtyp (irpl)  =  mosghd(nsgob,nhobtp)
    jcdtyp (irpl)  =  mosghd(nsgob,nhcode)
    zobhr  (irpl)  =  osghed(nsgob,nhtime)
    ystidl (irpl)  =  ' '
    ystidl (irpl)  =  yosghd(nsgob) (1:ilstid_chk)
  ENDDO

  IF (nrepsg >= 1) THEN
    ALLOCATE ( blk_loc (nrepsg) , STAT=istat )

    CALL obs_cdf_blacklist_local ( nrepsg , jobtyp , ilstid_chk , ystidl )
!   ============================

! determine which reports are missing on the whitelists
! -----------------------------------------------------

    CALL obs_cdf_whitelist_local ( nrepsg , jobtyp , jcdtyp , zobhr            &
                                 , ilstid_chk , ystidl , lblkw )
!   ============================

  ENDIF

!-------------------------------------------------------------------------------
! Section 2: Report body: Put the observations from the observation type
!            dependent buffers into standardised arrays
!-------------------------------------------------------------------------------

  maxlev = 0

  lrhtd   =  (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )        &
              .OR. (kcdftyp == ncdf_buoy ))    ! .AND. (kcdftyp /= ncdf_nat_dwd)
  ltd     =        (kcdftyp == ncdf_metar)
  lrhtdqx =  (     (kcdftyp == ncdf_amdar) .OR. (lacars)                       &
              .OR. (kcdftyp == ncdf_modes))
  lexitd  =  (lrhtdqx) .OR. (lrhtd) .OR. (ltd)
  lexirh  =  (lrhtdqx) .OR. (lrhtd)
  lexiqx  =  (lrhtdqx)

! allocate standardised arrays
! ----------------------------

  IF (nrepsg >= 1) THEN
    ALLOCATE ( lpract  (nrepsg) , STAT=istat )

    IF (lsurfob) THEN
      ALLOCATE ( npcode  (nrepsg) , STAT=istat )
      ALLOCATE ( landsy  (nrepsg) , STAT=istat )
      ALLOCATE ( kfog    (nrepsg) , STAT=istat )
      ALLOCATE ( icbh    (nrepsg) , STAT=istat )
    ENDIF

    IF ((lexitd) .OR. (lexirh) .OR. (lexiqx)) THEN
      ALLOCATE ( zrhw    (nrepsg) , STAT=istat )
      ALLOCATE ( ztd     (nrepsg) , STAT=istat )
      ALLOCATE ( zqx     (nrepsg) , STAT=istat )
      ALLOCATE ( zrhc    (nrepsg) , STAT=istat )
    ENDIF

! 'zpp' allocated for all types, although not available for AMDAR, METAR, SCATT
    ALLOCATE ( zpp     (nrepsg) , STAT=istat )
    ALLOCATE ( zdd     (nrepsg) , STAT=istat )
    ALLOCATE ( zff     (nrepsg) , STAT=istat )
    IF ((kcdftyp /= ncdf_satob) .AND. (.NOT. lscatto))                         &
      ALLOCATE ( ztt     (nrepsg) , STAT=istat )

    IF (lupair) THEN
      ALLOCATE ( zzz     (nrepsg) , STAT=istat )
      IF (kcdftyp /= ncdf_satob)  ALLOCATE ( iroll   (nrepsg) , STAT=istat )
      IF (kcdftyp == ncdf_satob)  ALLOCATE ( iqual   (nrepsg) , STAT=istat )
      IF (lacars)                 ALLOCATE ( zttq    (nrepsg) , STAT=istat )
      IF ((kcdftyp == ncdf_modes) .OR. (lacars))                               &
        ALLOCATE ( iqxq    (nrepsg) , STAT=istat )
      IF ((kcdftyp == ncdf_amdar) .OR. (lacars)) THEN
        ALLOCATE ( iturb   (nrepsg) , STAT=istat )
        ALLOCATE ( zvgust  (nrepsg) , STAT=istat )
      ENDIF
    ENDIF

    IF (lsurfob) THEN
      ALLOCATE ( zpmsl   (nrepsg) , STAT=istat )
!     ALLOCATE ( ztd     (nrepsg) , STAT=istat )
      ALLOCATE ( zgust   (nrepsg) , STAT=istat )
      ALLOCATE ( zzt     (nrepsg) , STAT=istat )
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        ALLOCATE ( ivdt    (nrepsg) , STAT=istat )
        ALLOCATE ( igudt   (nrepsg) , STAT=istat )
        ALLOCATE ( ivtisi  (nrepsg) , STAT=istat )
        ALLOCATE ( irrdt   (nrepsg) , STAT=istat )
        ALLOCATE ( zrr     (nrepsg) , STAT=istat )
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy )) THEN
        ALLOCATE ( zdpdt   (nrepsg) , STAT=istat )
!       ALLOCATE ( zrhw    (nrepsg) , STAT=istat )
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        ALLOCATE ( igudt2  (nrepsg) , STAT=istat )
        ALLOCATE ( irrdt2  (nrepsg) , STAT=istat )
        ALLOCATE ( zgust2  (nrepsg) , STAT=istat )
        ALLOCATE ( zrr2    (nrepsg) , STAT=istat )
      ENDIF
      IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
        ALLOCATE ( iclev   (nrepsg) , STAT=istat )
        ALLOCATE ( iclct   (nrepsg) , STAT=istat )
        ALLOCATE ( iclsig  (nrepsg) , STAT=istat )
        ALLOCATE ( iclclm  (nrepsg) , STAT=istat )
        ALLOCATE ( iccl    (nrepsg) , STAT=istat )
        ALLOCATE ( iccm    (nrepsg) , STAT=istat )
        ALLOCATE ( icch    (nrepsg) , STAT=istat )
        ALLOCATE ( iww     (nrepsg) , STAT=istat )
        ALLOCATE ( iwdt    (nrepsg) , STAT=istat )
        ALLOCATE ( iw1     (nrepsg) , STAT=istat )
        ALLOCATE ( iw2     (nrepsg) , STAT=istat )
        ALLOCATE ( itmxt   (nrepsg) , STAT=istat )
        ALLOCATE ( itmnt   (nrepsg) , STAT=istat )
        ALLOCATE ( zvis    (nrepsg) , STAT=istat )
        ALLOCATE ( zrr24   (nrepsg) , STAT=istat )
        ALLOCATE ( zcbase  (nrepsg) , STAT=istat )
        ALLOCATE ( ztmax   (nrepsg) , STAT=istat )
        ALLOCATE ( ztmin   (nrepsg) , STAT=istat )
      ENDIF
      IF     (kcdftyp == ncdf_synop) THEN
        ALLOCATE ( iee     (nrepsg) , STAT=istat )
        ALLOCATE ( ifi     (nrepsg) , STAT=istat )
        ALLOCATE ( zradgl  (nrepsg) , STAT=istat )
        ALLOCATE ( zraddf  (nrepsg) , STAT=istat )
        ALLOCATE ( zradlw  (nrepsg) , STAT=istat )
        ALLOCATE ( zhsnow  (nrepsg) , STAT=istat )
        ALLOCATE ( zpfi    (nrepsg) , STAT=istat )
      ELSEIF (kcdftyp == ncdf_metar) THEN
        ALLOCATE ( iclev   (nrepsg) , STAT=istat )
        ALLOCATE ( zvis    (nrepsg) , STAT=istat )
      ENDIF
    ENDIF
  ENDIF

! initialise standardised arrays
! ------------------------------

! PRINT *,'zrhrmdi ', rmdi, imdi, lexitd, lexirh, lexiqx, nrepsg
  DO irpl = 1 , nrepsg
    IF ((lexitd) .OR. (lexirh) .OR. (lexiqx)) THEN
      zrhw   (irpl) = rmdi
      ztd    (irpl) = rmdi
      zqx    (irpl) = rmdi
      zrhc   (irpl) = rmdi
    ENDIF

    zpp    (irpl) = rmdi
    zdd    (irpl) = rmdi
    zff    (irpl) = rmdi
    IF ((kcdftyp /= ncdf_satob) .AND. (.NOT. lscatto))                         &
      ztt    (irpl) = rmdi

    IF (lupair) THEN
      zzz    (irpl) = rmdi
      IF  (kcdftyp /= ncdf_satob)                 iroll  (irpl) = imdi
      IF  (kcdftyp == ncdf_satob)                 iqual  (irpl) = imdi
      IF  (lacars)                                zttq   (irpl) = rmdi
      IF ((kcdftyp == ncdf_modes) .OR. (lacars))  iqxq   (irpl) = imdi
      IF ((kcdftyp == ncdf_amdar) .OR. (lacars))  iturb  (irpl) = imdi
      IF ((kcdftyp == ncdf_amdar) .OR. (lacars))  zvgust (irpl) = rmdi
    ENDIF

    IF (lsurfob) THEN
      zpmsl  (irpl) = rmdi
!     ztd    (irpl) = rmdi
      zgust  (irpl) = rmdi
      zzt    (irpl) = rmdi
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        ivdt   (irpl) = imdi
        igudt  (irpl) = imdi
        ivtisi (irpl) = imdi
        irrdt  (irpl) = imdi
        zrr    (irpl) = rmdi
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy )) THEN
        zdpdt  (irpl) = rmdi
!       zrhw   (irpl) = rmdi
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        igudt2 (irpl) = imdi
        irrdt2 (irpl) = imdi
        zgust2 (irpl) = rmdi
        zrr2   (irpl) = rmdi
      ENDIF
      IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
        iclev  (irpl) = 0
        iclct  (irpl) = imdi
        iclsig (irpl) = imdi
        iclclm (irpl) = imdi
        iccl   (irpl) = imdi
        iccm   (irpl) = imdi
        icch   (irpl) = imdi
        iww    (irpl) = imdi
        iwdt   (irpl) = imdi
        iw1    (irpl) = imdi
        iw2    (irpl) = imdi
        itmxt  (irpl) = 0
        itmnt  (irpl) = 0
        zvis   (irpl) = rmdi
        zrr24  (irpl) = rmdi
        zcbase (irpl) = rmdi
        ztmax  (irpl) = rmdi
        ztmin  (irpl) = rmdi
      ENDIF
      IF     (kcdftyp == ncdf_synop) THEN
        iee    (irpl) = imdi
        ifi    (irpl) = imdi
        zradgl (irpl) = rmdi
        zraddf (irpl) = rmdi
        zradlw (irpl) = rmdi
        zhsnow (irpl) = rmdi
        zpfi   (irpl) = rmdi
      ELSEIF (kcdftyp == ncdf_metar) THEN
        iclev  (irpl) = 0
        zvis   (irpl) = rmdi
      ENDIF
    ENDIF
  ENDDO

! fill standardised arrays
! ------------------------

  IF (((lacars) .OR. (kcdftyp == ncdf_amdar)                                   &
                .OR. (kcdftyp == ncdf_modes)) .AND. (nrepsg >= 1)) THEN
    DO irpl = 1 , nrepsg
      IF (     ibuf(iioffs(irpl)+nicomh+ 2)  /= imiss   )                      &
        iroll (irpl) =       ibuf (iioffs(irpl)+nicomh+ 2)
      IF (     ibuf(iioffs(irpl)+nicomh+ 3)  /= imiss   ) THEN
        zdd   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 3) , wp )
      ENDIF
      IF (     ibuf(iioffs(irpl)+nicomh+ 4)  /= imiss   ) THEN
        zzz   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 4) , wp )
      ENDIF
!     PRINT *, 'zzz7  ',irpl, ibuf(iioffs(irpl)+nicomh+ 4), imiss, zzz(irpl)
!     PRINT *, 'zrhm1 ',irpl, ibuf(iioffs(irpl)+nicomh+ 5), imiss, zrhw(irpl)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        zff   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        ztt   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 3) ) < rmisschk)                      &
        ztd   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 3)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 4) ) < rmisschk)                      &
        zqx   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 4)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 5) ) < rmisschk)                      &
        zpp   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 5)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 6) ) < rmisschk)                      &
        zrhw  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 6) * 0.01_wp
      IF ((lacars) .OR. (kcdftyp == ncdf_amdar)) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 5)  /= imiss   )                    &
          iturb (irpl) =       ibuf (iioffs(irpl)+nicomh+ 5)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 7) ) < rmisschk)                    &
          zvgust(irpl) =       rbuf (iroffs(irpl)+nrcomh+ 7)
      ENDIF
      IF (lacars) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 6)  /= imiss   )                    &
          iqxq  (irpl) =       ibuf (iioffs(irpl)+nicomh+ 6)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 8) ) < rmisschk)                    &
          zttq  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 8)
      ENDIF
      IF (kcdftyp == ncdf_modes) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 5)  /= imiss   )                    &
          iqxq  (irpl) =       ibuf (iioffs(irpl)+nicomh+ 5)
      ENDIF
    ENDDO
  ENDIF

  IF ((kcdftyp == ncdf_satob) .AND. (nrepsg >= 1)) THEN
    DO irpl = 1 , nrepsg
      IF (     ibuf(iioffs(irpl)+nicomh+ 2)  /= imiss   )                      &
        iqual (irpl) =       ibuf (iioffs(irpl)+nicomh+ 2)
      IF (     ibuf(iioffs(irpl)+nicomh+ 3)  /= imiss   ) THEN
        zdd   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 3) , wp )
      ENDIF
!     IF (     ibuf(iioffs(irpl)+nicomh+ 4)  /= imiss   ) THEN
!       zzz   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 4) , wp )
!     ENDIF
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        zff   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        zpp   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
    ENDDO
  ENDIF

  IF ((lscatto) .AND. (nrepsg >= 1)) THEN
    DO irpl = 1 , nrepsg
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        zdd   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        zff   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
    ENDDO
  ENDIF

  IF ((lsurfob) .AND. (nrepsg >= 1) .AND. (kcdftyp /= ncdf_nat_dwd)) THEN
    DO irpl = 1 , nrepsg
      IF (     ibuf(iioffs(irpl)+nicomh+ 2)  /= imiss   ) THEN
        zdd   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 2) , wp )
      ENDIF
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        ztt   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        ztd   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 3) ) < rmisschk)                      &
        zff   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 3)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 4) ) < rmisschk)                      &
        zgust (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 4)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 5) ) < rmisschk)                      &
        zzt   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 5)
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy )) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 3)  /= imiss   ) THEN
          zrhw  (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 3) , wp ) * 0.01_wp
        ENDIF
        IF (     ibuf(iioffs(irpl)+nicomh+ 4)  /= imiss   )                    &
          ivdt  (irpl) =       ibuf (iioffs(irpl)+nicomh+ 4)
        IF (     ibuf(iioffs(irpl)+nicomh+ 5)  /= imiss   )                    &
          igudt (irpl) =       ibuf (iioffs(irpl)+nicomh+ 5)
        IF (     ibuf(iioffs(irpl)+nicomh+ 6)  /= imiss   )                    &
          ivtisi(irpl) =       ibuf (iioffs(irpl)+nicomh+ 6)
        IF (     ibuf(iioffs(irpl)+nicomh+ 7)  /= imiss   )                    &
          irrdt (irpl) =       ibuf (iioffs(irpl)+nicomh+ 7)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 6) ) < rmisschk)                    &
          zpp   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 6)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 7) ) < rmisschk)                    &
          zpmsl (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 7)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 8) ) < rmisschk)                    &
          zdpdt (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 8)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 9) ) < rmisschk)                    &
          zrr   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 9)
      ENDIF
      IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 8)  /= imiss   )                    &
          igudt2(irpl) =       ibuf (iioffs(irpl)+nicomh+ 8)
        IF (     ibuf(iioffs(irpl)+nicomh+ 9)  /= imiss   )                    &
          irrdt2(irpl) =       ibuf (iioffs(irpl)+nicomh+ 9)
        IF (     ibuf(iioffs(irpl)+nicomh+10)  /= imiss   )                    &
          iclev (irpl) =       ibuf (iioffs(irpl)+nicomh+10)
        IF (     ibuf(iioffs(irpl)+nicomh+11)  /= imiss   )                    &
          iclct (irpl) =       ibuf (iioffs(irpl)+nicomh+11)
        IF (     ibuf(iioffs(irpl)+nicomh+12)  /= imiss   )                    &
          iclsig(irpl) =       ibuf (iioffs(irpl)+nicomh+12)
        IF (     ibuf(iioffs(irpl)+nicomh+13)  /= imiss   )                    &
          iclclm(irpl) =       ibuf (iioffs(irpl)+nicomh+13)
        IF (     ibuf(iioffs(irpl)+nicomh+14)  /= imiss   )                    &
          iccl  (irpl) =       ibuf (iioffs(irpl)+nicomh+14)
        IF (     ibuf(iioffs(irpl)+nicomh+15)  /= imiss   )                    &
          iccm  (irpl) =       ibuf (iioffs(irpl)+nicomh+15)
        IF (     ibuf(iioffs(irpl)+nicomh+16)  /= imiss   )                    &
          icch  (irpl) =       ibuf (iioffs(irpl)+nicomh+16)
        IF (     ibuf(iioffs(irpl)+nicomh+17)  /= imiss   )                    &
          iww   (irpl) =       ibuf (iioffs(irpl)+nicomh+17)
        IF (     ibuf(iioffs(irpl)+nicomh+18)  /= imiss   )                    &
          iwdt  (irpl) =       ibuf (iioffs(irpl)+nicomh+18)
        IF (     ibuf(iioffs(irpl)+nicomh+19)  /= imiss   )                    &
          iw1   (irpl) =       ibuf (iioffs(irpl)+nicomh+19)
        IF (     ibuf(iioffs(irpl)+nicomh+20)  /= imiss   )                    &
          iw2   (irpl) =       ibuf (iioffs(irpl)+nicomh+20)
        IF (     ibuf(iioffs(irpl)+nicomh+21)  /= imiss   )                    &
          itmxt (irpl) =       ibuf (iioffs(irpl)+nicomh+21)
        IF (     ibuf(iioffs(irpl)+nicomh+22)  /= imiss   )                    &
          itmnt (irpl) =       ibuf (iioffs(irpl)+nicomh+22)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+10) ) < rmisschk)                    &
          zgust2(irpl) =       rbuf (iroffs(irpl)+nrcomh+10)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+11) ) < rmisschk)                    &
          zrr2  (irpl) =       rbuf (iroffs(irpl)+nrcomh+11)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+12) ) < rmisschk)                    &
          zvis  (irpl) =       rbuf (iroffs(irpl)+nrcomh+12)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+13) ) < rmisschk)                    &
          zrr24 (irpl) =       rbuf (iroffs(irpl)+nrcomh+13)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+14) ) < rmisschk)                    &
          zcbase(irpl) =       rbuf (iroffs(irpl)+nrcomh+14)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+15) ) < rmisschk)                    &
          ztmax (irpl) =       rbuf (iroffs(irpl)+nrcomh+15)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+16) ) < rmisschk)                    &
          ztmin (irpl) =       rbuf (iroffs(irpl)+nrcomh+16)
      ENDIF
      IF     (kcdftyp == ncdf_synop) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+23)  /= imiss   )                    &
          iee   (irpl) =       ibuf (iioffs(irpl)+nicomh+23)
        IF (     ibuf(iioffs(irpl)+nicomh+24)  /= imiss   )                    &
          ifi   (irpl) =       ibuf (iioffs(irpl)+nicomh+24)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+17) ) < rmisschk)                    &
          zradgl(irpl) =       rbuf (iroffs(irpl)+nrcomh+17)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+18) ) < rmisschk)                    &
          zraddf(irpl) =       rbuf (iroffs(irpl)+nrcomh+18)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+19) ) < rmisschk)                    &
          zradlw(irpl) =       rbuf (iroffs(irpl)+nrcomh+19)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+20) ) < rmisschk)                    &
          zhsnow(irpl) =       rbuf (iroffs(irpl)+nrcomh+20)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+21) ) < rmisschk)                    &
          zpfi  (irpl) =       rbuf (iroffs(irpl)+nrcomh+21)
      ELSEIF (kcdftyp == ncdf_metar) THEN
        IF (     ibuf(iioffs(irpl)+nicomh+ 3)  /= imiss   )                    &
          iclev (irpl) =       ibuf (iioffs(irpl)+nicomh+ 3)
        IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 6) ) < rmisschk)                    &
          zvis  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 6)
      ENDIF
    ENDDO

    IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
        .OR. (kcdftyp == ncdf_metar)) THEN
      DO irpl = 1 , nrepsg
        maxlev = MAX( maxlev, iclev(irpl) )
      ENDDO
    ENDIF

    IF ((     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )             &
         .OR. (kcdftyp == ncdf_metar))   .AND. (maxlev >= 1)) THEN
      ALLOCATE ( iclxsg  (MAX( maxlev,1 ), nrepsg) , STAT=istat )
      ALLOCATE ( iclxcl  (MAX( maxlev,1 ), nrepsg) , STAT=istat )
      ALLOCATE ( iclxct  (MAX( maxlev,1 ), nrepsg) , STAT=istat )
      ALLOCATE ( zclxbs  (MAX( maxlev,1 ), nrepsg) , STAT=istat )
      DO irpl = 1 , nrepsg
        DO ilev = 1 , maxlev
          iclxsg (ilev,irpl) = imdi
          iclxcl (ilev,irpl) = imdi
          iclxct (ilev,irpl) = imdi
          zclxbs (ilev,irpl) = rmdi
        ENDDO
      ENDDO

      DO irpl = 1 , nrepsg
        nlev = iclev(irpl)
        DO ilev = 1 , nlev
          IF (     ibuf(iioffs(irpl)+niscal       +ilev)  /= imiss   )         &
            iclxsg (ilev,irpl) = ibuf (iioffs(irpl)+niscal       +ilev)
          IF (     ibuf(iioffs(irpl)+niscal+  nlev+ilev)  /= imiss   )         &
            iclxcl (ilev,irpl) = ibuf (iioffs(irpl)+niscal+  nlev+ilev)
          IF (     ibuf(iioffs(irpl)+niscal+2*nlev+ilev)  /= imiss   )         &
            iclxct (ilev,irpl) = ibuf (iioffs(irpl)+niscal+2*nlev+ilev)
          IF (ABS( rbuf(iroffs(irpl)+nrscal       +ilev) ) < rmisschk)         &
            zclxbs (ilev,irpl) = rbuf (iroffs(irpl)+nrscal       +ilev)
        ENDDO
      ENDDO
    ENDIF
  ENDIF   !   ((lsurfob) .AND. (nrepsg >= 1))

  IF ((kcdftyp == ncdf_nat_dwd) .AND. (nrepsg >= 1)) THEN
    DO irpl = 1 , nrepsg
      IF (     ibuf(iioffs(irpl)+nicomh+ 2)  /= imiss   )                      &
        irrdt (irpl) =       ibuf (iioffs(irpl)+nicomh+ 2)
      IF (     ibuf(iioffs(irpl)+nicomh+ 3)  /= imiss   )                      &
        irrdt2(irpl) =       ibuf (iioffs(irpl)+nicomh+ 3)
      IF (     ibuf(iioffs(irpl)+nicomh+ 4)  /= imiss   )                      &
        igudt (irpl) =       ibuf (iioffs(irpl)+nicomh+ 4)
      IF (     ibuf(iioffs(irpl)+nicomh+ 5)  /= imiss   )                      &
        igudt2(irpl) =       ibuf (iioffs(irpl)+nicomh+ 5)
      IF (     ibuf(iioffs(irpl)+nicomh+ 6)  /= imiss   )                      &
        zdd   (irpl) = REAL( ibuf (iioffs(irpl)+nicomh+ 6) , wp )
      IF (     ibuf(iioffs(irpl)+nicomh+ 7)  /= imiss   )                      &
        ivdt  (irpl) =       ibuf (iioffs(irpl)+nicomh+ 7)
      IF (     ibuf(iioffs(irpl)+nicomh+ 8)  /= imiss   )                      &
        ivtisi(irpl) =       ibuf (iioffs(irpl)+nicomh+ 8)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 1) ) < rmisschk)                      &
        zrr   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 1)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 2) ) < rmisschk)                      &
        zrr2  (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 2)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 3) ) < rmisschk)                      &
        zgust (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 3)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 4) ) < rmisschk)                      &
        zgust2(irpl) =       rbuf (iroffs(irpl)+nrcomh+ 4)
      IF (ABS( rbuf(iroffs(irpl)+nrcomh+ 5) ) < rmisschk)                      &
        zff   (irpl) =       rbuf (iroffs(irpl)+nrcomh+ 5)
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 3: Process pressure and geopotential, derive local blacklist, and
!            fill pressure and height into ODR arrays
!-------------------------------------------------------------------------------

! ---------------
! surface reports: get pressure, reference height level, and pressure code
! ---------------

  IF (lsurfob) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
! (surface reports without station altitude have been discarded previously!)
      istalt  =  NINT( osghed(nsgob,nhalt) )
! use station pressure normally
      zpuse  (irpl) = zpp(irpl)
      zzuse  (irpl) = osghed(nsgob,nhalt)
      npcode (irpl) = 1
      lpract (irpl) = .FALSE.
! use msl pressure   if station pressure is missing
!                (or if station height is missing, but this does not occur)
      IF ((zpuse(irpl) < rmdich) .AND. (zpmsl(irpl) > rmdich)) THEN
        zpuse  (irpl) = zpmsl(irpl)
        zzuse  (irpl) = c0
        npcode (irpl) = 0
! use msl pressure   if msl pressure computed from reported station pressure and
!                       temperature deviates too much from reported msl pressure
      ELSEIF (      (zpp(irpl) > rmdich) .AND. (zpmsl(irpl) > rmdich)          &
              .AND. (ztt(irpl) > rmdich)) THEN
        ztmean = ztt(irpl) + c05* 0.0065_wp* osghed(nsgob,nhalt)
!   perform this check only if sta. height is below 800 m (reporting practice!)
!                          and input values are reasonable
        IF (     (ztmean >= 210._wp) .AND. (zpmsl(irpl) >=  85000._wp)         &
           .AND. (ztmean <= 330._wp) .AND. (zpmsl(irpl) <= 108000._wp)         &
           .AND. (istalt >= -400)    .AND. (istalt <= 800)) THEN
!   simple computation of reduced msl pressure
          IF (istalt /= 0) THEN
            zpred = zpp(irpl) * exp( r_g *osghed(nsgob,nhalt) / (r_d *ztmean) )
            zplim = 60._wp + 0.125_wp* ABS( osghed(nsgob,nhalt) )
          ELSE
            zpred = zpp(irpl)
            zplim = 60._wp
          ENDIF
          IF (ABS( zpmsl(irpl) - zpred ) > zplim) THEN
            zpuse  (irpl) = zpmsl(irpl)
            zzuse  (irpl) = c0
            npcode (irpl) = 0
          ENDIF
        ENDIF
      ENDIF
! use reduced pressure if neither station or msl pressure is reported for synop
      IF ((kcdftyp == ncdf_synop) .AND. (npcode(irpl) == 1)) THEN
        IF ((zpfi(irpl) > rmdich) .AND. (ifi(irpl) /= imdi)) THEN

! bad reporting practice assumed if:
!   station altitude  >  800m             , and mean sea level is reported
!   station altitude  >  800m             , and  1000hPa level is reported
!   station altitude  > 1500m             , and   925hPa level is reported
!   station altitude  > 1700m or  <  300m , and   900hPa level is reported
!   station altitude  > 2300m or  <  800m , and   850hPa level is reported
!   station altitude  > 3700m or  < 2300m , and   700hPa level is reported
!   station altitude              < 3700m , and   500hPa level is reported
          npfi = NINT( zpfi(irpl) *0.01_wp )
          IF     (npfi ==  850) THEN
            npcode (irpl) =  2
            IF ((istalt <  800) .OR. (istalt > 2300))  lpract (irpl) = .TRUE.
          ELSEIF (npfi ==  700) THEN
            npcode (irpl) =  3
            IF ((istalt < 2300) .OR. (istalt > 3700))  lpract (irpl) = .TRUE.
          ELSEIF (npfi ==  500) THEN
            npcode (irpl) = 11
            IF  (istalt < 3700)                        lpract (irpl) = .TRUE.
          ELSEIF (npfi ==  900) THEN
            npcode (irpl) =  9
            IF ((istalt <  300) .OR. (istalt > 1700))  lpract (irpl) = .TRUE.
          ELSEIF (npfi == 1000) THEN
            npcode (irpl) = 10
            IF                       (istalt >  800)   lpract (irpl) = .TRUE.
          ELSEIF (npfi ==  925) THEN
            npcode (irpl) = 12
            IF                       (istalt > 1500)   lpract (irpl) = .TRUE.
          ENDIF
          IF (zpuse(irpl) < rmdich) THEN
            zpuse  (irpl) = zpfi (irpl)
            zzuse  (irpl) = REAL( ifi(irpl) , wp )
! use reduced pressure if red. pressure computed from reported station pressure
!                 and temperature deviates too much from reported red. pressure
!   perform this check only if input values are available and reasonable
          ELSEIF (      (.NOT. lpract(irpl)) .AND. (zpp(irpl) > rmdich)        &
                                             .AND. (ztt(irpl) > rmdich)        &
                  .AND. (ztt(irpl) >= 210._wp) .AND. (npfi >=  450)            &
                  .AND. (ztt(irpl) <= 330._wp) .AND. (npfi <= 1080)            &
                  .AND. (istalt >= -400)       .AND. (istalt <= 6000)) THEN
!   simple computation of reduced pressure
            ztmean = ztt(irpl) + c05* 0.0065_wp* (istalt - ifi(irpl))
            zpred  = zpp(irpl) * EXP( r_g *(istalt - ifi(irpl)) /(r_d *ztmean) )
            zplim  = 60._wp + 0.125_wp * REAL(ABS( istalt - ifi(irpl) ), wp)
            IF (ABS( zpfi(irpl) - zpred ) > zplim) THEN
              zpuse  (irpl) = zpfi(irpl)
              zzuse  (irpl) = REAL( ifi(irpl) , wp )
            ELSE
              npcode (irpl) = 1
              lpract (irpl) = .FALSE.
            ENDIF
          ELSE
            npcode (irpl) = 1
            lpract (irpl) = .FALSE.
          ENDIF
        ENDIF
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .AND. (npcode(irpl) == 0)               &
         .AND. (istalt /= imdi) .AND. (istalt > 800))  lpract (irpl) = .TRUE.

! get 'station height minus model orography' + land/sea obs indicator, used
! further below to decide whether different variables are active or passive
      lseaobs       = BTEST( mosghd(nsgob,nhschr),nvsebp )
      landsy (irpl) = (mosghd(nsgob,nhobtp) == nsynop) .AND. (.NOT.lseaobs)
!     fisd   (irpl) = osghed(nsgob,nhsurf) - osghed(nsgob,nhalt)
      fisd   (irpl) = osghed(nsgob,nhalt)  - osghed(nsgob,nhsurf)

! get (model-independent!) 'reference pressure' for evaluating the blacklist
!   ( = zpuse , or pressure of standard atmosphere at station altitude)
      zpref  (irpl) = zpuse(irpl)
      IF (     (zpref(irpl) < rmdich) .OR. (npcode(irpl) >= 2)                 &
          .OR. ((npcode(irpl) == 0) .AND. (istalt /= 0))) THEN
!       zpref(irpl) = r_ps(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))
        ! safety: 'std_atmosphere' will always return a value if zstalt > -800.
        zstalt = MAX( osghed(nsgob,nhalt) , -400._wp )

        CALL std_atmosphere ( 'z', zstalt, 'p', zpref(irpl), r_g, r_d, irm )
!       ===================

      ENDIF
      ! if (ps < 1 Pa) then set ps to missing to avoid division by zero later on
      IF (zpuse(irpl) < 1._wp)  zpuse(irpl) = rmdi
    ENDDO
  ENDIF

! ------------------
! AMDAR / SATOB case:   derive pressure from height (when required)
! ------------------
! Note: complete reports are rejected if both pressure and height are missing,
!                        or if pressure or height is outside required limits

  IF (lupair) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      kobtyp = mosghd(nsgob,nhobtp)
      icma   = icmaa (irpl)
      lpract (irpl) = .FALSE.
      lrej  = .FALSE.
      iactr = 0
      IF (mosghd(nsgob,nhpass) == 0) iactr = 1
! reject upper-air obs without pressure and height
      IF ((zpp(irpl) < rmdich) .AND. (zzz(irpl) < rmdich)) THEN
        ilen = 2 + istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt1
          DO icl = 1 , istrej
            outbuf(nacout+2+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
        zpuse (irpl) = rmdi
        zpref (irpl) = rmdi
        lrej = .TRUE.
! use reported pressure if present
      ELSEIF (zpp(irpl) > rmdich) THEN
        zpuse (irpl)  =  zpp(irpl)
        zpref (irpl)  =  zpp(irpl)
      ELSE
! derive pressure from 'zzz': height or flight level;
! but always use (model state independent) the standard atmosphere
! to derive 'zpref' which is used to evaluate the blacklist !
!       rzzz_grb = REAL( zzz(irpl), irealgrib )
!       CALL atmhp ('single', irm, 'US1976', rzzz_grb, rppp_grb )  ! etc.

        CALL std_atmosphere ( 'z', zzz(irpl), 'p', zpref(irpl), r_g, r_d, irm )
!       ===================

        IF (irm /= 0)  zpref (irpl) = rmdi
        IF (irm /= 0)  lrej = .TRUE.
        IF (kobtyp == nairep) THEN
          zpuse (irpl)  =  zpref(irpl)
        ! SATOB: derive pressure from height according to model atmosphere
        ELSEIF (kobtyp == nsatob) THEN
          iob    = mosghd(nsgob,nhio)
          job    = mosghd(nsgob,nhjo)
          zpuse (irpl) = f_z2p ( zzz(irpl), ke, r_hhl(iob,job,:), r_p(iob,job,:) &
                               , r_t_ll(iob,job), r_g, r_d, rmdi, .FALSE. )
!                        =====
          ! indicate the model state was used to derive pressure
          lpract (irpl) = .TRUE.
        ENDIF
      ENDIF
! discard any reports above level 'rpplim' (if.not.lrej: 'zpuse' is not missing)
      IF ((.NOT. lrej) .AND. (zpref(irpl) < rpplim))  lrej = .TRUE.
! do not reject AMDAR obs below model orography (because 'zzz' is ficticious)
!     IF ( ... (zzz(irpl) < osghed(nsgob,nhsurf))) THEN
!       zpuse (irpl) = rmdi;   lrej = .TRUE.  ; ENDIF
      mpassiv  =  0
      IF (.NOT. lrej) THEN       ! (i.e. 'zpuse' /= missing value)
! set aircraft rep. to passive if less than ( 25m or)  3hPa  \  above (model)
! set  SATOB   rep. to passive if less than (200m or) 20hPa  /   'orography'
! --> do not use AMDAR height for the orography check because it is ficticious
!     and can be well below the real airport height within anticyclones
!   Note (!): Since there is no true height of AMDAR available, AMDAR pressure
!   ----      has to be compared to model surface pressure, and this renders
!             the report state dependent on the model state, i.e. the state
!             of the AMDAR report may depend on the ensemble member!!!
        zporo = r_ps(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))
        loro  = ((kobtyp == nairep) .AND. (zpuse(irpl) > zporo - 300._wp))
        IF ((kobtyp == nsatob) .AND. (zzz(irpl) >= rmdich)) THEN
          !  use of reported height for SATOB makes 'loro' model state independ.
          loro  = (zzz(irpl) < osghed(nsgob,nhsurf)+200._wp)
        ELSEIF (kobtyp == nsatob) THEN
          !  model state is used both to derive 'zpuse' and 'zporo', and
          !  this make 'loro' almost independent from the model state
          !  (exactly independent if '2000.' was replaced by '0.')
          loro  = (zpuse(irpl) > zporo -2000._wp)
        ENDIF
!       zporo = r_ps(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))
!       IF (    ((kobtyp == nairep) .AND. (zpuse(irpl) > zporo - 300._wp))     &
!          .OR. ((kobtyp == nsatob) .AND. (zpuse(irpl) > zporo -2000._wp))) THEN
        IF (loro) THEN
          mosghd (nsgob,nhschr)  =  IBSET ( mosghd(nsgob,nhschr), nvalbp )
          mosghd (nsgob,nhschr)  =  IBSET ( mosghd(nsgob,nhschr), nvpsbp )
          mosghd (nsgob,nhflag)  =  IBSET ( mosghd(nsgob,nhflag), FL_HEIGHT )
          mpassiv = 2
          IF ((iactr == 1) .AND. (mpassiv == 2)) THEN
            neventr (nezdif,icma) = neventr(nezdif,icma) + iactr
            cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
            cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
            mosghd (nsgob,nhpass)  =  2
            iactr                  =  0
          ENDIF
        ENDIF
      ENDIF

! update statistics and event report counters if complete reports are discarded
      IF (lrej)                                                                &
        neventr (nenops,icma) = neventr (nenops,icma) + iactr
      IF ((lrej) .OR. ((mpassiv == 2) .AND. (.NOT. lverpas))) THEN
        IF (iactr == 1) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        IF (iactr == 0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
        mosghd (nsgob,nhpass)  = -1
      ENDIF
      zzuse  (irpl) = zzz(irpl)
!     npcode (irpl) = imdi
!     (it is not necessary here to set zzuse=rmdi if (zpp < rmdich)
!      since ioberr=0 is set for aircrafts and satob further below anyway
!      in order to avoid assimilating height if pressure was not observed)
    ENDDO
  ENDIF

! -------------
! scatterometer
! -------------
  IF (lscatto) THEN
    DO irpl = 1 , nrepsg
      zzuse (irpl) = c0
      zpuse (irpl) = rmdi
      zpref (irpl) = rmdi
    ENDDO
  ENDIF

!                               ---------
! auxilliary: produce 2-dim obs blacklist field (report, variable v,p,T,q)
! ------------------------------------------------------------------------
! Note: For SATOB, the 'pressure' 'pref' derived from reported height using
!       the standard atmosphere is deployed; this can deviate significantly
!       from the 'true' pressure but makes the blacklisting independent from
!       the model state (ensemble member)
!       

  DO irpl = 1 , nrepsg
    kblk (irpl,1) = 0
    kblk (irpl,2) = 0
    kblk (irpl,3) = 0
    kblk (irpl,4) = 0
    DO intv = 1 , blk_loc(irpl)% ndim
! if pressure not defined: variables with a non-zero entry in the original
!                          global blacklist are always flagged here
      IF (zpref(irpl) <= rmdich) THEN
        kblk (irpl,blk_loc(irpl)%kvar(intv)) = 1
      ELSE
! if pressure defined: variables are flagged only if pressure is in the
!                      interval given in the original blacklist
        IF (      (zpref(irpl) <= blk_loc(irpl)%plow(intv) +epsy)              &
            .AND. (zpref(irpl) >= blk_loc(irpl)%pup (intv) -epsy))             &
          kblk (irpl,blk_loc(irpl)%kvar(intv)) = 1
      ENDIF
    ENDDO
  ENDDO
  DO irpl = 1 , nrepsg
! if missing on whitelist: all variables are flagged
    IF (lblkw(irpl)) THEN
      kblk (irpl,1) = 1
      kblk (irpl,2) = 1
      kblk (irpl,3) = 1
      kblk (irpl,4) = 1
    ENDIF
  ENDDO

  IF (nrepsg >= 1) THEN
    DEALLOCATE ( blk_loc , STAT=istat )
  ENDIF

! ---------------------------------------
! pressure observation error and flagging
! ---------------------------------------

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    kobtyp = mosghd(nsgob,nhobtp)
    icma   = icmaa (irpl)
    iactr  = 0
    IF (mosghd(nsgob,nhpass) == 0)  iactr  = 1
    nflgx  (irpl)  =   0
    zoberr (irpl)  =  c0
    ioberr (irpl)  =   1
    IF (     (mosghd(nsgob,nhpass) == -1)                                      &
        .OR. (.NOT. lsurfob)) THEN
      ioberr (irpl)  =  0
      fisdnu         =  rmdi
      IF (kobtyp == nairep) THEN
        ! aircraft: 'height' is always flight level and not geometrical height !
        nflgx (irpl)  =  IBSET( nflgx(irpl), nvfbps(5) )
      ELSEIF ((lupair) .AND. (zpp  (irpl) < rmdich)                            &
                       .AND. (zpuse(irpl) > rmdich)) THEN
        nflgx (irpl)  =  IBSET( nflgx(irpl), nvfbps(6) )
      ENDIF
! if pressure missing (no flag in 'nflgx' set)
    ELSEIF (zpuse(irpl) < rmdich) THEN
      ioberr (irpl)  =  0
      fisdnu         =  rmdi
      IF (kcdftyp /= ncdf_nat_dwd) THEN
        IF (nzex(kobtyp) >= 1)  neventd (nepmis,icma) =                        &
                                neventd (nepmis,icma) + iactr
      ENDIF
    ELSEIF (lsurfob) THEN
! if reporting practice for pressure bad
      IF (lpract(irpl)) THEN
        neventd (neprac,icma) = neventd(neprac,icma) + iactr
        nflgx  (irpl)  =  IBSET( nflgx(irpl), nvfbps(5) )
        ioberr (irpl)  =  0
      ENDIF
! if pressure blacklisted
      IF (kblk(irpl,2) == 1) THEN
        IF (ioberr(irpl) == 1)  neventd (nepflg,icma) =                        &
                                neventd (nepflg,icma) + iactr
        nflgx  (irpl)  =  IBSET( nflgx(irpl), nvfbps(2) )
        ioberr (irpl)  =  0
      ENDIF
      lseaonly = (ABS(altopsu(2)) < epsy)

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
!!!!  in the nudging, but still used in the LETKF. Then, the status bit in word
!!!!  'nbserr' is set to 1 (active), but at the same time, the invalid height
!!!!  range flag is set in 'nbsflg' !!!
! (for surface reports, 'zzuse' and station height are always well-defined)
      iob    = mosghd(nsgob,nhio)
      job    = mosghd(nsgob,nhjo)
      zzkeml = c05* (r_hhl(iob,job,ke+1)+r_hhl(iob,job,ke))
      IF (zzuse(irpl) >= zzkeml-epsy) THEN
        ! --> only extrapolate obs-value from station height to reported p-level
        !     (model-values are interpolated to reported p-level)
        fisdzz = ABS( osghed(nsgob,nhalt) - zzuse(irpl) )
        fisdnu = fisdzz + fdoro(2) *(zzuse(irpl) - zzkeml)
      ELSEIF (zzuse(irpl) >= osghed(nsgob,nhalt)-epsy) THEN
        ! i.e. zzkeml > zzuse >= sta_alt
        ! --> extrapolate mod-val from model orography down to reported p-level
        !            plus obs-val from station height  up   to reported p-level
        ! --> ABS(zzkeml-zzuse)+ABS(zzuse-sta_alt) = zzkeml-sta_alt
        fisdzz = zzkeml - osghed(nsgob,nhalt)
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
        fisdzz = MAX( osghed(nsgob,nhalt) - zzuse(irpl)                        &
                    , zzkeml              - zzuse(irpl) )
        fisdnu = fisdzz
      ENDIF
! additional observation error assigned due to extrapolation
      zadder = 0.04_wp * fisdzz
! if station height or its diff. to model orography too large for surf. pressure
      IF (     ((.NOT.lseaonly) .AND. (osghed(nsgob,nhalt) > altopsu(2)))      &
          .OR. ((     lseaonly) .AND. (landsy(irpl)))                          &
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
      zoberr (irpl) = zadder * zpref(irpl) * r_g / (r_d * tmelt)
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
    osgbdy (nsgob,nbsp  )  =  zpuse (irpl)
!     (for surface reports, 'zzuse' and station height are always well-defined)
    osgbdy (nsgob,nbsz  )  =  zzuse (irpl)
    osgbdy (nsgob,nbszer)  =  zoberr(irpl)
    osgbdy (nsgob,nbsplv)  =  zpref (irpl)
    osgbdy (nsgob,nbsviz)  =  fisdnu
!   IF (lsurfob)                                                               &
    osgbdy (nsgob,nbsvip)  =  r_ps(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))
    mosgbd (nsgob,nbserr)  =  insert( 0, ioberr(irpl), nvrz )
    mosgbd (nsgob,nbsflg)  =  insert( 0, nflgx(irpl), nvfzbp )
    mosgbd (nsgob,nbslid)  =  IBSET ( 0, nvlidp(7) )
    IF (kobtyp == nairep)  mosgbd (nsgob,nbslid)  =  0
    IF (kobtyp == nsatob)  mosgbd (nsgob,nbslid)  =  0
    IF (kobtyp == nsynop)  mosgbd (nsgob,nbslid)  =  npcode(irpl)
!   IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )              &
!       .OR. (kcdftyp == ncdf_buoy ))                                          &
!     osgbdy (nsgob,nbspsl)  =  zpmsl (irpl)

    mosgbd (nsgob,nbsqcf)  =  0
    ! fill ODR station characteristics: vertical levels reported as height
    IF ((lupair) .AND. (lpract(irpl)))                                         &
      CALL MVBITS( 1 , 0 , 1 , mosghd(nsgob,nhschr), nvlzbp )
  ENDDO

! ----------------------
! 3-hour pressure change (for surface pressure quality control threshold only)
! ----------------------

  IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )                &
      .OR. (kcdftyp == ncdf_buoy )) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF ((mosghd(nsgob,nhpass) /= -1) .AND. (zdpdt(irpl) > rmdich)) THEN
        IF (ABS( zdpdt(irpl) ) < 4000._wp) THEN
          osgbdy (nsgob,nbspst)  =  zdpdt(irpl)
! if pressure tendency gross error
        ELSEIF (mosghd(nsgob,nhpass) == 0) THEN
          neventd (nepsdt,icmaa(irpl)) = neventd (nepsdt,icmaa(irpl)) + 1
          ilen = 4 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt6
            outbuf(nacout+3) = 0
            outbuf(nacout+4) = NINT( zdpdt(irpl) * 0.1_wp )
            DO icl = 1 , istrej
              outbuf(nacout+4+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 4: Determine and store the main elements in the single-level ODR
!            (temperature, humidity, wind)
!-------------------------------------------------------------------------------

! -----------
! temperature
! -----------

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF ((mosghd(nsgob,nhpass) /= -1) .AND. (kcdftyp /= ncdf_satob)             &
                                     .AND. (kcdftyp /= ncdf_nat_dwd)           &
                                     .AND. (.NOT. lscatto)) THEN
      kobtyp = mosghd(nsgob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mosghd(nsgob,nhpass) == 0)  iactr  = 1
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
! if ACARS dataset temperature error flagged high
        IF ((kcdftyp == ncdf_amdar) .OR. (kcdftyp == ncdf_modes)) THEN
!   (missing value: 3)
          nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 3, 0 )
        ELSEIF (lacars) THEN
          IF ((zttq(irpl) <= epsy) .OR. (zttq(irpl) > 99._wp)) THEN
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 3, 0 )
          ELSEIF (zttq(irpl) < c2-epsy) THEN
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 0, 0 )
            zoberr (irpl) = zoberr(irpl) + zttq(irpl)
          ELSE
            IF (ioberr(irpl) == 1)  neventd (netflg,icma) =                    &
                                    neventd (netflg,icma) + iactr
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 1, 0 )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
! if height of sensor above ground is not close to 2 m (ok: 1.5m <= zzt <= 2.5m)
        IF (lsurfob) THEN
          IF ((zzt(irpl) > rmdich) .AND. (ABS( zzt(irpl)-c2 ) > c05+epsy)) THEN
            IF (ioberr(irpl) == 1)  neventd (netflg,icma) =                    &
                                    neventd (netflg,icma) + iactr
            nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(5) )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
! if temperature blacklisted
        IF (kblk(irpl,3) == 1) THEN
          IF (ioberr(irpl) == 1)  neventd (netflg,icma) =                      &
                                  neventd (netflg,icma) + iactr
          nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(2) )
          ioberr (irpl) = 0
        ENDIF
! if station height or its diff to model orography too large for 2-m temperature
        IF (lsurfob) THEN
          lseaonly = (ABS(altopsu(3)) < epsy)
          fisdtt =  ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd(irpl) )) &
                   * fisd(irpl)
          IF (     ((.NOT.lseaonly) .AND. (osghed(nsgob,nhalt) > altopsu(3)))  &
              .OR. ((     lseaonly) .AND. (landsy(irpl)))                      &
              .OR. (fisdtt > doromx(3))) THEN
            IF (ioberr(irpl) == 1)  neventd (netalt,icma) =                    &
                                    neventd (netalt,icma) + iactr
            nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(4) )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
! if temperature gross error
        IF (     (ztt(irpl) > tmelt +60._wp)                                   &
            .OR.((ztt(irpl) > tmelt +20._wp) .AND. (zpref(irpl) < 70000._wp))  &
            .OR.((ztt(irpl) > tmelt + 5._wp) .AND. (zpref(irpl) < 50000._wp))  &
            .OR.((ztt(irpl) > tmelt - 5._wp) .AND. (zpref(irpl) < 40000._wp))  &
            .OR. (ztt(irpl) < tmelt -90._wp)) THEN
          IF (ioberr(irpl) == 1)  neventd (netext,icma) =                      &
                                  neventd (netext,icma) + iactr
          nflgx  (irpl) = IBSET( nflgx(irpl), nvfbps(3) )
          ioberr (irpl) = 0
        ENDIF
      ENDIF

! observation error for 2-m temperature
      IF (lsurfob)  zoberr (irpl) = c1 + 0.003_wp *ABS( fisd(irpl) )

! fill ODR (temperature)
      osgbdy (nsgob,nbst  ) = ztt   (irpl)
      osgbdy (nsgob,nbster) = zoberr(irpl)
      mosgbd (nsgob,nbserr) = insert( mosgbd(nsgob,nbserr), ioberr(irpl), nvrt )
      mosgbd (nsgob,nbsflg) = insert( mosgbd(nsgob,nbsflg), nflgx(irpl), nvftbp)
      IF (kobtyp == nairep)                                                    &
        mosgbd (nsgob,nbslid) = IBSET( mosgbd(nsgob,nbslid), nvlidp(9) )
    ENDIF
  ENDDO


! --------
! humidity
! --------

IF ((lexitd) .OR. (lexirh) .OR. (lexiqx)) THEN

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
!   PRINT *, 'zrh0a ',irpl, mosghd(nsgob,nhpass), zrhw(irpl)
    IF (mosghd(nsgob,nhpass) /= -1) THEN
!     PRINT *, 'zqx0  ',irpl, ystidl(irpl), zrhw(irpl), zqx(irpl), mosghd(nsgob,nhhrmn)
      kobtyp = mosghd(nsgob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mosghd(nsgob,nhpass) == 0)  iactr  = 1
      nflgx  (irpl) =  0
      zoberr (irpl) = c0
      ioberr (irpl) =  1

! reported humidity is fully discarded if negative or exactly zero
!   (sometimes reported zero values denote undefined values)
      IF (zrhw(irpl) <= 1.E-10_wp) THEN
        zrhw   (irpl)  =  rmdi
!       ioberr (irpl)  =  0   ! is set further below if humidity is missing
!       IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =       & ! also set
!                               neventd (neqflg,icma) + iactr   ! below
      ENDIF
      IF (ztd(irpl) <= 1.E-10_wp)  ztd (irpl)  =  rmdi
      IF (zqx(irpl) <= 1.E-10_wp)  zqx (irpl)  =  rmdi

! as long as reported humidity of any kind is converted into relative humidity
! in order to be processed further (as is currently the case for nudging), then:
!   mixing ratio without temperature or pressure, and
!   dewpoint temperature without temperature      are completely discarded
!   (since these quantities are required for conversion into relative humidity)
      IF ((zqx(irpl) > rmdich) .AND. (     (ztt  (irpl) < rmdich)              &
                                      .OR. (zpuse(irpl) < rmdich))) THEN
        zqx (irpl)  =  rmdi
        IF ((zrhw(irpl) < rmdich) .AND. (ztd(irpl) < rmdich)) THEN
          IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
          ioberr (irpl) = 0
        ENDIF
      ENDIF
      IF ((ztd(irpl) > rmdich) .AND.       (ztt  (irpl) < rmdich) ) THEN
        ztd (irpl)  =  rmdi
        IF ((zrhw(irpl) < rmdich) .AND. (zqx(irpl) < rmdich)) THEN
          IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
          ioberr (irpl) = 0
        ENDIF
      ENDIF

! if (relative or dewpoint) humidity missing (no flag in 'nflgx' set)
!     IF (    ((lrhtdqx) .AND. (MAX(zrhw(irpl),ztd(irpl),zqx(irpl)) <= c0   )) &
!        .OR. ((lrhtd  ) .AND. (MAX(zrhw(irpl),ztd(irpl)          ) <= c0   )) &
!        .OR. ((ltd    ) .AND. (    ztd (irpl) <= epsy  ))) THEN
      IF (    ((lrhtdqx) .AND. (MAX(zrhw(irpl),ztd(irpl),zqx(irpl)) < rmdich)) &
         .OR. ((lrhtd  ) .AND. (MAX(zrhw(irpl),ztd(irpl)          ) < rmdich)) &
         .OR. ((ltd    ) .AND. (               ztd(irpl)        < rmdich))) THEN
        IF ((ioberr(irpl) == 1) .AND. (ntdex(kobtyp) >= 1))                    &
          neventd (neqmis,icma) = neventd(neqmis,icma) + iactr
        zrhw   (irpl) = rmdi
        ztd    (irpl) = rmdi
        zqx    (irpl) = rmdi
        zrhc   (irpl) = rmdi
        ioberr (irpl) = 0
      ELSE
! if ACARS humidity error / missing measurement flagged high
        IF (kcdftyp == ncdf_amdar) THEN
!   (missing value: 3)
          nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 3, 0 )
        ELSEIF ((lacars) .OR. (kcdftyp == ncdf_modes)) THEN
          IF ((iqxq(irpl) == imdi) .OR. (iqxq(irpl) >= 63)                     &
                                   .OR. (iqxq(irpl) <  0 )) THEN
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 3, 0 )
          ELSEIF (iqxq(irpl) == 0) THEN
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 0, 0 )
          ELSE
            IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
                                    neventd (neqflg,icma) + iactr
            nflgx  (irpl) = ireplace( nflgx(irpl), nvfbps(1), nvfboc(1), 1, 0 )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
! if humidity blacklisted
        IF (kblk(irpl,4) == 1) THEN
          IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(2) )
          ioberr (irpl) = 0
        ENDIF
! if station height or its differ. to model orography too large for 2-m humidity
        IF (lsurfob) THEN
          lseaonly = (ABS(altopsu(4)) < epsy)
          fisdrh =  ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd(irpl) )) &
                   * fisd(irpl)
          IF (     ((.NOT.lseaonly) .AND. (osghed(nsgob,nhalt) > altopsu(4)))  &
              .OR. ((     lseaonly) .AND. (landsy(irpl)))                      &
              .OR. (fisdrh > doromx(4))) THEN
            IF (ioberr(irpl) == 1)  neventd (neqalt,icma) =                    &
                                    neventd (neqalt,icma) + iactr
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(4) )
            ioberr (irpl) = 0
          ENDIF
! if (upper-air) humidity not in valid height range
        ELSEIF (osgbdy(nsgob,nbsplv) < pqmin) THEN
          IF (ioberr(irpl) == 1)  neventd (neq300,icma) =                      &
                                  neventd (neq300,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(4) )
          ioberr (irpl) = 0
        ENDIF
! derived rel. humidity is set passive if temperature passive (for other reasons
! than only not in valid height range) and 'zrhw' missing
        IF (      (      (.NOT. BTEST( mosgbd(nsgob,nbserr),nvrt ))            &
                   .AND. (IBCLR( IBITS( mosgbd(nsgob,nbsflg), nvftbp, nvfaoc ) &
                               , nvfbps(4) ) > 0))                             &
            .AND. (zrhw(irpl) < rmdich)                                        &
            .AND. ((ztd(irpl) > rmdich) .OR. (zqx(irpl) > rmdich))) THEN
          IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                      &
                                  neventd (neqflg,icma) + iactr
          nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
          ioberr (irpl) = 0
        ENDIF
! if dewpoint temperature gross error, or greater than temperature+2K
!       IF (lexitd) THEN    <-- not required (if 'ztd' is always initialised)
        IF (      (ztd(irpl) > rmdich)                                         &
            .AND. (      (ztd(irpl) < tmelt -150._wp)                          &
                   .OR. ((ztd(irpl) < tmelt - 90._wp) .AND. (lsurfob))         &
                   .OR.  (ztd(irpl) > tmelt + 40._wp)                          &
                   .OR.  (ztd(irpl) > ztt(irpl)+c2))) THEN
          IF (MAX( zrhw(irpl),zqx(irpl) ) < rmdich) THEN
            IF (ioberr(irpl) == 1)  neventd (neqlow,icma) =                    &
                                    neventd (neqlow,icma) + iactr
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
            ioberr (irpl) = 0
!    (minimum ztd must be > max( b4w,b4i ) to avoid division by zero)
            ztd  (irpl)  =  MAX( ztd(irpl), b4w+c1 )
            ztd  (irpl)  =  MIN( MIN( ztd(irpl), ztt(irpl)+c2 )                &
                               , tmelt + 40._wp )
          ELSE
!    (if 'ztd' has gross error, use 'zrhw' or 'zqx' unless both are missing)
            ztd  (irpl)  =  rmdi
          ENDIF
        ENDIF
! if mixing ratio gross error  , i.e. not in (1E-10 < zqx(w) <= 0.03)
        IF (      (zqx(irpl) > rmdich)                                         &
            .AND. (     (zqx(irpl) > 0.03_wp)                                  &
                   .OR. (zqx(irpl) <= 1.E-10_wp))) THEN
          IF (MAX( zrhw(irpl),ztd(irpl) ) < rmdich) THEN
            IF (ioberr(irpl) == 1) THEN
              IF (zqx(irpl) < 1.E-10_wp) neventd (neqlow,icma) =               &
                                         neventd (neqlow,icma) + iactr
              IF (zqx(irpl) > 0.03_wp)   neventd (neqbig,icma) =               &
                                         neventd (neqbig,icma) + iactr
            ENDIF
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(3) )
            ioberr (irpl) = 0
            zqx (irpl)  =  MAX( MIN( zqx(irpl) , 0.03_wp ) , 1.E-10_wp )
          ELSE
!    (if 'zqx' has gross error, use 'zrhw' or 'ztd' unless both are missing)
            zqx (irpl)  =  rmdi
          ENDIF
        ENDIF
! if relative humidity gross error  , i.e. not in (epsy < zrhw <= 102%)
        IF (      (zrhw(irpl) > rmdich)                                        &
            .AND. (     (zrhw(irpl) > rtshlm)                                  &
                   .OR. (zrhw(irpl) <= epsy))) THEN
          IF (MAX( ztd(irpl),zqx(irpl) ) < rmdich) THEN
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
          ELSE
!    (if 'zrh' has gross error, use 'ztd' or 'zqx' unless both are missing)
            zrhw  (irpl)  =  rmdi
          ENDIF
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

  IF ((lexirh) .AND. (nrepsg > 0)) THEN

    CALL obs_rhw2rh ( madj_hum, nrepsg, ztt, zrhw                              &
                    , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhc )
!   ===============
  ENDIF

  IF ((lexitd) .AND. (nrepsg > 0)) THEN
    ALLOCATE ( zrhw1 (nrepsg) , STAT=istat )
    ALLOCATE ( zrh1  (nrepsg) , STAT=istat )

    CALL obs_td2rh  ( madj_hum, nrepsg, ztt , ztd                              &
                    , rmdi, b1, b2w, b2i, b3, b4w, b4i, tmelt , zrhw1, zrh1 )
!   ==============
  ENDIF

  IF ((lexiqx) .AND. (nrepsg > 0)) THEN
    ALLOCATE ( zrhw2 (nrepsg) , STAT=istat )
    ALLOCATE ( zrh2  (nrepsg) , STAT=istat )
    ALLOCATE ( zqvw  (nrepsg) , STAT=istat )
    ALLOCATE ( zqv   (nrepsg) , STAT=istat )
    zqvw = rmdi

    CALL obs_qx2rh ( madj_hum, nrepsg, ztt, zpuse, zqx , zqvw                  &
                   , rmdi, b1,b2w,b2i, b3,b4w,b4i, rdv, tmelt, zqv, zrhw2, zrh2)
!   ==============
  ENDIF

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF (mosghd(nsgob,nhpass) /= -1) THEN
!     IF (lacars) PRINT *, 'zqx7a  ',irpl, ystidl(irpl), zqx(irpl), zrh2(irpl)
!     PRINT *, 'zqx7b  ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)
!     kobtyp = mosghd(nsgob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mosghd(nsgob,nhpass) == 0)  iactr  = 1

! check consistency between reported and derived relative humidity values
! -----------------------------------------------------------------------
      IF ((lexitd) .AND. (lexirh)) THEN
!   if dewpoint-derived rel. humidity and reported rel. humidity differ by > 4%
        IF ((zrh1(irpl) > rmdich) .AND. (zrhc(irpl) > rmdich)) THEN
          IF (ABS( zrh1(irpl) - zrhc(irpl) ) > 0.04_wp) THEN
            IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
                                    neventd (neqflg,icma) + iactr
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
      ENDIF
      IF ((lexiqx) .AND. (lexirh)) THEN
!   if mixing-ratio-derived rel. humidity and reported rel. humid. differ by >4%
        IF ((zrh2(irpl) > rmdich) .AND. (zrhc(irpl) > rmdich)) THEN
          IF (ABS( zrh2(irpl) - zrhc(irpl) ) > 0.04_wp) THEN
            IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
                                    neventd (neqflg,icma) + iactr
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
      ENDIF
      IF ((lexiqx) .AND. (lexitd)) THEN
!   if mixing-ratio-derived and dewpoint-derived rel. humidity differ by >4%
        IF ((zrh2(irpl) > rmdich) .AND. (zrh1(irpl) > rmdich)) THEN
          IF (ABS( zrh2(irpl) - zrh1(irpl) ) > 0.04_wp) THEN
            IF (ioberr(irpl) == 1)  neventd (neqflg,icma) =                    &
                                    neventd (neqflg,icma) + iactr
            nflgx  (irpl) = IBSET ( nflgx(irpl), nvfbps(5) )
            ioberr (irpl) = 0
          ENDIF
        ENDIF
      ENDIF
!     PRINT *, 'zqx8   ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)

! if relative humidity is not reported ('zrhc' = missing value)
! then use the dewpoint- or mixing-ratio-derived value 'zrh1' resp. 'zrh2'
! ------------------------------------------------------------------------
      IF (lexitd) THEN
!       PRINT *, 'zrhtd ',zrh1(irpl)
        IF ((zrh1(irpl) > rmdich) .AND. (zrhc(irpl) < rmdich)) THEN
          zrhc  (irpl)  =  zrh1 (irpl)
          zrhw  (irpl)  =  zrhw1(irpl)
!    (if 'zrhw1' has gross error, use 'zrh2' if well-defined)
          IF ((zrhw1(irpl) > rtshlm) .AND. (zqx(irpl) > rmdich))               &
            zrhc  (irpl)  =  rmdi
        ENDIF
      ENDIF
      IF (lexiqx) THEN
!       PRINT *, 'zrhqx ',zrh2(irpl)
        IF ((zrh2(irpl) > rmdich) .AND. (zrhc(irpl) < rmdich)) THEN
          zrhc  (irpl)  =  zrh2 (irpl)
          zrhw  (irpl)  =  zrhw2(irpl)
        ENDIF
      ENDIF
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
      IF ((zrhc(irpl) >= rhtsat(ixrhs)) .AND. (zrhc(irpl) < c1-epsy)) THEN
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
      IF (zrhc(irpl) >= rhtsat(ixrhs)) THEN
        zrhc   (irpl)  =  c1
      ENDIF
!     IF (lacars) PRINT *, 'zqx9a  ',irpl, ystidl(irpl), zqx(irpl), zrh2(irpl)
!     PRINT *, 'zqx9b  ',irpl, ystidl(irpl), zrhc(irpl), ioberr(irpl)

! observation error for 2-m humidity
      IF (lsurfob) THEN
        zoberr (irpl) =  0.01_wp *rherr1  + 0.0004_wp *ABS( fisd(irpl) )
      ENDIF
!     PRINT *, 'zrh4 ',irpl, zrhc(irpl)
!     PRINT *, 'zrh4  ',irpl, ystidl(irpl), zrhc(irpl)

! fill ODR (humidity)
! -------------------
      osgbdy (nsgob,nbsrh ) = zrhc  (irpl)
      osgbdy (nsgob,nbsqer) = zoberr(irpl)
      IF ((zrhc(irpl) > rmdich) .AND. (zrhw(irpl) > rmdich)) THEN
        osgbdy (nsgob,nbsdrh) = zrhc (irpl)  -  zrhw(irpl)
      ELSEIF (zrhc(irpl) > rmdich) THEN
        osgbdy (nsgob,nbsdrh) = c0
      ENDIF
      mosgbd (nsgob,nbserr) = insert( mosgbd(nsgob,nbserr), ioberr(irpl), nvrq )
      mosgbd (nsgob,nbsflg) = insert( mosgbd(nsgob,nbsflg), nflgx(irpl), nvfqbp)
!     IF (kobtyp == nairep)                                                    &
!       mosgbd (nsgob,nbslid) = IBSET( mosgbd(nsgob,nbslid), nvlidp(9) )
    ENDIF
  ENDDO

  !    fill dewpoint temperature into ODR
  IF (((lexitd) .OR. (lexirh) .OR. (lexiqx)) .AND. (nrepsg > 0)) THEN

    CALL obs_rh2td ( nrepsg, zrhc, ztt, rmdi, b1, b2w, b3, b4w , ztd )
!   ==============
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1)                                          &
        osgbdy (nsgob,nbstd ) = ztd (irpl)
    ENDDO
  ENDIF

  IF ((lexitd) .AND. (nrepsg > 0)) THEN
    DEALLOCATE ( zrhw1   , STAT=istat )
    DEALLOCATE ( zrh1    , STAT=istat )
  ENDIF
  IF ((lexiqx) .AND. (nrepsg > 0)) THEN
    DEALLOCATE ( zrhw2   , STAT=istat )
    DEALLOCATE ( zrh2    , STAT=istat )
    DEALLOCATE ( zqvw    , STAT=istat )
    DEALLOCATE ( zqv     , STAT=istat )
  ENDIF

ENDIF    !    ((lexitd) .OR. (lexirh) .OR. (lexiqx))

! ---------------
! horizontal wind
! ---------------

!       zpref(irpl) = r_ps(mosghd(nsgob,nhio),mosghd(nsgob,nhjo))

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF (mosghd(nsgob,nhpass) /= -1) THEN
      kobtyp = mosghd(nsgob,nhobtp)
      icma   = icmaa (irpl)
      iactr  = 0
      IF ((mosghd(nsgob,nhpass) == 0) .AND. (kcdftyp /= ncdf_nat_dwd)) iactr = 1
      nflgx  (irpl)  =   0
      zoberr (irpl)  =  c0
      ioberr (irpl)  =   1

! if wind speed or direction missing (no flag in 'nflgx' set)
      IF ((zff(irpl) < rmdich) .OR. (zdd(irpl) < rmdich)) THEN
        IF (nuex(kobtyp) >= 1) THEN
          IF (zff(irpl) <= rmdich) neventd (nefmis,icma) =                     &
                                   neventd (nefmis,icma) + iactr
          IF (zdd(irpl) <= rmdich) neventd (nedmis,icma) =                     &
                                   neventd (nedmis,icma) + iactr
        ENDIF
        zff    (irpl)  =  rmdi
        zdd    (irpl)  =  rmdi
        ioberr (irpl)  =  0
      ELSE
! if time significance not 'time averaged' or time period > 60 minutes
!                                   (usually, time period should be 10 minutes)
        IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )          &
            .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_nat_dwd)) THEN
          IF (     (      (ivtisi(irpl) /= imdi) .AND. (ivtisi(irpl) /= 31)    &
                    .AND. (ivtisi(irpl) /= 2))                                 &
              .OR. (      (ivdt  (irpl) /= imdi)                               &
                    .AND. (ABS(ivdt(irpl)) > 60))) THEN
            IF (ioberr(irpl) == 1)  neventd (nefflg,icma) =                    &
                                    neventd (nefflg,icma) + iactr
            nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(5))
            ioberr (irpl)  =  0
          ENDIF
        ENDIF
! if wind blacklisted
        IF (kblk(irpl,1) == 1) THEN
          IF (ioberr(irpl) == 1)  neventd (nedflg,icma) =                      &
                                  neventd (nedflg,icma) + iactr
          IF (ioberr(irpl) == 1)  neventd (nefflg,icma) =                      &
                                  neventd (nefflg,icma) + iactr
          nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(2) )
          ioberr (irpl)  =  0
        ENDIF
        IF (lsurfob) THEN
!AR:    Caution: in MEC, r_hhl is not a global field and r_z0 not read so far
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF (.FALSE.) THEN
          ! COSMO model: if roughness length and/or Laplacian of orography
          !              too large for 10-m wind
            iob    = mosghd(nsgob,nhio)
            job    = mosghd(nsgob,nhjo)
            zlaplac =   r_hhl(iob-1,job  ,ke+1) + r_hhl(iob  ,job-1,ke+1)      &
                      + r_hhl(iob+1,job  ,ke+1) + r_hhl(iob  ,job+1,ke+1)      &
                      - 4._wp* r_hhl(iob,job,ke+1)
            IF (r_z0(iob,job) < 0.1_wp)  zlaplac = zlaplac * 0.9_wp
            loro  =       ( zlaplac      > zlimv10(1))                         &
                     .OR. (-zlaplac      > zlimv10(2))                         &  
                     .OR. (r_z0(iob,job) > zlimv10(3))
!           IF (loro) PRINT '("UV10M oro ",A,1X,F7.0,F6.2,2F7.1,2F7.2)'        &
!                          , ystidl(irpl) , zlaplac, r_z0(iob,job)             &
!                          , r_hhl(iob,job,ke+1) , osghed(nsgob,nhalt)         &
!                          , osghed(nsgob,nhjlat), osghed(nsgob,nhilon)
          ENDIF
          ! DACE (MEC / ICON): if std. dev. of sub-grid scale orography
          !                    too large for 10-m wind
          loro  =  (osghed(nsgob,nhssos) > zlimv10(1))
!       ENDIF
! if station height or its diff to model orography too large for 10-m wind
!       IF (lsurfob) THEN
          lseaonly = (ABS(altopsu(1)) < epsy)
          fisduv =  ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd(irpl) )) &
                   * fisd(irpl)
          IF (     ((.NOT.lseaonly) .AND. (osghed(nsgob,nhalt) > altopsu(1)))  &
              .OR. ((     lseaonly) .AND. (landsy(irpl)))                      &
              .OR. (fisduv > doromx(1)) .OR. (loro)) THEN
            IF (ioberr(irpl) == 1)  neventd (nevalt,icma) =                    &
                                    neventd (nevalt,icma) + iactr
            nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(4) )
            iob    = mosghd(nsgob,nhio)
            job    = mosghd(nsgob,nhjo)
!           IF (ioberr(irpl) == 1)                                             &
!             WRITE( 0,'("UV10M zdz ",A,L2,3F7.0,2F7.2)' )  ystidl(irpl)       &
!                  , loro, osghed(nsgob,nhssos), osghed(nsgob,nhalt), fisduv   &
!                        , osghed(nsgob,nhjlat), osghed(nsgob,nhilon)
            ioberr (irpl)  =  0
          ENDIF
        ENDIF
! if wind speed or direction gross error
        IF (      (zdd(irpl) > 360.1_wp) .OR. (zdd(irpl) < c0-epsy)            &
            .OR.  (zff(irpl) < -epsy)    .OR. (zff(irpl) > 150._wp)            &
            .OR. ((zff(irpl) >  90.0_wp) .AND. (lsurfob))       ) THEN
          IF (ioberr(irpl) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(3) )
          ioberr (irpl)  =  0
        ENDIF
! zero wind speed (not only of buoy or scatterometer) as bad reporting practice
!       IF ((zff(irpl) <= epsy) .AND. (     (kobtyp == ndribu)                 &
!                                      .OR. (kobtyp == nscatt))) THEN
        IF (zff(irpl) <= epsy) THEN
          IF (ioberr(irpl) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(5) )
          ioberr (irpl)  =  0
        ENDIF
! for DWD national reports with character station ID, do not use wind actively
        IF ((kcdftyp == ncdf_nat_dwd) .AND. (ystidl(irpl) (1:2) /= '10')) THEN
          IF (ioberr(irpl) == 1)  neventd (nefneg,icma) =                      &
                                  neventd (nefneg,icma) + iactr
          nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(5) )
          ioberr (irpl)  =  0
        ENDIF
! aircraft roll angle as revised quality flag for wind (2 bits)
        IF (     (kcdftyp == ncdf_amdar) .OR. (lacars)                         &
            .OR. (kcdftyp == ncdf_modes)) THEN
!         IF (iroll(irpl) == 1) THEN
!           nflgx  (irpl)  =  IBSET ( nflgx(irpl), nvfbps(5) )
!           ioberr (irpl)  =  0
!         ENDIF
          IF ((iroll(irpl) < 0) .OR. (iroll(irpl) > nibits(nvfboc(1))))        &
            iroll (irpl) = nibits(nvfboc(1))
          nflgx (irpl)  =  ireplace( nflgx(irpl), nvfbps(1), nvfboc(1)         &
                                   , iroll(irpl), 0 )
          IF (iroll(irpl) == 1)  ioberr (irpl)  =  0
! aircraft roll angle in station characteristics header word ('nvaaoc' bits)
          IF (iroll(irpl) == nibits(nvfboc(1)))  iroll (irpl) = nibits(nvaaoc)
          mosghd (nsgob,nhschr) = ireplace( mosghd(nsgob,nhschr)               &
                                          , nvaabp, nvaaoc, iroll(irpl), 0 )
        ENDIF
!       IF ((lsurfob) .AND. (ioberr(irpl) >= 1))                               &
!           PRINT '("UV10M act ",A,1X,F7.0,F6.2,2F7.1,2F7.2)', ystidl(irpl)    &
!                 , zlaplac, r_z0(iob,job), r_hhl(iob,job,ke+1)                &
!              , osghed(nsgob,nhalt), osghed(nsgob,nhjlat), osghed(nsgob,nhilon)
      ENDIF
    ENDIF
  ENDDO

! Transformation of wind speed and direction to wind components,
! within COSMO (not DACE) also convert to the rotated model coordinate system
! ---------------------------------------------------------------------------
  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF ((zff(irpl) > rmdich) .AND. (mosghd(nsgob,nhpass) /= -1)) THEN
      zuu (irpl,1)  =  zff(irpl) * (-SIN( zdd(irpl) *r_degrad ))
      zvv (irpl,1)  =  zff(irpl) * (-COS( zdd(irpl) *r_degrad ))
    ELSE
      zuu (irpl,1)  =  c0
      zvv (irpl,1)  =  c0
    ENDIF
    zrlat (irpl,1)  =  osghed(nsgob,nhjlat)
    zrlon (irpl,1)  =  osghed(nsgob,nhilon)
  ENDDO

  IF (nrepsg > 0) THEN
#ifdef __COSMO__
    CALL uv2uvrot_vec ( zuu, zvv, zrlat, zrlon, r_pollat, r_pollon, nrepsg, 1 )
!   =================
#endif
  ENDIF

! fill ODR (horizontal wind)
! --------
!NEC$ nomove
  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF (mosghd(nsgob,nhpass) /= -1) THEN
      kobtyp = mosghd(nsgob,nhobtp)
      IF (zff(irpl) <= rmdich)  zuu (irpl,1)  =  rmdi
      IF (zff(irpl) <= rmdich)  zvv (irpl,1)  =  rmdi
      osgbdy (nsgob,nbsu  ) = zuu   (irpl,1)
      osgbdy (nsgob,nbsv  ) = zvv   (irpl,1)
      osgbdy (nsgob,nbsuer) = zoberr(irpl)
      mosgbd (nsgob,nbserr) = insert( mosgbd(nsgob,nbserr), ioberr(irpl), nvru )
      mosgbd (nsgob,nbsflg) = insert( mosgbd(nsgob,nbsflg), nflgx(irpl), nvfubp)
      osgbdy (nsgob,nbsff ) = zff   (irpl)
      osgbdy (nsgob,nbsdd ) = zdd   (irpl)
      IF (kobtyp == nairep)                                                    &
        mosgbd (nsgob,nbslid) = IBSET( mosgbd(nsgob,nbslid), nvlidp(8) )
      IF (kobtyp == nsatob)                                                    &
        mosgbd (nsgob,nbsqal) = iqual(irpl)
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Section 5: Determine and store the additional weather-related elements
!            (precipitation, visibility, cloud, weather, gusts, etc)
!-------------------------------------------------------------------------------

! initialise group words

  IF (lsurfob) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      mosgbd (nsgob,nbscwg)  =  imdi
      mosgbd (nsgob,nbswwe)  =  imdi
      mosgbd (nsgob,nbsclg)  =  imdi
      mosgbd (nsgob,nbscl1)  =  imdi
      mosgbd (nsgob,nbscl2)  =  imdi
      mosgbd (nsgob,nbscl3)  =  imdi
      mosgbd (nsgob,nbscl4)  =  imdi
    ENDDO
  ENDIF

! ---------------------------------------------------------------------
! aircraft: degree of turbulence, max. derived equivalent vertical gust
! ---------------------------------------------------------------------

  IF ((kcdftyp == ncdf_amdar) .OR. (lacars)) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
! degree of turbulence (WMO Table 011031)
      mosgbd (nsgob,nbstur)  =  imdi
      IF ((iturb(irpl) < 15) .AND. (iturb(irpl) >= 0))                         &
        mosgbd (nsgob,nbstur)  =  iturb(irpl)
! maximum derived equivalent vertical gust
      osgbdy (nsgob,nbsfgv)  =  zvgust(irpl)
!     PRINT *,'zvgust ', irpl, ystidl(irpl), mosgbd(nsgob,nbstur), osgbdy(nsgob,nbsfgv)
    ENDDO
  ENDIF

! -------------------------
! max. 10-m wind gust speed
! -------------------------

  IF (lsurfob) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
! if max. wind gust speed gross error
      IF ((zgust(irpl) > 90._wp) .OR. (zgust(irpl) < c0-epsy))                 &
        zgust (irpl) = rmdi
      IF ((mosghd(nsgob,nhpass) /= -1) .AND. (igudt(irpl) /= imdi)             &
                                       .AND. (zgust(irpl) > rmdich)) THEN
! max. wind gust speed over 1, 3, or 6 hours
        IF (ABS(igudt(irpl)) ==  60)  osgbdy (nsgob,nbsfg1) = zgust(irpl)
        IF (ABS(igudt(irpl)) == 180)  osgbdy (nsgob,nbsfg3) = zgust(irpl)
        IF (ABS(igudt(irpl)) == 360)  osgbdy (nsgob,nbsfg6) = zgust(irpl)
      ENDIF

      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_nat_dwd)) THEN
! if max. wind gust speed gross error
        IF ((zgust2(irpl) > 90._wp) .OR. (zgust2(irpl) < c0-epsy))             &
          zgust2 (irpl) = rmdi
! max. wind gust speed over 1 hour or 6 hours
        IF ((mosghd(nsgob,nhpass) /= -1) .AND. (igudt2(irpl) /= imdi)          &
                                         .AND. (zgust2(irpl) > rmdich)) THEN
          IF (ABS(igudt2(irpl)) ==  60)  osgbdy (nsgob,nbsfg1) = zgust2(irpl)
          IF (ABS(igudt2(irpl)) == 180)  osgbdy (nsgob,nbsfg3) = zgust2(irpl)
          IF (ABS(igudt2(irpl)) == 360)  osgbdy (nsgob,nbsfg6) = zgust2(irpl)
! shorter periods (e.g. over 1 hour) are preferred over longer periods
!  previously:  igudtl = ABS( igudt(irpl) )  (if MOD(ABS(igudt(irpl)),60) == 0))
!         IF (      (MOD( ABS(igudt2(irpl)), 60 ) == 0)                        &
!             .AND. ((ABS(igudt2(irpl)) < igudtl) .OR. (igudtl == imdi))) THEN
!           igudtl                = ABS( igudt2(irpl) )
!           osgbdy (nsgob,nbsfg1) =      zgust2(irpl)
!         ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! -------------
! precipitation : get sum over 1 hour, 6 hours, and 24 or preferably 12 hours
! -------------

  IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )                &
      .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_nat_dwd)) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1) THEN
        zrr1  = rmdi
        zrr3  = rmdi
        zrr6  = rmdi
        zr12  = rmdi
        zr24  = rmdi
! get 24-h precip (from MRR24)
        IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
          IF ((zrr24(irpl) > rmdich) .AND. (zrr24(irpl) > c0-epsy))            &
              zrr24    =      zrr24 (irpl)
        ENDIF
        IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )          &
            .OR. (kcdftyp == ncdf_nat_dwd)) THEN
          IF ((zrr2(irpl) > c0-epsy) .AND. (irrdt2(irpl) /= imdi)) THEN
            IF           (ABS(irrdt2(irpl)) ==  1) THEN
              zrr1     =      zrr2  (irpl)
            ELSEIF       (ABS(irrdt2(irpl)) ==  3) THEN
              zrr3     =      zrr2  (irpl)
            ELSEIF       (ABS(irrdt2(irpl)) ==  6) THEN
              zrr6     =      zrr2  (irpl)
            ELSEIF       (ABS(irrdt2(irpl)) == 12) THEN
              zr12     =      zrr2  (irpl)
            ELSEIF       (ABS(irrdt2(irpl)) == 24) THEN
              zr24     =      zrr2  (irpl)
            ENDIF
          ENDIF
        ENDIF
        IF ((zrr (irpl) > c0-epsy) .AND. (irrdt (irpl) /= imdi)) THEN
          IF           (ABS(irrdt (irpl)) ==  1) THEN
            zrr1     =      zrr   (irpl)
          ELSEIF       (ABS(irrdt (irpl)) ==  3) THEN
            zrr3     =      zrr   (irpl)
          ELSEIF       (ABS(irrdt (irpl)) ==  6) THEN
            zrr6     =      zrr   (irpl)
          ELSEIF       (ABS(irrdt (irpl)) == 12) THEN
            zr12     =      zrr   (irpl)
          ELSEIF       (ABS(irrdt (irpl)) == 24) THEN
            zr24     =      zrr   (irpl)
          ENDIF
        ENDIF

! if precipitation gross error
        IF (     (zrr1  >  99._wp)                                             &
            .OR. (zrr3  > 199._wp)                                             &
            .OR. (zrr6  > 399._wp)                                             &
            .OR. (zr12  > 599._wp)                                             &
            .OR. (zr24  > 699._wp)) THEN
          IF (mosghd(nsgob,nhpass) == 0)                                       &
            neventd (nerlim,icmaa(irpl)) = neventd (nerlim,icmaa(irpl)) + 1
          ilen = 4 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt2
            IF     (zrr1  >  99._wp) THEN
              outbuf(nacout+3) = INT(zrr1 *100._wp)
              outbuf(nacout+4) = 1
            ELSEIF (zrr3  > 199._wp) THEN
              outbuf(nacout+3) = INT(zrr3 *100._wp)
              outbuf(nacout+4) = 3
            ELSEIF (zrr6  > 399._wp) THEN
              outbuf(nacout+3) = INT(zrr6 *100._wp)
              outbuf(nacout+4) = 6
            ELSEIF (zr12  > 599._wp) THEN
              outbuf(nacout+3) = INT(zr12 *100._wp)
              outbuf(nacout+4) = 12
            ELSEIF (zr24  > 699._wp) THEN
              outbuf(nacout+3) = INT(zr24 *100._wp)
              outbuf(nacout+4) = 24
            ENDIF
            DO icl = 1 , istrej
              outbuf(nacout+4+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
          IF (zrr1  >  99._wp)  zrr1 = rmdi
          IF (zrr3  > 199._wp)  zrr3 = rmdi
          IF (zrr6  > 399._wp)  zrr6 = rmdi
          IF (zr12  > 599._wp)  zr12 = rmdi
          IF (zr24  > 699._wp)  zr24 = rmdi
        ENDIF

! write precipitation to ODR
        osgbdy (nsgob , nbsrr1) = zrr1
        osgbdy (nsgob , nbsrr3) = zrr3
        osgbdy (nsgob , nbsrr6) = zrr6
        osgbdy (nsgob , nbsr12) = zr12
        osgbdy (nsgob , nbsr24) = zr24
        nrtr  =  0
        IF (zr12 > rmdich)  nrtr = 2
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nrtrbp, nrtroc, nrtr , 0 )
      ENDIF
    ENDDO
  ENDIF

! ---------------------
! horizontal visibility --> fog
! -----------------------------

  IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )                &
      .OR. (kcdftyp == ncdf_metar)) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      kfog (irpl) = 0
      IF ((mosghd(nsgob,nhpass) /= -1) .AND. (zvis(irpl) > rmdich)) THEN
        osgbdy (nsgob,nbsvis) = MIN( zvis(irpl) , 99999._wp )
        ! pre-set fog indicator
        IF (zvis(irpl) < vfoglim)  kfog (irpl) = 2
      ENDIF
    ENDDO
  ENDIF

! ---------------
! present weather --> fog
! -----------------------

  IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN

    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1) THEN
        !   in entry 'nbscwg', WW has only nvwwoc = 7 bits
        IF ((iww(irpl) == imdi) .OR.(iww(irpl) >= 100) .OR.(iww(irpl) < 0)) THEN
          nww = nibits(nvwwoc)
        ELSE
          nww = iww(irpl)
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvwwbp, nvwwoc, nww , 0 )

        ! new weather group word: WW has nvw0oc = 9 bits
        IF (     (iww(irpl) == imdi) .OR. (iww(irpl) > nibits(nvw0oc))         &
            .OR. (iww(irpl) < 0)) THEN
          nww = nibits(nvw0oc)
        ELSE
          nww = iww(irpl)
        ENDIF
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nvw0bp, nvw0oc, nww , 0 )

        ! pre-set flag for fog
        !   present weather is fog (excluding sky visible)
        IF (     ((nww >=  43) .AND. (nww <=  49) .AND. (mod( nww,2 ) == 1))   &
            .OR. ((nww >= 132) .AND. (nww <= 135)) .OR. (nww == 130)           &
            .OR. ((nww >= 247) .AND. (nww <= 249))) THEN
           kfog (irpl) = 3
        !   dust(storm), sand(storm), smoke, (blowing) snow, haze, volcanic ash
        !     04-09 haze, dust,sand, smoke
        !     30-39 duststorm, sandstorm, drifting or blowing snow
        !     104, 105, 127-129 haze, smoke, dust, snow
        !     204, 206-208, 210 volcan. ash, dust haze, blowing spray, snow haze
        !     230, 239 duststorm, sandstorm, blowing snow
        ELSEIF (     ((nww >=   4) .AND. (nww <=   9))                         &
                .OR. ((nww >=  30) .AND. (nww <=  39))                         &
                .OR. ((nww >= 104) .AND. (nww <= 105))                         &
                .OR. ((nww >= 127) .AND. (nww <= 129)) .OR. (nww == 230)       &
                .OR. ((nww >= 204) .AND. (nww <= 210)) .OR. (nww == 239)) THEN
           kfog (irpl) = -1
        !   in case of heavy precip, it is unclear whether
        !   (moderately) low visibility is due to precip or fog
        ELSEIF (      (kfog(irpl)  == 2) .AND. (zvis(irpl) > 100._wp)          &
                .AND. (     ((nww >=  50) .AND. (nww <=  99))                  &
                       .OR. ((nww >= 150) .AND. (nww <= 199)))) THEN
           SELECT CASE (nww)
           CASE (54, 55, 57, 59, 64, 65, 67, 69, 74, 75, 81, 82, 84, 86, 88, 90&
            ,97 ,98, 99,153,156,158,163,166,168,173,176,183,184,187,189,195,196)
              kfog (irpl) = 1
           END SELECT
        ENDIF
        IF ((kfog(irpl) == -1) .OR. (kfog(irpl) == 1)) THEN
          ilen = 6 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt20
            outbuf(nacout+3) = nww
            outbuf(nacout+4) = -1
            IF (osgbdy(nsgob,nbsvis) >= -2._wp )                               &
              outbuf(nacout+4) = NINT( osgbdy(nsgob,nbsvis) )
            outbuf(nacout+5) = 0
            outbuf(nacout+6) = 0
            DO icl = 1 , istrej
              outbuf(nacout+6+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
        ENDIF

! ------------
! past weather
! ------------

        ! valid values are between 0 and 19
        IF (     (iw1(irpl) == imdi) .OR. (iw1(irpl) > nibits(nvw1oc))         &
            .OR. (iw1(irpl) < 0) .OR. (iw1(irpl) >= 20)) THEN
          nww = nibits(nvw1oc)
        ELSE
          nww = iw1(irpl)
        ENDIF
        IF (     (iw2(irpl) == imdi) .OR. (iw2(irpl) > nibits(nvw1oc))         &
            .OR. (iw2(irpl) < 0) .OR. (iw2(irpl) >= 20)) THEN
          nw2 = nibits(nvw1oc)
        ELSE
          nw2 = iw2(irpl)
        ENDIF
        IF (nww == nibits(nvw1oc)) THEN
          nww = nw2
        ELSEIF (nw2 /= nibits(nvw1oc)) THEN
        ! use 'nw2' if it has a higher value (modulo 10), or
        !           if 'nw1=nww' is between 0 and 9 and 'nw2'='nw1'+10
          IF ((MOD( nw2,10 ) > MOD( nww,10 )) .OR. (nw2 == nww+10)) nww = nw2
        ENDIF
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nvw1bp, nvw1oc, nww , 0 )
        ! time range of past weather
        nrtr    =  0
        IF (ABS(iwdt(irpl)) == 12)  nrtr = 2
        IF (ABS(iwdt(irpl)) == 24)  nrtr = 4
        IF (ABS(iwdt(irpl)) ==  6)  nrtr = 1
        IF (ABS(iwdt(irpl)) ==  1)  nrtr = 5
        IF (ABS(iwdt(irpl)) ==  2)  nrtr = 6
        IF (ABS(iwdt(irpl)) ==  3)  nrtr = 7
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nvwtbp, nvwtoc, nrtr , 0 )
      ENDIF
    ENDDO

! -------------------------
! general cloud group (GCG) --> cloud type, significance, base height, amount
! ---------------------------------------------------------------------------

    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1) THEN

  !       ----------
  ! GCG:  cloud type
  ! ----------------

        ! low cloud type
        ncc = iccl(irpl)
        IF ((ncc < 0) .OR. (ncc > nibits(nxctoc)))  ncc = nibits(nxctoc)
        mosgbd (nsgob,nbsclg)  =  ireplace ( mosgbd(nsgob,nbsclg)              &
                                           , nctlbp, nxctoc, ncc , 0 )
        IF ((iccl(irpl) >= 30) .AND. (iccl(irpl) <= 39)) THEN
          ncc = iccl(irpl) - 30
        !   clouds invisible due to darkness, fog, blowing dust or sand, etc.
        ELSEIF ((iccl(irpl) == 59) .OR. (iccl(irpl) == 62)) THEN
          ncc = 10
        ELSE
          ncc = nibits(nvcloc)
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvclbp, nvcloc, ncc , 0 )
        ! middle cloud type
        ncc = iccm(irpl)
        IF ((ncc < 0) .OR. (ncc > nibits(nxctoc)))  ncc = nibits(nxctoc)
        mosgbd (nsgob,nbsclg)  =  ireplace ( mosgbd(nsgob,nbsclg)              &
                                           , nctmbp, nxctoc, ncc , 0 )
        IF ((iccm(irpl) >= 20) .AND. (iccm(irpl) <= 29)) THEN
          ncc = iccm(irpl) - 20
        ELSEIF ((iccm(irpl) == 59) .OR. (iccm(irpl) == 61)) THEN
          ncc = 10
        ELSE
          ncc = nibits(nvcmoc)
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvcmbp, nvcmoc, ncc , 0 )
        ! high cloud type
        ncc = icch(irpl)
        IF ((ncc < 0) .OR. (ncc > nibits(nxctoc)))  ncc = nibits(nxctoc)
        mosgbd (nsgob,nbsclg)  =  ireplace ( mosgbd(nsgob,nbsclg)              &
                                           , ncthbp, nxctoc, ncc , 0 )
        IF ((icch(irpl) >= 10) .AND. (icch(irpl) <= 19)) THEN
          ncc = icch(irpl) - 10
        ELSEIF ((icch(irpl) == 59) .OR. (icch(irpl) == 60)) THEN
          ncc = 10
        ELSE
          ncc = nibits(nvchoc)
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvchbp, nvchoc, ncc , 0 )

  !       ------------------------------
  ! GCG:  vertical significance of cloud   (WMO Table 008002)
  ! ------------------------------------

        IF (iclsig(irpl) == 62) THEN
          ncsig = iclsig(irpl)
        ELSEIF (     (iclsig(irpl) == imdi)                                    &
                .OR. (iclsig(irpl) <  0) .OR. (iclsig(irpl) >= 12)) THEN
          ncsig = nibits(nxsgoc)
        ELSE
          ncsig = iclsig(irpl)
        ENDIF
        mosgbd (nsgob,nbsclg)  =  ireplace ( mosgbd(nsgob,nbsclg)              &
                                           , nxsgbp, nxsgoc, ncsig, 0 )
        iclsig(irpl) = ncsig

  !       -----------------
  ! GCG:  total cloud cover (amount)
  ! --------------------------------

        ! convert % into key code (octas), see WMO table 020011 (cloud amount)
        IF (     (iclct(irpl) == imdi)                                         &
            .OR. (iclct(irpl) < 0) .OR. (iclct(irpl) > 113)) THEN
          nclct = nibits(nvnoc)
          iclct(irpl) = 127
        ELSEIF (iclct(irpl) == 113) THEN
          nclct = 9
        ELSEIF (iclct(irpl) >= 101) THEN
          nclct = nibits(nvnoc)
        ELSEIF (iclct(irpl) == 100) THEN
          nclct = 8
        ELSEIF (iclct(irpl) ==   0) THEN
          nclct = 0
!       ELSEIF (iclct(irpl) <=  15) THEN
!         nclct = 1
!       ELSEIF (iclct(irpl) <=  35) THEN
!         nclct = 2
!       ELSEIF (iclct(irpl) >=  85) THEN
!         nclct = 7
!       ELSEIF (iclct(irpl) >=  65) THEN
!         nclct = 6
!       ELSE
!         nclct = NINT( iclct(irpl) * 0.1_wp ) - 1
        ELSE
          nclct = MIN( MAX( NINT( iclct(irpl) / 12.5_wp ) , 1 ) , 7 )
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvnbp, nvnoc, nclct , 0 )
        IF (nclct == nibits(nvnoc))  nclct = -1  ! more convenient missing value
        !   hereafter, -1 <= nclct <= 9

  !       -----------------
  ! GCG:  cloud base height (CBH)
  ! -----------------------------
        ! 1) convert into integer, missing value = -1
        !    WMO id. 020013 in BUFR: 11 bits, scale 0.1 m -> max. value: 20920
!       IF ((zcbase(irpl) > -0.1_wp) .AND. (zcbase(irpl) <= cbh_clr+0.5)) THEN
        IF ((zcbase(irpl) > -0.1_wp) .AND. (zcbase(irpl) <= 20920.1_wp)) THEN
          icbh (irpl)  = MIN( NINT( zcbase(irpl) ) , cbh_clr )
        ELSE
          icbh (irpl)  = -1
        ENDIF
        ! 2) consistency check / adjustment of CBH in cloud-free case
        IF (nclct == 0) THEN
          !   inconsistency:  cloud base height (CBH) < CBH_CLR
          !   (WMO Man. on Codes, Vol I.1, AnnexII, Alphanum. codes, Table 1600:
          !    CBH = 2500 m may indicate zero cloud, is hence not inconsistent)
          !   (mvtsu = 6 = 'clouds not detected below follow. height' -> France)
          IF (     ((icbh(irpl) > 0) .AND. (icbh(irpl) < 2500))                &
              .OR. (      (icbh(irpl) > 2500) .AND. (icbh(irpl) < cbh_clr)     &
                    .AND. (ncsig /= 6))) THEN
            icbh (irpl) = -1
            nclct       = -1
          !  CBH not reported or zero or >= 2500 m: assume cloud-free is correct
          ELSE
            icbh (irpl) = cbh_clr
          ENDIF
        ENDIF

        ! 3) convert [m] into key code Table 1600
        !    (WMO Manual on Codes, Vol I.1, AnnexII, Alphanumeric codes):
        !    lower limits of bins: 0, 50,100,200,300,600,1000,1500, 2000,2500m)
        IF (icbh(irpl) == -1) THEN
          nhkey = nibits(nvhoc)
        ELSEIF (icbh(irpl) <=   99) THEN
          nhkey =  icbh(irpl)        /  50
        ELSEIF (icbh(irpl) <=  399) THEN
          nhkey =  icbh(irpl)        / 100 + 1
        ELSEIF (icbh(irpl) <=  999) THEN
          nhkey = (icbh(irpl) - 200) / 400 + 4
        ELSEIF (icbh(irpl) <= 2499) THEN
          nhkey =  icbh(irpl)        / 500 + 4
        ELSE               ! >= 2500 m
          nhkey = 9
        ENDIF
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvhbp, nvhoc, nhkey , 0 )
!       IF (icbh(irpl) /= -1)  osgbdy (nsgob,nbscbs) = REAL( icbh(irpl), wp )

  !       ----------------------------
  ! GCG:  low or mid-level cloud cover  ( <-- 'MNH' )
  ! ----------------------------------

        IF (     (iclclm(irpl) == imdi) .OR. (iclclm(irpl) > nibits(nxcloc))   &
            .OR. (iclclm(irpl) < 0))                                           &
          iclclm (irpl)  =  nibits(nxcloc)
        mosgbd (nsgob,nbsclg)  =  ireplace ( mosgbd(nsgob,nbsclg)              &
                                           , nxclbp, nxcloc, iclclm(irpl), 0 )
        mosgbd (nsgob,nbscwg)  =  ireplace ( mosgbd(nsgob,nbscwg)              &
                                           , nvnhbp, nxcloc, iclclm(irpl), 0 )

  ! GCG:  derive low, middle, and high cloud cover (cloud amount)
  ! -------------------------------------------------------------
  !   WMO Manual on Codes, Vol I.1, Annex II, Alphanumeric codes: 12.2.7.3:
  !       The coding of CL, CM and CH clouds (for 'MHN' + cloud types) shall be
  !       as specified in the International Cloud Atlas (WMO-No. 407), Vol I.
  !       Note:  This is according to cloud type, not definite height levels.
  !              Hence, cloud base height is given preference over vertical
  !              significance (or cloud type).
  !   In ICON, layer of low  cloud is for p > 0.8 ps,
  !            layer of high cloud is for p < 0.4 ps, mid-level cloud in between
  !   In standard atmosphere, p = 0.8 ps is at ~ 1940 m,
  !                           p = 0.4 ps is at ~ 7180 m;
  !   hence height limits of <= 1999 m and <= 6999 m are used.
  !   Keep in mind that for FM12, cloud base height = 2500 m means >= 2500 m
  !   (WMO Manual on Codes, Vol I.1, AnnexII, Alphanumeric codes: Table 1600)
  !   Note: Only 1 cloud layer (low, mid, or high) derived from GCG can be > 0.

        nclqf    = 0
        clcd     = -1._wp
        !   re-set missing value of iclclm to -1 for convenience,
        !     incl. original value 10 = sky partially obscured (e.g. by fog)
        IF (     (iclclm(irpl) < 0) .OR. (iclclm(irpl) > 13)                   &
            .OR. (iclclm(irpl) == 10))  iclclm (irpl) = -1
        !   inconsistency: zero total, but non-zero 'MNH' reported: discard all
        IF ((nclct == 0) .AND. (iclclm(irpl) >= 1)) THEN
          nclct        = -1
          iclclm(irpl) = -1
          icbh  (irpl) = -1
        ENDIF

        !   for scattered / few cloud, type 'real' is needed
        clcd (1) = REAL( iclclm(irpl) , wp )
        ! re-assign cloud amount (WMO Table 020011) from element 'MNH'
        !    FM 15 Metar, regulation 15.9.1.1:  = 12: broken    --> 5 - 7 octas
        !                                       = 11: scattered --> 3 - 4 octas
        !                                       = 13: few       --> 1 - 2 octas
        IF ((iclclm(irpl) >= 11) .AND. (iclclm(irpl) <= 13)) nclqf = 1
        IF (iclclm(irpl) == 12)  clcd(1) = 6.0_wp
        IF (iclclm(irpl) == 11)  clcd(1) = 3.5_wp
        IF (iclclm(irpl) == 13)  clcd(1) = 1.5_wp
        iclclm (irpl) = NINT( clcd(1) )
        clcd (2) = clcd(1)
        clcd (3) = -1._wp
        clcd (4) = REAL( nclct , wp )
        ! at this point, -1 <= iclclm, clcd(1:4) <= 9

        ! derive low / mid-level / high cloud amount where possible;
        ! total cloud is used for this only if it indicates zero cloud or fog,
        !                     or for high cloud if low + mid-level cloud is zero
        ! zero total cloud: (MNH > 0 excluded, see above --> MNH = 0 or -1)
        IF (nclct == 0) THEN
          clcd (:) = 0._wp
        ! non-zero total cloud and zero MNH --> high cloud is zero
        ! (note: this may lead to cloud layers slightly inconsistent with model
        !        equivalents because reported low / mid-level cloud is based on
        !        cloud type, not on cloud base / top height as in the model)
        ELSEIF (     (nclct >= 1) .AND. (nclct /= 9) .AND. (iclclm(irpl) == 0) &
               .AND. ((icbh(irpl) >= 2500) .OR. (icbh(irpl) == -1))) THEN
          clcd (3) = clcd(4)
        !   cloud amount = 9 : 'sky obscured by fog or other meteorol. phenom.'
        !   assume fog / low cloud with fully covered sky (8 octas) if both
        !    - sky obscured ('MNH'=9, or 'N'=9 (then 'MNH', 'NH' not reported))
        !    - fog flag is set or vertical significance = ceiling, and
        ELSEIF (      ((iclclm(irpl) == 9) .OR. (nclct == 9))                  &
                .AND. (      (kfog(irpl) >=  3)                                &
                       .OR. ((kfog(irpl) /= -1) .AND. (ncsig == 5)))) THEN
!               .AND. (lfog(irpl)) THEN  ! hitherto
          nclqf    =  1
          clcd (4) =  8._wp
          clcd (1) = -1._wp
          clcd (2) = -1._wp
          IF (     ((icbh(irpl) >= 0) .AND. (icbh(irpl) <= h_low))             &
              .OR. (kfog(irpl) >= 2)) THEN
            clcd (1) =  8._wp
          ELSEIF ((icbh(irpl) > h_low) .AND. (icbh(irpl) <= h_mid)) THEN
            clcd (1) =  0._wp
            clcd (2) =  8._wp
          ENDIF
        ! otherwise if sky obscured then do not trust CBH or vert. significance
        !   (and hope to get info from individual cloud layer group)
        ELSEIF (iclclm(irpl) == 9) THEN
          clcd (1) = -1._wp
          clcd (2) = -1._wp
        ! assign cloud layer from cloud base height (CBH)
        ELSEIF (icbh(irpl) /= -1) THEN
          !   inconsistency: ('MNH' = 0 .and. CBH < 2500 m)
          !           (note: CBH = 2500 m in Table 1600 for zero cloud is ok)
          IF ((icbh(irpl) < 2500) .AND. (iclclm(irpl) == 0)) THEN
            clcd (1) = -1._wp
            clcd (2) = -1._wp
            icbh (irpl) = -1           ! discard CBH if inconsistent with MNH ?
          ELSE
            IF (icbh(irpl) <= h_low)  clcd (2) = -1._wp  ! --> low cloud
            IF (icbh(irpl) >  h_low)  clcd (1) =  0._wp  ! --> no low cloud
            IF (icbh(irpl) >  h_mid)  clcd (2) =  0._wp  ! --> no mid-lev cloud,
            IF (icbh(irpl) >  h_mid)  clcd (3) = clcd(4) !    ... but high cloud
!           IF (zcbase(irpl) <  ccl_lim)  nclcm = nibits(nxcloc)
!           IF (zcbase(irpl) >= ccl_lim)  nclcl = 0
!           IF (zcbase(irpl) >  cch_lim)  nclcm = 0
!           IF (zcbase(irpl) >  cch_lim)  nclch = nclct
          ENDIF
        ! vertical significance ('MVTSU'; ncsig): low, middle, or high cloud
        ! -> very few reports exist with missing CBH and existing 'MHN' +'MVTSU'
        !    from which cloud layer could be derived (if trusting 'ncsig')
        ELSEIF ((ncsig >= 7) .AND. (ncsig <= 9)) THEN
          IF (ncsig == 7) clcd(2) = -1._wp       ! low  cloud (WMO Table 008002)
          IF (ncsig >= 8) clcd(1) =  0._wp       ! mid-level or high cloud
          IF (ncsig == 9) clcd(2) =  0._wp       ! high cloud
          IF (ncsig == 9) clcd(3) = clcd(4)      ! high cloud
        ! (ncsig = 62): vertical significance = value not applicable: in many
        !   countries (DE, CH, AT, ...) this denotes 'clear sky (no clouds)';
        !   however in DK, (ncsig == 62) and ('N' > 0) occurs --> do not infer
        !   (zero) low/mid/high cloud from (ncsig == 62) .and. ('N' missing)
!       ELSEIF ((ncsig == 62) .AND. (nclct == -1)) THEN
!         clcd (1:3) = 0._wp
        ! it cannot be determined whether 'MNH' is low or mid-level cloud
        ELSEIF (iclclm(irpl) >= 1) THEN
          clcd (1) = -1._wp
          clcd (2) = -1._wp
        ENDIF
        ! total cloud and derived high cloud: obscured or undefined, or missing
        IF (NINT( clcd(4) ) == 9)  clcd (4) = -1._wp
        IF (NINT( clcd(3) ) == 9)  clcd (3) = -1._wp
        IF (clcd(1) > epsy)        clcd (2) = -1._wp
        ! (hereafter -1. <= clcd(1:3) <= 8.)

        ! consistency checks:
        ! ------------------
        ! - discard cloud levels if CBH and vert. significance not consistent
        !   (emphasis on low cloud, with tolerances):
        !   (note: tolerances for low / middle cloud bounds are asymmetric since
        !          middle cloud (type) with CBH < H_LOW can reach into mid-level
        !          cloud layer, whereas low cloud (type) with CBH >= H_LOW
        !          cannot reach down into low cloud layer)
        !   (note: for high cloud (ncsig == 9), reported CBH can be 2500 m)
        IF (     ((ncsig == 7) .AND. (icbh(irpl) > h_low+ 400))                &
            .OR. ((ncsig == 8) .AND. (icbh(irpl) < h_low- 600))                &
            .OR. ((ncsig == 9) .AND. (icbh(irpl) < h_mid-2000)                 &
                               .AND. (icbh(irpl) /= 2500))) THEN
          clcd (1)    = -1._wp
          clcd (2)    = -1._wp
          clcd (3)    = -1._wp
          icbh (irpl) = -1
        ENDIF
        ! - discard cloud levels if 'MNH' and vert. significance not consistent
        IF (     ((iclclm(irpl) == 0) .AND. ((ncsig == 7) .OR. (ncsig == 8)))  &
            .OR. ((iclclm(irpl) >= 1) .AND.  (ncsig == 9))) THEN
          clcd (1)    = -1._wp
          clcd (2)    = -1._wp
          clcd (3)    = -1._wp
          iclclm (irpl) = -1
        ENDIF
        ! - total cloud cover < low or mid-level cloud cover:
        !   - allow tolerance of 1 octa and adjust total cloud cover
        IF (clcd(4) >= -epsy) THEN
          IF (clcd(4) >= clcd(1)-1.01_wp)  clcd (4)  = MAX( clcd(4), clcd(1) )
          IF (clcd(4) >= clcd(2)-1.01_wp)  clcd (4)  = MAX( clcd(4), clcd(2) )
        !   - otherwise reject cloud amounts
          IF (clcd(4) < MAX( clcd(1), clcd(2) )-epsy)  clcd (1:4) = -1._wp
        ENDIF
        !   if high cloud cover, adjust (TAC) reported CBH of 2500 m to high CBH
        !                        to make it consistent
        IF ((clcd(3) > epsy) .AND. (icbh(irpl) == 2500)) icbh(irpl) = h_mid+ 501
        ! (former creation or adjustment of total cloud cover obs removed:
        !  --> clcd(4) = MAXVAL( clcd(1:4) )  .)

!       IF (     (clcd(4) >= -epsy)        .OR. (nhkey /= nibits(nvhoc))       &
!           .OR. (ncsig /= nibits(nxsgoc))) THEN
!         nclcl        = NINT( clcd(1) )   ;   IF (nclcl == -1) nclcl = 15
!         nclcm        = NINT( clcd(2) )   ;   IF (nclcm == -1) nclcm = 15
!         nclch        = NINT( clcd(3) )   ;   IF (nclch == -1) nclch = 15
!         nclct        = NINT( clcd(4) )   ;   IF (nclct == -1) nclct = 15
!         IF (iclclm(irpl) == -1)  iclclm(irpl) = 15
!         WRITE( 0,'("ZCL1 ",A8,I4,I5,6I3,I6,I3,I2,I5,F7.0)' )                 &
!                ystidl(irpl), irpl, iclct(irpl), nclct, iclclm(irpl)          &
!                     , nclcl, nclcm, nclch, ncsig, icbh(irpl), nhkey          &
!                     , kfog(irpl), iww(irpl), osgbdy(nsgob,nbsvis)
!       ENDIF

        ! write to ODR: low, middle, high, total cloud cover;
        ! ------------  include low cloud flag in combined flag word
        IF (clcd(1) >= -epsy)  osgbdy (nsgob,nbscl) = clcd(1)
        IF (clcd(2) >= -epsy)  osgbdy (nsgob,nbscm) = clcd(2)
        IF (clcd(3) >= -epsy)  osgbdy (nsgob,nbsch) = clcd(3)
        IF (clcd(4) >= -epsy)  osgbdy (nsgob,nbsct) = clcd(4)
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nvcqbp, nvcqoc, nclqf , 0 )

  ! GCG:  ceiling   (note: only either low .xor. mid-level .xor. high cloud
  ! -------------          can be non-zero according to the above derivation)
        IF ((MAXVAL( clcd(1:3) ) > 4.1_wp) .AND. (icbh(irpl) /= -1)) THEN
          osgbdy (nsgob,nbscil) = REAL( icbh(irpl) , wp )
        !   if total cloud cover is > 0 (but <= 4 octa), CBH is requested to be
        !   reported here to avoid the situation that stations never reporting
        !   CBH can only yield ceiling = cbh_clr (which leads to biased obs)
        !   (if total cloud cover == 0, CBH is not requested to be reported)
        ELSEIF (      (clcd(4) >= -epsy) .AND. (clcd(4) <= 4.1_wp)             &
                .AND. ((clcd(4) <=  epsy) .OR. (icbh(irpl) /= -1))) THEN
          osgbdy (nsgob,nbscil) = REAL( cbh_clr   , wp )
        ENDIF

  ! GCG:  cloud base height (observed CBH used only if some cloud cover observed
  ! -----------------------                      or if total cloud is zero)
        IF (MAXVAL( clcd ) > epsy) THEN
          IF (icbh(irpl) /= -1)  osgbdy (nsgob,nbscbs) = REAL( icbh(irpl) , wp )
        ELSEIF ((NINT( clcd(4) ) == 0) .AND. (icbh(irpl) == cbh_clr)) THEN
                                 osgbdy (nsgob,nbscbs) = REAL( icbh(irpl) , wp )
        ELSE
          icbh(irpl) = -1
        ENDIF
      ENDIF
    ENDDO
  ENDIF

! -----------------------------------
! individual cloud layer group (ICLG)
! -----------------------------------
! IF ((my_cart_id == 0) .AND. (lfirst)) THEN
!   WRITE( 0, '("Cloud -                 reported     derived  "               &
!              ," derived (incl individual cloud layers)")' )
!   WRITE( 0, '("Cloud --          general cloud group   |     "               &
!              ," and   saved   as   observations")' )
!   WRITE( 0, '("Cloud ---          20003 oct 8002 [m] index   "               &
!              ," cld-cov-octa   [m]   [m]    [m]")' )
!   WRITE( 0, '("Cloud -sta ID   UTC   WW MNH vsig CBH  fog    "               &
!              ," N_L _M _H  N   CBH  CEIL    VIS")' )
! ENDIF

  IF ((     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )               &
       .OR. (kcdftyp == ncdf_metar))   .AND. (maxlev >= 1)) THEN
    ALLOCATE ( icbase (MAX( maxlev,1 )) , STAT=istat )
    ALLOCATE ( zclxc  (MAX( maxlev,1 )) , STAT=istat )
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      nlev  = iclev(irpl)
      IF (mosghd(nsgob,nhpass) == -1)  nlev = 0
      ! derive and store individual cloud layer group entries
      ! -----------------------------------------------------
      mosgbd (nsgob,nbscl1)  =  imdi
      mosgbd (nsgob,nbscl2)  =  imdi
      mosgbd (nsgob,nbscl3)  =  imdi
      mosgbd (nsgob,nbscl4)  =  imdi
      DO ilev = 1 , maxlev
        nicl  = imdi
        icbase(ilev) = nibits(nxbsoc)
        IF (ilev <= nlev) THEN
          iclxsg (ilev,irpl)  =  MIN( iclxsg(ilev,irpl) , nibits(nxsgoc) )
          iclxcl (ilev,irpl)  =  MIN( iclxcl(ilev,irpl) , nibits(nxcloc) )
          iclxct (ilev,irpl)  =  MIN( iclxct(ilev,irpl) , nibits(nxctoc) )
          IF (iclxsg(ilev,irpl) < 0)  iclxsg (ilev,irpl)  =  nibits(nxsgoc)
          IF (iclxcl(ilev,irpl) < 0)  iclxcl (ilev,irpl)  =  nibits(nxcloc)
          IF (iclxct(ilev,irpl) < 0)  iclxct (ilev,irpl)  =  nibits(nxctoc)
          IF (      (zclxbs(ilev,irpl) >= c0)                                  &
              .AND. (zclxbs(ilev,irpl) <= REAL( nibits(nxbsoc) - 2, wp)))      &
            icbase(ilev) = NINT( zclxbs(ilev,irpl) )
          nicl = ireplace ( nicl, nxsibp, nxsgoc, iclxsg(ilev,irpl), 0 )
          nicl = ireplace ( nicl, nxclbp, nxcloc, iclxcl(ilev,irpl), 0 )
          nicl = ireplace ( nicl, nxctbp, nxctoc, iclxct(ilev,irpl), 0 )
          nicl = ireplace ( nicl, nxbsbp, nxbsoc, icbase(ilev)     , 0 )
        ENDIF
        IF (ilev == 1)  mosgbd (nsgob,nbscl1)  =  nicl
        IF (ilev == 2)  mosgbd (nsgob,nbscl2)  =  nicl
        IF (ilev == 3)  mosgbd (nsgob,nbscl3)  =  nicl
        IF (ilev == 4)  mosgbd (nsgob,nbscl4)  =  nicl
      ENDDO

! derive / update low, middle, and high cloud cover, cloud base height, ceiling
! -----------------------------------------------------------------------------
! - cloud base height, ceiling: use ICLG to replace info from GCG (general cloud
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~  group) .... ======= .... where possible, because
!   - cloud base height in GCG does no need to be exact but should only indicate
!     some value within the correct layer according to WMO Table 1600 (Vol I.1)
!   - e.g. in case of several low cloud layers and cloud cover <= 4 octas for
!     the lowest layer and >= 5 octas for another low layer, the use of GCG info
!     will result in either (both usages of 'MHN' are found):
!     - too little low cloud cover if 'MHN' refers to the lowest layer only, or
!     - ceiling being erroneously assigned the base of the lowest layer if 'MNH'
!       refers to all low cloud layers
! - cloud amounts: use ICLG to cross-check and complement info from GCG
!   ~~~~~~~~~~~~~              ===========     ==========

      ! remove 'cloud' layers with layer index > 1 and zero cloud base height
      !   (reported e.g. by some Belgian stations)
      nlv = nlev
      DO ilev = nlv , 2 , -1
        IF (icbase(ilev) == 0)  nlev = MAX( nlev , ilev-1 )
      ENDDO

      ! re-assignment of cloud amount (WMO Table 020011)
      IF (nlev >= 1)  zclxc (1:nlev) = -1._wp
      DO ilev = nlev , 1 , -1
        IF ((iclxcl(ilev,irpl) >= 0) .AND. (iclxcl(ilev,irpl) <= 13))          &
          zclxc (ilev) = real( iclxcl(ilev,irpl) , wp )
        !   WMO Manual on Codes, Vol I.1, Annex II, Alphanumeric codes:
        !    FM 15 Metar, regulation 15.9.1.1:  = 12: broken    --> 5 - 7 octas
        !                                       = 11: scattered --> 3 - 4 octas
        !                                       = 13: few       --> 1 - 2 octas
        IF (iclxcl(ilev,irpl) == 12)  zclxc (ilev) = 6.0_wp
        IF (iclxcl(ilev,irpl) == 11)  zclxc (ilev) = 3.5_wp
        IF (iclxcl(ilev,irpl) == 13)  zclxc (ilev) = 1.5_wp
        !    = 10: sky partially obscured (e.g. by fog) --> missing value
        IF (iclxcl(ilev,irpl) == 10)  zclxc (ilev) = -1._wp
        !    =  9: 'sky obscured by fog or other meteorol. phenomena'
        IF (iclxcl(ilev,irpl) ==  9)  zclxc (ilev) = -1._wp
        IF ((iclxcl(ilev,irpl) == 9) .AND. (ilev == 1)) THEN
          IF (      (kfog(irpl) >=  3)                                         &
              .OR. ((kfog(irpl) /= -1) .AND. (     (iclxsg(ilev,irpl) == 5)    &
                                              .OR. (icbase(ilev) <= 99))))     &
            zclxc (ilev) = 8._wp
        ENDIF
      ENDDO
      !    hereafter, -1 <= zclxc <= 8 is used instead of 'iclxcl'
      !    (note: cross-check with vert. significance is not done for individual
      !           cloud layers as values 7 - 9 (low, middle, high) are not
      !           reported (only values 1 - 4, 21 - 24, 5 (and 0, 63) are found)

      ! derive low / mid-level / high cloud amounts, cloud base height, ceiling
      ! based on individual cloud layer group info only
      ! -----------------------------------------------------------------------
      clg_cl  = -1._wp
      clg_cm  = -1._wp
      clg_ch  = -1._wp
      clg_cml = 15._wp ! lowest mid-level cloud layer
      iclg_cb = -1
      iclg_ci = -1
      DO ilev = nlev , 1 , -1
        !   use only layers with both valid + non-zero cloud amount + valid CBH
        IF ((NINT(zclxc(ilev)) >= 1) .AND.(icbase(ilev) /= nibits(nxbsoc))) THEN
          !   derive values in top down loop
          !     cloud base height, ceiling
                                      iclg_cb = icbase(ilev)
          IF (zclxc(ilev) >= 4.1_wp)  iclg_ci = icbase(ilev)
          !     cloud amount of low / mid-level / high cloud layers
          IF     (icbase(ilev) <= h_low) THEN
            clg_cl  = MAX( clg_cl , zclxc(ilev) )
          ELSEIF (icbase(ilev) <= h_mid) THEN
            clg_cm  = MAX( clg_cm , zclxc(ilev) )
            !     lowest mid-level layer
            clg_cml =               zclxc(ilev)
          ELSE
            clg_ch  = MAX( clg_ch , zclxc(ilev) )
          ENDIF
        ENDIF
      ENDDO
      !   complement cloud-free layers, if cloud layers exist only further above
      IF ((clg_ch > epsy) .AND. (clg_cl <= epsy))  clg_cm = MAX( clg_cm , 0._wp)
      IF (MAX( clg_ch, clg_cm ) > epsy)            clg_cl = MAX( clg_cl , 0._wp)

      ! consistency check on ceiling from GCG:
      ! -------------------------------------
      !   ceiling based on 'MNH'-based cloud cover >= 5 octa (i.e. if (lexcil))
      !   assumed suspect if ICLG-derived cloud layer further above is < 5 octa
      lexcil =       (osgbdy(nsgob,nbscil) > rmdich)                           &
               .and. (osgbdy(nsgob,nbscil) < REAL( cbh_clr, wp ) - 0.5_wp)
      IF (lexcil) THEN     !   'nbscil' = icbh(irpl)
        IF (     ((icbh(irpl) <= h_mid) .AND. (clg_ch < 4.1_wp)                &
                                        .AND. (clg_ch >= -epsy))               &
            .OR. ((icbh(irpl) <= h_low) .AND. (clg_cm < 4.1_wp)                &
                                        .AND. (clg_cm >= -epsy)))              &
          osgbdy (nsgob,nbscil) = rmdi
      ENDIF

      !   for some British stations, mid-level or high cloud < low cloud cover
      !   --> false reporting of mid-level resp. high cloud cover
      IF ((MIN( clg_cl, clg_cm ) >= -epsy) .AND. (clg_cm < clg_cl-epsy))       &
        clg_cm = -1._wp
      IF ((MIN( clg_cl, clg_ch ) >= -epsy) .AND. (clg_ch < clg_cl-epsy))       &
        clg_ch = -1._wp
      IF ((MIN( clg_cm, clg_ch ) >= -epsy) .AND. (clg_ch < clg_cm-epsy))       &
        clg_ch = -1._wp

      ! consistency checks on cloud cover from GCG:
      ! ------------------------------------------
      !   discard ('MNH'-derived) cloud amounts if suspect (overestimated)
      !   (i.e. if contradicting ICLG which is given preference;
      !    after these checks, amounts from ICLG or from GCG may still be
      !    imcomplete, i.e. underestimated, since not all cloud layers are
      !    always reported resp. accounted)
      clcd (1) = -1._wp
      clcd (2) = -1._wp
      clcd (3) = -1._wp
      IF (osgbdy (nsgob,nbscl) > rmdich)  clcd (1) = osgbdy(nsgob,nbscl)
      IF (osgbdy (nsgob,nbscm) > rmdich)  clcd (2) = osgbdy(nsgob,nbscm)
      IF (osgbdy (nsgob,nbsch) > rmdich)  clcd (3) = osgbdy(nsgob,nbsch)

      !   high cloud:
      !   - 'MNH' sometimes reports zero by wrongly assigning (mid-level) cloud
      !     (mainly above 6000m) -> high cloud derived from total cloud suspect
      !     - unless its cloud amount is clcd(3) > clg_cm: in this case assume
      !       that 'clcd(3)' was based on a different (higher) cloud layer); or
      !     - if mid-level (or low) cloud with base below ~ 5000 m exists
      !   - high cloud derived from total cloud is suspect if low cloud exists
      IF (      (clcd(3) > epsy)                                               &
          .AND. (     (clcd(3) <= clg_cm+epsy)                                 &
                 .OR. ((iclg_cb < 5000) .AND. (iclg_cb /= -1))                 &
                 .OR. (clg_cl > epsy)))                                        &
        clcd (3) = -1._wp
      !   mid-level cloud:
      !   - high cloud being misplaced as mid-level cloud in 'MNH' is assumed
      !     if individual cloud layer group contains only high cloud;
      !   - low cloud with base height <= 1999 in individual cloud layer group
      !     is occasionally misplaced as mid-level cloud in 'MNH' (this must be
      !     the lowest cloud layer - if an even lower cloud layer existed, it
      !     would have been interpreted as low cloud in 'MNH')
      IF (      (clcd(2) > epsy)                                               &
          .AND. (     ((MAX( clg_cl, clg_cm ) <= epsy) .AND. (clg_ch > epsy))  &
                 .OR. (clg_cl > epsy)))                                        &
        clcd (2) = -1._wp
      !   low cloud: suspect
      !   - if individual cloud group reports cloud layers, but no low cloud; or
      !   - if 'MNH' (>)== clg_cml (= cloud amount of lowest mid-level layer)
      !     ('MNH'-derived low cloud occasionally includes 1 cloud layer with
      !      2000 m <= CBH < 2500 m (above other low cloud layers))
      IF (      (clcd(1) > epsy)                                               &
          .AND. (     ((MAX( clg_ch, clg_cm ) > epsy) .AND. (clg_cl <= epsy))  &
                 .OR. (clcd(1) >= clg_cml-epsy)))                              &
        clcd (1) = -1._wp
      IF (clcd(1) < -epsy)  osgbdy (nsgob,nbscl) = rmdi
      IF (clcd(2) < -epsy)  osgbdy (nsgob,nbscm) = rmdi
      IF (clcd(3) < -epsy)  osgbdy (nsgob,nbsch) = rmdi

      ! assign values to ODR
      ! --------------------
      ! cloud amounts: take max of general and individual cloud group info
      ! ~~~~~~~~~~~~~  (info from any cloud group may still be incomplete and
      !                 not include all layers e.g. of low cloud)
      !                (zero amount is also allowed; missing is < 0)
      IF (clg_ch  >= 0)  osgbdy (nsgob,nbsch ) = MAX( clg_ch , clcd(3) )
      IF (clg_cm  >= 0)  osgbdy (nsgob,nbscm ) = MAX( clg_cm , clcd(2) )
      IF (clg_cl  >= 0)  osgbdy (nsgob,nbscl ) = MAX( clg_cl , clcd(1) )
      ! cloud base height CBH, ceiling: prefer individual cloud group info
      IF (iclg_cb >= 0)  osgbdy (nsgob,nbscbs) = REAL( iclg_cb , wp )
      IF (iclg_ci >= 0)  osgbdy (nsgob,nbscil) = REAL( iclg_ci , wp )
      !    (note: iclg_cb == -1 if CBH missing or cloud cover missing or zero,
      !                         or cloud cover == 9 and not fog diagnosed)

      !   control output
!     IF (     (iww(irpl) == imdi) .OR. (iww(irpl) >= nibits(nvw0oc))          &
!         .OR. (iww(irpl) < 0)) iww(irpl) = -1
!     ivis = -1
!     IF ((zvis(irpl) >= -epsy) .AND. (zvis(irpl) < 99999.9_wp))               &
!       ivis = NINT( zvis(irpl) )
!     DO ilev = 1 , nlev
!       IF (icbase(ilev) == nibits(nxbsoc))  icbase(ilev) = -1
!       WRITE( 0,'("ZCL9 ",A8, I5 ,2I3,3F4.0,I6,I3," GCG",I4,I3,I5,I6,I3       &
!                         ,2I6," ODR ",4I3,2I6,I7)' )                          &
!              ystidl(irpl), mosghd(nsgob,nhhrmn), ilev, iclxsg(ilev,irpl)     &
!                          , clg_cl, clg_cm, clg_ch                            &
!                          , icbase(ilev), iclxcl(ilev,irpl)                   &
!                          , iww(irpl), kfog(irpl), iclclm(irpl), icbh(irpl)   &
!                          , iclsig(irpl), iclg_cb, iclg_ci                    &
!            , MOD( NINT( MAX( osgbdy(nsgob,nbscl) , -1._wp ) ) +100 , 100 )   &
!            , MOD( NINT( MAX( osgbdy(nsgob,nbscm) , -1._wp ) ) +100 , 100 )   &
!            , MOD( NINT( MAX( osgbdy(nsgob,nbsch) , -1._wp ) ) +100 , 100 )   &
!            , MOD( NINT( MAX( osgbdy(nsgob,nbsct) , -1._wp ) ) +100 , 100 )   &
!                 , NINT( MAX( osgbdy(nsgob,nbscbs), -1._wp ) )                &
!                 , NINT( MAX( osgbdy(nsgob,nbscil), -1._wp ) ), ivis
!     ENDDO
!     WRITE( 0,'("Cloud ",A8,I2,":"I2.2,I5,I4,I3,I6,I3," saved:",4I3,2I6,I7)') &
!            ystidl(irpl), mosghd(nsgob,nhhrmn) /100,                          &
!                     MOD( mosghd(nsgob,nhhrmn), 100 ) , iww(irpl)             &
!          , MOD( iclclm(irpl) +100 , 100 ), iclsig(irpl), icbh(irpl)          &
!          ,                                               kfog(irpl)          &
!          , MOD( NINT( MAX( osgbdy(nsgob,nbscl) , -1._wp ) ) +100 , 100 )     &
!          , MOD( NINT( MAX( osgbdy(nsgob,nbscm) , -1._wp ) ) +100 , 100 )     &
!          , MOD( NINT( MAX( osgbdy(nsgob,nbsch) , -1._wp ) ) +100 , 100 )     &
!          , MOD( NINT( MAX( osgbdy(nsgob,nbsct) , -1._wp ) ) +100 , 100 )     &
!               , NINT( MAX( osgbdy(nsgob,nbscbs), -1._wp ) )                  &
!               , NINT( MAX( osgbdy(nsgob,nbscil), -1._wp ) ), ivis
    ENDDO
    DEALLOCATE ( icbase  , STAT=istat )
    DEALLOCATE ( zclxc   , STAT=istat )
    DEALLOCATE ( iclxsg  , STAT=istat )
    DEALLOCATE ( iclxcl  , STAT=istat )
    DEALLOCATE ( iclxct  , STAT=istat )
    DEALLOCATE ( zclxbs  , STAT=istat )
  ENDIF

! ---------------------------
! max. / min. 2-m temperature
! ---------------------------

  IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1) THEN
        osgbdy (nsgob,nbstx ) =  ztmax(irpl)
        osgbdy (nsgob,nbstn ) =  ztmin(irpl)
        mosgbd (nsgob,nbsttr) =   insert( 0 , ABS( itmxt(irpl) ) , ntxbp )     &
                                + insert( 0 , ABS( itmnt(irpl) ) , ntnbp )
      ENDIF
    ENDDO
  ENDIF

! ---------
! radiation (global, diffuse, long-wave)
! --------------------------------------

  IF (kcdftyp == ncdf_synop) THEN
    DO irpl = 1 , nrepsg
      nsgob = ntotsg + irpl
      IF (mosghd(nsgob,nhpass) /= -1) THEN
        osgbdy (nsgob,nbsrad) = zradgl(irpl)
        osgbdy (nsgob,nbsrdd) = zraddf(irpl)
        osgbdy (nsgob,nbsrdl) = zradlw(irpl)

! ----------------
! total snow depth
! ----------------

        osgbdy (nsgob,nbshsw) = zhsnow(irpl)

! ---------------
! state of ground
! ---------------

        IF ((iee(irpl) < 0) .OR. (iee(irpl) > nibits(nveeoc)))                 &
          iee (irpl) = nibits(nveeoc)
        mosgbd (nsgob,nbswwe)  =  ireplace ( mosgbd(nsgob,nbswwe)              &
                                           , nveebp, nveeoc, iee(irpl), 0 )
      ENDIF

!     PRINT '("OSG 12-20 ",A,8F8.2)',ystidl(irpl),(osgbdy(nsgob,ncc),ncc=12,19)
!     PRINT '("OSG 20-26 ",A,7F8.2)',ystidl(irpl),(osgbdy(nsgob,ncc),ncc=20,26)
!     PRINT '("MOS  6- 9 ",A,4I12 )',ystidl(irpl),(mosgbd(nsgob,ncc),ncc= 6, 9)

    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Section 6: Get observation errors
!-------------------------------------------------------------------------------

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl
    IF ((mosghd(nsgob,nhpass) /= -1) .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
      zlop = LOG( osgbdy(nsgob,nbsp) )

      CALL obs_find_level ( nerlev, rolnlv, zlop , ilverr, fiperr )
!     ===================
    ENDIF

! height obs error (only for surface obs)
    IF (      (     ((.NOT. lsurfob) .AND. (osgbdy(nsgob,nbsz) > rmdich))      &
               .OR. ((      lsurfob) .AND. (osgbdy(nsgob,nbsp) > rmdich)))     &
        .AND. (mosghd(nsgob,nhpass) /= -1)) THEN
      IF (lsurfob)  zz2p = osgbdy(nsgob,nbsp) * r_g / (r_d * tmelt)
      IF     (kcdftyp == ncdf_buoy ) THEN
        osgbdy(nsgob,nbszer) =   osgbdy(nsgob,nbszer) + oezdrib *zz2p
      ELSEIF (kcdftyp == ncdf_ship ) THEN
        osgbdy(nsgob,nbszer) =   osgbdy(nsgob,nbszer) + oezship *zz2p
      ELSEIF (      ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_metar))     &
              .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
        osgbdy(nsgob,nbszer) =   osgbdy(nsgob,nbszer)                          &
                               + (      fiperr * oezsynp(ilverr)               &
                                  + (c1-fiperr)* oezsynp(ilverr+1)) *zz2p
      ELSE
        IF (.NOT. lsurfob)  zz2p = c1
        osgbdy(nsgob,nbszer) = oezsynp(nerlev) *zz2p
        mosgbd(nsgob,nbserr) = IBCLR ( mosgbd(nsgob,nbserr), nvrz )
      ENDIF
    ENDIF

! temperature obs error (only for aircrafts; 2-m T error is already specified)
    IF ((osgbdy(nsgob,nbst ) >= rmdich) .AND. (mosghd(nsgob,nhpass) /= -1)) THEN
      IF (      (mosghd(nsgob,nhobtp) == nairep)                               &
          .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
        osgbdy(nsgob,nbster) =                        fiperr * oetairp(ilverr) &
                                                + (c1-fiperr)* oetairp(ilverr+1)
!     (also consider entry 'MPOTO' --> zttq for ACARS)
        IF (lacars) THEN
          IF ((zttq(irpl) > 0.2_wp) .AND. (zttq(irpl) < 10._wp))               &
            osgbdy(nsgob,nbster) = MAX( osgbdy(nsgob,nbster) , zttq(irpl) )
        ENDIF
      ENDIF
    ENDIF

! rel. humidity obs error (only for aircrafts; 2-m RH error already specified)
    IF ((osgbdy(nsgob,nbsrh) >= rmdich) .AND. (mosghd(nsgob,nhpass) /= -1)) THEN
      IF (mosghd(nsgob,nhobtp) == nairep) THEN
                                           osgbdy(nsgob,nbsqer) = rherr1 * c100r
        IF (osgbdy(nsgob,nbsrh) < rrhlim)  osgbdy(nsgob,nbsqer) = rherr2 * c100r
!   increased error if very cold or if temperature not defined
        IF (osgbdy(nsgob,nbst ) < rttlim)  osgbdy(nsgob,nbsqer) = rherr3 * c100r
      ENDIF
    ENDIF

! horizontal wind obs error
    IF ((osgbdy(nsgob,nbsu ) >= rmdich) .AND. (mosghd(nsgob,nhpass) /= -1)) THEN
      IF     (kcdftyp == ncdf_buoy ) THEN
        osgbdy(nsgob,nbsuer) = oevdrib
      ELSEIF (lscatto) THEN
        osgbdy(nsgob,nbsuer) = oevscat
      ELSEIF (      (mosghd(nsgob,nhobtp) == nsynop)                           &
              .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
        osgbdy(nsgob,nbsuer) =                        fiperr * oevsynp(ilverr) &
                                                + (c1-fiperr)* oevsynp(ilverr+1)
      ELSEIF (      (mosghd(nsgob,nhobtp) == nairep)                           &
              .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
        osgbdy(nsgob,nbsuer) =                        fiperr * oeairep(ilverr) &
                                                + (c1-fiperr)* oeairep(ilverr+1)
      ELSEIF (      (mosghd(nsgob,nhobtp) == nsatob)                           &
              .AND. (osgbdy(nsgob,nbsp) > rmdich)) THEN
        osgbdy(nsgob,nbsuer) =                        fiperr * oesatob(ilverr) &
                                                + (c1-fiperr)* oesatob(ilverr+1)
      ELSEIF (      (mosghd(nsgob,nhobtp) == nsynop)                           &
              .AND. (osghed(nsgob,nhalt) > -400._wp)                           &
              .AND. (osghed(nsgob,nhalt) < 1000._wp)) THEN
        osgbdy(nsgob,nbsuer) = oevsynp(1)
      ELSEIF        (mosghd(nsgob,nhobtp) == nsynop) THEN
        osgbdy(nsgob,nbsuer) = oevsynp(4)
      ELSE
        osgbdy(nsgob,nbsuer) = oevsynp(nerlev)
        mosgbd(nsgob,nbserr) = IBCLR ( mosgbd(nsgob,nbserr), nvru )
      ENDIF
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Section 7: Check existence of data in report
!-------------------------------------------------------------------------------

  DO irpl = 1 , nrepsg
    nsgob = ntotsg + irpl

    IF (mosghd(nsgob,nhpass) /= -1) THEN
      icma   = icmaa (irpl)
      iactr  = 0
      IF (mosghd(nsgob,nhpass) == 0)  iactr  = 1

      iactx = -1
      ipasx = -1
      IF (BTEST( mosgbd(nsgob,nbserr),nvru )) iactx (1) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrz )) iactx (2) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrt )) iactx (3) = 1
      IF (BTEST( mosgbd(nsgob,nbserr),nvrq )) iactx (4) = 1
      IF (osgbdy(nsgob,nbsu  ) > rmdich) ipasx (1) = 1
      IF ((.NOT. lsurfob) .AND. (osgbdy(nsgob,nbsz) > rmdich)) ipasx (2) = 1
      IF ((      lsurfob) .AND. (osgbdy(nsgob,nbsp) > rmdich)) ipasx (2) = 1
      IF (osgbdy(nsgob,nbst  ) > rmdich) ipasx (3) = 1
      IF (osgbdy(nsgob,nbsrh ) > rmdich) ipasx (4) = 1
      iactx (5)             = MAX( iactx(1), iactx(2), iactx(3), iactx(4) )
      ipasx (5)             = MAX( ipasx(1), iactx(2), ipasx(3), ipasx(4) )
      IF (iactx(5) == -1) THEN
        IF (     (osgbdy(nsgob,nbscbs) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbscil) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbscl ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbscm ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsch ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsct ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsvis) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsrr1) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsrr3) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsrr6) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsr12) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsr24) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsfg1) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsfg3) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbsfg6) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbstn ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbstx ) > rmdich)                               &
            .OR. (osgbdy(nsgob,nbshsw) > rmdich)                               &
            .OR. (mosgbd(nsgob,nbscwg) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbswwe) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbsclg) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbscl1) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbscl2) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbscl3) /= imdi)                                &
            .OR. (mosgbd(nsgob,nbscl4) /= imdi))  iactx(5) = 2
      ENDIF
      nzaexi                = MAX( MIN( 0 , -ipasx(5) ) , iactx(5) )
!     PRINT   '("precip-1h ",A,I6,F4.1,2F5.1,4I4,F6.1,2X,5I2,2X,5I2,3I5)'      &
!           , yosghd(nsgob), mosghd(nsgob,nhkz)                                &
!           , osghed(nsgob,nhtime), osghed(nsgob,nhilon), osghed(nsgob,nhjlat) &
!           , mosghd(nsgob,nhobtp), mosghd(nsgob,nhcode), mosghd(nsgob,nhpass) &
!           , nzaexi, osgbdy(nsgob,nbsrr1), iactx, ipasx                       &
!           , mosghd(nsgob,nhitot), mosghd(nsgob,nhjtot), my_cart_id
! (nzaexi ==  1): active data present
! (nzaexi ==  2): additional data present (for MEC/verification)
! (nzaexi == -1): only passive data present
! (nzaexi ==  0): no data at all in total report
      IF (nzaexi <= 0)                                                         &
        neventr (nenoda,icma) = neventr(nenoda,icma) + iactr
      IF ((nzaexi == -1) .AND. (iactr == 1) .AND. (lverpas)) THEN
        mosghd (nsgob,nhpass) =  2
        mosghd (nsgob,nhflag) =  IBSET ( mosghd(nsgob,nhflag), FL_NO_OBS )
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1

! flag reports which are to be discarded completely
      ELSEIF ((nzaexi == 0) .OR. ((nzaexi == -1) .AND. (.NOT. lverpas))) THEN
        mosghd (nsgob,nhpass) = -1
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
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt3
          DO icl = 1 , istrej
            outbuf(nacout+2+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
    ENDIF
  ENDDO

! this is done after obs_cdf_redundancy
! TODO : delete reports with (moxxhd(nxxob,nhpass) == -1)
! CALL obs_del_old_reports
! ========================

! clean up: de-allocate standardised arrays
! -----------------------------------------

  IF (nrepsg >= 1) THEN
    DEALLOCATE ( lpract  , STAT=istat )
    IF (lsurfob) THEN
      DEALLOCATE ( npcode  , STAT=istat )
      DEALLOCATE ( landsy  , STAT=istat )
      DEALLOCATE ( kfog    , STAT=istat )
      DEALLOCATE ( icbh    , STAT=istat )
    ENDIF

    IF ((lexitd) .OR. (lexirh) .OR. (lexiqx)) THEN
      DEALLOCATE ( zrhw   , STAT=istat )
      DEALLOCATE ( ztd    , STAT=istat )
      DEALLOCATE ( zqx    , STAT=istat )
      DEALLOCATE ( zrhc   , STAT=istat )
    ENDIF

    DEALLOCATE ( zpp    , STAT=istat )
    DEALLOCATE ( zdd    , STAT=istat )
    DEALLOCATE ( zff    , STAT=istat )
    IF ((kcdftyp /= ncdf_satob) .AND. (.NOT. lscatto))                         &
      DEALLOCATE ( ztt    , STAT=istat )

    IF (lupair) THEN
      DEALLOCATE ( zzz     , STAT=istat )
      IF (kcdftyp /= ncdf_satob)  DEALLOCATE ( iroll   , STAT=istat )
      IF (kcdftyp == ncdf_satob)  DEALLOCATE ( iqual   , STAT=istat )
      IF (lacars)                 DEALLOCATE ( zttq    , STAT=istat )
      IF ((kcdftyp == ncdf_modes) .OR. (lacars))                               &
        DEALLOCATE ( iqxq    , STAT=istat )
      IF ((kcdftyp == ncdf_amdar) .OR. (lacars)) THEN
        DEALLOCATE ( iturb   , STAT=istat )
        DEALLOCATE ( zvgust  , STAT=istat )
      ENDIF
    ENDIF

    IF (lsurfob) THEN
      DEALLOCATE ( zpmsl  , STAT=istat )
!     DEALLOCATE ( ztd    , STAT=istat )
      DEALLOCATE ( zgust  , STAT=istat )
      DEALLOCATE ( zzt    , STAT=istat )
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy ) .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        DEALLOCATE ( ivdt   , STAT=istat )
        DEALLOCATE ( igudt  , STAT=istat )
        DEALLOCATE ( ivtisi , STAT=istat )
        DEALLOCATE ( irrdt  , STAT=istat )
        DEALLOCATE ( zrr    , STAT=istat )
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_buoy )) THEN
        DEALLOCATE ( zdpdt  , STAT=istat )
!       DEALLOCATE ( zrhw   , STAT=istat )
      ENDIF
      IF (     (kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )            &
          .OR. (kcdftyp == ncdf_nat_dwd)) THEN
        DEALLOCATE ( igudt2 , STAT=istat )
        DEALLOCATE ( irrdt2 , STAT=istat )
        DEALLOCATE ( zgust2 , STAT=istat )
        DEALLOCATE ( zrr2   , STAT=istat )
      ENDIF
      IF ((kcdftyp == ncdf_synop) .OR. (kcdftyp == ncdf_ship )) THEN
        DEALLOCATE ( iclev  , STAT=istat )
        DEALLOCATE ( iclct  , STAT=istat )
        DEALLOCATE ( iclsig , STAT=istat )
        DEALLOCATE ( iclclm , STAT=istat )
        DEALLOCATE ( iccl   , STAT=istat )
        DEALLOCATE ( iccm   , STAT=istat )
        DEALLOCATE ( icch   , STAT=istat )
        DEALLOCATE ( iww    , STAT=istat )
        DEALLOCATE ( iwdt   , STAT=istat )
        DEALLOCATE ( iw1    , STAT=istat )
        DEALLOCATE ( iw2    , STAT=istat )
        DEALLOCATE ( itmxt  , STAT=istat )
        DEALLOCATE ( itmnt  , STAT=istat )
        DEALLOCATE ( zvis   , STAT=istat )
        DEALLOCATE ( zrr24  , STAT=istat )
        DEALLOCATE ( zcbase , STAT=istat )
        DEALLOCATE ( ztmax  , STAT=istat )
        DEALLOCATE ( ztmin  , STAT=istat )
      ENDIF
      IF     (kcdftyp == ncdf_synop) THEN
        DEALLOCATE ( iee    , STAT=istat )
        DEALLOCATE ( ifi    , STAT=istat )
        DEALLOCATE ( zradgl , STAT=istat )
        DEALLOCATE ( zraddf , STAT=istat )
        DEALLOCATE ( zradlw , STAT=istat )
        DEALLOCATE ( zhsnow , STAT=istat )
        DEALLOCATE ( zpfi   , STAT=istat )
      ELSEIF (kcdftyp == ncdf_metar) THEN
        DEALLOCATE ( iclev  , STAT=istat )
        DEALLOCATE ( zvis   , STAT=istat )
      ENDIF
    ENDIF
  ENDIF
  lfirst = .FALSE.

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_store_singlev
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_store_singlev


!===============================================================================

ELEMENTAL INTEGER FUNCTION ireplace   ( invar, ipos, iboc, irepl, ipsr )
  !---------------------------------------------------------------------
  INTEGER         , INTENT (IN)  ::  invar, ipos, iboc, irepl, ipsr
  !-----------------------------------------------------------------------------
  ! replaces 'iboc' bits starting at bit position 'ipos' of integer word 'invar'
  ! by the 'iboc' bits starting at bit position 'ipsr' from integer word 'irepl'
  !-----------------------------------------------------------------------------
  !
  ireplace = IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )            &
                , ISHFT( IAND( ISHFT( irepl,-ipsr ), nibits(iboc) ), ipos ) )
  !
END FUNCTION ireplace

!-------------------------------------------------------------------------------

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

END MODULE src_obs_cdfin_sing
