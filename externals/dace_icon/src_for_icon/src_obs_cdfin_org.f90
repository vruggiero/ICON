!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------
#if __NEC_VERSION__ >= 30005    /* options directive supported since 3.0.5 */
!NEC$ options "-fno-move-loop-invariants"
#endif

MODULE src_obs_cdfin_org

!-------------------------------------------------------------------------------
! Description:
!   This module organizes the reading, pre-selection, distribution, and storage
!   in ODR of observation reports to be read from NetCDF observation input
!   files.
!   Special tasks are:
!    - reading reports from NetCDF files, selecting, assigning them to model
!      grid points and distributing them to the nodes (precessors) according
!      to the sub-domain which they lie on
!    - pre-processing, including some gross error checking, blacklist checking,
!      assigning observation errors, bias correction, redundancy checking,
!      multi-level gross error checking (wind shear, stability), etc.
!    - storing the reports in the ODR internal arrays, clean up ODR arrays from
!      reports which are not used any more
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files.
!         It is shared between the COSMO model and the DACE program package !
!
!   Known bug: Redundancy check is not guaranteed to be independent from domain
!   =========  decomposition unless option 'lredn_repro' is set which allows for
!              redundancy only between reports assigned to the same grid point.
!
! Method:
!   This module contains the following public module procedures:
!    - obs_cdf_read_org       : called by obs_org_cdf_proc, read_cosmo_obs,
!                               main routine calling obs_cdf_read_temp_pilot etc
!    - obs_cdf_interface      : called by obs_org_cdf_proc, read_cosmo_obs
!    - obs_cdf_proc_init      : called by obs_org_cdf_proc, read_cosmo_obs
!    - obs_cdfin_open_files   : called by obs_org_cdf_proc, read_cosmo_obs
!    - obs_fof_check_files    : called by obs_org_cdf_proc, read_cosmo_obs
!    - obs_cdf_tower_copy     : called by read_cosmo_obs in DACE
!    - obs_cdf_bcorr_rass     : called by read_cosmo_obs in DACE
!    - obs_cdf_mult_qualicheck: called by obs_org_cdf_proc, read_cosmo_obs
!    - obs_cdf_raso_rh_bias   : called by obs_org_cdf_proc, organize_obs_proc
!
!   This module contains the following private module procedures:
!    - obs_cdf_redundancy     : called by obs_cdf_read_org
!    - obs_cdf_del_old_rep    : called by obs_cdf_read_org
!    - obs_cdf_supob_levs     : called by obs_cdf_proc_init
!
!   This module also contains elemental functions, formerly statement functions:
!   - rmod         : MOD function for positive REALS
!   - fpvsw        : Magnus formula for water: saturation vapour pressure from T
!   - fpvsi        : Magnus formula for ice  : saturation vapour pressure from T
!   - ireplace     : replaces a bit pattern in an integer by another bit pattern
!   - ibflip       : reverts (flips) a given bit
!   - ileap        : detects leap year (=0: no leap year, =1: leap year)
!
!   It uses from:
!    - src_obs_cdfin_mult:    - obs_cdf_read_temp_pilot
!                             - obs_cdf_read_profiler
!    - src_obs_cdfin_sing:    - obs_cdf_read_surface
!                             - obs_cdf_read_aircraft
!    - src_obs_cdfin_gps :    - obs_cdf_read_gps
!    - src_obs_fdbk_in   :    - obs_fof_read
!    - src_obs_cdfin_util:    - get_global_surface
!                             - obs_solar_zenith_angle
!                             - obs_rm_duplicate_levels
!    - parallel_utilities:    - global_values
!                             - distribute_values
!    - environment:           - model_abort
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
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012-01-31 Christoph Schraff
!  Initial release, extracted from modules 'src_obs_proc_cdf' and
!    'src_obs_processing' and adapted (e.g. modified routine interfaces and
!    diagnostic arrays, 'obs_pointrs' --> 'i_cma', modified ODR (observation
!    status word introduced)).
!  - NetCDF observation file input made more flexible: several files for the
!    same observation type may be read.
!  - After the pure redundancy check, the redundancy check loop is re-run with
!    different settings to put TEMP or PILOT parts A, B, C, D together.
!  - No redundancy between a wind profiler and a RASS report allowed.
!  - No redundancy between reports from two different ships allowed.
!  - Correction of a bug in the redundancy check which could lead to vertically
!    collocated levels within a (composite TEMP radiosonde) profile.
!  - Computation of solar zenith angle moved out of routine for bias correction
!    of Vaisala RS92 radiosonde humidity 'obs_cdf_raso_rh_bias' into routine
!    'obs_solar_zenith_angle' (module 'src_obs_cdfin_util.f90'),
!    and solar zenith angle is written to the ODR now.
!  - New routine for bias correction of Vaisala RS92 radiosonde humidity.
!  - Karolin Eichler: for any GPS station processed by several centres,
!    preference is given in the redundancy check according to the order of
!    centres in array 'i_gpscen'.
! V4_23        2012/05/10 Oliver Fuhrer
!  Bug Fix in sending character variable yncannex according to usage of
!   distribute_values
! V4_27        2013/03/19 Christoph Schraff, Ulrich Schaettler
!  Character length of (y)ydate_ref increased from 10 to 14.
!  Call to global_values only for num_compute > 1 (around line 960) (US)
! V4_28        2013/07/12 Christoph Schraff
!  Optional capability prepared for reading obs from feedobs files;   interface
!  in obs_cdf_interface extended by 'irun_osse', 'losse_fg', 'fperturb','iseed'.
!  Bug fix when updating main flag word in lapse rate check.
!  Statement functions replaced by elemental or intrinsic functions.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  If (losse_fg) then, in order to obtain the same set of obs irrespective of
!  applying random perturbations, the lapse rate check is replaced by a lapse
!  rate adjustment of observed values, and discarding a report in the redundancy
!  check depends only the levels and observed variables rather than observed
!  values being a subset of another report. (CS)
!  Extensions for reading / processing Mode-S aircraft observations, from NetCD
!  files with 2 different templates: (i) converted from BUFR as received from
!  Mode-S processing centre KNMI, (ii) converted into DWD-ACARS template. (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_1a        2015-03-30 Christoph Schraff
!  Enlarged horizontal and vertical check limits for redundancy check of ship
!  and buoy to account for differences in BUFR and TAC (alphanumeric) reports.
! V5_3         2015-10-09 Christoph Schraff
!  - For verification (model equivalent calculator) purposes, additional
!    quantities of single-level ODR written to active report in redundancy check
!  - Enlarged horizontal, vertical and temporal check limits in redundancy check
!    of SYNOP to account for differences in BUFR and TAC (alphanumeric) reports.
!  - Use of obs_fof_read conditioned on ifdef __COSMO__ to make it compatible
!    with LETKF / MEC code.
!  - For DWD: prefer SYNOP with small KZ (TAC) over SYNOP with large KZ (BUFR)
!    in redundancy check.
! V5_4         2016-03-10 Christoph Schraff
!  - Dimension of 'neventr' and 'neventd' reduced from 3 to 2.
!  - Variables related to AOF interface removed.
! V5_4d        2016-12-12 Christoph Schraff
!  Interface for modules reading observations extended by 'zlimv10'.
! V5_4f        2017-09-01 Christoph Schraff
!  - Minor change to the redundancy check to help to select the TEMP parts B, A
!    (in particular B, and not C, D) as active report, in order to get the obs
!    and decoding times from Part B in the header of the composite radiosonde
!    report.
!  - Update from DACE V1_51: In redundancy check, horiz. distances between obs
!    made compatible with DACE source + applicable also to ICON, GME, IFS grids.
! V5_4h        2017-12-15 Christoph Schraff
!  Call to 'obs_cdf_read_amv' added for reading SATOB / AMV reports.
! V5_5b        2018-10-29 Christoph Schraff
!  Fix to ensure that nbsu and nbsv never have missing values
! V5_6         2019-02-27 Christoph Schraff
!  Option 'lredn_repro' introduced to ensure reproducibility of redundancy check
!  irrespective of domain decomposition by allowing for redundancy only between
!  reports assigned to the same grid point.
! V5_6b        2019-10-16 Christoph Schraff
!  - Obs input file names 'cdfin_*.x' are now inquired and read also for x = 1.
!  - For processing steps (e.g. redundancy checking) affecting the state (active,
!    passive...) of a report or observation, 'observed pressure' (nbtp, nbsp) is
!    replaced by 'pressure level independent from model state' (nbtplv, nbsplv)
!    in order to render the number of reports and observations independent
!    from the model state and thus from the member of an ensemble (for safety,
!    this is required for 'online LETKF').
!    Because of this change, redundancy checking between 2 (upper-air) reports
!    is excluded where exactly 1 of the 2 reports has height as level info
!    (and uses standard atmosphere to derive a model state independent, but
!    highly inaccurate pressure).
!  - Bug fix in the redundancy checking (not autochecking) to avoid 2 levels
!    with the same pressure: the pressure and height of a level in the active
!    report is replaced by that of a mandatory level from the rejected report
!    only if the newly assigned level will be below (i.e. at higher pressure
!    than) the subsequent level of the active report.
!  - Redundancy check adjusted to avoid supplementing a superobbed profile by
!    data from the original report (with state 'merged').
!  - Subroutine 'obs_cdf_supob_levs' added for transforming namelist input
!    'av_levs', 'av_incr' into 'nref' superobbing reference levels 'pref'.
!  - 'maxav', 'av_reso' added to interface of routine 'obs_cdf_interface' for
!    superobbing of high-resolution (radiosonde and other) profile reports
!    (based on 'r_av_levs', 'r_av_incr'), 'imaxtvl' removed.
!  - Introduction of reading tower profile reports, high-resolution BUFR
!    radiosonde profile reports and descending TEMP reports.
!  - Minor bug fix in the use of flags to decide on replacing single data
!    from the active (multi-level) report by data from the passive report.
!  - Revised interface of 'obs_cdf_proc_init' ('nuspecif', 'yuspecif' added).
!  - Set subroutines used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! V5_7a        2020-05-11 Christoph Schraff
!  - Processing of ground-based wind lidar added.
!  - Exclude merging (redundancy) of a radiosonde descent report with a report
!    of any other code type (e.g. with a radiosonde ascent report).
!  - The criterion for 'lmerge' (idea: merging TAC TEMP parts A-D) is tightened
!    in order to allow more easily for ordinary redundany checking between
!    TAC and BUFR reports.
!  - New criteria of preferences in redundancy check for multi-level reports
!    introduced:
!    - prefer high-resolution superobbed report over other report,
!    - prefer larger vertical range of superobbed levels,
!    - prefer BUFR TEMP report over TAC TEMP report.
!  - The criterion for 'lsub' (i.e. for discarding a multi-level report because
!    it is subset of another report) is relaxed to still apply
!    - if mandatory levels are added to a superobbed report;
!    - if only the top level of a superobbed redundant profile is absent in
!      (but within the range of) the active profile;
!    - even though observed values differ by some tolerances (1.1 m/s, 0.4 K,
!      RH=0.5) (and the pressure levels differ by <= rprlim.
!  - Bug fix at redundancy check with national DWD reports.
!  - Set losse_fg = .false. if there are no 'fof' input files.
! V5_7b        2020-06-19 Ulrich Schaettler
!  Added a NEC NOMOVE directive
! @VERSION@  @DATE@     Christoph Schraff
!  - Some code clean-up.
!  - Bug fix in redundancy check, avoiding missing obs error for valid RH obs.
! 2021-10-27 Christoph Schraff
!  Opening of files only for obs types not processed by dace obs operators.
! 2021-11-03 Christoph Schraff
!  - Minor bug fix in redundancy check to give preference to complete TEMPHIRS
!    reports over early reports with top level at 100 hPa.
! 2021-11-03 CS: Redundancy check: minor bug fix to give preference to complete
!                TEMPHIRS reports over early reports with top level at 100 hPa.
! 2022-08-01 CS: Replacement of cloud ceiling added in redundancy check.
! 2023-04-17 CS: Special redundancy limits for synop vs. radiosonde surface rep.
! 2023-04-18 CS: 'rhtsat' and redundancy time limits are namelist parameters
!                instead of fixed values and moved to data_obs_lib_cosmo.
! 2023-04-24 CS: Processing of ICOS tower reports added. 
!                Subr. 'obs_cdf_tower_copy' added for duplicating tower reports
!                to allow for processing with an alternative obs operator.
! 2023-05-09 CS: In redundancy check, add condition to prefer later obs time
!                (assuming to be closer to analysis time) over earlier obs time.
!                Discriminate between reports set passive solely due to obstype
!                from reports set passive ('rejected') for other reasons
! 2023-05-16 CS: Namelist parameters 'itim_wp', 'icdt_tws', 'rtmltow' added.
! 2023-05-22 CS: For wind profilers in redundancy check, prefer high mode 
!                (larger pulse length and signal-to-noise ratio) over low mode,
!                also to prevent higher vertical resolution of obs than of model
!                in mid-to-high troposphere.
! 2023-06-05 CS: In redundancy check, prefer Synop reports valid at full hour,
!                and avoid adding temporal non-local obs to (any) reports not
!                valid at full hour.
! 2024-06-14 CS: Re-write of 'obs_cdf_tower_copy' in order to allow for
!                exceptions from tower obs operator rule given by 'dhosag'.
! 2024-07-17 CS: Subr. 'obs_cdf_bcorr_rass' added for bias correction of RASS
!                virtual temperature (fixed values according to namelist input).
! 2024-10-07 CS: Namelist parameters 'icdt_rss' and 'lsytac' added.
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
!   c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
!   c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_wp : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_wp : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_wp : commonly used very small value > 0

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    ! horizontal and vertical sizes of the model fields
    ie_loc         ,& ! number of grid pts in zonal direction (local sub-domain)
    je_loc         ,& ! number of grid pts in meridional dir. (local sub-domain)
    ke             ,& ! number of grid pts in vertical direction (--> 'f_z2p')
    ie_tot         ,& ! number of grid pts in zonal direction (in total domain)
    je_tot         ,& ! number of grid pts in meridional dir. (in total domain)

    ! variables related to parallelisation / domain decomposition
    num_compute    ,& ! number of compute PEs
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers   ,& ! INTEGER   type used for MPI
    imp_character  ,& ! CHARACTER type used for MPI

    ! report dimensions of ODR (Observation Data Record, for observation storage
    ! on local sub-domains), and other variables related to namelist parameters
    maxmlv         ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    maxmll         ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl         ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl         ,& ! size (report dimension) of the (ground-based) GPS  ODR
!   maxtvl         ,& ! size (report dimension) of the satellite retrieval ODR
                      !
    nolbc          ,& ! number of grid rows at lateral boundaries
                      !   where obs are neglected
    itim_wp        ,& !   0 : mode of correction of obs time for wind profiler
    icdt_tws       ,& !   0 : mode for obs/code type of tower surface-level rep.
    icdt_rss       ,& !   0 : mode for obs/code type of surface report from TEMP
    madj_hum       ,& ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable
    lredn_repro    ,& ! ensure reproducibility of redundancy check irrespective
                      !   of domain decomposition by allowing for redundancy
                      !   only between reports assigned to the same grid point
    lsytac         ,& ! TAC preferred over BUFR Synop reports (or vice versa)
                      !
    mxav           ,& ! max. length of level definition list for superobbing
                      !   of high-resolution radisondes
    mxgpc          ,& ! max. number of GPS processing centres used
    irun_osse      ,& ! model run to derive obs values from file yfofin='fof'
    losse_fg       ,& ! f.g. check flag from 'fof' converted to 'dataset flag'
    iseed             ! external seed for random number generator

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
    acthr          ,& ! actual model time [hours] with respect to 'ydate_ref'
    fperturb       ,& ! factor to obs error variances to define size of random
                      ! perturbations added to the obs (only from yfofin='fof')
    av_reso        ,& ! apply superobbing if the averaged resolution of the obs
                      !   profile exceeds 'av_reso' times the model resolution
    rhtsat         ,& ! relative humidity threshold above which obs is set =100%
    rtmlrs         ,& ! redundancy time limit for radiosondes              [hrs]
    rtmlsy         ,& ! redundancy time limit for SYNOP                    [hrs]
    rtmlair        ,& ! redundancy time limit for AIREP                    [hrs]
    rtmltow        ,& ! redundancy time limit for tower                    [hrs]
    rtmlrsy        ,& ! redundancy time limit for raso sfc. level vs Synop [hrs]
    rtmlim         ,& ! redundancy time limit for other obs                [hrs]

    ! physical constants
    r_g            ,& ! acceleration due to gravity
    tmelt          ,& ! melting temperature of ice
    r_d            ,& ! gas constant for dry air
    rdv            ,& ! r_d / r_v
    o_m_rdv        ,& ! 1 - r_d/r_v
    rdocp          ,& ! r_d / cp_d
    b1             ,& ! variables for computing the saturation vapour pressure
    b2w            ,& ! over water (w) and ice (i)
    b2i            ,& !               -- " --
    b3             ,& !               -- " --
    b4w            ,& !               -- " --
    b4i            ,& !               -- " --

    ! switches related to namelist parameters and other
    lverpas        ,& ! write also passive reports on feedback files
    lwonl          ,& ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open
!   lgpsbias          ! .t. ==> bias correction to GPS IWV applied
    ydate_ref      ,& ! reference date (e.g. start of the forecast)
                      ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! 2b. Variables / fields as pointers (which need to be targeted)
! --------------------------------------------------------------

    i_gpscen       ,& ! array of processing centres of GPS reports used actively
    r_av_levs      ,& ! level definition list  \  for superobbing layers of
    r_av_incr      ,& ! level increment  list  /  high-res. radiosonde reports
                      ! -- model fields: --
!   r_p            ,& ! pressure (at main levels)
    r_hhl          ,& ! geometrical height (at half levels)
!   r_t_ll         ,& ! temperature at lowest model (main) level
!   r_ps           ,& ! surface pressure
    r_frland       ,& ! land fraction

! 2d: Quantities related to above pointers
!-----------------------------------------

    pref           ,& ! reference interpolation pressure levels for
                      ! superobbing of high-resolution radiosonde reports
    nref           ,& ! actual number of interpolation reference levels

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    yucautn    ,& ! caution messages if too many obs for ODR size
    nucautn       ! caution messages if too many obs for current ODR size

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP report
    nairep     ,& ! AIREP report (all aircraft reports)
    nsatob     ,& ! SATOB report
    ndribu     ,& ! DRIBU report
    ntemp      ,& ! TEMP  report
    npilot     ,& ! PILOT report
!   nsatem     ,& ! SATEM report
!   nsattv     ,& ! SATEM report
    ngps       ,& ! GPS   report
    nscatt     ,& ! SCATT report
    mxobtp     ,& ! number of observation types
!   nsrscd     ,& !   synop surface report
!   natscd     ,& !   automatic synop surface report
    nshscd     ,& !   ship synop report
!   nabscd     ,& !   ship synop abbreviated report
!   nshred     ,& !   shred report
    natshs     ,& !   automatic ship synop report
    nmetar     ,& !   Metar
    nsytmp        !   TEMP-derived surface report
!   naircd     ,& !   aircraft report
!   ncodar     ,& !   codar report
!   ncolba     ,& !   colba report
!   namdar     ,& !   amdar report
!   nacar      ,& !   acar  report
!   ndrbcd     ,& !   dribu report
!   nbathy     ,& !   bathy report
!   ntesac        !   tesac report

USE data_obs_lib_cosmo, ONLY :  &
    ntdesc     ,& !   temp high-res. descending report (HA 02.04.19: as EC)
    nldpcd     ,& !   pilot land   report
    nshpcd     ,& !   pilot ship   report
!   nmopcd     ,& !   pilot mobile report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    npr_us     ,& !   US Wind Profiler/RASS report
    ntower     ,& !   tower profile report
    ntowic     ,& !   icos tower profile report
    nwlidr        !   ground-based wind lidar report
!   nascat     ,& !   ASCAT scatterometer report
!   nqscat        !   QuickScat scatterometer report
!   nstmcd     ,& !   satem report
!   nstovs     ,& !   high resolution ATOVS satellite data
!   nsmsg1     ,& !   MSG_1  satellite retrieval  (nsmsg1 == 200 + idmsg1  )
!   nnoa15     ,& !   NOAA15 satellite retrieval  (nnoa15 == 200 + idnoaa15)
!   nnoa16     ,& !   NOAA16 satellite retrieval  (nnoa16 == 200 + idnoaa16)
!   nnoa17     ,& !   NOAA17 satellite retrieval  (nnoa17 == 200 + idnoaa17)
!   nnoa18     ,& !   NOAA18 satellite retrieval

! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &
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

    ntype_cdfin    ,& ! number of NetCDF observation input file types
    mxcdfin        ,& ! max. number of NetCDF observation input files
    mxfofin        ,& ! max. number of NetCDF feedobs (feedback) input files
    icdfinlen      ,& ! maximum length of NetCDF observation input file name
    iannexlen      ,& ! maximum length of annex of NetCDF obs input file name
    ncdf_temp      ,& ! indicator for processing of NetCDF TEMP       input
    ncdf_tempship  ,& ! indicator for processing of NetCDF TEMPSHIP   input
    ncdf_temphirs  ,& ! indicator for proc. NetCDF TEMP high-res BUFR input
!   ncdf_tempdrop  ,& ! indicator for proc. NetCDF TEMP Dropsonde     input
    ncdf_tempdesc  ,& ! indicator for proc. NetCDF descending TEMP    input
    ncdf_pilot     ,& ! indicator for proc. NetCDF PILOT (z-levels)   input
    ncdf_pilot_p   ,& ! indicator for proc. NetCDF PILOT (p-levels)   input
    ncdf_amdar_ml  ,& ! indicator for proc. NetCDF AMDAR multi-level  input
!   ncdf_amdar_vp  ,& ! indicator for proc. NetCDF AMDAR vert.profile input
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
    ncdf_nat_dwd   ,& ! indicator for proc. NetCDF DWD national surf. input
    ncdf_tower     ,& ! indicator for proc. NetCDF tower profile      input
    ncdf_tower_icos,& ! indicator for proc. NetCDF ICOS tower profile input
    ncdf_wlidar_wp ,& ! indicator for proc. NetCDF ground-b wind lidar input
    ncdf_synop_tst ,& ! indicator for proc. NetCDF test SYNOP         input
    ncdf_fof          ! indicator for proc. 'fof' feedobs (feedback) file input

USE data_obs_cdfin, ONLY :  &

    ycdfin         ,& ! file names of NetCDF observation input files
    yfofin         ,& ! file name of feedobs (feedback) input file(s) (='fof')
    n_cdfin        ,& ! number of existing NetCDF observation input files
    n_fofin        ,& ! number of existing NetCDF feedobs (feedback) input files
    icdfin         ,& ! obs file type of NetCDF observation input files
    ncinid         ,& ! unit numbers of NetCDF observation input files

    yncannex       ,& ! annex of NetCDF observation input file names
    yfofannex         ! annex of NetCDF feedobs (feedback) input file names
!   dimids            ! dimension IDs in NetCDF files

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
!-------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
    mxreve     ,& ! length of report event counter array
    neredn     ,& ! redundancy between 2 multi-level, or 2 single-level reports
    neslml     ,& ! one multi-level report made from other reports (either
                  !   from single-level aircraft, or from TEMP parts A,B,C,D)
    neslps     ,& ! report (either single-level aircraft, or TEMP part A,B,C,
                  !   or D) put in multi-level report and set passive
    neredx     ,& ! redundancy between 1 multi- and 1 single-level report
    nenoml     ,& ! multi-levl report not made due to ODR array size

!         3.1.2   Data event counter array format
!                 -------------------------------
    mxdeve     ,& ! length of data event counter array
    netlps     ,& ! lapse rate of multi-level temperature too large
    nefshr     ,& ! wind speed shear too large
    nedshr     ,& ! directional wind shear too large

!         3.2    Event counter arrays
!                --------------------
    neventr    ,& ! counter array of report events
    neventd    ,& ! counter array of data events

!         4.1    Observation error levels
!                ------------------------
!   nerlev     ,& ! number of standard error levels
    rlevel     ,& ! error levels
    rolnlv     ,& ! ln(rlevel(15))

!         5.1    Redundancy check limits
!                -----------------------
!   rtmlim     ,& ! time limit for all reports except AIREP      [hrs]
!   rtmlrs     ,& ! time limit for radiosondes (TEMP, PILOT)     [hrs]
!   rtmlsy     ,& ! time limit for SYNOP (> 20 min)              [hrs]
!   rtmlair    ,& ! time limit for reports of obs type 'AIREP'   [hrs]
!   rtmlrsy    ,& ! time limit for raso surface level vs Synop   [hrs]
                  !  (time of lowest level of multi-level ODR)
    rhzlim     ,& ! horiz.  dist. limit for obs sta. (\ AIREP)    [km]
    rhzlshp    ,& ! horiz.  dist. limit for ship/buoy/SYNOP       [km]
    rhzlair    ,& ! horizont. distance limit for AIREP reports    [km]
    rvtlim     ,& ! vertic. dist. limit for obs sta. (\ AIREP)     [m]
    rvtlsy     ,& ! vertic. dist. limit (sta.ht.) for SYNOP        [m]
    rvtlshp    ,& ! vertic. dist. limit for ship/buoy              [m]
    rvtlair    ,& ! vertical  distance limit for AIREP reports   [hpa]
    rdplim     ,& ! vertic. dist. limit within multi-level rep.  [hpa]
    rprlim     ,& ! vertic. dist. limit for replacing 'missing
                  !   data' within multi-level reports           [hpa]

!         5.2    Temperature / humidity / pressure / height / fog limits
!                -------------------------------------------------------
!   rhtsat     ,& ! rel. humidity threshold for saturation with real obs
    picoshift  ,& ! shift [pa] applied to colocated ICOS tower wind obs

!         5.4    Limits for the directional wind shear check
!                -------------------------------------------
    nnqcdd     ,& ! number of levels in the quality control threshold tables
    nqcdd      ,& ! thresholds for directional shear, as funct. of speed
    nqcddff    ,& ! limit values for sum of wind speeds (for 'nqcdd')
    qcfddff       ! wind speed limits for directional and speed shear checks

USE data_obs_cdfin, ONLY :  &

! Section 7 :  For reporting rejection of data: Output buffer, size and formats
!-------------------------------------------------------------------------------

    outbuf     ,& ! buffer containing output for a single node
    nacout     ,& ! actual number of records stored in the output buffer
    nmxoln     ,& ! maximum length of output buffer
    istrej     ,& ! length of strings (station id) in output buffer
    nfmt7      ,& ! excess of lapse rate
    nfmt8      ,& ! excess of wind speed shear
    nfmt9      ,& ! excess of directional shear
    nfmt10     ,& ! redundancy of surface-level report
    nfmt11     ,& ! redundancy of multi-level report
!   nfmt12     ,& ! redundancy of aircraft report
    nfmt13     ,& ! redundancy of wind
    nfmt14     ,& ! redundancy of temperature
    nfmt15     ,& ! redundancy of humidity
    nfmt16     ,& ! redundancy of pressure / height

! Section 8 : Temporary global model fields
!-------------------------------------------------------------------------------

    hsurf_tot  ,& ! total array of model surface height
    fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_cdfin

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &

!       1.1    ODR header format
!       ------------------------
!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
    mxrhed     ,& ! header length of multi-level reports
    mxshed     ,& ! header length of single-level reports
    mxghed     ,& ! header length of GPS reports
    nhilon     ,& ! longitude of observing station
    nhjlat     ,& ! latitude  of observing station
    nhalt      ,& ! station altitude [m]
    nhtime     ,& ! (exact) time of observation in forecast hours
!   nhsurf     ,& ! height of model grid pt. to which obs. is assigned
!   nhzio      ,& ! longitude of obs. station (or lowest datum) in grid pt. unit
!   nhzjo      ,& ! latitude  of obs. station in grid pt. units
!   nhsynt     ,& ! nominal (synoptic) time of observation in forecast hours
    nhsolz        ! solar zenith angle [deg]

USE data_obs_record, ONLY :   &

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
    mxrhdf     ,& ! header length of multi-level reports
    mxshdf     ,& ! header length of single-level reports
    mxghdf     ,& ! header length of GPS reports
!   nhio       ,& ! (local) x-coord. of grid pt. assigned to obs
!   nhjo       ,& ! (local) y-coord. of grid pt. assigned to obs
    nhitot     ,& ! global x-coord. of grid pt. assigned to obs
    nhjtot     ,& ! global y-coord. of grid pt. assigned to obs
    nhobtp     ,& ! observation type
    nhcode     ,& ! code type
    nhschr     ,& ! station characteristics                      (see 1.1.4)
    nhpass     ,& ! flag for report being set to 'passive'       (see 1.1.4)
!   nhqcfw     ,& ! threshold quality control flags for pressure, status of
                  ! verification output
    nhflag     ,& ! report flags (obs type, surface, altitude, station ID)
    nhcorr     ,& ! update sequence number (station correction indicator)
!   nhcat      ,& ! data     category (from BUFR Section 1)
!   nhcats     ,& ! data sub-category (from BUFR Section 1)
    nhkz       ,& ! DWD internal classification number (observation type)
    nhcent     ,& ! originating centre
!   nhstid     ,& ! station identity number
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
!   nhsyhr     ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    nhnlev     ,& ! number of obs. levels (for multi-level reports)
!   nhvqcf     ,& ! for satellite retrieval: threshold quality control flags
!   nhaexi     ,& ! flag for exist. of wind or temperature in multi-level report
    nhuexi     ,& ! flag for existence of wind data        in multi-level report
    nhtexi     ,& ! flag for existence of temperature data in multi-level report
    nhqexi     ,& ! flag for existence of humidity data    in multi-level report
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)

!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------
    ilstid     ,& ! character length of the station identity
    ilstidp    ,& ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!              -----------------------------------------------------------
    nvpsbp     ,& ! bit pos. for report set passive since 1 of next   nhschr
                  !              5 flags or flight track flag applies
    nvrdbp     ,& ! bit pos. for flag: 'redundant report'               "
    nvlzbp        ! bit pos. for vertical levels reported as height     "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------
!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd'
!              ------------------------------------------------------
!   maxrsl     ,& ! max. number of levels in multi-level ODR


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
!   nbtzio     ,& ! longitude in geographical coordinates
!   nbtzjo     ,& ! latitude  in geographical coordinates
!   nbttim     ,& ! observation time relative to report (header) time
    nbtlop     ,& ! LOG( pressure )
    nbtplv     ,& ! pressure [Pa] level independent from model state
                  !   (if not reported, then derived from height by std. atm.)
    nbtdrh     ,& ! bias correction for relative humidity [/] or RASS Tv
    nbtw       ,& ! vertical velocity [m/s]
    nbtsnr     ,& ! signal to noise ratio
    nbtuac     ,& ! accuracy (std dev from data provider) of horiz. wind [m/s]


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
    mxrbdf     ,& ! body length of multi-level reports
    nbtflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
    nbterr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
!   nbtqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbtlsg     ,& ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
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
    nbscbs     ,& ! (lowest) cloud base height CBH above surface       [m]
    nbscil     ,& ! ceiling (CBH of lowest cloud layer >= 5 octas)     [m]
    nbscl      ,& ! low       cloud cover        (BUFR Table 020011)   [octas]
    nbscm      ,& ! mid-level cloud cover        (BUFR Table 020011)   [octas]
    nbsch      ,& ! high      cloud cover        (BUFR Table 020011)   [octas]
    nbsct      ,& ! total     cloud cover        (BUFR Table 020011)   [octas]
    nbsvis     ,& ! (horizontal) visibility                            [m]
    nbsff      ,& ! wind speed                                         [m/s]
    nbsdd      ,& ! wind direction                                     [deg]
    nbstd      ,& ! dewpoint temperature                               [K]
    nbsrr1     ,& ! precipitation amount over 1 hour                   [mm]
    nbsrr3     ,& ! precipitation amount over 3 hours                  [mm]
    nbsrr6     ,& ! precipitation amount over 6 hours                  [mm]
    nbsr12     ,& ! precipitation amount over 12 hours                 [mm]
    nbsr24     ,& ! precipitation amount over 24 hours                 [mm]
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
    mxsbdf     ,& ! body length of single-level reports
    nbsflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
!   nbsqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
    nbserr     ,& ! status flag word        (bit pattern, see below: 'nb?err')
!   nbslid     ,& ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
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
    nbsttr        ! time periods    (bit pattern of values, see below: nbsttr)

USE data_obs_record, ONLY :   &

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
    mxgbdy     ,& ! body length of GPS reports
    nbgtze     ,& ! error in total zenith delay [mm]
    nbgzpd     ,& ! zenith path delay (total zenith delay)             [mm]
    nbgzwd     ,& ! zenith wet delay [mm]
    nbgiwv     ,& ! integrated water vapour [mm]
    nbgp       ,& ! pressure [Pa]
!   nbgt       ,& ! temperature [K]
!   nbgrh      ,& ! relative humidity [/]
!   nbgbia     ,& ! bias correction to integrated water vapour [mm]
!   nbgiwa     ,& ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
    mxgbdf     ,& ! body length of GPS reports
!   nbgflg     ,& ! main flag word          (bit pattern, see below: 'nb?flg')
!   nbgqcf     ,& ! threshold quality control flags      (see below: 'nb?qcf')
!   nbglid     ,& ! level identity          (bit pattern, see below: 'nb?lid')
    nbgerr        ! pre-processing status flag

USE data_obs_record, ONLY :   &

!       1.4    Bit patterns for packed information in ODR (and VOF) body
!       ------------------------------------------------------------------------
!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru       ,& ! bit pos. for status 'active' for horiz. wind    nb?err
    nvrt       ,& ! bit pos. for status 'active' for temperature      "
    nvrq       ,& ! bit pos. for status 'active' for humidity         "
    nvrz       ,& ! bit pos. for status 'active' for pressure/height  "
    nvriwv     ,& ! bit pos. for status 'active' for ZPD
!   nvrw       ,& ! bit pos. for status 'active' for vertical wind    "
    nvfubp     ,& ! bit pos. for main flag on wind                  nb?flg
    nvftbp     ,& ! bit pos. for main flag on temperature             "
    nvfqbp     ,& ! bit pos. for main flag on humidity                "
    nvfzbp     ,& ! bit pos. for main flag on pressure / geopot.      "
    nvfaoc     ,& ! no. of bits occ. by each main flag                "
    nvfbps     ,& ! bit pattern for main flags:                       "
    nvfboc     ,& ! no. of bits occ. for main flags                   "
!   nvflbp     ,& ! bit pos. for level flag: level below surface      "
    nvffbp     ,& ! bit pos. for flag: tower obs operator for wind    "
                  !                    disagrees with 'dhosag'
    nvfmbp     ,& ! bit pos. for flag: tower obs operator for mass    "
                  !                    (T, qv) disagrees with 'dhosag'
    ntxbp      ,& ! bit pos. for time period for T-max (1,..,24 hr)   "
    ntnbp      ,& ! bit pos. for time period for T-min (1,..,24 hr)   "
    ntxoc      ,& ! no. of bits occ. by each time period              "
    nvlidp        ! level id. bit pattern                           nb?lid
!   nvlido        ! no. bits occ. by each indicator in level id.      "

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------
    imdi       ,& ! missing data indicator for ODR integers (2^31-1)
    ntotml     ,& ! tot. number of stored multi-level reports
    ntotsg     ,& ! tot. number of stored single-level reports
    ntotgp     ,& ! tot. number of stored GPS reports

!       2.     Observation data records (ODR) and associated arrays
!       ----------------------------------------------------------------------

!       2.1    Formats of observation data records
!              -----------------------------------
    omlbdy     ,& ! body   of multi-level ODR
    omlhed     ,& ! header of multi-level ODR
    osgbdy     ,& ! body   of single-level ODR
    osghed     ,& ! header of single-level ODR
    ogpbdy     ,& ! body   of GPS ODR
    ogphed     ,& ! header of GPS ODR
    momlbd     ,& ! body   of multi-level ODR
    momlhd     ,& ! header of multi-level ODR
    mosgbd     ,& ! body   of single-level ODR
    mosghd     ,& ! header of single-level ODR
    mogpbd     ,& ! body   of GPS ODR
    mogphd     ,& ! header of GPS ODR
    yomlhd     ,& ! header of multi-level ODR
    yosghd     ,& ! header of single-level ODR
    yogphd     ,& ! header of GPS ODR

!       3.     Masking constants
!       ------------------------

    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables,          ONLY :  &

    FL_OBSTYPE    ,& ! passive report type (at obs. location)
    FL_REDUNDANT  ,& ! redundant report (COSMO: also TEMP ABCD parts merged
                     !                          in a complete profile)
    FL_THIN       ,& ! thinning
    RS_92_DIG12   ,& ! Vaisala RS 92 Digicora I,II or Marvin
    RS_92_DIG3    ,& ! Vaisala RS 92 Digicora III
    RS_92_AUTO    ,& ! Vaisala RS 92 Autosonde
    LS_SURFACE    ,& ! surface
    LS_STANDARD   ,& ! standard level
    LS_TROPO      ,& ! tropopause level
    LS_MAX        ,& ! maximum wind level
!   LS_SIGN       ,& ! significant level
    LS_SUPEROBS      ! superobservation (layer average)

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values   ,& ! computes global values by operating on local arrays
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    get_global_surface     ,& ! get orography / land fraction on global domain
    obs_solar_zenith_angle ,& ! get solar zenith angle for an array of reports
    obs_rm_duplicate_levels   ! remove duplicate levels

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_mult,          ONLY :  &
    obs_cdf_read_temp_pilot,& ! read, pre-process, store radiosonde (TEMP/PILOT)
    obs_cdf_read_profiler     ! read, pre-process, store remote-sensing profiler

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_sing,          ONLY :  &
    obs_cdf_read_surface   ,& ! read, pre-process, store surface-level reports
    obs_cdf_read_aircraft  ,& ! read, pre-process, store aircraft reports
    obs_cdf_read_amv          ! read, pre-process, store SATOB / AMV reports

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_gps,            ONLY :  &
    obs_cdf_read_gps          ! read, pre-process, store GPS reports

!-------------------------------------------------------------------------------

#ifdef __COSMO__
 USE src_obs_fdbk_in,              ONLY :  &
    obs_fof_read              ! read, pre-process, store reports from 'fof' file
#endif

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE
PRIVATE
          ! public subroutines
PUBLIC :: obs_cdf_read_org,      obs_cdf_interface,   obs_cdf_proc_init,       &
          obs_cdfin_open_files,  obs_cdf_tower_copy,  obs_cdf_mult_qualicheck, &
          obs_cdf_bcorr_rass,    obs_cdf_raso_rh_bias

! for nudging in COSMO model code only:
PUBLIC :: obs_fof_check_files

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_read_org ( lneedob, min_sta, min_end, ycdfdir, icdfdirlen   &
                            , nodrold, nexceed, ldo_airmult, ldo_pr_out )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_cdfin_org" organizes the reading,
!   pre-selection, distribution, and storage in ODR of observation reports.
!
! Method:
!   According to the input parameters 'lneedob', 'min_sta', 'min_end',
!   new reports from certain periods are read from certain NetCDF observation
!   input files.
!   The observations are assigned to model grid points, distributed to the
!   appropriate nodes (sub-domains, which contain the assigned grid points),
!   pre-processed and stored in the internal storage arrays ODR.
!   The pre-processing of observations includes redundancy checking
!   (including putting TEMP parts together), gross error checking, assignment
!   of observation errors, etc. Furthermore, the ODR is cleaned from reports
!   which are not used any more, and writing of a range of monitoring /
!   diagnostic (ASCII) output is prepared.
!   The routine however does not include initialisations (such as reading
!   the blacklist file), thinning and aggregation of single-level aircraft
!   reports to piecewise vertical profiles, and flight track checking.
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

  INTEGER        , INTENT (IN)    ::  &
    min_sta  (:),& ! start of reading time interval [min]
    min_end  (:),& ! end of reading time interval [min]
    icdfdirlen     ! max. length of name of directory where
                   !   NetCDF observation input files reside

  LOGICAL        , INTENT (IN)    ::  &
    lneedob  (:)   ! new observation need to be read at current timestep

  CHARACTER (LEN= *)  , INTENT (IN)    ::  &
    ycdfdir          ! directory where NetCDF obs input + feedobs files reside

  INTEGER        , INTENT (INOUT) ::  &
    nodrold  (3)   ! number of multi-level / single-level / GPS reports
                   ! before having read new data at the current timestep
                   !   (can be modified in 'obs_cdf_del_old_rep')

  INTEGER        , INTENT (OUT  ) ::  &
    nexceed  (3)   ! number of reports in excess of ODR array size

  LOGICAL        , INTENT (INOUT) ::  &
    ldo_airmult ,& ! after leaving this routine do call 'obs_air_org_mult'
                   ! (to produce multi-level aircraft reports)
                   ! because new aircraft reports have been read
    ldo_pr_out     ! after leaving this routine do call 'obs_cdf_print_odr'
                   ! (to write control output on rejected/processed reps.)
                   ! because new reports have been read

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    kcdf        ,& ! obs file type of NetCDF observation input file
    nodrnew (4) ,& ! number of new (multi-level/ single-level/ GPS) reports
    ntotnew (2) ,& ! (total/aircraft) number of new reports (sum over loop)
    ilcf        ,& ! loop index
    istat , ierr   ! error status variables

  CHARACTER (LEN=20)       ::  &
    yroutine       ! name of this subroutine
  CHARACTER (LEN=30)       ::  &
    yerr           ! error message

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_read_org
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_read_org'

  IF (     (SIZE( lneedob ) /= SIZE( min_sta ))                                &
      .OR. (SIZE( lneedob ) /= SIZE( min_end )))                               &
    PRINT *, 'ERROR in obs_cdf_read_org, incorrect SIZE of arrays'

!-------------------------------------------------------------------------------
! Section 3: Read from NetCDF files, select, assign to grid points, and
!            distribute observational reports to nodes (processors) according
!            to the sub-domains which contain the assigned grid points;
!            then pre-process and store the reports in the ODR arrays;
!            finally check for redundancy of reports; also clean up the ODR
!            arrays from reports which are not used any more
!-------------------------------------------------------------------------------

! compose full fields of model surface height and land fraction
! (to prepare the assignment of observation to grid points)
! -------------------------------------------------------------
  ALLOCATE ( hsurf_tot (ie_tot,je_tot)      , STAT=istat )
  ALLOCATE ( fland_tot (ie_tot,je_tot)      , STAT=istat )

  CALL get_global_surface ( num_compute     , my_cart_id                       &
                          , ie_loc , je_loc , r_hhl(:,:,ke+1) , r_frland       &
                          , ie_tot , je_tot , hsurf_tot       , fland_tot )
! =======================

! delete observations which have been marked to be too old to be used any more
! ----------------------------------------------------------------------------

! (note: the numbers of reports 'nodrold' which have been read prior
!        to the current timestep, and the total number of reports
!        (ntotml, ntotsg, ntotgp) are updated in 'obs_cdf_del_old_rep')

  CALL obs_cdf_del_old_rep ( nodrold )
! ========================

! read the observations from different types of NetCDF files
! (and pre-process and store them as outlined above)
! ----------------------------------------------------------

  nexceed (:) = 0
  ntotnew (:) = 0

  DO ilcf = 1 , n_cdfin + n_fofin
    IF (ilcf <= n_cdfin)  kcdf  =  icdfin(ilcf)
    IF (ilcf >  n_cdfin)  kcdf  =  ncdf_fof
    IF (lneedob(ilcf)) THEN
!     IF (my_cart_id == 0) PRINT *,'need obs_cdf_read ', ilcf, kcdf            &
!                                 , lneedob(ilcf), min_end(ilcf), min_sta(ilcf)
      nodrnew (:) =  0

      IF (     (kcdf == ncdf_temp    ) .OR. (kcdf == ncdf_tempship)            &
          .OR. (kcdf == ncdf_temphirs) .OR. (kcdf == ncdf_tempdesc)            &
          .OR. (kcdf == ncdf_pilot   ) .OR. (kcdf == ncdf_pilot_p )) THEN
!                                      .OR. (kcdf == ncdf_tempdrop)            &

        CALL obs_cdf_read_temp_pilot ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(1:2) , nexceed(1:2) )
!       ============================

      ELSEIF (     (kcdf == ncdf_synop  ) .OR. (kcdf == ncdf_synop_mob)        &
              .OR. (kcdf == ncdf_ship   ) .OR. (kcdf == ncdf_buoy )            &
              .OR. (kcdf == ncdf_ascat  ) .OR. (kcdf == ncdf_qscat)            &
              .OR. (kcdf == ncdf_nat_dwd) .OR. (kcdf == ncdf_synop_tst)) THEN
!             .OR. (kcdf == ncdf_metar)) THEN

        CALL obs_cdf_read_surface    ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(  2) , nexceed(  2) )
!       =========================

      ELSEIF (     (kcdf == ncdf_amdar)    .OR. (kcdf == ncdf_acars)           &
              .OR. (kcdf == ncdf_acars_uk) .OR. (kcdf == ncdf_acars_us)        &
              .OR. (kcdf == ncdf_modes)    .OR. (kcdf == ncdf_modes_acr)       &
              .OR. (kcdf == ncdf_amdar_ml)) THEN
!             .OR. (kcdf == ncdf_amdar_ml) .OR. (kcdf == ncdf_amdar_vp)) THEN

        CALL obs_cdf_read_aircraft   ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(1:2) , nexceed(1:2) )
!       ==========================
        ntotnew (2)  =  ntotnew(2)  +  nodrnew(2)

      ELSEIF (     (kcdf == ncdf_wprof) .OR. (kcdf == ncdf_radar_vad )         &
              .OR. (kcdf == ncdf_rass ) .OR. (kcdf == ncdf_wlidar_wp )         &
              .OR. (kcdf == ncdf_tower) .OR. (kcdf == ncdf_tower_icos)) THEN

        CALL obs_cdf_read_profiler   ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(1:2) , nexceed(1:2) )
!       ==========================

      ELSEIF       (kcdf == ncdf_satob)  THEN

        CALL obs_cdf_read_amv        ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(  2) , nexceed(  2) )
!       =====================

      ELSEIF       (kcdf == ncdf_gps_zenith)  THEN

        CALL obs_cdf_read_gps        ( min_sta(ilcf) , min_end(ilcf) , ilcf    &
                                     , nodrnew(  3) , nexceed(  3) )
!       =====================

#ifdef __COSMO__
      ELSEIF       (kcdf == ncdf_fof)  THEN

        CALL obs_fof_read      ( min_sta(ilcf) , min_end(ilcf) , ilcf-n_cdfin  &
                               , ycdfdir , icdfdirlen , nodrnew(1:4) , nexceed(1:3) )
!       =================
        ntotnew (2)  =  ntotnew(2)  +  nodrnew(4)
#endif

      ELSE
        IF (my_cart_id == 0)                                                   &
          PRINT *,'NOTE: obs file type ', kcdf, 'not available'
      ENDIF

! 'nexceed'   is determined by the preceeding routines
! 'nodrnew'   is determined in 'obs_cdf_store_comhead' (called by routines above)
! 'ntotml', 'ntotsg', 'ntotgp'  are updated below and in 'obs_cdf_del_old_rep'
! 'nodrold'   is updated in 'obs_cdf_del_old_rep' only

      ntotml      =  ntotml      +  nodrnew(1)
      ntotsg      =  ntotsg      +  nodrnew(2)
      ntotgp      =  ntotgp      +  nodrnew(3)
      ntotnew (1) =  ntotnew(1)  +  nodrnew(1) + nodrnew(2) + nodrnew(3)
!     PRINT *,'after obs_cdf_read ', nodrnew, ntotml, ntotsg, ntotgp, nodrold

! redundancy checking (flag redundant reports)
! -------------------
!     IF (kcdf /= ncdf_fof)                                                    &

      CALL obs_cdf_redundancy ( nodrnew(1), nodrnew(2), nodrnew(3), nodrold(2) &
                              , maxmlv )
!     =======================

! clean up ODR (delete rejected reports, defragment ODR)
! ------------  ('ntotml', 'ntotsg', 'ntotgp', 'nodrold' are also updated)

      CALL obs_cdf_del_old_rep ( nodrold )
!     ========================

!     IF (MAXVAL( nodrnew ) > 0)                                               &
!       PRINT *,'after obs_del_oldr ', nodrnew, ntotml, ntotsg, ntotgp         &
!                                    , ntotnew, nodrold(2)
    ENDIF
!   IF (ilcf == n_cdfin) THEN
!     nodrcdf (1)  =  ntotml
!     nodrcdf (2)  =  ntotsg
!     nodrcdf (3)  =  ntotgp
!   ENDIF
  ENDDO

  DEALLOCATE ( hsurf_tot , STAT=istat )
  DEALLOCATE ( fland_tot , STAT=istat )

  IF (num_compute > 1) THEN

    CALL global_values ( ntotnew, 2, 'SUM',imp_integers, icomm_cart, -1,yerr,ierr)
    ! ------------------
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11015, yerr, yroutine)
  ENDIF

  IF (ntotnew(1) > 0)  ldo_pr_out  = .TRUE.
  IF (ntotnew(2) > 0)  ldo_airmult = .TRUE.

! compute solar zenith angle for new reports read from 'cdfin' files
! ------------------------------------------------------------------

! IF (nodrcdf(1) > nodrold(1))                                                 &
  IF (ntotml > nodrold(1))                                                     &

    CALL obs_solar_zenith_angle ( ntotml - nodrold(1) , rmdich                 &
                                , momlhd(nodrold(1)+1:ntotml,nhdate)           &
                                , momlhd(nodrold(1)+1:ntotml,nhhrmn)           &
                                , omlhed(nodrold(1)+1:ntotml,nhjlat)           &
                                , omlhed(nodrold(1)+1:ntotml,nhilon)           &
                                , omlhed(nodrold(1)+1:ntotml,nhsolz) )
!   ===========================

! IF (nodrcdf(2) > nodrold(2))                                                 &
  IF (ntotsg > nodrold(2))                                                     &

    CALL obs_solar_zenith_angle ( ntotsg - nodrold(2) , rmdich                 &
                                , mosghd(nodrold(2)+1:ntotsg,nhdate)           &
                                , mosghd(nodrold(2)+1:ntotsg,nhhrmn)           &
                                , osghed(nodrold(2)+1:ntotsg,nhjlat)           &
                                , osghed(nodrold(2)+1:ntotsg,nhilon)           &
                                , osghed(nodrold(2)+1:ntotsg,nhsolz) )
!   ===========================

! IF (nodrcdf(3) > nodrold(3))                                                 &
  IF (ntotgp > nodrold(3))                                                     &

    CALL obs_solar_zenith_angle ( ntotgp - nodrold(3) , rmdich                 &
                                , mogphd(nodrold(3)+1:ntotgp,nhdate)           &
                                , mogphd(nodrold(3)+1:ntotgp,nhhrmn)           &
                                , ogphed(nodrold(3)+1:ntotgp,nhjlat)           &
                                , ogphed(nodrold(3)+1:ntotgp,nhilon)           &
                                , ogphed(nodrold(3)+1:ntotgp,nhsolz) )
!   ===========================

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_read_org
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_read_org


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org.f90" for interface to model environmt.
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_interface ( iie, ije, ike, iie_tot, ije_tot                 &
                             , npe, nlbc_lines, my_id, icomm                   &
                             , mpi_reals, mpi_integers, mpi_character          &
                             , zpollon, zpollat, zpolgam, zdlon, zdlat         &
                             , zstartlat_tot, zstartlon_tot, zdegrad           &
                             , imaxmll, imaxsgl, imaxgpl, imaxmlv              &
                             , inolbc, imadj_hum, zdoromx, zaltopsu, zzlimv10  &
                             , zav_reso, zrhtsat, iitim_wp,iicdt_tws,iicdt_rss &
                             , zrtmlrs, zrtmlsy, zrtmlair, zrtmltow            &
                             , zrtmlrsy, zrtmlim, llredn_repro, llsytac, zacthr&
                             , krun_osse, llosse_fg, zfperturb, iseed_ex       &
                             , z_g, z_t0_melt, z_r_d, z_rdv, z_rdocp           &
                             , z_b1, z_b2w, z_b2i, z_b3, z_b4w, z_b4i          &
                             , llverpas, llwonl, yydate_ref, maxgpc, maxav )

!-------------------------------------------------------------------------------
! Description:
!   This procedure of module "src_obs_cdfin_org.f90" makes the values of those
!   input variables (and fields and parameters)
!   - (1) which must be commonly available in the observation pre-processing
!         library and the observation operator library related to COSMO model
!         applications, but
!   - (2) which are set outside of these libraries in the 'model' environment
!   commonly available within these libraries.
!   These variables are stored in module 'data_obs_lib_cosmo'.
!
!   Remarks:
!   - The values of all the other variables which are set outside the libraries
!     but which must be available within a certain subroutine of the libraries
!     are transfered directly from the 'model' environment through the parameter
!     list of the respective subroutines.
!   - All the other variables which must be commonly available within the
!     libraries but can be set within the libraries are stored in one of the
!     data modules that are part of the libraries.
!   Note: Only the model fields and 'acthr' are variable in time.
!
! Method:
!   This subroutine is called from outside the obs processing / obs operator
!   libraries. In the call, the parameter list consists of variables from the
!   'model' environment. Within this routine, the values in the parameter list
!   are assigned to the variables stored in module 'data_obs_lib_cosmo'.
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

  INTEGER        , INTENT (IN)   ::  &

! horizontal and vertical sizes of the model fields
    ! iie    , ije    : used in 'get_global_surface', 'obs_cdf_interface'
    ! ike             : used in 'f_z2p'
    ! iie_tot, ije_tot: used in 'get_global_surface', 'obs_cdf_read_org',
    !                           'obs_assign_gridpt'
    iie            ,& ! number of grid pts in zonal direction (local sub-domain)
    ije            ,& ! number of grid pts in meridional dir. (local sub-domain)
    ike            ,& ! number of grid pts in vertical direction (--> 'f_z2p')
    iie_tot        ,& ! number of grid pts in zonal direction (in total domain)
    ije_tot        ,& ! number of grid pts in meridional dir. (in total domain)

! variables related to parallelisation / domain decomposition
    npe            ,& ! number of compute PEs
    nlbc_lines     ,& ! number of overlapping boundary lines of the subdomains
    my_id          ,& ! rank of this subdomain in the cartesian communicator
    icomm          ,& ! communicator for the virtual cartesian topology
    mpi_reals      ,& ! REAL      type used for MPI
    mpi_integers   ,& ! INTEGER   type used for MPI
    mpi_character  ,& ! CHARACTER type used for MPI

! report dimensions of ODR (Observation Data Record, for observation storage on
! local sub-domains), and other variables related to namelist parameters
    imaxmlv        ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    imaxmll        ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    imaxsgl        ,& ! size (report dimension) of the single-level (s-l)  ODR
    imaxgpl        ,& ! size (report dimension) of the (ground-based) GPS  ODR
!   imaxtvl        ,& ! size (report dimension) of the satellite retrieval ODR
                      !
    inolbc         ,& ! number of grid rows at lateral boundaries
                      !   where obs are neglected
    iitim_wp       ,& ! mode of correction of obs time for wind profiler
    iicdt_tws      ,& ! mode for obs/code type of tower surface-level reports
    iicdt_rss      ,& ! mode for obs/code type of surface report from TEMP
    imadj_hum      ,& ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable
    maxav          ,& ! max. length of level definition list for superobbing
                      !       of high-resolution radisondes
    maxgpc         ,& ! max. number of GPS processing centres used
    krun_osse      ,& ! model run to derive obs values from file yfofin='fof'
    iseed_ex          ! external seed for random number generator

  REAL (KIND=wp) , INTENT (IN)   ::  &

! constants for the horizontal rotated grid and related variables
    zpollon        ,& ! longitude of the rotated north pole (in degrees, E>0)
    zpollat        ,& ! latitude of the rotated north pole (in degrees, N>0)
    zpolgam        ,& ! angle between the north poles of the systems
    zdlon          ,& ! grid point distance in zonal direction (in degrees)
    zdlat          ,& ! grid point distance in meridional direction (in degrees)
    zstartlat_tot  ,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    zstartlon_tot  ,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    zdegrad        ,& ! factor for transforming degree to rad

! variables related to namelist parameters
    zdoromx    (4) ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    zaltopsu   (4) ,& ! SYNOP obs. above height 'altopsu' are set passive
    zzlimv10   (3) ,& ! additional limits for the use of 10-m wind observations:
                      !   pos./neg. Laplacian of orography, surface roughness
    zfperturb      ,& !  factor to obs error variances to define size of random
                      !  perturbations added to the obs (only from yfofin='fof')
    zav_reso       ,& ! apply superobbing if the averaged resolution of the obs
                      !   profile exceeds 'av_reso' times the model resolution
    zrhtsat    (5) ,& ! relative humidity threshold above which obs is set =100%
    zrtmlrs        ,& ! redundancy time limit for radiosondes              [hrs]
    zrtmlsy        ,& ! redundancy time limit for SYNOP                    [hrs]
    zrtmlair       ,& ! redundancy time limit for AIREP                    [hrs]
    zrtmltow       ,& ! redundancy time limit for tower                    [hrs]
    zrtmlrsy       ,& ! redundancy time limit for raso sfc. level vs Synop [hrs]
    zrtmlim        ,& ! redundancy time limit for other obs                [hrs]

! other variables
    zacthr         ,& ! actual model time [hours] with respect to 'yydate_ref'

! physical constants
    z_g            ,& ! acceleration due to gravity
    z_t0_melt      ,& ! melting temperature of ice
    z_r_d          ,& ! gas constant for dry air
    z_rdv          ,& ! r_d / r_v
!   z_o_m_rdv      ,& ! 1 - r_d/r_v
    z_rdocp        ,& ! r_d / cp_d
    z_b1           ,& ! variables for computing the saturation vapour pressure
    z_b2w          ,& ! over water (w) and ice (i)
    z_b2i          ,& !               -- " --
    z_b3           ,& !               -- " --
    z_b4w          ,& !               -- " --
    z_b4i             !               -- " --

  LOGICAL        , INTENT (IN)   ::  &
    llverpas       ,& ! write also passive reports on feedback files
    llosse_fg      ,& ! f.g. check flag from 'fof' converted to 'dataset flag'
    llredn_repro   ,& ! ensure reproducibility of redundancy check irrespective
                      !   of domain decomposition by allowing for redundancy
                      !   only between reports assigned to the same grid point
    llsytac        ,& ! TAC preferred over BUFR Synop reports (or vice versa)
    llwonl            ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open

  CHARACTER (LEN= 14)  , INTENT (IN)   ::  &
    yydate_ref        ! reference date (e.g. start of the forecast)
                      ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! Local parameters: None
! ----------------

! Local scalars:
! -------------

! INTEGER        ::  &
!   i , j , k        ,& ! loop indices
!   istat               ! error status variables

! CHARACTER (LEN=20)       ::  &
!   yroutine            ! name of this subroutine

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_interface
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_interface'

!-------------------------------------------------------------------------------
! Section 1: Copy values to the common variables available in the observation
!            processing library
!-------------------------------------------------------------------------------

! horizontal and vertical sizes of the model fields
! -------------------------------------------------

  ie_loc  =  iie
  je_loc  =  ije
  ke      =  ike
  ie_tot  =  iie_tot
  je_tot  =  ije_tot

! constants for the horizontal rotated grid and related variables
! ---------------------------------------------------------------

  r_pollon       =  zpollon
  r_pollat       =  zpollat
  r_polgam       =  zpolgam
  r_dlon         =  zdlon
  r_dlat         =  zdlat
  r_startlat_tot =  zstartlat_tot
  r_startlon_tot =  zstartlon_tot
  r_degrad       =  zdegrad

! variables related to parallelisation / domain decomposition
! -----------------------------------------------------------

  num_compute    =  npe
  nboundlines    =  nlbc_lines
  my_cart_id     =  my_id
  icomm_cart     =  icomm
  imp_reals      =  mpi_reals
  imp_integers   =  mpi_integers
  imp_character  =  mpi_character

! report dimensions of ODR (related to namelist parameters)
! ------------------------

  maxmlv  =  imaxmlv
  maxmll  =  imaxmll
  maxsgl  =  imaxsgl
  maxgpl  =  imaxgpl
! maxtvl  =  imaxtvl

! other variables related to namelist parameters
! ----------------------------------------------

  nolbc       =  inolbc
  itim_wp     =  iitim_wp 
  icdt_tws    =  iicdt_tws
  icdt_rss    =  iicdt_rss
  madj_hum    =  imadj_hum
  doromx      =  zdoromx
  altopsu     =  zaltopsu
  zlimv10     =  zzlimv10
  av_reso     =  zav_reso
  rhtsat      =  zrhtsat
  rtmlrs      =  zrtmlrs
  rtmlsy      =  zrtmlsy
  rtmlair     =  zrtmlair
  rtmltow     =  zrtmltow
  rtmlrsy     =  zrtmlrsy
  rtmlim      =  zrtmlim
  lverpas     =  llverpas
  lredn_repro =  llredn_repro
  lsytac      =  llsytac
  ydate_ref   =  yydate_ref
  mxav        =  maxav
  mxgpc       =  maxgpc
  irun_osse   =  krun_osse
  losse_fg    =  llosse_fg
  fperturb    =  zfperturb
  iseed       =  iseed_ex

! other variables
! ---------------

  acthr   =  zacthr
  lwonl   =  llwonl

! physical constants
! ------------------

  r_g     =  z_g
  tmelt   =  z_t0_melt
  r_d     =  z_r_d
  rdv     =  z_rdv
  o_m_rdv =  1.0_wp - rdv
  rdocp   =  z_rdocp
  b1      =  z_b1
  b2w     =  z_b2w
  b2i     =  z_b2i
  b3      =  z_b3
  b4w     =  z_b4w
  b4i     =  z_b4i

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_interface
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_interface


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for initialising the obs processing
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_proc_init ( lobpfrs , ndimev2 , nuout_in, yuout_in )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" initialises the
!   observation processing.
!
! Method:
!   Allocation and initialization of long term storage arrays.
!   Opening of files for control output.
!
! Written by        :  Christoph Schraff, DWD  (original version: 23.08.04, part
!                                               of previous 'organize_obs_proc')
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
  LOGICAL           , INTENT (IN)              :: &
    lobpfrs          ! variable for 'first time of obs. processing'

  INTEGER           , INTENT (IN)              :: &
    ndimev2          ! length of 2nd dimenstion of event counter arrays

  INTEGER           , INTENT (IN)  , OPTIONAL  :: &
    nuout_in         ! unit number of output file

  CHARACTER (LEN=8) , INTENT (IN)  , OPTIONAL  :: &
    yuout_in         ! file name for unit number 'nuout_in'

! Local parameters: None
! ----------------

! Local variables:
! ----------------
  INTEGER             ::  &
    istat                  ! error status variable

  CHARACTER (LEN=20)  ::  &
    yroutine               ! name of this subroutine
! CHARACTER (LEN=40)  ::  &
!   yerrmsg                ! error message
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_proc_init
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdf_proc_init'

!-------------------------------------------------------------------------------
!  Section  1.   Allocation and initialization of long term storage arrays
!-------------------------------------------------------------------------------

  IF (lobpfrs)  THEN

    ! report and data event counters
    ! ------------------------------
    IF (ALLOCATED( neventr )) THEN
      PRINT '("CAUTION in src_obs_cdfin_org: neventr is already allocated "    &
            &,"at time / PE",F6.3,I5)', acthr, my_cart_id
      DEALLOCATE ( neventr , neventd , STAT=istat )
    ENDIF
    ALLOCATE ( neventr( mxreve, ndimev2 ), STAT=istat )
    ALLOCATE ( neventd( mxdeve, ndimev2 ), STAT=istat )

    neventr (:,:) = 0
    neventd (:,:) = 0

    ! Observation data record (ODR)
    ! -----------------------------

    IF (ALLOCATED( omlbdy )) THEN
      PRINT '("CAUTION in src_obs_cdfin_org: omlbdy is already allocated "     &
            &,"at time / PE",F6.3,I5)', acthr, my_cart_id
      DEALLOCATE ( omlbdy , omlhed , momlbd , momlhd , yomlhd , STAT=istat )
      DEALLOCATE ( osgbdy , osghed , mosgbd , mosghd , yosghd , STAT=istat )
      DEALLOCATE ( ogpbdy , ogphed , mogpbd , mogphd , yogphd , STAT=istat )
    ENDIF
    ALLOCATE ( omlbdy( maxmll, maxmlv, mxrbdy), STAT=istat)
    ALLOCATE ( omlhed( maxmll,         mxrhed), STAT=istat)
    ALLOCATE ( momlbd( maxmll, maxmlv, mxrbdf), STAT=istat)
    ALLOCATE ( momlhd( maxmll,         mxrhdf), STAT=istat)
    ALLOCATE ( yomlhd( maxmll                ), STAT=istat)
    ALLOCATE ( osgbdy( maxsgl,         mxsbdy), STAT=istat)
    ALLOCATE ( osghed( maxsgl,         mxshed), STAT=istat)
    ALLOCATE ( mosgbd( maxsgl,         mxsbdf), STAT=istat)
    ALLOCATE ( mosghd( maxsgl,         mxshdf), STAT=istat)
    ALLOCATE ( yosghd( maxsgl                ), STAT=istat)
    ALLOCATE ( ogpbdy( maxgpl,         mxgbdy), STAT=istat)
    ALLOCATE ( ogphed( maxgpl,         mxghed), STAT=istat)
    ALLOCATE ( mogpbd( maxgpl,         mxgbdf), STAT=istat)
    ALLOCATE ( mogphd( maxgpl,         mxghdf), STAT=istat)
    ALLOCATE ( yogphd( maxgpl                ), STAT=istat)

    omlbdy (:,:,:) = rmdi
    omlhed (:,:)   = rmdi
    momlbd (:,:,:) = imdi
    momlhd (:,:)   = imdi
    yomlhd (:)     = '        '
    osgbdy (:,:)   = rmdi
    osghed (:,:)   = rmdi
    mosgbd (:,:)   = imdi
    mosghd (:,:)   = imdi
    yosghd (:)     = '        '
    ogpbdy (:,:)   = rmdi
    ogphed (:,:)   = rmdi
    mogpbd (:,:)   = imdi
    mogphd (:,:)   = imdi
    yogphd (:)     = '        '

    ! Calculate log of pressure of the observation error levels
    ! ---------------------------------------------------------
    rolnlv(:) = LOG (rlevel(:))

    ! derive list of reference interpolation levels 'pref' from namelist input
    ! ------------------------------------------------------------------------
    IF ((PRESENT( nuout_in )) .AND. (PRESENT( yuout_in ))) THEN

      CALL obs_cdf_supob_levs ( nuout_in, yuout_in )
  !   =======================
    ELSE
      CALL obs_cdf_supob_levs
  !   =======================
    ENDIF

  ENDIF

!-------------------------------------------------------------------------------
!  Section  2.   Initialising buffer array for report rejection messages
!-------------------------------------------------------------------------------

! Find maximum length of output buffer
  istrej  = ilstidp
  nmxoln  = maxmll*(100+istrej) + (maxsgl+maxgpl)*(5+(istrej+1)/2)

  IF (ALLOCATED( outbuf )) THEN
    PRINT '("CAUTION in src_obs_cdfin_org: outbuf is already allocated "       &
          &,"at time / PE",F6.3,I5)', acthr, my_cart_id
    DEALLOCATE ( outbuf , STAT=istat )
  ENDIF
  ALLOCATE ( outbuf(nmxoln)             , STAT=istat )
  outbuf (:)   = 0

! Set initial number of output records to 0
  nacout       = 0

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_proc_init
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_proc_init


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for preparing radiosonde superobbing
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_supob_levs ( nuout_in, yuout_in )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" derives the list of
!   reference interpolation levels 'pref' from the namelist input for
!   the superobbing of high-resolution radiosonde profile data.
!
! Method:
!   Straightforward.
!   Input : arrays 'r_av_levs', 'r_av_incr', which are copied from namelist
!           variables 'av_levs', 'av_incr'.
!   Output: number of reference interpolation levels 'nref' and the
!           levels 'pref(1:nref)' themselves; 'pref' is allocated here.
!   (Further input (argument list) is for control output only.)
!
! Written by        :  Christoph Schraff, DWD  (original version: 22.04.19)
! Current Code Owner:  Christoph Schraff, DWD
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------

  INTEGER           , INTENT (IN)  , OPTIONAL  :: &
    nuout_in         ! unit number of output file

  CHARACTER (LEN=8) , INTENT (IN)  , OPTIONAL  :: &
    yuout_in         ! file name for unit number 'nuout_in'

! Local parameters: None
! ----------------

! Local variables:
! ----------------

  LOGICAL        ::  &
    lzopen           ! .true. if file unit 'nuout' is already open
                     !           before calling this routine

  INTEGER        ::  &
    navl          ,& ! length of non-empty part of level definition list
    jav           ,& ! index over interpolation level definition list
    kmax          ,& ! max. number of interpolation reference levels
    kk            ,& ! loop index
    itmp  , icl   ,& ! auxilliary buffer
    nuout            ! unit number of output file

  REAL (KIND=wp) ::  &
    pmin             ! limit for reference interpolation levels

  CHARACTER (LEN=8) ::  &
    yform            ! write format

  REAL (KIND=wp) , ALLOCATABLE ::  &
    ptmp   (:)    ,& ! temporary interpolation levels for superobbing
    ptln   (:)    ,& ! log( ptmp )
    pmean  (:)    ,& ! target pressure for the superobbing layer
    pmlnp  (:)       ! approx. pressure obtained by interpolation in ln(p)
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_supob_levs
!-------------------------------------------------------------------------------

    ! derive list of reference interpolation levels 'pref' from namelist input
    ! ------------------------------------------------------------------------
!!  av_levs  (1:7) = [ 1075., 755., 710., 610. , 587.5, 87.5, 60., 0. ]
!!  av_incr  (1:7) = [   10.,  15.,  20.,  22.5,  25. , 27.5, 20., 0. ]
!   av_levs_d(1:7) = [ 1075., 755., 710., 90., 75., 5. ]
!   av_incr_d(1:7) = [   10.,  15.,  20., 15., 10., 0. ]
!   dp_surf_d      = 200.   ! is not a namelist parameter currently
    navl    = COUNT( r_av_levs > 0._wp )
    kmax    = NINT( SUM( (r_av_levs(1:navl-1)-r_av_levs(2:navl))             &
                        / r_av_incr(1:navl-1)                  ) ) + navl
    ALLOCATE( ptmp (kmax) )
!   ptmp    = 0._wp
    ptmp(1) = r_av_levs(1)                             ! first avg.layer bound
    pmin    = 0.1_wp                                   ! Pa
    nref    = 1
    jav     = 1
    DO kk = 2, kmax
      IF (     (r_av_incr(jav) <= 0._wp)             & ! no increments left
          .OR. (ptmp(kk-1)-r_av_incr(jav) <= pmin))  & ! p>0 for ln(p) be real
        EXIT
      ptmp (kk)  =  ptmp(kk-1) - r_av_incr(jav)        ! next layer bound
      IF (ptmp(kk) <= r_av_levs(jav+1))  jav = jav + 1 ! take next increment
      nref  = kk
    ENDDO
    IF (nref > 1) THEN
      ! assuming that the target pressure for the superobbing layers is the mean
      ! of the pressure between the layer bounds as given by the namelist:
      ! apply a correction to the layer bounds such that averaging in log(p)
      ! within the layers results in (approx.) the target pressure
      ! ------------------------------------------------------------------------
      ALLOCATE( pmean (nref-1) )
      ALLOCATE( pmlnp (nref-1) )
      ALLOCATE( ptln  (nref) )
      DO kk = 1 , nref
        ptln (kk)  =  LOG( ptmp(kk) )
      ENDDO
      DO kk = 1 , nref-1
        !   this is the target pressure for the superobbing layer
        pmean (kk)  =       c05*( ptmp(kk) + ptmp(kk+1) )
        !   this is the approx. pressure obtained by interpolation in ln(p)
        pmlnp (kk)  =  EXP( c05*( ptln(kk) + ptln(kk+1) ) )
      ENDDO
      !   apply the correction
      DO kk = 2 , nref-1
        ptmp (kk)  =  ptmp(kk) + c05*(  pmean(kk) + pmean(kk-1)                &
                                      - pmlnp(kk) - pmlnp(kk-1))
      ENDDO
      ptmp (nref)  =  ptmp(nref) + c05*(pmean(nref-1) - pmlnp(nref-1))
      ptmp (1)     =  ptmp(1)    + c05*(pmean(     1) - pmlnp(     1))
      DEALLOCATE( pmean )
      DEALLOCATE( pmlnp )
      DEALLOCATE( ptln  )
    ENDIF
    ! 'pref' has save atttribute, is allocated once and never de-allocated
    ALLOCATE( pref (nref) )
    pref (1:nref) = ptmp(1:nref)
    DEALLOCATE( ptmp )

    ! writing superobbing reference levels as control output
    ! ------------------------------------------------------
    IF (my_cart_id == 0) THEN
      ! write to file as indicated by optional subroutine arguments if given,
      ! otherwise write to standard output
      nuout  = 0
      lzopen = .TRUE.
      IF ((PRESENT( nuout_in )) .AND. (PRESENT( yuout_in )))  nuout = nuout_in
      IF (nuout /= 0) THEN
        INQUIRE( UNIT=nuout, OPENED=lzopen )
        IF (.NOT. lzopen) OPEN( nuout, FILE=yuout_in, FORM='FORMATTED'         &
                                     , STATUS='UNKNOWN', POSITION='APPEND')
      ENDIF
      WRITE( nuout,'(" ")' )
      WRITE( nuout,'(" Layer bounds for superobbing of high-res radiosondes:", &
                    &I4,", as follows (in hPa):")') nref
      itmp = nref / 10
      DO kk = 1 , itmp
        WRITE( nuout,'(10F8.1)' ) (pref(icl)*0.01_wp, icl=10*(kk-1)+1,10*kk)
      ENDDO
      itmp = MOD( nref, 10 )
      IF (itmp > 0) THEN
        WRITE( yform,'("(",I2.2,"F8.1)")' ) itmp
        WRITE( nuout,yform ) (pref(icl)*0.01_wp, icl=nref-itmp+1,nref)
      ENDIF
      IF (.NOT. lzopen) CLOSE( nuout )
    ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_supob_levs
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_supob_levs


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdfin_open_files  ( icdfdirlen , ycdfdir , ydace_op )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" opens all the
!   NetCDF observation input files.
!   The corresponding file unit number is also used as indicator of
!   file existence and is distributed to all nodes.
!
! Method:
!   The file names given by 'ycdfin' in 'data_obs_cdfin.f90' can used as is or
!   with an annex '.nc'. For any of the observation file types, none, one, or
!   several files may be present. If there are several files, the annex '.2',
!   '.3', etc. is used (in front of the annex '.nc').
!   Here, the full path names of the files are created and the existence of
!   the files are checked first. Existing files are opened.
!   Output values of this routine are the number of existing files, the obs type
!   of the files, the file unit numbers, and the annex for file names. This
!   information is distributed to all nodes.
!   Note that the model DOES NOT ABORT if the file does not exist.
!   It does abort only if opening of an existing file fails.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER            , INTENT (IN)  ::  &
    icdfdirlen       ! max. length of name of directory where
                     !   NetCDF observation input files reside

  CHARACTER (LEN=*)  , INTENT (IN)  ::  &
    ycdfdir          ! directory where the NetCDF observation input files reside

  CHARACTER (LEN= *) , INTENT (IN)  , OPTIONAL  ::  &
    ydace_op (:)   ! obs types processed by dace obs operators

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL        ::  &
    lexicdf       ,& ! indicates whether NetCDF observation input file exists
    lproc (mxobtp+1) ! links obs types to choice of obs operators (COSMO/dace)

  INTEGER        ::  &
    kcdf          ,& ! obs file type of NetCDF observation input file
    ncid          ,& ! file unit number of NetCDF observation input file
    iannex        ,& ! number of NetCDF obs input file of current obs file type
    status        ,& ! NetCDF status variable
    ilenpath      ,& ! length of file name / path name of NetCDF obs input file
    ierr  , istat ,& ! error status variable
    iot_cdfin (ntype_cdfin)  ! links file units to obs types

  CHARACTER (LEN=MIN( icdfdirlen+icdfinlen, 255 ))  ::  &
    yncfile       ,& ! total path of NetCDF observation input file, w/out annex
    yncfilea         ! total path of NetCDF observation input file, with annex
  CHARACTER (LEN=icdfinlen+iannexlen)  ::  &
    yfn              ! file name of NetCDF observation input file, with annex
  CHARACTER (LEN=iannexlen)  ::  &
    yannex           ! total path of NetCDF observation input file
  CHARACTER (LEN=20)       ::  &
    yroutine         ! name of this subroutine
  CHARACTER (LEN=60)       ::  &
    yerrmsg          ! error message
  CHARACTER (LEN=30)       ::  &
    yerr             ! error message
  CHARACTER (LEN=100)        ::  &
    ytmp (mxcdfin)   ! temporal file for message passing with distribute_values
                     !  (which expects characters with length LEN=100 !)

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdfin_open_files
!-------------------------------------------------------------------------------

  yroutine = 'obs_cdfin_open_files'

  n_cdfin  = 0
  icdfin   (:)  =  0
  yncannex (:)  =  '      '

  ierr  =  0
  yncfile = ' '
  IF (my_cart_id == 0) THEN

    ! prepare check which 'cdfin' files are not opened here
    ! because they are processed by the dace obs operators
    ! ------------------------------------------------------
    ! 1. link file units to obs types
    !    (ensure that iot_cdfin(:) is well defined for all indices
    !     with values >= 1 and <= mxobtp+1 = size of lproc)
    iot_cdfin (:              ) = mxobtp + 1
    iot_cdfin (ncdf_temp      ) = ntemp
    iot_cdfin (ncdf_tempship  ) = ntemp
    iot_cdfin (ncdf_temphirs  ) = ntemp
!   iot_cdfin (ncdf_tempdrop  ) = ntemp
    iot_cdfin (ncdf_tempdesc  ) = ntemp
    iot_cdfin (ncdf_pilot     ) = npilot
    iot_cdfin (ncdf_pilot_p   ) = npilot
    iot_cdfin (ncdf_amdar_ml  ) = nairep
!   iot_cdfin (ncdf_amdar_vp  ) = nairep
    iot_cdfin (ncdf_amdar     ) = nairep
    iot_cdfin (ncdf_acars     ) = nairep
    iot_cdfin (ncdf_acars_uk  ) = nairep
    iot_cdfin (ncdf_acars_us  ) = nairep
    iot_cdfin (ncdf_modes     ) = nairep
    iot_cdfin (ncdf_modes_acr ) = nairep
    iot_cdfin (ncdf_wprof     ) = npilot
    iot_cdfin (ncdf_rass      ) = npilot
    iot_cdfin (ncdf_radar_vad ) = npilot
    iot_cdfin (ncdf_synop     ) = nsynop
    iot_cdfin (ncdf_synop_mob ) = nsynop
    iot_cdfin (ncdf_ship      ) = nsynop
    iot_cdfin (ncdf_buoy      ) = ndribu
    iot_cdfin (ncdf_metar     ) = nsynop
    iot_cdfin (ncdf_gps_zenith) = ngps
    iot_cdfin (ncdf_ascat     ) = nscatt
    iot_cdfin (ncdf_qscat     ) = nscatt
    iot_cdfin (ncdf_satob     ) = nsatob
    iot_cdfin (ncdf_nat_dwd   ) = nsynop
!   iot_cdfin (ncdf_gnss_ztd  ) = ngps
!   iot_cdfin (ncdf_gnss_std  ) = ngps
    iot_cdfin (ncdf_tower     ) = npilot
    iot_cdfin (ncdf_tower_icos) = npilot
    iot_cdfin (ncdf_wlidar_wp ) = npilot
    iot_cdfin (ncdf_synop_tst ) = nsynop
    ! 2. link obs types to choice of obs operators (COSMO or dace)
    lproc (:)  =  .TRUE.
    IF (PRESENT( ydace_op )) THEN
      lproc (nsynop) = (.NOT. ANY( ydace_op == 'SYNOP' ))
      lproc (nairep) = (.NOT. ANY( ydace_op == 'AIREP' ))
      lproc (nsatob) = (.NOT. ANY( ydace_op == 'SATOB' ))
      lproc (ndribu) = (.NOT. ANY( ydace_op == 'DRIBU' ))
      lproc (ntemp ) = (.NOT. ANY( ydace_op == 'TEMP'  ))
      lproc (npilot) = (.NOT. ANY( ydace_op == 'PILOT' ))
      lproc (nscatt) = (.NOT. ANY( ydace_op == 'SCATT' ))
      lproc (ngps  ) = (.NOT. ANY( ydace_op == 'GPS'   ))
    ENDIF

! loop over types of NetCDF observation input files
! -------------------------------------------------
    DO kcdf = 1 , ntype_cdfin
      iannex  = 0
      lexicdf = .TRUE.

      ! get full path of NetCDF observation input file
      yncfile  = ycdfin(kcdf)
      ilenpath = LEN_TRIM(ycdfdir)
      IF (ilenpath >= 1) THEN
        IF (ycdfdir(ilenpath:ilenpath) == '/')  ilenpath = ilenpath - 1
        yncfile = ycdfdir(1:ilenpath) // '/' // yncfile(1:icdfinlen)
      ENDIF
      ilenpath = LEN_TRIM(yncfile)

! loop: search for file of the current obs file type
! --------------------------------------------------
      DO WHILE (lexicdf)
        iannex = iannex + 1

! check existence of NetCDF file (without or with '.nc' annex)
! ------------------------------
        IF (iannex == 1) THEN
          WRITE( yannex, '(A1,I1,A4)' )  '.', iannex, '    '
          yncfilea = yncfile(1:ilenpath) // yannex
          INQUIRE (FILE=yncfilea, EXIST=lexicdf)
          IF (.NOT. lexicdf) yannex = '      '
        ELSEIF (iannex <= 9) THEN
          WRITE( yannex, '(A1,I1,A4)' )  '.', iannex, '    '
        ELSEIF (iannex <= 99) THEN
          WRITE( yannex, '(A1,I2,A3)' )  '.', iannex, '   '
        ENDIF
        yncfilea = yncfile(1:ilenpath) // yannex
        INQUIRE (FILE=yncfilea, EXIST=lexicdf)

        IF (.NOT. lexicdf) THEN
          yannex = yannex(1:LEN_TRIM(yannex)) // '.nc'
          yncfilea = yncfile(1:ilenpath) // yannex
          INQUIRE (FILE=yncfilea, EXIST=lexicdf)
        ENDIF
        IF ((iannex == 1) .AND. (.NOT. lexicdf))                               &
          PRINT '(3A)' , '  NOTE: NO FILE ', yncfile(1:ilenpath), ' (.nc)'

! check if dace obs operator to be used for this obs type
! ---------------------------------------------------------
        IF ((lexicdf) .AND. (.NOT. lproc(iot_cdfin(kcdf)))) THEN
          yfn     = ycdfin(kcdf) (1:LEN_TRIM(ycdfin(kcdf))) // yannex
          PRINT          '("NOTICE ! ",A ," not opened by COSMO routines"      &
                         &," because processed by dace operators")' ,  yfn

! open NetCDF file (if it exists)
! ----------------
        ELSEIF ((lexicdf) .AND. (n_cdfin < mxcdfin)) THEN
          status   =   nf90_open ( yncfilea, NF90_NOWRITE, ncid )

          IF (status /= nf90_noerr) THEN
            PRINT '(2A)' , 'ERROR OPENING FILE ', yncfilea
            yfn     = ycdfin(kcdf) (1:LEN_TRIM(ycdfin(kcdf))) // yannex
            yerrmsg = 'OPENING OF ' // yfn(1:LEN_TRIM(yfn)) // ' FAILED'
            CALL model_abort (my_cart_id, 1025, yerrmsg, yroutine)
          ENDIF

! update file counter, obs type of files,
! file unit numbers, and annex for file names
! -------------------------------------------
          n_cdfin  =  n_cdfin + 1
          icdfin   (n_cdfin)  =  kcdf
          ncinid   (n_cdfin)  =  ncid
          yncannex (n_cdfin)  =  yannex

        ELSEIF (lexicdf) THEN

! write message if file exists and (n_cdfin >= mxcdfin)
! -----------------------------------------------------
          OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'      &
                                     , POSITION='APPEND', IOSTAT=istat)
          IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
          IF (istat /= 0) CALL model_abort (my_cart_id, 1419, yerr, yroutine)
          yfn     = ycdfin(kcdf) (1:LEN_TRIM(ycdfin(kcdf))) // yannex
          WRITE( nucautn,'("CAUTION !!!!! ",A ," not opened because # files >" &
                         &," mxcdfin =",I4)' )  yfn, mxcdfin
          WRITE( nucautn,'("   -->  decrease # files or increase hard-coded "  &
                         &,"parameter mxcdfin =",I4)' )  mxcdfin
          PRINT          '("CAUTION !!!!! ",A ," not opened because # files >" &
                         &," mxcdfin =",I4)' ,  yfn, mxcdfin
          PRINT          '("   -->  decrease # files or increase hard-coded "  &
                         &,"parameter mxcdfin =",I4)',  mxcdfin
          CLOSE (nucautn)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

! distribute unit number of NetCDF input file   (unit number = 0 indicates
! -------------------------------------------    that file does not exist)
  IF (num_compute > 1)                                                         &

    CALL distribute_values ( n_cdfin, 1, 0, imp_integers, icomm_cart, ierr )
!   ======================

  IF ((num_compute > 1) .AND. (n_cdfin >= 1)) THEN

    CALL distribute_values ( icdfin  , n_cdfin, 0,imp_integers ,icomm_cart,ierr)
!   ======================
    CALL distribute_values ( ncinid  , n_cdfin, 0,imp_integers ,icomm_cart,ierr)
!   ======================

    ! careful with that axe Eugene !
    ! arrays of characters in distribute_values must have (LEN=100)
    ytmp             = ' '
    ytmp (1:n_cdfin) = yncannex(1:n_cdfin)

    CALL distribute_values ( ytmp, n_cdfin , 0, imp_character, icomm_cart, ierr)
!   ======================
    yncannex (1:n_cdfin) = ytmp(1:n_cdfin) (1:iannexlen)

  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_cdfin_open_files
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdfin_open_files


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_fof_check_files ( icdfdirlen , ycdfdir )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" checks the existence of
!   one or several NetCDF feedobs (or feedback) files as observational input.
!   The corresponding file unit number is also used as indicator of
!   file existence and is distributed to all nodes.
!   Note that in contrast to the 'cdfin' NetCDF observation input files in
!   routine 'obs_cdfin_open_files', the feedobs files are not opened here.
!   Instead, they will be opened in the call of 'read_fdbk_head'.
!
! Method:
!   The file name 'fof' (given by 'yfofin' in 'data_obs_cdfin.f90') can used
!   as is or with an annex '.nc'. None, one, or several files may be present.
!   If there are several files, the annex '.2', '.3', etc. is used (in front of
!   the annex '.nc').
!   Here, the full path names of the files are created and the existence of
!   the files are checked first. Existing files are not opened here.
!   Output values of this routine are the number of existing files, and the
!   annex for file names. This information is distributed to all nodes.
!   Note that the model DOES NOT ABORT if the file does not exist.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)  ::  &
    icdfdirlen       ! max. length of name of directory where
                     !   NetCDF observation input files reside

  CHARACTER (LEN=*)  , INTENT (IN)  ::  &
    ycdfdir          ! directory where the NetCDF observation input files reside

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  LOGICAL        ::  &
    lexifof          ! indicates whether NetCDF observation input file exists

  INTEGER        ::  &
    iannex        ,& ! number of NetCDF obs input file of current obs file type
    ilenpath      ,& ! length of file name / path name of NetCDF obs input file
    ierr   , istat   ! error status variable

  CHARACTER (LEN=MIN( icdfdirlen+icdfinlen, 255 ))  ::  &
    yfofpath      ,& ! total path of NetCDF input feedobs file, w/out annex
    yfofpatha        ! total path of NetCDF input feedobs file, with annex
  CHARACTER (LEN=icdfinlen+iannexlen)  ::  &
    yfn              ! file name of NetCDF observation input file, with annex
  CHARACTER (LEN=iannexlen)  ::  &
    yannex           ! total path of NetCDF observation input file
  CHARACTER (LEN=20)       ::  &
    yroutine         ! name of this subroutine
  CHARACTER (LEN=30)       ::  &
    yerr             ! error message
  CHARACTER (LEN=100)        ::  &
    ytmp (mxcdfin)   ! temporal file for message passing with distribute_values
                     !  (which expects characters with length LEN=100 !)

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_fof_check_files
!-------------------------------------------------------------------------------

  yroutine = 'obs_fof_check_files'

  n_fofin  = 0
  yfofannex (:)  =  '      '
  ierr  =  0

  ! get full path of input feedobs (feedback) file(s)
  yfofpath = yfofin
  ilenpath = LEN_TRIM(ycdfdir)
  IF (ilenpath >= 1) THEN
    IF (ycdfdir(ilenpath:ilenpath) == '/')  ilenpath = ilenpath - 1
    yfofpath = ycdfdir(1:ilenpath) // '/' // yfofpath(1:icdfinlen)
  ENDIF
  ilenpath = LEN_TRIM(yfofpath)

  IF (my_cart_id == 0) THEN
    iannex  = 0
    lexifof = .TRUE.

! loop: search for feedobs files
! ------------------------------
    DO WHILE (lexifof)
      iannex = iannex + 1

! check existence of NetCDF file (without or with '.nc' annex)
! ------------------------------
      IF (iannex == 1) THEN
        yannex = '      '
      ELSEIF (iannex <= 9) THEN
        WRITE( yannex, '(A1,I1,A4)' )  '.', iannex, '    '
      ELSEIF (iannex <= 99) THEN
        WRITE( yannex, '(A1,I2,A3)' )  '.', iannex, '   '
      ENDIF
      yfofpatha = yfofpath(1:ilenpath) // yannex
      INQUIRE (FILE=yfofpatha, EXIST=lexifof)

      IF (.NOT. lexifof) THEN
        yannex = yannex(1:LEN_TRIM(yannex)) // '.nc'
        yfofpatha = yfofpath(1:ilenpath) // yannex
        INQUIRE (FILE=yfofpatha, EXIST=lexifof)
      ENDIF
      IF ((iannex == 1) .AND. (.NOT. lexifof))                                 &
        PRINT '(3A)' , '  NOTE: NO FILE ', yfofpath(1:ilenpath), ' (.nc)'

! update file counter, obs type of files,
! file unit numbers, and annex for file names
! -------------------------------------------
      IF ((lexifof) .AND. (n_fofin <= mxfofin)) THEN
        n_fofin  =  n_fofin + 1
        yfofannex (n_fofin)  =  yannex

      ELSEIF (lexifof) THEN

! write message if file exists and (n_fofin >= mxfofin)
! -----------------------------------------------------
        OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'      &
                                   , POSITION='APPEND', IOSTAT=istat)
        IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
        IF (istat /= 0) CALL model_abort (my_cart_id, 1420, yerr, yroutine)
        yfn     = yfofin (1:LEN_TRIM(yfofin)) // yannex
        WRITE( nucautn,'("CAUTION !!!!! ",A ," not opened because # files >" &
                       &," mxfofin =",I4)' )  yfn, mxfofin
        WRITE( nucautn,'("   -->  decrease # files or increase hard-coded "  &
                       &,"parameter mxfofin =",I4)' )  mxfofin
        PRINT          '("CAUTION !!!!! ",A ," not opened because # files >" &
                       &," mxfofin =",I4)' ,  yfn, mxfofin
        PRINT          '("   -->  decrease # files or increase hard-coded "  &
                       &,"parameter mxfofin =",I4)',  mxfofin
        CLOSE (nucautn)
      ENDIF
    ENDDO
  ENDIF

! distribute unit number of NetCDF input file   (unit number = 0 indicates
! -------------------------------------------    that file does not exist)
  IF (num_compute > 1)                                                         &

    CALL distribute_values ( n_fofin, 1, 0, imp_integers, icomm_cart, ierr )
!   ======================

  IF ((num_compute > 1) .AND. (n_fofin >= 1)) THEN

    ! careful with that axe Eugene !
    ! arrays of characters in distribute_values must have (LEN=100)
    ytmp             = ' '
    ytmp (1:n_fofin) = yncannex(1:n_fofin)

    CALL distribute_values ( ytmp, n_fofin , 0, imp_character, icomm_cart, ierr)
!   ======================
    yncannex (1:n_fofin) = ytmp(1:n_fofin) (1:iannexlen)

  ENDIF

! set losse_fg = .false. if there are no 'fof' input files
  IF (n_fofin == 0)  losse_fg = .FALSE.

!-------------------------------------------------------------------------------
! End Subroutine obs_fof_check_files
!-------------------------------------------------------------------------------
END SUBROUTINE obs_fof_check_files


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for flagging redundant reports/data
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_redundancy ( nnmlnew, nnsgnew, nngpnew, ntotsgo, maxmlv )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" flags redundant reports
!   resp. data.
!
! Method:
!   Separate treatment for redundancy checking for each type of ODR arrays
!   (i.e. between single-level, multi-level, resp. GPS-IWV reports.
!   Each of those reports which have been read from file and put to the ODR
!   immediately before calling this routine are checked against all the other
!   reports (which have a smaller report index value in the ODR). If several
!   reports meet the redundancy criteria (such as 4-D quasi-colocation) for
!   a report to be checked then only one (preferably with the same status
!   (active / passive)) is chosen.
!   If appropriate, the active report is complemented by data of the passive
!   report, e.g. in order to merge TEMP A, B, C, and D parts into one report
!   or to merge separate temperature and wind profiles.
!   In addition, multi-level reports are first 'auto-checked' against themselves
!   in order to merge or discard quasi-colocated vertical levels.
!   Redundant reports are flagged, but not discarded unless a redundant
!   multi-level report is fully contained in another report.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! (Some changes in previous routine 'obs_redundancy':
!  1.31 1999/07/01 Christoph Schraff: Bug corr. at merging multi-level reps:
!                                     - above top of active report
!                                     - when merging T-profile with wind profile
!                                     - removal of approx. colocated levels.
!                                     Complementation with mandatory levels
!                                     from rejected reports.
!  1.36 2000/02/24 Michael Buchhold : Introduction of ACAR aircraft reports.
!  2.5  2001/06/01 Christoph Schraff: Correction to redundancy checking between
!                                     more than 2 colocated reports.
!  3.6  2003/12/11 Christoph Schraff: Inclusion of observation type SATOB.
!-------------------------------------------------------------------------------

#ifndef __COSMO__
  !
  ! acces to DACE modules
  !
  USE data_obs_lib_cosmo, ONLY: gen_grid  ! generic grid meta date (for ICON)
  USE mo_physics,         ONLY: r2d       ! factor rad -> degree
#endif

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER        , INTENT (IN)     ::  &
    nnmlnew       ,& ! number of single-level  \  reports which are new and have
    nnsgnew       ,& ! number of  multi-level   > to be checked for redundancy
    nngpnew       ,& ! number of  GPS  (IWV)   /  (against all other reports)
    ntotsgo       ,& ! number of 'old' single-level reports (i.e. before having
                     !   read new data at the current timestep)
    maxmlv           ! max. number of vertical levels in multi-level report

! Local parameters
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c0e = 0.0001_wp  ! small constant less than accuracy of values

! Local scalars
! -------------

  LOGICAL        ::  &
    lairc  , lair2   ,& ! .TRUE., if aircraft report
    lsatbc , lsatb2  ,& ! .TRUE., if satob report
    lship  , lship2  ,& ! .TRUE., if ship report
    lcoloc           ,& ! .TRUE. for spatial colocation
    lraso  , lraso2  ,& ! .TRUE., if obs type is radiosonde (TEMP or true PILOT)
    lwprofl, lwprof2 ,& ! .TRUE., if obs type is wind profiler, RASS or radarVAD
    ltower , ltower2 ,& ! .TRUE., if obs type is tower (or ICOS tower)
    lnatdwd, lnatdwd2,& ! .TRUE., if obs type is DWD national surface report
    lrbufr1, lrbufr2 ,& ! .TRUE., if radiosonde BUFR report
    lsupoba, lsupobr ,& ! .TRUE., if superobbed report
    ltop             ,& ! .TRUE., if top level
    lrepla           ,& ! .TRUE., if data from the redundant report may replace
                        !            missing data of the active report
    lreplanl         ,& ! .TRUE., if lrepla unless temporal non-local variable
                        !            and active report not valid at full hour
    lmand            ,& ! .TRUE., if rejected level is mandatory level
    leqsfc           ,& ! .TRUE., if both levels are surface or are not surface
    lactsurf         ,& ! .TRUE., if level from active report is surface level
    lmerge           ,& ! .TRUE., if active and redundant reports are two TEMP
                        !            (PILOT) parts to be merged together
    lsupoac          ,& ! .TRUE., if active layer is superobs
    lsuporj          ,& ! .TRUE., if rejected layer is superobs
    lsgexch          ,& ! .TRUE., if level significance of superob level is to
                        !            be replaced (by lev. sign. of rejected l.)
    lsub             ,& ! .TRUE., if current report is contained in active rep.
                        !            (and hence is discarded)
    laddu, laddt, laddq ! .TRUE., if current quantity is to be added

  INTEGER        ::  &
    nsgob  , nmlob   ,& !  \  ODR outer loop indices over single-level / multi-
    ngpob            ,& !  /  level / GPS reports to be checked for redundancy
    nsglo  , nmllo   ,& !  \  ODR inner loop indices over single-level / multi-
    ngplo            ,& !  /  level / GPS reports to be checked against
    nsgo2  , nmlo2   ,& !  \  ODR indices of single-level / multi-level / GPS
    ngpo2            ,& !  /  report which meets redundancy criteria
    nmlza            ,& ! ODR index used to store the (active) temporary array
    nmlrj            ,& ! ODR index used to store the rejected report
    nchkend          ,& ! number of reports to check against the current report
    nodract          ,& ! index of active ODR report
    nodrrej          ,& ! index of rejected ODR report
    iactio           ,& ! =1 : redundancy check , =2 : merge TEMP parts A,B,C,D
    ktypchk          ,& ! =1 : auto checking    , =2 : ordinary redundancy check
    icma             ,& ! pointer index for 'cma' / diagnostic arrays
    kobtyp , kobtyp2 ,& ! observation type
    kcdtyp , kcdtyp2 ,& ! code type
    nstcor           ,& ! station correction bit
    npassiv, npassv2 ,& ! passive (bit) flag
    ksubcen, ksubce2 ,& ! data (sub-)centre which has processed the GPS data
    igplon , igplon2 ,& ! gridpoint in zonal direction assigned to observ.
    igplat , igplat2 ,& ! gridpoint in meridional direction assig. to observ.
    nreplpr          ,& ! number of replaced data
    nactpr              ! active report

  INTEGER        ::  &
    kobtypr, kcdtypr ,& ! observation/code type of rejected report
    klev, krej, kact ,& ! level indices
    krjc, kclo, klvi ,& ! level indices
    nflaga , nflagr  ,& ! variable flags
    ilen             ,& ! length of output record
    ilvsg            ,& ! level significance (level identity)
    isdel            ,& ! levels below surface level to be discarded
    nrejlv           ,& ! number of levels to be removed
!   kz1    , kz2     ,& ! DWD internal classification number
    nlev   , nlv2    ,& ! number of vertical levels
    nlvcom           ,& ! number of common pressure levels in 2 reports
    nactgpc          ,& ! number of actively used processing centers in i_gpscen
    ixcen  , ixce2   ,& ! index of processing center in i_gpscen
    iiv    , icl     ,& ! loop indices
    iis    , iie     ,& ! start / end of characters of station id to be checked
    istat               ! error status variable

  REAL (KIND=wp) ::  &
    rdegkm2          ,& ! factor for conversion of degree to km (**2)
    rscale           ,& ! scaling factor
    obtime , obtime2 ,& ! observation date/time
    obbalt , obbalt2 ,& ! observation altitude
    zptop1 , zptop2  ,& ! pressure of uppermost superobbed obs level
    zpbase1, zpbase2 ,& ! pressure of lowermost superobbed obs level
    rtmlmt , rtmlmte ,& ! colocation threshold for time
    rhzlmt , rhzlmte ,& ! colocation threshold for horizontal distance
    rvtlmt , rvtlmte ,& ! colocation threshold for vertical distance
    rdplmt , rprlmt  ,& ! colocation thresholds for vert. dist. in pres. units
    rdplmt2          ,& ! colocation threshold (temporary)
    dprjc  , dprjca  ,& ! pressure differences
    pnexlmt          ,& ! pressure limit imposed by next active level above
    rhzob            ,& ! horizontal distance in km
    epsy2               ! 2 * epsy

  CHARACTER (LEN=ilstid)   :: &
    ystid , ystidz      ! station identity

! Local (automatic) arrays:
! -------------------------

  REAL (KIND=wp) ::  &
    zmlbdy (2*maxmlv,mxrbdy),& ! temporary real array for multi-level rep. body
                               ! (array length '2*maxmlv' prevents from writing
                               !  levels beyond array limits)
    zmlhed          (mxrhed),& ! temporary array for one multi-level rep. header
    zsgbdy          (mxsbdy),& ! temporary real array for single-level rep. body
    zsghed          (mxshed)   ! temporary array for a single-level rep. header

  INTEGER        ::  &
    iomlbd (2*maxmlv,mxrbdf),& ! temporary int array for multi-level report body
    iomlhd          (mxrhdf),& ! temporary array for one multi-level rep. header
    mzsgbd          (mxsbdf),& ! temporary int array for sing.-level report body
    mzsghd          (mxshdf)   ! temporary array for one sing.-level rep. header

! Local arrays:
! ------------

  REAL (KIND=wp) , ALLOCATABLE :: &
    cosolat2     (:)    ! ( cos(lat(obs)) ) **2

  LOGICAL        , ALLOCATABLE :: &
    lredn        (:)    ! report meet redundancy criteria with current report
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_redundancy
!-------------------------------------------------------------------------------

! preset some variables
! ---------------------

  rdegkm2 = 111._wp **2
  rscale  = 100._wp*(1._wp+epsy)
  epsy2   = 2* epsy

#ifdef __COSMO__
  ALLOCATE ( cosolat2 (MAX( nnsgnew, nnmlnew, nngpnew, 1 )) , STAT=istat )
#endif
  ALLOCATE ( lredn    (MAX( ntotsg , ntotml , ntotgp , 1 )) , STAT=istat )

!-------------------------------------------------------------------------------
! Section 1.0  Redundancy check for single-level reports
!-------------------------------------------------------------------------------

#ifdef __COSMO__
  DO nsgob = 1 , nnsgnew
    igplat = mosghd (nsgob+ntotsg-nnsgnew,nhjtot)
    cosolat2 (nsgob) = COS( (r_startlat_tot+(igplat-1)*r_dlat) * r_degrad )
    cosolat2 (nsgob) = cosolat2(nsgob) * cosolat2(nsgob)
  ENDDO
#endif

  !  check limits for radiosonde-derived surface level report vs. Synop
  rtmlmte = rtmlrsy
  rhzlmte = rhzlshp
  rvtlmte = rvtlsy
  rhzlmte = rhzlmte * rhzlmte

DO nsgob = ntotsg - nnsgnew + 1 , ntotsg

! get some information of current report
! --------------------------------------
!  Observation and code type

  kobtyp   =  mosghd (nsgob,nhobtp)
  kcdtyp   =  mosghd (nsgob,nhcode)
! lship    =       (kcdtyp == nshscd) .OR. (kcdtyp == nabscd)                  &
!             .OR. (kcdtyp == nshred) .OR. (kcdtyp == natshs)
  lship    =       (kcdtyp == nshscd) .OR. (kcdtyp == natshs)
  lnatdwd  =       (mosghd(nsgob,nhkz) == 40005) .OR.(mosghd(nsgob,nhkz) == 5) &
              .OR. (mosghd(nsgob,nhkz) == 10005)
  lraso    =  (kobtyp == ntemp) .OR. (kobtyp == npilot)
  lairc    =  (kobtyp == nairep)
  lsatbc   =  (kobtyp == nsatob)
  ltower   =  (kcdtyp == ntower) .OR. (kcdtyp == ntowic)

!  Get WMO station identity, station correction bit, and passive flag

  ystid    =  yosghd (nsgob)
  nstcor   =  mosghd(nsgob,nhcorr)
  npassiv  =  mosghd(nsgob,nhpass)
  IF ((npassiv == 2) .AND. (IBCLR( mosghd(nsgob,nhflag), FL_OBSTYPE ) > 0))    &
    npassiv = 3

!  Get 4-d location, and cos(latitude) for distances

  obtime = osghed (nsgob,nhtime)
  igplon = mosghd (nsgob,nhitot)
  igplat = mosghd (nsgob,nhjtot)
  obbalt = osghed (nsgob,nhalt )
  IF (lairc .OR. lsatbc) obbalt = - osgbdy (nsgob,nbsplv)

!  Preset thresholds for 4-d colocation

  IF ((kobtyp == nairep) .OR. (kobtyp == nsatob)) THEN
    rtmlmt = rtmlair
    rhzlmt = rhzlair
    rvtlmt = rvtlair * 100._wp
  ELSEIF ((kobtyp == ndribu) .OR. (lship)) THEN
    rtmlmt = rtmlim
    !  enlarged check limits to account for BUFR and TAC (alphanumeric) reports
    rhzlmt = rhzlshp
    rvtlmt = rvtlshp
  ELSEIF (kobtyp == nsynop) THEN
    !  enlarged check limits to account for BUFR and TAC (alphanumeric) reports
    rtmlmt = rtmlsy
    rhzlmt = rhzlshp
    rvtlmt = rvtlsy
  ELSE
    rtmlmt = rtmlim
    rhzlmt = rhzlim
    rvtlmt = rvtlim
    IF (lraso)  rtmlmt = rtmlrs
    IF (ltower) rtmlmt = rtmltow
  ENDIF
  rhzlmt = rhzlmt * rhzlmt

! check if current report meets the redundancy criteria with another report
! -------------------------------------------------------------------------
! (do this with a vectorisable loop)

  nchkend = nsgob - 1
  IF (npassiv == -1) nchkend = 0

  DO nsglo = 1 , nchkend
    lredn (nsglo)  = .FALSE.

!   Get 4-d location and pressure of second report
    obtime2 = osghed (nsglo,nhtime)
    igplon2 = mosghd (nsglo,nhitot)
    igplat2 = mosghd (nsglo,nhjtot)
    obbalt2 = osghed (nsglo,nhalt )
    kobtyp2 = mosghd (nsglo,nhobtp)
    lair2    =  (kobtyp2 == nairep)
    lsatb2   =  (kobtyp2 == nsatob)
    lraso2   =  (kobtyp2 == ntemp ) .OR. (kobtyp2 == npilot)
!   ltower2  =  (kcdtyp2 == ntower) .OR. (kcdtyp2 == ntowic)
    IF (lair2 .OR. lsatb2) obbalt2 = - osgbdy (nsglo,nbsplv)
    lnatdwd2 =     (mosghd(nsglo,nhkz) == 40005) .OR.(mosghd(nsglo,nhkz) == 5) &
              .OR. (mosghd(nsglo,nhkz) == 10005)

!   Check for horizontal and vertical (quasi-) colocation and simultaneity

#ifdef __COSMO__
    rhzob  =    cosolat2(nsgob-(ntotsg-nnsgnew))                               &
               *((igplon2-igplon) *r_dlon) *((igplon2-igplon) *r_dlon)         &
             +  ((igplat2-igplat) *r_dlat) *((igplat2-igplat) *r_dlat)
#else
    rhzob  = r2d **2 * SUM ((  gen_grid% xnglob (igplon ,igplat ,:,1)          &
                             - gen_grid% xnglob (igplon2,igplat2,:,1)) **2)
#endif
    IF (     ((lraso ) .AND. (kobtyp2 == nsynop))                              &
        .OR. ((lraso2) .AND. (kobtyp  == nsynop))) THEN
      lcoloc =      (      rhzob * rdegkm2     <=  rhzlmte )                   &
              .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmte )                   &
              .AND. ( ABS( obtime2 - obtime )  <=  rtmlmte )
    ELSE
      lcoloc =      (      rhzob * rdegkm2     <=  rhzlmt )                    &
              .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmt )                    &
              .AND. ( ABS( obtime2 - obtime )  <=  rtmlmt )
    ENDIF
    ! reproducibility irrespective of domain decomposition ensured by allowing
    ! for redundancy only between reports assigned to the same grid point
    IF (lredn_repro)  lcoloc = lcoloc .AND. (igplon2 == igplon)                &
                                      .AND. (igplat2 == igplat)
    IF ((lcoloc) .AND. (mosghd(nsglo,nhpass) /= -1)) THEN

!     check for station / report ID

      lredn (nsglo)  =  (yosghd(nsglo) == ystid)
      !   initial 2 characters may be different for TAC and BUFR DRIBU reports
      IF ((kobtyp == ndribu) .OR. (lship)) THEN
        iis = 1
        !   include redundancy between e.g. 'DRIBU 44505' and 'SYNOP 00505'
        IF ((yosghd(nsglo)(1:1) >= '0') .AND. (yosghd(nsglo)(1:1) <= '9')) iis=3
        iie = LEN_TRIM( yosghd(nsglo) )
        IF (yosghd(nsglo)(iis:iie) == ystid(iis:iie))  lredn (nsglo) =  .TRUE.
      ENDIF
!     IF ((lnatdwd) .OR. (lnatdwd2))  lredn (nsglo) =  .TRUE.
!     kcdtyp2  =  mosghd (nsglo ,nhcode)
!     lship2   =       (kcdtyp2 == nshscd) .OR. (kcdtyp2 == natshs)
!     lraso2   =  (kobtyp2 == ntemp) .OR. (kobtyp2 == npilot)
!     lredn (nsglo)  =        (     (yosghd(nsglo) == ystid)                   &
!                              .OR. (lraso .NEQV. lraso2)                      &
!                              .OR. ((lship) .AND. (lship2)))                  &
!                       .AND. (lairc   .EQV. lair2  )                          &
!                       .AND. (lsatbc  .EQV. lsatb2 )
      ! exclude merging (redundancy) between 2 (upper-air) reports where exactly
      ! 1 of the 2 reports has height as level info (-> uses standard atmosphere
      ! to derive a model state independent, but highly inaccurate pressure).
      IF (       BTEST( mosghd(nsglo,nhschr), nvlzbp )                         &
          .NEQV. BTEST( mosghd(nsgob,nhschr), nvlzbp ))  lredn (nsglo) = .FALSE.
    ENDIF
  ENDDO
  nsgo2 = 0
! first find redundancy between 2 rejected or 2 passive reports
  DO nsglo = 1 , nchkend
    IF ((lredn(nsglo)) .AND. (nsgo2 == 0) .AND. (npassiv >= 2)) THEN
      npassv2 = mosghd(nsglo,nhpass)
      IF ((npassv2 == 2) .AND. (IBCLR( mosghd(nsglo,nhflag), FL_OBSTYPE) > 0)) &
        npassv2 = 3
      IF (npassv2 == npassiv)  nsgo2 = nsglo
    ENDIF
  ENDDO
! otherwise find redundancy between 2 active reports (or 1 passive + 1 rejected)
  DO nsglo = 1 , nchkend
    IF (      (lredn(nsglo)) .AND. (nsgo2 == 0)                                &
        .AND. (mosghd(nsglo,nhpass)/2 == npassiv/2))  nsgo2 = nsglo
  ENDDO
! otherwise find redundancy between active and passive report
  DO nsglo = 1 , nchkend
    IF (      (lredn(nsglo)) .AND. (nsgo2 == 0))  nsgo2 = nsglo
  ENDDO
  IF (nsgo2 > 0) THEN

! decide which report is redundant
! --------------------------------

    npassv2  =  mosghd(nsgo2,nhpass)
    IF ((npassv2 == 2) .AND. (IBCLR( mosghd(nsgo2,nhflag), FL_OBSTYPE) > 0))   &
      npassv2 = 3
    lraso2   =       (mosghd(nsgo2,nhobtp) == ntemp)                           &
                .OR. (mosghd(nsgo2,nhobtp) == npilot)
    lnatdwd2 =     (mosghd(nsgo2,nhkz) == 40005) .OR.(mosghd(nsgo2,nhkz) == 5) &
              .OR. (mosghd(nsgo2,nhkz) == 10005)
    lrepla   = .TRUE.
! a single-level aircraft report used in multi-level report is treated as active
    IF     ((npassiv <  2) .AND. (npassv2 >= 2)) THEN
      nodract  = nsgob
      lrepla   = .FALSE.
    ELSEIF ((npassiv >= 2) .AND. (npassv2 <  2)) THEN
      nodract  = nsgo2
      lrepla   = .FALSE.
    !   prefer non-METAR
    ELSEIF ((kcdtyp /= nmetar) .AND. (mosghd(nsgo2,nhcode) == nmetar)) THEN
      nodract  = nsgob
    ELSEIF ((kcdtyp == nmetar) .AND. (mosghd(nsgo2,nhcode) /= nmetar)) THEN
      nodract  = nsgo2
    !   prefer non-natdwd
    ELSEIF ((.NOT. lnatdwd) .AND. (lnatdwd2)) THEN
      nodract  = nsgob
    ELSEIF ((lnatdwd) .AND. (.NOT. lnatdwd2)) THEN
      nodract  = nsgo2
    !   prefer Synop over non-Synop (e.g. TEMP-derived surface report)
    ELSEIF ((kobtyp == nsynop) .AND. (mosghd(nsgo2,nhobtp) /= nsynop)) THEN
      nodract  = nsgob
    ELSEIF ((kobtyp /= nsynop) .AND. (mosghd(nsgo2,nhobtp) == nsynop)) THEN
      nodract  = nsgo2
    !   prefer non-TEMP-derived
    ELSEIF ((kcdtyp /= nsytmp) .AND. (mosghd(nsgo2,nhcode) == nsytmp)) THEN
      nodract  = nsgob
    ELSEIF ((kcdtyp == nsytmp) .AND. (mosghd(nsgo2,nhcode) /= nsytmp)) THEN
      nodract  = nsgo2
    !   prefer obs time not later than 'ydate_ref' (= end of period of fof/ekf)
    ELSEIF ((obtime <= epsy) .AND. (osghed(nsgo2,nhtime) >  epsy)) THEN
      nodract  = nsgob
    ELSEIF ((obtime >  epsy) .AND. (osghed(nsgo2,nhtime) <= epsy)) THEN
      nodract  = nsgo2
    !   Synop: prefer report valid at full hours
    ELSEIF ((kobtyp == nsynop).AND. (MOD(mosghd(nsgob,nhhrmn), 100) == 0)     &
                              .AND. (MOD(mosghd(nsgo2,nhhrmn), 100) /= 0)) THEN
      nodract  = nsgob
    ELSEIF ((kobtyp == nsynop).AND. (MOD(mosghd(nsgob,nhhrmn), 100) /= 0)     &
                              .AND. (MOD(mosghd(nsgo2,nhhrmn), 100) == 0)) THEN
      nodract  = nsgo2
    !   prefer obs time closer to 'ydate_ref' (analysis time)
    ELSEIF (obtime > osghed(nsgo2,nhtime)+0.01_wp) THEN
      nodract  = nsgob
    ELSEIF (obtime < osghed(nsgo2,nhtime)-0.01_wp) THEN
      nodract  = nsgo2
    ELSEIF (nstcor >  mosghd(nsgo2,nhcorr)) THEN
      nodract  = nsgob
    ELSEIF (nstcor <  mosghd(nsgo2,nhcorr)) THEN
      nodract  = nsgo2
!   ELSEIF ((kobtyp == nsynop) .AND. (mosghd(nsgob,nhkz) >= 10000)             &
!                              .AND. (mosghd(nsgo2,nhkz) <  10000)) THEN
!                 nodract  = nsgob
!     IF (lsytac) nodract  = nsgo2
!   ELSEIF ((kobtyp == nsynop) .AND. (mosghd(nsgob,nhkz) <  10000)             &
!                              .AND. (mosghd(nsgo2,nhkz) >= 10000)) THEN
!                 nodract  = nsgo2
!     IF (lsytac) nodract  = nsgob
    ELSEIF (      (kobtyp == nsynop)                                           &
            .AND. (mosghd(nsgob,nhkz) < mosghd(nsgo2,nhkz))) THEN
                  nodract  = nsgo2
      IF (lsytac) nodract  = nsgob
    ELSEIF (      (kobtyp == nsynop)                                           &
            .AND. (mosghd(nsgob,nhkz) > mosghd(nsgo2,nhkz))) THEN
                  nodract  = nsgob
      IF (lsytac) nodract  = nsgo2
    ELSE
      ! the first report is active, the second report (processed now) is passive
      nodract  = nsgo2
    ENDIF
    IF (nodract == nsgo2) THEN
      nodrrej  = nsgob
      nactpr   = 1
    ELSEIF (nodract == nsgob) THEN
      nodrrej  = nsgo2
      nactpr   = 2
    ENDIF
    lreplanl = lrepla
    IF (MOD(mosghd(nodract,nhhrmn), 100) /= 0)  lreplanl = .FALSE.
!   PRINT *,'ztt5 ', nodract, nodrrej, yosghd(nodract), osgbdy(nodract,nbst)   &
!                  , mosgbd(nodract,nbserr), nsgob

! replace missing data of active report by available data of rejected report
! --------------------------------------------------------------------------

    nreplpr = 0
    IF (     (.NOT. ((lraso) .OR. (lraso2) .OR. (lairc) .OR. (lsatbc)))        &
        .AND.(.NOT. BTEST( mosgbd(nodract,nbserr), nvrz ))                     &
        .AND.(      BTEST( mosgbd(nodrrej,nbserr), nvrz )) .AND. (lrepla)) THEN
      osgbdy (nodract,nbsp  ) = osgbdy (nodrrej,nbsp  )
      osgbdy (nodract,nbsplv) = osgbdy (nodrrej,nbsplv)
      osgbdy (nodract,nbsz  ) = osgbdy (nodrrej,nbsz  )
      osgbdy (nodract,nbszer) = osgbdy (nodrrej,nbszer)
      mosgbd (nodract,nbserr) = IBSET( mosgbd(nodract,nbserr), nvrz )
      CALL MVBITS( mosgbd(nodrrej,nbsflg), nvfzbp, nvfaoc                      &
                 , mosgbd(nodract,nbsflg), nvfzbp )
      nreplpr = nreplpr + 1
    ENDIF
    IF (     (.NOT. BTEST( mosgbd(nodract,nbserr), nvru ))                     &
        .AND.(      BTEST( mosgbd(nodrrej,nbserr), nvru )) .AND. (lrepla)) THEN
      osgbdy (nodract,nbsu  ) = osgbdy (nodrrej,nbsu  )
      osgbdy (nodract,nbsv  ) = osgbdy (nodrrej,nbsv  )
      osgbdy (nodract,nbsuer) = osgbdy (nodrrej,nbsuer)
      mosgbd (nodract,nbserr) = IBSET( mosgbd(nodract,nbserr), nvru )
      CALL MVBITS( mosgbd(nodrrej,nbsflg), nvfubp, nvfaoc                      &
                 , mosgbd(nodract,nbsflg), nvfubp )
      nreplpr = nreplpr + 2
    ENDIF
    IF (     (.NOT. BTEST( mosgbd(nodract,nbserr), nvrt ))                     &
        .AND.(      BTEST( mosgbd(nodrrej,nbserr), nvrt )) .AND. (lrepla)) THEN
      osgbdy (nodract,nbst  ) = osgbdy (nodrrej,nbst  )
      osgbdy (nodract,nbster) = osgbdy (nodrrej,nbster)
      mosgbd (nodract,nbserr) = IBSET( mosgbd(nodract,nbserr), nvrt )
      CALL MVBITS( mosgbd(nodrrej,nbsflg), nvftbp, nvfaoc                      &
                 , mosgbd(nodract,nbsflg), nvftbp )
      nreplpr = nreplpr + 4
    ENDIF
    IF (     (.NOT. BTEST( mosgbd(nodract,nbserr), nvrq ))                     &
        .AND.(      BTEST( mosgbd(nodrrej,nbserr), nvrq )) .AND. (lrepla)) THEN
      osgbdy (nodract,nbsrh ) = osgbdy (nodrrej,nbsrh )
      osgbdy (nodract,nbsqer) = osgbdy (nodrrej,nbsqer)
      osgbdy (nodract,nbsdrh) = osgbdy (nodrrej,nbsdrh)
      mosgbd (nodract,nbserr) = IBSET( mosgbd(nodract,nbserr), nvrq )
      CALL MVBITS( mosgbd(nodrrej,nbsflg), nvfqbp, nvfaoc                      &
                 , mosgbd(nodract,nbsflg), nvfqbp )
      nreplpr = nreplpr + 8
    ENDIF
    IF (     (osgbdy(nodract,nbsr12) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsr12) >  rmdich) .AND. (lreplanl)) THEN
      osgbdy (nodract,nbsr12) = osgbdy (nodrrej,nbsr12)
      nreplpr = nreplpr + 16
    ENDIF
    IF (     (osgbdy(nodract,nbscl ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbscl ) >  rmdich) .AND. (lrepla)) THEN
      osgbdy (nodract,nbscl ) = osgbdy (nodrrej,nbscl )
      nreplpr = nreplpr + 32
    ENDIF
    IF (     (osgbdy(nodract,nbscm ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbscm ) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbscm ) = osgbdy (nodrrej,nbscm )
    IF (     (osgbdy(nodract,nbsch ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsch ) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbsch ) = osgbdy (nodrrej,nbsch )
    IF (     (osgbdy(nodract,nbsct ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsct ) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbsct ) = osgbdy (nodrrej,nbsct )
    IF (     (osgbdy(nodract,nbscbs) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbscbs) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbscbs) = osgbdy (nodrrej,nbscbs)
    IF (     (osgbdy(nodract,nbscil) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbscil) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbscil) = osgbdy (nodrrej,nbscil)
    IF (     (osgbdy(nodract,nbsvis) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsvis) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbsvis) = osgbdy (nodrrej,nbsvis)
    IF (     (osgbdy(nodract,nbspst) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbspst) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbspst) = osgbdy (nodrrej,nbspst)
    IF (     (osgbdy(nodract,nbsff ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsff ) >  rmdich) .AND. (lrepla)                 &
        .AND.(osgbdy(nodract,nbsu  ) >  rmdich))                               &
      osgbdy (nodract,nbsff ) = osgbdy (nodrrej,nbsff )
    IF (     (osgbdy(nodract,nbsdd ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsdd ) >  rmdich) .AND. (lrepla)                 &
        .AND.(osgbdy(nodract,nbsu  ) >  rmdich))                               &
      osgbdy (nodract,nbsdd ) = osgbdy (nodrrej,nbsdd )
    IF (     (osgbdy(nodract,nbstd ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbstd ) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbstd ) = osgbdy (nodrrej,nbstd )
    IF (     (osgbdy(nodract,nbsrr1) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrr1) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrr1) = osgbdy (nodrrej,nbsrr1)
    IF (     (osgbdy(nodract,nbsrr3) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrr3) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrr3) = osgbdy (nodrrej,nbsrr3)
    IF (     (osgbdy(nodract,nbsrr6) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrr6) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrr6) = osgbdy (nodrrej,nbsrr6)
    IF (     (osgbdy(nodract,nbsr24) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsr24) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsr24) = osgbdy (nodrrej,nbsr24)
    IF (     (osgbdy(nodract,nbsfg1) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsfg1) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsfg1) = osgbdy (nodrrej,nbsfg1)
    IF (     (osgbdy(nodract,nbsfg3) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsfg3) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsfg3) = osgbdy (nodrrej,nbsfg3)
    IF (     (osgbdy(nodract,nbsfg6) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsfg6) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsfg6) = osgbdy (nodrrej,nbsfg6)
    IF (     (osgbdy(nodract,nbstn ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbstn ) >  rmdich) .AND. (lreplanl)) THEN
      osgbdy (nodract,nbstn ) = osgbdy (nodrrej,nbstn )
      mosgbd (nodract,nbsttr) = ireplace( mosgbd(nodract,nbsttr), ntnbp, ntxoc &
                                        , mosgbd(nodrrej,nbsttr), ntnbp )
    ENDIF
    IF (     (osgbdy(nodract,nbstx ) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbstx ) >  rmdich) .AND. (lreplanl)) THEN
      osgbdy (nodract,nbstx ) = osgbdy (nodrrej,nbstx )
      mosgbd (nodract,nbsttr) = ireplace( mosgbd(nodract,nbsttr), ntxbp, ntxoc &
                                        , mosgbd(nodrrej,nbsttr), ntxbp )
    ENDIF
    IF (     (osgbdy(nodract,nbsrad) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrad) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrad) = osgbdy (nodrrej,nbsrad)
    IF (     (osgbdy(nodract,nbsrdd) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrdd) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrdd) = osgbdy (nodrrej,nbsrdd)
    IF (     (osgbdy(nodract,nbsrdl) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbsrdl) >  rmdich) .AND. (lreplanl))              &
      osgbdy (nodract,nbsrdl) = osgbdy (nodrrej,nbsrdl)
    IF (     (osgbdy(nodract,nbshsw) <  rmdich)                                &
        .AND.(osgbdy(nodrrej,nbshsw) >  rmdich) .AND. (lrepla))                &
      osgbdy (nodract,nbshsw) = osgbdy (nodrrej,nbshsw)
    IF (     (osgbdy(nodract,nbscwg) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbscwg) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbscwg) = mosgbd (nodrrej,nbscwg)
    IF (     (osgbdy(nodract,nbswwe) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbswwe) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbswwe) = mosgbd (nodrrej,nbswwe)
    IF (     (osgbdy(nodract,nbstur) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbstur) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbstur) = mosgbd (nodrrej,nbstur)
    IF (     (osgbdy(nodract,nbsclg) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbsclg) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbsclg) = mosgbd (nodrrej,nbsclg)
    IF (     (osgbdy(nodract,nbscl1) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbscl1) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbscl1) = mosgbd (nodrrej,nbscl1)
    IF (     (osgbdy(nodract,nbscl2) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbscl2) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbscl2) = mosgbd (nodrrej,nbscl2)
    IF (     (osgbdy(nodract,nbscl3) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbscl3) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbscl3) = mosgbd (nodrrej,nbscl3)
    IF (     (osgbdy(nodract,nbscl4) ==   imdi)                                &
        .AND.(osgbdy(nodrrej,nbscl4) /=   imdi) .AND. (lrepla))                &
      mosgbd (nodract,nbscl4) = mosgbd (nodrrej,nbscl4)

! flag redundant report
! ---------------------

! If data of a current redundant aircraft report replace missing data of an
! 'old' (read at a previous timestep) active report then replace (the order of)
! these two, so that the extended active report gets a report index which
! suggests that it has been read at the current timestep.
! This will facilite the inclusion of this aircraft identifier in the list for
! which multi-level aircraft reports are (re-)created (see src_obs_air_org.f90).
! (And it avoids possible needs for redundancy checking of single-level against
!  multi-level aircraft reports.)
    IF ((lairc) .AND. (nreplpr > 0) .AND. (nodract <= ntotsgo)) THEN
      zsgbdy (      1:mxsbdy) = osgbdy (nsgob,1:mxsbdy)
      osgbdy (nsgob,1:mxsbdy) = osgbdy (nsgo2,1:mxsbdy)
      osgbdy (nsgo2,1:mxsbdy) = zsgbdy (      1:mxsbdy)
      mzsgbd (      1:mxsbdf) = mosgbd (nsgob,1:mxsbdf)
      mosgbd (nsgob,1:mxsbdf) = mosgbd (nsgo2,1:mxsbdf)
      mosgbd (nsgo2,1:mxsbdf) = mzsgbd (      1:mxsbdf)
      zsghed (      1:mxshed) = osghed (nsgob,1:mxshed)
      osghed (nsgob,1:mxshed) = osghed (nsgo2,1:mxshed)
      osghed (nsgo2,1:mxshed) = zsghed (      1:mxshed)
      mzsghd (      1:mxshdf) = mosghd (nsgob,1:mxshdf)
      mosghd (nsgob,1:mxshdf) = mosghd (nsgo2,1:mxshdf)
      mosghd (nsgo2,1:mxshdf) = mzsghd (      1:mxshdf)
      yosghd (nsgob         ) = yosghd (nsgo2)
      yosghd (nsgo2         ) = ystid
      nodract                 = nsgob
      nodrrej                 = nsgo2
    ENDIF
    IF ((lverpas) .AND. (.NOT. lnatdwd) .AND. (.NOT. lnatdwd2)) THEN
! ordinary flagging
      mosghd (nodrrej,nhpass)  = 2
      mosghd (nodrrej,nhschr)  = IBSET( mosghd(nodrrej,nhschr) , nvrdbp )
      mosghd (nodrrej,nhschr)  = IBSET( mosghd(nodrrej,nhschr) , nvpsbp )
      mosghd (nodrrej,nhflag)  = IBSET( mosghd(nodrrej,nhflag) , FL_REDUNDANT )
    ELSE
! this report will be discarded in 'obs_cdf_del_old_rep'
      mosghd (nodrrej,nhpass)  = -1
    ENDIF
!   PRINT *,'ztt6 ', nodract, nodrrej, yosghd(nodract), osgbdy(nodract,nbst)   &
!                  , mosgbd(nodract,nbserr), nsgob, mosghd(nodract,nhpass)     &
!                  , mosghd(nodrrej,nhpass)

! update statistics and report events, and prepare control output
! ---------------------------------------------------------------

    kobtypr = mosghd(nodrrej,nhobtp)
    kcdtypr = mosghd(nodrrej,nhcode)
    icma    = i_cma ( kobtypr , kcdtypr )
!             =====
    IF ((kobtypr /= ntemp) .AND. (kobtypr /= npilot)) THEN
      IF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
        cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        neventr (neredn,icma) = neventr(neredn,icma) + 1
      ELSE
        cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      ENDIF
      IF (      lverpas) cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      IF (.NOT. lverpas) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ELSEIF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
      neventr (neredx,icma) = neventr(neredx,icma) + 1
    ENDIF

    ilen = 15 + 2*istrej
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt10
      outbuf(nacout+ 3) = kcdtypr
      DO iiv = 1,5
        IF (osgbdy(nodract,iiv) > rmdich) THEN
          outbuf(nacout+3+iiv) = NINT(osgbdy(nodract,iiv)*rscale)
        ELSE
          outbuf(nacout+3+iiv) = imdi
        ENDIF
      ENDDO
      outbuf(nacout+ 9) = mosghd(nodract,nhcode)
      outbuf(nacout+10) = NINT(osghed(nodract,nhilon)*rscale)
      outbuf(nacout+11) = NINT(osghed(nodract,nhjlat)*rscale)
      outbuf(nacout+12) = nactpr
      outbuf(nacout+13) = 0
      IF (mosghd(nodract,nhcorr) /= mosghd(nodrrej,nhcorr))                    &
        outbuf(nacout+13) = 1
      outbuf(nacout+14) = nreplpr
      outbuf(nacout+15) = NINT(osghed(nodract,nhtime)*rscale)
      DO icl = 1 , istrej
        outbuf(nacout+15       +icl) = ICHAR( yosghd(nodrrej) (icl:icl) )
        outbuf(nacout+15+istrej+icl) = ICHAR( yosghd(nodract) (icl:icl) )
      ENDDO
      nacout  = nacout + ilen
    ENDIF  ! lprodr
  ENDIF    ! (nsgo2 > 0)
ENDDO


!-------------------------------------------------------------------------------
! Section 2.0  Redundancy check for multi-level reports
!-------------------------------------------------------------------------------

#ifdef __COSMO__
  DO nmlob = 1 , nnmlnew
    igplat = momlhd (ntotml-nnmlnew+nmlob,nhjtot)
    cosolat2 (nmlob) = COS( (r_startlat_tot+(igplat-1)*r_dlat) * r_degrad )
    cosolat2 (nmlob) = cosolat2(nmlob) * cosolat2(nmlob)
  ENDDO
#endif

! iactio = 1 : redundancy check
! iactio = 2 : merging TEMP / PILOT parts A, B, C, D

DO iactio = 1 , 2
!~~~~~~~~~~~~~~~~

 loop_over_new_reports:  DO nmlob = ntotml - nnmlnew + 1 , ntotml
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  (For redundancy of multi-level AIREPs (not active),
!   lowest data have to be (quasi-) colocated)

! get some information of current report
! --------------------------------------

!  observation and code type

  kobtyp   =  momlhd (nmlob,nhobtp)
  kcdtyp   =  momlhd (nmlob,nhcode)
  lraso    =  (kobtyp == npilot) .AND. (     (kcdtyp == nldpcd)                &
                                        .OR. (kcdtyp == nshpcd))
  lraso    =  (kobtyp == ntemp) .OR. (lraso)
  lwprofl  =  (kobtyp == npilot) .AND. (     (kcdtyp == nwp_eu)                &
                                        .OR. (kcdtyp == nra_eu)                &
                                        .OR. (kcdtyp == npr_us)                &
                                        .OR. (kcdtyp == nravad)                &
                                        .OR. (kcdtyp == nwlidr))
  lairc    =  (kobtyp == nairep)
  ltower   =  (kcdtyp == ntower) .OR. (kcdtyp == ntowic)

  ! perform (iactio = 2) only for radiosondes
  IF ((iactio == 2) .AND. (.NOT. lraso))             CYCLE loop_over_new_reports

!  get WMO station identity, station correction bit, and passive flag

  ystid    =  yomlhd (nmlob)
  nstcor   =  momlhd(nmlob,nhcorr)
  npassiv  =  momlhd(nmlob,nhpass)
  IF (npassiv == -1)                                 CYCLE loop_over_new_reports
  IF ((npassiv == 2) .AND. (IBCLR( momlhd(nmlob,nhflag), FL_OBSTYPE ) > 0))    &
    npassiv = 3

!  get 4-d location, and cos(latitude) for distances

  obtime = omlhed (nmlob,nhtime)
! IF (lraso) obtime = omlhed (nmlob,nhsynt)
  igplon = momlhd (nmlob,nhitot)
  igplat = momlhd (nmlob,nhjtot)
  obbalt = omlhed (nmlob,nhalt )
  IF (lairc) obbalt = - omlbdy (nmlob,1,nbtplv)

!  preset thresholds for 4-d colocation

  IF (kobtyp == nairep) THEN
    rtmlmt = rtmlair
    rhzlmt = rhzlair
    rvtlmt = rvtlair * 100._wp
    rdplmt = rvtlair * 100._wp + epsy
    rprlmt = rprlim  * 100._wp + epsy
  ELSE
    rtmlmt = rtmlim
    rhzlmt = rhzlim
    rvtlmt = rvtlim
    rdplmt = rdplim * 100._wp + epsy
    rprlmt = rprlim * 100._wp + epsy
    IF (lraso)        rtmlmt = rtmlrs
    IF (ltower)       rtmlmt = rtmltow
    ! apply a small value for 'rdplmt' for putting TEMP parts A,B,C,D together
    IF (iactio == 2)  rdplmt = rprlmt
    ! apply a very small value for 'rdplmt' for towers
    IF (ltower)       rprlmt = MIN( rprlmt ,.5_wp* picoshift )
  ENDIF
  rhzlmt = rhzlmt * rhzlmt

! check if current report meets the redundancy criteria with other reports
! ------------------------------------------------------------------------
! (do this with a vectorisable loop)

  nchkend = nmlob - 1

  DO nmllo = 1 , nchkend
    lredn (nmllo)  = .FALSE.
    kobtyp2  =  momlhd (nmllo ,nhobtp)
    kcdtyp2  =  momlhd (nmllo ,nhcode)
    lraso2   =  (kobtyp2 == npilot) .AND. (     (kcdtyp2 == nldpcd)            &
                                           .OR. (kcdtyp2 == nshpcd))
    lraso2   =  (kobtyp2 == ntemp) .OR. (lraso2)
    lair2    =  (kobtyp2 == nairep)
    ltower2  =  (kcdtyp2 == ntower) .OR. (kcdtyp2 == ntowic)

    ! get 4-d location and pressure of second report
    obtime2 = omlhed (nmllo,nhtime)
!   IF (lraso2) obtime2 = omlhed (nmllo,nhsynt)
    igplon2 = momlhd (nmllo,nhitot)
    igplat2 = momlhd (nmllo,nhjtot)
    obbalt2 = omlhed (nmllo,nhalt )
    IF (lair2)  obbalt2 = - omlbdy (nmllo,1,nbtplv)

!   check for horizontal and vertical (quasi-) colocation and simultaneity

#ifdef __COSMO__
!   rhzob  =    cosolat2(ntotml-nnmlnew+nmlob)                                 &
    rhzob  =    cosolat2(nmlob-(ntotml-nnmlnew))                               &
               *((igplon2-igplon) *r_dlon) *((igplon2-igplon) *r_dlon)         &
             +  ((igplat2-igplat) *r_dlat) *((igplat2-igplat) *r_dlat)
#else
    rhzob  = r2d **2 * SUM ((  gen_grid% xnglob (igplon ,igplat ,:,1)          &
                             - gen_grid% xnglob (igplon2,igplat2,:,1)) **2)
#endif
    lcoloc =      (      rhzob * rdegkm2     <=  rhzlmt )                      &
            .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmt )                      &
            .AND. ( ABS( obtime2 - obtime )  <=  rtmlmt )
    ! reproducibility irrespective of domain decomposition ensured by allowing
    ! for redundancy only between reports assigned to the same grid point
    IF (lredn_repro)  lcoloc = lcoloc .AND. (igplon2 == igplon)                &
                                      .AND. (igplat2 == igplat)
    IF ((lcoloc) .AND. (momlhd(nmllo,nhpass) /= -1)) THEN

!   check for 'similar' observation types (but do not request equal station ID!)
!                                          ------------------------------------
      lwprof2  =  (kobtyp2 == npilot) .AND. (     (kcdtyp2 == nwp_eu)          &
                                             .OR. (kcdtyp2 == nra_eu)          &
                                             .OR. (kcdtyp2 == npr_us)          &
                                             .OR. (kcdtyp2 == nravad)          &
                                             .OR. (kcdtyp2 == nwlidr))
      lredn (nmllo)  =        (lraso   .EQV. lraso2 )                          &
                        .AND. (lwprofl .EQV. lwprof2)                          &
                        .AND. (lairc   .EQV. lair2  )                          &
                        .AND. (ltower  .EQV. ltower2)                          &
                        .AND. ((kobtyp == kobtyp2) .OR. (lraso))
      ! (allow merging of TEMP with PILOT, but ... )
      ! exclude merging (redundancy) of RASS with Wind Profiler
      IF (lwprofl)  lredn (nmllo)  =  (lredn(nmllo)) .AND. (kcdtyp == kcdtyp2)
      ! exclude merging (redundancy) of radiosonde descent w. another code type
      lredn (nmllo)  =  (lredn(nmllo)) .AND. (      (kcdtyp  == ntdesc)        &
                                              .EQV. (kcdtyp2 == ntdesc))
      ! exclude merging (redundancy) between 2 reports where exactly 1 of the 2
      ! reports has height as level information (and uses standard atmosphere
      ! to derive a model state independent, but highly inaccurate pressure).
      IF (       BTEST( momlhd(nmllo,nhschr), nvlzbp )                         &
          .NEQV. BTEST( momlhd(nmlob,nhschr), nvlzbp ))  lredn (nmllo) = .FALSE.
!     IF (       BTEST( momlhd(nmllo,nhschr), nvlzbp )                         &
!         .NEQV. BTEST( momlhd(nmlob,nhschr), nvlzbp ))                        &
!       PRINT *,'coloc1 ',ystid, nmlob, nmllo, lredn(nmllo), obtime, obtime2   &
!              , BTEST( momlhd(nmlob,nhschr), nvlzbp )                         &
!              , BTEST( momlhd(nmllo,nhschr), nvlzbp ), kcdtyp, kcdtyp2, lcoloc
    ENDIF
  ENDDO

!   for radiosonde: check whether redundancy or different TEMP parts

  lmerge = .FALSE.
  IF (lraso) THEN
    DO nmllo = 1 , nchkend
      IF (lredn(nmllo)) THEN
!       kz1 = momlhd(nmlob,nhkz)
!       kz2 = momlhd(nmllo,nhkz)
!       lmerge = (      (kz1 /= imdi) .AND. (kz2 /= imdi) .AND. (kz1 /= kz2)   &
!                 .AND. (MIN(kz1,kz2) > 500) .AND. (MAX(kz1,kz2) < 800))
        nlev = momlhd(nmlob,nhnlev)
        nlv2 = momlhd(nmllo,nhnlev)
        nlvcom = -1
        IF (     (omlbdy(nmlob,nlev,nbtplv) >= omlbdy(nmllo,1,nbtplv))         &
            .OR. (omlbdy(nmllo,nlv2,nbtplv) >= omlbdy(nmlob,1,nbtplv))) THEN
          ! if 1 report is fully above the other then assume TEMP parts
          lmerge = .TRUE.
          ! if only 1 report is superobbed then do not merge
        ELSEIF (       (momlhd(nmlob,nhcorr) >= 10)                            &
                .NEQV. (momlhd(nmllo,nhcorr) >= 10)) THEN
          lmerge = .FALSE.
        ELSE
          ! if <= 1/2 of the levels of the smaller report exist in the larger
          ! report (with a margin of 0.99 hPa)
          ! then assume that these reports are different TEMP parts
!         ! report then assume that these reports are different TEMP parts
          nlvcom = 0
          DO klev = 1 , nlev
            DO kact = 1 , nlv2
              IF (ABS( omlbdy(nmlob,klev,nbtplv)                               &
                      -omlbdy(nmllo,kact,nbtplv) ) < 99._wp) nlvcom = nlvcom + 1
!                     -omlbdy(nmllo,kact,nbtplv) ) < epsy)  nlvcom = nlvcom + 1
            ENDDO
          ENDDO
          lmerge = (     (MIN( nlev,nlv2 ) >= 2* nlvcom)                       &
                    .OR. (MIN( nlev,nlv2 ) <= 2) )
        ENDIF
        lredn (nmllo)  =  (lmerge  .EQV.  (iactio == 2))
      ENDIF
    ENDDO
  ENDIF

! select the most appropriate of the other reports that meet redundancy criteria
! ------------------------------------------------------------------------------

  nmlo2 = 0
! first find redundancy between 2 rejected or 2 passive reports
  DO nmllo = 1 , nchkend
    IF ((lredn(nmllo)) .AND. (nmlo2 == 0) .AND. (npassiv >= 2)) THEN
      npassv2 = momlhd(nmllo,nhpass)
      IF ((npassv2 == 2) .AND. (IBCLR( momlhd(nmllo,nhflag), FL_OBSTYPE) > 0)) &
        npassv2 = 3
      IF (npassv2 == npassiv)  nmlo2 = nmllo
    ENDIF
  ENDDO
! otherwise find redundancy between 2 active reports (or 1 passive + 1 rejected)
  DO nmllo = 1 , nchkend
    IF (      (lredn(nmllo)) .AND. (nmlo2 == 0)                                &
        .AND. (momlhd(nmllo,nhpass) == momlhd(nmlob,nhpass)))  nmlo2 = nmllo
  ENDDO
! otherwise find redundancy between active and passive report
  DO nmllo = 1 , nchkend
    IF (      (lredn(nmllo)) .AND. (nmlo2 == 0))  nmlo2 = nmllo
  ENDDO
! PRINT *,'coloc2 ',ystid, nmlob, nmlo2, nchkend, momlhd(nmllo,nhpass), npassiv

! loop of length 1 or 2 only:
! auto-checking and possibly ordinary redundancy checking
! -------------------------------------------------------

!NEC$ nomove
  DO ktypchk = iactio , 1 + MIN( nmlo2, 1 )
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  at first, the report is auto-checked against itself to get rid of (approx.)
!  colocated levels within the report
    IF (ktypchk == 1) THEN
      nodract  = nmlob
      nodrrej  = nmlob
      lrepla   = .TRUE.

! ordinary redundancy checking: decide which report is redundant
! --------------------------------------------------------------

    ELSEIF (ktypchk == 2) THEN
      npassv2  =  momlhd(nmlo2,nhpass)
      IF ((npassv2 == 2) .AND. (IBCLR( momlhd(nmlo2,nhflag), FL_OBSTYPE) > 0)) &
        npassv2 = 3
      lrepla   = .TRUE.
      lrbufr1 = .FALSE.
      lrbufr2 = .FALSE.
      zptop1  = -1._wp
      zptop2  = -1._wp
      zpbase1 = -2._wp
      zpbase2 = -2._wp
!     IF ((lraso) .AND. (momlhd(nmlo2,nhobtp) == ntemp)) THEN
        DO klev = 1 , momlhd(nmlob,nhnlev)
          !  TEMP pressure level accuracy: hPa for TAC, 0.1 hPa for BUFR
          IF (MOD( NINT( omlbdy(nmlob,klev,nbtplv) ), 100 ) > 0) lrbufr1 =.TRUE.
          IF (BTEST( momlbd(nmlob,klev,nbtlsg), LS_SUPEROBS ) ) THEN
            zptop1 = omlbdy(nmlob,klev,nbtplv)
            IF (zpbase1 < 0._wp)  zpbase1 = omlbdy(nmlob,klev,nbtplv)
          ENDIF
        ENDDO
        DO klev = 1 , momlhd(nmlo2,nhnlev)
          IF (MOD( NINT( omlbdy(nmlo2,klev,nbtplv) ), 100 ) > 0) lrbufr2 =.TRUE.
          IF (BTEST( momlbd(nmlo2,klev,nbtlsg), LS_SUPEROBS ) ) THEN
            zptop2 = omlbdy(nmlo2,klev,nbtplv)
            IF (zpbase2 < 0._wp)  zpbase2 = omlbdy(nmlo2,klev,nbtplv)
          ENDIF
        ENDDO
!     ENDIF
!  determine indices of active resp. redundant report
      !   prefer active over merged over passive (rejected) report
      IF     (npassiv < npassv2) THEN
        nodract  = nmlob
        lrepla   = .FALSE.
      ELSEIF (npassiv > npassv2) THEN
        nodract  = nmlo2
        lrepla   = .FALSE.
      !   prefer TEMP over PILOT sonde report
      ELSEIF ((kobtyp == ntemp) .AND. (momlhd(nmlo2,nhobtp) /= ntemp)) THEN
        nodract  = nmlob
      ELSEIF ((kobtyp /= ntemp) .AND. (momlhd(nmlo2,nhobtp) == ntemp)) THEN
        nodract  = nmlo2
      !   multi-level aircraft: prefer report with lower lowest level
      ELSEIF ((lairc) .AND. (obbalt2-obbalt >  rvtlmt)) THEN
        nodract  = nmlob
      ELSEIF (lairc) THEN
        nodract  = nmlo2
      !   prefer high-resolution superobbed report over other report
      ELSEIF (nstcor >  momlhd(nmlo2,nhcorr)) THEN
        nodract  = nmlob
      ELSEIF (nstcor <  momlhd(nmlo2,nhcorr)) THEN
        nodract  = nmlo2
      !   prefer larger vertical range of superobbed levels
      ELSEIF (zpbase1-zptop1 > zpbase2-zptop2+.5_wp) THEN
        nodract  = nmlob
      ELSEIF (zpbase2-zptop2 > zpbase1-zptop1+.5_wp) THEN
        nodract  = nmlo2
      !   prefer BUFR over TAC TEMP report
      ELSEIF ((lrbufr1) .AND. (.NOT. lrbufr2)) THEN
        nodract  = nmlob
      ELSEIF ((lrbufr2) .AND. (.NOT. lrbufr1)) THEN
        nodract  = nmlo2
      !   prefer TEMP parts A, B over parts C, D (-> obs time, decoding time)
      ELSEIF ((lraso) .AND. (omlbdy(nmlob,1,nbtplv) >= 15000._wp)              &
                      .AND. (omlbdy(nmlo2,1,nbtplv) <  15000._wp)) THEN
        nodract  = nmlob
      ELSEIF ((lraso) .AND. (omlbdy(nmlo2,1,nbtplv) >= 15000._wp)              &
                      .AND. (omlbdy(nmlob,1,nbtplv) <  15000._wp)) THEN
        nodract  = nmlo2
      ELSEIF ((lraso) .AND. (momlhd(nmlob,nhnlev) > momlhd(nmlo2,nhnlev))) THEN
        !   usually, TEMP part B has most levels   (-> obs time, decoding time)
        nodract  = nmlob
      ELSEIF ((lwprofl) .AND. (  omlbdy(nmlob,1,nbtplv)                        &
                               < omlbdy(nmlo2,1,nbtplv))) THEN
        !   prefer wind profiler high mode (larger pulse length + hence signal-
        !   to-noise ratio) over low mode (too high resolution at upper levels)
        nodract  = nmlob
      ELSEIF ((lwprofl) .AND. (  omlbdy(nmlob,1,nbtplv)                        &
                               > omlbdy(nmlo2,1,nbtplv))) THEN
        nodract  = nmlo2
      !   prefer obs time not later than 'ydate_ref' (= end of fof/ekf period)
      ELSEIF ((ltower) .AND. (obtime <= epsy) .AND. (omlhed(nmlo2,nhtime) >  epsy)) THEN
        nodract  = nmlob
      ELSEIF ((ltower) .AND. (obtime >  epsy) .AND. (omlhed(nmlo2,nhtime) <= epsy)) THEN
        nodract  = nmlo2
      !   prefer obs time closer to 'ydate_ref' (analysis time)
      ELSEIF ((ltower) .AND. (obtime > omlhed(nmlo2,nhtime)+0.01_wp)) THEN
        nodract  = nmlob
      ELSEIF ((ltower) .AND. (obtime < omlhed(nmlo2,nhtime)-0.01_wp)) THEN
        nodract  = nmlo2
      ELSE
        nodract  = nmlo2
      ENDIF
      IF (nodract == nmlo2)  nodrrej = nmlob
      IF (nodract == nmlob)  nodrrej = nmlo2
    ENDIF
    lsupoba  =  (momlhd(nodract,nhcorr) >= 10)
    lsupobr  =  (momlhd(nodrrej,nhcorr) >= 10)

!   IF (yomlhd(nodract) (1:5) == '03743')     &
!     WRITE(0,*) 'nodrrej ', ystid, nodrrej, nodract, nmlob, ntotml, nnmlnew   &
!                          , omlhed(nodrrej,nhtime), omlhed(nodract,nhtime)    &
!                          , momlhd(nodrrej,nhnlev), momlhd(nodract,nhnlev)    &
!                          , momlhd(nodrrej,nhnlev), omlbdy(nodrrej,1,nbtplv)  &
!                          , omlbdy(nodrrej,momlhd(nodrrej,nhnlev),nbtplv)     &
!                          , momlhd(nodract,nhnlev), omlbdy(nodract,1,nbtplv)  &
!                          , omlbdy(nodract,momlhd(nodrrej,nhnlev),nbtplv)

! Fill data of active report and supplementary data
! of redundant report into temporary field
! -------------------------------------------------

!  Fill header of active report into temporary field;
!  fill body of temporary field with missing values

    zmlhed (1:mxrhed) = omlhed (nodract,1:mxrhed)
    iomlhd (1:mxrhdf) = momlhd (nodract,1:mxrhdf)
    ystidz            = yomlhd (nodract)
    zmlbdy (1:2*maxmlv,1:mxrbdy) = rmdi
    iomlbd (1:2*maxmlv,1:mxrbdf) = imdi

!  Loop over levels of body of active report

    lsub = .TRUE.
    !   do not discard any tower reports (with different obs times)
    IF ((ltower) .AND. (ktypchk == 2)) THEN
      IF (ABS( obtime - omlhed(nmlo2,nhtime) ) > 0.01_wp) lsub = .FALSE.
    ENDIF
    ltop = .FALSE.
    klev = 0
    krej = 1

    Levels_acr: DO  kact = 1 , momlhd(nodract,nhnlev)
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      lactsurf = (      (BTEST( momlbd(nodract,kact,nbtlsg), LS_SURFACE ))     &
                  .AND. (momlbd(nodract,kact,nbtlsg) /= imdi))

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'LAPSEred1 ', yomlhd(nodract),kact, omlbdy(nodract,kact,nbtp)

! 'while' loop: Get all the rejected levels which are more than 'rdplmt' [Pa]
!               (e.g. 25 hPa) below the active level 'kact'.
!               From these levels, fill those levels, which are not within
!               'rdplmt' [Pa] of any level of the active report,
!               into the temporary field
!               (for wind profiler, 'rdplmt' is reduced to add (more) lower
!                levels of rejected low mode to the active high mode report)
!     IF (yomlhd(nodract) (1:5) == '03743')                                    &
!       WRITE(0,*) 'krejm   ',kact,krej, klev, omlbdy(nodrrej,krej,nbtp), rdplmt
      rdplmt2 = rdplmt
      IF ((lwprofl) .AND. (kact == 1) .AND. (omlbdy(nodract,1,nbtz) > rmdich)) &
        rdplmt2 = MIN( rdplmt , MAX( 500._wp , .5_wp* omlbdy(nodract,1,nbtz) ) )

      DO WHILE (      (   omlbdy(nodrrej,krej,nbtplv)                          &
                       >= omlbdy(nodract,kact,nbtplv)+rdplmt2)                 &
                .AND. (.NOT. ltop) .AND. (klev < 2*maxmlv))

!       IF (yomlhd(nodract) (1:5) == '03743')                                  &
!         WRITE(0,*) 'krej0   ', kact, krej, klev, omlbdy(nodrrej,krej,nbtp)   &
!                                                , omlbdy(nodract,kact,nbtp)
        lmand = .FALSE.
        IF (lrepla) THEN
          klev = klev + 1
          zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
          iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
          IF (BTEST( iomlbd(klev,nbterr), nvru ))  iomlhd (nhuexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrt ))  iomlhd (nhtexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrq ))  iomlhd (nhqexi) = 1
          lmand =   (rmod( omlbdy(nodrrej,krej,nbtplv),10000._wp ) < epsy2)    &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 5000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 7000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-15000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-25000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-85000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-92500._wp ) < epsy)
        ENDIF
        IF (krej <  momlhd(nodrrej,nhnlev)) THEN
          krej = krej + 1
        ELSE
          ltop = .TRUE.
        ENDIF
        ! set 'lsub' to .false. unless
        !   -  a mandatory level is added to a superobbed report, or
        !   -  only the top level of a superobbed redundant report (with nlev>2)
        !      is absent in (but within the range of) the active profile
        IF (.NOT. (     (lsupoba .AND. lmand)                                  &
                   .OR. (lsupobr .AND. ltop .AND. (krej > 2))))  lsub = .FALSE.
      ENDDO
!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'krej1   ', kact,krej,klev, omlbdy(nodrrej,krej,nbtp), rprlmt

! 'while' loop: Get all the rejected levels which are less than 'rdplmt' [pa],
!               but more than 'rprlmt' [pa] below the active level 'kact'.
!               From these levels, fill those observed quantities, which are not
!               present in the active level 'kact', into the temporary field.
!               (Fill all the observed quantities if the rejected level is
!                the surface level or a mandatory level.)
!        (Note: This may be important if e.g. a TEMP report contains only
!               temperature and humidity and a PILOT report has been derived
!               and disseminated from one and the same radiosonde ascent.)

      DO WHILE (      (   omlbdy(nodrrej,krej,nbtplv)                          &
                       >= omlbdy(nodract,kact,nbtplv)+rprlmt)                  &
                .AND. (.NOT. ltop) .AND. (klev < 2*maxmlv))
!       IF (yomlhd(nodract) (1:5) == '03743')                                  &
!         WRITE(0,*) 'krej0b  ', kact, krej, klev, omlbdy(nodrrej,krej,nbtp)   &
!                                                , omlbdy(nodract,kact,nbtp)
        lmand = .FALSE.
        laddu = .FALSE.
        laddt = .FALSE.
        laddq = .FALSE.
        IF (lrepla) THEN
          lmand =   (rmod( omlbdy(nodrrej,krej,nbtplv),10000._wp ) < epsy2)    &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 5000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 7000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-15000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-25000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-85000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-92500._wp ) < epsy)
          laddu =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvru)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvru))
          laddt =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvrt)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvrt))
          laddq =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvrq)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvrq))
        ENDIF
        IF (      (BTEST( momlbd(nodrrej,krej,nbtlsg), LS_SURFACE ))           &
            .AND. (momlbd(nodrrej,krej,nbtlsg) /= imdi)) THEN
          laddu = .TRUE.
          laddt = .TRUE.
          laddq = .TRUE.
        ENDIF
        IF ((laddu) .OR. (laddt) .OR. (laddq)) THEN
          klev = klev + 1
          zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
          iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
          IF (.NOT. laddu) zmlbdy (klev,nbtu  ) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtv  ) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtuer) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtw  ) = rmdi
          IF (.NOT. laddu) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvru )
          IF (.NOT. laddu) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvfubp)
          IF (.NOT. laddt) zmlbdy (klev,nbtt  ) = rmdi
          IF (.NOT. laddt) zmlbdy (klev,nbtter) = rmdi
          IF (.NOT. laddt) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvrt )
          IF (.NOT. laddt) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvftbp)
          IF (.NOT. laddq) zmlbdy (klev,nbtrh ) = rmdi
          IF (.NOT. laddq) zmlbdy (klev,nbtqer) = rmdi
          IF (.NOT. laddq) zmlbdy (klev,nbtdrh) = rmdi
          IF (.NOT. laddq) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvrq )
          IF (.NOT. laddq) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvfqbp)
          IF (BTEST( iomlbd(klev,nbterr), nvru ))  iomlhd (nhuexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrt ))  iomlhd (nhtexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrq ))  iomlhd (nhqexi) = 1
        ENDIF
        IF (krej <  momlhd(nodrrej,nhnlev)) THEN
          krej = krej + 1
        ELSE
          ltop = .TRUE.
        ENDIF
        ! set 'lsub' to .false. unless
        !   -  a mandatory level is added to a superobbed report, or
        !   -  only the top level of a superobbed redundant report (with nlev>2)
        !      is absent in (but within the range of) the active profile
        IF (.NOT. (     (lsupoba .AND. lmand)                                  &
                   .OR. (lsupobr .AND. ltop .AND. (krej > 2))))  lsub = .FALSE.
      ENDDO
!     (no need to check whether (klev > 2*maxmlv)

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'krej2   ',kact,krej, klev, omlbdy(nodrrej,krej,nbtp), rprlmt

! 'while' loop: if a rejected level which is less than 'rprlmt' [pa] below the
!               active level 'kact' is the surface level then add this level
!               and update 'krej' accordingly

      krjc = krej
      DO WHILE (      (  omlbdy(nodrrej,krjc,nbtplv)                           &
                       > omlbdy(nodract,kact,nbtplv)+epsy)                     &
                .AND. (.NOT. ltop) .AND. (klev < 2*maxmlv))
        IF (      (      BTEST( momlbd(nodrrej,krjc,nbtlsg), LS_SURFACE ))     &
            .AND. (.NOT. BTEST( momlbd(nodract,kact,nbtlsg), LS_SURFACE ))     &
            .AND. (momlbd(nodrrej,krjc,nbtlsg) /= imdi)) THEN
          klev = klev + 1
          zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krjc,1:mxrbdy)
          iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krjc,1:mxrbdf)
          krej = krjc + 1
          lsub = .FALSE.
        ENDIF
        IF (krjc <  momlhd(nodrrej,nhnlev)) THEN
          krjc = krjc + 1
        ELSE
          ltop = .TRUE.
        ENDIF
      ENDDO

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'krej3   ',kact,krej, klev, omlbdy(nodrrej,krej,nbtp), rprlmt

!  'nearly' colocated levels within the active report are removed
!  (except if the lower level is the surface level)
!  (however, the active level is complemented by colocated level data
!   in the auto-checking)

      IF ((kact > 1) .AND. (  omlbdy(nodract,MAX(kact-1,1),nbtplv)-rprlmt)     &
                            < omlbdy(nodract,kact,nbtplv)) THEN
        IF (BTEST( momlbd(nodract,kact-1,nbtlsg), LS_SURFACE ) .EQV.           &
            BTEST( momlbd(nodract,kact  ,nbtlsg), LS_SURFACE )      )          &
                                                                CYCLE Levels_acr
      ENDIF

!     IF ((kact == 1) .OR. (   omlbdy(nodract,MAX(kact-1,1),nbtplv)-rprlmt)    &
!                           >= omlbdy(nodract,kact,nbtplv)) THEN

!  Fill next level of active report into temporary field

      klev = klev + 1
      zmlbdy (klev,1:mxrbdy) = omlbdy (nodract,kact,1:mxrbdy)
      iomlbd (klev,1:mxrbdf) = momlbd (nodract,kact,1:mxrbdf)
      krjc  = krej
      kclo  = 0
      dprjc = - rdplmt

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'kclo1  ',kact, krej, klev, kclo, omlbdy(nodract,kact,nbtp) &
!                                                  , omlbdy(nodrrej,krej,nbtp)

! 'while' loop: get all rejected levels which are within 'rprlmt' [pa] of the
!               active level 'kact';
!               find the rejected level 'kclo' most distant to the act. level
!               'kact' (since that level is most likely to contain additional
!               relevant information, particularly when auto-checking).
!               (This search for rejected levels is ceased when a surface level
!                is found among the rejected levels.)

      DO WHILE (      (   omlbdy(nodrrej,krej,nbtplv)                          &
                       >  omlbdy(nodract,kact,nbtplv)-rprlmt)                  &
                .AND. (.NOT. ltop))
        dprjca = ABS(  omlbdy(nodrrej,krej,nbtplv)                             &
                     - omlbdy(nodract,kact,nbtplv) )
        leqsfc =  (BTEST( momlbd(nodrrej,krej,nbtlsg), LS_SURFACE ) .EQV.      &
                   BTEST( momlbd(nodract,kact,nbtlsg), LS_SURFACE )      )
!       IF (yomlhd(nodract) (1:5) == '03743')     &
!         WRITE(0,*) 'kclo1a ',kact,krej,klev, kclo, omlbdy(nodract,kact,nbtp) &
!                                                  , omlbdy(nodrrej,krej,nbtp) &
!             , momlbd(nodrrej,krej,nbtlsg), momlbd(nodract,kact,nbtlsg), leqsfc
        IF ((dprjca >= dprjc) .AND. (leqsfc)) THEN
          dprjc = dprjca
          kclo  = krej
        ENDIF
        IF (krej <  momlhd(nodrrej,nhnlev)) THEN
          IF ((leqsfc) .AND. (krjc == krej)) krjc = krjc + 1
          krej = krej + 1
        ELSE
          ltop = .TRUE.
        ENDIF
      ENDDO
      krej = krjc
!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'kclo2   ', kact, krej, klev, kclo                          &
!                            , omlbdy(nodrrej,MAX(kclo,1),nbtp)

!  If level 'kclo' is close enough to level 'kact', complement missing data
!  and replace data if override bit is higher or flag is lower.

      IF (      (ABS( omlbdy(nodrrej,MAX(kclo,1),nbtplv)                       &
                     -omlbdy(nodract,    kact   ,nbtplv)) <= rprlmt)           &
          .AND. (kclo > 0) .AND. (lrepla)) THEN

!       IF (yomlhd(nodract) (1:5) == '03743')     &
!         WRITE(0,*) 'kclo3   ', kact, krej, klev, kclo                        &
!                              , omlbdy(nodrrej,MAX(kclo,1),nbtp), rprlmt

        !   If level of rejected report is mandatory, replace pressure + height
        !   (provided that the newly assigned level is below (i.e. at higher
        !    pressure) than the subsequent level of the active report, or
        !    that autochecking is done)

        lmand  =   (rmod( omlbdy(nodrrej,kclo,nbtplv),10000.0_wp ) < epsy2)    &
              .OR. ( ABS( omlbdy(nodrrej,kclo,nbtplv)- 5000.0_wp ) < epsy)     &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 7000.0_wp ) < epsy)     &
              .OR. ( ABS( omlbdy(nodrrej,kclo,nbtplv)-15000.0_wp ) < epsy)     &
              .OR. ( ABS( omlbdy(nodrrej,kclo,nbtplv)-25000.0_wp ) < epsy)     &
              .OR. ( ABS( omlbdy(nodrrej,kclo,nbtplv)-85000.0_wp ) < epsy)     &
              .OR. ( ABS( omlbdy(nodrrej,kclo,nbtplv)-92500.0_wp ) < epsy)
        lsupoac = BTEST( momlbd(nodract,kact,nbtlsg), LS_SUPEROBS )
        lsuporj = BTEST( momlbd(nodrrej,kclo,nbtlsg), LS_SUPEROBS )
        lsgexch = .FALSE.
        IF (    (ABS( omlbdy(nodrrej,kclo,nbtplv)                              &
                     -omlbdy(nodract,kact,nbtplv)) > epsy)                     &
           .AND.(    ( omlbdy(nodrrej,     kclo  ,nbtplv)                      &
                      -omlbdy(nodract,MIN( kact+1, momlhd(nodract,nhnlev) )    &
                                                 ,nbtplv) > c05)               &
                 .OR.(ktypchk == 1))                                           &
           .AND.(    (    (.NOT. BTEST( momlbd(nodract,kact,nbterr), nvrz ))   &
                     .AND.(      BTEST( momlbd(nodrrej,kclo,nbterr), nvrz )))  &
                 .OR.(    (    (.NOT.BTEST(momlbd(nodract,kact,nbterr),nvrz))  &
                           .OR.(     BTEST(momlbd(nodrrej,kclo,nbterr),nvrz))) &
                     .AND.((lmand) .OR. (lsuporj))                             &
                     .AND.(.NOT. lsupoac)))) THEN
          ilen = 8 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+ 1) = ilen
            outbuf(nacout+ 2) = nfmt16
            outbuf(nacout+ 3) = NINT(omlbdy(nodract,kact,nbtp  )*rscale)
            outbuf(nacout+ 4) = NINT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
            outbuf(nacout+ 5) = imdi
            outbuf(nacout+ 6) = imdi
            outbuf(nacout+ 7) = imdi
            outbuf(nacout+ 8) = imdi
            IF (BTEST( momlbd(nodract,kact,nbterr), nvrz ))                    &
              outbuf(nacout+ 5) = NINT(omlbdy(nodract,kact,nbtzer)*rscale)
            IF (BTEST( momlbd(nodrrej,kclo,nbterr), nvrz ))                    &
              outbuf(nacout+ 6) = NINT(omlbdy(nodrrej,kclo,nbtzer)*rscale)
            IF (omlbdy(nodract,kact,nbtz  ) > rmdich)                          &
              outbuf(nacout+ 7) = NINT(omlbdy(nodract,kact,nbtz  )*rscale)
            IF (omlbdy(nodrrej,kclo,nbtz  ) > rmdich)                          &
              outbuf(nacout+ 8) = NINT(omlbdy(nodrrej,kclo,nbtz  )*rscale)
            DO icl = 1 , istrej
              outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF

          zmlbdy (klev,nbtp  ) = omlbdy(nodrrej,kclo,nbtp)
          zmlbdy (klev,nbtz  ) = omlbdy(nodrrej,kclo,nbtz)
          zmlbdy (klev,nbtzer) = omlbdy(nodrrej,kclo,nbtzer)
          zmlbdy (klev,nbtlop) = omlbdy(nodrrej,kclo,nbtlop)
          zmlbdy (klev,nbtplv) = omlbdy(nodrrej,kclo,nbtplv)
          CALL MVBITS( momlbd(nodrrej,kclo,nbterr), nvrz  , 1                  &
                     , iomlbd(        klev,nbterr), nvrz           )
          CALL MVBITS( momlbd(nodrrej,kclo,nbtflg), nvfzbp, nvfaoc             &
                     , iomlbd(        klev,nbtflg), nvfzbp         )
        ENDIF

        !   a: wind components

        nflaga = IBITS( momlbd(nodract,kact,nbtflg), nvfubp+nvfbps(1),nvfboc(1))
        nflagr = IBITS( momlbd(nodract,kclo,nbtflg), nvfubp+nvfbps(1),nvfboc(1))
        IF (     (      (.NOT. BTEST( momlbd(nodract,kact,nbterr), nvru ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvru )))    &
            .OR. (      (omlbdy(nodract,kact,nbtu  ) < rmdich)                 &
                  .AND. (omlbdy(nodrrej,kclo,nbtu  ) > rmdich))                &
            .OR. (      (      BTEST( momlbd(nodract,kact,nbterr), nvru ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvru ))     &
                  .AND. (nflagr == 0) .AND. (nflaga /= 0))) THEN
          ilen = 8 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+ 1) = ilen
            outbuf(nacout+ 2) = nfmt13
            outbuf(nacout+ 3) = NINT(omlbdy(nodract,kact,nbtp  )*rscale)
            outbuf(nacout+ 4) = NINT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
            outbuf(nacout+ 5) = imdi
            outbuf(nacout+ 6) = imdi
            IF (BTEST( momlbd(nodract,kact,nbterr), nvru ))                    &
              outbuf(nacout+ 5) = NINT(omlbdy(nodract,kact,nbtuer)*rscale)
            IF (BTEST( momlbd(nodrrej,kclo,nbterr), nvru ))                    &
              outbuf(nacout+ 6) = NINT(omlbdy(nodrrej,kclo,nbtuer)*rscale)
            outbuf(nacout+ 7) = nflaga
            outbuf(nacout+ 8) = nflagr
            DO icl = 1 , istrej
              outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF

          zmlbdy (klev,nbtu  ) = omlbdy (nodrrej,kclo,nbtu  )
          zmlbdy (klev,nbtv  ) = omlbdy (nodrrej,kclo,nbtv  )
          zmlbdy (klev,nbtuer) = omlbdy (nodrrej,kclo,nbtuer)
          zmlbdy (klev,nbtw  ) = omlbdy (nodrrej,kclo,nbtw  )
          zmlbdy (klev,nbtsnr) = omlbdy (nodrrej,kclo,nbtsnr)
          zmlbdy (klev,nbtuac) = omlbdy (nodrrej,kclo,nbtuac)
          CALL MVBITS( momlbd(nodrrej,kclo,nbterr), nvru  , 1                  &
                     , iomlbd(        klev,nbterr), nvru           )
          CALL MVBITS( momlbd(nodrrej,kclo,nbtflg), nvfubp, nvfaoc             &
                     , iomlbd(        klev,nbtflg), nvfubp         )
          iomlhd (     nhuexi) = 1
          lsgexch = .TRUE.
        ENDIF

        !   b: temperature

        nflaga = IBITS( momlbd(nodract,kact,nbtflg), nvftbp+nvfbps(1),nvfboc(1))
        nflagr = IBITS( momlbd(nodract,kclo,nbtflg), nvftbp+nvfbps(1),nvfboc(1))
        IF (     (      (.NOT. BTEST( momlbd(nodract,kact,nbterr), nvrt ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvrt )))    &
            .OR. (      (omlbdy(nodract,kact,nbtt  ) < rmdich)                 &
                  .AND. (omlbdy(nodrrej,kclo,nbtt  ) > rmdich))                &
            .OR. (      (      BTEST( momlbd(nodract,kact,nbterr), nvrt ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvrt ))     &
                  .AND. (nflagr == 0) .AND. (nflaga /= 0))) THEN
          ilen = 8 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+ 1) = ilen
            outbuf(nacout+ 2) = nfmt14
            outbuf(nacout+ 3) = NINT(omlbdy(nodract,kact,nbtp  )*rscale)
            outbuf(nacout+ 4) = NINT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
            outbuf(nacout+ 5) = imdi
            outbuf(nacout+ 6) = imdi
            IF (BTEST( momlbd(nodract,kact,nbterr), nvrt ))                    &
              outbuf(nacout+ 5) = NINT(omlbdy(nodract,kact,nbtter)*rscale)
            IF (BTEST( momlbd(nodrrej,kclo,nbterr), nvrt ))                    &
              outbuf(nacout+ 6) = NINT(omlbdy(nodrrej,kclo,nbtter)*rscale)
            outbuf(nacout+ 7) = nflaga
            outbuf(nacout+ 8) = nflagr
            DO icl = 1 , istrej
              outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF

          zmlbdy (klev,nbtt  ) = omlbdy (nodrrej,kclo,nbtt  )
          zmlbdy (klev,nbtter) = omlbdy (nodrrej,kclo,nbtter)
          zmlbdy (klev,nbtsnr) = omlbdy (nodrrej,kclo,nbtsnr)
          zmlbdy (klev,nbtuac) = omlbdy (nodrrej,kclo,nbtuac)
          CALL MVBITS( momlbd(nodrrej,kclo,nbterr), nvrt  , 1                  &
                     , iomlbd(        klev,nbterr), nvrt           )
          CALL MVBITS( momlbd(nodrrej,kclo,nbtflg), nvftbp, nvfaoc             &
                     , iomlbd(        klev,nbtflg), nvftbp         )
          iomlhd (     nhtexi) = 1
          lsgexch = .TRUE.
        ENDIF

        !   c: humidity

        nflaga = IBITS( momlbd(nodract,kact,nbtflg), nvfqbp+nvfbps(1),nvfboc(1))
        nflagr = IBITS( momlbd(nodract,kclo,nbtflg), nvfqbp+nvfbps(1),nvfboc(1))
        IF (     (      (.NOT. BTEST( momlbd(nodract,kact,nbterr), nvrq ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvrq )))    &
            .OR. (      (omlbdy(nodract,kact,nbtrh ) < rmdich)                 &
                  .AND. (omlbdy(nodrrej,kclo,nbtrh ) > rmdich))                &
            .OR. (      (      BTEST( momlbd(nodract,kact,nbterr), nvrq ))     &
                  .AND. (      BTEST( momlbd(nodrrej,kclo,nbterr), nvrq ))     &
                  .AND. (nflagr == 0) .AND. (nflaga /= 0))) THEN
          ilen = 8 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+ 1) = ilen
            outbuf(nacout+ 2) = nfmt15
            outbuf(nacout+ 3) = NINT(omlbdy(nodract,kact,nbtp  )*rscale)
            outbuf(nacout+ 4) = NINT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
            outbuf(nacout+ 5) = imdi
            outbuf(nacout+ 6) = imdi
            IF (BTEST( momlbd(nodract,kact,nbterr), nvrq ))                    &
              outbuf(nacout+ 5) = NINT(omlbdy(nodract,kact,nbtqer)*rscale)
            IF (BTEST( momlbd(nodrrej,kclo,nbterr), nvrq ))                    &
              outbuf(nacout+ 6) = NINT(omlbdy(nodrrej,kclo,nbtqer)*rscale)
            outbuf(nacout+ 7) = nflaga
            outbuf(nacout+ 8) = nflagr
            DO icl = 1 , istrej
              outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF

          zmlbdy (klev,nbtrh ) = omlbdy (nodrrej,kclo,nbtrh )
          zmlbdy (klev,nbtqer) = omlbdy (nodrrej,kclo,nbtqer)
          zmlbdy (klev,nbtdrh) = omlbdy (nodrrej,kclo,nbtdrh)
          CALL MVBITS( momlbd(nodrrej,kclo,nbterr), nvrq  , 1                  &
                     , iomlbd(        klev,nbterr), nvrq           )
          CALL MVBITS( momlbd(nodrrej,kclo,nbtflg), nvfqbp, nvfaoc             &
                     , iomlbd(        klev,nbtflg), nvfqbp         )
          iomlhd (     nhqexi) = 1
        ENDIF

        !   d: level identity (surface, tropopause, max wind id complement
        !                     standard level id or replace significant level id)

        ilvsg  =  momlbd(nodrrej,kclo,nbtlsg)
        IF (     (BTEST( ilvsg, LS_SURFACE ))                                  &
            .OR. (BTEST( ilvsg, LS_TROPO   ))                                  &
            .OR. (BTEST( ilvsg, LS_MAX     ))) THEN
          IF (BTEST( momlbd(nodract,kact,nbtlsg), LS_SUPEROBS )) THEN
            IF (lsgexch) iomlbd (klev,nbtlsg) = ilvsg
            IF (lsgexch) iomlbd (klev,nbtlid) = momlbd(nodrrej,kclo,nbtlid)
          ELSEIF (BTEST( momlbd(nodract,kact,nbtlsg), LS_STANDARD )) THEN
            iomlbd (klev,nbtlsg) = IOR( ilvsg , iomlbd(klev,nbtlsg) )
            iomlbd (klev,nbtlid) = IOR( momlbd(nodrrej,kclo,nbtlid)            &
                                      ,         iomlbd(klev,nbtlid) )
          ELSE
            iomlbd (klev,nbtlsg) = ilvsg
            iomlbd (klev,nbtlid) = momlbd(nodrrej,kclo,nbtlid)
          ENDIF
        ENDIF
      ENDIF

      IF ((ABS( omlbdy(nodrrej,MAX(kclo,1),nbtplv)                             &
               -omlbdy(nodract,kact,nbtplv)) <= rprlmt) .AND. (kclo > 0)) THEN
        IF (.NOT. losse_fg) THEN
          !  tolerances for being subset: 1.1 m/s, 0.4 K, 0.05 RH
          IF (     (      (omlbdy(nodrrej,kclo,nbtu ) >  rmdich)               &
                    .AND. (ABS( omlbdy(nodract,kact,nbtu )                     &
                               -omlbdy(nodrrej,kclo,nbtu ) ) > 1.1_wp))        &
              .OR. (      (omlbdy(nodrrej,kclo,nbtv ) >  rmdich)               &
                    .AND. (ABS( omlbdy(nodract,kact,nbtv )                     &
                               -omlbdy(nodrrej,kclo,nbtv ) ) > 1.1_wp))        &
              .OR. (      (omlbdy(nodrrej,kclo,nbtt ) >  rmdich)               &
                    .AND. (ABS( omlbdy(nodract,kact,nbtt )                     &
                               -omlbdy(nodrrej,kclo,nbtt ) ) > .4_wp))         &
              .OR. (      (omlbdy(nodrrej,kclo,nbtrh) >  rmdich)               &
                    .AND. (ABS( omlbdy(nodract,kact,nbtrh)                     &
                               -omlbdy(nodrrej,kclo,nbtrh) ) > .05_wp)))       &
            lsub =.FALSE.
        ELSE
          IF (MAX( ABS( ABS( omlbdy(nodract,kact,nbtu ) )                      &
                       -ABS( omlbdy(nodrrej,kclo,nbtu ) ) )                    &
                 , ABS( ABS( omlbdy(nodract,kact,nbtt ) )                      &
                       -ABS( omlbdy(nodrrej,kclo,nbtt ) ) )                    &
                 , ABS( ABS( omlbdy(nodract,kact,nbtrh) )                      &
                       -ABS( omlbdy(nodrrej,kclo,nbtrh) ) ) ) >= ABS(rmdich))  &
            lsub =.FALSE.
        ENDIF
      ENDIF

! update 'krej': if 'kact' is surface level then 'krej' must be above pressure
!                                                level of 'kact' (minus 0.5 Pa),
!                if 'kact' is not sfc level then 'krej' must be above pressure
!                                                level of 'kact' minus 'rprlmt'

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'kclo7   ',kact,krej, klev, kclo, omlbdy(nodrrej,krej,nbtp) &
!                            ,omlbdy(nodract,kact,nbtp)
      IF (.NOT. lactsurf)  pnexlmt = omlbdy(nodract,kact,nbtplv) - rprlmt
      IF (      lactsurf)  pnexlmt = omlbdy(nodract,kact,nbtplv) - c05
      DO WHILE ((omlbdy(nodrrej,krej,nbtplv) >= pnexlmt) .AND. (.NOT. ltop))
        IF (krej <  momlhd(nodrrej,nhnlev)) THEN
          krej = krej + 1
        ELSE
          ltop = .TRUE.
        ENDIF
!       IF (yomlhd(nodract) (1:5) == '03743')     &
!         WRITE(0,*) 'kclo7b  ',kact,krej,klev,kclo, omlbdy(nodrrej,krej,nbtp) &
!                              ,omlbdy(nodract,kact,nbtp)
      ENDDO

! 'while' loop: Get all the rejected levels which are more than 'rprlmt' [Pa],
!               but less than 'rdplmt' [Pa] above the active level 'kact',
!               and more than 'rprlmt' [Pa] below the level 'kact+1'.
!               From these levels, fill those observed quantities, which are not
!               present in either of the 2 active levels, into the temporary
!               field.
!               (If 'kact' is the surface level, then 'kact' should not have
!                an effect on the treatment of rejected levels further above.
!                This is accomplished by jumping over this while loop.)
!        (Note: This may be important if e.g. a TEMP report containing only
!               temperature and humidity and a PILOT report have been derived
!               and disseminated from one and the same radiosonde ascent.)

      pnexlmt  =  c0
      IF (kact < momlhd(nodract,nhnlev))                                       &
        pnexlmt  =  omlbdy(nodract,kact+1,nbtplv) + rprlmt

!     IF (yomlhd(nodract) (1:5) == '03743')     &
!       WRITE(0,*) 'kclo8   ',kact,krej, klev, kclo, omlbdy(nodrrej,krej,nbtp) &
!                            ,omlbdy(nodract,kact,nbtp), pnexlmt
      DO WHILE (      (   omlbdy(nodrrej,krej,nbtplv)                          &
                       >= omlbdy(nodract,kact,nbtplv)-rdplmt)                  &
                .AND. (   omlbdy(nodrrej,krej,nbtplv) >= pnexlmt)              &
                .AND. (.NOT. ltop) .AND. (klev < 2*maxmlv)                     &
                .AND. (.NOT. lactsurf))
        lmand = .FALSE.
        laddu = .FALSE.
        laddt = .FALSE.
        laddq = .FALSE.
        IF (lrepla) THEN
          lmand =   (rmod( omlbdy(nodrrej,krej,nbtplv),10000._wp ) < epsy2)    &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 5000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)- 7000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-15000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-25000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-85000._wp ) < epsy)     &
               .OR. ( ABS( omlbdy(nodrrej,krej,nbtplv)-92500._wp ) < epsy)
          laddu =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvru)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvru))
          laddt =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvrt)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvrt))
          laddq =    ((.NOT.BTEST(momlbd(nodract,kact,nbterr),nvrq)).OR.lmand) &
                .AND. (     BTEST(momlbd(nodrrej,krej,nbterr),nvrq))
        ENDIF
        IF ((lrepla) .AND. (kact < momlhd(nodract,nhnlev))) THEN
          IF (   omlbdy(nodrrej,krej  ,nbtplv)                                 &
              <= omlbdy(nodract,kact+1,nbtplv)+rdplmt) THEN
            IF (laddu) laddu= ( .NOT.BTEST(momlbd(nodract,kact+1,nbterr),nvru) &
                             .OR. lmand)
            IF (laddt) laddt= ( .NOT.BTEST(momlbd(nodract,kact+1,nbterr),nvrt) &
                             .OR. lmand)
            IF (laddq) laddq= ( .NOT.BTEST(momlbd(nodract,kact+1,nbterr),nvrq) &
                             .OR. lmand)
          ENDIF
        ENDIF
        IF (      (BTEST( momlbd(nodrrej,krej,nbtlsg), LS_SURFACE ))           &
            .AND. (momlbd(nodrrej,krej,nbtlsg) /= imdi)) THEN
          laddu = .TRUE.
          laddt = .TRUE.
          laddq = .TRUE.
        ENDIF
!       IF (yomlhd(nodract) (1:5) == '03743')     &
!         WRITE(0,*) 'kclo6   ',kact, krej, klev, kclo, lmand, laddu,laddt,laddq
        IF ((laddu) .OR. (laddt) .OR. (laddq)) THEN
          klev = klev + 1
          zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
          iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
          IF (.NOT. laddu) zmlbdy (klev,nbtu  ) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtv  ) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtuer) = rmdi
          IF (.NOT. laddu) zmlbdy (klev,nbtw  ) = rmdi
          IF (.NOT. laddu) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvru )
          IF (.NOT. laddu) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvfubp)
          IF (.NOT. laddt) zmlbdy (klev,nbtt  ) = rmdi
          IF (.NOT. laddt) zmlbdy (klev,nbtter) = rmdi
          IF (.NOT. laddt) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvrt )
          IF (.NOT. laddt) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvftbp)
          IF (.NOT. laddq) zmlbdy (klev,nbtrh ) = rmdi
          IF (.NOT. laddq) zmlbdy (klev,nbtqer) = rmdi
          IF (.NOT. laddq) zmlbdy (klev,nbtdrh) = rmdi
          IF (.NOT. laddq) iomlbd (klev,nbterr) = IBCLR ( iomlbd(klev,nbterr)  &
                                                        , nvrq )
          IF (.NOT. laddq) CALL MVBITS( 0,0,nvfaoc, iomlbd(klev,nbtflg), nvfqbp)
          IF (BTEST( iomlbd(klev,nbterr), nvru ))  iomlhd (nhuexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrt ))  iomlhd (nhtexi) = 1
          IF (BTEST( iomlbd(klev,nbterr), nvrq ))  iomlhd (nhqexi) = 1
        ENDIF
        IF (krej <  momlhd(nodrrej,nhnlev)) THEN
          krej = krej + 1
        ELSE
          ltop = .TRUE.
        ENDIF
        ! set 'lsub' to .false. unless
        !   -  a mandatory level is added to a superobbed report, or
        !   -  only the top level of a superobbed redundant report (with nlev>2)
        !      is absent in (but within the range of) the active profile
        IF (.NOT. (     (lsupoba .AND. lmand)                                  &
                   .OR. (lsupobr .AND. ltop .AND. (krej > 2))))  lsub = .FALSE.
      ENDDO

!  Close the loop over levels of body of active report

    ENDDO Levels_acr
!   ~~~~~~~~~~~~~~~~

!  Get all rejected levels which are more than 'rdplmt' [pa] above
!  the highest active level. Fill those levels into the temporary field.

    DO WHILE ((.NOT. ltop) .AND. (klev < 2*maxmlv))
      IF (lrepla) THEN
        klev = klev + 1
        zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
        iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
      ENDIF
      IF (krej <  momlhd(nodrrej,nhnlev)) THEN
        krej = krej + 1
      ELSE
        ltop = .TRUE.
      ENDIF
      ! set 'lsub' to .false. unless a level is added on top of a radiosonde rep
      IF (.NOT. lraso)  lsub = .FALSE.
    ENDDO

!  determine number of levels

    iomlhd (nhnlev) = klev

!  prepare discarding any levels below the surface level
    isdel = 0
    DO iiv = klev , 1 , -1
      IF (      (BTEST( iomlbd(iiv,nbtlsg), LS_SURFACE ))                      &
          .AND. (iomlbd(iiv,nbtlsg) /= imdi))  isdel = iiv - 1
    ENDDO
    iomlhd (nhnlev) = iomlhd(nhnlev) - isdel

!   DO iiv = 1 , iomlhd(nhnlev) + isdel
!     WRITE(0,*) 'nodract2 ', ystid, ktypchk, iomlhd(nhnlev), isdel, iiv       &
!                                  , zmlbdy(iiv,nbtp), iomlbd(iiv,nbtlsg)
!   ENDDO

!  limit nlev = momlhd(.,nhnlev) to maxmlv
!    (to prevent crashes due to access beyond array limits !!!)
    iomlhd (nhnlev) = MIN( iomlhd(nhnlev) , maxmlv )

! flag the redundant report and store the active report (from the temporary
! array) in the ODR; if the active report contains any new data (i.e. from
! the report with index 'nmlob'), it has to be stored with index 'nmlob'
! in order to facilitate its further processing (e.g. mult_quality_check)
! -------------------------------------------------------------------------

    IF ((ktypchk == 1) .OR. (nodract == nmlob)                                 &
                       .OR. ((lrepla) .AND. (.NOT. lsub))) THEN
      nmlza  =  nmlob
    ELSE
      nmlza  =  nmlo2
    ENDIF

!   IF (yomlhd(nodract) (1:5) == '03743')     &
!     WRITE(0,*) 'nodrrej2 ', ystid, ktypchk, nmlza, nmlob, nmlo2, nodrrej     &
!            , momlhd(nmlza,nhnlev), momlhd(nmlo2,nhnlev), iomlhd(nhnlev)      &
!            , lsub, lrepla
!            , momlhd(nmlza,nhnlev), momlhd(nmlo2,nhnlev), iomlhd(nhnlev), isdel

!  if necessary, move redundant report within the ODR to index 'nmlo2'
    IF ((ktypchk == 2) .AND. (nodrrej == nmlob) .AND. (nmlza == nmlob)) THEN
      omlbdy (nmlo2 ,1:maxmlv,1:mxrbdy) = omlbdy (nodrrej,1:maxmlv,1:mxrbdy)
      momlbd (nmlo2 ,1:maxmlv,1:mxrbdf) = momlbd (nodrrej,1:maxmlv,1:mxrbdf)
      omlhed (nmlo2          ,1:mxrhed) = omlhed (nodrrej         ,1:mxrhed)
      momlhd (nmlo2          ,1:mxrhdf) = momlhd (nodrrej         ,1:mxrhdf)
      yomlhd (nmlo2)                    = yomlhd (nodrrej)
    ENDIF

!  copy temporary field into ODR (for active report)
    omlbdy (nmlza,1:maxmlv,1:mxrbdy) = zmlbdy (1+isdel:maxmlv+isdel,1:mxrbdy)
    momlbd (nmlza,1:maxmlv,1:mxrbdf) = iomlbd (1+isdel:maxmlv+isdel,1:mxrbdf)
    omlhed (nmlza         ,1:mxrhed) = zmlhed (                     1:mxrhed)
    momlhd (nmlza         ,1:mxrhdf) = iomlhd (                     1:mxrhdf)
    yomlhd (nmlza)                   = ystidz

    ! for safety: remove duplicate levels in active report
    nlev = momlhd(nmlza,nhnlev)
    CALL obs_rm_duplicate_levels ( nlev, mxrbdy, mxrbdf, nbtplv, rmdich        &
                                 , omlbdy(nmlza,1:nlev,:)                      &
                                 , momlbd(nmlza,1:nlev,:), nrejlv )
  ! ============================
    momlhd (nmlza,nhnlev) = nlev - nrejlv

  ENDDO   !   DO ktypchk = 1 , 1 + MIN( nmlo2, 1 )
! ~~~~~

! IF (yomlhd(nodract) (1:5) == '03743')     &
!   WRITE(0,*) 'nodrrej3 ', ystid, ktypchk, nmlza, nmlob, nmlo2, nodrrej       &
!                         , momlhd(nmlza,nhnlev), momlhd(nmlo2,nhnlev)

! if a redundant report is found (not in auto-checking)
  IF (nmlo2 > 0) THEN
    IF (nmlza == nmlob)  nmlrj = nmlo2
    IF (nmlza /= nmlob)  nmlrj = nmlob
    IF ((lverpas) .AND. (.NOT. lsub) .AND. (iactio == 1)) THEN
! ordinary flagging of redundant report
      momlhd (nmlrj ,nhpass)  = 2
      momlhd (nmlrj ,nhschr)  = IBSET( momlhd(nmlrj,nhschr) , nvrdbp )
      momlhd (nmlrj ,nhschr)  = IBSET( momlhd(nmlrj,nhschr) , nvpsbp )
      momlhd (nmlrj ,nhflag)  = IBSET( momlhd(nmlrj,nhflag) , FL_REDUNDANT )
    ELSE
! redundant report will be discarded in 'obs_cdf_del_old_rep'
      momlhd (nmlrj ,nhpass)  = -1
    ENDIF

!   IF (yomlhd(nodract) (1:5) == '03743')     &
!     WRITE(0,*) 'drrej4 ', ystid, iactio, ktypchk, nmlza, nmlob, nmlo2, nmlrj &
!                         , momlhd(nmlrj,nhpass), lsub, npassiv, npassv2

! update statistics and report events, and prepare control output
! ---------------------------------------------------------------

    kobtypr = momlhd(nmlrj,nhobtp)
    kcdtypr = momlhd(nmlrj,nhcode)
    icma    = i_cma ( kobtypr , kcdtypr )
!             =====
    IF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
      cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      IF (iactio == 1)  neventr (neredn,icma) = neventr(neredn,icma) + 1
      IF (iactio == 2)  neventr (neslps,icma) = neventr(neslps,icma) + 1
    ELSE
      cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
    ENDIF
    IF ((lverpas) .AND. (.NOT. lsub) .AND. (iactio == 1)) THEN
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
    ELSE
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF

    ilen = 12 + 2*istrej + (momlhd(nmlza,nhnlev) + momlhd(nmlrj,nhnlev)) *11
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt11
      outbuf(nacout+ 3) = kcdtypr
      outbuf(nacout+ 4) = momlhd(nmlza,nhitot)
      outbuf(nacout+ 5) = momlhd(nmlza,nhjtot)
      outbuf(nacout+ 6) = NINT(omlhed(nmlza,nhilon)*rscale)
      outbuf(nacout+ 7) = NINT(omlhed(nmlza,nhjlat)*rscale)
      IF (omlhed(nmlza,nhalt) > rmdich) THEN
        outbuf(nacout+8) = NINT(omlhed(nmlza,nhalt))
      ELSE
        outbuf(nacout+8) = imdi
      ENDIF
      outbuf(nacout+ 9) = NINT(omlhed(nmlza,nhtime)*rscale)
!     IF (lraso) outbuf(nacout+ 9) = NINT(omlhed(nmlza,nhsynt)*rscale)
      outbuf(nacout+10) = momlhd(nmlza,nhcode)
      outbuf(nacout+11) = momlhd(nmlza,nhnlev)
      outbuf(nacout+12) = momlhd(nmlrj,nhnlev)
      DO icl = 1 , istrej
        outbuf(nacout+12       +icl) = ICHAR( yomlhd(nmlrj)(icl:icl) )
        outbuf(nacout+12+istrej+icl) = ICHAR( yomlhd(nmlza)(icl:icl) )
      ENDDO
      nacout = nacout + 12 + 2*istrej
      DO klvi = 1, momlhd(nmlrj,nhnlev) + momlhd(nmlza,nhnlev)
        DO iiv = 1, 10
          outbuf (nacout+ iiv) = imdi
        ENDDO
        nactpr = nmlrj
        klev   = klvi
        IF (klvi > momlhd(nmlrj,nhnlev)) THEN
          nactpr = nmlza
          klev   = klvi - momlhd(nmlrj,nhnlev)
        ENDIF
        IF (omlbdy(nactpr,klev,nbtu  ) > rmdich)                               &
          outbuf (nacout+ 1) = NINT(omlbdy(nactpr,klev,nbtu  )*rscale)
        IF (omlbdy(nactpr,klev,nbtv  ) > rmdich)                               &
          outbuf (nacout+ 2) = NINT(omlbdy(nactpr,klev,nbtv  )*rscale)
        IF (omlbdy(nactpr,klev,nbtt  ) > rmdich)                               &
          outbuf (nacout+ 3) = NINT(omlbdy(nactpr,klev,nbtt  )*rscale)
        IF (omlbdy(nactpr,klev,nbtrh ) > rmdich)                               &
          outbuf (nacout+ 4) = NINT(omlbdy(nactpr,klev,nbtrh )*rscale)
        IF (omlbdy(nactpr,klev,nbtp  ) > rmdich)                               &
          outbuf (nacout+ 5) = NINT(omlbdy(nactpr,klev,nbtp  )*rscale)
        IF (omlbdy(nactpr,klev,nbtz  ) > rmdich)                               &
          outbuf (nacout+ 6) = NINT(omlbdy(nactpr,klev,nbtz  )*rscale)
        IF (BTEST( momlbd(nactpr,klev,nbterr), nvru ))                         &
          outbuf (nacout+ 7) = NINT(omlbdy(nactpr,klev,nbtuer)*rscale)
        IF (BTEST( momlbd(nactpr,klev,nbterr), nvrt ))                         &
          outbuf (nacout+ 8) = NINT(omlbdy(nactpr,klev,nbtter)*rscale)
        IF (BTEST( momlbd(nactpr,klev,nbterr), nvrq ))                         &
          outbuf (nacout+ 9) = NINT(omlbdy(nactpr,klev,nbtqer)*rscale)
        IF (BTEST( momlbd(nactpr,klev,nbterr), nvrz ))                         &
          outbuf (nacout+10) = NINT(omlbdy(nactpr,klev,nbtzer)*rscale)
        outbuf(nacout+11) = momlbd(nactpr,klev,nbtlid)
        nacout  = nacout + 11
      ENDDO
    ENDIF  ! lprodr
  ENDIF    ! (nmlo2 > 0)

 ENDDO loop_over_new_reports
!~~~~~~~~~~~~~~~~~~~~~~~~~~~
ENDDO      ! iactio
!~~~~

!-------------------------------------------------------------------------------
! Section 3.0  Redundancy check for GPS (IWV) reports
!-------------------------------------------------------------------------------

#ifdef __COSMO__
  DO ngpob = 1 , nngpnew
    igplat = mogphd (ngpob+ntotgp-nngpnew,nhjtot)
    cosolat2 (ngpob) = COS( (r_startlat_tot+(igplat-1)*r_dlat) * r_degrad )
    cosolat2 (ngpob) = cosolat2(ngpob) * cosolat2(ngpob)
  ENDDO
#endif

  IF (nngpnew > 0) THEN
    DO iiv = 1 , mxgpc
      IF (i_gpscen(iiv) >= 0)  nactgpc = iiv
    ENDDO
  ENDIF

DO ngpob = ntotgp - nngpnew + 1 , ntotgp

! get some information of current report
! --------------------------------------

  npassiv = mogphd (ngpob,nhpass)
  ksubcen = mogphd (ngpob,nhcent) /1000
  obtime  = ogphed (ngpob,nhtime)
  igplon  = mogphd (ngpob,nhitot)
  igplat  = mogphd (ngpob,nhjtot)
  obbalt  = ogphed (ngpob,nhalt )
!  Preset thresholds for 4-d colocation
  rtmlmt = rtmlim
  rhzlmt = rhzlim * rhzlim
  rvtlmt = rvtlim


! check if current report meets the redundancy criteria with another report
! -------------------------------------------------------------------------
! (do this with a vectorisable loop)

  nchkend = ngpob - 1
  IF (npassiv == -1) nchkend = 0

  DO ngplo = 1 , nchkend
    igplon2 = mogphd (ngplo,nhitot)
    igplat2 = mogphd (ngplo,nhjtot)

!  check for horiz/vert. (quasi-)colocation, simultaneity, identical station ID
#ifdef __COSMO__
!   rhzob  =    cosolat2(ngpob+ntotgp-nngpnew)                                 &
    rhzob  =    cosolat2(ngpob-(ntotgp-nngpnew))                               &
               *((igplon2-igplon) *r_dlon) *((igplon2-igplon) *r_dlon)         &
             +  ((igplat2-igplat) *r_dlat) *((igplat2-igplat) *r_dlat)
#else
    rhzob  = r2d **2 * sum ((  gen_grid% xnglob (igplon ,igplat ,:,1)&
                             - gen_grid% xnglob (igplon2,igplat2,:,1)) **2)
#endif
    lcoloc =      (      rhzob * rdegkm2                  <=  rhzlmt )         &
            .AND. ( ABS( ogphed(ngplo,nhalt ) - obbalt )  <=  rvtlmt )         &
            .AND. ( ABS( ogphed(ngplo,nhtime) - obtime )  <=  rtmlmt )
    ! reproducibility irrespective of domain decomposition ensured by allowing
    ! for redundancy only between reports assigned to the same grid point
    IF (lredn_repro)  lcoloc = lcoloc .AND. (igplon2 == igplon)                &
                                      .AND. (igplat2 == igplat)
!   lredn (ngplo) =       (lcoloc) .AND. (yogphd(ngplo) == yogphd(ngpob))      &
!                   .AND. (mogphd(ngplo,nhpass) /= -1)
    ! reports from same station may have different station id's
    ! if processed from different centres
    lredn (ngplo) =       (lcoloc) .AND. (mogphd(ngplo,nhpass) /= -1)
  ENDDO
  ngpo2 = 0
! prefer to find redundancy between active reports (or between passive reports)
  DO ngplo = 1 , nchkend
    IF (      (lredn(ngplo)) .AND. (ngpo2 == 0)                                &
        .AND. (mogphd(ngplo,nhpass) == npassiv))  ngpo2 = ngplo
  ENDDO
! otherwise find redundancy between active and passive report
  DO ngplo = 1 , nchkend
    IF (      (lredn(ngplo)) .AND. (ngpo2 == 0))  ngpo2 = ngplo
  ENDDO

  IF (ngpo2 > 0) THEN

! decide which report is redundant
! --------------------------------
    npassv2  =  mogphd(ngpo2,nhpass)
    ksubce2  =  mogphd(ngpo2,nhcent) /1000

    ! get indices of reported processing centers within 'i_gpscen'
    ! (the report with the smaller index will be prefered)
    ixcen = nactgpc + 1
    ixce2 = nactgpc + 1
    DO iiv = nactgpc , 1 , -1
      IF (ksubcen == i_gpscen(iiv))  ixcen = iiv
      IF (ksubce2 == i_gpscen(iiv))  ixce2 = iiv
    ENDDO

! prefer active report over passive report
    IF ((mogphd(ngpob,nhpass) < 2) .AND. (mogphd(ngpo2,nhpass) == 2)) THEN
      nodract  =  ngpob
    ELSEIF ((mogphd(ngpob,nhpass) == 2) .AND. (mogphd(ngpo2,nhpass) < 2)) THEN
      nodract  =  ngpo2
    ELSEIF (      (      BTEST( mogpbd(ngpob,nbgerr), nvriwv ))                &
            .AND. (.NOT. BTEST( mogpbd(ngpo2,nbgerr), nvriwv ))) THEN
      nodract  =  ngpob
    ELSEIF (      (.NOT. BTEST( mogpbd(ngpob,nbgerr), nvriwv ))                &
            .AND. (      BTEST( mogpbd(ngpo2,nbgerr), nvriwv ))) THEN
      nodract  =  ngpo2
! prefer processing centre according to the order in 'i_gpscen'
! (BKG, then GFZ, LPT, GOP, SGN, KNMI, NGAA, ROB, MET ...)
    ELSEIF (ixcen < ixce2) THEN
      nodract  =  ngpob
    ELSEIF (ixce2 < ixcen) THEN
      nodract  =  ngpo2
! prefer higher update sequence number
    ELSEIF (mogphd(ngpob,nhcorr) > mogphd(ngpo2,nhcorr)) THEN
      nodract  =  ngpob
    ELSEIF (mogphd(ngpob,nhcorr) < mogphd(ngpo2,nhcorr)) THEN
      nodract  =  ngpo2
! prefer smaller observation error, if obs error exist for both reports
    ELSEIF (ogpbdy(ngpo2,nbgtze) > MAX(ogpbdy(ngpob,nbgtze),c0)) THEN
      nodract  =  ngpob
    ELSEIF (ogpbdy(ngpob,nbgtze) > MAX(ogpbdy(ngpo2,nbgtze),c0)) THEN
      nodract  =  ngpo2
    ELSE
      nodract  =  ngpo2
    ENDIF
    nodrrej  =  ngpob + ngpo2 - nodract

! flag redundant report
! ---------------------

    IF (lverpas) THEN
! ordinary flagging (if used for quality control of IWV or verif of passive rep)
      mogphd (nodrrej,nhpass)  = 2
      mogphd (nodrrej,nhschr)  = IBSET( mogphd(nodrrej,nhschr) , nvrdbp )
      mogphd (nodrrej,nhschr)  = IBSET( mogphd(nodrrej,nhschr) , nvpsbp )
      mogphd (nodrrej,nhflag)  = IBSET( mogphd(nodrrej,nhflag) , FL_REDUNDANT )
    ELSE
! this report will be discarded in 'obs_cdf_del_old_rep'
      mogphd (nodrrej,nhpass)  = -1
    ENDIF

! update statistics and report events, and prepare control output
! ---------------------------------------------------------------

    kobtypr = mogphd(nodrrej,nhobtp)
    kcdtypr = mogphd(nodrrej,nhcode)
    icma    = i_cma ( kobtypr , kcdtypr )
!             =====
    IF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
      cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
      neventr (neredn,icma) = neventr(neredn,icma) + 1
    ELSE
      cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
    ENDIF
    IF (lverpas) THEN
      cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
    ELSE
      cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
    ENDIF

    ilen = 15 + 2*istrej
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt10
      outbuf(nacout+ 3) = kcdtypr
      DO iiv = 1, 5
        outbuf(nacout+ 3+iiv) = imdi
      ENDDO
      IF (ogpbdy(nodract,nbgiwv) > rmdich)                                     &
        outbuf(nacout+ 4) = INT(ogpbdy(nodract,nbgiwv)*rscale)
      IF (ogpbdy(nodract,nbgzwd) > rmdich)                                     &
        outbuf(nacout+ 5) = INT(ogpbdy(nodract,nbgzwd)*rscale)
      IF (ogpbdy(nodract,nbgzpd) > rmdich)                                     &
        outbuf(nacout+ 6) = INT(ogpbdy(nodract,nbgzpd)*rscale)
      IF (ogpbdy(nodract,nbgtze) > rmdich)                                     &
        outbuf(nacout+ 7) = INT(ogpbdy(nodract,nbgtze)*rscale)
      IF (ogpbdy(nodract,nbgp  ) > rmdich)                                     &
        outbuf(nacout+ 8) = INT(ogpbdy(nodract,nbgp  )*rscale)
      outbuf(nacout+ 9) = mogphd(nodract,nhcode)
      outbuf(nacout+10) = INT(ogphed(nodract,nhilon)*rscale)
      outbuf(nacout+11) = INT(ogphed(nodract,nhjlat)*rscale)
      outbuf(nacout+12) = 1 + MAX( 0 , MIN( 1, nodract - nodrrej ) )
!     outbuf(nacout+12) = nactpr
      outbuf(nacout+13) = 0
      IF (mogphd(nodract,nhcorr) /= mogphd(nodrrej,nhcorr))                    &
        outbuf(nacout+13) = 1
      outbuf(nacout+14) = 0
      outbuf(nacout+15) = INT(ogphed(nodract,nhtime)*rscale)
      DO icl = 1 , istrej
        outbuf(nacout+15       +icl) = ICHAR( yogphd(nodrrej) (icl:icl) )
        outbuf(nacout+15+istrej+icl) = ICHAR( yogphd(nodract) (icl:icl) )
      ENDDO
      nacout  = nacout + ilen
    ENDIF  ! lprodr
  ENDIF    ! (ngpo2 > 0)
ENDDO

! clean up
! --------

#ifdef __COSMO__
  DEALLOCATE ( cosolat2 , STAT=istat )
#endif
  DEALLOCATE ( lredn    , STAT=istat )

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_redundancy
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_redundancy


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for deleting 'old' reports in the ODR
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_del_old_rep ( nodrold )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" deletes reports
!   which are not used any more (for nudging or verification), and it puts
!   forward the current reports in the data records ODR.
!   The procedure also updates the number of (remaining) 'old' reports in the
!   ODR (i.e. reports which have been read before this timestep) and the total
!   number of reports stored in the ODR.
!
! Method:
!   Reports have been previously been marked (by setting 'mo??hd(,.nhpass) =-1')
!   to be deleted for two reasons:
!   - either they are too old to be used any more,
!   - or they are flagged (ODR status entry 'nhpass' set to -1) in the
!     observation pre-processing within the current module, e.g. in the
!     redundancy checking.
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

  INTEGER        , INTENT (INOUT)  ::  &
    nodrold  (3)     ! number of multi-level / single-level / GPS reports
                     ! before having read new data at the current timestep

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    nmlob , nsgob ,& ! loop indices
    ngpob         ,& ! loop indices
    narest        ,& ! ODR report index of first missing 'report'
    ndiff         ,& ! total number of deleted reports
    ndiffo           ! number of old (read at past timesteps) deleted reports

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_del_old_rep
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Multi-level reports
!-------------------------------------------------------------------------------

  ndiff  = 0
  ndiffo = 0
  DO nmlob = 1 , ntotml
    IF (momlhd(nmlob,nhpass) == -1) THEN
      ndiff = ndiff + 1
      IF (nmlob <= nodrold(1)) ndiffo = ndiffo + 1
    ELSE
      omlbdy (nmlob-ndiff,1:maxmlv,1:mxrbdy) = omlbdy (nmlob,1:maxmlv,1:mxrbdy)
      momlbd (nmlob-ndiff,1:maxmlv,1:mxrbdf) = momlbd (nmlob,1:maxmlv,1:mxrbdf)
      omlhed (nmlob-ndiff,1:mxrhed)          = omlhed (nmlob,1:mxrhed)
      momlhd (nmlob-ndiff,1:mxrhdf)          = momlhd (nmlob,1:mxrhdf)
      yomlhd (nmlob-ndiff)                   = yomlhd (nmlob)
    ENDIF
  ENDDO

  narest  = ntotml  - ndiff + 1
  omlbdy (narest:ntotml,1:maxmlv,1:mxrbdy) = rmdi
  momlbd (narest:ntotml,1:maxmlv,1:mxrbdf) = imdi
  omlhed (narest:ntotml,1:mxrhed)          = rmdi
  momlhd (narest:ntotml,1:mxrhdf)          = imdi
  yomlhd (narest:ntotml)                   = '        '

  ntotml     = ntotml     - ndiff
  nodrold(1) = nodrold(1) - ndiffo

!-------------------------------------------------------------------------------
!  Section 2: Single-level reports
!-------------------------------------------------------------------------------

  ndiff  = 0
  ndiffo = 0
  DO  nsgob = 1 , ntotsg
    IF (mosghd(nsgob,nhpass) == -1) THEN
      ndiff = ndiff + 1
      IF (nsgob <= nodrold(2)) ndiffo = ndiffo + 1
    ELSE
      osgbdy (nsgob-ndiff,1:mxsbdy) = osgbdy (nsgob,1:mxsbdy)
      mosgbd (nsgob-ndiff,1:mxsbdf) = mosgbd (nsgob,1:mxsbdf)
      osghed (nsgob-ndiff,1:mxshed) = osghed (nsgob,1:mxshed)
      mosghd (nsgob-ndiff,1:mxshdf) = mosghd (nsgob,1:mxshdf)
      yosghd (nsgob-ndiff)          = yosghd (nsgob)
    ENDIF
  ENDDO
  narest  = ntotsg  - ndiff + 1
  osgbdy (narest:ntotsg,1:mxsbdy) = rmdi
  mosgbd (narest:ntotsg,1:mxsbdf) = imdi
  osghed (narest:ntotsg,1:mxshed) = rmdi
  mosghd (narest:ntotsg,1:mxshdf) = imdi
  yosghd (narest:ntotsg)          = '        '

  ntotsg     = ntotsg     - ndiff
  nodrold(2) = nodrold(2) - ndiffo

!-------------------------------------------------------------------------------
!  Section 3: GPS reports
!-------------------------------------------------------------------------------

  ndiff  = 0
  ndiffo = 0
  DO  ngpob = 1 , ntotgp
    IF (mogphd(ngpob,nhpass) == -1) THEN
      ndiff = ndiff + 1
      IF (ngpob <= nodrold(3)) ndiffo = ndiffo + 1
    ELSE
      ogpbdy (ngpob-ndiff,1:mxgbdy) = ogpbdy (ngpob,1:mxgbdy)
      mogpbd (ngpob-ndiff,1:mxgbdf) = mogpbd (ngpob,1:mxgbdf)
      ogphed (ngpob-ndiff,1:mxghed) = ogphed (ngpob,1:mxghed)
      mogphd (ngpob-ndiff,1:mxghdf) = mogphd (ngpob,1:mxghdf)
      yogphd (ngpob-ndiff)          = yogphd (ngpob)
    ENDIF
  ENDDO
  narest  = ntotgp  - ndiff + 1
  ogpbdy (narest:ntotgp,1:mxgbdy) = rmdi
  mogpbd (narest:ntotgp,1:mxgbdf) = imdi
  ogphed (narest:ntotgp,1:mxghed) = rmdi
  mogphd (narest:ntotgp,1:mxghdf) = imdi
  yogphd (narest:ntotgp)          = '        '

  ntotgp     = ntotgp     - ndiff
  nodrold(3) = nodrold(3) - ndiffo

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_del_old_rep
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_del_old_rep


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for douplicating tower reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_tower_copy ( dhosag, ntotmlo, ntwex, ilstidtw               &
                              , htwex, ivtwex, ytwex , ntotml, nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" douplicates multi-level
!   tower reports (i.e. reports stored in multi-level ODR).
!   The purpose is the following:
!   - The obs in the original reports have the mandatory-level flag set, and the
!     observation operators will compute the model equivalents at the correct
!     height of sensor above ground (possibly incl. a height correction).
!   - The obs in the copied reports do not have the mandatory-level flag set,
!     and the observation operators will compute the model equivalents at the
!     correct altitude (height above mean sea level) of the observations. 
!   Flags are set in the original and / or the duplicated report depending on
!   subr. argument 'dhosag' which is set by namelist and determines whether the
!   whole report or specific levels should be processed at reported height of 
!   sensor above ground or at correct altitude.
!   Whole reports are flagged with FL_THIN, and single levels with 'nvfsbp'.
!   Later on, these reports / observations are set passive (not: rejected)
!   in the conversion to the DACE data structures (for 'nvfsbp', this needs
!   special treatment in subr. 'check_store_conv' in mo_cosmo_conv.f90).
!
! Method:
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

  INTEGER        , INTENT (IN)     ::  &
    ntotmlo     ,& ! number of (multi-level) reports before having read
                   !   new data at the current timestep
    ntwex       ,& ! number of exceptions for tower processing
    ilstidtw    ,& ! max. length of tower station ID's
    ivtwex(MAX(1,ntwex)) ! variable indicators related to exception levels:
                         !   (O: none) 1: wind; 2: T + humidity; 3: all

  REAL (KIND=wp) , INTENT (IN)     ::  &
    dhosag      ,& ! limit in terms of height of sensor above ground [m]
                   !   below which model equivalents are computed at the
                   !               observed height of sensor above ground, and
                   !   above which model equivalents are computed at the
                   !               observed height (altitude)
    htwex (MAX(1,ntwex)) ! sensor heights a. ground [m] of exception levels;
                         !   if = 0 then apply to all levels (of that tower)

  CHARACTER (LEN=ilstidtw) , INTENT (IN)  :: &
    ytwex (MAX(1,ntwex)) ! tower station ID's related to exception levels


  INTEGER        , INTENT (INOUT)  ::  &
    ntotml      ,& ! total number of multi-level reports
    nexceed        ! number of (multi-level) reports in excess of ODR array size

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    itw_s (MAX(1,ntotml-ntotmlo)),& ! report indices of original tower reports
    itw_t (MAX(1,ntotml-ntotmlo)),& ! report indices of duplicated tower reports
    iob , io2   ,& ! loop indices over multi-level reports (in ODR)
    itw         ,& ! loop index   over new tower reports
    iex         ,& ! loop index   over exceptions from 'dhosag' rule
    ilev        ,& ! loop index   over vertical levels
    nadd        ,& ! number of duplicated reports
    nexcetw     ,& ! number of tower reports in excess of ODR array size
    kobtyp      ,& ! observation type
    kcdtyp      ,& ! observation code type
    icma           ! indices for statistics ('cma', events)

  LOGICAL        ::  &
    lrepass (MAX(1,ntotml-ntotmlo)) ! set report to passive if 'dhosag'
                                    !   criterion is not met (i.e. for all obs)

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_tower_copy
!-------------------------------------------------------------------------------

  nadd    = 0
  nexcetw = 0
  !   get source + target indices of reports to be duplicated
  DO iob = ntotmlo+1 , ntotml
    IF (      (momlhd(iob,nhpass) >= 0)                                        &
        .AND. (momlhd(iob,nhobtp) == npilot)                                   &
        .AND. (     (momlhd(iob,nhcode) == ntower)                             &
               .OR. (momlhd(iob,nhcode) == ntowic))) THEN
      IF (ntotml+nadd < maxmll) THEN
        nadd = nadd + 1
        ! store report indices
        itw_s (nadd) = iob
        itw_t (nadd) = ntotml + nadd
      ELSE
        nexcetw = nexcetw + 1
        kobtyp = momlhd (ntotml,nhobtp)
        kcdtyp = momlhd (ntotml,nhcode)
        icma  = i_cma ( kobtyp , kcdtyp )
        neventr (nenoml,icma) = neventr (nenoml,icma) + 1
      ENDIF
    ENDIF
  ENDDO
  !   duplicate reports
  DO itw = 1 , nadd
    omlbdy (itw_t(itw),1:maxmlv,1:mxrbdy) = omlbdy(itw_s(itw),1:maxmlv,1:mxrbdy)
    momlbd (itw_t(itw),1:maxmlv,1:mxrbdf) = momlbd(itw_s(itw),1:maxmlv,1:mxrbdf)
    omlhed (itw_t(itw)         ,1:mxrhed) = omlhed(itw_s(itw)         ,1:mxrhed)
    momlhd (itw_t(itw)         ,1:mxrhdf) = momlhd(itw_s(itw)         ,1:mxrhdf)
    yomlhd (itw_t(itw))                   = yomlhd(itw_s(itw))
  ENDDO
  DO itw = 1 , nadd
    kobtyp = momlhd (itw_t(itw),nhobtp)
    kcdtyp = momlhd (itw_t(itw),nhcode)
    icma  = i_cma ( kobtyp , kcdtyp )
    neventr (neslml,icma) = neventr (neslml,icma) + 1
  ENDDO
  !   flag levels as passive according to 'dhosag'
  DO itw = 1 , nadd
    DO ilev = 1 , momlhd(itw_s(itw),nhnlev)
      ! flag HOSAG levels as passive if above 'dhosag'
      iob = itw_s(itw)
      IF (omlbdy(iob,ilev,nbtz)-omlhed(iob,nhalt) > dhosag+0.1_wp) THEN
!       momlbd(iob,ilev,nbtflg) = IBSET( momlbd(iob,ilev,nbtflg), nvfsbp )
        momlbd(iob,ilev,nbtflg) = IBSET( momlbd(iob,ilev,nbtflg), nvffbp )
        momlbd(iob,ilev,nbtflg) = IBSET( momlbd(iob,ilev,nbtflg), nvfmbp )
      ENDIF
      ! mark height levels, and flag them as passive if below 'dhosag'
      io2 = itw_t(itw)
      momlbd (io2,ilev,nbtlsg) = IBCLR( momlbd(io2,ilev,nbtlsg), LS_STANDARD )
      momlbd (io2,ilev,nbtlid) = IBCLR( momlbd(io2,ilev,nbtlid), nvlidp(6) )
      IF (omlbdy(io2,ilev,nbtz)-omlhed(io2,nhalt) <= dhosag+0.1_wp) THEN
!       momlbd(io2,ilev,nbtflg) = IBSET( momlbd(io2,ilev,nbtflg), nvfsbp )
        momlbd(io2,ilev,nbtflg) = IBSET( momlbd(io2,ilev,nbtflg), nvffbp )
        momlbd(io2,ilev,nbtflg) = IBSET( momlbd(io2,ilev,nbtflg), nvfmbp )
      ENDIF
    ENDDO
    lrepass (itw) = .TRUE.
  ENDDO
  !   execute exceptions from 'dhosag' rule
  DO iex = 1 , ntwex
    ! go through list of exceptions
    DO itw = 1 , nadd
      ! check if the current exception refers to this station ID
      IF (      (TRIM( yomlhd(itw_s(itw)) ) == TRIM( ytwex(iex) ))) THEN
!         .AND. (htwex(iex) > -0.1) .AND. (ivtwex(iex) > 0) ! done at reading NL
        iob = itw_s(itw)
        io2 = itw_t(itw)
        DO ilev = 1 , momlhd(iob,nhnlev)
          ! check if this level is listed for exceptions
          IF ((   ABS( omlbdy(iob,ilev,nbtz)-omlhed(iob,nhalt) -htwex(iex) )   &
              <= 1.1_wp) .OR. (ABS( htwex(iex) ) < 0.1_wp)) THEN
            IF (MOD( ivtwex(iex),2 ) == 1) THEN
              ! change state of wind obs
              momlbd(iob,ilev,nbtflg) = ibflip( momlbd(iob,ilev,nbtflg), nvffbp)
              momlbd(io2,ilev,nbtflg) = ibflip( momlbd(io2,ilev,nbtflg), nvffbp)
            ENDIF
            IF (ivtwex(iex) >= 2) THEN
              ! change state of temperature / humidity obs
              momlbd(iob,ilev,nbtflg) = ibflip( momlbd(iob,ilev,nbtflg), nvfmbp)
              momlbd(io2,ilev,nbtflg) = ibflip( momlbd(io2,ilev,nbtflg), nvfmbp)
            ENDIF
            lrepass (itw) = .FALSE.
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  !   set report to passive if all obs do not meet 'dhosag' criterion
  DO itw = 1 , nadd
    IF (lrepass(itw)) THEN
      iob = itw_s(itw)
      io2 = itw_t(itw)
      ! if the whole report is to be processed at reported altitude
      ! then set the original report to passive
      IF (dhosag <= omlbdy(iob,1,nbtz) -omlhed(iob,nhalt)+0.1_wp) THEN
        momlhd (iob,nhpass) = 2
        momlhd (iob,nhflag) = IBSET( momlhd(iob,nhflag) , FL_THIN )
!       momlhd (iob,nhschr) = IBSET( momlhd(iob,nhschr) , nvexbp )
      ENDIF
      ! if the whole report is to processed at height of sensor above ground
      ! then set duplicated report to passive
      IF (dhosag >  omlbdy(io2,MAX(momlhd(io2,nhnlev),1),nbtz)                   &
                   -omlhed(io2,nhalt)+0.1_wp) THEN
        momlhd (io2,nhpass) = 2
        momlhd (io2,nhflag) = IBSET( momlhd(io2,nhflag) , FL_THIN )
!       momlhd (io2,nhschr) = IBSET( momlhd(io2,nhschr) , nvexbp )
      ENDIF
    ENDIF
  ENDDO

  ntotml  = ntotml  + nadd
  nexceed = nexceed + nexcetw

  IF (nexcetw > 0) PRINT '(''CAUTION !!!!! multi-level ODR size'',I5           &
                         &,'' too small to duplicate'',I5,'' tower reports'')' &
                         , maxmll, nexcetw

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_tower_copy
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_tower_copy


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for douplicating tower reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_bcorr_rass ( ntotmlo, ntotml, nbcrr, ilstidtw               &
                              , bcrrt, bcrrhl, bcrrhu, ybcrr )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_cdfin_org" applies bias corrections
!   to RASS (virtual) temperature according to 'RASS bias corrections rules'
!   (see input arguments 'ybcrr', 'bcrrt', 'bcrrhl', 'bcrrhu').
!   (The correction rules are determined from namelist input. Also, the wind
!    profiler reports are checked for temperature obs to be bias corrected.)
!
! Method:
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

  INTEGER        , INTENT (IN)     ::  &
    ntotmlo     ,& ! number of (multi-level) reports before having read
                   !   new data at the current timestep
    ntotml      ,& ! total number of multi-level reports
    nbcrr       ,& ! number of RASS bias correction rules
    ilstidtw       ! max. length of tower station ID's

  REAL (KIND=wp) , INTENT (IN)     ::  &
    bcrrt  (MAX(1,nbcrr)) ,& ! bias correction of RASS (virtual) temperature [K]
    bcrrhl (MAX(1,nbcrr)) ,& ! lower \ height limit [m] in RASS profile to apply
    bcrrhu (MAX(1,nbcrr))    ! upper / bias correct. 'bcrrt' for station 'ybcrr'

  CHARACTER (LEN=ilstidtw) , INTENT (IN)  :: &
    ybcrr  (MAX(1,nbcrr))    ! station ID's for RASS bias correction rules

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER        ::  &
    irrl (MAX(1,ntotml-ntotmlo)),& ! list of ODR report indices of RASS reports
    irr         ,& ! ODR  index of RASS report
    iob         ,& ! loop index over multi-level reports
    irep        ,& ! loop index over new RASS reports
    ibc         ,& ! loop index over 'RASS bias correction rules'
    ilev        ,& ! loop index over vertical levels
    nrass          ! number of new RASS reports

! Local arrays: None
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_bcorr_rass
!-------------------------------------------------------------------------------

  nrass = 0
  DO iob = ntotmlo+1 , ntotml
    IF (      (momlhd(iob,nhpass) >= 0)                                        &
        .AND. (     (momlhd(iob,nhcode) == nra_eu)                             &
               .OR. (momlhd(iob,nhcode) == nwp_eu))) THEN
      nrass = nrass + 1
      irrl (nrass) = iob
    ENDIF
  ENDDO

  !   go through list of RASS bias correction rules
  DO ibc = 1 , nbcrr
    ! loop over RASS (and wind profiler) reports
    DO irep = 1 , nrass
      irr = irrl(irep)
      ! check if the current BC rule refers to this station ID
      IF (TRIM( yomlhd(irr) ) == TRIM( ybcrr(ibc) )) THEN
        DO ilev = 1 , momlhd(irr,nhnlev)
          IF (      (omlbdy(irr,ilev,nbtz) >= bcrrhl(ibc))                     &
              .AND. (omlbdy(irr,ilev,nbtz) <  bcrrhu(ibc))                     &
              .AND. (omlbdy(irr,ilev,nbtt) >  rmdich)) THEN
            ! if obs already bias corrected, then reconstruct original value
            IF (omlbdy(irr,ilev,nbtdrh) > rmdich)                              &
              omlbdy (irr,ilev,nbtt) =   omlbdy(irr,ilev,nbtt)                 &
                                       - omlbdy(irr,ilev,nbtdrh)
            ! apply bias correction according to current rule
            omlbdy (irr,ilev,nbtdrh) = bcrrt(ibc)
            omlbdy (irr,ilev,nbtt  ) = bcrrt(ibc) + omlbdy(irr,ilev,nbtt)
            ! commented out: no not set flag
!           mflgc  =  NINT( c05 + SIGN( c05, ABS(bcrrt(ibc)) -epsy ) )
!           CALL MVBITS( mflgc, 0, 1, momlbd(irr,ilev,nbtflg), nvftbp+nvfbps(6))
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_bcorr_rass
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_bcorr_rass


!===============================================================================
!+ Module procedure in "src_obs_cdfin_org" for shear and lapse rate checks
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_mult_qualicheck ( ntotmlo, ntotml, maxmlv, ievnt, rdocp )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "src_obs_cdfin_org" performs
!   model-independent observation quality control checks for multi-level data,
!   namely a superadiabatic lapse rate check, a wind speed shear check, and a
!   directional shear check.
!
! Method:
!   The checks follow the ECMWF pre-processing (Met. Bull. M1.4/3).
!   17.02.01
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------

  INTEGER        , INTENT(IN) :: &
    ntotmlo            ,& ! number of multi-level reports before having
                          ! read new data at current timestep
    ntotml             ,& ! total number of multi-level reports
    maxmlv                ! max. number of vertical levels in multi-level report

  INTEGER        , INTENT(IN) :: &
    ievnt  (ntotml+1)     ! index for data event counters

  REAL (KIND=wp) , INTENT(IN) :: &
    rdocp                 ! r_d / cp_d

! Local parameters: None
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c90     =  90.0_wp      ,& !
    c100    = 100.0_wp      ,& !
    zlop0   =  11.5129_wp   ,& ! LOG( reference pressure for pot. temp.)
    zrpi180 =  57.2957795_wp   ! conversion from radians to degrees

! Local scalars:
! -------------

  INTEGER        ::  &
    iob              ,& ! index of report in the ODR (obs data record)
    kobtyp , kcdtyp  ,& ! observation type , code type
    nlev             ,& ! number of obs. levels in report 'iob'
    ilev             ,& ! vertical loop index over observation levels
    klev             ,& ! loop index over levels to be checked
    nclev            ,& ! number of obs. levels containing checked quantity
    ilvqdd           ,& ! index determining threshold for directional shear
    ilvq             ,& ! loop index over standard pressure levels
    ilva   , ilvb    ,& ! indices of obs. levels rejected by wind shear check
    nphpa            ,& ! pressure [hPa] of current obs. level
    nfdbp            ,& ! bit position of rejection flag in main flag
    ilen             ,& ! length of output record
    icl                 ! loop index

  REAL (KIND=wp) ::  &
!   zdth             ,& ! potential temperature difference
    zthcorr          ,& ! -0.5*(pot. temperature decrease with height minus
                        !       allowed pot. temperature decrease with height)
                        ! = +/- pot. temperature correction to upper/lower layer
    shearff          ,& ! wind speed shear
    sheardd          ,& ! directional shear
    sumff               ! sum of wind speeds of 2 adjacent mandatory levels

  LOGICAL        ::  &
    lactrep             ! .TRUE. if the report is active

! Local (automatic) arrays:
! -------------------------

  INTEGER        ::  &
    iclev   (maxmlv) ,& ! level indices of levels containing checked quantity
    krej    (maxmlv) ,& ! indicator for levels (data) to be rejected
    kqcdd   (15)        ! type of limit values used for directional shear check

  REAL (KIND=wp) ::  &
    zthob   (maxmlv) ,& ! potential temperature
    zsupcor (maxmlv) ,& ! superadiabatic correction limit
    zdd     (15)     ,& ! wind direction                   \   at mandatory
    zff     (15)     ,& ! wind speed                        >  pressure
    zshear  (15)        ! shear of direction or speed      /   levels
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_mult_qualicheck
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Superadiabatic lapse rate check
!-------------------------------------------------------------------------------

! ___________________________
  DO iob = ntotmlo+1 , ntotml
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    kobtyp  =  momlhd(iob,nhobtp)
    kcdtyp  =  momlhd(iob,nhcode)
    nlev    =  momlhd(iob,nhnlev)
    lactrep = (momlhd(iob,nhpass) == 0)

! compute potential temperature and specify superadiabatic correction
! -------------------------------------------------------------------

    nclev = 0
    DO ilev = 1 , nlev
      IF (BTEST( momlbd(iob,ilev,nbterr), nvrt )) THEN
        nclev          =  nclev + 1
        iclev (nclev)  =  ilev
        krej  (nclev)  =  0
        !   use 'nbtplv' instead of 'nbtp' to ensure check is independent from
        !   model state (ensemble member) even though the potential temperature
        !   and hence the lapse rate is not quite correct (relevant for RASS)
        zthob (nclev)  =  omlbdy(iob,ilev,nbtt)                                &
                        * EXP( rdocp *(zlop0 - LOG( omlbdy(iob,ilev,nbtplv) )) )
!                       * EXP( rdocp *(zlop0 - omlbdy(iob,ilev,nbtlop)) )
        IF (omlbdy(iob,ilev,nbtplv) <  40000.0_wp) THEN
          zsupcor (nclev) = - 0.5_wp
        ELSEIF (omlbdy(iob,ilev,nbtplv) > 100000.0_wp) THEN
          zsupcor (nclev) = - 4.5_wp
        ELSEIF (omlbdy(iob,ilev,nbtplv) >  85000.0_wp) THEN
          zsupcor (nclev) = - 3.5_wp
        ELSEIF (omlbdy(iob,ilev,nbtplv) >  70000.0_wp) THEN
          zsupcor (nclev) = - 2.5_wp
        ELSEIF (omlbdy(iob,ilev,nbtplv) >  50000.0_wp) THEN
          zsupcor (nclev) = - 1.5_wp
        ELSE
          zsupcor (nclev) = - 1.0_wp
        ENDIF
      ENDIF
    ENDDO

!   IF ((losse_fg) .AND. (fperturb > epsy)) THEN
    IF (losse_fg) THEN
      DO klev = 2 , nclev - 2
        IF (zthob(klev+1)-zthob(klev) < zsupcor(klev)) THEN
          zthcorr =  0.5_wp* (zsupcor(klev) - zthob(klev+1) + zthob(klev))
          ilev    =  iclev(klev)
          omlbdy(iob,ilev,nbtt) = omlbdy(iob,ilev,nbtt)                        &
              - zthcorr *EXP( -rdocp *(zlop0 - LOG( omlbdy(iob,ilev,nbtplv) )) )
          ilev    =  iclev(klev+1)
          omlbdy(iob,ilev,nbtt) = omlbdy(iob,ilev,nbtt)                        &
              + zthcorr *EXP( -rdocp *(zlop0 - LOG( omlbdy(iob,ilev,nbtplv) )) )
        ENDIF
      ENDDO
      nclev = 0
    ENDIF

! determine levels to be rejected
! -------------------------------

    DO klev = 2 , nclev - 2
      IF (zthob(klev+1)-zthob(klev) < zsupcor(klev)) THEN
        IF (      (zthob(klev+1)-zthob(klev-1) <  zsupcor(klev-1))             &
            .AND. (zthob(klev+2)-zthob(klev  ) >= zsupcor(klev  ))) THEN
          krej (klev+1) = -1
        ELSEIF (      (zthob(klev+1)-zthob(klev-1) >= zsupcor(klev-1))         &
                .AND. (zthob(klev+2)-zthob(klev  ) <  zsupcor(klev  ))) THEN
          krej (klev  ) =  1
        ELSE
          krej (klev  ) =  1
          krej (klev+1) =  2
        ENDIF
      ENDIF
    ENDDO
    IF (nclev >= 2) THEN
      IF (zthob(2)-zthob(1) < zsupcor(1)) THEN
        krej (1) =  1
        krej (2) =  2
      ELSEIF (zthob(nclev)-zthob(nclev-1) < zsupcor(nclev-1)) THEN
        krej (nclev-1) =  1
        krej (nclev  ) =  2
      ENDIF
    ENDIF

! flagging of data which are rejected by the lapse-rate check
! -----------------------------------------------------------
! (rejection of temperature implies rejection of humidity)

    DO klev = 1 , nclev
      IF (krej(klev) /= 0) THEN
        ilev  =  iclev(klev)
        momlbd (iob,ilev,nbterr) = IBCLR ( momlbd(iob,ilev,nbterr), nvrt )
        momlbd (iob,ilev,nbterr) = IBCLR ( momlbd(iob,ilev,nbterr), nvrq )
        IF (.NOT. BTEST( momlbd(iob,ilev,nbtflg), nvftbp+nvfbps(4) ))          &
          momlbd (iob,ilev,nbtflg) = IBSET( momlbd(iob,ilev,nbtflg)            &
                                          , nvftbp +nvfbps(4) )
        IF (lactrep) THEN
          neventd (netlps,ievnt(iob)) = neventd (netlps,ievnt(iob)) + 1
          ilen = 7 + istrej
          IF ((nacout+ilen <= nmxoln) .AND. (krej(klev) /= 2)) THEN
            outbuf(nacout+1) = ilen
            outbuf(nacout+2) = nfmt7
            outbuf(nacout+3) = NINT( zsupcor(klev) * c100 )
            outbuf(nacout+4) = NINT( (zthob(klev+1) - zthob(klev)) *c100 )
            outbuf(nacout+5) = NINT( omlbdy(iob,ilev,nbtp) )
            outbuf(nacout+6) = NINT( omlbdy(iob,ilev,nbtp) )
            outbuf(nacout+7) = NINT( omlhed(iob,nhtime) * c100 )
            IF (krej(klev+1) == 2) THEN
              outbuf(nacout+6) = NINT( omlbdy(iob,iclev(klev+1),nbtp) )
            ELSEIF (krej(klev) == -1) THEN
              outbuf(nacout+3) = NINT( zsupcor(klev-1) * c100 )
              outbuf(nacout+4) = NINT( (zthob(klev  ) - zthob(klev-1)) *c100 )
            ENDIF
            DO icl = 1 , istrej
              outbuf(nacout+7+icl) = ICHAR( yomlhd(iob) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
        ENDIF
      ENDIF
    ENDDO

!-------------------------------------------------------------------------------
!  Section 3: Wind shear checks  (for TEMP , PILOT only)
!-------------------------------------------------------------------------------

    IF (.NOT. ((kobtyp == ntemp) .OR. (kobtyp == npilot)))  nlev = -1

! compute wind direction and speed of mandatory levels (in rotated coordinates)
! -----------------------------------------------------------------------------
! (if 'zff < qcfddff', then direction is not required since the speed check will
!  not be passed whenever there is a chance not to pass the direction check)

    nclev = 0
    DO ilev = 1 , nlev
      IF (BTEST( momlbd(iob,ilev,nbterr), nvru )) THEN
        nphpa = NINT( omlbdy(iob,ilev,nbtplv) * 0.01_wp )
        IF (     (nphpa == 1000) .OR. (nphpa == 850) .OR. (nphpa == 700)       &
            .OR. ((MOD( nphpa,100 ) == 0) .AND. (nphpa <= 500))                &
            .OR. ((MOD( nphpa, 50 ) == 0) .AND. (nphpa <= 250))                &
            .OR. ((MOD( nphpa, 10 ) == 0) .AND. (nphpa <=  30))) THEN
          nclev          =  nclev + 1
          iclev (nclev)  =  ilev
          krej  (nclev)  =  0
          zdd   (nclev)  =  c0
          zshear(nclev)  =  c0
          kqcdd (nclev)  =  2
          IF ((nphpa >  849) .OR. (nphpa < 151))  kqcdd (nclev) = 1
          IF  (nphpa == 700)                      kqcdd (nclev) = 0
        ENDIF
      ENDIF
    ENDDO
    DO klev = 1 , nclev
      ilev = iclev(klev)
      IF (ABS( omlbdy(iob,ilev,nbtu) ) > epsy) THEN
        zff (klev) = SQRT(  omlbdy(iob,ilev,nbtu) * omlbdy(iob,ilev,nbtu)      &
                          + omlbdy(iob,ilev,nbtv) * omlbdy(iob,ilev,nbtv) )
        IF (zff(klev) > qcfddff(MAX( kqcdd(klev),1 )))                         &
          zdd (klev) = 180._wp + SIGN( c90 , omlbdy(iob,ilev,nbtu) )           &
                                   - ATAN(  omlbdy(iob,ilev,nbtv)              &
                                          / omlbdy(iob,ilev,nbtu) ) *zrpi180
      ELSE
        zff (klev) = ABS( omlbdy(iob,ilev,nbtv) )
        IF (zff(klev) > qcfddff(MAX( kqcdd(klev),1 )))                         &
          zdd (klev) = 270._wp - SIGN( c90 , omlbdy(iob,ilev,nbtv) )
      ENDIF
      IF (kqcdd(klev) == 0)  kqcdd (klev) = 2
    ENDDO

! determine levels to be rejected
! -------------------------------

    DO klev = 1 , nclev - 1
! wind speed shear check
      shearff = ABS( zff(klev) - zff(klev+1) )
      sumff   =      zff(klev) + zff(klev+1)
      IF (shearff > 20.6_wp+0.275_wp*sumff) THEN
        krej   (klev)  =  1
        zshear (klev)  =  shearff
! compute wind direction shear (if needed),
! and indicate if check needs to be done
      ELSEIF (      (NINT( sumff ) >= nqcddff(nnqcdd,kqcdd(klev)))             &
              .AND. (MIN( zdd(klev),zdd(klev+1) ) > epsy)) THEN
        sheardd = ABS( zdd(klev) - zdd(klev+1) )
        IF (sheardd > 180._wp)  sheardd = 360._wp - sheardd
        IF (sheardd > REAL ( nqcdd(1), wp )+epsy) THEN
          krej   (klev)  = -1
          zshear (klev)  =  sheardd
        ENDIF
      ENDIF
    ENDDO

! wind direction shear check
    DO klev = 1 , nclev - 1
      IF (krej(klev) == -1) THEN
        krej (klev) = 0
        sumff   = zff(klev) + zff(klev+1)
        ilvqdd  = nnqcdd
        DO ilvq = nnqcdd - 1 , 1 , -1
          IF (NINT( sumff ) >= nqcddff(ilvq,kqcdd(klev)))  ilvqdd = ilvq
        ENDDO
        IF (zshear(klev) > REAL (nqcdd(ilvqdd), wp)+epsy)   krej (klev) = 1 + ilvqdd
!       PRINT '(''sheardd: '',A ,F5.1, 2I4, F8.0, F7.1,I2, F6.1,I2)'           &
!              , yomlhd(iob), omlhed(iob,nhtime), klev, iclev(klev)            &
!              , omlbdy(iob,iclev(klev),nbtp), sumff, ilvqdd                   &
!              , zshear(klev), krej(klev)
      ENDIF
    ENDDO

! flagging of data not to be accepted by the wind shear check
! -----------------------------------------------------------

    DO klev = 1 , nclev
      IF (krej(klev) /= 0) THEN
        ilvb  =  iclev(klev)
        ilva  =  iclev(klev+1)
        DO ilev = ilvb , ilva
          IF (BTEST( momlbd(iob,ilev,nbterr), nvru )) THEN
            momlbd (iob,ilev,nbterr) = IBCLR ( momlbd(iob,ilev,nbterr), nvru )
            IF (krej(klev) == 1)  nfdbp = nvfubp + nvfbps(4)
            IF (krej(klev) >= 2)  nfdbp = nvfubp + nvfbps(5)
            IF (.NOT. BTEST( momlbd(iob,ilev,nbtflg), nfdbp ))                 &
              momlbd (iob,ilev,nbtflg) = IBSET( momlbd(iob,ilev,nbtflg), nfdbp )
          ENDIF
        ENDDO
        IF ((lactrep) .AND. (krej(klev) == 1))                                 &
          neventd (nefshr,ievnt(iob)) = neventd (nefshr,ievnt(iob)) + 1
        IF ((lactrep) .AND. (krej(klev) >= 2))                                 &
          neventd (nedshr,ievnt(iob)) = neventd (nedshr,ievnt(iob)) + 1
        IF (lactrep) THEN
          ilen = 7 + istrej
          IF (nacout+ilen <= nmxoln) THEN
            outbuf(nacout+1) = ilen
            IF (krej(klev) == 1) THEN
              outbuf(nacout+2) = nfmt8
              outbuf(nacout+3) = NINT( 2060.0_wp                               &
                                      +  27.5_wp *(zff(klev) + zff(klev+1)))
              outbuf(nacout+4) = NINT( zshear(klev) * c100 )
            ELSE
              outbuf(nacout+2) = nfmt9
              outbuf(nacout+3) = NINT( nqcdd(krej(klev)-1) * c100 )
              outbuf(nacout+4) = NINT( zshear(klev) * c100 )
            ENDIF
            outbuf(nacout+5) = NINT( omlbdy(iob,ilvb,nbtp) )
            outbuf(nacout+6) = NINT( omlbdy(iob,ilva,nbtp) )
            outbuf(nacout+7) = NINT( omlhed(iob,nhtime) * c100 )
            DO icl = 1 , istrej
              outbuf(nacout+7+icl) = ICHAR( yomlhd(iob) (icl:icl) )
            ENDDO
            nacout  = nacout + ilen
          ENDIF
        ENDIF
      ENDIF
    ENDDO
! _____
  ENDDO
! ~~~~~

!-------------------------------------------------------------------------------
! End of module procedure obs_cdf_mult_qualicheck
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_mult_qualicheck


!===============================================================================
!+ Module procedure in "src_obs_proc_cdf" for bias correction of RS92 humidity
!-------------------------------------------------------------------------------

SUBROUTINE obs_cdf_raso_rh_bias ( mqcorr92_in, ntotmlo, ntotml )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_cdf" computes a bias
!   correction for Vaisala RS92 humidity observations.
!
! Method:
!   Correction computed for clear-sky conditions approximately according to:
!   Miloshevich et al., 2009, J. Geophys. Res. (Atm.), D11305 :
!   The mean percentage bias is a function pressure and given for different
!   relative humidity values; here, dry conditions are extrapolated.
!   Only levels below 100 hPa are corrected, further above, different
!   equations not implemented here do apply.
!   The solar elevation angle needed for the solar radiation bias error is
!   computed using the same equations as in the radiation parameterisation.
!   An additional empirical reduction of the solar radiation bias error is
!   added to account for cloudy layers.
!
! Initial release: Christoph Schraff, 25.08.09
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
    mqcorr92_in      ,& ! = 0 : no correction for humidity
                        ! = 1 : correct only solar radiation bias for humidity
                        ! = 2 : correct total humidity bias (incl. nighttime bias)
    ntotmlo          ,& ! number of multi-level reports read before current
                        !   timestep
    ntotml              ! total number of  multi-level reports

! Local parameters:
! ----------------

  REAL (KIND=wp) , PARAMETER  :: &
    c_n_rh (4) = (/ 0.50_wp, 0.42_wp, 0.25_wp, 0.12_wp /)                     ,&
    ! coefficients for RS92 mean percentage bias as a function of pressure and
    ! for a given relative humidity at night
    c_n (4,7) = RESHAPE( (/                                                    &
      ! RH >= 50            = 42            = 25 (30,20)          = 12  %      &
   3.7854E+01_wp,-7.5226E+00_wp,-8.4463E+00_wp,-8.6609E+00_wp,                 &
  -4.9026E-01_wp,-9.4387E-02_wp,-6.7739E-02_wp,-2.3153E-01_wp,                 &
   2.0313E-03_wp, 5.6012E-04_wp, 2.1850E-04_wp, 1.1601E-03_wp,                 &
  -3.9299E-06_wp,-1.0285E-06_wp, 2.4128E-07_wp,-1.6559E-06_wp,                 &
   3.9439E-09_wp, 8.1621E-10_wp,-1.1680E-09_wp, 4.7114E-10_wp,                 &
  -1.9776E-12_wp,-2.4513E-13_wp, 1.1593E-12_wp, 6.4842E-13_wp,                 &
   3.8808E-16_wp, 3.3189E-18_wp,-3.6948E-16_wp,-3.7600E-16_wp/)                &
                       , (/4,7/) )

  REAL (KIND=wp) , PARAMETER  :: &
    c_d_rh (4) = (/ 0.34_wp, 0.22_wp, 0.11_wp, 0.05_wp /)                     ,&
    ! coefficients for RS92 mean percentage bias as a function of pressure and
    ! for a given relative humidity at daytime (solar elevation angle = 66 deg.)
    c_d (4,6) = RESHAPE( (/                                                    &
      ! RH >= 34            = 22               = 11               =  5  %      &
  -6.0024E+01_wp,-6.6938E+01_wp,-6.7112E+01_wp,-6.6681E+01_wp,                 &
   1.4726E-01_wp, 1.1812E-01_wp, 1.1009E-01_wp, 1.4741E-01_wp,                 &
  -6.9462E-05_wp, 2.8349E-04_wp, 3.7366E-04_wp, 1.6426E-05_wp,                 &
  -2.0216E-07_wp,-1.0166E-06_wp,-1.2284E-06_wp,-1.4146E-07_wp,                 &
   3.1579E-10_wp, 1.0377E-09_wp, 1.2520E-09_wp, 8.9222E-12_wp,                 &
  -1.3450E-13_wp,-3.5797E-13_wp,-4.3857E-13_wp, 4.0390E-14_wp/)                &
                       , (/4,6/) )

  REAL (KIND=wp) , PARAMETER  :: &
    ! coefficients for solar radiation error
    f_sre  (4) = (/ -1.6061E-3_wp, 3.7746E-2_wp                                &
                  , -4.7402E-4_wp, 2.0018E-6_wp /)

  REAL (KIND=wp) , PARAMETER  :: &
    ! parameters for mixed water/ice phase (set as in ECMWF IFS (CY31r1))
    tmpmin = 250.16_wp ,& ! minimum temperature of mixed-phase T-range [K]
    tmpmax = 273.16_wp ,& ! maximum temperature of mixed-phase T-range [K]
    exp_mp = 2._wp        ! exponent in the interpolation formula
                              ! for the mixed-phase water fraction [-]

  REAL (KIND=wp) , PARAMETER  ::  &
    zlwpcal = 60.0_wp  ,& ! e-folding decay scale [g/m2] to reduce (f_cloud)
                              !   the solar radiation error due to cloud layers
    c100r   = 0.01_wp     ! maximum number of elements to be received

  INTEGER        , PARAMETER  :: &
    mdd_offs (13) = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334   &
                                                                       , 365/)

  REAL (KIND=wp) ::  &
    pi                  ! circle constant

! Local scalars:
! -------------

  INTEGER        ::  &
    mqcorr92         ,& ! mqcorr92_in
    nob              ,& ! loop index over reports
    klev             ,& ! loop index over report levels
!   iyyyy            ,& ! year (incl. century) of observation
!   imm              ,& ! month                of observation
!   idd              ,& ! day                  of observation
!   iddjul           ,& ! Julian day (i.e. days since the beginning of the year)
!   nlev             ,& ! number of vertical levels in current report
    ilva   , ilvb    ,& ! levels of c_n or c_d next to the obs level
    mflgc               ! flag indicating non-zero bias correction

  REAL (KIND=wp) ::  &
    zrhorig          ,& ! originally reported relative humidity
    zrhcorr          ,& ! bias corrected      relative humidity
!   zhour            ,& ! hour of observation
!   ztwo   , ztho    ,& ! factors in time equation (day (angle) of solar year)
!   zdtzgl , zdek    ,& ! factors in time equation (day (angle) of solar year)
!   zlonoon          ,& ! longitude where it is exactly noon at the obs time
!   ztimrad          ,& ! longitudinal radians related to the time
!   zdeksin, zdekcos ,& ! factors from time equation
!   zsinlat, zcoslat ,& ! factors depending on latitude
!   zcosthi, zcoszen ,& ! COS( solar zenith angle )  / above horizon
!   zsolzen          ,& ! solar zenith angle    [rad]
    zsolelv          ,& ! solar elevation angle [deg]
    zsolerr          ,& ! fraction of the solar radiation error SRE
    r_fsol           ,& ! solar angle dependent factor for cloud reflectivity
    f_cloud          ,& ! reduction of SRE due to cloud layers
    fbias            ,& ! total percentage bias
    fbias_s          ,& ! mean percentage solar radiation bias error
    fbias_n, fbias_d ,& ! nighttime / daytime mean percentage bias at obs level
    fb_n2  , fb_d2   ,& ! nighttime / daytime mean percentage bias at c_n level
    zfrh             ,& ! interpolation factor
    zqc    , zqc1    ,& ! estimated specific water content [g/kg] /previous lev.
    zlwp             ,& ! vertical liquid water path
    zpp              ,& ! pressure [hPa] of observation level
    zppcl            ,& ! pressure [hPa] of next layer above if cloudy
    zpp2   , zpp3    ,& ! square( zpp )  ;  zpp **3
    ztt              ,& ! temperature [K] at observation level
!   ztbias           ,& ! temperature bias [K]
    fr_wat           ,& ! water fraction in mixed phase
    es_w   , es_m    ,& ! saturation vapour pressure over water /for mixed phase
    es_o                ! bias-corrected observed vapour pressure

! CHARACTER (LEN=25)       :: &
!   yroutine            ! name of this subroutine
! CHARACTER (LEN=ilstidp)     :: &
!   ystid               ! observation type (input array for local blacklist)

! Local arrays:
! ------------
!
!============ End of header ====================================================

!-------------------------------------------------------------------------------
! Begin Subroutine obs_cdf_raso_rh_bias
!-------------------------------------------------------------------------------

! yroutine = 'obs_cdf_raso_rh_bias'

  mqcorr92  =  mqcorr92_in
  IF ((mqcorr92 >= 3) .OR. (mqcorr92 < 0))  mqcorr92 = 0

! IF (ntotml > ntotmlo)  pi = 4.0_wp * ATAN( c1 )

!-------------------------------------------------------------------------------
! Section 1: solar zenith angle at observation time + location (multi-level obs)
!-------------------------------------------------------------------------------

  DO nob = ntotmlo + 1 , ntotml

    zsolelv = c0
    IF (mqcorr92 >= 1) THEN
!     imm   = MOD( momlhd(nob,nhdate) , 10000 ) / 100
!     idd   = MOD( momlhd(nob,nhdate) , 100 )
!     iyyyy =      momlhd(nob,nhdate)           / 10000
!     zhour =   REAL(      momlhd(nob,nhhrmn) / 100   )                        &
!             + REAL( MOD( momlhd(nob,nhhrmn) , 100 ) ) / 60._wp
!     ! get Julian day (i.e. days since the beginning of the years)
!     iddjul  =  mdd_offs(imm) + idd
!     IF (imm > 2)  iddjul  =  iddjul + ileap(iyyyy)
!     ! use time equation as in src_radiation
!     ztwo    = 0.681 + 0.2422*(iyyyy-1949)-(iyyyy-1949)/4
!     ztho    = 2.*pi*( REAL(iddjul, wp) -1.0 + ztwo )/365.2422
!     zdtzgl  = 0.000075 + 0.001868*COS(   ztho) - 0.032077*SIN(   ztho)       &
!                        - 0.014615*COS(2.*ztho) - 0.040849*SIN(2.*ztho)
!     zdek    = 0.006918 - 0.399912*COS(   ztho) + 0.070257*SIN(   ztho)       &
!                        - 0.006758*COS(2.*ztho) + 0.000907*SIN(2.*ztho)       &
!                        - 0.002697*COS(3.*ztho) + 0.001480*SIN(3.*ztho)
!     ! longitude-dependent part
!       ! get the longitude where it is exactly noon at the obs time
!     zlonoon = pi*(zhour-12._wp)/12._wp + zdtzgl
!     ztimrad = zlonoon + r_degrad* omlhed(nob,nhilon)

!!    WRITE(0,'("SOLAR ",A ,6F7.2,F7.4,2I11,2I4)' )                            &
!!          yomlhd(nob), zhour, omlhed(nob,nhjlat), omlhed(nob,nhilon)         &
!!        , zdtzgl, zdek, ztimrad, r_degrad, momlhd(nob,nhdate)                &
!!        , momlhd(nob,nhhrmn), ntotmlo, ntotml
!!    WRITE(0,'("SOLAR1 ",I3,1X,A ,5I6,F7.2,2I11,2F8.4)' )  my_cart_id,        &
!!          yomlhd(nob), ntotmlo, ntotml, nob, iddjul, iyyyy, zhour            &
!!                     , momlhd(nob,nhdate), momlhd(nob,nhhrmn)                &
!!                     , omlhed(nob,nhjlat), omlhed(nob,nhilon)
!!    WRITE(0,'("SOLAR2 ",I3,1X,A ,4(F9.2,1X)         )' )  my_cart_id,        &
!!          yomlhd(nob), zdtzgl, zdek, ztimrad, r_degrad

!     ! latitude-dependent part
!     zdeksin = SIN( zdek )
!     zdekcos = COS( zdek )
!     zsinlat = SIN( r_degrad* omlhed(nob,nhjlat) )
!     zcoslat = SQRT( c1 - zsinlat*zsinlat )
!     ! get COS( solar zenith angle )  (above horizon)
!     zcosthi = zdeksin * zsinlat + zdekcos * zcoslat * COS( ztimrad )
!     zcoszen = MIN( MAX( zcosthi, c0 ) , c1 )
!     ! get solar zenith angle
!     zsolzen = ACOS( zcoszen )
!     zsolelv = 90.0_wp - zsolzen /r_degrad
      zsolelv = 90.0_wp - omlhed(nob,nhsolz)
    ENDIF
!   IF ((zsolelv < -epsy) .OR. (zsolelv > 90.0_wp+epsy))                       &
!     WRITE(0,'("SOLAR ",A ,I5,4F7.2,F8.3)' )                                  &
!            yomlhd(nob), nob, zhour, omlhed(nob,nhjlat), omlhed(nob,nhilon)   &
!                       , zsolelv, zsolzen

!-------------------------------------------------------------------------------
! Section 2: Bias correction for RS92 humidity
!-------------------------------------------------------------------------------

    ! if Vaisala RS92 report and bias correction to be done
    IF ((mqcorr92 >= 1) .AND. (momlhd(nob,nhobtp) == ntemp)                    &
                        .AND. (     (momlhd(nob,nhrtyp) == RS_92_DIG12)        &
                               .OR. (momlhd(nob,nhrtyp) == RS_92_DIG3 )        &
                               .OR. (momlhd(nob,nhrtyp) == RS_92_AUTO ))) THEN

! fraction of the solar radiation error SRE
! (ratio of SRE at current solar elevation angle to SRE at 66 deg.)
! (Fig.10 in Miloshevich et al., 2009, J. Geophys. Res. (Atm.), D11305)
! ---------------------------------------------------------------------
      IF (zsolelv <= c0+epsy) THEN
        zsolerr = c0
      ELSE
        zsolerr = f_sre(1) + f_sre(2) *zsolelv + f_sre(3)*zsolelv*zsolelv      &
                           + f_sre(4) *zsolelv*zsolelv*zsolelv
      ENDIF

      ! get solar angle dependent correction for cloud reflectivity
!     r_fsol = COS( c05* zsolzen )
      r_fsol = COS( c05* r_degrad* omlhed(nob,nhsolz) )
      r_fsol = r_fsol * r_fsol

!     WRITE(0,'("SOLAR ",A ,6F7.2,F9.2,F6.3)' )                                &
!           yomlhd(nob), zhour, omlhed(nob,nhjlat), omlhed(nob,nhilon),        &
!           zdtzgl, zdek, zcosthi, zsolelv, zsolerr

! initialise bias correction
! --------------------------
      f_cloud =  c1
      zppcl   = -c1
      zqc1    =  c0

      DO klev = momlhd(nob,nhnlev) , 1 , -1

        ! uncorrected reported relative humidity value
        zrhorig  = omlbdy(nob,klev,nbtrh )
        IF (      (omlbdy(nob,klev,nbtrh ) > rmdich)                           &
            .AND. (omlbdy(nob,klev,nbtdrh) > rmdich))                          &
          zrhorig = omlbdy(nob,klev,nbtrh ) - omlbdy(nob,klev,nbtdrh)
        ! preset zrhcorr to a defined value
        zrhcorr = zrhorig

! compute bias correction (if RH exists and < 100 % over water and p >= 100 hPa)
! -----------------------
! set RH > rhtsat = 96 % to 100 % (because values > 99 % are hardly reported,
!                          and reported values > 96 % often indicate saturation)
        IF (zrhorig > rhtsat(3)) THEN
          zrhcorr = c1
! Miloshevich correction if RH <= rhtsat = 96 % and p >= 100 hPa
!       ELSEIF (      (zrhorig > rmdich) .AND. (zrhorig <= c1-epsy)            &
        ELSEIF (      (zrhorig > rmdich)                                       &
                .AND. (omlbdy(nob,klev,nbtp) >= 10000._wp-epsy)) THEN

!       ! solar radiation bias correction using factor from Voemel 2007
!         IF (zsolerr > epsy) THEN
!         ! divisional factor according to Voemel et al., 2007
!           ! pressure dependent part (temperature dependent part is neglected)
!           zlnhpa  = LOG( omlbdy(nob,klev,nbtp) * c100r )
!           zfacvoe = - 0.12158_wp* zlnhpa* zlnhpa                             &
!                     + 1.66400_wp* zlnhpa         - 4.7855_wp
!         ! scale deviation of Voemel factor from 1 by solar factor
!           zfacvoe = c1 + (zfacvoe - c1)* zsolerr
!         ! apply divisional Voemal factor (zfacvoe <= 1)
!           zrhcorr = zrhorig / zfacvoe
!         ENDIF

        ! bias correction approx. according to Miloshevich et al., 2009
          zpp  = omlbdy(nob,klev,nbtp ) * c100r
          zpp2 = zpp * zpp
          zpp3 = zpp2* zpp
          ! nighttime mean percentage bias
          ilvb = 4
          IF (zrhorig >= c_n_rh(3))  ilvb = 3
          IF (zrhorig >= c_n_rh(2))  ilvb = 2
          IF (zrhorig >= c_n_rh(1))  ilvb = 1
          fbias_n = c_n(ilvb,1) + c_n(ilvb,2)*zpp       + c_n(ilvb,3)*zpp2     &
                                + c_n(ilvb,4)*zpp3      + c_n(ilvb,5)*zpp3*zpp &
                                + c_n(ilvb,6)*zpp3*zpp2 + c_n(ilvb,7)*zpp3*zpp3
          IF (ilvb >= 2) THEN
            ! rather dry conditions
            ilva  = ilvb - 1
            zfrh  = (c_n_rh(ilva) - zrhorig) / (c_n_rh(ilva) - c_n_rh(ilvb))
            fb_n2 = c_n(ilva,1) + c_n(ilva,2)*zpp       + c_n(ilva,3)*zpp2     &
                                + c_n(ilva,4)*zpp3      + c_n(ilva,5)*zpp3*zpp &
                                + c_n(ilva,6)*zpp3*zpp2 + c_n(ilva,7)*zpp3*zpp3
            fbias_n = zfrh *fbias_n  +  (c1-zfrh) *fb_n2
          ENDIF
          ! solar radiation percentage bias for clear-sky conditions
          IF (zsolerr <= epsy) THEN
            fbias_s = c0
          ELSE
            ! daytime mean percentage bias (valid for solar elevation = 66 deg.)
            ilvb = 4
            IF (zrhorig >= c_d_rh(3))  ilvb = 3
            IF (zrhorig >= c_d_rh(2))  ilvb = 2
            IF (zrhorig >= c_d_rh(1))  ilvb = 1
            fbias_d = c_d(ilvb,1) + c_d(ilvb,2)*zpp     + c_d(ilvb,3)*zpp2     &
                                  + c_d(ilvb,4)*zpp3    + c_d(ilvb,5)*zpp3*zpp &
                                  + c_d(ilvb,6)*zpp3*zpp2
            zfrh  = c1
            fb_d2 = fbias_d
            IF (ilvb >= 2) THEN
              ! rather dry conditions
              ilva  = ilvb - 1
              zfrh  = (c_d_rh(ilva) - zrhorig) / (c_d_rh(ilva) - c_d_rh(ilvb))
              fb_d2 = c_d(ilva,1) + c_d(ilva,2)*zpp     + c_d(ilva,3)*zpp2     &
                                  + c_d(ilva,4)*zpp3    + c_d(ilva,5)*zpp3*zpp &
                                  + c_d(ilva,6)*zpp3*zpp2
              fbias_d = zfrh *fbias_d  +  (c1-zfrh) *fb_d2
            ENDIF
            fbias_s = (fbias_d - fbias_n) * zsolerr
          ! account for cloudy layers further above
            fbias_s = fbias_s * f_cloud

!           WRITE(0,'("FB_ND ",A ,I4,F6.0,F8.4,I3,F6.2,4F8.3,F9.0,F12.0)' )    &
!                    yomlhd(nob), klev, zpp, zrhorig, ilvb, zfrh               &
!                  , fb_d2, fbias_d, fbias_n, fbias_s, zpp2, zpp3
          ENDIF
          ! total percentage bias
          IF (mqcorr92 == 1) fbias =           fbias_s !corr only solar rad bias
          IF (mqcorr92 == 2) fbias = fbias_n + fbias_s !correct total bias
          ! bias-corrected relative humidity
          zrhcorr = MAX( epsy , MIN( c1 , zrhorig /(c1 + c100r* fbias) ) )

!           WRITE(0,'("FBIAS ",A ,I4,F6.0,F8.4,3F9.3)' )                       &
!                 yomlhd(nob), klev, zpp, zrhcorr, fbias_s, fbias, f_cloud
        ENDIF

! compute cloud correction 'f_cloud' to solar radiation error (if RH,T,p exist)
! -----------------------------------------------------------
! - liquid water path in the vertical direction LWP:
!     LWP [in g/m2] = qc_mean [g/kg] * density [kg/m3] * delta_z [m]
!                   = qc_mean [g/kg] * delta_p [Pa=kg/(m*s2)] / g [m/s2]
!                   = 10 * qc_mean [in g/kg] * delta_p [in hPa]
!   - assuming qc_mean = 0.1 g/kg for water clouds, it follows:
!       LWP [in g/m2] = cloud thickness delta_p [in hPa]
!   - for ice clouds, qc_mean = 0.1 g/kg for convective or lower troposheric
!                 and qc_mean = 0.01 g/kg for cirrus cloud are reasonable
!     --> let qc_mean vary linearly in p from 0.1 at 700 hPa to 0.01 at 300 hPa
!   ==> qc_mean = fr_wat*0.1 + (1-fr_wat)*(.01 +.09* MIN(1,MAX(0,(p-300)/400 )))
! - required: diffuse transmissivity T:
!             T  =  1 - cloud reflectivity R  -  absorptivity
!   from Slingo 1989 (J.Atm.Sci.), Fig.1: absorptivity can be nearly neglected,
!     and approx.: T(zsolzen=0) = 1 - R(zsolzen=0) = exp( - LWP / 60 )
!                  T = 1 - R    = 1 - ( R(0) + (1-R(0))*(sin(zsolzen/2))^2 )
!                               = (1-R(0)) *(cos(zsolzen/2))^2
!   ==> T = exp( - LWP / 60 ) * (cos(zsolzen/2))^2
!   since several obs within the profile can indicate cloud, the total LWP is
!   split up into several pieces:
!       T = exp(-(LWP1+LWP2+...+LWPn)/60) * (cos(zsolzen/2))^2
!         = exp(-LWP1/60) * (cos(zsolzen/2))^2 * exp(-LWP2/60)*...*exp(-LWPn/60)
        IF (      (zrhcorr                 > rmdich)                           &
            .AND. (omlbdy(nob,klev,nbtt) > rmdich)                             &
            .AND. (omlbdy(nob,klev,nbtp) > rmdich)) THEN
          ztt  = omlbdy(nob,klev,nbtt )
          zpp  = omlbdy(nob,klev,nbtp ) * c100r
          ! check if current layer is cloudy (obs vapour pressure > saturation)
            ! use saturation for mixed phase region as in Tiedtke convection:
            !   water fraction for mixed water-ice phase as dep. on temperature
          IF (ztt <= tmpmin+epsy) THEN
            fr_wat = c0
          ELSEIF (ztt >= tmpmax-epsy) THEN
            fr_wat = c1
          ELSE
            fr_wat = ((ztt - tmpmin) /(tmpmax - tmpmin)) **exp_mp
          ENDIF
            ! saturation vapour pressure over water
          es_w = fpvsw( ztt, b1, b2w, b3, b4w )
            ! saturation vapour pressure for water / mixed / ice phase
          es_m =       fr_wat  *fpvsw( ztt, b1, b2w, b3, b4w )                 &
                 + (c1-fr_wat) *fpvsi( ztt, b1, b2i, b3, b4i )
            ! bias-corrected observed vapour pressure
          es_o = zrhcorr *es_w
          IF (es_o >= es_m -epsy) THEN
          ! if cloudy level: if previous level also cloudy then adjust 'f_cloud'
          !              and store pressure of cloudy level for next cycle
            zqc =      fr_wat  * 0.1_wp                                        &
                  +(c1-fr_wat) *(0.01_wp + 0.09_wp*                            &
                             MIN( c1, MAX( c0,(zpp-300._wp)/400._wp ) ))
            IF (zppcl > epsy) THEN
            ! to approximate qc_mean, just take mean of qc at the 2 obs levels
              zqc1 = c05 *(zqc1 + zqc)
            ! liquid water path = qc_mean [g/kg] * 100* delta_p [hPa] / g [m/s2]
              zlwp = 100._wp * zqc1 * (zpp - zppcl) / r_g
            ! Tnew = Told *exp( - LWP / 60 ) * (cos(zsolzen/2))^2
            !        but the factor (cos(zsolzen/2))^2 may be applied only once
              IF (f_cloud >= c1-epsy)  f_cloud = f_cloud * r_fsol
              f_cloud = f_cloud *EXP( - zlwp / zlwpcal )
            ENDIF
            zppcl = zpp
            zqc1  = zqc
          ELSE
            zppcl = -c1
            zqc1  =  c0
          ENDIF

!         WRITE(0,'("F_CLD ",A ,I4,F6.0,F6.1,F5.2,3F8.3,F5.2,F7.1,F6.0,F7.3)') &
!               yomlhd(nob), klev, zpp, ztt, fr_wat, es_w, es_m, es_o          &
!                          , zqc1, zlwp, zppcl, f_cloud
        ENDIF

! fill ODR
! --------
        omlbdy (nob,klev,nbtrh ) = zrhcorr
        IF ((zrhcorr > rmdich) .AND. (zrhorig > rmdich)) THEN
          IF (omlbdy(nob,klev,nbtdrh) > rmdich) THEN
            omlbdy (nob,klev,nbtdrh) =   omlbdy(nob,klev,nbtdrh)               &
                                       + zrhcorr - zrhorig
          ELSE
            omlbdy (nob,klev,nbtdrh) =   zrhcorr - zrhorig
          ENDIF
          mflgc  =  NINT( c05 + SIGN( c05, ABS(zrhcorr-zrhorig) -epsy ) )
          CALL MVBITS( mflgc, 0, 1, momlbd(nob,klev,nbtflg), nvfqbp+nvfbps(6) )
        ENDIF
      ENDDO  ! loop over levels
    ENDIF    ! Vaisala RS92 report

! set instrument type to correct value
! ------------------------------------
    IF (momlhd(nob,nhrtyp) /= imdi)                                            &
      momlhd (nob,nhrtyp) = ABS( momlhd(nob,nhrtyp) )

  ENDDO      ! loop over multi-level reports

!-------------------------------------------------------------------------------
! omitted from experimental version:
! Section 3: Rejection of daytime temperature below 800 hPa
!            and of daytime surface pressure                (multi-level obs)
! Section 4: Solar zenith angle at obs time + location (single-level obs)
! Section 5: Rejection of daytime temperature below 800 hPa
!            and of daytime surface pressure                (single-level obs)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! End Subroutine obs_cdf_raso_rh_bias
!-------------------------------------------------------------------------------
END SUBROUTINE obs_cdf_raso_rh_bias


!===============================================================================

ELEMENTAL REAL (KIND=wp) FUNCTION rmod  ( ztti, ziv )
  !----------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  ztti, ziv
  !---------------------------------------------------------------------------
  ! MOD function for positive REALS
  !---------------------------------------------------------------------------
  !
  rmod  =  ztti  -  ziv * REAL(INT( ztti/ziv + epsy ), wp) + epsy
  !
END FUNCTION rmod

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsw  ( zt, b1, b2w, b3, b4w )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2w, b3, b4w
  !---------------------------------------------------------------------------
  ! Magnus formula for water:  input  'zt'   : temperature
  !                            output 'fpvsw': saturation water vapour pressure
  !---------------------------------------------------------------------------
  !
  fpvsw  =  b1 * EXP( b2w *(zt-b3) /(zt-b4w) )
  !
END FUNCTION fpvsw

!-------------------------------------------------------------------------------

ELEMENTAL REAL (KIND=wp) FUNCTION fpvsi  ( zt, b1, b2i, b3, b4i )
  !----------------------------------------------------
  REAL    (KIND=wp)         , INTENT (IN)  ::  zt, b1, b2i, b3, b4i
  !---------------------------------------------------------------------------
  ! Magnus formula for ice:  input  'zt'   : temperature
  !                          output 'fpvsi': saturation water vapour pressure
  !---------------------------------------------------------------------------
  !
  fpvsi  =  b1 * EXP( b2i *(zt-b3) /(zt-b4i) )
  !
END FUNCTION fpvsi

!-------------------------------------------------------------------------------

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

ELEMENTAL INTEGER FUNCTION ibflip  ( invar, ibp )
  !---------------------------------------------------------------------
  INTEGER         , INTENT (IN)  ::  invar, ibp
  !-----------------------------------------------------------------------------
  ! reverts (i.e. flips, toggles) the bit at position 'ipb' of int word 'invar'
  !-----------------------------------------------------------------------------
  !
  ibflip = IEOR( invar , ISHFT( 1 , ibp ) )
  !
END FUNCTION ibflip

!-------------------------------------------------------------------------------

ELEMENTAL INTEGER FUNCTION ileap  ( iyy )
  !--------------------------------------
  INTEGER         , INTENT (IN)  ::  iyy    ! year [YYYY]
  !---------------------------------------------------------------------------
  ! detects leap year:  ileap(iyy) = 0:  no leap year,
  !                     ileap(iyy) = 1:  leap year
  !---------------------------------------------------------------------------
  !
  ileap =   IABS( MOD(iyy,  4) -  4) /  4   & ! every     4 years is a leapyear
          - IABS( MOD(iyy,100) -100) /100   & ! but every 100 ys. is no leapyear
          + IABS( MOD(iyy,400) -400) /400     ! but every 400 ys. is a leapyear
  !
END FUNCTION ileap

!-------------------------------------------------------------------------------

END MODULE src_obs_cdfin_org

