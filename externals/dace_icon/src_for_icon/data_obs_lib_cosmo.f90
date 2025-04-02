!+ Data module for common variables for obs process. and obs operator libraries
!-------------------------------------------------------------------------------

MODULE data_obs_lib_cosmo

!-------------------------------------------------------------------------------
!
! Description:
!   This module contains those variables (scalars and arrays) which must be
!   commonly available in the observation pre-processing library and the
!   observation operator library for application related to the COSMO model.
!   This includes:
!    - some general parameters
!    - scalar variables and arrays originally defined in data modules outside of
!      these libraries, obtained by calling 'obs_cdf_interface' (scalars) or by
!      assigning pointers to target arrays (arrays) in the model environment
!    - tables for pressure dependent scales and geometry of horiz. correlations
!    - I/O device numbers for obs processing / nudging and file names
!    - CMA observation type and code type numbers
!    - data type with meta information / rules for CMA obs and code types
!    - a function on this data type (to obtain array index for given CMA types)
!
!   Note: This module belongs to a group of COSMO-related modules for reading
!   ----  conventional data from ('cdfin') NetCDF observation input files and
!         for observation operators including quality control for conventional
!         observations.
!         It is shared between the COSMO model and the DACE program package !
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
! V4_22        2012/01/31 Christoph Schraff
!  Initial version, partly extracted from module 'data_nudge_all' and extended,
!  e.g. by data types with meta information / rules for CMA obs and code types.
!  A separate CMA code type is introduced for each processing type from each
!  GPS processing centre.
! V4_27        2013/03/19 Christoph Schraff
!  Character length of (y)ydate_ref increased from 10 to 14.
! V4_28        2013/07/12 Christoph Schraff
!  'luse_mlz', 'irun_osse', 'losse_fg', 'fperturb', and 'iseed' added.
! V5_1         2014-11-28 Christoph Schraff, Oliver Fuhrer
!  Observation code type introduced for Mode-S aircraft data. (CS)
!  New value for 'nacar' (145) to be compatible with OC_ACARS in 3DVar. (CS)
!  Replaced ireals by wp (working precision) (OF)
! V5_4d        2016-12-12 Christoph Schraff
!  Additional limits 'zlimv10' for the use of 10-m wind observations.
!  Surface roughness length added to interface of observation reading.
! V5_4f        2017-09-01 Christoph Schraff
!  Update frmo DACE V1_51: 'ifndef' conditioned introduction of 't_grid' and
!  'gen_grid' to make it compatible for ICON grid / for use in ICON.
! V5_6         2019-02-27 Christoph Schraff
!  Switch 'lredn_repro' introduced for constraining the redundancy check
!    in order to make it reproducible.
! V5_6b        2019-10-16 Christoph Schraff
!  - Lists 'r_av_levs', 'r_av_incr', 'av_reso', 'nref', 'pref' added for
!    superobbing of high-resolution radiosonde profile data.
!  - Observation code types 'ntdesc' (high-res. BUFR descending TEMP),
!    'nbtemp' (high-res. BUFR land temp), 'nbship' (high-res. BUFR ship temp),
!    'nbdrop' (high-res. BUFR temp drop), and 'ntower' (tower) added.
!  - Unused code types 'nabscd', 'nshred' (abbrev. ship), 'ncolba', 'nsst',
!    and 'nbathy' removed. 
!  - GPS centre 'IGN_' added.
!  - Set parameters and variables used by other modules explicitly to public
!    (required by ICON Compiler directives when modules are ported to ICON).
!  - All 'iintegers' removed and replaced by standard integers.
!  - Use of 'kind_parameters' instead of 'data_parameters'.
! V5_7a        2020-05-11 Christoph Schraff
!  Observation code types 'nwlidr' (ground-based wind lidar) and 'nsytst'
!  (test synop surface report) added.
!  - 16.05.22: 'r_sso_sd' added (std. dev. of sub-grid scale model orography).
!  - 17.04.23: fixed parameters converted into namelist variables added here: 
!              'rhtsat' (limit above which rel. humidity obs set to 100%), and
!              redundancy time limits 'rtmlrs(y)', 'rtmlsy', 'rtmlair','rtmlim'.
!  - 24.04.23: Observation code type 'ntowic' (ICOS tower profiles)  added.
!  - 16.05.23: Namelist parameters 'itim_wp', 'icdt_tws', 'rtmltow' added.
!  - 07.10.24: Namelist parameters 'icdt_rss', 'lsytac' and 
!              code types 'nsytmp', 'nsytow' added.
!
! ==============================================================================
! DACE:
! ----
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
! Modules used:
!
!-------------------------------------------------------------------------------

USE kind_parameters, ONLY :   &
    wp           ! KIND-type parameter for real variables

USE data_constants,  ONLY :   &
    rprecision   ! precision in current floating point format

#ifndef __COSMO__

  USE mo_atm_grid, ONLY :  &
    t_grid             ! generic grid meta data type (also for ICON)

#endif

!-------------------------------------------------------------------------------

IMPLICIT NONE

PRIVATE
          ! 1    General parameters
PUBLIC :: c0, c1, c2, c05, c3600, rmdi, rmdich, epsy, i1, luse_mlz,            &
          ! 2    Variables related to cdf_interface
          ie_loc, je_loc, ke, ie_tot, je_tot, num_compute, nboundlines,        &
          my_cart_id, icomm_cart, imp_reals, imp_integers, imp_character,      &
          maxmlv, maxmll, maxsgl, maxgpl, nolbc, itim_wp, icdt_tws, icdt_rss,  &
          madj_hum, irun_osse, iseed, mxgpc, lverpas, losse_fg, lredn_repro,   &
          lsytac,                                                              &
          r_pollon, r_pollat, r_polgam, r_dlon, r_dlat,                        &
          r_startlat_tot, r_startlon_tot, r_degrad,                            &
          rhtsat, rtmlrs, rtmlsy, rtmlair, rtmltow, rtmlrsy, rtmlim,           &
          doromx, altopsu, zlimv10, fperturb, acthr, lwonl, ydate_ref,         &
          r_g, tmelt, r_d, rdv, o_m_rdv, rdocp, b1, b2w, b2i, b3, b4w, b4i,    &
          i_subpos, i_gpscen,                                                  &
          r_p, r_hhl, r_t_ll, r_ps, r_frland, r_sso_sd, r_z0,                  &
          r_av_levs, r_av_incr, mxav, av_reso, pref, nref,                     &
          ! 3    Variables related to cdf_interface
          ncolev, tabcolp, tabcop, rhvsond, rhtsond, rhqsond,                  &
          ! 4    I/O files
          yucautn, yustats, yurejct, yuobsdr, yuprint,                         &
          nucautn, nustat , nurej  , nuodr  , nupr   , lopen_odr, lopen_rej,   &
          ! 5    CMA observation and code types
          nsynop, nairep, nsatob, ndribu, ntemp , npilot, nsatem, nsattv,      &
          nscatt, ngps  ,                                                      &
          nsrscd, natscd, nshscd, natshs, nmetar, nsytst, nsytmp, nsytow,      &
          naircd, ncodar, namdar, nacar , nmodes,                              &
          nstbcd, nhrvis, namv  ,         ndrbcd, ntesac, nrocob, nrocsh,      &
          nldtcd, nshtcd, nmotcd, nbtemp, nbship, nbdrop, ntdrop, ntdesc,      &
          nldpcd, nshpcd, nmopcd, nwp_eu, nra_eu, npr_us, nravad, ntower,      &
          ntowic, nwlidr,         nascat, nqscat,         nstmcd, nstovs,      &
          ! 6    Data types with meta info / rules
          t_gpscen, n_gpc, gpc, ngpgfz, t_cmatyp, cs, cn, cw, ce, cma, n_cma,  &
          n_gpc_offs, i_cma

! for nudging in COSMO model code only:
PUBLIC :: fnodivc, yuquctl, yuverif, nuqc, nuverif, mxobtp, nacar_vof

! for nudging in 1DVAR only:
! PUBLIC :: r_t, r_qv, r_t_g, r_t_2m, r_qv_2m, r_u_10m, r_v_10m
! PUBLIC :: nsmsg1, nsmsg2, nnoa15, nnoa16, nnoa17, nnoa18, nusatin , maxtvl

#ifndef __COSMO__
PUBLIC :: gen_grid
#endif

!===============================================================================

! Local Declarations:

!-------------------------------------------------------------------------------
! Global (i.e. public) Declarations:
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Section 1:  General parameters
!-------------------------------------------------------------------------------

  REAL (KIND=wp) , PARAMETER  :: &
    c0     = 0.0_wp     ,& ! standard real constant 0.0
    c1     = 1.0_wp     ,& ! standard real constant 1.0
    c2     = 2.0_wp     ,& ! standard real constant 2.0
    c05    = 0.5_wp     ,& ! standard real constant 0.5
    c3600  = 3600.0_wp  ,& ! standard real constant 3600.0
    rmdi   = -1.E31_wp  ,& ! commonly used missing data indicator
    rmdich = -1.E30_wp  ,& ! commonly used check value for missing data
    epsy   = MAX(1.0E-8_wp ,& ! commonly used very small value > 0
                 rprecision)  ! (limited to rprecision in single precision)

  INTEGER        , PARAMETER  :: &
    i1     = 1             ! standard integer constant 1

  LOGICAL        , PARAMETER  :: &
    luse_mlz  =  .FALSE.   ! prepare for possible future namelist parameter:
                           ! if .false. then multi-level temperature rather than
                           !   height is used, and height obs are set passive
                           !   except for the lowest z-p-obs

!-------------------------------------------------------------------------------
! Section 2:  Scalar variables originally defined in other data modules
!             (which belong to the model environment),
!             values are obtained by calling 'obs_cdf_interface'
!-------------------------------------------------------------------------------

! horizontal and vertical sizes of the model fields
    ! ie_loc, je_loc: used in 'get_global_surface', 'obs_cdf_interface'
    ! ke            : used in 'f_z2p'
    ! ie_tot, je_tot: used in 'get_global_surface', 'obs_cdf_read_org',
    !                         'obs_assign_gridpt'

  INTEGER        :: &
    ie_loc         ,& ! number of grid pts in zonal direction (local sub-domain)
    je_loc         ,& ! number of grid pts in meridional dir. (local sub-domain)
    ke             ,& ! number of grid pts in vertical direction (--> 'f_z2p')
    ie_tot         ,& ! number of grid pts in zonal direction (in total domain)
    je_tot            ! number of grid pts in meridional dir. (in total domain)

! variables related to parallelisation / domain decomposition

  INTEGER        :: &
    num_compute    ,& ! number of compute PEs
    nboundlines    ,& ! number of overlapping boundary lines of the subdomains
    my_cart_id     ,& ! rank of this subdomain in the cartesian communicator
    icomm_cart     ,& ! communicator for the virtual cartesian topology
    imp_reals      ,& ! REAL      type used for MPI
    imp_integers   ,& ! INTEGER   type used for MPI
    imp_character     ! CHARACTER type used for MPI

! report dimensions of ODR (Observation Data Record, for observation storage on
! local sub-domains), and other variables related to namelist parameters

  INTEGER        :: &
                      !  (guarantee initialization of report dimensions)
    maxmlv = 1     ,& ! size (level  dimension) of the  multi-level (m-l)  ODR
    maxmll = 1     ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl = 1     ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl = 1     ,& ! size (report dimension) of the (ground-based) GPS  ODR
!   maxtvl = 1     ,& ! size (report dimension) of the satellite retrieval ODR
                      !
    nolbc          ,& ! number of grid rows at lateral boundaries
                      !   where obs are neglected
    itim_wp        ,& !   0 : mode of correction of obs time for wind profiler
                      !       = 0 : no correction
                      !       = 1 : correct by half of the obs averaging period
                      !       = 2 : correct to the end of the obs aver. period
                      !       = 3 : correct to the beginning of the obs period
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
    madj_hum       ,& ! = 1 : adjust observed humidity (by ratio of saturation
                      !       vapour pressure over water to the one over ice,
                      !       to be applied if cloud ice is not a state variable
    irun_osse      ,& !   0 : for data from input feedback file yfofin='fof':
                      !       == 0 : obs from yfofin are used as obs
                      !       >= 1 : simulated obs from yfofin with model run
                      !              index = irun_osse are used as obs (thus,
                      !              run 'irun_osse' is nature run for OSSE)
    iseed          ,& !   0 : external seed for random number generator
                      !       (this seed is combined with another seed that
                      !        depends on the initial date and time of the
                      !        model run;  note that for LETKF OSSE, 'iseed'
                      !        should be identical for all ensemble members
                      !        so that the perturbed obs will be identical)
                      !
    mxav           ,& ! max. length of level definition list for superobbing
                      !       of high-resolution radisondes
    mxgpc             ! max. number of GPS processing centres used

  LOGICAL        :: &
    lverpas        ,& ! write also passive reports on feedback files
    losse_fg       ,& ! .f. : if true then first guess check flag from 'fof' is
                      !             converted into 'dataset' pre-processing flag
                      !       if false then first guess check flag is discarded
                      !             and obs may be used actively
    lredn_repro    ,& ! .f. : if true then reproducibility of redundancy check
                      !       ensured irrespective of domain decomposition by
                      !       allowing for redundancy only between reports
                      !       assigned to the same grid point
    lsytac            ! .t. : prefer TAC over BUFR Synop reports (or vice versa)

! constants for the horizontal rotated grid and related variables

  REAL (KIND=wp) :: &
    r_pollon       ,& ! longitude of the rotated north pole (in degrees, E>0)
    r_pollat       ,& ! latitude of the rotated north pole (in degrees, N>0)
    r_polgam       ,& ! angle between the north poles of the systems
    r_dlon         ,& ! grid point distance in zonal direction (in degrees)
    r_dlat         ,& ! grid point distance in meridional direction (in degrees)
    r_startlat_tot ,& ! transformed latitude of the lower left grid point
                      ! of the total domain (in degrees, N>0)
    r_startlon_tot ,& ! transformed longitude of the lower left grid point
                      ! of the total domain (in degrees, E>0)
    r_degrad          ! factor for transforming degree to rad

! generic grid meta data for ICON and other models

#ifndef __COSMO__

  TYPE (t_grid) ,pointer :: &
    gen_grid          ! grid meta data

#endif

! variables related to namelist parameters

  REAL (KIND=wp) :: &
    doromx     (4) ,& ! SYNOP obs. with height differences betw. model orography
                      !  and station height larger than 'doromx' are set passive
    altopsu    (4) ,& ! SYNOP obs. above height 'altopsu' are set passive
    zlimv10    (3) ,& ! additional limits for the use of 10-m wind observations:
                      !   pos. resp. neg. scaled Laplacian of orography [m],
                      !   surface roughness length [m]
    fperturb       ,& ! factor to obs error variances to define size of random
                      ! perturbations added to the obs (only from 'fof' file)
    av_reso        ,& ! apply superobbing if the averaged resolution of the obs
                      !   profile exceeds 'av_reso' times the model resolution
    rhtsat     (5) ,& ! relative humidity threshold above which obs is set =100%
                      !   (1: all, 2: surface obs, 3: TEMP, 4: PILOT, 5: AIREP;
                      !    non-default values 2 - 5 overrule value at index 1)
    rtmlrs         ,& ! redundancy time limit for radiosondes              [hrs]
    rtmlsy         ,& ! redundancy time limit for SYNOP                    [hrs]
    rtmlair        ,& ! redundancy time limit for AIREP                    [hrs]
    rtmltow        ,& ! redundancy time limit for tower                    [hrs]
    rtmlrsy        ,& ! redundancy time limit for raso sfc. level vs Synop [hrs]
    rtmlim         ,& ! redundancy time limit for other obs                [hrs]

! time-dependent variable

    acthr             ! actual model time [hours] with respect to 'ydate_ref'

! physical constants

  REAL (KIND=wp) :: &
    r_g            ,& ! acceleration due to gravity    (  9.80665)
    tmelt          ,& ! melting temperature of ice     (273.15)
    r_d            ,& ! gas constant for dry air       (287.05)
!   r_v            ,& ! gas constant for water vapor   (461,51)
    rdv            ,& ! r_d / r_v
    o_m_rdv        ,& ! 1 - r_d/r_v
!   cp_d           ,& ! specific heat of dry air at constant pressure  (1005.0)
    rdocp          ,& ! r_d / cp_d
    b1             ,& ! variables to compute saturation vapour pressure (610.78)
    b2w            ,& ! over water (w) and ice (i)                 ( 17.2693882)
    b2i            ,& !               -- " --                      ( 21.8745584)
    b3             ,& !               -- " --                      (273.16)
    b4w            ,& !               -- " --                      ( 35.86)
    b4i               !               -- " --                      (  7.66)

! switches related to namelist parameters and other

  LOGICAL        :: &
    lwonl             ! .true. for the node (sub-domain) at which file with
                      !        the unit number 'nupr' is open
!   lgpsbias          ! .t. ==> bias correction to GPS IWV applied

  CHARACTER (LEN= 14)  ::  &
    ydate_ref         ! reference date (e.g. start of the forecast)
                      ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

!-------------------------------------------------------------------------------
! Section 2b: Pointers for arrays originally defined in other data modules
!             (which belong to the model environment),
!             values are obtained by assigning pointers to target arrays
!             and takes place in the model environment
!-------------------------------------------------------------------------------

! arrays related to parallelisation / domain decomposition

  INTEGER        , POINTER    :: &
    i_subpos (:,:) => NULL() ! positions of the subdomains in the total domain
                             ! (i- and j-indices of the lower left + upper right
                             !  grid point in the order: i_ll, j_ll, i_ur, j_ur;
                             !  only the domain interior is considered, not the
                             !  boundary lines)

! arrays related to namelist parameters

  INTEGER        , POINTER    :: &
    i_gpscen (:) => NULL()   ! array of processing centres of GPS reports used
                             ! actively: order of centres determines preference
                             ! in redundancy check; '-1' means no active centre

  REAL (KIND=wp) , POINTER    :: &
    r_av_levs(:) => NULL(),& ! level definition list      \  for superobbing
                             !  (1075,755,710, 90, 75, 5)  \ layers of
    r_av_incr(:) => NULL()   ! level increment  list       / high-resolution
                             !  (  10, 15, 20, 15, 10, 0) /  radiosonde reports

! model fields

  REAL (KIND=wp) , POINTER    :: &
    r_p    (:,:,:) => NULL()  ,& ! pressure (at main levels)
    r_hhl  (:,:,:) => NULL()  ,& ! height   (at half levels)
    r_t_ll   (:,:) => NULL()  ,& ! temperature at lowest model (main) level
    r_ps     (:,:) => NULL()  ,& ! surface pressure
    r_frland (:,:) => NULL()  ,& ! land fraction
    r_sso_sd (:,:) => NULL()  ,& ! std. dev. of sub-grid scale model orography
    r_z0     (:,:) => NULL()     ! surface roughness length

!-------------------------------------------------------------------------------
! Section 2c: Pointers for arrays originally defined in other data modules,
!             values obtained by assignment to target arrays (1DVAR only)
!-------------------------------------------------------------------------------

! REAL (KIND=wp) , POINTER    :: &
!!                               !   model fields:
!   r_t    (:,:,:) => NULL()  ,& ! temperature                           (  K  )
!   r_qv   (:,:,:) => NULL()  ,& ! specific water vapor content          (kg/kg)
!   r_t_g    (:,:) => NULL()  ,& ! weighted surface temperature          (  K  )
!   r_t_2m   (:,:) => NULL()  ,& ! temperature in 2m                     (  K  )
!   r_qv_2m  (:,:) => NULL()  ,& ! specific water vapor content in 2m    (kg/kg)
!   r_u_10m  (:,:) => NULL()  ,& ! zonal wind in 10m                     ( m/s )
!   r_v_10m  (:,:) => NULL()     ! meridional wind in 10m                ( m/s )

!-------------------------------------------------------------------------------
! Section 2d: Quantities related to above pointers
!-------------------------------------------------------------------------------

  REAL (KIND=wp) , ALLOCATABLE :: &
    pref       (:)    ! reference interpolation pressure levels for
                      ! superobbing of high-resolution radiosonde reports

  INTEGER        :: &
    nref              ! actual number of interpolation reference levels

!-------------------------------------------------------------------------------
! Section 3:  Tables for pressure dependent scales and geometry
!             of horizontal correlations
!             (used in 'obs_air_correl' and for spreading)
!-------------------------------------------------------------------------------

!      3.1  Pressure levels used for the tables
!           -----------------------------------

  INTEGER        , PARAMETER  :: &
    ncolev = 11    ! number of levels in the correlation scale tables

  REAL (KIND=wp) :: &
    tabcolp (11)   ! ln(tabcop(11))

  REAL (KIND=wp) , PARAMETER  :: &
    tabcop  (11) & ! levels of the correlation scale tables
              = (/ 100000._wp, 85000._wp, 70000._wp, 50000._wp,&
                    40000._wp, 30000._wp, 25000._wp, 20000._wp,&
                    15000._wp, 10000._wp,  5000._wp /)

!      3.2  Pressure dependent (part of) horizontal correlation scales
!      ---------------------------------------------------------------

!           defined at the standard levels tabcop(1:ncolev)

  REAL (KIND=wp) , PARAMETER  :: &
                   ! upper-air     wind      horizontal correlation scales
    rhvsond (11) = (/ 70.0_wp,   80.0_wp,   90.0_wp,  100.0_wp,&
                     100.0_wp,  110.0_wp,  115.0_wp,  120.0_wp,&
                     125.0_wp,  125.0_wp,  125.0_wp /)            ,&
                   ! upper-air  temperature  horizontal correlation scales
    rhtsond (11) = (/ 70.0_wp,   80.0_wp,   90.0_wp,  100.0_wp,&
                     100.0_wp,  100.0_wp,  100.0_wp,  110.0_wp,&
                     120.0_wp,  120.0_wp,  120.0_wp /)            ,&
                   ! upper-air   humidity    horizontal correlation scales
    rhqsond (11) = (/ 70.0_wp,   80.0_wp,   90.0_wp,  100.0_wp,&
                     100.0_wp,  100.0_wp,  100.0_wp,  110.0_wp,&
                     120.0_wp,  120.0_wp,  120.0_wp /)

!      3.3  Pressure dependent factors to the non-divergence
!           correction in the 2-dim wind correlations
!           ------------------------------------------------

!           defined at the standard levels tabcop(1:ncolev)

  REAL (KIND=wp) , PARAMETER  :: &
                   ! factor to non-divergence correction in the 2-d wind correl.
    fnodivc (11) = (/ 0.40_wp,   0.50_wp,   0.50_wp,   0.50_wp,&
                      0.50_wp,   0.60_wp,   0.65_wp,   0.70_wp,&
                      0.75_wp,   0.75_wp,   0.75_wp /)

!-------------------------------------------------------------------------------
! Section 4:  I/O device numbers for obs processing / nudging and file names
!-------------------------------------------------------------------------------

!      4.1  File names
!           ----------

  CHARACTER (LEN=7)    , PARAMETER  :: &
    yucautn = 'YUCAUTN'  ,& ! caution messages if too many obs for ODR size
    yuquctl = 'YUQUCTL'  ,& ! data rejected by threshold QC at obs time
    yustats = 'YUSTATS'  ,& ! statistics of processed reports
    yurejct = 'YUREJCT'  ,& ! direct reporting of rejected obs. reports
    yuobsdr = 'YUOBSDR'  ,& ! observations stored in the observation data record
    yuverif = 'YUVERIF'  ,& ! VOF (output): verification observation file (obs.
                            !      incl. quality control flag for verification)
    yuprint = 'YUPRINT'     ! all the remaining information

!      4.2  Device numbers   (obtained by call of 'get_free_unit'
!           --------------    outside of the library)

  INTEGER        :: &
!   nusatin ,& ! satellite meta data input files (bias correction,
               !                                  channel selection)
    nucautn ,& ! caution messages if too many obs for current ODR size
    nuqc    ,& ! data rejected by threshold quality control at obs time
    nustat  ,& ! statistics of processed reports
    nurej   ,& ! direct reporting of rejected obs. reports
    nuodr   ,& ! observations stored in the observation data record
    nuverif ,& ! VOF (output): verification observation file (observations
               !   incl. quality control flag for verification)
    nupr       ! all the remaining information

!      4.3  File status
!           -----------

  LOGICAL        :: &
    lopen_odr = .FALSE. ,& ! .true. if yuobsdr is open
    lopen_rej = .FALSE.    ! .true. if yurejct is open

!-------------------------------------------------------------------------------
! Section 5:  CMA observation type and code type numbers
!             Note: These tables should correspond to the tables 'obstype' and
!                   'codetype' in module 'mo_fdbk_tables.f90' !!!
!-------------------------------------------------------------------------------

!         5.1   Observation type numbers
!               ------------------------

  INTEGER        , PARAMETER  :: &
    nsynop =  1  ,& !   SYNOP report
    nairep =  2  ,& !   AIREP report (all aircraft reports)
    nsatob =  3  ,& !   SATOB report
    ndribu =  4  ,& !   DRIBU report
    ntemp  =  5  ,& !   TEMP  report
    npilot =  6  ,& !   PILOT report
    nsatem =  7  ,& !   SATEM report
    nsattv =  7  ,& !   SATEM report
    nscatt =  9  ,& !   SCATT report (from NetCDF only)
    ngps   =  12 ,& !   GPS   report
    mxobtp =  12    ! number of observation types (should be ' = mxot_blk ' in
                    !                              module 'data_obs_cdfin.f90')

!         5.2   Observation code type numbers
!               -----------------------------

  INTEGER        , PARAMETER  :: &
    nsrscd =  11 ,& !   synop surface report
    natscd =  14 ,& !   automatic synop surface report
    nshscd =  21 ,& !   ship synop report
!   nabscd =  22 ,& !   ship synop abbreviated report
!   nshred =  23 ,& !   shred report
    natshs =  24 ,& !   automatic ship synop report
    nmetar = 140 ,& !   Metar
    nsytst = 811 ,& !   test synop surface report
    nsytmp = 835 ,& !   surface report from TEMP
    nsytow = 839 ,& !   surface report from TEMP
    naircd = 141 ,& !   aircraft report
    ncodar =  41 ,& !   codar report
!   ncolba = 241 ,& !   colba report
    namdar = 144 ,& !   amdar report
    nacar  = 145 ,& !   acar  report
!   nacar  = 244 ,& !   acar  report
    nmodes = 146 ,& !   mode-s report
    nacar_vof=244,& !   acar  report (only in VOF)
    nstbcd =  88 ,& !   satob report
    nhrvis =  89 ,& !   high-res. VIS wind report
    namv   =  90 ,& !   AMV   report
!   nsst   = 188 ,& !   sst report
    ndrbcd = 165 ,& !   dribu report
!   nbathy =  63 ,& !   bathy report
    ntesac =  64 ,& !   tesac report
    nldtcd =  35 ,& !   temp land   report
    nshtcd =  36 ,& !   temp ship   report
    nmotcd =  37 ,& !   temp mobile report
    ntdrop = 135 ,& !   temp drop   report
    nbtemp = 109 ,& !   temp land   report (high-res. BUFR)
    nbship = 111 ,& !   temp ship   report (high-res. BUFR)
    nbdrop = 230 ,& !   temp drop   report (high-res. BUFR)
    ntdesc = 231 ,& !   temp high-res. descending report (HA 02.04.19: as EC)
    nrocob =  39 ,& !   rocob      report
    nrocsh =  40 ,& !   rocob ship report
    nldpcd =  32 ,& !   pilot land   report
    nshpcd =  33 ,& !   pilot ship   report
    nmopcd =  38 ,& !   pilot mobile report
    nwp_eu = 132 ,& !   wind profiler report (European)
    nra_eu = 133 ,& !   sodar/rass report (European)
    npr_us = 136 ,& !   wind profiler/rass report (USA)
    nravad = 137 ,& !   radar vad wind profile report
    ntower = 139 ,& !   tower profile report
    ntowic = 159 ,& !   icos tower profile report
    nwlidr = 187 ,& !   ground-based wind lidar report
    nascat = 123 ,& !   ASCAT scatterometer report
    nqscat = 122    !   QuickScat scatterometer report
!   ndrone = 147 ,& !   drone report
!   nraman = 190    !   Raman lidar report (> 190)

  INTEGER        , PARAMETER  :: &
    nstmcd =  86 ,& !   satem report
    nstovs = 186    !   high resolution ATOVS satellite data
!HE    nsmsg1 =  55 ,& !   MSG_1 (METEOSAT-8) satellite retrieval
!HE    nsmsg2 =  56 ,& !   MSG_2 (METEOSAT-9) satellite retrieval
!   nsmsg1 =  71 ,& !   MSG_1 (METEOSAT-8) satellite retrieval
!   nsmsg2 =  72 ,& !   MSG_2 (METEOSAT-9) satellite retrieval
!   nnoa15 = 206 ,& !   NOAA15 satellite retrieval
!   nnoa16 = 207 ,& !   NOAA16 satellite retrieval
!   nnoa17 = 208 ,& !   NOAA17 satellite retrieval
!   nnoa18 = 209    !   NOAA18 satellite retrieval
!   ngpgfz  --> is set further below (required only for reading GPS data
!                                     from COST-716 ASCII files)
!   ngpgfz =  96 ,& !   GPS report processed by GFZ  (Geoforschungszentrum, Germany)
!   ngpasi =  97 ,& !   GPS report processed by ASI  (Agenzia Spaziale Italiana, Italy)
!   ngpbkg =  98 ,& !   GPS report processed by BKG  (Bundesamt fuer Kartographie und Geodaesie, Germany)
!   ngpgop =  99 ,& !   GPS report processed by GOP  (Geodetic Observatory Pecny, Czech Republic
!   ngpieec= 100 ,& !   GPS report processed by IEEC (Institut d'Estudis Espacials de Catalunya, Spain)
!   ngpsgn = 101 ,& !   GPS report processed by SGN  (Institute de Geodesie National, France)
!   ngplpt = 102 ,& !   GPS report processed by LPT  (Swiss Federal Office of Topography, Swiss)
!   ngpmet = 103 ,& !   GPS report processed by MET  (UK Met Office, UK)
!   ngprob = 104 ,& !   GPS report processed by ROB  (Royal Observatory of Belgium, Belgium)
!   ngpige = 105 ,& !   GPS report processed by IGE  (Instituto Geografico Nacional de Espana, Spain)
!   ngpknmi= 106 ,& !   GPS report processed by KNMI (Koninklijk Nederlands Meteorologisch Instituut, Netherlands)
!   ngpnga = 107    !   GPS report processed by NGA  (Nordic GPS Atmospheric Analysis Centre, Sweden)

!-------------------------------------------------------------------------------
! Section 6:  Data type with meta information / rules for CMA obs and code types
!-------------------------------------------------------------------------------

! data type to relate number and name of GPS processing center and CMA code type
! ------------------------------------------------------------------------------
! (note: in order to add the processing / use of further GPS processing centres,
!        only (!) the table 'gpc' and its dimension 'n_gpc' need to be updated)

  TYPE t_gpscen
    CHARACTER (LEN=4)  :: name      ! 4-character centre name, contained
                                    !   in station identifier (char. 6-9)
                                    !   of BUFR reports obtained via UKMO
    INTEGER            :: gpscen    ! (original) processing centre
                                    !   (WMO Common Code Table C-12)
    INTEGER            :: cdtyp     ! CMA code type
  END TYPE t_gpscen

  INTEGER        , PARAMETER  ::  n_gpc = 21

  TYPE (t_gpscen), PARAMETER  ::  gpc (n_gpc) =   (/&
      ! GPS:
    t_gpscen( 'METO' ,  0, 800 ) ,& ! UK: UK Met Office
    t_gpscen( 'MET_' ,  0, 900 ) ,& ! UK: UK Met Office
    t_gpscen( 'ASI_' , 21, 821 ) ,& ! I : Agenzia Spaziale Italiana
    t_gpscen( 'GFZ_' , 23, 823 ) ,& ! D : Geoforschungszentrum
    t_gpscen( 'GOP_' , 24, 824 ) ,& ! CZ: Geodetic Observatory Pecny
    t_gpscen( 'GOPE' , 24, 924 ) ,& ! CZ: Geodetic Observatory Pecny
    t_gpscen( 'IEEC' , 25, 825 ) ,& ! E : Inst. d'Estudis Espacials de Catalunya
    t_gpscen( 'LPT_' , 26, 826 ) ,& ! CH: Swiss Federal Office of Topography
    t_gpscen( 'LPTR' , 26, 926 ) ,& ! CH: Swiss Federal Office of Topography
    t_gpscen( 'SGN_' , 29, 829 ) ,& ! F : Institute de Geodesie National
    t_gpscen( 'SGN1' , 29, 929 ) ,& ! F : Institute de Geodesie National
    t_gpscen( 'IGN_' , 29, 829 ) ,& ! F : Institute de Geodesie National
    t_gpscen( 'BKG_' , 30, 830 ) ,& ! D : Bundesamt Kartographie und Geodaesie
    t_gpscen( 'BKGH' , 30, 930 ) ,& ! D : Bundesamt Kartographie und Geodaesie
    t_gpscen( 'KNMI' , 33, 833 ) ,& ! NL: Koninklijk Nederlands Meteorol. Inst.
    t_gpscen( 'KNM1' , 33, 933 ) ,& ! NL: Koninklijk Nederlands Meteorol. Inst.
    t_gpscen( 'NGAA' , 34, 834 ) ,& ! S : Nordic GPS Atmospheric Analysis Centre
    t_gpscen( 'NGA_' , 34, 934 ) ,& ! S : Nordic GPS Atmospheric Analysis Centre
    t_gpscen( 'IGE_' , 35, 835 ) ,& ! E : Instit. Geografico Nacional de Espana
    t_gpscen( 'ROB_' , 37, 837 ) ,& ! B : Royal Observatory of Belgium
    t_gpscen( 'XXX_' , 99, 899 ) /) ! -   yet unknown processing centres

  INTEGER        , PARAMETER  :: &
    ngpgfz = gpc(n_gpc) %cdtyp      ! Geoforschungszentrum, required for reading
                                    !   from ASCII files in COST-716 format

! data type for information on CMA observation and code types
! -----------------------------------------------------------

  TYPE t_cmatyp
    INTEGER            :: obtyp     ! CMA observation type
    INTEGER            :: cdtyp     ! CMA code type
    CHARACTER (LEN=20) :: name      ! description
    REAL  (KIND=wp)    :: obslat    ! southern latitude limit for active use
    REAL  (KIND=wp)    :: obnlat    ! northern latitude limit for active use
    REAL  (KIND=wp)    :: obwlon    ! western longitude limit for active use
    REAL  (KIND=wp)    :: obelon    ! eastern longitude limit for active use
    REAL  (KIND=wp)    :: exslat    ! southern latitude limit of exclusion area
    REAL  (KIND=wp)    :: exnlat    ! northern latitude limit of exclusion area
    REAL  (KIND=wp)    :: exwlon    ! western longitude limit of exclusion area
    REAL  (KIND=wp)    :: exelon    ! eastern longitude limit of exclusion area
    INTEGER            :: cnt_pr    ! counter: processed reports
    INTEGER            :: cnt_ac    ! counter: active    reports
    INTEGER            :: cnt_ps    ! counter: passive   reports
    INTEGER            :: cnt_rj    ! counter: rejected  reports
    INTEGER            :: attr      ! attribute:
                                    ! -SATEM: =0: do not process this type
                                    !         =1: process SEVIRI radiances
                                    !         =2: process ATOVS radiances
  END TYPE t_cmatyp

  REAL (KIND=wp) , PARAMETER  :: &
    cs  =  -90.0_wp ,& ! southern latitude limit
    cn  =   90.0_wp ,& ! northern latitude limit
    cw  = -180.0_wp ,& ! western longitude limit
    ce  =  180.0_wp    ! eastern longitude limit

! meta data on CMA observation and code types

! split 'cma' into 'cma_p1' and 'cma_p2' to avoid more than 39 continuation
! lines  (construction requires 'paramater' attribute)

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_p0 ( 1) =         (/&
    t_cmatyp( 0     , 0     , 'unknown obs/code typ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 )/)

  INTEGER        , PARAMETER  , PRIVATE  ::  n_cma_p1 = 33

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_p1 (n_cma_p1) =   (/&
      ! SYNOP:
    t_cmatyp( nsynop, 0     , 'SYNOP               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nsrscd, 'SYNOP Manual Land   ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, natscd, 'SYNOP Automatic Land',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nshscd, 'SHIP                ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsynop, nabscd, 'SHIP Abbreviated    ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsynop, nshred, 'SHIP Reduced SHRED  ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, natshs, 'SHIP Automatic      ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nmetar, 'METAR               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nsytst, 'SYNOP TEST          ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nsytmp, 'surface TEMP        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsynop, nsytow, 'surface tower       ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
      ! AIREP:
    t_cmatyp( nairep, 0     , 'AIREP               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nairep, ncodar, 'CODAR               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nairep, naircd, 'AIREP Aircraft      ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nairep, namdar, 'AMDAR               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nairep, nacar , 'ACARS               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nairep, nmodes, 'MODE-S              ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nairep, ndrone, 'DRONE               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nairep, ncolba, 'COLBA Const Lev Ball',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
      ! SATOB
    t_cmatyp( nsatob, 0     , 'SATOB               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsatob, nstbcd, 'SATOB               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsatob, nhrvis, 'High-res VIS Wind   ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsatob, namv  , 'AMV                 ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsatob, nsst  , 'SST as DRIBU        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
      ! DRIBU:
    t_cmatyp( ndribu, 0     , 'DRIBU               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ndribu, nbathy, 'BATHY               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ndribu, ntesac, 'TESAC               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ndribu, ndrbcd, 'DRIBU Drifting Buoy ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
      ! TEMP:
    t_cmatyp( ntemp , 0     , 'TEMP                ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nldtcd, 'TEMP Land           ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nshtcd, 'TEMP Ship           ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nmotcd, 'TEMP Mobile         ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nrocob, 'ROCOB Land          ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nrocsh, 'ROCOB SHIP          ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nbtemp, 'TEMP Land high-res  ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nbship, 'TEMP Ship high-res  ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , nbdrop, 'TEMP DROP high-res  ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , ntdrop, 'TEMP DROP           ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( ntemp , ntdesc, 'TEMP DESCENT hi-res.',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 )/)

  INTEGER        , PARAMETER  , PRIVATE  ::  n_cma_p2 = 17
! INTEGER        , PARAMETER  , PRIVATE  ::  n_cma_p2 = 21

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_p2 (n_cma_p2) =   (/&
      ! PILOT:
    t_cmatyp( npilot, 0     , 'PILOT               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nldpcd, 'PILOT Land          ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nshpcd, 'PILOT SHIP          ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nmopcd, 'PILOT Mobile        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nwp_eu, 'Wind Profiler (Eur) ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nra_eu, 'RASS / SODAR  (Eur) ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, npr_us, 'Wind Prof/RASS (US) ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nravad, 'RADAR VAD Wind Prof.',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, ntower, 'TOWER / MAST        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, ntowic, 'ICOS TOWER / MAST   ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( npilot, nwlidr, 'Wind Lidar ground-b ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
      ! SATEM:
    t_cmatyp( nsattv, 0     , 'SATEM               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsattv, nstmcd, 'GTS SATEM (500km)   ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nsattv, nstovs, 'Hi-res SATEM (250km)',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nsmsg1, 'SEVIRI MSG-1        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nsmsg2, 'SEVIRI MSG-2        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nnoa15, 'ATOVS NOAA15        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nnoa16, 'ATOVS NOAA16        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nnoa17, 'ATOVS NOAA17        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( nsattv, nnoa18, 'ATOVS NOAA18        ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp(  7,200,'GTS BUFR SATEM 250km                     '),&
!   t_cmatyp(  7,201,'GTS BUFR SATEM Clear Radiance            '),&
!   t_cmatyp(  7,202,'GTS BUFR retr. profiles/clear radiances  '),&
!   t_cmatyp(  7,210,'ATOVS                                    '),&
!   t_cmatyp(  7,211,'RTOVS                                    '),&
!   t_cmatyp(  7,212,'TOVS                                     '),&
!   t_cmatyp(  7,215,'SSMI                                     '),&
    ! SCATT:
    t_cmatyp( nscatt, 0     , 'Scatterometer       ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nscatt, nascat, 'ASCAT scatterometer ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
    t_cmatyp( nscatt, nqscat, 'QuickScat scatterom.',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 )/)

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_gps           =     &
    t_cmatyp( ngps  , 0     , 'GPS                 ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 )

  INTEGER        , PARAMETER  , PRIVATE  ::  n_cma_p3 = n_gpc + 1

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_p3 (n_cma_p3) = cma_gps

! INTEGER        , PARAMETER  , PRIVATE  ::  n_cma_p3 = 13
! TYPE (t_cmatyp), PARAMETER  , PRIVATE  ::  cma_p3 (n_cma_p3) =   (/&
!     ! GPS:
!   t_cmatyp( ngps, 0     , 'GPS                  ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpgfz, 'by GFZ               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpasi, 'by ASI               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpbkg, 'by BKG               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpgop, 'by GOP               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpieec,'by IEEC              ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpsgn, 'by SGN               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngplpt, 'by LPT               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpmet, 'by MET               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngprob, 'by ROB               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpige, 'by IGE               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpknmi,'by KNMI              ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 ),&
!   t_cmatyp( ngps, ngpnga, 'by NGA               ',cs,cn,cw,ce,cn,cs,ce,cw,0,0,0,0,0 )/)

  TYPE (t_cmatyp), PARAMETER  , PRIVATE  :: &
    cma_p ( 0 : SIZE(cma_p1)+SIZE(cma_p2)+SIZE(cma_p3) )  =  (/ cma_p0, cma_p1, cma_p2, cma_p3 /)

  TYPE (t_cmatyp)                        :: &
    cma   ( 0 : SIZE(cma_p1)+SIZE(cma_p2)+SIZE(cma_p3) )  =  cma_p

  INTEGER        , PARAMETER  :: &
    n_cma  =  n_cma_p1 + n_cma_p2 + n_cma_p3 ,& ! number of CMA obs + code types
    n_gpc_offs  =  n_cma - n_gpc                ! index offset of GPS code type
                                                !   elements in array 'cma'
!-------------------------------------------------------------------------------
! Section 7:  Functions
!-------------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------------
FUNCTION i_cma ( kobtyp , kcdtyp )

!---------------------------------------------------
! function to determine the index of cma
! referring to a given CMA observation and code type
!---------------------------------------------------
  INTEGER        , INTENT (IN) ::  kobtyp  ,& ! CMA observation type
                                   kcdtyp     ! CMA code type
  INTEGER        :: i_cma , ii
!-------------------------------------------------------------------------------

  i_cma = 0
  DO ii = 1, SIZE( cma ) - 1
    IF ((cma(ii)%obtyp == kobtyp) .AND. (cma(ii)%cdtyp == kcdtyp))  i_cma = ii
  ENDDO

END FUNCTION i_cma
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------

END MODULE data_obs_lib_cosmo
