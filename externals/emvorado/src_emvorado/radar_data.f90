! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2005-2024, DWD, KIT
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_data

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares all namelist variables and fields needed for the
!  radar operators. It also contains functions for conversion of some
!  new derived TYPEs.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

#ifndef NOMPI
  USE mpi
#endif

  ! Kind parameters of radar variables:
  USE radar_kind, ONLY : dp, wp, sp, i8, wpfwo
  
  USE radar_dbzcalc_params_type, ONLY : t_dbzcalc_params
  
  IMPLICIT NONE

  !==============================================================================

  ! INCLUDE statements

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

  !==============================================================================

  PUBLIC

  !==============================================================================

  INTERFACE rsm_init_strings_blanks
    MODULE PROCEDURE              &
         rsm_init_strings_blanks_multime, &
         rsm_init_strings_blanks_onetime, &
         rsm_init_strings_blanks_multimv, &
         rsm_init_strings_blanks_onetimv
  END INTERFACE rsm_init_strings_blanks

  !==============================================================================

  ! Data types for radar specific MPI:
  INTEGER, PARAMETER ::   &
    imp_fwo_double   = MPI_DOUBLE_PRECISION, & ! determines the REAL type for MPI
    imp_fwo_integers = MPI_INTEGER             ! determines the INTEGER type for MPI

  !.. Data types for radar_meta- and dbzparams-structures for MPI-exchange:
  INTEGER            ::          &
    mpi_dbzcalc_params_typ,      &
    mpi_polmp_typ,               &
    mpi_radar_meta_typ_alltimes, &
    mpi_radar_meta_typ_onetime,  &
    mpi_comp_meta_typ,           &
    mpi_voldata_ostream_typ

  !.. Basic constants in wpfwo precision:
  REAL(KIND=wpfwo), PARAMETER :: pi     = 4.0_wpfwo * ATAN(1.0_wpfwo)
  REAL(KIND=wpfwo), PARAMETER :: degrad = pi / 180.0_wpfwo
  REAL(KIND=wpfwo), PARAMETER :: raddeg = 180.0_wpfwo / pi

  ! Length of character strings for filenames, directories and error messages:
  INTEGER, PARAMETER :: cmaxlen = 300

  ! Length of character strings for obs filenames in rs_meta_type:
  INTEGER, PARAMETER :: cobsflen = 300

  ! Length of character strings for varnames in radar_meta_type:
  INTEGER, PARAMETER :: cvarlen = 10

  ! Length of character strings for namelist printing:
  INTEGER, PARAMETER :: cnmlstlen = 2000

  ! max. length of elevation-vector and obstimes-vector and other vectors:
  INTEGER, PARAMETER :: nel_max = 25
  INTEGER, PARAMETER :: nel_composite_max = nel_max
  INTEGER, PARAMETER :: nscanstrategies_max = 10
  INTEGER, PARAMETER :: nobstimes_max = 500
  INTEGER, PARAMETER :: ngpsm_max = 15
  INTEGER, PARAMETER :: nradsta_max = 140
  INTEGER, PARAMETER :: noutput_fields_max = 50  ! length of list of output vars for volume data
  INTEGER, PARAMETER :: ndoms_max = 3  ! Maximum number of radar-active model domains supported by EMVORADO
  INTEGER, PARAMETER :: nautobubbles_max = 30  ! Maximum allowed number of automatically detected bubbles for the bubble generator
  INTEGER, PARAMETER :: noutstreams_max = 5
  
  ! Setup of internal identifiers used as indices into data arrays
  !  to address certain radar parameters:
  ! - Index 1 = vr
  ! - Index 2 = qv
  ! - Index 3 = z
  ! - Index 4 = qz
  ! - Index 5 = zdr
  ! ...
  INTEGER, PARAMETER              :: ndatakind  = 10
  CHARACTER(len=*), PARAMETER     :: cndatakind = '10' ! character form of ndatakind

  ! Internal identifiers for radar moments, have to be > 0 and <= ndatakind (increase ndatakind if needed):
  INTEGER, PARAMETER :: i_vrad     = 1    ! internal unique identifier for radial wind
  INTEGER, PARAMETER :: i_qualvrad = 2    !    - " -  for quality flag (DWD only) of vrad
  INTEGER, PARAMETER :: i_dbzh     = 3    !    - " -  for horizontal polar. reflectivity
  INTEGER, PARAMETER :: i_qualdbzh = 4    !    - " -  for quality flag (DWD only) of dbzh
  INTEGER, PARAMETER :: i_zdr      = 5    !    - " -  for differential reflectivity
  INTEGER, PARAMETER :: i_kdp      = 6    !    - " -  for spec. differential phase shift [deg/km from -15 to 15]
  INTEGER, PARAMETER :: i_phidp    = 7    !    - " -  for total differential propagation phase [deg from -180 to 180]
  INTEGER, PARAMETER :: i_rhv      = 8    !    - " -  for h/v correlation coefficient
  INTEGER, PARAMETER :: i_ldr      = 9    !    - " -  for linear depolarization ratio
  INTEGER, PARAMETER :: i_cflags   = 10   !    - " -  for quality flags of polarization data

  ! Max. number of countries:
  INTEGER, PARAMETER :: ncountry_max = 11

  ! Flags for different "countries" (= collection of radars of similar type, properties and scan strategy)
  !  (must be a number between 0 and ncountry_max-1 )
  INTEGER, PARAMETER :: i_fakeobs     = 0  ! special flag for faking obs data from simulated data for testsuite;
                                           ! will get properties of dwd radars
  INTEGER, PARAMETER :: i_dwd         = 1
  INTEGER, PARAMETER :: i_meteoswiss  = 2
  INTEGER, PARAMETER :: i_arpasim     = 3
  INTEGER, PARAMETER :: i_belgium     = 4
  INTEGER, PARAMETER :: i_denmark     = 5
  INTEGER, PARAMETER :: i_france      = 6
  INTEGER, PARAMETER :: i_poland      = 7
  INTEGER, PARAMETER :: i_czech       = 8
  INTEGER, PARAMETER :: i_netherlands = 9
  INTEGER, PARAMETER :: i_slovakia    = 10

  ! Number of radar stations of the German radar network (icountry = 1 = i_dwd):
  INTEGER, PARAMETER :: nradsta_dwd         = 64 ! (17 new + 5 old stations + 5 "fail safe radars" + 4 OSSEs radars
                                                 !  + precip scans for each of these 31 radar
                                                 !  + 2 KIT radars)
  INTEGER, PARAMETER :: nradsta_swiss       = 10 ! (5 stations, both with full 20 elev and reduced 5 elev)
  INTEGER, PARAMETER :: nradsta_italy       = 17 ! Only Emilia-Romagna radar network
  INTEGER, PARAMETER :: nradsta_belgium     = 2  ! Only Wideumont and Jabekke; others do not deliver all necessary moments on common elevations
  INTEGER, PARAMETER :: nradsta_denmark     = 5
  INTEGER, PARAMETER :: nradsta_france      = 20 ! Only the 20 C-Band radars; there are also 5 S-Bands
  INTEGER, PARAMETER :: nradsta_poland      = 8
  INTEGER, PARAMETER :: nradsta_czech       = 2
  INTEGER, PARAMETER :: nradsta_netherlands = 2
  INTEGER, PARAMETER :: nradsta_slovakia    = 3
  
  ! Number of all radar stations in the entire network:
  INTEGER, PARAMETER :: nradsta_all = nradsta_dwd + nradsta_italy + nradsta_swiss + &
       &                nradsta_belgium + nradsta_denmark + nradsta_france + nradsta_poland + nradsta_czech + &
       &                nradsta_netherlands + nradsta_slovakia

  ! Internal format names for different obsfile formats:
  CHARACTER(len=15), PARAMETER    :: c_format_cdfin           = 'cdfin          ' ! CDFIN format, single-moment multi-sweep, either single-time or multi-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_native_smss  = 'h5-native-smss ' ! native hdf5 single-moment single-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_native_smms  = 'h5-native-smms ' ! native hdf5 single-moment multi-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_native_mmss  = 'h5-native-mmss ' ! native hdf5 multi-moment single-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_native_mmms  = 'h5-native-mmms ' ! native hdf5 multi-moment multi-sweep single-time 
  CHARACTER(len=15), PARAMETER    :: c_format_h5_2_opera_smss = 'h5-2-opera-smss' ! native hdf5 single-moment single-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_2_opera_smms = 'h5-2-opera-smms' ! native hdf5 single-moment multi-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_2_opera_mmss = 'h5-2-opera-mmss' ! native hdf5 multi-moment single-sweep single-time
  CHARACTER(len=15), PARAMETER    :: c_format_h5_2_opera_mmms = 'h5-2-opera-mmms' ! native hdf5 multi-moment multi-sweep single-time
  
  ! Dummy filename for obs/fdbk files:
  CHARACTER(len=*), PARAMETER     :: obsfile_missingname = 'nofile_obs'
  CHARACTER(len=*), PARAMETER     :: fdbkfile_missingname = 'nofile_fdbk'
  
  ! FO-internal identifiers type of grid-scale precipitation physics:
  !  Numbers in parentheses at the end of the comment give the equivalent COSMO-like itype_gscp previously used
  INTEGER, PARAMETER :: fwo_gscp_allowed(7) = [120,&      ! 1-mom, incl. cr,     ie. frozen hydrometeors at all (1)
                                               130,&      ! 1-mom, incl. cri,    ie. no snow (and graupel)      (2)
                                               140,&      ! 1-mom, incl. cris,   ie. no graupel                 (3)
                                               150,&      ! 1-mom, incl. crisg,  ie. with graupel               (4) 
                                               151,&      ! 1-mom, incl. crisg,  ie. with graupel, and with modified snow microphysics (5-99)
                                               250,&      ! 2-mom, incl. crisg,  ie. no hail                    (100:1999)
                                               260&       ! 2-mom, incl. crisgh, ie. with hail                  (2000:)
                                              ]


  ! Recurring thresholds and markers for initializations, missing values, etc.:
  REAL(dp), PARAMETER :: Z_crit_radar    = 1d-9, &
                         ext_crit_radar  = 1d-12, &
                         zero_value      = -99.99_dp, &   ! correct 0's ("no signal") in log10 space (dBZ)
                         shield_value    = -988.88_dp, &  ! flag value indicating terrain shielded radar bins
                         shield_value_rhv= -HUGE(1.0_dp)*0.98_dp, &  ! flag value indicating terrain shielded radar bins for rhohv
                         miss_threshold  = -900.0_dp, &   ! threshold to distinguish between correct 0's and missing values
                         miss_value      = -999.99_dp, &  ! flag value for missing values ("no data")
                         missing_obs     = -999.99_dp, &  ! flag value for missing values ("no data") from obs files (could be replaced by miss_value)
                         miss_thresh_rhv = -HUGE(1.0_dp)*0.9_dp, &   ! threshold to distinguish between correct 0's and missing values for rhohv
                         miss_value_rhv  = -HUGE(1.0_dp), &  ! flag value for missing values ("no data") for rhohv
                         reject_value    = -955.55_dp, &  ! flag value for rejected observations (due to quality checks)
                         reject_threshold= -960.0_dp, &   ! threshold to distinguish rejected observations from missing obs
                         unused_value    = -999.9_dp, &   ! general "neutral" initialization value, could be replaced by miss_value, if some formats in namelist ctrl outputs are adjusted to 2 decimals.
                         eps_dtobs       = 1e-6_dp        ! for the dt_obs-triplet: if dt_obs(1) and/or dt_obs(2) differ by eps_dtobs from model init and/or end time, the series of obs_times computed from this triplet will be synced to time t=0.0 s.

  REAL(dp), PARAMETER :: dBZ_crit_radar       = 10.0d0*LOG10(Z_crit_radar), &       ! -90
                         logext_crit_radar    = LOG10(ext_crit_radar), &            ! -12
                         shield_low_threshold = FLOOR(shield_value*0.1d0)*10.0d0, & ! -990
                         shield_up_threshold  = shield_low_threshold + 10.0d0, &    ! -980
                         shield_low_thresh_rhv= shield_value_rhv*0.99_dp, &         ! < shield_value_rhv
                         shield_up_thresh_rhv = shield_value_rhv*0.97_dp            ! > shield_value_rhv, but < miss_thresh_rhv

  ! type conversions
  REAL(wp), PARAMETER :: missval_wp = REAL(miss_value, kind=wp), &
                         missthr_wp = REAL(miss_threshold, kind=wp)

  REAL(sp), PARAMETER :: missval_sp        = REAL(miss_value, kind=sp), &
                         missthr_sp        = REAL(miss_threshold, kind=sp), &
                         zeroval_sp        = REAL(zero_value, kind=sp), &
                         Z_crit_radar_sp   = REAL(Z_crit_radar, kind=sp), &
                         dBZ_crit_radar_sp = REAL(dBZ_crit_radar, kind=sp)

  INTEGER,  PARAMETER :: missval_int = INT(miss_value), &
                         missthr_int = INT(miss_threshold)


!======================================================================================
!
! Variables from the driving model that have to be also present in EMVORADO
!
!======================================================================================

  ! Some physical constants:
  REAL    (KIND=dp)   ::   &
       time_mod_sec        ! model forecast time in [s] since model start

  ! Some flags for timing the runtime of different code parts. These are provided by the atmospheric model
  !  and are used to identify where to add the respective computation time in the timing output:
  INTEGER            ::   &
       i_fwo_prep_compute,&! Timing flag for preparations of the forward operator
       i_fwo_bubbles,   &  ! Timing flag for the bubble generator
       i_fwo_composites,&  ! Timing flag for the composite generation
       i_fwo_ini,       &  ! Time for the initialization of the forward operator
       i_fwo_compgrid,  &  ! Time for computations on the model grid
       i_fwo_comm,      &  ! Time for MPI-communications
       i_fwo_ongeom,    &  ! Time for the ray tracing in online beam propagation
       i_fwo_comppolar, &  ! Time for interpolation of reflectivity and radial wind
                           !  from model grid points to the radar bins/auxiliary azi slice grid
       i_fwo_out,       &  ! Time for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                           !  reading obs data, producing feedback files)
       i_fwo_barrier       ! Time for barrier waiting in MPI-communications (measure for load imbalance)

!!$ UB: re-organization in the future: put to model-specific file radar_interface.f90 or introduce #ifdef COSMO, #ifdef ICON etc.
  ! Field dimensions of model fields used in the radar operator:
  INTEGER            :: &
       ie_fwo, je_fwo, ke_fwo

  INTEGER            :: &
       ! variables from the driving model
       num_compute_fwo, & ! number of compute PEs
       my_cart_id_fwo,  & ! rank of this PE (=subdomain) in the cartesian communicator
       icomm_cart_fwo,  & ! communicator for the virtual cartesian topology

       ! only used in radar-routines
       num_radar,       & ! number of radar PEs (num_compute + num_radario)
       num_radario,     & ! total number of radar-IO PEs
       radar_master,    & ! root-PE of the radar PE group                      (in the world comm., normally = 0)
       radario_master,  & ! root-PE of the total radario group for all domains (in the world comm.!!!)
       my_radar_id,     & ! rank of this PE in the radar communicator (cart+radario)
       my_radario_id,   & ! rank of this PE in the (asynchroneous) radario communicator (radario)
       icomm_radar,     & ! communicator for the group of radar-IO PEs (all domains) + compute PEs
       icomm_radario      ! communicator for the group of radar-IO PEs (all domains)

  INTEGER         :: &
       ndoms_radar,  &  ! number of model domains with active radar simulation, is determined in prep_domains_radar() below.
       idom             ! Internal storage in EMVORADO of the domain index idom_model of the hosting model for the actual EMVORADO call

  ! Dictionaries to relate a domain index list of the radar-active domains (1:idom_radar) to the original model domain indices (1:idom_model)
  INTEGER, ALLOCATABLE :: &
       list_domains_for_model(:), & ! list_domains_for_model(idom_radar) returns the original model domain index of the internal
                                    ! idom_radar'th EMVORADO domain. This list is created in prep_domains_radar() below.
       list_domains_for_radar(:)    ! list_domains_for_radar(idom_model) returns the internal EMVORADO domain index of the
                                    ! hosting model's idom_model'th domain. This list is created in prep_domains_radar() below.

  INTEGER, ALLOCATABLE :: &
       num_radar_dom(:),       & ! number of radar PEs (num_compute + num_radario_dom) per radar-active model domain
       num_radario_dom(:),     & ! number of radar-IO PEs per radar-active model domain
       radario_master_dom(:),  & ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       icomm_radar_dom(:),     & ! communicator for the group of radar-IO PEs of each domain + compute PEs
       icomm_radario_dom(:),   & ! communicator for the group of radar-IO PEs for each domain
       my_radar_id_dom(:),     & ! rank of this PE in the radar communicator (cart+radario_dom)
       my_radario_id_dom(:)      ! rank of this PE in the (asynchroneous) radario communicator icomm_radario_dom


! INTEGER   , ALLOCATABLE , TARGET   ::     &
!   isubpos(:,:)       ! positions of the subdomains in the total domain. Given
!                      ! are the i- and the j-indices of the lower left and the
!                      ! upper right grid point in the order
!                      !                  i_ll, j_ll, i_ur, j_ur.
!                      ! Only the interior of the domains are considered, not
!                      ! the boundary lines.
!                      ! ('target' attribute is required because a pointer will
!                      !  point to this in the assimilation (obs processing))

  LOGICAL      ::         &
       lcompute_pe_fwo,   & ! indicates whether this is a compute PE or not
       lradar_pe,         & ! indicates whether this is a radar PE for any domain or not (compute or radar-IO)
       lradario_pe          ! indicates whether this is a radar-IO PE for any domain or not

  LOGICAL, ALLOCATABLE :: &
       lradar_pe_dom(:),  & ! indicates whether this is a radar PE for a certain domain or not (compute or radar-IO)
       lradario_pe_dom(:)   ! indicates whether this is a radar-IO PE for a certain domain or not

  ! Some physical constants from the model converted to dp:
  REAL    (KIND=dp)  ::   &
       r_earth_dp         ! mean radius of the earth

  ! Some physical constants in model wp, directly overtaken from the model:
  REAL(KIND=wp), PARAMETER :: pi_model  = 4.0_wp * ATAN(1.0_wp)
  REAL    (KIND=wp)  ::   &
       t0_melt_model,     & ! melting temperature (K)
       rho_w_model,       & ! density of liquid water (kg/m^3)
       rho_ice_model,     & ! density of ice          (kg/m^3)
       K_w_model,         & ! dielectric constant for water
       K_ice_model          ! dielectric constant for ice
  INTEGER :: itype_gscp_model ! the original microphysics choice flag of the hosting model
  LOGICAL :: luse_muD_relation_rain_model ! the original mu-D-relation flag of the hosting model
  
  ! year, month, day, hour, minute and second of model start (YYYYMMDDhhmmss):
  CHARACTER (len=14)    :: ydate_ini_mod

!======================================================================================
!
! Some derived TYPES for data and metadata:
!
!======================================================================================

  !------------------------------------------------------------------------------

  ! NOTE 1: all arrays in "radar_meta_type" are fixed-sized, because
  !         only then the SEQUENCE statement assures a continuous
  !         block in memory and enables simple MPI distribution.
  !         In addition, the data have to be naturally aligned, i.e.,
  !         first declare the REAL, then the INTEGER, LOGICAL and CHARACTER
  !         TYPE components in that order.
  !
  ! NOTE 2: There is a corresponding type below "radar_meta_type_onetime" which is for the
  !         case of 1 observation time instead of "nobstimes_max".
  !         If you change something in "radar_meta_type", do the same
  !         in "radar_meta_type_onetime", but replacing "nobstimes_max" by "1".
  !
  ! NOTE 3: Add new elements also to the below copy-function "rsm_multitime2onetime()"
  !
  ! NOTE 4: One is tempted to make use of the new class inheritance
  !         mechanism in F2008 to extend "radar_meta_type" and "radar_meta_type_onetime"
  !         from a base TYPE. However, this won't work because of the SEQUENCE
  !         keyword!

  TYPE radar_meta_type
    SEQUENCE   ! Important: ensures that the following parameters are one block in memory
               !       and that an MPI-distribution by a derived MPI-datatype is possible:
               !       - parallel_utilities.f90, def_mpi_radar_meta_type()
               !       - CALL def_mpi_radar_meta_type (mpi_radar_meta_typ)
               !         MPI_BCAST( <rsmeta>, 1, ..., mpi_radar_meta_typ, ...)
               !       For this, the following data have to be naturally aligned, i.e.,
               !       first the REAL, then the INTEGER, LOGICAL and CHARACTER at the end.
    REAL    (KIND=wpfwo)      :: lambda         ! Station radar wavelength (m)
    REAL    (KIND=wpfwo)      :: lat            ! Station latitude  (deg)
    REAL    (KIND=wpfwo)      :: lon            ! Station longitude (deg)
    REAL    (KIND=wpfwo)      :: alt_agl_mod    ! Station height AGL model orogr.   (m)
    REAL    (KIND=wpfwo)      :: alt_msl        ! Station height MSL used in the model, derived from namelist input   (m)
    REAL    (KIND=wpfwo)      :: alt_msl_true   ! True station height MSL, taken from obs files  (m)
    REAL    (KIND=wpfwo)      :: msl_mod        ! Model oro height MSL at station   (m)
    REAL    (KIND=wpfwo)      :: az_inc         ! azimuth increment  (deg)
    REAL    (KIND=wpfwo)      :: az_start       ! start azimut of first az bin rel. to true north  (deg)
    REAL    (KIND=wpfwo)      :: ra_inc         ! range increment that is actually used for simulations (m)
    REAL    (KIND=wpfwo)      :: ra_inc_obs     ! range increment of observations which can be smaller than ra_inc, in which case aggretation to ra_inc is done on input (m)
    REAL    (KIND=wpfwo)      :: ra_inc_coarse  ! approximate range increment for range aggretation [m]. ra_inc will be set to the nearest value which is a multiple of ra_inc_obs
    REAL    (KIND=wpfwo)      :: Theta3         ! vertical 3-dB-oneway beam width in the beam coord. system  [degrees]
                                                !   (half width of the one-way beamfunction)
    REAL    (KIND=wpfwo)      :: Phi3           ! horizontal 3-dB-oneway beam width in the beam coord. system [degrees]
    REAL    (KIND=wpfwo)      :: dalpha         ! azimutal averaging interval for calculation of one
                                             !   averaged pulse (averaging over the number
                                             !   of statistically independent pulses)
    REAL    (KIND=wpfwo)      :: alpha3_eff_0   ! effective horizontal 3-dB-oneway beam width at
                                             !   elevation = 0.0, depending on the ratio (dalpha/phi3)
                                             !   --> will be automatically determined from a lookup table
    REAL    (KIND=wpfwo)      :: smth_interv_fact ! Factor to determine the azimutal and elevational integration range
                                             !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                             !   factor to the effective 3-dB-oneway beamwidths. See Blahak, JAOTECH, 2008.
                                             !   (a value of 1.29 leads to the 90-%-weight-range of the beam function)
    REAL    (KIND=wpfwo)      :: dt_obs(3)             ! triplet for time increment of observations in seconds: notation either <incr>,<miss>,<miss> or <from>,<to>,<incr>
    REAL    (KIND=wpfwo)      :: dt_obs_fdbk(3)        ! triplet for time increment of observations to write to the feedback file: notation either <incr>,<miss>,<miss> or <from>,<to>,<incr>
    REAL    (KIND=wpfwo)      :: dt_obs_voldata(3)     ! triplet for time increment of observations to write to the volume data files: notation either <incr>,<miss>,<miss> or <from>,<to>,<incr>
    REAL    (KIND=wpfwo)      :: ext_nyq(nel_max,nobstimes_max)  = -999.9_wpfwo      ! Extended Nyquist velocity [m/s]
    REAL    (KIND=wpfwo)      :: high_nyq(nel_max) = -999.9_wpfwo      ! high Nyquist velocity [m/s]
    REAL    (KIND=wpfwo)      :: prf(nel_max)      = -999.9_wpfwo      ! Low PRF        [Hz]
    REAL    (KIND=wpfwo)      :: dualprf_ratio(nel_max) = -999.9_wpfwo ! dual PRF ratio [-]
    REAL    (KIND=wpfwo)      :: rngate_len = -999.9_wpfwo ! range gate length [m]
    REAL    (KIND=wpfwo)      :: mds_Z0             ! Minimum detectable signal at a reference range [dBZ]
    REAL    (KIND=wpfwo)      :: mds_r0             ! Reference range for minimum detectable signal [m]
    REAL    (KIND=wpfwo)      :: obs_times(nobstimes_max) ! array of observation times (s) as difference to model start in seconds
    REAL    (KIND=wpfwo)      :: obs_times_obs(nobstimes_max) ! array of observation times (s) coming from the obs files
    REAL    (KIND=wpfwo)      :: obs_times_fdbk(nobstimes_max) ! array of observation times (s) to write to the feedback file
    REAL    (KIND=wpfwo)      :: obs_times_voldata(nobstimes_max) ! array of observation times (s) to write to the volume data files
    REAL    (KIND=wpfwo)      :: el_arr(nel_max)         ! array of nominal elevations (deg) (these are used for all computations)
    REAL    (KIND=wpfwo)      :: el_arr_default(nel_max,nscanstrategies_max) ! array of nominal elevations (deg) (these are used for all computations)
    REAL    (KIND=wpfwo)      :: xabscsm_v(ngpsm_max)  ! array of abscissa values for vertical Gauss-Legendre quadrature [-]
    REAL    (KIND=wpfwo)      :: xabscsm_h(ngpsm_max)  ! array of abscissa values for horizonal Gauss-Legendre quadrature [-]
    REAL    (KIND=wpfwo)      :: weigsm_v(ngpsm_max)   ! array of weights for vertical Gauss-Legendre quadrature [-]
    REAL    (KIND=wpfwo)      :: weigsm_h(ngpsm_max)   ! array of weights for horiztontal Gauss-Legendre quadrature [-]
    REAL    (KIND=wpfwo)      :: vnyq_min_for_vr_active_fdbk ! Filter for writing radial winds to fdbk-files: only data with larger Nyquist-Veloc. are written
    INTEGER                   :: icountry           ! Internal Country identifier
    INTEGER                   :: station_id         ! WMO Station ident code
    INTEGER                   :: ista               ! internal station index  (only for internal bookkeeping of cdfin-output; do not set in namelist!)
    INTEGER                   :: nel                ! actual number of elevations in the nominal scan strategy
    INTEGER                   :: nel_present        ! actual number of elevations which are present at the actual time (<= nel)
    INTEGER                   :: nel_fdbk           ! actual number of elevations written into feedback file
    INTEGER                   :: nel_voldata        ! actual number of elevations written into volume data files
    INTEGER                   :: nel_default(nscanstrategies_max) ! number of elevations for each possible default scan strategy
    INTEGER                   :: eleind_for_composite_bub  ! elevation index for construction of the composite
                                                           !  in the warm bubble generator (comp_dbzsim_bub, comp_dbzobs_bub)
    INTEGER                   :: naz                ! number of nominal azimuths
    INTEGER                   :: naz_ncdf(nobstimes_max,ndatakind)! number of azimuths in CDFIN NetCDF files (can be > naz)
    INTEGER                   :: nra                ! max. number of range bins occuring in a volume scan
    INTEGER                   :: nra_obs            !   "-" from observation data
    INTEGER                   :: n_aggr_ra_obs      ! number of range bins to aggregate/average when reading obs data ( = nint(ra_inc/ra_inc_obs))
    INTEGER                   :: nobs_times         ! number of observation times
    INTEGER                   :: nobs_times_obs     ! number of observation times in obs files
    INTEGER                   :: nobs_times_fdbk    ! number of observation times in feedback files
    INTEGER                   :: nobs_times_voldata ! number of observation times in volume data files
    INTEGER                   :: i_nearest_mod      ! i-Index of nearest model grid point to radar station
    INTEGER                   :: j_nearest_mod      ! j-Index of nearest model grid point to radar station
    INTEGER                   :: ngpsm_v            ! number of vertical smoothing points for Gauss-Legendre-quadrature
    INTEGER                   :: ngpsm_h            ! number of horizontal smoothing points for Gauss-Legendre-quadrature
    INTEGER                   :: nrep_ncdf(ndatakind) = missval_int  ! actual number of reports in fdkbfile; for bookkeeping
    INTEGER                   :: num_gates = missval_int   ! number of gates averaged
    INTEGER                   :: num_pulses= missval_int   ! number of pulses integrated
    INTEGER                   :: obs_startrec(nobstimes_max,ndatakind) ! array of observation start records, second index for file 1=vr, 2=qv, 3=z, 4=qz
    INTEGER                   :: obs_endrec(nobstimes_max,ndatakind)   ! array of observation end records, second index for file 1=vr, 2=qv, 3=z, 4=qz
    INTEGER                   :: ind_ele_present(nel_max)  ! List of the indices of nominal elevations which are present for the actual time; will be determined only when reading the data for the actual timestep
    INTEGER                   :: ind_ele_fdbk(nel_max) ! array of indices which elevations of the list of present elevations are written into feedback file
    INTEGER                   :: ind_ele_voldata(nel_max) ! array of indices which elevations of the list of present elevations are written into volume data files
    INTEGER                   :: eleindlist_for_composite(nel_composite_max)  ! elevation index list of the list of present elevations for construction of composite(s) (comp_dbzsim, comp_dbzobs)
    LOGICAL                   :: lobs_avail(nobstimes_max,ndatakind)     ! observations available
    LOGICAL                   :: lfdbkfile_exist    ! feedback file exists
    LOGICAL                   :: lobstimes_ovwrt_recalc  ! whether or not it should be enabled that, for station meta data re-definition after
                                                         ! metadata reading, the parameters rs_meta%dt_obs and rs_meta%nobs_times are used
                                                         ! to compute the series of desired radar output times. By default (.false.), only
                                                         ! rs_meta%obs_times takes effect when overwritten. If .true., it has the side effect that
                                                         ! it is possible to add further stations in the re-definition step that are not present
                                                         ! in the observations. The added stations however must be defined in the internal background meta data list.
    LOGICAL                   :: lvrad_to_fdbk      ! write radial wind of this stationt to feedback files
    LOGICAL                   :: ldbzh_to_fdbk      ! write horizontal reflectivity of this station to feedback files
    CHARACTER (LEN=3)         :: station_name       ! short station name
    CHARACTER (LEN=10)        :: scanname           ! short name to denote the scan strategy: e.g. PPI1230, RHI3456, PREC
    CHARACTER (LEN=14)        :: obs_cdate(nobstimes_max)   ! date and time of observation
    CHARACTER (LEN=cobsflen)  :: obsfile(nobstimes_max,ndatakind)      ! name of input files, for vr, qv, z, qz, ...
    CHARACTER (len=15)        :: obsfile_format(nobstimes_max)
    CHARACTER (LEN=60)        :: fdbkfile         ! name of NetCDF feedback file for vr, qv, z and qz
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_vrad ! HDF5 shortnames of the desired radar moments from hdf5 obs files
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_dbzh !  Have to be set in the background meta data lists in radar_obs_meta_list.f90
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_zdr  !  Changing/Overwriting them by namelist has no effect!
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_rhv  ! These names in the type simplify the handling of different names for different icountries.
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_kdp
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_phidp
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_ldr
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_cflags
    CHARACTER (len=cobsflen)  :: dummy_buf_for_memalign  ! dummy buffer to compensate for possible memory alignment in mpi transfers
  END TYPE radar_meta_type

  TYPE radar_meta_type_onetime
    SEQUENCE
    REAL    (KIND=wpfwo)      :: lambda
    REAL    (KIND=wpfwo)      :: lat
    REAL    (KIND=wpfwo)      :: lon
    REAL    (KIND=wpfwo)      :: alt_agl_mod
    REAL    (KIND=wpfwo)      :: alt_msl
    REAL    (KIND=wpfwo)      :: alt_msl_true
    REAL    (KIND=wpfwo)      :: msl_mod
    REAL    (KIND=wpfwo)      :: az_inc
    REAL    (KIND=wpfwo)      :: az_start
    REAL    (KIND=wpfwo)      :: ra_inc
    REAL    (KIND=wpfwo)      :: ra_inc_obs
    REAL    (KIND=wpfwo)      :: ra_inc_coarse
    REAL    (KIND=wpfwo)      :: Theta3
    REAL    (KIND=wpfwo)      :: Phi3
    REAL    (KIND=wpfwo)      :: dalpha
    REAL    (KIND=wpfwo)      :: alpha3_eff_0
    REAL    (KIND=wpfwo)      :: smth_interv_fact
    REAL    (KIND=wpfwo)      :: dt_obs(3)
    REAL    (KIND=wpfwo)      :: dt_obs_fdbk(3)
    REAL    (KIND=wpfwo)      :: dt_obs_voldata(3)
    REAL    (KIND=wpfwo)      :: ext_nyq(nel_max,1)     = -999.9_wpfwo
    REAL    (KIND=wpfwo)      :: high_nyq(nel_max)      = -999.9_wpfwo
    REAL    (KIND=wpfwo)      :: prf(nel_max)           = -999.9_wpfwo
    REAL    (KIND=wpfwo)      :: dualprf_ratio(nel_max) = -999.9_wpfwo
    REAL    (KIND=wpfwo)      :: rngate_len             = -999.9_wpfwo
    REAL    (KIND=wpfwo)      :: mds_Z0
    REAL    (KIND=wpfwo)      :: mds_r0
    REAL    (KIND=wpfwo)      :: obs_times(1)
    REAL    (KIND=wpfwo)      :: obs_times_obs(1)
    REAL    (KIND=wpfwo)      :: obs_times_fdbk(1)
    REAL    (KIND=wpfwo)      :: obs_times_voldata(1)
    REAL    (KIND=wpfwo)      :: el_arr(nel_max)
    REAL    (KIND=wpfwo)      :: el_arr_default(nel_max,nscanstrategies_max)
    REAL    (KIND=wpfwo)      :: xabscsm_v(ngpsm_max)
    REAL    (KIND=wpfwo)      :: xabscsm_h(ngpsm_max)
    REAL    (KIND=wpfwo)      :: weigsm_v(ngpsm_max)
    REAL    (KIND=wpfwo)      :: weigsm_h(ngpsm_max)
    REAL    (KIND=wpfwo)      :: vnyq_min_for_vr_active_fdbk
    INTEGER                   :: icountry
    INTEGER                   :: station_id
    INTEGER                   :: ista
    INTEGER                   :: nel
    INTEGER                   :: nel_present
    INTEGER                   :: nel_fdbk
    INTEGER                   :: nel_voldata
    INTEGER                   :: nel_default(nscanstrategies_max)
    INTEGER                   :: eleind_for_composite_bub
    INTEGER                   :: naz
    INTEGER                   :: naz_ncdf(1,ndatakind)
    INTEGER                   :: nra
    INTEGER                   :: nra_obs
    INTEGER                   :: n_aggr_ra_obs
    INTEGER                   :: nobs_times
    INTEGER                   :: nobs_times_obs
    INTEGER                   :: nobs_times_fdbk
    INTEGER                   :: nobs_times_voldata
    INTEGER                   :: i_nearest_mod
    INTEGER                   :: j_nearest_mod
    INTEGER                   :: ngpsm_v
    INTEGER                   :: ngpsm_h
    INTEGER                   :: nrep_ncdf(ndatakind) = missval_int
    INTEGER                   :: num_gates = missval_int
    INTEGER                   :: num_pulses= missval_int
    INTEGER                   :: obs_startrec(1,ndatakind)
    INTEGER                   :: obs_endrec(1,ndatakind)
    INTEGER                   :: ind_ele_present(nel_max)
    INTEGER                   :: ind_ele_fdbk(nel_max)
    INTEGER                   :: ind_ele_voldata(nel_max)
    INTEGER                   :: eleindlist_for_composite(nel_composite_max)
    LOGICAL                   :: lobs_avail(1,ndatakind)
    LOGICAL                   :: lfdbkfile_exist
    LOGICAL                   :: lobstimes_ovwrt_recalc
    LOGICAL                   :: lvrad_to_fdbk
    LOGICAL                   :: ldbzh_to_fdbk
    CHARACTER (LEN=3)         :: station_name
    CHARACTER (LEN=10)        :: scanname
    CHARACTER (LEN=14)        :: obs_cdate(1)
    CHARACTER (LEN=cobsflen)  :: obsfile(1,ndatakind)
    CHARACTER (LEN=15)        :: obsfile_format(1)
    CHARACTER (LEN=60)        :: fdbkfile
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_vrad
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_dbzh
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_zdr
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_rhv
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_kdp
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_phidp
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_ldr
    CHARACTER (len=cvarlen)   :: obs_hdf5_varname_cflags
    CHARACTER (len=cobsflen)  :: dummy_buf_for_memalign  ! dummy buffer to compensate for possible memory alignment in mpi transfers
  END TYPE radar_meta_type_onetime


  !------------------------------------------------------------------------------

  TYPE radar_data_type
    REAL    (KIND=wpfwo), POINTER :: w_intp(:,:) => NULL()   ! array of horizontal interpolation
                                                             ! weights for each observation (first dimension)
                                                             ! and each spatial direction i,j,k (second dimension)
    REAL    (KIND=wpfwo), POINTER :: w_intp_smth(:,:) => NULL()
    INTEGER             , POINTER :: ind_intp(:,:) => NULL() ! grid indices for each observation (first dimension)
                                                             ! the second dimension contains of
                                                             ! 1:     the continuous number of the model grid
                                                             !        cell associated with the observation
                                                             ! 2,3,4: the observation indices in azimuthal,
                                                             !        radial and elevation direction (m,n,o)
                                                             ! or 2 : the continuous observation index combined
                                                             !        from radial, azimutal and elevation index
    INTEGER             , POINTER :: ind_intp_smth(:,:) => NULL()
    INTEGER             , POINTER :: radpos_all(:) => NULL() ! used for station data output: continuous index constructed from ra, az, el
    INTEGER             , POINTER :: radpos_all_smth(:) => NULL() ! used for station data output: continuous index constructed from ra, az, el
    REAL    (KIND=wpfwo), POINTER :: vt_mod(:) => NULL()    ! model aequivalents of fall speed on local PE
    REAL    (KIND=wpfwo), POINTER :: vt_mod_smth(:) => NULL()     ! model aequivalents of radial winds on local PE
    REAL    (KIND=wpfwo), POINTER :: radwind_mod(:) => NULL()     ! model aequivalents of radial winds on local PE
    REAL    (KIND=wpfwo), POINTER :: radwind_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: zh_radar_mod(:) => NULL()    ! model aequivalents of horizontal radar reflectivity on local PE
    REAL    (KIND=wpfwo), POINTER :: zh_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: ah_radar_mod(:) => NULL()    ! model aequivalents of horizontal extinction on local PE
    REAL    (KIND=wpfwo), POINTER :: ah_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: zv_radar_mod(:) => NULL()    ! model aequivalents of vertical radar reflectivity on local PE
    REAL    (KIND=wpfwo), POINTER :: zv_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: rrhv_radar_mod(:) => NULL()  ! model aequivalents of real correlation coeff. on local PE
    REAL    (KIND=wpfwo), POINTER :: rrhv_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: irhv_radar_mod(:) => NULL()  ! model aequivalents of imag correlation coeff. on local PE
    REAL    (KIND=wpfwo), POINTER :: irhv_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: kdp_radar_mod(:) => NULL()   ! model aequivalents of spec. diff. phase on local PE
    REAL    (KIND=wpfwo), POINTER :: kdp_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: adp_radar_mod(:) => NULL()   ! model aequivalents of spec. diff. attenuation on local PE
    REAL    (KIND=wpfwo), POINTER :: adp_radar_mod_smth(:) => NULL()
    REAL    (KIND=wpfwo), POINTER :: zvh_radar_mod(:) => NULL()   ! model aequivalents of hor-turned-vert radar reflectivity (for LDR) on local PE
    REAL    (KIND=wpfwo), POINTER :: zvh_radar_mod_smth(:) => NULL()

    REAL    (KIND=wpfwo), POINTER :: hl_loc(:) => NULL()  ! beam height above MSL
    REAL    (KIND=wpfwo), POINTER :: el_loc(:) => NULL()  ! local elevation
    REAL    (KIND=wpfwo), POINTER :: s_loc (:) => NULL()  ! arc distance to radar station at MSL

    INTEGER             , POINTER :: n_sh(:,:,:,:) => NULL()  ! first range index which is blocked by orography for output_radar_smth()

    INTEGER                       :: nobs               ! number of simulated observations on local PE
    INTEGER                       :: nobs_above_sfc     ! number of simulated observations on local PE above SFC for 4/3 earth model
    INTEGER                       :: nsmth              ! like nobs, but for smth version
    INTEGER                       :: nsmth_above_sfc    ! like nobs_above_sfc, but for smth version with 4/3 earth model

    ! Data vectors for observations:
    INTEGER                       :: nobs_obs(ndatakind)  ! number of observations read from NetCDF on local PE
    INTEGER                       :: check(4)
    REAL    (KIND=wpfwo), POINTER :: radwind_obs(:) => NULL()      ! observed radial winds on local PE
    REAL    (KIND=wpfwo), POINTER :: radwind_obs_q(:) => NULL()    ! quality index of observed radial winds
    REAL    (KIND=wpfwo), POINTER :: zh_radar_obs(:) => NULL()     ! observed reflectivity
    REAL    (KIND=wpfwo), POINTER :: zh_radar_obs_q(:) => NULL()   ! quality index of observed reflectivity
    REAL    (KIND=wpfwo), POINTER :: zdr_radar_obs(:) => NULL()    ! observed differential reflectivity
    REAL    (KIND=wpfwo), POINTER :: kdp_radar_obs(:) => NULL()    ! observed spec. differential phase shift
    REAL    (KIND=wpfwo), POINTER :: phidp_radar_obs(:) => NULL()  ! observed differential propagation phase
    REAL    (KIND=wpfwo), POINTER :: rhv_radar_obs(:) => NULL()    ! observed H/V correllation coefficient

    INTEGER             , POINTER :: ind_intp_obs(:) => NULL()  ! continuous index constructed from ra, az, el
    INTEGER             , POINTER :: ind_intp_obs_vr(:) => NULL()
    INTEGER             , POINTER :: ind_intp_obs_qv(:) => NULL()
    INTEGER             , POINTER :: ind_intp_obs_z(:) => NULL()
    INTEGER             , POINTER :: ind_intp_obs_qz(:) => NULL()

  END TYPE radar_data_type

  !===============================================================================

  TYPE cart_data_type

     REAL    (KIND=wpfwo)          :: pollon
     REAL    (KIND=wpfwo)          :: pollat
     REAL    (KIND=wpfwo)          :: polgam
     REAL    (KIND=wpfwo)          :: startlon
     REAL    (KIND=wpfwo)          :: startlat
     REAL    (KIND=wpfwo)          :: endlon
     REAL    (KIND=wpfwo)          :: endlat
     REAL    (KIND=wpfwo)          :: dlat ! latitude resolution of cartesian coordinates [deg]
     REAL    (KIND=wpfwo)          :: dlon ! longitude resolution of cartesian coordinates [deg]
     REAL    (KIND=wpfwo)          :: resolution !  resolution of cartesian coordinates [m]
     REAL    (KIND=wpfwo)          :: width     ! half total averaging width [m]
     REAL    (KIND=wpfwo)          :: width_deg ! half total averaging width [deg]
     REAL    (KIND=wpfwo)          :: minra_vr ! minimal range to radar station to do averaging for radial wind [m]
     REAL    (KIND=wpfwo)          :: minra_z  ! minimal range to radar station to do averaging for reflectivity [m]
     REAL    (KIND=wpfwo), POINTER :: dis(:,:) => NULL() ! distance to the closest radar point
     REAL    (KIND=wpfwo), POINTER :: calt(:,:) => NULL()! altitude of cartesian grid point
     INTEGER                    :: ni_tot  ! total number in i-direction
     INTEGER                    :: nj_tot  ! total number in j-direction
     INTEGER                    :: ni  ! number in i-direction within a rectangle surrounding the radar scan radius
     INTEGER                    :: nj  ! number in j-direction within a rectangle surrounding the radar scan radius
     INTEGER                    :: aw_max ! aw should not be bigger than aw_max [degrees]
     INTEGER                    :: naz_average
     INTEGER                    :: nra_average
     INTEGER          , POINTER :: rdata_ind(:,:) => NULL()   ! indices of azimuth, elevation and range of the closest radar point
     INTEGER          , POINTER :: cartdata_ind(:,:) => NULL() ! combined i/j index of cartesian point
     INTEGER          , POINTER :: aw(:,:)  => NULL() ! half averaging azimuthal width

  END TYPE cart_data_type

  !------------------------------------------------------------------------------

  TYPE radar_grid_type

    REAL    (KIND=wpfwo)   , POINTER :: w_intp(:,:) => NULL()        ! array of horizontal interpolation
                                                            ! weights for each observation (first dimension)
                                                            ! and each spatial direction i,j,k (second dimension)
    INTEGER                , POINTER :: ind_intp(:,:) => NULL()      ! grid indices for each observation (first dimension)
                                                            ! the second dimension contains of
                                                            ! 1:     the continuous number of the model grid
                                                            !        cell associated with the observation
                                                            ! 2,3,4: the observation indices in azimuthal,
    INTEGER                , POINTER :: ind_azgrd(:) => NULL()
    REAL    (KIND=wpfwo)   , POINTER :: hl_grd(:) => NULL()
    REAL    (KIND=wpfwo)   , POINTER :: hl_azgrd(:) => NULL()
    REAL    (KIND=wpfwo)             :: al_inc             ! arc length increment
    INTEGER                          :: naz_nbl            ! number of azimuths including nbl_az extra points (for lsmooth=.true.)
    INTEGER                          :: nal                ! number of arc length
    REAL    (KIND=wpfwo)   , POINTER :: vt_grd(:) => NULL()          ! model aequivalents of fall speed on local PE
    REAL    (KIND=wpfwo)   , POINTER :: u_grd(:) => NULL()           ! model aequivalents of u on local PE
    REAL    (KIND=wpfwo)   , POINTER :: v_grd(:) => NULL()           ! model aequivalents of v on local PE
    REAL    (KIND=wpfwo)   , POINTER :: w_grd(:) => NULL()           ! model aequivalents of w on local PE
    REAL    (KIND=wpfwo)   , POINTER :: zh_radar_grd(:) => NULL()    ! model aequivalents of horizontal radar reflectivity on local PE
    REAL    (KIND=wpfwo)   , POINTER :: ah_radar_grd(:) => NULL()    ! model aequivalents of horizontal extinction on local PE
    REAL    (KIND=wpfwo)   , POINTER :: zv_radar_grd(:) => NULL()    ! model aequivalents of vertical radar reflectivity on local PE
    REAL    (KIND=wpfwo)   , POINTER :: rrhv_radar_grd(:) => NULL()  ! model aequivalents of real RhoHV on local PE
    REAL    (KIND=wpfwo)   , POINTER :: irhv_radar_grd(:) => NULL()  ! model aequivalents of imag RhoHV on local PE
    REAL    (KIND=wpfwo)   , POINTER :: kdp_radar_grd(:) => NULL()   ! model aequivalents of KDP on local PE
    REAL    (KIND=wpfwo)   , POINTER :: adp_radar_grd(:) => NULL()   ! model aequivalents of ADP on local PE
    REAL    (KIND=wpfwo)   , POINTER :: zvh_radar_grd(:) => NULL()   ! model aequivalents of hor-turned-vert radar reflectivity (for LDR) on local PE
    REAL    (KIND=wpfwo)   , POINTER :: rfridx_grd(:) => NULL()
    REAL    (KIND=wpfwo)             :: rfridx_sta         ! refractive index at radar station
    REAL    (KIND=wpfwo)             :: rfridxgrad_sta     ! gradient of refractive index at radar station
    INTEGER                          :: ngrd               ! total number of radar bins in vectors *_grd

    REAL    (KIND=wpfwo)   , POINTER :: vt_azgrd(:) => NULL()          ! model aequivalents of fall speed on local PE after collection of entire azimut slices
    REAL    (KIND=wpfwo)   , POINTER :: u_azgrd(:) => NULL()           ! model aequivalents of u on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: v_azgrd(:) => NULL()           ! model aequivalents of v on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: w_azgrd(:) => NULL()           ! model aequivalents of w on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: zh_radar_azgrd(:) => NULL()    ! model aequivalents of horizontal radar reflectivity on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: ah_radar_azgrd(:) => NULL()    ! model aequivalents of horizontal extinction on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: zv_radar_azgrd(:) => NULL()    ! model aequivalents of vertical radar reflectivity on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: rrhv_radar_azgrd(:) => NULL()  ! model aequivalents of real RhoHV on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: irhv_radar_azgrd(:) => NULL()  ! model aequivalents of imag RhoHV on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: kdp_radar_azgrd(:) => NULL()   ! model aequivalents of KDP on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: adp_radar_azgrd(:) => NULL()   ! model aequivalents of ADP on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: zvh_radar_azgrd(:) => NULL()   ! model aequivalents of vertical radar reflectivity on local PE      - " -
    REAL    (KIND=wpfwo)   , POINTER :: rfridx_azgrd(:) => NULL()
    INTEGER                          :: nazgrd              ! total number of radar bins in vectors *_azgrd

  END TYPE radar_grid_type

  !------------------------------------------------------------------------------

  ! Type to hold a pointer to a vector of reals. Can be used to construct "vectors" of pointers:
  TYPE rpvect
    REAL    (KIND=wpfwo)   , POINTER :: p(:)  => NULL()
    INTEGER                       :: n
  END TYPE rpvect

  !------------------------------------------------------------------------------

  ! Type to hold a pointer to a vector of integers. Can be used to construct "vectors" of pointers:
  TYPE ipvect
    INTEGER              , POINTER :: p(:)  => NULL()
    INTEGER                        :: n
  END TYPE ipvect

  !------------------------------------------------------------------------------

  ! .. Type for meta data of feedback files:
  TYPE fdbk_meta_type
    INTEGER                   :: &
         ie_tot      ,& ! number of grid points in zonal direction
         je_tot      ,& ! number of grid points in meridional direction
         ke_tot         ! number of grid points in vertical direction
    REAL    (KIND=wpfwo)  :: &
         pollat      ,& ! latitude  of the rotated north pole (in degrees, N>0)
         pollon      ,& ! longitude of the rotated north pole (in degrees, E>0)
         polgam      ,& ! angle between the north poles of the systems
         dlat        ,& ! grid point distance in meridional direction (in degrees)
         dlon        ,& ! grid point distance in zonal      direction (in degrees)
         startlat_tot,& ! rotated latitude   \  of the lower left grid point of the
         startlon_tot   ! rotated longitude  /  total domain (in degrees, N,E>0)
    INTEGER                   :: &
         iveri_ens_member  ,& ! ensemble member ( -1: deterministic)
         nvers                ! exp. id + 16384*(class of model run)
    CHARACTER (LEN=cmaxlen)  :: &
                                ! NetCDF global attributes:
         yglatt_institution,& ! originating center name
         yglatt_source        ! program name and version
    REAL    (KIND=wpfwo)  :: &
         hversta     ,& ! start of verification period in 'model integration hours'
         hverend        ! end of verification period in 'model integration hours'
  END TYPE fdbk_meta_type

  !------------------------------------------------------------------------------

  ! .. Type for meta data of rotated lat-lon-grid for radar composites:
  TYPE composite_meta_type
    SEQUENCE   ! Important: ensures that the following parameters are one block in memory
               !       and that an MPI-distribution by a derived MPI-datatype is possible:
               !       - parallel_utilities.f90, def_mpi_compmeta_type()
               !       - CALL def_mpi_compmeta_type (mpi_compmeta_typ)
               !         MPI_BCAST( <comp_meta>, 1, ..., mpi_compmeta_typ, ...)
               !       For this, the following data have to be naturally aligned, i.e.,
               !       first the REAL, then the INTEGER, LOGICAL and CHARACTER at the end.
    REAL    (KIND=wpfwo)  :: &
         r_earth     ,& ! Earth radius for computing rotated lat/lon coordinates
         pollat      ,& ! latitude  of the rotated north pole (in degrees, N>0)
         pollon      ,& ! longitude of the rotated north pole (in degrees, E>0)
         dlat        ,& ! grid point distance in meridional direction (in degrees)
         dlon        ,& ! grid point distance in zonal      direction (in degrees)
         polgam      ,& ! angle between the north poles of the systems
         startlat    ,& ! rotated latitude   \  of the lower left grid point of the
         startlon       ! rotated longitude  /  total domain (in degrees, N,E>0)
    INTEGER                   :: &
         ni          ,& ! total number of grid points in zonal direction
         nj             ! total number of grid points in meridional direction
  END TYPE composite_meta_type

  !------------------------------------------------------------------------------

  ! .. Type to hold buffers for bubble positions and times detected by the
  !     bubble generator for each domain with "neutral" initialisation:

  TYPE bubble_list_type
    INTEGER :: nbubbles = 0, nbubbles_reject = 0
    REAL (kind=wpfwo), DIMENSION(nautobubbles_max) :: &
         bub_centlon   = -HUGE(1.0_dp), &
         bub_centlat   = -HUGE(1.0_dp), &
         bub_timestamp = -HUGE(1.0_dp)
  END TYPE bubble_list_type

!======================================================================================
!
! Instances of the above TYPES:
!
!======================================================================================

  ! array of radar station meta data (local PE)
  TYPE(radar_meta_type), ALLOCATABLE, TARGET, DIMENSION(:,:) :: rs_meta_container
  TYPE(radar_meta_type), POINTER, DIMENSION(:) :: rs_meta => NULL()

  ! array of radar data (local PE)
  TYPE(radar_data_type), ALLOCATABLE, TARGET, DIMENSION(:,:) :: rs_data_container
  TYPE(radar_data_type), POINTER, DIMENSION(:) :: rs_data => NULL()

  ! array of auxiliary azimut grid data (local PE)
  TYPE(radar_grid_type), ALLOCATABLE, TARGET, DIMENSION(:,:) :: rs_grid_container
  TYPE(radar_grid_type), POINTER, DIMENSION(:) :: rs_grid => NULL()

  ! cartesian coordinates for superobing (the whole model domain)
  TYPE(cart_data_type), ALLOCATABLE, TARGET, DIMENSION(:,:)  :: cart_data_container
  TYPE(cart_data_type), POINTER, DIMENSION(:)  :: cart_data => NULL()

  ! meta data for configuration of grid point reflectivity calculation
  ! for each radar
  TYPE(t_dbzcalc_params), ALLOCATABLE, TARGET, DIMENSION(:,:)  :: dbz_meta_container
  TYPE(t_dbzcalc_params), POINTER, DIMENSION(:)  :: dbz_meta => NULL()

  ! meta data for configuration of the rotated lat/lon grid for radar composites
  TYPE(composite_meta_type), ALLOCATABLE, TARGET, DIMENSION(:)  ::   comp_meta_container
  TYPE(composite_meta_type), POINTER                            ::   comp_meta => NULL()

  ! meta data for configuration of the rotated lat/lon grid for the special radar composite for bubble generator
  TYPE(composite_meta_type), ALLOCATABLE, TARGET, DIMENSION(:)  ::   comp_meta_bub_container
  TYPE(composite_meta_type), POINTER                            ::   comp_meta_bub => NULL()

  TYPE(bubble_list_type), ALLOCATABLE, TARGET, DIMENSION(:)  ::   bubble_list_container
  TYPE(bubble_list_type), POINTER                            ::   bubble_list => NULL()

  ! time levels used in src_radar.f90 for the model fields u,v,w,T,p,qX,qnX
  !    (default initialistion "1", so that init_lookup_mie() can be called also
  !     from other modules like organize_data.f90 with a defined itlrad. The exact
  !     value does not matter, because those calls are only for generating lookup tables,
  !     not computing any meaningful reflectivity)
  INTEGER                     :: itlrad_dyn = 1
  INTEGER                     :: itlrad_qx  = 1

  ! number of radar stations (will be set by the program, either to nradsta_namelist (runs without observation files) or to
  !  the number of radar stations found in observation files:
  INTEGER                     :: nradsta

  ! overlap for the auxiliary azimut-slice grid at the lower and upper azimut ends
  INTEGER                     :: nbl_az


!======================================================================================
!
! Module procedures
!
!======================================================================================

CONTAINS

  !======================================================================================
  !
  ! Routines to initialize the character variables with blanks in the radar_meta_type
  !
  !======================================================================================

  SUBROUTINE rsm_init_strings_blanks_multime (rsmm)

    TYPE(radar_meta_type), INTENT(inout) :: rsmm

    rsmm%station_name(:) = ' '
    rsmm%scanname(:) = ' '
    rsmm%obs_cdate(:)(:) = ' '
    rsmm%obsfile(:,:)(:) = ' '
    rsmm%obsfile_format(:)(:) = ' '
    rsmm%fdbkfile(:) = ' ' 
    rsmm%obs_hdf5_varname_vrad(:) = ' '
    rsmm%obs_hdf5_varname_dbzh(:) = ' '
    rsmm%obs_hdf5_varname_zdr(:) = ' '
    rsmm%obs_hdf5_varname_rhv(:) = ' '
    rsmm%obs_hdf5_varname_kdp(:) = ' '
    rsmm%obs_hdf5_varname_phidp(:) = ' '
    rsmm%obs_hdf5_varname_ldr(:) = ' '
    rsmm%obs_hdf5_varname_cflags(:) = ' '
    rsmm%dummy_buf_for_memalign(:) = ' '
    
  END SUBROUTINE rsm_init_strings_blanks_multime
  
  SUBROUTINE rsm_init_strings_blanks_onetime (rsmo)

    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsmo

    rsmo%station_name(:) = ' '
    rsmo%scanname(:) = ' '
    rsmo%obs_cdate(:)(:) = ' '
    rsmo%obsfile(:,:)(:) = ' '
    rsmo%obsfile_format(:)(:) = ' '
    rsmo%fdbkfile(:) = ' ' 
    rsmo%obs_hdf5_varname_vrad(:) = ' '
    rsmo%obs_hdf5_varname_dbzh(:) = ' '
    rsmo%obs_hdf5_varname_zdr(:) = ' '
    rsmo%obs_hdf5_varname_rhv(:) = ' '
    rsmo%obs_hdf5_varname_kdp(:) = ' '
    rsmo%obs_hdf5_varname_phidp(:) = ' '
    rsmo%obs_hdf5_varname_ldr(:) = ' '
    rsmo%obs_hdf5_varname_cflags(:) = ' '
    rsmo%dummy_buf_for_memalign(:) = ' '

  END SUBROUTINE rsm_init_strings_blanks_onetime
  
  SUBROUTINE rsm_init_strings_blanks_multimv (rsmm)

    TYPE(radar_meta_type), INTENT(inout) :: rsmm(:)
    INTEGER :: i

    DO i=LBOUND(rsmm, dim=1), UBOUND(rsmm, dim=1)
      rsmm(i)%station_name(:) = ' '
      rsmm(i)%scanname(:) = ' '
      rsmm(i)%obs_cdate(:)(:) = ' '
      rsmm(i)%obsfile(:,:)(:) = ' '
      rsmm(i)%obsfile_format(:)(:) = ' '
      rsmm(i)%fdbkfile(:) = ' '
      rsmm(i)%obs_hdf5_varname_vrad(:) = ' '
      rsmm(i)%obs_hdf5_varname_dbzh(:) = ' '
      rsmm(i)%obs_hdf5_varname_zdr(:) = ' '
      rsmm(i)%obs_hdf5_varname_rhv(:) = ' '
      rsmm(i)%obs_hdf5_varname_kdp(:) = ' '
      rsmm(i)%obs_hdf5_varname_phidp(:) = ' '
      rsmm(i)%obs_hdf5_varname_ldr(:) = ' '
      rsmm(i)%obs_hdf5_varname_cflags(:) = ' '
      rsmm(i)%dummy_buf_for_memalign(:) = ' '
    END DO

  END SUBROUTINE rsm_init_strings_blanks_multimv
  
  SUBROUTINE rsm_init_strings_blanks_onetimv (rsmo)

    TYPE(radar_meta_type_onetime), INTENT(inout) :: rsmo(:)
    INTEGER :: i

    DO i=LBOUND(rsmo, dim=1), UBOUND(rsmo, dim=1)
      rsmo(i)%station_name(:) = ' '
      rsmo(i)%scanname(:) = ' '
      rsmo(i)%obs_cdate(:)(:) = ' '
      rsmo(i)%obsfile(:,:)(:) = ' '
      rsmo(i)%obsfile_format(:)(:) = ' '
      rsmo(i)%fdbkfile(:) = ' '
      rsmo(i)%obs_hdf5_varname_vrad(:) = ' '
      rsmo(i)%obs_hdf5_varname_dbzh(:) = ' '
      rsmo(i)%obs_hdf5_varname_zdr(:) = ' '
      rsmo(i)%obs_hdf5_varname_rhv(:) = ' '
      rsmo(i)%obs_hdf5_varname_kdp(:) = ' '
      rsmo(i)%obs_hdf5_varname_phidp(:) = ' '
      rsmo(i)%obs_hdf5_varname_ldr(:) = ' '
      rsmo(i)%obs_hdf5_varname_cflags(:) = ' '
      rsmo(i)%dummy_buf_for_memalign(:) = ' '
    END DO

  END SUBROUTINE rsm_init_strings_blanks_onetimv
  
  !======================================================================================
  !
  ! Functions to convert back and forth an instance of type radar_meta_type "rsmm" (multiple output time steps)
  !  to radar_meta_type_onetime "rsmo" (only one output time step).
  !
  ! All type components are simply copied, with the following exceptions:
  !  - nobs_times is set to 1 if it is > 0.
  !  - only the first element of all fields having a dimension nobstimes_max
  !    is retained along this dimension.
  !
  !======================================================================================

  FUNCTION rsm_multitime2onetime (rsmm) RESULT (rsmo)

    TYPE(radar_meta_type), INTENT(in) :: rsmm   ! Input
    TYPE(radar_meta_type_onetime)     :: rsmo   ! Result

    rsmo%lambda                       = rsmm%lambda
    rsmo%lat                          = rsmm%lat
    rsmo%lon                          = rsmm%lon
    rsmo%alt_agl_mod                  = rsmm%alt_agl_mod
    rsmo%alt_msl                      = rsmm%alt_msl
    rsmo%alt_msl_true                 = rsmm%alt_msl_true
    rsmo%msl_mod                      = rsmm%msl_mod
    rsmo%az_inc                       = rsmm%az_inc
    rsmo%az_start                     = rsmm%az_start
    rsmo%ra_inc                       = rsmm%ra_inc
    rsmo%ra_inc_obs                   = rsmm%ra_inc_obs
    rsmo%ra_inc_coarse                = rsmm%ra_inc_coarse
    rsmo%Theta3                       = rsmm%Theta3
    rsmo%Phi3                         = rsmm%Phi3
    rsmo%dalpha                       = rsmm%dalpha
    rsmo%alpha3_eff_0                 = rsmm%alpha3_eff_0
    rsmo%smth_interv_fact             = rsmm%smth_interv_fact
    rsmo%dt_obs                       = rsmm%dt_obs
    rsmo%dt_obs_fdbk                  = rsmm%dt_obs_fdbk
    rsmo%dt_obs_voldata               = rsmm%dt_obs_voldata
    rsmo%ext_nyq(:,1)                 = rsmm%ext_nyq(:,1)
    rsmo%high_nyq                     = rsmm%high_nyq
    rsmo%prf                          = rsmm%prf
    rsmo%dualprf_ratio                = rsmm%dualprf_ratio
    rsmo%rngate_len                   = rsmm%rngate_len
    rsmo%mds_Z0                       = rsmm%mds_Z0
    rsmo%mds_r0                       = rsmm%mds_r0
    rsmo%obs_times(1)                 = rsmm%obs_times(1)
    rsmo%obs_times_obs(1)             = rsmm%obs_times_obs(1)
    rsmo%obs_times_fdbk(1)            = rsmm%obs_times_fdbk(1)
    rsmo%obs_times_voldata(1)         = rsmm%obs_times_voldata(1)
    rsmo%el_arr                       = rsmm%el_arr
    rsmo%el_arr_default               = rsmm%el_arr_default
    rsmo%xabscsm_v                    = rsmm%xabscsm_v
    rsmo%xabscsm_h                    = rsmm%xabscsm_h
    rsmo%weigsm_v                     = rsmm%weigsm_v
    rsmo%weigsm_h                     = rsmm%weigsm_h
    rsmo%vnyq_min_for_vr_active_fdbk  = rsmm%vnyq_min_for_vr_active_fdbk
    rsmo%icountry                     = rsmm%icountry
    rsmo%station_id                   = rsmm%station_id
    rsmo%ista                         = rsmm%ista
    rsmo%nel                          = rsmm%nel
    rsmo%nel_present                  = rsmm%nel_present
    rsmo%nel_fdbk                     = rsmm%nel_fdbk
    rsmo%nel_voldata                  = rsmm%nel_voldata
    rsmo%nel_default                  = rsmm%nel_default
    rsmo%eleind_for_composite_bub     = rsmm%eleind_for_composite_bub
    rsmo%naz                          = rsmm%naz
    rsmo%naz_ncdf(1,:)                = rsmm%naz_ncdf(1,:)
    rsmo%nra                          = rsmm%nra
    rsmo%nra_obs                      = rsmm%nra_obs
    rsmo%n_aggr_ra_obs                = rsmm%n_aggr_ra_obs
    rsmo%nobs_times                   = MIN(rsmm%nobs_times, 1)
    rsmo%nobs_times_obs               = MIN(rsmm%nobs_times_obs, 1)
    rsmo%nobs_times_fdbk              = MIN(rsmm%nobs_times_fdbk, 1)
    rsmo%nobs_times_voldata           = MIN(rsmm%nobs_times_voldata, 1)
    rsmo%i_nearest_mod                = rsmm%i_nearest_mod
    rsmo%j_nearest_mod                = rsmm%j_nearest_mod
    rsmo%ngpsm_v                      = rsmm%ngpsm_v
    rsmo%ngpsm_h                      = rsmm%ngpsm_h
    rsmo%nrep_ncdf                    = rsmm%nrep_ncdf
    rsmo%num_gates                    = rsmm%num_gates
    rsmo%num_pulses                   = rsmm%num_pulses
    rsmo%obs_startrec(1,:)            = rsmm%obs_startrec(1,:)
    rsmo%obs_endrec(1,:)              = rsmm%obs_endrec(1,:)
    rsmo%ind_ele_present              = rsmm%ind_ele_present
    rsmo%ind_ele_fdbk                 = rsmm%ind_ele_fdbk
    rsmo%ind_ele_voldata              = rsmm%ind_ele_voldata
    rsmo%eleindlist_for_composite     = rsmm%eleindlist_for_composite
    rsmo%lobs_avail(1,:)              = rsmm%lobs_avail(1,:)
    rsmo%lfdbkfile_exist              = rsmm%lfdbkfile_exist
    rsmo%lobstimes_ovwrt_recalc       = rsmm%lobstimes_ovwrt_recalc
    rsmo%lvrad_to_fdbk                = rsmm%lvrad_to_fdbk
    rsmo%ldbzh_to_fdbk                = rsmm%ldbzh_to_fdbk
    rsmo%station_name                 = rsmm%station_name
    rsmo%scanname                     = rsmm%scanname
    rsmo%obs_cdate(1)                 = rsmm%obs_cdate(1)
    rsmo%obsfile(1,:)                 = rsmm%obsfile(1,:)
    rsmo%obsfile_format(1)            = rsmm%obsfile_format(1)
    rsmo%fdbkfile                     = rsmm%fdbkfile
    rsmo%obs_hdf5_varname_vrad        = rsmm%obs_hdf5_varname_vrad 
    rsmo%obs_hdf5_varname_dbzh        = rsmm%obs_hdf5_varname_dbzh 
    rsmo%obs_hdf5_varname_zdr         = rsmm%obs_hdf5_varname_zdr  
    rsmo%obs_hdf5_varname_rhv         = rsmm%obs_hdf5_varname_rhv  
    rsmo%obs_hdf5_varname_kdp         = rsmm%obs_hdf5_varname_kdp  
    rsmo%obs_hdf5_varname_phidp       = rsmm%obs_hdf5_varname_phidp
    rsmo%obs_hdf5_varname_ldr         = rsmm%obs_hdf5_varname_ldr  
    rsmo%obs_hdf5_varname_cflags      = rsmm%obs_hdf5_varname_cflags

  END FUNCTION rsm_multitime2onetime

  FUNCTION rsm_onetime2multitime (rsmo) RESULT (rsmm)

    TYPE(radar_meta_type_onetime), INTENT(in) :: rsmo   ! Input
    TYPE(radar_meta_type)                     :: rsmm   ! Result

    rsmm%lambda                       = rsmo%lambda
    rsmm%lat                          = rsmo%lat
    rsmm%lon                          = rsmo%lon
    rsmm%alt_agl_mod                  = rsmo%alt_agl_mod
    rsmm%alt_msl                      = rsmo%alt_msl
    rsmm%alt_msl_true                 = rsmo%alt_msl_true
    rsmm%msl_mod                      = rsmo%msl_mod
    rsmm%az_inc                       = rsmo%az_inc
    rsmm%az_start                     = rsmo%az_start
    rsmm%ra_inc                       = rsmo%ra_inc
    rsmm%ra_inc_obs                   = rsmo%ra_inc_obs
    rsmm%ra_inc_coarse                = rsmo%ra_inc_coarse
    rsmm%Theta3                       = rsmo%Theta3
    rsmm%Phi3                         = rsmo%Phi3
    rsmm%dalpha                       = rsmo%dalpha
    rsmm%alpha3_eff_0                 = rsmo%alpha3_eff_0
    rsmm%smth_interv_fact             = rsmo%smth_interv_fact
    rsmm%dt_obs                       = rsmo%dt_obs
    rsmm%dt_obs_fdbk                  = rsmo%dt_obs_fdbk
    rsmm%dt_obs_voldata               = rsmo%dt_obs_voldata
    rsmm%ext_nyq(:,:)                 = -999.9_wpfwo
    rsmm%ext_nyq(:,1)                 = rsmo%ext_nyq(:,1) 
    rsmm%high_nyq                     = rsmo%high_nyq
    rsmm%prf                          = rsmo%prf
    rsmm%dualprf_ratio                = rsmo%dualprf_ratio
    rsmm%rngate_len                   = rsmo%rngate_len
    rsmm%mds_Z0                       = rsmo%mds_Z0
    rsmm%mds_r0                       = rsmo%mds_r0
    rsmm%obs_times(:)                 = -999.9_wpfwo
    rsmm%obs_times(1)                 = rsmo%obs_times(1)
    rsmm%obs_times_obs(:)             = -999.9_wpfwo
    rsmm%obs_times_obs(1)             = rsmo%obs_times_obs(1)
    rsmm%obs_times_fdbk(:)            = -999.9_wpfwo
    rsmm%obs_times_fdbk(1)            = rsmo%obs_times_fdbk(1)
    rsmm%obs_times_voldata(:)         = -999.9_wpfwo
    rsmm%obs_times_voldata(1)         = rsmo%obs_times_voldata(1)
    rsmm%el_arr                       = rsmo%el_arr
    rsmm%el_arr_default               = rsmo%el_arr_default
    rsmm%xabscsm_v                    = rsmo%xabscsm_v
    rsmm%xabscsm_h                    = rsmo%xabscsm_h
    rsmm%weigsm_v                     = rsmo%weigsm_v
    rsmm%weigsm_h                     = rsmo%weigsm_h
    rsmm%vnyq_min_for_vr_active_fdbk  = rsmo%vnyq_min_for_vr_active_fdbk
    rsmm%icountry                     = rsmo%icountry
    rsmm%station_id                   = rsmo%station_id
    rsmm%ista                         = rsmo%ista
    rsmm%nel                          = rsmo%nel
    rsmm%nel_present                  = rsmo%nel_present
    rsmm%nel_fdbk                     = rsmo%nel_fdbk
    rsmm%nel_voldata                  = rsmo%nel_voldata
    rsmm%nel_default                  = rsmo%nel_default
    rsmm%eleind_for_composite_bub     = rsmo%eleind_for_composite_bub
    rsmm%naz                          = rsmo%naz
    rsmm%naz_ncdf(:,:)                = 0
    rsmm%naz_ncdf(1,:)                = rsmo%naz_ncdf(1,:)
    rsmm%nra                          = rsmo%nra
    rsmm%nra_obs                      = rsmo%nra_obs
    rsmm%n_aggr_ra_obs                = rsmo%n_aggr_ra_obs
    rsmm%nobs_times                   = rsmo%nobs_times
    rsmm%nobs_times_obs               = rsmo%nobs_times_obs
    rsmm%nobs_times_fdbk              = rsmo%nobs_times_fdbk
    rsmm%nobs_times_voldata           = rsmo%nobs_times_voldata
    rsmm%i_nearest_mod                = rsmo%i_nearest_mod
    rsmm%j_nearest_mod                = rsmo%j_nearest_mod
    rsmm%ngpsm_v                      = rsmo%ngpsm_v
    rsmm%ngpsm_h                      = rsmo%ngpsm_h
    rsmm%nrep_ncdf                    = rsmo%nrep_ncdf
    rsmm%num_gates                    = rsmo%num_gates
    rsmm%num_pulses                   = rsmo%num_pulses
    rsmm%obs_startrec(:,:)            = missval_int
    rsmm%obs_startrec(1,:)            = rsmo%obs_startrec(1,:)
    rsmm%obs_endrec(:,:)              = missval_int
    rsmm%obs_endrec(1,:)              = rsmo%obs_endrec(1,:)
    rsmm%ind_ele_present              = rsmo%ind_ele_present
    rsmm%ind_ele_fdbk                 = rsmo%ind_ele_fdbk
    rsmm%ind_ele_voldata              = rsmo%ind_ele_voldata
    rsmm%eleindlist_for_composite     = rsmo%eleindlist_for_composite
    rsmm%lobs_avail(:,:)              = .FALSE.
    rsmm%lobs_avail(1,:)              = rsmo%lobs_avail(1,:)
    rsmm%lfdbkfile_exist              = rsmo%lfdbkfile_exist
    rsmm%lobstimes_ovwrt_recalc       = rsmo%lobstimes_ovwrt_recalc
    rsmm%lvrad_to_fdbk                = rsmo%lvrad_to_fdbk
    rsmm%ldbzh_to_fdbk                = rsmo%ldbzh_to_fdbk
    rsmm%station_name                 = rsmo%station_name
    rsmm%scanname                     = rsmo%scanname
    rsmm%obs_cdate(:)                 = ' '
    rsmm%obs_cdate(1)                 = rsmo%obs_cdate(1)
    rsmm%obsfile(:,:)                 = TRIM(obsfile_missingname)
    rsmm%obsfile(1,:)                 = rsmo%obsfile(1,:)
    rsmm%obsfile_format(:)            = ' '
    rsmm%obsfile_format(1)            = rsmo%obsfile_format(1)
    rsmm%fdbkfile                     = rsmo%fdbkfile
    rsmm%obs_hdf5_varname_vrad        = rsmo%obs_hdf5_varname_vrad 
    rsmm%obs_hdf5_varname_dbzh        = rsmo%obs_hdf5_varname_dbzh 
    rsmm%obs_hdf5_varname_zdr         = rsmo%obs_hdf5_varname_zdr  
    rsmm%obs_hdf5_varname_rhv         = rsmo%obs_hdf5_varname_rhv  
    rsmm%obs_hdf5_varname_kdp         = rsmo%obs_hdf5_varname_kdp  
    rsmm%obs_hdf5_varname_phidp       = rsmo%obs_hdf5_varname_phidp
    rsmm%obs_hdf5_varname_ldr         = rsmo%obs_hdf5_varname_ldr  
    rsmm%obs_hdf5_varname_cflags      = rsmo%obs_hdf5_varname_cflags

  END FUNCTION rsm_onetime2multitime


  SUBROUTINE prep_domains_radar ( ndoms_model, radar_flag_doms_model)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ndoms_model
    LOGICAL, INTENT(in) :: radar_flag_doms_model(ndoms_model)  ! Flag for each domain whether EMVORADO is applied to this domain or not

    INTEGER :: i

    !=======================================================================================
    ! Initialize the information on model domains with active radar simulator:
    !---------------------------------------------------------------------------------------

    ALLOCATE(list_domains_for_radar(ndoms_model))
    ALLOCATE(list_domains_for_model(ndoms_max  ))

    list_domains_for_radar (:) = -1
    list_domains_for_model (:) = -1
    ndoms_radar = 0
    DO i=1, ndoms_model
      IF (radar_flag_doms_model(i) .AND. ndoms_radar+1 <= ndoms_max) THEN
        ndoms_radar = ndoms_radar + 1
        list_domains_for_radar(i)           = ndoms_radar
        list_domains_for_model(ndoms_radar) = i
      ELSE IF (ndoms_radar+1 > ndoms_max) THEN
        WRITE (*, '(a,i4)') 'ERROR prep_domains_radar(): too many model domains with active radar simulation, '// &
             'allowed are ndoms_max = ', ndoms_max
        STOP
      END IF
    END DO

    !=======================================================================================
    ! Initialize the information on model domains with active radar simulator:
    !---------------------------------------------------------------------------------------

    IF (ndoms_radar > 0) THEN
      ALLOCATE(rs_meta_container       (ndoms_radar,nradsta_max))
      ALLOCATE(rs_data_container       (ndoms_radar,nradsta_max))
      ALLOCATE(rs_grid_container       (ndoms_radar,nradsta_max))
      ALLOCATE(cart_data_container     (ndoms_radar,nradsta_max))
      ALLOCATE(dbz_meta_container      (ndoms_radar,nradsta_max))
      ALLOCATE(comp_meta_container     (ndoms_radar            ))
      ALLOCATE(comp_meta_bub_container (ndoms_radar            ))
      ALLOCATE(bubble_list_container   (ndoms_radar            ))
    END IF

    IF (ndoms_model > 0) THEN
      ALLOCATE(num_radar_dom       (ndoms_model))
      ALLOCATE(num_radario_dom     (ndoms_model))
      ALLOCATE(radario_master_dom  (ndoms_model))
      ALLOCATE(icomm_radar_dom     (ndoms_model))
      ALLOCATE(icomm_radario_dom   (ndoms_model))
      ALLOCATE(lradar_pe_dom       (ndoms_model))
      ALLOCATE(lradario_pe_dom     (ndoms_model))
      ALLOCATE(my_radar_id_dom     (ndoms_model))
      ALLOCATE(my_radario_id_dom   (ndoms_model))
    END IF

  END SUBROUTINE prep_domains_radar

  SUBROUTINE switch_to_domain_radar ( idom_model )
    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom_model

    INTEGER             :: idom_radar

    ! .. Set the EMVORADO-wide domain index storage variable
    !     for the actual model domain on which EMVORADO works at the moment:
    idom       = idom_model

    idom_radar = list_domains_for_radar(idom_model)

    rs_data       => rs_data_container       (idom_radar,:)
    rs_meta       => rs_meta_container       (idom_radar,:)
    rs_grid       => rs_grid_container       (idom_radar,:)
    cart_data     => cart_data_container     (idom_radar,:)
    dbz_meta      => dbz_meta_container      (idom_radar,:)
    comp_meta     => comp_meta_container     (idom_radar  )
    comp_meta_bub => comp_meta_bub_container (idom_radar  )
    bubble_list   => bubble_list_container   (idom_radar  )

  END SUBROUTINE switch_to_domain_radar


END MODULE radar_data

