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

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_data_namelist

!------------------------------------------------------------------------------
!
! Description:
!   This module provides global variables for the namelist parameters
!   of the radar forward operator EMVORADO.
!
! Method:
!   Declarations of global variables
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp, sp, wpfwo

  USE radar_data, ONLY :   &
       cmaxlen, &
       radar_meta_type, t_dbzcalc_params, &
       nel_max, nel_composite_max, noutstreams_max, &
       nobstimes_max, noutput_fields_max, &
       nradsta, nbl_az, ndoms_radar, ndoms_max, &
       list_domains_for_radar


  IMPLICIT NONE

  PUBLIC

  ! Default path and name for the EMVORADO namelist file(s):
#ifdef __COSMO__
  CHARACTER(len=cmaxlen)  :: radarnmlfile     = 'INPUT_RADARSIM' // REPEAT(' ',cmaxlen)
#endif
#ifdef __ICON__
  CHARACTER(len=cmaxlen)  :: radarnmlfile     = 'NAMELIST_EMVORADO' // REPEAT(' ',cmaxlen)
#endif
  CHARACTER(len=cmaxlen)  :: radarnmlfilepath = REPEAT(' ', cmaxlen)

  ! Default path and name for the EMVORADO namelist control output:
#ifdef __COSMO__
  CHARACTER(len=cmaxlen)  :: nmlctrlfile      = 'YUSPECIF_RADAR' // REPEAT(' ',cmaxlen)
#endif
#ifdef __ICON__
  CHARACTER(len=cmaxlen)  :: nmlctrlfile      = 'nml.emvorado.log' // REPEAT(' ',cmaxlen)
#endif
  CHARACTER(len=cmaxlen)  :: nmlctrlfilepath  = REPEAT(' ',cmaxlen)


  ! .. number of radar stations in namelist file:
  INTEGER                               :: nradsta_namelist

  ! .. switch to enable ASCII-output (default = .true.):
  !     (disabling just prevents the output-step in subroutine output3d_ascii_radar(),
  !       not the computations and the output in the diagnosis files!)
  LOGICAL                    :: lvoldata_output

  ! .. Derived type to hold informations for one voldata output stream:
  TYPE t_voldata_ostream
    SEQUENCE
    ! .. Format of the volume data output:  'ascii', 'f90-binary', 'ncdf'
    !     If 'f90-binary' is used, the resulting Fortran binary files (.bin) need
    !     a conversion by the PROGRAM "bin2ascii_convrates" from U. Blahak
    CHARACTER(LEN=12)              :: format
    CHARACTER(len=12)              :: grib2_packingtype    ! 'grid_simple', 'grid_ccsds' (default), 'png'
    ! .. Filename pattern for output in cdfin or grib2 formats. Ascii and bin have predefined filenames.
    !     May consist of the following:
    !     - Any valid characters for filenames
    !     - Patterns for radar-ID, time and parameter name identification:
    !       - <stationid>    : will result in e.g. 'id-123456'
    !       - <varname>      : will result in e.g. 'zrsim' or 'all', depending on format 
    !       - <tstart>       : will result in YYYYMMDDhhmmss  (absolute datetime)
    !       - <tend>         : will result in YYYYMMDDhhmmss  (absolute datetime)
    !       - <tvvz>         : will result in DDhhmmss  (forecast lead time)
    !       - <scantype>     : will result in 'volscan' or 'precipscan' by default, but can be defined by pat_scantype_XXX below
    !     The extention will be chosen according to the data format:
    !     - '.nc' for 'cdfin' and 'cdfin-mulmom'
    !     - no extention for 'grib2'
    !    Some patterns are mandatory
    CHARACTER(LEN=cmaxlen)         :: file_pattern
    CHARACTER(len=12)              :: pat_scantype_volscan, pat_scantype_precipscan
    CHARACTER(LEN=cmaxlen)         :: output_subdir
    ! .. List of output fields:  possible names are e.g.
    !     'zrsim', 'zrobs', 'vrsim', 'vrobs',
    !     'ersim', 'losim', 'lasim', 'hrsim'
    !    For complete list see radar_namlist_read.f90
    CHARACTER(LEN=12)              :: output_list(noutput_fields_max)
    ! .. namelist parameters for defining the time batches contained in output "cdfin" or "grib2" files
    REAL (KIND=wpfwo)              :: content_dt    ! length of time batches in output-files [s]
    REAL (KIND=wpfwo)              :: content_tref  ! reference time for synchronization of the
                                                    !  batches in output-files [seconds sice model start]
  END TYPE t_voldata_ostream
  ! .. Vector of output streams
  TYPE(t_voldata_ostream)          :: voldata_ostream(noutstreams_max)

  !.. Switch to enable output of READY files for each EMVORADO output time step:
  LOGICAL                    :: lwrite_ready
  CHARACTER (len=cmaxlen)    :: ready_file_pattern    ! file pattern containing (optional and mandatory) <keys>

  ! .. switches to choose non-blocking mpi_isend - mpi_irecv - mpi_wait communication
  !    (=more efficient, but needs more memory)for:
  !    a) collecting radar data on radar-IO-PEs (for voldata, fdbk, composites, bubble generator)
  !         (comm. from workers -> radar-IO-PEs)
  !    b) collecting interpolated model data on the azimuthal slices for online beam propagation
  !         (comm. from workers <-> workers)
  LOGICAL                    :: lcomm_nonblocking_output   ! a)
  LOGICAL                    :: lcomm_nonblocking_online   ! b)

  !.. switch to enable radar-param output (and nwp-field output) on model grid (default = .false.):
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
  LOGICAL                    :: lmodfield_output
  LOGICAL                    :: loutnwp
  LOGICAL                    :: loutvolaux

  !.. switch to enable output of netcdf feedback files (only possible with -DNETCDF and -DNUDGING)
  LOGICAL                    :: lfdbk_output
  !.. enable debug info output:
  LOGICAL                    :: ldebug_radsim
  !.. enable calculation and output of radial wind for each radar:
  LOGICAL                    :: loutradwind
  !.. fill areas with simulated missing radial wind by model background radial wind from model U, V, W:
  !    (may happen if lsmooth=.TRUE., so this switch only effective if lsmooth=.TRUE.):
  LOGICAL                    :: lfill_vr_backgroundwind
  !.. enable output of beam geometry (bin heights, local elevations) for each radar:
  LOGICAL                    :: lout_geom
  !.. enable calculation and output of reflectivity for each radar:
  LOGICAL                    :: loutdbz
  !.. if radar moments such as dbz or polar. param. should be computed on radar bins or not.
  !    If .true., the model state is first interpolated to radarbins, before moments are derived from them.
  !    Otherwise, moments are computed on model grid and then interpolated to radarbins (old default method).
  LOGICAL                    :: lcalc_dbz_on_radarbins
  !.. enable calculation and output of standard polarimetric
  !     parameters (ZDR, KDP, RHV) for each radar:
  LOGICAL                    :: loutpolstd
  !.. enable calculation and output of additional polarimetric
  !     parameters (LDR) for each radar:
  LOGICAL                    :: loutpolall
  !.. enable calculation and consideration of extinction effects on
  !     reflectivity and polarimetric parameters for each radar:
  LOGICAL                    :: lextdbz
  !.. enable to take the reflectivity as weight into
  !     account when calculating the terminal fall speed of hydrometeors (lfall=.true.) as well as if lsmooth=.true.:
  LOGICAL                    :: lweightdbz
  !.. enable to take the falling velocity of hydrometeors into account
  !     while caculating the radial wind for each radar:
  LOGICAL                    :: lfall
  !.. enable to change the strategy for computing of radio propagation path:
  !     static 4/3-earth ("offline") or
  !     dynamic ("online") depending on actual refractive index of air
  LOGICAL                    :: lonline, lsode
  !.. enable spatial smoothing (weighted spatial mean over measuring
  !     volume increasing with distance):
  LOGICAL                    :: lsmooth
  !.. enable to honor the minimum detectable signal (quadratic function of range), i.e.,
  !     if Z < Z0(r), then set zero_value ("correct 0") for Ze and miss_value for vr ("no valid measurement"):
  LOGICAL                    :: lmds_z
  LOGICAL                    :: lmds_vr
  !.. meta data of radar station will be read from netCDF data files, overriding
  !   any radar meta data from the namelist:
  LOGICAL                    :: lreadmeta_from_netcdf
  !.. if lreadmeta_from_netcdf=.true., check the record numbers in files for
  !   vr, qv, z, qz for each station on correspondence in time. If not,
  !   set record numbers for corresponding times for each data type file
  !   separately (VERY SLOW ON VECTOR MACHINES!)
  !   This might help, if some records are missing in an input file (e.g., for qv) and
  !   are present in another file (e.g., vr) for the same station.
  LOGICAL                    :: lcheck_inputrecords
  LOGICAL                    :: lequal_azi_alldatasets

  !.. whether or not a model run should abort if serious problems with
  !    required observation files or meta data occur, i.e., no obs files at all,
  !    errors in file content, missing variables, missing or wrong station ID or scan strategy, etc.
  !   In operational runs, it should be up to the user if the run should abort or continue, therefore the switch:
  LOGICAL                    :: labort_if_problems_obsfiles

  
  !.. whether or not a model run should abort if serious problems with
  !   output gribs are encountered
  LOGICAL                    :: labort_if_problems_gribout

  !.. flag to take into account quality control flags for the observations from NetCDF files:
  LOGICAL                    :: lqc_flag

  !.. flag for dealiasing observations from NetCDF files:
  LOGICAL                    :: ldealiase_vr_obs

  ! .. Testmode for reflectivity calculation in the first time step for idealized cases / testsuite:
  LOGICAL                    :: ltestpattern_hydrometeors

  ! .. Compute global composites of obs and sim reflectivity on the model grid
  !     for later use in the NWP-model, using data of one specific elevation of each station:
  LOGICAL                    ::   ldo_composite                  ! Create radar composites

  ! .. Subswitch of ldo_composite, to trigger own grib2 output of the composites
  !    aside from the model-internal grib output facilities. This makes it independent
  !    of the model grid and is important for asynchroneous radar IO:
  LOGICAL                    :: lcomposite_output
  CHARACTER(len=12)          :: comp_grib2_packingtype    ! 'grid_simple', 'grid_ccsds' (default), 'png'
  CHARACTER (len=cmaxlen)    :: composite_file_pattern    ! file pattern containing (optional and mandatory) <keys>

  ! .. Automatic and manual warm bubbles:
  LOGICAL             ::   &
       ldo_bubbles,        & ! If .true., automatically detect the need for warm bubbles from composites; requires ldo_composite=.TRUE. in radar namelist and works together with ldo_bubbles in the driving model
       ldo_bubbles_manual    ! If .true., allows manual specification of bubbles in real case runs (lartif_data=.false.) using the bubble parameters from the ARTIFDATA namelist
  LOGICAL                    :: lcomposite_output_bub
  CHARACTER (len=cmaxlen)    :: composite_file_pattern_bub    ! file pattern containing (optional and mandatory) <keys>
  ! Parameters for detection and downstream advection of automatic bubbles
  REAL (kind=wpfwo)          :: dt_bubble_search             ! Time interval from one automatic bubble search to the next [seconds] (SYNCHRONIZE WITH COMPOSITE TIMES!!!)
  REAL (kind=wpfwo)          :: tstart_bubble_search         ! Start time for the bubble search, and convection trigger [seconds]   (SYNCHRONIZE WITH COMPOSITE TIMES!!!)
  REAL (kind=wpfwo)          :: tend_bubble_search           ! End time for the bubble search and convection trigger [seconds]      (SYNCHRONIZE WITH COMPOSITE TIMES!!!)

  REAL (kind=wpfwo)          :: t_offset_bubble_trigger_async! In case of asynchr. radar IO for runtime optimization:
                                                             !   Time delay of bubble triggering after detection step.
                                                             !   It adds to the bubble advection time dt_bubble_advect below.
                                                             !   The worker PEs of the model may continue in parallel to the async IO PEs
                                                             !   for this amount of model time, until they synchronize with the IO PEs to
                                                             !   to gather information on new bubbles. An offset of 0.0 would be physically ideal but
                                                             !   but destroys the runtime advantage of asynchr. IO.
  REAL (kind=wpfwo)          :: prob_bubble                  ! Probability of triggering a bubble when it is found [0-1]
  REAL (kind=wpfwo)          :: maxdim_obs                   ! Maximum dimension of observed cells that can trigger a bubble (very large objects are not targeted) [meters]
  LOGICAL                    :: lbub_isolated                       ! Check that the bubble candidate objects are isolated from other objects in observations
  REAL (kind=wpfwo)          :: threshold_obs(2), threshold_mod(2)  ! Thresholds that define the minimum dBZ in an object- and the high intensity region [dBZ]
  REAL (kind=wpfwo)          :: areamin_mod(2), areamin_obs(2)      ! Minimum area of the object- and high intensity regions for detecting an object for detecting an object [m^2]
  REAL (kind=wpfwo)          :: mult_dist_obs                       ! Multiplicative axis factor for the obs object [-]. If it is high, separation- and isolation checks are more restrictive
  REAL (kind=wpfwo)          :: mult_dist_mod                       ! Multiplicative axis factor for the sim object [-]. If it is high, separation- and isolation checks are more restrictive
  REAL (kind=wpfwo)          :: add_dist_obs                        ! Additive axis increase for the obs object [meters]. If it is high, separation- and isolation checks are more restrictive
  REAL (kind=wpfwo)          :: add_dist_mod                        ! Additive axis increase for the sim object [meters]. If it is high, separation- and isolation checks are more restrictive
  REAL (kind=wpfwo)          :: dt_bubble_advect             ! Time scale for downstream advection of automatic bubbles to compensate for transport effects during initial bubble formation [seconds]
  REAL (kind=wpfwo)          :: zlow_meanwind_bubble_advect  ! The lower bound of averaging height interval for bubble advection speed [meters AMSL]
  REAL (kind=wpfwo)          :: zup_meanwind_bubble_advect   ! The upper bound of averaging height interval for bubble advection speed [meters AMSL]
  CHARACTER (len=12)         :: bubble_type         ! Type of perturbation 'cos-hrd', 'cos-instant'
  REAL (kind=wpfwo)          :: bubble_heatingrate  ! Constant heating rate for the bubbles of type 'cos-hrd' [K/s]
  REAL (kind=wpfwo)          :: bubble_timespan     ! Total timespan for heating the bubbles of type 'cos-hrd' [s]
  REAL (kind=wpfwo)          :: bubble_dT           ! Temperature disturbance for the bubbles of type 'cos-instant' [K]
  REAL (kind=wpfwo)          :: bubble_centz        ! Center height (Z) of temperature disturbances [m AGL]
  REAL (kind=wpfwo)          :: bubble_radx         ! Horizontal main axis (radius) in X (longitudinal) direction of temperature disturbances [m]
  REAL (kind=wpfwo)          :: bubble_rady         ! Horizontal main axis (radius) in Y (latitudinal) direction of temperature disturbances [m]
  REAL (kind=wpfwo)          :: bubble_radz         ! Vertical main axis (radius) of temperature disturbances [m]
  REAL (kind=wpfwo)          :: bubble_rotangle     ! Rotation angle of the horizontal main axes of temperature disturbances [degrees]
  REAL (kind=wpfwo)          :: bubble_dT_noise     ! In case of ladd_bubblenoise_t=.true., relative noise level, such that
                                                    !   dT          = dT          * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-instant')
                                                    !   heatingrate = heatingrate * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-hrd')
                                                    ! Whether or not to keep the relative humidity constant during heating

  LOGICAL                    :: bubble_holdrhconst  ! Switch to choose if RH should kept constant during heating or not
  LOGICAL                    :: bubble_addnoise_T   ! Switch to activate some random noise on the bubbles with a
                                                    !  relative amplitude of "bubble_dT_noise" (FEATURE NOT YET IMPLEMENTED IN COSMO/ICON)

  ! .. Namelist parameter to define all reflectivty related parameters
  !    in a global way for all radars.
  !      (THIS REPLACES THE PREVIOUS NAMELIST PARAMETERS
  !       ITYPE_REFL_GLOB AND LLOOKUP_MIE_GLOB. BOTH ARE COMPONENTS
  !       OF DBZ_META_GLOB)
  TYPE(t_dbzcalc_params)      :: dbz_meta_glob

  ! .. Namelist parameters for global ngpsm_h and ngpsm_v default values.
  !    These are the background values for all radars, which can
  !    later be overwritten for the single stations via namelist.
  !    If no values or -99 is given, the background default values
  !    1 resp. 1 will win, which means "no smoothing"!
  INTEGER                     :: ngpsm_h_glob
  INTEGER                     :: ngpsm_v_glob


  ! namelist setting for superobing: 0: no superobing; 1: horizontal averaging;  2: horizontal median
  INTEGER                        :: itype_supobing
  ! lower threshold for number of radar bins used for superobbing:
  INTEGER                        :: supob_nrb
  REAL (KIND=wpfwo)              :: supob_azi_maxsector_vr  ! maximal azimut sector (symetrical to its center) for v_r superobing [deg]
  REAL (KIND=wpfwo)              :: supob_cart_resolution   ! resoltion of the cartesian grid for superobing [m]
  REAL (KIND=wpfwo)              :: supob_ave_width         ! width of averaging area for superobing [m]
  REAL (KIND=wpfwo)              :: supob_minrange_vr       ! min. range of superobing points for vr [m]
  REAL (KIND=wpfwo)              :: supob_minrange_z        ! min. range of superobing points for z  [m]
  REAL (KIND=wpfwo)              :: supob_vrw   ! threshold for stddev of radial wind values within a superobed bin [m/s]
  REAL (KIND=wpfwo)              :: supob_rfl   ! threshold for stddev of reflectivity values within a superobed bin [mm^6/m^3]
  REAL (KIND=wpfwo)              :: supob_lowthresh_z_obs   ! lower threshold for valid Ze values, to which no rain data are
  REAL (KIND=wpfwo)              :: supob_lowthresh_z_sim   !   set before average_superobing
  INTEGER                        :: itype_obserr_vr         ! Mode to specify observation errors for (superobe'd or original) radial wind in feedback files (fdbk%e_o):
                                                            !   = 0: constant observation error (=baseval_obserr_vr)
                                                            !   = 1: e_o is a linear ramp as function of (superobe'd or original) reflectivity (see sketch below)
                                                            !   = 2: superobe'd e_o + e_o-weighted superobe'd radial wind. For this,
                                                            !        e_o is computed from reflectivity before superobing (same linear ramp as for option 1) and is used:
                                                            !        a) as second weight (inverse e_o) additionally to spatial Cressman-weight in superobbing of radial wind
                                                            !        b) for superobbing of e_o itself with the same weighting
                                                            ! / sub-switches for itype_obserr_vr (see sketch down below):
  REAL (KIND=wpfwo)              :: ramp_lowdbz_obserr_vr,& ! for options   1/2: lower ramp dBZ-threshold for increasing observation error in radial wind (miss_value = neutral value)
                                    ramp_highdbz_obserr_vr  ! for options   1/2: upper ramp dBZ-threshold for increasing observation error in radial wind (miss_value = neutral value)
  REAL (kind=wpfwo)              :: maxval_obserr_vr        ! for options   1/2: obs error for dBZ-values <  ramp_highdbz_obserr_vr
  REAL (kind=wpfwo)              :: baseval_obserr_vr       ! for options 0/1/2: obs error for dBZ-values >= ramp_highdbz_obserr_vr (options 1/2) or general value (option 0)
                                                            !   in general: 1.0 means that the obserr for vr is generally defined as a relative error.
                                                            !               2.5 would be a good choice for an absolute error
  REAL (kind=wpfwo), PARAMETER   :: baseval_obserr_dbz = 1.0_wpfwo  ! base value for reflectivity obs error (fdbk%e_o)

  ! Sketch of observation error for radial wind (e_o) as function of dBZ:
  ! =====================================================================
  !
  ! Linear ramp function:   e_o(dBZ) = e_o_0 + (e_o_1 - e_o_0) * min( max( (dBZ - dBZ_0) / (dBZ_1 - dBZ_0) , 0.0) , 1.0)
  !
  !   e_o ^
  !       |               |             |   with:  dBZ_0 = ramp_lowdbz_obserr_vr  [dBZ] (default: -1000.99 neutral value)
  ! e_o_0 |---------------\             |          dBZ_1 = ramp_highdbz_obserr_vr [dBZ] (default: -999.99  neutral value)
  !       |               | \           |          e_o_0 = maxval_obserr_vr        [-]  (default: 10.0)
  !       |               |   \         |          e_0_1 = baseval_obserr_vr       [-]  (default: 1.0)
  !       |               |     \       |
  !       |               |       \     |
  !       |               |         \   |
  !       |               |           \ |
  ! e_o_1 | - - - - - - - | - - - - - - \------------------------------
  !       |               |             |
  !     -------------------------------------------------------------------->
  !       |             dBZ_0          dBZ_1                             dBZ

  ! namelist parameter for data thinning in NetCDF feedback files:
  !   thinning step withs in azi, range and elevation directions.
  ! NOTE: - range- and azi-thinning only effective if itype_supobing == 0
  !          (no horizontal superobing)
  !       - elevational thinning is always effective
  INTEGER                         :: thin_step_azi, thin_step_range, thin_step_ele

  ! namelist parameter for selecting elevations in feedback files:
  INTEGER                         :: ind_ele_fdbk_glob(nel_max)

  ! namelist parameter for selecting elevations in volume data files:
  INTEGER                         :: ind_ele_voldata_glob(nel_max)

  ! namelist parameters for selecting observation times in feedback files (seconds since model start):
  REAL (KIND=wpfwo)              :: obs_times_fdbk_glob(nobstimes_max)
  REAL (KIND=wpfwo)              :: dt_obs_fdbk_glob(3)  ! notation either <incr>,<miss>,<miss> or <from>,<to>,<incr>

  ! namelist parameters for selecting observation times in volume data files (seconds since model start):
  REAL (KIND=wpfwo)              :: obs_times_voldata_glob(nobstimes_max)
  REAL (KIND=wpfwo)              :: dt_obs_voldata_glob(3)  ! notation either <incr>,<miss>,<miss> or <from>,<to>,<incr>


  ! namelist setting to decide about the metric for reflectivity in the feedback files:
  !  itype_metric_refl_fdbk = 1:  write dBZ to feedback files
  !  itype_metric_refl_fdbk = 2:  convert to effective LWC = 0.004*Z^0.55, but leave obs error factor at 1.0
  !  itype_metric_refl_fdbk = 3:  convert to effective LWC as for (2) and set obs-dependent obs error factor
  !                               with an asymptotic value of <minval_obserr_lwc> for LWC -> 0.0
  INTEGER                        :: itype_metric_refl_fdbk
  REAL (KIND=wpfwo)              :: minval_obserr_lwc           ! [g/m^3]

  !.. Directory where NetCDF obs input resides (from namelist):
  CHARACTER(LEN=cmaxlen)          :: ydirradarin
  !.. Filename of existing file with a listing of the contents of ydirradarin:
  CHARACTER(LEN=cmaxlen)          :: ydirlistfile

  !.. Directory where radar output files are written to:
  CHARACTER(LEN=cmaxlen)          :: ydirradarout   ! main output dir
  CHARACTER(LEN=cmaxlen)          :: ysubdirfof     ! subdir for fof's herein
  CHARACTER(LEN=cmaxlen)          :: ysubdircomp    ! subdir for composites herein

  !.. Directory for reading Mie lookup tables:
  CHARACTER(LEN=cmaxlen)          :: ydir_mielookup_read
  !.. Directory for writing new Mie lookup tables (should normally be equal to ydir_mielookup_read):
  CHARACTER(LEN=cmaxlen)          :: ydir_mielookup_write
    !.. Directory for writing ready files (should normally be equal to ydirradarout):
  CHARACTER(LEN=cmaxlen)          :: ydir_ready_write
  !.. Parallelization strategy for lookup table generation:
  INTEGER                         :: itype_mpipar_lookupgen  ! 1 = parallelization over dbzmeta-sets and hydrometeors
                                                             ! 2 = parallelization over table elements
  !.. Processor number range in icomm_compute_fwo for parallel lookup table generation for all itype_mpipar_lookupgen methods:
  INTEGER                         :: pe_start_lookupgen  ! start PE  [0...num_compute_fwo-1]
  INTEGER                         :: pe_end_lookupgen    ! end PE    [0...num_compute_fwo-1]
  LOGICAL                         :: llookup_interp_mode_dualpol ! if to use exact same lookup table interpolation method
                                                                 ! for all radar moments (cubic in lin-lin-space).
                                                                 ! The default is the (historically grown) previous method:
                                                                 ! - zh, zv, zvh, ah linear in log-log space
                                                                 ! - rrhv, irhv, kdp, adp cubic in lin-lin space
  
  ! Country flag for the available observational radar data and/or
  !  desired background radar meta data list:
  !  0 = DWD German radar network       (default)
  !  1 = MeteoSwiss Swiss radar network
  INTEGER                 :: icountry

  ! .. Compute global composites of obs and sim reflectivity on the model grid
  !     for later use in the NWP-model, using data of one specific elevation of each station:
  ! ( EFFECTIVE ONLY IF NAMELIST PARAMETER ldo_bubbles=.TRUE. )
  LOGICAL                    ::   lsmooth_composite_bub_glob          ! If composite is smoothed by binomial filter
  INTEGER                    ::   eleind_for_composite_bub_glob,   &  ! elevation index to be used for compositing for warm bubble generator (comp_dbzsim_bub, comp_dbzobs_bub)
                                  nsmoothpoints_for_comp_bub_glob, &  ! width of symetric 2D binomial smoother in grid points
                                  nfilt_for_comp_bub_glob                 ! number of consecutive applications of the smoother

  ! ( EFFECTIVE ONLY IF NAMELIST PARAMETER ldo_composite=.TRUE. )
  LOGICAL                    ::   lsmooth_composite_glob          ! If composite is smoothed by binomial filter
  INTEGER                    ::   eleindlist_for_composite_glob(1:nel_composite_max), &    ! elevation index list to be used for compositing (comp_dbzsim, comp_dbzobs)
                                  levelidlist_for_composite_glob(1:nel_composite_max), &   ! level identifier of the composite in grib2-output
                                  nel_composite,  & ! Number of actual composites computed in the operator (must be <= nel_composite_max!)
                                  nsmoothpoints_for_comp_glob, &  ! width of symetric 2D binomial smoother in grid points
                                  nfilt_for_comp_glob             ! number of consecutive applications of the smoother

  ! (not yet a namelist parameter, but may become one in the future):
  REAL    (KIND=wpfwo)             ::   &
       htop            ! Maximum height MSL for radar computations ( <= model domain top height )
                       !  Has to be set in radar_interface.f90, get_model_config_for_radar()

  ! Meta data for attributes for writing cdfin files:
  ! Not yet a namelist parameter ...
  CHARACTER (len=cmaxlen)    :: cdfin_creator_model

  ! Approximate range resolution [m] (Global background value) to which the input observation data are aggregated in range on input and
  !  which is internally used for simulating synthetic radar volume data:
  !  (only effective if lreadmeta_from_netcdf=.true.)
  REAL    (KIND=wpfwo)             :: ra_inc_coarse_glob
  
  
  !==========================================================================================================
  !
  ! .. Type to hold the above namelist parameters (e.g., to store different sets for different model domains)
  !
  !==========================================================================================================

  TYPE glob_nml_type
    ! .. nradsta is not an active namelist parameter (because it is automatically determined
    !     by the program), but it also has to be stored and recovered for each domain:
    INTEGER                    :: nradsta
    INTEGER                    :: nbl_az
    REAL(kind=wpfwo)           :: htop
    ! .. and here are the "true" namelist parameters:
    INTEGER                    :: nradsta_namelist
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
    LOGICAL                    :: lmodfield_output
    LOGICAL                    :: loutnwp
    LOGICAL                    :: loutvolaux
  !
    LOGICAL                    :: lvoldata_output
    TYPE(t_voldata_ostream)    :: voldata_ostream(noutstreams_max)
    LOGICAL                    :: lfdbk_output
    LOGICAL                    :: lwrite_ready
    LOGICAL                    :: ldebug_radsim
    LOGICAL                    :: loutradwind
    LOGICAL                    :: lfill_vr_backgroundwind
    LOGICAL                    :: lout_geom
    LOGICAL                    :: loutdbz
    TYPE(t_dbzcalc_params)       :: dbz_meta_glob
    INTEGER                    :: itype_mpipar_lookupgen
    INTEGER                    :: pe_start_lookupgen
    INTEGER                    :: pe_end_lookupgen
    LOGICAL                    :: llookup_interp_mode_dualpol
    LOGICAL                    :: lcalc_dbz_on_radarbins
    LOGICAL                    :: loutpolstd
    LOGICAL                    :: loutpolall
    LOGICAL                    :: lextdbz
    LOGICAL                    :: lweightdbz
    LOGICAL                    :: lfall
    LOGICAL                    :: lonline
    LOGICAL                    :: lcomm_nonblocking_online
    LOGICAL                    :: lsode
    LOGICAL                    :: lsmooth
    LOGICAL                    :: lmds_z
    LOGICAL                    :: lmds_vr
    LOGICAL                    :: lreadmeta_from_netcdf
    LOGICAL                    :: lcheck_inputrecords
    LOGICAL                    :: lequal_azi_alldatasets
    LOGICAL                    :: lqc_flag
    LOGICAL                    :: ldealiase_vr_obs
    LOGICAL                    :: ltestpattern_hydrometeors
    LOGICAL                    :: ldo_composite
    LOGICAL                    :: lcomposite_output
    LOGICAL                    :: lcomposite_output_bub
    CHARACTER(LEN=12)          :: comp_grib2_packingtype
    CHARACTER(LEN=cmaxlen)     :: composite_file_pattern
    CHARACTER(LEN=cmaxlen)     :: composite_file_pattern_bub
    CHARACTER(LEN=cmaxlen)     :: ready_file_pattern
    LOGICAL                    :: ldo_bubbles
    LOGICAL                    :: ldo_bubbles_manual
    LOGICAL                    :: lcomm_nonblocking_output
    LOGICAL                    :: labort_if_problems_obsfiles
    REAL (kind=wpfwo)          :: dt_bubble_search
    REAL (kind=wpfwo)          :: t_offset_bubble_trigger_async
    REAL (kind=wpfwo)          :: prob_bubble
    REAL (kind=wpfwo)          :: maxdim_obs
    LOGICAL                    :: lbub_isolated
    REAL (kind=wpfwo)          :: threshold_obs(2)
    REAL (kind=wpfwo)          :: threshold_mod(2)
    REAL (kind=wpfwo)          :: areamin_mod(2)
    REAL (kind=wpfwo)          :: areamin_obs(2)
    REAL (kind=wpfwo)          :: mult_dist_obs
    REAL (kind=wpfwo)          :: mult_dist_mod
    REAL (kind=wpfwo)          :: add_dist_obs
    REAL (kind=wpfwo)          :: add_dist_mod
    REAL (kind=wpfwo)          :: dt_bubble_advect
    REAL (kind=wpfwo)          :: zlow_meanwind_bubble_advect
    REAL (kind=wpfwo)          :: zup_meanwind_bubble_advect
    CHARACTER(len=12)          :: bubble_type
    REAL (kind=wpfwo)          :: bubble_heatingrate
    REAL (kind=wpfwo)          :: bubble_timespan
    REAL (kind=wpfwo)          :: bubble_dT
    REAL (kind=wpfwo)          :: bubble_centz
    REAL (kind=wpfwo)          :: bubble_radx
    REAL (kind=wpfwo)          :: bubble_rady
    REAL (kind=wpfwo)          :: bubble_radz
    REAL (kind=wpfwo)          :: bubble_rotangle
    REAL (kind=wpfwo)          :: bubble_dT_noise
    LOGICAL                    :: bubble_holdrhconst
    LOGICAL                    :: bubble_addnoise_T
    INTEGER                    :: ngpsm_h_glob
    INTEGER                    :: ngpsm_v_glob
    INTEGER                    :: itype_supobing
    INTEGER                    :: supob_nrb
    INTEGER                    :: itype_obserr_vr
    INTEGER                    :: itype_metric_refl_fdbk
    REAL (KIND=wpfwo)          :: supob_azi_maxsector_vr
    REAL (KIND=wpfwo)          :: supob_cart_resolution
    REAL (KIND=wpfwo)          :: supob_ave_width
    REAL (KIND=wpfwo)          :: supob_minrange_vr
    REAL (KIND=wpfwo)          :: supob_minrange_z
    REAL (KIND=wpfwo)          :: supob_vrw
    REAL (KIND=wpfwo)          :: supob_rfl
    REAL (KIND=wpfwo)          :: supob_lowthresh_z_obs
    REAL (KIND=wpfwo)          :: supob_lowthresh_z_sim
    REAL (KIND=wpfwo)          :: ramp_lowdbz_obserr_vr
    REAL (KIND=wpfwo)          :: ramp_highdbz_obserr_vr
    REAL (KIND=wpfwo)          :: maxval_obserr_vr
    REAL (KIND=wpfwo)          :: baseval_obserr_vr
    REAL (KIND=wpfwo)          :: minval_obserr_lwc
    INTEGER                    :: thin_step_azi
    INTEGER                    :: thin_step_range
    INTEGER                    :: thin_step_ele
    INTEGER                    :: ind_ele_fdbk_glob(nel_max)
    INTEGER                    :: ind_ele_voldata_glob(nel_max)
    REAL (KIND=wpfwo)          :: obs_times_fdbk_glob(nobstimes_max)
    REAL (KIND=wpfwo)          :: dt_obs_fdbk_glob(3)
    REAL (KIND=wpfwo)          :: obs_times_voldata_glob(nobstimes_max)
    REAL (KIND=wpfwo)          :: dt_obs_voldata_glob(3)
    REAL (KIND=wpfwo)          :: ra_inc_coarse_glob
    CHARACTER(LEN=cmaxlen)     :: ydirradarin
    CHARACTER(LEN=cmaxlen)     :: ydirlistfile
    CHARACTER(LEN=cmaxlen)     :: ydirradarout
    CHARACTER(LEN=cmaxlen)     :: ysubdirfof
    CHARACTER(LEN=cmaxlen)     :: ysubdircomp
    CHARACTER(LEN=cmaxlen)     :: ydir_mielookup_read
    CHARACTER(LEN=cmaxlen)     :: ydir_mielookup_write
    CHARACTER(LEN=cmaxlen)     :: ydir_ready_write
    INTEGER                    :: icountry
    LOGICAL                    :: lsmooth_composite_bub_glob
    INTEGER                    :: eleind_for_composite_bub_glob
    INTEGER                    :: nsmoothpoints_for_comp_bub_glob
    INTEGER                    :: nfilt_for_comp_bub_glob
    LOGICAL                    :: lsmooth_composite_glob
    INTEGER                    :: eleindlist_for_composite_glob(1:nel_composite_max)
    INTEGER                    :: levelidlist_for_composite_glob(1:nel_composite_max)
    INTEGER                    :: nel_composite
    INTEGER                    :: nsmoothpoints_for_comp_glob
    INTEGER                    :: nfilt_for_comp_glob
  END TYPE glob_nml_type

  TYPE(glob_nml_type), ALLOCATABLE, DIMENSION(:), TARGET :: glob_nml_container
  TYPE(glob_nml_type), POINTER                           :: glob_nml

CONTAINS

  SUBROUTINE prep_domains_radar_nml ( )

    IMPLICIT NONE

    IF (ndoms_radar > ndoms_max) THEN
      WRITE (*,'(a,i3,a,i3,a)') 'ERROR prep_domains_radar_nml(): ndoms_radar = ', ndoms_radar, &
           ' exceeds allowed ndoms_max = ', ndoms_max, '!'
      STOP
    END IF

    ALLOCATE(glob_nml_container(ndoms_radar))

  END SUBROUTINE prep_domains_radar_nml

  SUBROUTINE crosscheck_domains_radar_nml ( )

    IMPLICIT NONE

    INTEGER                             :: i, j

    IF (ndoms_radar > ndoms_max) THEN
      WRITE (*,'(a,i3,a,i3,a)') 'ERROR crosscheck_domains_radar_nml(): ndoms_radar = ', ndoms_radar, &
           ' exceeds allowed ndoms_max = ', ndoms_max, '!'
      STOP
    END IF

    ! .. Some cross-checks among the global namelist parameters for different domains:
    DO i=1, ndoms_radar
      DO j=1, ndoms_radar
        IF (i /= j) THEN

          ! .. It is forbidden that EMVORADO writes his output for different domains into the same output directory:
          IF (TRIM(glob_nml_container(i)%ydirradarout) == TRIM(glob_nml_container(j)%ydirradarout)) THEN
            WRITE (*,'(a,/,a,/,a)') 'ERROR crosscheck_domains_radar_nml: same radar output directory (ydirradarout) ', &
                 TRIM(glob_nml_container(i)%ydirradarout), &
                 'chosen for different domains, but must be different for each domain!'
            STOP
          END IF

        END IF
      END DO
    END DO

  END SUBROUTINE crosscheck_domains_radar_nml

  SUBROUTINE store_domain_radar_nml ( idom_model )
    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom_model

    glob_nml   => glob_nml_container(list_domains_for_radar(idom_model))

    ! .. store all above namelist parameters in the container structure:

    glob_nml%nradsta                          = nradsta
    glob_nml%nbl_az                           = nbl_az
    glob_nml%htop                             = htop

    glob_nml%nradsta_namelist                 = nradsta_namelist
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
    glob_nml%lmodfield_output                 = lmodfield_output
    glob_nml%loutnwp                          = loutnwp
    glob_nml%loutvolaux                       = loutvolaux
  !
    glob_nml%lvoldata_output                  = lvoldata_output
    glob_nml%voldata_ostream                  = voldata_ostream
    glob_nml%lfdbk_output                     = lfdbk_output
    glob_nml%lwrite_ready                     = lwrite_ready
    glob_nml%ldebug_radsim                    = ldebug_radsim
    glob_nml%loutradwind                      = loutradwind
    glob_nml%lfill_vr_backgroundwind          = lfill_vr_backgroundwind
    glob_nml%lout_geom                        = lout_geom
    glob_nml%loutdbz                          = loutdbz
    glob_nml%dbz_meta_glob                    = dbz_meta_glob
    glob_nml%itype_mpipar_lookupgen           = itype_mpipar_lookupgen
    glob_nml%pe_start_lookupgen               = pe_start_lookupgen
    glob_nml%pe_end_lookupgen                 = pe_end_lookupgen
    glob_nml%llookup_interp_mode_dualpol      = llookup_interp_mode_dualpol
    glob_nml%lcalc_dbz_on_radarbins           = lcalc_dbz_on_radarbins
    glob_nml%loutpolstd                       = loutpolstd
    glob_nml%loutpolall                       = loutpolall
    glob_nml%lextdbz                          = lextdbz
    glob_nml%lweightdbz                       = lweightdbz
    glob_nml%lfall                            = lfall
    glob_nml%lonline                          = lonline
    glob_nml%lcomm_nonblocking_online         = lcomm_nonblocking_online
    glob_nml%lsode                            = lsode
    glob_nml%lsmooth                          = lsmooth
    glob_nml%lmds_z                           = lmds_z
    glob_nml%lmds_vr                          = lmds_vr
    glob_nml%lreadmeta_from_netcdf            = lreadmeta_from_netcdf
    glob_nml%lcheck_inputrecords              = lcheck_inputrecords
    glob_nml%lequal_azi_alldatasets           = lequal_azi_alldatasets
    glob_nml%lqc_flag                         = lqc_flag
    glob_nml%ldealiase_vr_obs                 = ldealiase_vr_obs
    glob_nml%ltestpattern_hydrometeors        = ltestpattern_hydrometeors
    glob_nml%ldo_composite                    = ldo_composite
    glob_nml%lcomposite_output                = lcomposite_output
    glob_nml%lcomposite_output_bub            = lcomposite_output_bub
    glob_nml%comp_grib2_packingtype           = comp_grib2_packingtype
    glob_nml%composite_file_pattern           = composite_file_pattern
    glob_nml%composite_file_pattern_bub       = composite_file_pattern_bub
    glob_nml%ready_file_pattern               = ready_file_pattern
    glob_nml%ldo_bubbles                      = ldo_bubbles
    glob_nml%ldo_bubbles_manual               = ldo_bubbles_manual
    glob_nml%lcomm_nonblocking_output         = lcomm_nonblocking_output
    glob_nml%labort_if_problems_obsfiles      = labort_if_problems_obsfiles
    glob_nml%dt_bubble_search                 = dt_bubble_search
    glob_nml%t_offset_bubble_trigger_async    = t_offset_bubble_trigger_async
    glob_nml%prob_bubble                      = prob_bubble
    glob_nml%maxdim_obs                       = maxdim_obs
    glob_nml%lbub_isolated                    = lbub_isolated
    glob_nml%threshold_obs                    = threshold_obs
    glob_nml%threshold_mod                    = threshold_mod
    glob_nml%areamin_mod                      = areamin_mod
    glob_nml%areamin_obs                      = areamin_obs
    glob_nml%mult_dist_obs                    = mult_dist_obs
    glob_nml%mult_dist_mod                    = mult_dist_mod
    glob_nml%add_dist_obs                     = add_dist_obs
    glob_nml%add_dist_mod                     = add_dist_mod
    glob_nml%dt_bubble_advect                 = dt_bubble_advect
    glob_nml%zlow_meanwind_bubble_advect      = zlow_meanwind_bubble_advect
    glob_nml%zup_meanwind_bubble_advect       = zup_meanwind_bubble_advect
    glob_nml%bubble_type                      = bubble_type
    glob_nml%bubble_heatingrate               = bubble_heatingrate
    glob_nml%bubble_timespan                  = bubble_timespan
    glob_nml%bubble_dT                        = bubble_dT
    glob_nml%bubble_centz                     = bubble_centz
    glob_nml%bubble_radx                      = bubble_radx
    glob_nml%bubble_rady                      = bubble_rady
    glob_nml%bubble_radz                      = bubble_radz
    glob_nml%bubble_rotangle                  = bubble_rotangle
    glob_nml%bubble_dT_noise                  = bubble_dT_noise
    glob_nml%bubble_holdrhconst               = bubble_holdrhconst
    glob_nml%bubble_addnoise_T                = bubble_addnoise_T
    glob_nml%ngpsm_h_glob                     = ngpsm_h_glob
    glob_nml%ngpsm_v_glob                     = ngpsm_v_glob
    glob_nml%itype_supobing                   = itype_supobing
    glob_nml%itype_obserr_vr                  = itype_obserr_vr
    glob_nml%itype_metric_refl_fdbk           = itype_metric_refl_fdbk
    glob_nml%supob_nrb                        = supob_nrb
    glob_nml%supob_azi_maxsector_vr           = supob_azi_maxsector_vr
    glob_nml%supob_cart_resolution            = supob_cart_resolution
    glob_nml%supob_ave_width                  = supob_ave_width
    glob_nml%supob_minrange_vr                = supob_minrange_vr
    glob_nml%supob_minrange_z                 = supob_minrange_z
    glob_nml%supob_vrw                        = supob_vrw
    glob_nml%supob_rfl                        = supob_rfl
    glob_nml%supob_lowthresh_z_obs            = supob_lowthresh_z_obs
    glob_nml%supob_lowthresh_z_sim            = supob_lowthresh_z_sim
    glob_nml%ramp_lowdbz_obserr_vr            = ramp_lowdbz_obserr_vr
    glob_nml%ramp_highdbz_obserr_vr           = ramp_highdbz_obserr_vr
    glob_nml%maxval_obserr_vr                 = maxval_obserr_vr
    glob_nml%baseval_obserr_vr                = baseval_obserr_vr
    glob_nml%minval_obserr_lwc                = minval_obserr_lwc
    glob_nml%thin_step_azi                    = thin_step_azi
    glob_nml%thin_step_range                  = thin_step_range
    glob_nml%thin_step_ele                    = thin_step_ele
    glob_nml%ind_ele_fdbk_glob                = ind_ele_fdbk_glob
    glob_nml%ind_ele_voldata_glob             = ind_ele_voldata_glob
    glob_nml%obs_times_fdbk_glob              = obs_times_fdbk_glob
    glob_nml%dt_obs_fdbk_glob                 = dt_obs_fdbk_glob
    glob_nml%obs_times_voldata_glob           = obs_times_voldata_glob
    glob_nml%dt_obs_voldata_glob              = dt_obs_voldata_glob
    glob_nml%ra_inc_coarse_glob               = ra_inc_coarse_glob
    glob_nml%ydirradarin                      = ydirradarin
    glob_nml%ydirlistfile                     = ydirlistfile
    glob_nml%ydirradarout                     = ydirradarout
    glob_nml%ysubdirfof                       = ysubdirfof
    glob_nml%ysubdircomp                      = ysubdircomp
    glob_nml%ydir_mielookup_read              = ydir_mielookup_read
    glob_nml%ydir_mielookup_write             = ydir_mielookup_write
    glob_nml%ydir_ready_write                 = ydir_ready_write
    glob_nml%icountry                         = icountry
    glob_nml%lsmooth_composite_bub_glob       = lsmooth_composite_bub_glob
    glob_nml%eleind_for_composite_bub_glob    = eleind_for_composite_bub_glob
    glob_nml%nsmoothpoints_for_comp_bub_glob  = nsmoothpoints_for_comp_bub_glob
    glob_nml%nfilt_for_comp_bub_glob          = nfilt_for_comp_bub_glob
    glob_nml%lsmooth_composite_glob           = lsmooth_composite_glob
    glob_nml%eleindlist_for_composite_glob    = eleindlist_for_composite_glob
    glob_nml%levelidlist_for_composite_glob   = levelidlist_for_composite_glob
    glob_nml%nel_composite                    = nel_composite
    glob_nml%nsmoothpoints_for_comp_glob      = nsmoothpoints_for_comp_glob
    glob_nml%nfilt_for_comp_glob              = nfilt_for_comp_glob



  END SUBROUTINE store_domain_radar_nml

  SUBROUTINE switch_to_domain_radar_nml ( idom_model )
    IMPLICIT NONE

    INTEGER, INTENT(in) :: idom_model

    glob_nml   => glob_nml_container(list_domains_for_radar(idom_model))

    ! .. restore all above namelist parameters from the container structure:

    nradsta                             =   glob_nml%nradsta
    nbl_az                              =   glob_nml%nbl_az
    htop                                =   glob_nml%htop

    nradsta_namelist                    =   glob_nml%nradsta_namelist
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
    lmodfield_output                    =   glob_nml%lmodfield_output
    loutnwp                             =   glob_nml%loutnwp
    loutvolaux                          =   glob_nml%loutvolaux
  !
    lvoldata_output                     =   glob_nml%lvoldata_output
    voldata_ostream                     =   glob_nml%voldata_ostream
    lfdbk_output                        =   glob_nml%lfdbk_output
    lwrite_ready                        =   glob_nml%lwrite_ready
    ldebug_radsim                       =   glob_nml%ldebug_radsim
    loutradwind                         =   glob_nml%loutradwind
    lfill_vr_backgroundwind             =   glob_nml%lfill_vr_backgroundwind
    lout_geom                           =   glob_nml%lout_geom
    loutdbz                             =   glob_nml%loutdbz
    dbz_meta_glob                       =   glob_nml%dbz_meta_glob
    itype_mpipar_lookupgen              =   glob_nml%itype_mpipar_lookupgen
    pe_start_lookupgen                  =   glob_nml%pe_start_lookupgen
    pe_end_lookupgen                    =   glob_nml%pe_end_lookupgen
    llookup_interp_mode_dualpol         =   glob_nml%llookup_interp_mode_dualpol
    lcalc_dbz_on_radarbins              =   glob_nml%lcalc_dbz_on_radarbins
    loutpolstd                          =   glob_nml%loutpolstd
    loutpolall                          =   glob_nml%loutpolall
    lextdbz                             =   glob_nml%lextdbz
    lweightdbz                          =   glob_nml%lweightdbz
    lfall                               =   glob_nml%lfall
    lonline                             =   glob_nml%lonline
    lcomm_nonblocking_online            =   glob_nml%lcomm_nonblocking_online
    lsode                               =   glob_nml%lsode
    lsmooth                             =   glob_nml%lsmooth
    lmds_z                              =   glob_nml%lmds_z
    lmds_vr                             =   glob_nml%lmds_vr
    lreadmeta_from_netcdf               =   glob_nml%lreadmeta_from_netcdf
    lcheck_inputrecords                 =   glob_nml%lcheck_inputrecords
    lequal_azi_alldatasets              =   glob_nml%lequal_azi_alldatasets
    lqc_flag                            =   glob_nml%lqc_flag
    ldealiase_vr_obs                    =   glob_nml%ldealiase_vr_obs
    ltestpattern_hydrometeors           =   glob_nml%ltestpattern_hydrometeors
    ldo_composite                       =   glob_nml%ldo_composite
    lcomposite_output                   =   glob_nml%lcomposite_output
    lcomposite_output_bub               =   glob_nml%lcomposite_output_bub
    comp_grib2_packingtype              =   glob_nml%comp_grib2_packingtype
    composite_file_pattern              =   glob_nml%composite_file_pattern
    composite_file_pattern_bub          =   glob_nml%composite_file_pattern_bub
    ready_file_pattern                  =   glob_nml%ready_file_pattern
    ldo_bubbles                         =   glob_nml%ldo_bubbles
    ldo_bubbles_manual                  =   glob_nml%ldo_bubbles_manual
    lcomm_nonblocking_output            =   glob_nml%lcomm_nonblocking_output
    labort_if_problems_obsfiles         =   glob_nml%labort_if_problems_obsfiles
    dt_bubble_search                    =   glob_nml%dt_bubble_search
    t_offset_bubble_trigger_async       =   glob_nml%t_offset_bubble_trigger_async
    prob_bubble                         =   glob_nml%prob_bubble
    maxdim_obs                          =   glob_nml%maxdim_obs
    lbub_isolated                       =   glob_nml%lbub_isolated
    threshold_obs                       =   glob_nml%threshold_obs
    threshold_mod                       =   glob_nml%threshold_mod
    areamin_mod                         =   glob_nml%areamin_mod
    areamin_obs                         =   glob_nml%areamin_obs
    mult_dist_obs                       =   glob_nml%mult_dist_obs
    mult_dist_mod                       =   glob_nml%mult_dist_mod
    add_dist_obs                        =   glob_nml%add_dist_obs
    add_dist_mod                        =   glob_nml%add_dist_mod
    dt_bubble_advect                    =   glob_nml%dt_bubble_advect
    zlow_meanwind_bubble_advect         =   glob_nml%zlow_meanwind_bubble_advect
    zup_meanwind_bubble_advect          =   glob_nml%zup_meanwind_bubble_advect
    bubble_type                         =   glob_nml%bubble_type
    bubble_heatingrate                  =   glob_nml%bubble_heatingrate
    bubble_timespan                     =   glob_nml%bubble_timespan
    bubble_dT                           =   glob_nml%bubble_dT
    bubble_centz                        =   glob_nml%bubble_centz
    bubble_radx                         =   glob_nml%bubble_radx
    bubble_rady                         =   glob_nml%bubble_rady
    bubble_radz                         =   glob_nml%bubble_radz
    bubble_rotangle                     =   glob_nml%bubble_rotangle
    bubble_dT_noise                     =   glob_nml%bubble_dT_noise
    bubble_holdrhconst                  =   glob_nml%bubble_holdrhconst
    bubble_addnoise_T                   =   glob_nml%bubble_addnoise_T
    ngpsm_h_glob                        =   glob_nml%ngpsm_h_glob
    ngpsm_v_glob                        =   glob_nml%ngpsm_v_glob
    itype_supobing                      =   glob_nml%itype_supobing
    itype_obserr_vr                     =   glob_nml%itype_obserr_vr
    itype_metric_refl_fdbk              =   glob_nml%itype_metric_refl_fdbk
    supob_nrb                           =   glob_nml%supob_nrb
    supob_azi_maxsector_vr              =   glob_nml%supob_azi_maxsector_vr
    supob_cart_resolution               =   glob_nml%supob_cart_resolution
    supob_ave_width                     =   glob_nml%supob_ave_width
    supob_minrange_vr                   =   glob_nml%supob_minrange_vr
    supob_minrange_z                    =   glob_nml%supob_minrange_z
    supob_vrw                           =   glob_nml%supob_vrw
    supob_rfl                           =   glob_nml%supob_rfl
    supob_lowthresh_z_obs               =   glob_nml%supob_lowthresh_z_obs
    supob_lowthresh_z_sim               =   glob_nml%supob_lowthresh_z_sim
    ramp_lowdbz_obserr_vr               =   glob_nml%ramp_lowdbz_obserr_vr
    ramp_highdbz_obserr_vr              =   glob_nml%ramp_highdbz_obserr_vr
    maxval_obserr_vr                    =   glob_nml%maxval_obserr_vr
    baseval_obserr_vr                   =   glob_nml%baseval_obserr_vr
    minval_obserr_lwc                   =   glob_nml%minval_obserr_lwc
    thin_step_azi                       =   glob_nml%thin_step_azi
    thin_step_range                     =   glob_nml%thin_step_range
    thin_step_ele                       =   glob_nml%thin_step_ele
    ind_ele_fdbk_glob                   =   glob_nml%ind_ele_fdbk_glob
    ind_ele_voldata_glob                =   glob_nml%ind_ele_voldata_glob
    obs_times_fdbk_glob                 =   glob_nml%obs_times_fdbk_glob
    dt_obs_fdbk_glob                    =   glob_nml%dt_obs_fdbk_glob
    obs_times_voldata_glob              =   glob_nml%obs_times_voldata_glob
    dt_obs_voldata_glob                 =   glob_nml%dt_obs_voldata_glob
    ra_inc_coarse_glob                  =   glob_nml%ra_inc_coarse_glob
    ydirradarin                         =   glob_nml%ydirradarin
    ydirlistfile                        =   glob_nml%ydirlistfile
    ydirradarout                        =   glob_nml%ydirradarout
    ysubdirfof                          =   glob_nml%ysubdirfof
    ysubdircomp                         =   glob_nml%ysubdircomp
    ydir_mielookup_read                 =   glob_nml%ydir_mielookup_read
    ydir_mielookup_write                =   glob_nml%ydir_mielookup_write
    ydir_ready_write                    =   glob_nml%ydir_ready_write
    icountry                            =   glob_nml%icountry
    lsmooth_composite_bub_glob          =   glob_nml%lsmooth_composite_bub_glob
    eleind_for_composite_bub_glob       =   glob_nml%eleind_for_composite_bub_glob
    nsmoothpoints_for_comp_bub_glob     =   glob_nml%nsmoothpoints_for_comp_bub_glob
    nfilt_for_comp_bub_glob             =   glob_nml%nfilt_for_comp_bub_glob
    lsmooth_composite_glob              =   glob_nml%lsmooth_composite_glob
    eleindlist_for_composite_glob       =   glob_nml%eleindlist_for_composite_glob
    levelidlist_for_composite_glob      =   glob_nml%levelidlist_for_composite_glob
    nel_composite                       =   glob_nml%nel_composite
    nsmoothpoints_for_comp_glob         =   glob_nml%nsmoothpoints_for_comp_glob
    nfilt_for_comp_glob                 =   glob_nml%nfilt_for_comp_glob

  END SUBROUTINE switch_to_domain_radar_nml


END MODULE radar_data_namelist
