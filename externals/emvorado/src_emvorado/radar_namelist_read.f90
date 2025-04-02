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


MODULE radar_namelist_read

!------------------------------------------------------------------------------
!
! Description:
!   This module provides methods to read the namelist parameters
!   of the radar forward operator EMVORADO and to print their values as
!   diagnostic output to some files.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

  USE radar_kind, ONLY : dp

  USE radar_dbzcalc_params_type, ONLY : t_dbzcalc_params, dbz_namlst_d

  USE radar_data, ONLY : &
       unused_value, miss_value, missval_int, miss_threshold, eps_dtobs, obsfile_missingname, &
       nradsta, &
       radar_meta_type, rs_meta, rsm_init_strings_blanks, &
       dbz_meta, &
       composite_meta_type, comp_meta, comp_meta_bub, mpi_comp_meta_typ, mpi_voldata_ostream_typ, &
       nradsta_max, cmaxlen, noutput_fields_max, nel_max, nobstimes_max, &
       nel_composite_max, ncountry_max, &
       my_radar_id, num_radar, num_radario, icomm_radar, lcompute_pe_fwo, num_compute_fwo, ydate_ini_mod, &
       mpi_dbzcalc_params_typ, mpi_radar_meta_typ_alltimes, &
       ndatakind, cndatakind, r_earth_dp, &
       i_fakeobs, i_dwd, &
       nradsta_all

  USE radar_data_namelist, ONLY:       &
       radarnmlfile, radarnmlfilepath, &
       nmlctrlfile,  nmlctrlfilepath,  &
       ldebug_radsim, loutradwind, &
       loutdbz, lcalc_dbz_on_radarbins, loutpolstd, loutpolall, lextdbz, &
       dbz_meta_glob, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       lcheck_inputrecords, lvoldata_output, t_voldata_ostream, noutstreams_max, voldata_ostream, &
       lfdbk_output, lwrite_ready, ready_file_pattern, &
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
       lmodfield_output, loutnwp, loutvolaux, &
  !
       nradsta_namelist, ydirradarin, ydirradarout, ysubdirfof, ysubdircomp, lequal_azi_alldatasets, &
       ydir_mielookup_read, ydir_mielookup_write, ydir_ready_write, icountry, itype_mpipar_lookupgen, &
       pe_start_lookupgen, pe_end_lookupgen, lqc_flag, llookup_interp_mode_dualpol, ldealiase_vr_obs, &
       itype_supobing, supob_nrb, &
       thin_step_azi, thin_step_range, thin_step_ele, ind_ele_fdbk_glob, ind_ele_voldata_glob, &
       supob_azi_maxsector_vr, supob_cart_resolution, &
       supob_ave_width, supob_minrange_vr, supob_minrange_z, supob_vrw, supob_rfl, &
       ngpsm_h_glob, ngpsm_v_glob, lmds_z, lmds_vr, &
       itype_obserr_vr, baseval_obserr_vr, maxval_obserr_vr, ramp_lowdbz_obserr_vr, ramp_highdbz_obserr_vr, &
       supob_lowthresh_z_obs, supob_lowthresh_z_sim, ltestpattern_hydrometeors, &
       ldo_composite, lcomposite_output, composite_file_pattern, obs_times_fdbk_glob, dt_obs_fdbk_glob, &
       obs_times_voldata_glob, dt_obs_voldata_glob, ra_inc_coarse_glob, &
       lfill_vr_backgroundwind, &
       lsmooth_composite_bub_glob, eleind_for_composite_bub_glob, &
       nsmoothpoints_for_comp_bub_glob, nfilt_for_comp_bub_glob, &
       lsmooth_composite_glob, eleindlist_for_composite_glob, levelidlist_for_composite_glob, &
       nsmoothpoints_for_comp_glob, nfilt_for_comp_glob, nel_composite, &
       ldo_bubbles, ldo_bubbles_manual, ydirlistfile, lcomposite_output_bub, composite_file_pattern_bub, comp_grib2_packingtype, &
       dt_bubble_search, tstart_bubble_search, tend_bubble_search, &
       prob_bubble, maxdim_obs, lbub_isolated, &
       t_offset_bubble_trigger_async, &
       threshold_obs, threshold_mod, areamin_mod, areamin_obs, &
       mult_dist_obs, add_dist_obs, mult_dist_mod, add_dist_mod, &
       dt_bubble_advect, zlow_meanwind_bubble_advect, zup_meanwind_bubble_advect, &
       bubble_type, bubble_heatingrate, bubble_timespan, bubble_dT, &
       bubble_centz, bubble_radx, bubble_rady, bubble_radz, bubble_rotangle, &
       bubble_dT_noise, bubble_holdrhconst, bubble_addnoise_T, &
       htop, lcomm_nonblocking_output, lcomm_nonblocking_online, &
       cdfin_creator_model, labort_if_problems_obsfiles, itype_metric_refl_fdbk, minval_obserr_lwc, &
       labort_if_problems_gribout

  USE radar_obs_meta_list, ONLY : set_scanname, get_meta_proto
#ifdef NETCDF
  USE radar_obs_meta_list, ONLY : get_meta_network_all
  USE radar_obs_meta_read, ONLY : read_meta_info_all
#endif

  USE radar_utilities, ONLY : new_datetime, diff_seconds, get_alpha3_eff_0, get_free_funit, &
                              split_string, tolower, toupper

  USE radar_parallel_utilities, ONLY : distribute_values_radar, distribute_path_radar, global_values_radar

  USE radar_interface, ONLY : abort_run, hhl, get_model_inputdir, &
       num_regular_obstimes, &
       check_obstime_within_modelrun, get_domain_starttime_in_sec, get_domain_endtime_in_sec, &
       get_model_outputdir, bottomlevel_stag, &
       geo2model_coord_domaincheck, geo2model_cellindex, &
       interp2d_model2geo_horiz_scalar, get_model_top_height, &
       get_composite_metadata, run_is_restart

#ifdef __ICON__
  USE radar_interface, ONLY : setup_auxgrid_for_cellindex
#endif

#ifndef NOMPI
  USE mpi
#endif

  !================================================================================
  !================================================================================

  IMPLICIT NONE

  PRIVATE

  !==============================================================================
  ! Public Subroutines:

  PUBLIC :: input_radarnamelist, &
            ctrl_output_ista_nuspec_fwo, ctrl_output_dbzmeta_nuspec_fwo

  !==============================================================================
  ! Local variables:

  ! Buffers for MPI distributing namelists:
  INTEGER, PARAMETER :: ibuflen = 200
  INTEGER            :: intbuf  (ibuflen)
  REAL(KIND=dp)      :: realbuf (ibuflen)
  LOGICAL            :: logbuf  (ibuflen)
  CHARACTER(LEN=100) :: charbuf (ibuflen)

  REAL(kind=dp), PARAMETER           :: eps = 1e-20_dp

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !=========================================================================
  !
  ! Subroutine for input of the namelist RADARSIM_PARAMS
  !
  !=========================================================================

  SUBROUTINE input_radarnamelist (izdom)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Description: Opens and reads the namelist on proc 0 for domain izdom, distributes the
    !              values/structs to all other nodes, does some
    !              cross-checks and computes some derived parameters within the
    !              structs rs_meta(ista) and dbz_meta(ista).
    !
    !              The global namelist parameters are stored on the global namelist variables
    !              which means, that the values from a possible previous call for a different
    !              domain are lost (overwritten). However, there is the possibility to
    !              store the parameters for a specific domain in a container by afterwards calling
    !              the subroutine
    !                  CALL store_domain_radar_nml ( idom )
    !              from module radar_data_namelist.f90.
    !              The set of global namelist parameters for a specific domain can be
    !              later restored to the global namelist variables by
    !                  CALL switch_to_domain_radar_nml ( idom )
    !              which has to be done in each call to organize_radar('compute', itimelevel, idom)
    !
    !              Should there be separate namelists for different domains, they should
    !              be in different files, ending in suffix "_DOM01", "_DOM02" etc.
    !              This is necessary because each namelist is actually read several
    !              times before all namelist parameters are properly defined.
    !
    ! Method:
    !
    ! Input files: Namelist-file INPUT_RADARSIM, containing namelist RADARSIM_PARAMS
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in/out):

    INTEGER, INTENT(in) :: izdom

    !------------------------------------------------------------------------------
    !


    !------------------------------------------------------------------------------
    ! .. Namelist parameter dom to identify the domain for which a
    !    namelist group RADARSIM_PARAMS from the INPUT_RADARSIM file belongs to.
    !    This is a mandatory namelist parameter and will be checked and matched
    !    in subroutine read_radarnamelist() against izdom:

    INTEGER :: dom

    !------------------------------------------------------------------------------
    ! .. Defaults for namelist parameters:

    TYPE(radar_meta_type)     :: rs_meta_d(1)
    TYPE(radar_meta_type), ALLOCATABLE :: rs_meta_ncdf(:)
    TYPE(t_dbzcalc_params)    :: dbz_meta_d(1)
    TYPE(composite_meta_type) :: comp_meta_d, comp_meta_bub_d
    TYPE(t_voldata_ostream)   :: voldata_ostream_d(noutstreams_max)
    TYPE(t_voldata_ostream)   :: voldata_ostream_save(noutstreams_max)
    INTEGER :: nradsta_namelist_d, nradsta_ncdf, nradsta_changed
    INTEGER :: ngpsm_h_glob_d      ! global setting for ngpsm_h value
    INTEGER :: ngpsm_v_glob_d      ! global setting for ngpsm_v value
    LOGICAL :: ldebug_radsim_d,  & ! turn on/off debug mode
               lvoldata_output_d,  & ! enable output of ASCII files
               lfdbk_output_d,   & ! enable output of NetCDF feedback files
               lwrite_ready_d,   & ! enable output of READY files after each EMVORADO output timestep
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
               lmodfield_output_d, & ! enable radar param output on model grid
               loutnwp_d,          & ! enable nwp field output on model grid
               loutvolaux_d,       & ! enable additional volume scan radar param output
  !
               lout_geom_d,      & ! enable output of radar bin heights and local elevations
               loutradwind_d,    & ! enable calculation and output of
                                   !   radial wind for each radar
               loutdbz_d,        & ! enable calculation and output of
                                   !   reflectivity for each radar
  
               lcalc_dbz_on_radarbins_d, & !.. if radar moments such as dbz or polar. param. should be computed on radar bins or not.
                                   !    If .true., the model state is first interpolated to radarbins, before moments are derived from them.
                                   !    Otherwise, moments are computed on model grid and then interpolated to radarbins (old default method).
               loutpolstd_d,     & ! enable calculation and output of
                                   !   standard polarimetric parameters (ZDR, KDP, RHV)
                                   !   for each radar
               loutpolall_d,     & ! enable calculation and output of
                                   !   additional polarimetric parameters (LDR)
                                   !   for each radar
               lextdbz_d,        & ! enable calculation and consideration of
                                   !   extinction effects on reflectivity and polarimetric
                                   !   parameters for each radar
               lweightdbz_d,&      ! enable to take the (attenuated)
                                   !   reflectivity as weight into account
                                   !   while caculating the radial wind for
                                   !   each radar
               lfall_d,          & ! enable to take the falling velocity of
                                   !   hydrometeors into account while
                                   !   caculating the radial wind for each radar
               lfill_vr_backgroundwind_d, & ! fill simulated missing radial winds with background wind from model U, V, W
               lonline_d,        & ! enable to change the strategy for
                                   !   computing of radio propagation path: static 4/3-earth ("offline") or
                                   !   dynamic ("online") depending on actual refractive index of air
               lcomm_nonblocking_online_d, & ! non-blocking (=more efficient, but more memory) communication for azim. slices
               lcomm_nonblocking_output_d, & ! non-blocking (=more efficient, but more memory) communication for radar-IO-PEs

               lsode_d,          & ! Mehode of SODE
               lsmooth_d,        & ! enable spatial smoothing (weighted spatial
                                   !   mean over measuring volume increasing with distance)
               lmds_z_d,         & ! enable minimum detectable signal (quadratic function of range), i.e.,
               lmds_vr_d,        & !   if Z < Z0(r), then set zero_value for Ze and miss_value for vr
               lreadmeta_from_netcdf_d, & ! read meta data from netCDF files
               lcheck_inputrecords_d, &   ! check input records in netCDF files and set them individually
                                          !   for each input data type file
               lequal_azi_alldatasets_d,& ! indicate if all azimuts in the different files for z, vr, qz, qv are the same
                                          !   to save time during reading of DWD netCDF files
               labort_if_problems_obsfiles_d,& ! whether or not problems with obs data files and/or radar meta data
               labort_if_problems_gribout_d, &   !  and grib output files should be fatal or not
               lqc_flag_d           , & ! quality control flags evaluation for the obs
               ldealiase_vr_obs_d,    & ! dealiasing radial wind for the vr obs
               ltestpattern_hydrometeors_d,  &  ! Testmode for reflectivity calculation in the
                                                !   first time step for idealized cases / testsuite:
               ldo_composite_d, lcomposite_output_d, &
               ldo_bubbles_d,   &    ! If .true., automatically detect the need for warm bubbles from composites; requires ldo_composite=.TRUE. in radar namelist and works together with ldo_bubbles in the driving model
               ldo_bubbles_manual_d,&! If .true., allows manual specification of bubbles in real case runs (lartif_data=.false.) using the bubble parameters from the ARTIFDATA namelist
               lcomposite_output_bub_d, &
               lbub_isolated_d       ! Bubble generator: check that the bubbles are isolated from other objects in observations



    CHARACTER(len=cmaxlen)      :: ydirradarin_d, ydirlistfile_d, ydirradarout_d,         &
         &                         ysubdirfof_d, ysubdircomp_d, composite_file_pattern_d, &
         &                                                      composite_file_pattern_bub_d, &
         &                         ydir_mielookup_read_d, ydir_mielookup_write_d,         &
         &                         ydir_ready_write_d, ready_file_pattern_d
    !.. Parallelization strategy for lookup table generation:
    INTEGER                         :: itype_mpipar_lookupgen_d  ! 1 = parallelization over dbzmeta-sets and hydrometeors
                                                                 ! 2 = parallelization over table elements
    !.. Processor number range in icomm_compute_fwo for parallel lookup table generation for all itype_mpipar_lookupgen methods:
    INTEGER                         :: pe_start_lookupgen_d  ! start PE  [0...num_compute_fwo-1]
    INTEGER                         :: pe_end_lookupgen_d    ! end PE    [0...num_compute_fwo-1]
    LOGICAL                         :: llookup_interp_mode_dualpol_d  ! if to use exact same lookup table interpolation method
                                                                      ! for all radar moments (cubic in lin-lin-space).
                                                                      ! The default is the (historically grown) previous method:
                                                                      ! - zh, zv, zvh, ah linear in log-log space
                                                                      ! - rrhv, irhv, kdp, adp cubic in lin-lin space
    CHARACTER(LEN=12)       :: comp_grib2_packingtype_d

    INTEGER :: &
         icountry_d, & ! Country flag for the available observational radar data and/or
                       !  desired background radar meta data list:
                       !  1 = DWD German radar network       (default)
                       !  2 = MeteoSwiss Swiss radar network
                       !  3 = Italy (Emilia-Romagna network)
         itype_supobing_d , &  ! option for superobing: 0: no superobing; 1: averaging;  2: median
         itype_obserr_vr_d, &  ! option for observation error for (superobe'd or original) radial wind (fdbk%e_o):
                               !   = 0: constant observation error (=baseval_obserr_vr)
                               !   = 1: e_o is a linear ramp as function of (superobe'd or original) reflectivity
                               !   = 2: superobe'd e_o + e_o-weighted superobe'd radial wind. For this,
                               !        e_o is computed from reflectivity before superobing (same linear ramp as for option 1) and is used:
                               !        a) as second weight (inverse e_o) additionally to spatial Cressman-weight in superobbing of radial wind
                               !        b) for superobbing of e_o itself with the same weighting
         itype_metric_refl_fdbk_d, & ! namelist setting to decide about the metric for reflectivity in the feedback files:
                                     !  = 1: write dBZ to feedback files
                                     !  = 2: convert to effective LWC = 0.004*Z^0.55, but leave obs error factor at 1.0
                                     !  = 3: convert to effective LWC as for (2) and set obs-dependent obs error factor
         supob_nrb_d      , &  ! lower threshold for number of radar bins used for superobbing
         thin_step_azi_d  , &  ! if no superobing: step width of data thinning in azi direction
         thin_step_range_d, &  ! if no superobing: step width of data thinning in range direction
         thin_step_ele_d,   &  ! always effective: step width of data thinning for elevations
         ind_ele_fdbk_glob_d(nel_max),    & ! array of indices which elevations are written into feedback file (global for all radar stations)
         ind_ele_voldata_glob_d(nel_max), & ! array of indices which elevations are written into volume data files (global for all radar stations)
         eleind_for_composite_bub_glob_d,      & ! elevation index to be used for compositing for warm bubble generator (comp_dbzsim_bub, comp_dbzobs_bub)
         eleindlist_for_composite_glob_d(nel_composite_max), & ! elevation index to be used for the other composites (comp_dbzsim, comp_dbzobs)
         levelidlist_for_composite_glob_d(nel_composite_max), & ! level identifier of the composite in grib2-output
         nel_composite_d

    REAL (KIND=dp)  :: &
         supob_azi_maxsector_vr_d , &  ! maximal azimut sector (symetrical to its center) for v_r superobing [deg]
         supob_cart_resolution_d  , &  ! resoltion of the cartesian grid for superobing [m]
         supob_ave_width_d        , &  ! width of averaging area for superobing [m]
         supob_minrange_vr_d      , &  ! min. range of superobing points for vr [m]
         supob_minrange_z_d       , &  ! min. range of superobing points for z  [m]
         supob_vrw_d              , &  ! threshold for stddev of radial wind values within a superobed bin [m/s]
         supob_rfl_d              , &  ! threshold for stddev of reflectivity values within a superobed bin [mm^6/m^3]
         supob_lowthresh_z_obs_d  , &  ! lower threshold for valid Ze values, to which no rain data are
         supob_lowthresh_z_sim_d  , &  !   set before average_superobing
                                       ! / sub-switches for itype_obserr_vr (see sketch down below):
         ramp_lowdbz_obserr_vr_d  , &  ! for options   1/2: lower ramp dBZ-threshold for increasing observation error in radial wind (miss_value = neutral value)
         ramp_highdbz_obserr_vr_d , &  ! for options   1/2: upper ramp dBZ-threshold for increasing observation error in radial wind (miss_value = neutral value)
         maxval_obserr_vr_d       , &  ! for options   1/2: obs error for dBZ-values <  ramp_highdbz_obserr_vr
         baseval_obserr_vr_d      , &  ! for options 0/1/2: obs error for dBZ-values >= ramp_highdbz_obserr_vr (options 1/2) or general value (option 0)
                                       !   in general: 1.0 means that the obserr for vr is generally defined as a relative error
         minval_obserr_lwc_d      , &  ! for itype_metric_refl_fdbk = 3: asymptotic obs error for LWC -> 0.0
         obs_times_fdbk_glob_d(nobstimes_max), & ! global list of obs times to write into the feedback files
         dt_obs_fdbk_glob_d(3)       , &  ! global value of triplet for dt_obs to construct obs times to write into the feedback files
         obs_times_voldata_glob_d(nobstimes_max), & ! global list of obs times to write into the volume data files
         dt_obs_voldata_glob_d(3)    , &  ! global value of triplet for dt_obs to construct obs times to write into the volume data files
         ra_inc_coarse_glob_d          , &  ! global value of ra_inc to definde the range resolution for coarsening of input data
         dt_bubble_search_d       , &  ! Bubble generator: time interval from one automatic bubble search to the next [seconds] (SYNCHRONIZE WITH COMPOSITE TIMES!!!)
         tstart_bubble_search_d, tend_bubble_search_d, &
         t_offset_bubble_trigger_async_d, & ! Bubble generator: in case of asynchr. radar IO for runtime optimization:
                                            !  Time delay of (advection corrected) bubble triggering after detection step.
                                            !  The worker PEs of the model may continue in parallel to the async IO PEs
                                            !  for this amount of model time, until they synchronize with the IO PEs to
                                            !  to gather information on new bubbles. An offset of 0.0 would be physically ideal but
                                            !  but destroys the runtime advantage of asynchr. IO.
         prob_bubble_d, &      ! Bubble generator: Probability of triggering a bubble when it is found [0-1]
         threshold_obs_d(2), & ! Bubble generator: thresholds for obs that define the minimum dBZ in an object- and the high intensity region [dBZ]
         threshold_mod_d(2), & ! Bubble generator: thresholds for sim that define the minimum dBZ in an object- and the high intensity region [dBZ]
         maxdim_obs_d, &       ! Bubble generator: maximum dimension of not founded cell that can trigger a bubble (very large objects are not targeted) [meters]
         areamin_mod_d(2), &   ! Bubble generator: minimum area in obs of the object- and high intensity region for detecting an object  [m^2]
         areamin_obs_d(2), &   ! Bubble generator: minimum area in sim of the object- and high intensity region for detecting an object  [m^2]
         mult_dist_obs_d, &    ! Bubble generator: multiplicative axis factor for the obs object. If it is high, bubbles are more difficult [-]
         mult_dist_mod_d, &    ! Bubble generator: multiplicative axis factor for the sim object. If it is high, bubbles are more difficult [-]
         add_dist_obs_d, &     ! Bubble generator: additive axis increase for the obs object. If it is high, bubbles are more difficult [meters]
         add_dist_mod_d, &     ! Bubble generator: additive axis increase for the sim object. If it is high, bubbles are more difficult [meters]
         dt_bubble_advect_d, & ! Bubble generator: time scale for downstream advection of automatic bubbles [seconds]
         zlow_meanwind_bubble_advect_d,& ! Bubble generator: the lower bound of averaging height interval for bubble advection speed [meters AMSL]
         zup_meanwind_bubble_advect_d, & ! Bubble generator: the upper bound of averaging height interval for bubble advection speed [meters AMSL]
         bubble_heatingrate_d, & ! Bubble generator: heating rate for the bubbles of type 'cos-hrd' [K/s]
         bubble_timespan_d, &    ! Bubble generator: timespan for heating the bubbles of type 'cos-hrd' [s]
         bubble_dT_d, &          ! Bubble generator: temperature disturbance for the bubbles of type 'cos-instant' [K]
         bubble_centz_d, &       ! Bubble generator: center height MSL (Z) of the bubbles [m]
         bubble_radx_d, &        ! Bubble generator: horizontal radius (main axis) in X-dir of the bubbles [m]
         bubble_rady_d, &        ! Bubble generator: horizontal radius (main axis) in Y-dir of the bubbles [m]
         bubble_radz_d, &        ! Bubble generator: vertical radius in Z-dir of the bubbles [m]
         bubble_rotangle_d, &    ! Bubble generator: Rotation angle of the main axes of bubbles [degrees]
         bubble_dT_noise_d       ! Bubble generator: In case of bubble_addnoise_T=.true., relative noise level, such that
                                 !  dT_bubble = dT_bubble * (1 + bub_dT_bubblenoise * random_noise[-1,1])
    LOGICAL :: &
         bubble_holdrhconst_d, & ! Bubble generator: Switch to choose if RH should kept constant during heating or not
         bubble_addnoise_T_d     ! Bubble generator: Switch to activate some random noise on the bubbles with a
                                 !  relative amplitude of "bubble_dT_noise"
    CHARACTER (LEN=12) :: bubble_type_d  ! Bubble generator: Type of the bubble ('cos-hrd', 'cos-instant')

    ! .. Namelist parameter to define all reflectivty related parameters
    !    in a global way for all radars.
    !      (IN PRINCIPLE, THIS COULD REPLACE THE EXISTING
    !       ITYPE_REFL_GLOB AND LLOOKUP_MIE_GLOB NAMELIST PARAMETERS,
    !       BUT THE LATTER ARE KEPT FOR COMPATIBILITY WITH OLDER NAMELIST VERSIONS.
    !       THE LATTER WILL TAKE PRECENDENCE OVER THEIR DBZ_GLOB% COUNTERPARTS,
    !       IF PRESENT IN THE NAMELIST)
    TYPE(t_dbzcalc_params)          :: dbz_meta_glob_d

    !------------------------------------------------------------------------------
    ! .. Local variables:

    CHARACTER (LEN=10)              :: cnobs
    CHARACTER (LEN=*), PARAMETER    :: yzroutine = 'emvorado::input_radarnamelist()'
    CHARACTER(len=cmaxlen)          :: errstring, cstation
    CHARACTER(len=2*cmaxlen)        :: nlradarsimfile, ctrlfile, io_errmsg
    LOGICAL                         :: nlfile_exists, found
    INTEGER                         :: ctrlfile_unit

    TYPE(radar_meta_type)           :: rs_meta_buf(1)
    TYPE(radar_meta_type), ALLOCATABLE :: rs_meta_for_country(:)
    TYPE(t_dbzcalc_params)          :: dbz_meta_buf(1)
    TYPE(t_dbzcalc_params)          :: dbz_meta_for_country(nradsta_max)
    TYPE(t_dbzcalc_params)          :: dbz_meta_ncdf(nradsta_max)
    TYPE(radar_meta_type)           :: rs_meta_from_nmlst(nradsta_max)
    INTEGER                         :: i, ii, it, iend, ista, &
                                       nuin, nibuf, nrbuf, nlbuf, &
                                       ierr, err_meta(nradsta_max), err_miss(nradsta_max), &
                                       ista_namelist(nradsta_max), &
                                       ista_ncdf(nradsta_max), nradsta_valid
    INTEGER                         :: iu, io, ju, jo, &
                                       nel_fdbk_glob, nel_fdbk_tmp, nobs_times_fdbk_glob, &
                                       nel_voldata_glob, nel_voldata_tmp, nobs_times_voldata_glob, &
                                       station_id_tmp
    REAL    (KIND=dp)               :: rlon_r, rlat_r, hsurf_r, wi, wj, &
                                       rlon_min, rlon_max, rlat_min, rlat_max, tmptime
    LOGICAL                         :: is_inside
    LOGICAL, SAVE                   :: lfirst_ctrl_output = .TRUE.

    CHARACTER(len=*), PARAMETER       :: format_char       = '(T8,A,T40,A,T55,A,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_path       = '(T8,A,T30,A,/,T30,A,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_int        = '(T8,A,T40,I12,T55,I12,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_int_vec    = '(T8,A,i3.3,A,T40,I12,T55,I12,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_real       = '(T8,A,T40,es12.5,T55,es12.5,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f     = '(T8,A,T40,f12.1,T55,f12.1,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f_vec = '(T8,A,i3.3,A,T40,f12.1,T55,f12.1,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f_vec2= '(T4,A,i2.2,A,T40,f12.1,T55,f12.1,T70,A)'
    CHARACTER(len=*), PARAMETER       :: format_logical    = '(T8,A,T40,L12,T55,L12,T70,A)'

    
    !------------------------------------------------------------------------------
    ! .. Definition of namelist:

    NAMELIST /RADARSIM_PARAMS/ dom, ldebug_radsim, loutradwind, &
                               loutdbz, lcalc_dbz_on_radarbins, loutpolstd, loutpolall, lextdbz, &
                               dbz_meta_glob, lout_geom, &
                               lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
                               lcheck_inputrecords, lvoldata_output, voldata_ostream, lfdbk_output, &
                               lwrite_ready, ready_file_pattern, &
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
                               lmodfield_output, loutnwp, loutvolaux, &
  !
                               nradsta_namelist, rs_meta, dbz_meta, &
                               ydirradarin, ydirradarout, ysubdirfof, ysubdircomp, lequal_azi_alldatasets, &
                               ydir_mielookup_read, ydir_mielookup_write, ydir_ready_write, icountry, itype_mpipar_lookupgen, &
                               pe_start_lookupgen, pe_end_lookupgen, llookup_interp_mode_dualpol, lqc_flag, ldealiase_vr_obs, &
                               itype_supobing, supob_nrb, &
                               thin_step_azi, thin_step_range, thin_step_ele, ind_ele_fdbk_glob, ind_ele_voldata_glob, &
                               supob_azi_maxsector_vr, supob_cart_resolution, &
                               supob_ave_width, supob_minrange_vr, supob_minrange_z, supob_vrw, supob_rfl, &
                               ngpsm_h_glob, ngpsm_v_glob, lmds_z, lmds_vr, &
                               itype_obserr_vr, baseval_obserr_vr, maxval_obserr_vr, &
                               ramp_lowdbz_obserr_vr, ramp_highdbz_obserr_vr, &
                               supob_lowthresh_z_obs, supob_lowthresh_z_sim, ltestpattern_hydrometeors, &
                               ldo_composite, obs_times_fdbk_glob, dt_obs_fdbk_glob, &
                               obs_times_voldata_glob, dt_obs_voldata_glob, ra_inc_coarse_glob, &
                               lfill_vr_backgroundwind, itype_metric_refl_fdbk, minval_obserr_lwc, &
                               eleind_for_composite_bub_glob, nel_composite, eleindlist_for_composite_glob, &
                               levelidlist_for_composite_glob, comp_meta_bub, &
                               ldo_bubbles, ldo_bubbles_manual, ydirlistfile, &
                               lcomposite_output_bub, composite_file_pattern_bub, &
                               lcomposite_output, composite_file_pattern, comp_grib2_packingtype, &
                               dt_bubble_search, tstart_bubble_search, tend_bubble_search, &
                               prob_bubble, maxdim_obs, lbub_isolated, &
                               threshold_obs, threshold_mod, areamin_mod, areamin_obs, &
                               mult_dist_obs, add_dist_obs, mult_dist_mod, add_dist_mod, &
                               dt_bubble_advect, zlow_meanwind_bubble_advect, zup_meanwind_bubble_advect, &
                               t_offset_bubble_trigger_async, &
                               bubble_type, bubble_heatingrate, bubble_timespan, bubble_dT, &
                               bubble_centz, bubble_radx, bubble_rady, bubble_radz, bubble_rotangle, &
                               bubble_dT_noise, bubble_holdrhconst, bubble_addnoise_T, &
                               lcomm_nonblocking_output, lcomm_nonblocking_online, comp_meta, labort_if_problems_obsfiles,&
                               labort_if_problems_gribout


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    IF (ldebug_radsim .OR. my_radar_id == 0) WRITE (*,'(a,a,i2,a,i5)') &
         TRIM(yzroutine), ' for domain ', izdom,' on proc ', my_radar_id

    ! name of namelist file (one file for different domain specific namelists):
    nlradarsimfile(:) = ' '
    IF (LEN_TRIM(radarnmlfilepath) > 0) THEN
      nlradarsimfile = TRIM(radarnmlfilepath) // '/' // TRIM(radarnmlfile)
    ELSE
      nlradarsimfile = TRIM(radarnmlfile)
    ENDIF

    !-------------------------------------------------------------------------------
    ! Section 1:  Set defaults for the namelist parameters
    !-------------------------------------------------------------------------------

    ! Not yet a namelist parameter ...
    cdfin_creator_model(:) = ' '
#ifdef __COSMO__
    cdfin_creator_model = 'COSMO+EMVORADO'
#endif
#ifdef __DACE__
    cdfin_creator_model = 'DACE+EMVORADO'
#endif
#ifdef __ICON__
    cdfin_creator_model = 'ICON+EMVORADO'
#endif

    ldebug_radsim_d = .FALSE.
    lout_geom_d     = .FALSE.
    loutdbz_d       = .TRUE.
      lcalc_dbz_on_radarbins_d = .FALSE.
      loutpolstd_d  = .FALSE.  ! only effective if dbz_meta(ista)%itype_refl > 4
      loutpolall_d  = .FALSE.  ! only effective if dbz_meta(ista)%itype_refl > 4
      lextdbz_d     = .FALSE.  ! only effective if dbz_meta(ista)%itype_refl = 1 or > 4
      lmds_z_d      = .FALSE.
    loutradwind_d   = .TRUE.
      lweightdbz_d  = .FALSE.
      lfall_d       = .FALSE.
      lfill_vr_backgroundwind_d = .FALSE.
      lmds_vr_d     = .FALSE.
    lonline_d       = .FALSE.
      lsode_d       = .FALSE.
      ! Non-blocking comm. for collecting interpolated model data on the azimuthal slices
      !  (=more efficient, but needs more memory) for online beam propagation
      !  (comm. from workers <-> workers):
      lcomm_nonblocking_online_d = .TRUE.
    lsmooth_d       = .FALSE.

    ! .. get the top height "htop" up to which radar computations should be performed,
    !     normally this is vertical model top height (MSL). Must be called on all
    !     PEs in icomm_radar. NOT YET A NAMELIST PARAMETER, BUT STORED FOR EACH DOMAIN
    !     IN THE NAMELIST CONTAINER STRUCTURE IN SUBROUTINE SWITCH_TO_DOMAIN_NML().
    CALL get_model_top_height (izdom, htop)

    ! For the composites that are output for different underlying elevations:
    ldo_composite_d = .FALSE.
      CALL get_composite_metadata ( izdom, comp_meta_d )     ! Default from radar_interface.f90
      lcomposite_output_d = .FALSE.
      composite_file_pattern_d(:) = ' '   ! empty pattern, leads to default file name
      nel_composite_d                    = 2
      eleindlist_for_composite_glob_d(:) = -99
      eleindlist_for_composite_glob_d(1:nel_composite) = (/ 1, 2 /)
      levelidlist_for_composite_glob_d(:) = -99
      comp_grib2_packingtype_d(:) = ' '
      comp_grib2_packingtype_d    = 'grid_ccsds'  ! 'grid_simple', 'grid_ccsds' (default), 'png'
      ! Not yet in the namelist regarding the composites:
      lsmooth_composite_bub_glob = .FALSE.    ! If composite is smoothed by binomial filter
      nsmoothpoints_for_comp_bub_glob = 9     ! width of symetric 2D binomial smoother in grid points
      nfilt_for_comp_bub_glob = 1             ! number of consecutive applications of the smoother

    ! For automatic warm bubble generator:
    ldo_bubbles_d   = .FALSE.
      CALL get_composite_metadata ( izdom, comp_meta_bub_d ) ! Default from radar_interface.f90
      lcomposite_output_bub_d = .FALSE.
      composite_file_pattern_bub_d(:) = ' '   ! empty pattern, leads to default file name
      eleind_for_composite_bub_glob_d = -99  ! Elevation index to be used for the composite to detect the need for bubbles
      ! Not yet in the namelist regarding the composite:
      lsmooth_composite_bub_glob = .FALSE.    ! If composite is smoothed by binomial filter
      nsmoothpoints_for_comp_bub_glob = 9     ! width of symetric 2D binomial smoother in grid points
      nfilt_for_comp_bub_glob = 1             ! number of consecutive applications of the smoother
      ! For the bubble detection and advection:
      dt_bubble_search_d     = 900.0_dp          ! Time interval from one automatic bubble search to the next [s]
      tstart_bubble_search_d = 0.0_dp
      tend_bubble_search_d   = HUGE(1.0_dp)
      t_offset_bubble_trigger_async_d = 0.0_dp  ! In case of asynchr. radar IO for runtime optimization:
                                                !   Time delay of bubble triggering after detection step.
                                                !   It adds to the bubble advection time dt_bubble_advect below.
                                                !   The worker PEs of the model may continue in parallel to the async IO PEs
                                                !   for this amount of model time, until they synchronize with the IO PEs to
                                                !   to gather information on new bubbles. An offset of 0.0 would be physically ideal but
                                                !   but destroys the runtime advantage of asynchr. IO.
      prob_bubble_d        = 1.0_dp             ! Probability of triggering a bubble when it is found [0-1]
      lbub_isolated_d      = .FALSE.            ! Check that the bubble candidate objects are isolated from other objects in observations
      maxdim_obs_d         = 75000.0_dp         ! Maximum dimension of observed cells that can trigger a bubble (very large objects are not targeted) [meters]
      threshold_obs_d(1:2) = (/ 25.0, 30.0/)    ! Thresholds for obs that define the minimum dBZ in an object- and the high intensity region [dBZ]
      threshold_mod_d(1:2) = (/ 25.0, 30.0/)    ! Thresholds for sim that define the minimum dBZ in an object- and the high intensity region [dBZ]
      areamin_mod_d(1:2)   = (/ 25e6, 9e6/)     ! Minimum area in obs of the object- and high intensity region for detecting an object  [m^2]
      areamin_obs_d(1:2)   = (/ 25e6, 9e6/)     ! Minimum area in sim of the object- and high intensity region for detecting an object  [m^2]
      mult_dist_obs_d      = 1.0_dp             ! Multiplicative axis factor for the obs object [-]. If it is high, separation- and isolation checks are more restrictive
      mult_dist_mod_d      = 1.0_dp             ! Multiplicative axis factor for the sim object [-]. If it is high, separation- and isolation checks are more restrictive
      add_dist_obs_d       = 10000.0_dp         ! Additive axis increase for the obs object [meters]. If it is high, separation- and isolation checks are more restrictive
      add_dist_mod_d       = 10000.0_dp         ! Additive axis increase for the sim object [meters]. If it is high, separation- and isolation checks are more restrictive
      dt_bubble_advect_d   = 300.0_dp           ! Time scale for downstream advection of automatic bubbles to compensate for transport effects during initial bubble formation [seconds]
      zlow_meanwind_bubble_advect_d = 3000.0_dp ! The lower bound of averaging height interval for bubble advection speed [meters AMSL]
      zup_meanwind_bubble_advect_d  = 6000.0_dp ! The upper bound of averaging height interval for bubble advection speed [meters AMSL]
      bubble_type_d(:)     = ' '; bubble_type_d = 'cos-hrd'  ! Type of perturbation 'cos-hrd', 'cos-instant'
      bubble_heatingrate_d = 3.0_dp / 200.0_dp  ! Constant heating rate for 'cos-hrd' bubbles [K/s]
      bubble_timespan_d    = 200.0_dp           ! Total duration for release of temperature disturbance for 'cos-hrd' bubble [s]
      bubble_dT_d          = 3.0_dp             ! Temperature increment for 'cos-instant' bubbles [K]
      bubble_centz_d       = 2000.0_dp          ! Center height (Z) of temperature disturbances [m]
      bubble_radx_d        = 7500.0_dp          ! Horizontal main axis (radius) in X (longitudinal) direction of temperature disturbances [m]
      bubble_rady_d        = 7500.0_dp          ! Horizontal main axis (radius) in Y (latitudinal) direction of temperature disturbances [m]
      bubble_radz_d        = 1400.0_dp          ! Vertical main axis (radius) of temperature disturbances [m]
      bubble_rotangle_d    = 0.0_dp             ! Rotation angle of main axes of temperature disturbances [degrees]
      bubble_addnoise_T_d  = .FALSE.            ! .TRUE. is not yet implemented in COSMO/ICON for automatic bubbles!
        bubble_dT_noise_d  = 0.1_dp             !           not yet implemented in COSMO/ICON for automatic bubbles!
      bubble_holdrhconst_d = .TRUE.             ! Whether or not to keep the relative humidity constant during heating

    ! For doing manual bubbles, which are determined by namelist /ARTIFCTL/ (ONLY COSMO!):
    ldo_bubbles_manual_d   = .FALSE.

    lreadmeta_from_netcdf_d     = .FALSE.
      lcheck_inputrecords_d     = .FALSE.    ! sub-switch
      lequal_azi_alldatasets_d  = .TRUE.     ! sub-switch
      labort_if_problems_obsfiles_d = .TRUE.
      lvoldata_output_d = .TRUE.
    DO i=1, noutstreams_max
      voldata_ostream_d(i)%format(:) = ' ';  voldata_ostream_d(:)%format = 'ascii'
      voldata_ostream_d(i)%grib2_packingtype(:) = ' ';  voldata_ostream_d(:)%grib2_packingtype = 'grid_ccsds'
      voldata_ostream_d(i)%output_list(:)(:) = ' '! per default, there is no active output stream!
                                                  ! An active stream CONTAINS at least one non-empty element in output_list(:):
                                                  !   - 'zrsim', 'ahsim', 'zrobs', etc.
                                                  !   - If you set 'all', then you get all available vars.
      voldata_ostream_d(i)%file_pattern(:) = ' '  ! If empty, the default name will be chosen, depending on the data format
      voldata_ostream_d(i)%output_subdir(:) = ' ' ! If empty, store the output files in ydirradarout, otherwise in the respective subdirectory therein
      voldata_ostream_d(i)%pat_scantype_volscan(:) = ' '    ! Defaults for <scantype> placeholder
      voldata_ostream_d(i)%pat_scantype_volscan = 'volscan'
      voldata_ostream_d(i)%pat_scantype_precipscan(:) = ' '
      voldata_ostream_d(i)%pat_scantype_precipscan = 'precipscan'
      voldata_ostream_d(i)%content_tref = 0.0_dp  ! Ref.time for time batches for output-files = 0.0 [s], i.e., batches synchr. to model start time
      voldata_ostream_d(i)%content_dt   = 1e9_dp  ! a very large time in seconds to indicate that an output-file should
                                                  !  contain the whole forecast time
    END DO
    
    labort_if_problems_gribout_d = .TRUE.
    lfdbk_output_d  = .FALSE.
    lwrite_ready_d  = .FALSE.
      ready_file_pattern_d(:) = ' '   ! empty pattern, leads to default file name
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
    lmodfield_output_d = .FALSE.
      loutnwp_d = .FALSE.
    loutvolaux_d = .FALSE.
  !
    lqc_flag_d            = .TRUE.
    ldealiase_vr_obs_d    = .TRUE.
    ltestpattern_hydrometeors_d = .FALSE.

    ! .. switch to choose non-blocking mpi_isend - mpi_irecv - mpi_wait communication
    !    (=more efficient, but needs more memory) for collecting radar data on radar-IO-PEs
    !    (comm. from workers -> radar-IO-PEs)
    lcomm_nonblocking_output_d = .TRUE.

    nradsta_namelist_d = 0

    ydirradarin_d(:) = ' '
    ydirlistfile_d(:) = ' '
    ydirradarout_d(:) = ' '
    ysubdirfof_d(:) = ' '
    ysubdircomp_d(:) = ' '
    ydir_mielookup_read_d(:) = ' '
    ydir_mielookup_write_d(:) = ' '
    ydir_ready_write_d(:) = ' '

    itype_mpipar_lookupgen_d = 2  ! 1 = parallelization over dbzmeta-sets and hydrometeors
                                  ! 2 = parallelization over table elements (more efficient on parallel computers)
    pe_start_lookupgen_d     = 0  ! Start PE of the gang from icomm_cart_fwo which computes the lookup tables
    pe_end_lookupgen_d       = num_compute_fwo-1 ! End PE of the gang
    llookup_interp_mode_dualpol_d = .FALSE.  ! if to use exact same lookup table interpolation method
                                             ! for all radar moments (cubic in lin-lin-space).
                                             ! The default is the (historically grown) previous method:
                                             ! - zh, zv, zvh, ah linear in log-log space
                                             ! - rrhv, irhv, kdp, adp cubic in lin-lin space
    
    ! .. Parameters in struct dbz_meta_glob for reflectivity calculations:
    !      Initialize default-values structure dbz_meta_glob_d with global defaults
    !      from dbz_namlst_d and adjust some entries later depending on
    !      the radar meta data structure
    dbz_meta_glob_d = dbz_namlst_d
    dbz_meta_glob_d%itype_refl = 3   ! dbz_namlst_d%itype_refl is 4 and needs to be 4 because of the COSMO implementation of standalone-calls to calc_dbz_vec() for DBZ model grid output

    
    ! .. dbz_namlst_d is copied to a vector with size 1, to compensate for a possible compiler bug
    !    encountered in the past with some compilers when copying dbz_namlst_d to vector elements of the derived type:
    dbz_meta_d(1) = dbz_namlst_d

    
    ngpsm_h_glob_d = -99
    ngpsm_v_glob_d = -99

    icountry_d = i_dwd   ! Country flag for the radar meta data background list:
                         !  i_dwd = Germany
                         !  i_meteoswiss = Switzerland
                         !  i_arpasim = Italy

    itype_supobing_d         = 0              ! 0: no superobing; 1: averaging;  2: median
    supob_nrb_d              = 3              ! lower threshold for number of radar bins used for superobbing
    supob_azi_maxsector_vr_d = 40.0_dp        ! maximal azimut sector (symetrical to its center) for v_r superobing [deg]
    supob_cart_resolution_d  = 20000.0_dp     ! resoltion of the cartesian grid for superobing [m]
    supob_ave_width_d        = SQRT(2.0_dp) * supob_cart_resolution_d ! width of averaging area for superobing [m]
    supob_minrange_vr_d      = 0.75_dp*supob_ave_width_d  ! min. range of superobing points for vr [m]
    supob_minrange_z_d       = 0.75_dp*supob_ave_width_d  ! min. range of superobing points for z  [m]
    supob_vrw_d              = 10.0_dp                    ! threshold for stddev of radial wind values within a superobed bin [m/s]
    supob_rfl_d              = SQRT(HUGE(1.0_dp))  ! threshold for stddev of reflectivity values within a superobed bin [mm^6/m^3]
    supob_lowthresh_z_obs_d  = miss_value    ! dBZ
    supob_lowthresh_z_sim_d  = miss_value    ! dBZ
    itype_metric_refl_fdbk_d = 1  ! namelist setting to decide about the metric for reflectivity in the feedback files:
                                  !  = 1: write dBZ to feedback files
                                  !  = 2: convert to an effective LWC (0.004*Z)^(0.55)
                                  !  = 3: convert to effective LWC as for (2) and set obs-dependent obs error factor
    minval_obserr_lwc_d = 5e-4    ! for itype_metric_refl_fdbk = 3: asymptotic obs error for LWC -> 0.0

    itype_obserr_vr_d        = 0  ! option for observation error for radial wind:
                                  !   = 0: constant observation error (=baseval_obserr_vr)
                                  !   = 1: e_o is a linear ramp as function of superobe'd reflectivity (see sketch below)
                                  !   = 2: superobe'd e_o + e_o-weighted superobe'd radial wind. For this,
                                  !        e_o is computed from reflectivity before superobing (same linear ramp as for option 1) and is used:
                                  !        a) as second weight (inverse e_o) additionally to spatial Cressman-weight in superobbing of radial wind
                                  !        b) for superobbing of e_o itself with the same weighting
    ramp_lowdbz_obserr_vr_d  = 0.0_dp             ! for options   1/2: lower ramp dBZ-threshold for increasing observation error in radial wind
    ramp_highdbz_obserr_vr_d = 10.0_dp            ! for options   1/2: upper ramp dBZ-threshold for increasing observation error in radial wind
    maxval_obserr_vr_d       = 10.0_dp            ! for options   1/2: obs error for dBZ-values <  ramp_highdbz_obserr_vr
    baseval_obserr_vr_d      = 1.0_dp             ! for options 0/1/2: obs error for dBZ-values >= ramp_highdbz_obserr_vr (options 1/2) or general value (option 0)
                                                  !   in general: 1.0 means that the obserr for vr is generally defined as a relative error
                                                  !               2.5 would be a good choice for an absolute error

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


    thin_step_azi_d   = 1   ! if no superobing: step widths of data thinning
    thin_step_range_d = 1
    thin_step_ele_d   = 1

    ind_ele_fdbk_glob_d (1:nel_max) = missval_int    ! default: neutral value
    ind_ele_voldata_glob_d (1:nel_max) = missval_int ! default: neutral value

    obs_times_fdbk_glob_d (1:nobstimes_max) = unused_value    ! default: missing value to indicate that all obs times should be written to feedback
    obs_times_voldata_glob_d (1:nobstimes_max) = unused_value ! default: missing value to indicate that all obs times should be written to feedback
    dt_obs_fdbk_glob_d = unused_value                         ! default: missing value to indicate that all obs times should be written to feedback
    dt_obs_voldata_glob_d = unused_value                      ! default: missing value to indicate that all obs times should be written to feedback

    ra_inc_coarse_glob_d = unused_value  ! default: missing value to indicate that the individual ra_inc of each station is used
    
    ! .. Initialize all character strings in the rs_meta type with blanks:
    !     (otherwise there might be problems with filenames after mpi_bcast of this type)
    CALL rsm_init_strings_blanks(rs_meta)
      
    !-------------------------------------------------------------------------------
    ! Section 2: A First pass of reading the namelist.
    !   This is only for determining the country code of the background meta data list
    !   respectively the available observational radar data, the dbz_meta_glob,
    !   the ngpsm_h and ngpsm_v
    !-------------------------------------------------------------------------------

    IF (my_radar_id == 0) THEN

      DO ista = 1, nradsta_max
        rs_meta(ista)%icountry = icountry_d
      END DO

      nradsta_namelist = nradsta_namelist_d
      icountry         = icountry_d
      dbz_meta_glob    = dbz_meta_glob_d
      ngpsm_h_glob     = ngpsm_h_glob_d
      ngpsm_v_glob     = ngpsm_v_glob_d

      ind_ele_fdbk_glob = ind_ele_fdbk_glob_d
      ind_ele_voldata_glob = ind_ele_voldata_glob_d

      labort_if_problems_gribout = labort_if_problems_gribout_d

      ldo_composite = ldo_composite_d
        comp_meta = comp_meta_d
        lcomposite_output = lcomposite_output_d
        nel_composite = nel_composite_d
        eleindlist_for_composite_glob = eleindlist_for_composite_glob_d
        levelidlist_for_composite_glob = levelidlist_for_composite_glob_d
        comp_grib2_packingtype = comp_grib2_packingtype_d

      ldo_bubbles = ldo_bubbles_d
        comp_meta_bub = comp_meta_bub_d
        lcomposite_output_bub = lcomposite_output_bub_d
        eleind_for_composite_bub_glob = eleind_for_composite_bub_glob_d
        dt_bubble_search            =   dt_bubble_search_d
        tstart_bubble_search        =   tstart_bubble_search_d
        tend_bubble_search          =   tend_bubble_search_d
        t_offset_bubble_trigger_async = t_offset_bubble_trigger_async_d
        prob_bubble                 =   prob_bubble_d
        lbub_isolated               =   lbub_isolated_d
        maxdim_obs                  =   maxdim_obs_d
        threshold_obs(1:2)          =   threshold_obs_d(1:2)
        threshold_mod(1:2)          =   threshold_mod_d(1:2)
        areamin_mod(1:2)            =   areamin_mod_d(1:2)
        areamin_obs(1:2)            =   areamin_obs_d(1:2)
        mult_dist_obs               =   mult_dist_obs_d
        mult_dist_mod               =   mult_dist_mod_d
        add_dist_obs                =   add_dist_obs_d
        add_dist_mod                =   add_dist_mod_d
        dt_bubble_advect            =   dt_bubble_advect_d
        zlow_meanwind_bubble_advect =   zlow_meanwind_bubble_advect_d
        zup_meanwind_bubble_advect  =   zup_meanwind_bubble_advect_d
        bubble_type                 =   bubble_type_d
        bubble_heatingrate          =   bubble_heatingrate_d
        bubble_timespan             =   bubble_timespan_d
        bubble_dT                   =   bubble_dT_d
        bubble_centz                =   bubble_centz_d
        bubble_radx                 =   bubble_radx_d
        bubble_rady                 =   bubble_rady_d
        bubble_radz                 =   bubble_radz_d
        bubble_rotangle             =   bubble_rotangle_d
        bubble_addnoise_T           =   bubble_addnoise_T_d
        bubble_dT_noise             =   bubble_dT_noise_d
        bubble_holdrhconst          =   bubble_holdrhconst_d
      ldo_bubbles_manual = ldo_bubbles_manual_d

      obs_times_fdbk_glob = obs_times_fdbk_glob_d
      obs_times_voldata_glob = obs_times_voldata_glob_d
      dt_obs_fdbk_glob    = dt_obs_fdbk_glob_d
      dt_obs_voldata_glob = dt_obs_voldata_glob_d

      ra_inc_coarse_glob         = ra_inc_coarse_glob_d

      !***********************************************************************************
      ! .. Open, read and close the namelist:
      !
      CALL read_radarnamelist (izdom)
      !
      !***********************************************************************************

      !-------------------------------------------------------------------------------
      ! .. Check nradsta_namelist:

      IF (nradsta_namelist > nradsta_max) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a,i0,a)') 'nradsta_namelist > nradsta_max (',nradsta_max,') !'
        CALL abort_run (my_radar_id, 10065, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      ! .. Parameters in struct rs_meta_d for default radar station parameters:

      ! .. Check icountry:
      IF (icountry < 0 .OR. icountry >= ncountry_max) THEN
        ! A valid ncountry is in the range [0,ncountry_max-1]
        errstring(:) = ' '
        WRITE (errstring,'(a,i0,a,i0,a)') 'icountry=',icountry,' not in valid range [0,',ncountry_max-1,'] !'
        CALL abort_run (my_radar_id, 10075, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')        
      END IF
      
      ! .. Take default values for the radar meta data structure from prototype for icountry:
      rs_meta_d(1) = get_meta_proto ( icountry )

      ! .. Adjust some entries of struct dbz_meta_glob:
      dbz_meta_glob%station_id = rs_meta_d(1)%station_id   ! Any spurious station_id in dbz_meta_glob% should be eliminated
      ! will later be overwritten in dbz_meta(ista) from the values in rs_meta(ista):
      dbz_meta_glob%lambda_radar = rs_meta_d(1)%lambda

      ! .. Pre-set rs_meta and dbz_meta with the default values:
      DO ista = 1, nradsta_max
        rs_meta(ista)  = rs_meta_d(1)
        dbz_meta(ista) = dbz_meta_glob
      END DO

      !***********************************************************************************
      ! .. Open, read and close the namelist a second time, this time for the presetting of
      !     parameters in struct rs_meta for actual radar station parameters:
      !
      CALL read_radarnamelist (izdom)
      !
      !***********************************************************************************

      ! .. Adjust some entries of struct dbz_meta_glob again:
      dbz_meta_glob%station_id = rs_meta_d(1)%station_id   ! Any spurious station_id in dbz_glob% should be eliminated
      ! will later be overwritten in dbz_meta(ista) from the values in rs_meta(ista):
      dbz_meta_glob%lambda_radar = rs_meta_d(1)%lambda

      ! .. Check some entries in dbz_meta_glob:
      IF (dbz_meta_glob%itype_refl < 1 .OR. dbz_meta_glob%itype_refl > 6) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value for dbz_meta_glob%itype_refl in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 and <= 6!'
        CALL abort_run (my_radar_id, 98077, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF (dbz_meta_glob%igraupel_type < 1 .OR. dbz_meta_glob%igraupel_type > 3) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value for dbz_meta_glob%igraupel_type in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 and <= 3!'
        CALL abort_run (my_radar_id, 98078, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF (dbz_meta_glob%itype_Dref_fmelt < 1 .OR. dbz_meta_glob%itype_Dref_fmelt > 2) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value for dbz_meta_glob%itype_Dref_fmelt in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 and <= 2!'
        CALL abort_run (my_radar_id, 98079, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      ! .. Adjust melting scheme params if (seemingly) contradictory:
      dbz_meta_glob%Tmax_min_i = MIN(dbz_meta_glob%Tmax_max_i,dbz_meta_glob%Tmax_min_i)
      dbz_meta_glob%Tmax_min_s = MIN(dbz_meta_glob%Tmax_max_s,dbz_meta_glob%Tmax_min_s)
      dbz_meta_glob%Tmax_min_g = MIN(dbz_meta_glob%Tmax_max_g,dbz_meta_glob%Tmax_min_g)
      dbz_meta_glob%Tmax_min_h = MIN(dbz_meta_glob%Tmax_max_h,dbz_meta_glob%Tmax_min_h)


      ! .. Pre-set rs_meta for each radar in the namelist depending on its rs_meta%icountry
      !     (either from namelist or from the background icountry):
      DO ista = 1, nradsta_namelist
        rs_meta(ista)               = get_meta_proto ( rs_meta(ista)%icountry )
        dbz_meta(ista)              = dbz_meta_glob
        dbz_meta(ista)%station_id   = rs_meta(ista)%station_id
        dbz_meta(ista)%lambda_radar = rs_meta(ista)%lambda
      END DO

      DO ista = nradsta_namelist+1, nradsta_max
        rs_meta(ista)  = rs_meta_d(1)
        dbz_meta(ista) = dbz_meta_glob
      END DO

      
      !-------------------------------------------------------------------------------
      ! .. Set rs_meta(ista)ngpsm_h/v = ngpsm_h_glob resp. ngpsm_v_glob for all stations
      !    if ngpsm_h_glob and/or ngpsm_v_glob have been specified in the namelist.
      !    This is then the background value for all radars, which can
      !    later be overwritten for the single stations via namelist.
      !    If no value or -99 for ngpsm_h_glob and/or ngpsm_v_glob is given, the default
      !    values (1 resp. 1) will win (no smoothing!!!).
      !-------------------------------------------------------------------------------

      IF (ngpsm_h_glob >= 1) THEN
        DO ista = 1, nradsta_max
          rs_meta(ista)%ngpsm_h = ngpsm_h_glob
        END DO
      ELSE IF (ngpsm_h_glob /= ngpsm_h_glob_d) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value ngpsm_h_glob in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 or -99!'
        CALL abort_run (my_radar_id, 10079, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF (ngpsm_v_glob >= 1) THEN
        DO ista = 1, nradsta_max
          rs_meta(ista)%ngpsm_v = ngpsm_v_glob
        END DO
      ELSE IF (ngpsm_v_glob /= ngpsm_v_glob_d) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value ngpsm_v_glob in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 or -99!'
        CALL abort_run (my_radar_id, 10079, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      !-------------------------------------------------------------------------------
      ! .. Similarly for rs_meta(ista)%ra_inc_coarse = ra_inc_coarse_glob for all stations
      !    if ra_inc_coarse_glob has been specified in the namelist.
      !    This is then the background value for all radars, which can
      !    later be overwritten for the single stations via namelist.
      !    If no value for ra_inc_coarse_glob is given, the default
      !    value will win (= -999.99).
      !-------------------------------------------------------------------------------

      IF (ra_inc_coarse_glob > 0.0_dp) THEN
        rs_meta(:)%ra_inc_coarse = ra_inc_coarse_glob
      END IF
      
      !-------------------------------------------------------------------------------
      ! .. Similarly for rs_meta(ista)%eleindlist_for_composite = eleindlist_for_composite_glob for all stations.
      !    This is then the background values for all radars, which can
      !    later be overwritten for the single stations via namelist.
      !-------------------------------------------------------------------------------

      IF (ldo_composite .AND. nel_composite < 1) THEN
        ldo_composite = .FALSE.
        IF (my_radar_id == 0) THEN
          WRITE (*,'(a)') 'WARNING input_radarnamelist(): ldo_composite = .TRUE. but nel_composite < 1!'// &
               ' ldo_composite is set to .FALSE.!'
        END IF
      END IF

      IF (ldo_composite) THEN
        DO i=1, nel_composite
          IF (eleindlist_for_composite_glob(i) >= 1) THEN
            DO ista = 1, nradsta_max
              rs_meta(ista)%eleindlist_for_composite(i) = eleindlist_for_composite_glob(i)
            END DO
          ELSE
            errstring(:) = ' '
            WRITE (errstring,'(a,i2.2,a,i0,a,i0,a)') 'Wrong value eleindlist_for_composite_glob(',i,')=', &
                 eleindlist_for_composite_glob(i),' in '//trim(radarnmlfile)//'.' // &
                 ' Must be >= 1 and <= nel_max = ',nel_max,'!'
            CALL abort_run (my_radar_id, 10079, &
                 'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                 'radar_namelist_read.f90, input_radarnamelist()')
          END IF
        END DO
      END IF

      !-------------------------------------------------------------------------------
      ! .. Similarly for rs_meta(ista)%eleind_for_composite_bub = eleind_for_composite_bub_glob for all stations
      !    if eleind_for_composite_bub_glob has been specified in the namelist.
      !    This is then the background value for all radars, which can
      !    later be overwritten for the single stations via namelist.
      !    If no value for eleind_for_composite_bub_glob is given, the default
      !    value will win (= -99).
      !-------------------------------------------------------------------------------

      IF (ldo_bubbles .AND. eleind_for_composite_bub_glob >= 1) THEN
        DO ista = 1, nradsta_max
          rs_meta(ista)%eleind_for_composite_bub = eleind_for_composite_bub_glob
        END DO
      ELSE IF (ldo_bubbles .AND. eleind_for_composite_bub_glob == eleind_for_composite_bub_glob_d) THEN
        WRITE (*,'(a)') 'INFO: No value for eleind_for_composite_bub_glob in '//trim(radarnmlfile)//'.' // &
             ' We use the value 1 (first elevation).'
        DO ista = 1, nradsta_max
          rs_meta(ista)%eleind_for_composite_bub = 1
        END DO
      ELSE IF (ldo_bubbles .AND. eleind_for_composite_bub_glob /= eleind_for_composite_bub_glob_d) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value eleind_for_composite_bub_glob in '//trim(radarnmlfile)//'.' // &
             ' Must be >= 1 or -99!'
        CALL abort_run (my_radar_id, 10079, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

#ifdef __COSMO__
      !-------------------------------------------------------------------------------
      ! .. For automatic bubble generator:
      !    Noise on automatic bubbles is not yet implemented in COSMO src_artifdata.f90,
      !    SR artif_heatrate_dist()

      IF (ldo_bubbles .AND. bubble_addnoise_T) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value bubble_addnoise_T in '//trim(radarnmlfile)//'.' // &
             ' Has to be .FALSE. (default) because this option is not yet implemented in COSMO for automatic bubbles!'
        CALL abort_run (my_radar_id, 10979, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
#endif

      !-------------------------------------------------------------------------------
      ! .. Override any user setting of the earth radius of the composite grid:

      comp_meta    %r_earth = r_earth_dp
      comp_meta_bub%r_earth = r_earth_dp

      !-------------------------------------------------------------------------------
      ! .. Check user choice of the grib2 packing type for composite output:

      SELECT CASE (TRIM(comp_grib2_packingtype))
      CASE ('grid_simple','grid_ccsds','png')
      CASE default
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value comp_grib2_packingtype='//TRIM(comp_grib2_packingtype)// &
             ' in '//TRIM(radarnmlfile)//'.' // &
             ' Possible values are ''grid_simple'', ''grid_ccsds'' or ''png''!'
        CALL abort_run (my_radar_id, 10980, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')

      END SELECT

    END IF    !  my_radar_id == 0

    IF (num_radar > 1) THEN
      CALL distribute_values_radar (nradsta_namelist,               1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (icountry,                       1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (ngpsm_h_glob,                   1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (ngpsm_v_glob,                   1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (nel_composite,                  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (eleind_for_composite_bub_glob,  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (eleindlist_for_composite_glob,  nel_composite_max, 0, icomm_radar, ierr)
      CALL distribute_values_radar (levelidlist_for_composite_glob, nel_composite_max, 0, icomm_radar, ierr)
      CALL distribute_values_radar (obs_times_fdbk_glob,            nobstimes_max,     0, icomm_radar, ierr)
      CALL distribute_values_radar (obs_times_voldata_glob,         nobstimes_max,     0, icomm_radar, ierr)
      CALL distribute_values_radar (dt_obs_fdbk_glob,               SIZE(dt_obs_fdbk_glob),    0, icomm_radar, ierr)
      CALL distribute_values_radar (dt_obs_voldata_glob,            SIZE(dt_obs_voldata_glob), 0, icomm_radar, ierr)
      CALL distribute_values_radar (ra_inc_coarse_glob,             1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (labort_if_problems_gribout,     1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (ldo_composite,                  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (lcomposite_output,              1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (comp_grib2_packingtype, LEN(comp_grib2_packingtype), 0, icomm_radar, ierr)
      CALL distribute_values_radar (ldo_bubbles,                    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (ldo_bubbles_manual,             1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (lcomposite_output_bub,          1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (dt_bubble_search,               1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (tstart_bubble_search,           1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (tend_bubble_search,             1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (t_offset_bubble_trigger_async,  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (prob_bubble,                    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (lbub_isolated,                  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (maxdim_obs,                     1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (threshold_obs(1:2),             2,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (threshold_mod(1:2),             2,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (areamin_mod(1:2),               2,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (areamin_obs(1:2),               2,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (mult_dist_obs,                  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (mult_dist_mod,                  1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (add_dist_obs,                   1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (add_dist_mod,                   1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (dt_bubble_advect,               1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (zlow_meanwind_bubble_advect,    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (zup_meanwind_bubble_advect,     1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_type,                    LEN(bubble_type),  0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_heatingrate,             1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_timespan,                1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_dT,                      1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_centz,                   1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_radx,                    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_rady,                    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_radz,                    1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_rotangle,                1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_addnoise_T,              1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_dT_noise,                1,                 0, icomm_radar, ierr)
      CALL distribute_values_radar (bubble_holdrhconst,             1,                 0, icomm_radar, ierr)

      CALL distribute_values_radar (ind_ele_fdbk_glob,              nel_max,           0, icomm_radar, ierr)
      CALL distribute_values_radar (ind_ele_voldata_glob,           nel_max,           0, icomm_radar, ierr)
      CALL mpi_bcast               (comp_meta,                   1, mpi_comp_meta_typ, 0, icomm_radar, ierr)
      CALL mpi_bcast               (comp_meta_bub,               1, mpi_comp_meta_typ, 0, icomm_radar, ierr)

    END IF


    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    ! Section 3: read the namelist on proc 0, third pass.
    !            This time it is for all the parameters.
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------

    IF (my_radar_id == 0) THEN

      itype_supobing   = itype_supobing_d
      itype_obserr_vr  = itype_obserr_vr_d
      itype_metric_refl_fdbk = itype_metric_refl_fdbk_d
      supob_nrb        = supob_nrb_d
      thin_step_azi    = thin_step_azi_d
      thin_step_range  = thin_step_range_d
      thin_step_ele    = thin_step_ele_d

      supob_azi_maxsector_vr = supob_azi_maxsector_vr_d
      supob_cart_resolution  = supob_cart_resolution_d
      supob_ave_width        = supob_ave_width_d
      supob_minrange_vr      = supob_minrange_vr_d
      supob_minrange_z       = supob_minrange_z_d
      supob_vrw              = supob_vrw_d
      supob_rfl              = supob_rfl_d
      supob_lowthresh_z_obs  = supob_lowthresh_z_obs_d
      supob_lowthresh_z_sim  = supob_lowthresh_z_sim_d
      ramp_highdbz_obserr_vr = ramp_highdbz_obserr_vr_d
      ramp_lowdbz_obserr_vr  = ramp_lowdbz_obserr_vr_d
      baseval_obserr_vr      = baseval_obserr_vr_d
      maxval_obserr_vr       = maxval_obserr_vr_d
      minval_obserr_lwc      = minval_obserr_lwc_d

      ldebug_radsim = ldebug_radsim_d
      lout_geom     = lout_geom_d
      loutradwind   = loutradwind_d
      loutdbz       = loutdbz_d
      lcalc_dbz_on_radarbins = lcalc_dbz_on_radarbins_d
      loutpolstd    = loutpolstd_d
      loutpolall    = loutpolall_d
      lextdbz       = lextdbz_d
      lweightdbz = lweightdbz_d
      lfall         = lfall_d
      lfill_vr_backgroundwind = lfill_vr_backgroundwind_d
      lonline       = lonline_d
      lcomm_nonblocking_online = lcomm_nonblocking_online_d
      lsode         = lsode_d
      lsmooth       = lsmooth_d
      lmds_z        = lmds_z_d
      lmds_vr       = lmds_vr_d
      lreadmeta_from_netcdf = lreadmeta_from_netcdf_d
      lcheck_inputrecords   = lcheck_inputrecords_d
      lequal_azi_alldatasets = lequal_azi_alldatasets_d
      labort_if_problems_obsfiles = labort_if_problems_obsfiles_d
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
      lmodfield_output = lmodfield_output_d
      loutnwp          = loutnwp_d
      loutvolaux       = loutvolaux_d
  !
      lvoldata_output = lvoldata_output_d
      lfdbk_output  = lfdbk_output_d
      lwrite_ready  = lwrite_ready_d
      lqc_flag      = lqc_flag_d
      ldealiase_vr_obs = ldealiase_vr_obs_d
      ltestpattern_hydrometeors_d = ltestpattern_hydrometeors
      lcomm_nonblocking_output = lcomm_nonblocking_output_d
      llookup_interp_mode_dualpol = llookup_interp_mode_dualpol_d

      ydirradarin   = ydirradarin_d
      ydirlistfile  = ydirlistfile_d
      ydirradarout  = ydirradarout_d
      ysubdirfof    = ysubdirfof_d
      ysubdircomp   = ysubdircomp_d
      ydir_mielookup_read = ydir_mielookup_read_d
      ydir_mielookup_write= ydir_mielookup_write_d
      ydir_ready_write=ydir_ready_write_d
      composite_file_pattern = composite_file_pattern_d
      composite_file_pattern_bub = composite_file_pattern_bub_d
      ready_file_pattern = ready_file_pattern_d

      itype_mpipar_lookupgen = itype_mpipar_lookupgen_d
      pe_start_lookupgen     = pe_start_lookupgen_d
      pe_end_lookupgen       = pe_end_lookupgen_d

      voldata_ostream        = voldata_ostream_d
      
      !***********************************************************************************
      ! .. Open, read and close the namelist:
      !
      CALL read_radarnamelist (izdom)
      !
      !***********************************************************************************

      !-------------------------------------------------------------------------------
      !- Section 4a: Check values for errors and consistency
      !-------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      ! .. Check loutpolstd and loutpolall. Should only be .TRUE. if
      !    itype_refl = 1 .OR. itype_refl > 4

      IF (loutpolstd .AND. dbz_meta_glob%itype_refl > 1 .AND. dbz_meta_glob%itype_refl <= 4) THEN
        WRITE (errstring,'(a,i0,a)') 'loutpolstd=.TRUE. but dbz_meta_glob%itype_refl=', &
             dbz_meta_glob%itype_refl,' in '//TRIM(radarnmlfile)//'.' // &
             ' Can only be .TRUE. if dbz_meta_glob%itype_refl=1,5,6!'
        CALL abort_run (my_radar_id, 10979, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF (loutpolall .AND. dbz_meta_glob%itype_refl > 1 .AND. dbz_meta_glob%itype_refl <= 4) THEN
        WRITE (errstring,'(a,i0,a)') 'loutpolall=.TRUE. but dbz_meta_glob%itype_refl=', &
             dbz_meta_glob%itype_refl,' in '//TRIM(radarnmlfile)//'.' // &
             ' Can only be .TRUE. if dbz_meta_glob%itype_refl=1,5,6!'
        CALL abort_run (my_radar_id, 10979, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      
      !-------------------------------------------------------------------------------
      ! .. Check that every radar in the namelist has a valid station_id and
      !     and set the scanname according to the given elevations:

      DO ista = 1, nradsta_namelist
        IF (rs_meta(ista)%station_id == rs_meta_d(1)%station_id .OR. &
             rs_meta(ista)%station_id <= 0) THEN
          WRITE (*,'(a,i0,a)') 'WARNING input_radarnamelist(): No station_id given for station ', &
               ista,'. Might become critical later if operator not in obs data mode!'
        END IF
        ! .. set the scanname according to the given elevations if not in obs data mode.
        !     Recognizes also DWD precip scans:
        IF (.NOT.lreadmeta_from_netcdf) THEN
          CALL set_scanname ( rs_meta(ista) )
        END IF
        
        ! .. Adjust melting scheme params of indiv station if (seemingly) contradictory:
        ! ?? Do we actually need this? Done already on dbz_meta_glob.
        ! ?? But indiv. stations might have different settings.
        ! ?? Hence, rather do here and NOT on dbz_meta_glob?
        dbz_meta(ista)%Tmax_min_i = MIN(dbz_meta(ista)%Tmax_max_i,dbz_meta(ista)%Tmax_min_i)
        dbz_meta(ista)%Tmax_min_s = MIN(dbz_meta(ista)%Tmax_max_s,dbz_meta(ista)%Tmax_min_s)
        dbz_meta(ista)%Tmax_min_g = MIN(dbz_meta(ista)%Tmax_max_g,dbz_meta(ista)%Tmax_min_g)
        dbz_meta(ista)%Tmax_min_h = MIN(dbz_meta(ista)%Tmax_max_h,dbz_meta(ista)%Tmax_min_h)
      END DO

      !-------------------------------------------------------------------------------
      ! .. For automatic bubble generator:
      !    Time delay for bubble triggering after detection for asynchr. IO must
      !    be smaller or equal than the time from one detection to the next:

      IF (ldo_bubbles .AND. lreadmeta_from_netcdf .AND. loutdbz) THEN
        IF (tstart_bubble_search < 0.0_dp) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a)') 'Wrong value tstart_bubble_search in '//TRIM(radarnmlfile)//'.' // &
               ' Must be >= 0.0 s!'
          CALL abort_run (my_radar_id, 10979, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
        IF (dt_bubble_search <= 0.0_dp) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a)') 'Wrong value dt_bubble_search in '//TRIM(radarnmlfile)//'.' // &
               ' Must be > 0.0 s!'
          CALL abort_run (my_radar_id, 10980, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
        IF (tend_bubble_search <= 0.0_dp) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a)') 'Wrong value tend_bubble_search in '//TRIM(radarnmlfile)//'.' // &
               ' Must be > 0.0 s!'
          CALL abort_run (my_radar_id, 10981, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
        IF (num_radario > 0) THEN
          IF (t_offset_bubble_trigger_async < 0.0_dp) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a)') 'Wrong value t_offset_bubble_trigger_async in '//TRIM(radarnmlfile)//'.' // &
                 ' Must be >= 0.0 s!'
            CALL abort_run (my_radar_id, 10982, &
                 'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                 'radar_namelist_read.f90, input_radarnamelist()')
          END IF
          IF (t_offset_bubble_trigger_async > dt_bubble_search) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a,f0.1,a)') 'Wrong value t_offset_bubble_trigger_async in '//TRIM(radarnmlfile)//'.' // &
                 ' Must be <= dt_bubble_search (',dt_bubble_search ,' s)!'
            CALL abort_run (my_radar_id, 10983, &
                 'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                 'radar_namelist_read.f90, input_radarnamelist()')
          END IF
        END IF
      END IF

      !-------------------------------------------------------------------------------
      ! .. Check voldata_output_list and enforce full list of all possible fields, if none
      !     is given in the namelist:

      DO i=1, noutstreams_max
        IF ( ANY( INDEX( voldata_ostream(i)%output_list(:), 'all') > 0 ) ) THEN

          ! reset output list to empty:
          voldata_ostream(i)%output_list(:)(:) = ' '

          ! and fill it with the list of all possible fields:
          voldata_ostream(i)%output_list(1)  = 'zrsim'
          voldata_ostream(i)%output_list(2)  = 'zdrsim'
          voldata_ostream(i)%output_list(3)  = 'rhvsim'
          voldata_ostream(i)%output_list(4)  = 'kdpsim'
          voldata_ostream(i)%output_list(5)  = 'phidpsim'
          voldata_ostream(i)%output_list(6)  = 'ldrsim'
          voldata_ostream(i)%output_list(7)  = 'vrsim'

          voldata_ostream(i)%output_list(8)  = 'zrobs'
          voldata_ostream(i)%output_list(9)  = 'zdrobs'
          voldata_ostream(i)%output_list(10) = 'rhvobs'
          voldata_ostream(i)%output_list(11) = 'kdpobs'
          voldata_ostream(i)%output_list(12) = 'phidpobs'
          voldata_ostream(i)%output_list(13) = 'ldrobs'
          voldata_ostream(i)%output_list(14) = 'vrobs'

          voldata_ostream(i)%output_list(15) = 'qzobs'
          voldata_ostream(i)%output_list(16) = 'qzdrobs'
          voldata_ostream(i)%output_list(17) = 'qrhvobs'
          voldata_ostream(i)%output_list(18) = 'qkdpobs'
          voldata_ostream(i)%output_list(19) = 'qphidpobs'
          voldata_ostream(i)%output_list(20) = 'qldrobs'
          voldata_ostream(i)%output_list(21) = 'qvobs'

          voldata_ostream(i)%output_list(22) = 'losim'
          voldata_ostream(i)%output_list(23) = 'lasim'
          voldata_ostream(i)%output_list(24) = 'hrsim'
          voldata_ostream(i)%output_list(25) = 'ersim'
          voldata_ostream(i)%output_list(26) = 'adsim'

          voldata_ostream(i)%output_list(27) = 'ahpisim'
          voldata_ostream(i)%output_list(28) = 'ahsim'
          voldata_ostream(i)%output_list(29) = 'adpsim'
          voldata_ostream(i)%output_list(30) = 'adppisim'

          voldata_ostream(i)%output_list(31) = 'vrsupsim'
          voldata_ostream(i)%output_list(32) = 'vrsupobs'
          voldata_ostream(i)%output_list(33) = 'zrsupsim'
          voldata_ostream(i)%output_list(34) = 'zrsupobs'
          voldata_ostream(i)%output_list(35) = 'zdrsupsim'
          voldata_ostream(i)%output_list(36) = 'zdrsupobs'
          voldata_ostream(i)%output_list(37) = 'rhvsupsim'
          voldata_ostream(i)%output_list(38) = 'rhvsupobs'
          voldata_ostream(i)%output_list(39) = 'kdpsupsim'
          voldata_ostream(i)%output_list(40) = 'kdpsupobs'
          !voldata_ostream(i)%output_list(  ) = 'phidpsupsim'
          !voldata_ostream(i)%output_list(  ) = 'phidpsupobs'
          !voldata_ostream(i)%output_list(  ) = 'ldrsupsim'
          !voldata_ostream(i)%output_list(  ) = 'ldrsupobs'

          voldata_ostream(i)%output_list(41) = 'losupsim'
          voldata_ostream(i)%output_list(42) = 'lasupsim'
          voldata_ostream(i)%output_list(43) = 'vasim'

          voldata_ostream(i)%output_list(44) = 'vrobserr'
          voldata_ostream(i)%output_list(45) = 'vrsupobserr'

          ! Note: currently noutput_fields_max = 50 in radar_data.f90!
          !       Increase there if needed.
        END IF
      END DO

      !-------------------------------------------------------------------------------
      ! .. Check namelist parameters for parallel lookup table generation:

      IF ( itype_mpipar_lookupgen < 1 .OR. itype_mpipar_lookupgen > 2 ) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a,i0,a)') 'Wrong value itype_mpipar_lookupgen in '//TRIM(radarnmlfile)//'.' // &
             ' Must be 1 or 2, but is ', itype_mpipar_lookupgen, ' !'
        CALL abort_run (my_radar_id, 10984, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF ( pe_start_lookupgen < 0 .OR. pe_start_lookupgen > num_compute_fwo-1 ) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a,2(i0,a))') 'Wrong value pe_start_lookupgen in '//TRIM(radarnmlfile)//'.' // &
             ' Must be >= 0 and < ', num_compute_fwo,' (num_compute_fwo), but is ', pe_start_lookupgen,'!'
        CALL abort_run (my_radar_id, 10985, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      IF ( pe_end_lookupgen < pe_start_lookupgen .OR. pe_end_lookupgen > num_compute_fwo-1 ) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a,3(i0,a))') 'Wrong value pe_end_lookupgen in '//TRIM(radarnmlfile)//'.' // &
             ' Must be >= ', pe_start_lookupgen,' and < ', num_compute_fwo,' (num_compute_fwo), but is ', pe_end_lookupgen,'!'
        CALL abort_run (my_radar_id, 10986, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      !-------------------------------------------------------------------------------
      ! .. Check lcalc_dbz_on_radarbins (option to compute radar moments on radarbins
      !     after having interpolated the model state to radar bins)

      IF (lcalc_dbz_on_radarbins .AND. lonline) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Wrong value lcalc_dbz_on_radarbins in '//TRIM(radarnmlfile)//': ' // &
             'both lonline=.TRUE. and lcalc_dbz_on_radarbins=.TRUE., ' // &
             'but the latter is not yet implemented for online beam propagation!'
        CALL abort_run (my_radar_id, 10987, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF
      
      !-------------------------------------------------------------------------------
      ! .. Check namelist parameter combinations, which are not yet implemented:

      ! NONE at the moment :-)

    END IF  ! if (my_radar_id == 0)

    !-------------------------------------------------------------------------------
    !- Section 5: Communication Part 1:
    !             Distribute global namelist parameters to all other nodes
    !-------------------------------------------------------------------------------

    IF (num_radar > 1) THEN

      ! distribute global namelist parameters
      IF (my_radar_id == 0) THEN

        nibuf = 0
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = itype_supobing
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = itype_obserr_vr
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = itype_metric_refl_fdbk
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = supob_nrb
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = thin_step_azi
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = thin_step_range
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = thin_step_ele
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = itype_mpipar_lookupgen
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = pe_start_lookupgen
        nibuf = nibuf + 1 ;   intbuf(nibuf)  = pe_end_lookupgen

        nrbuf = 0
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_azi_maxsector_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_cart_resolution
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_ave_width
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_minrange_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_minrange_z
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_vrw
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_rfl
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_lowthresh_z_obs
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = supob_lowthresh_z_sim
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = ramp_highdbz_obserr_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = ramp_lowdbz_obserr_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = baseval_obserr_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = maxval_obserr_vr
        nrbuf = nrbuf + 1 ;   realbuf(nrbuf)  = minval_obserr_lwc

        nlbuf = 0
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = ldebug_radsim
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lout_geom
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutradwind
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutdbz
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lcalc_dbz_on_radarbins
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutpolstd
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutpolall
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lextdbz
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lweightdbz
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lfall
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lfill_vr_backgroundwind
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lonline
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lcomm_nonblocking_online
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lsode
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lsmooth
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lmds_z
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lmds_vr
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lreadmeta_from_netcdf
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lcheck_inputrecords
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lequal_azi_alldatasets
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = labort_if_problems_obsfiles
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lvoldata_output
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lfdbk_output
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lwrite_ready
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = llookup_interp_mode_dualpol
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lmodfield_output
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutnwp
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = loutvolaux
  !
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lqc_flag
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = ldealiase_vr_obs
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = ltestpattern_hydrometeors
        nlbuf = nlbuf + 1 ;   logbuf(nlbuf)  = lcomm_nonblocking_output
      END IF

      ! First, distribute the buffer sizes:
      CALL distribute_values_radar (nibuf, 1, 0,       icomm_radar, ierr)
      CALL distribute_values_radar (nrbuf, 1, 0,       icomm_radar, ierr)
      CALL distribute_values_radar (nlbuf, 1, 0,       icomm_radar, ierr)

      ! Then distribute the actual buffers:
      CALL distribute_values_radar (intbuf,  nibuf, 0, icomm_radar, ierr)
      CALL distribute_values_radar (realbuf, nrbuf, 0, icomm_radar, ierr)
      CALL distribute_values_radar (logbuf,  nlbuf, 0, icomm_radar, ierr)

      IF (my_radar_id /= 0) THEN

        nibuf = 0
        nibuf = nibuf + 1 ;   itype_supobing         = intbuf(nibuf)
        nibuf = nibuf + 1 ;   itype_obserr_vr        = intbuf(nibuf)
        nibuf = nibuf + 1 ;   itype_metric_refl_fdbk = intbuf(nibuf)
        nibuf = nibuf + 1 ;   supob_nrb              = intbuf(nibuf)
        nibuf = nibuf + 1 ;   thin_step_azi          = intbuf(nibuf)
        nibuf = nibuf + 1 ;   thin_step_range        = intbuf(nibuf)
        nibuf = nibuf + 1 ;   thin_step_ele          = intbuf(nibuf)
        nibuf = nibuf + 1 ;   itype_mpipar_lookupgen = intbuf(nibuf)
        nibuf = nibuf + 1 ;   pe_start_lookupgen     = intbuf(nibuf)
        nibuf = nibuf + 1 ;   pe_end_lookupgen       = intbuf(nibuf)

        nrbuf = 0
        nrbuf = nrbuf + 1 ;   supob_azi_maxsector_vr = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_cart_resolution  = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_ave_width        = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_minrange_vr      = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_minrange_z       = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_vrw              = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_rfl              = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_lowthresh_z_obs  = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   supob_lowthresh_z_sim  = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   ramp_highdbz_obserr_vr = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   ramp_lowdbz_obserr_vr  = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   baseval_obserr_vr      = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   maxval_obserr_vr       = realbuf(nrbuf)
        nrbuf = nrbuf + 1 ;   minval_obserr_lwc      = realbuf(nrbuf)

        nlbuf = 0
        nlbuf = nlbuf + 1 ;   ldebug_radsim    = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lout_geom        = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutradwind      = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutdbz          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lcalc_dbz_on_radarbins = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutpolstd       = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutpolall       = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lextdbz          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lweightdbz       = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lfall            = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lfill_vr_backgroundwind = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lonline          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lcomm_nonblocking_online = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lsode            = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lsmooth          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lmds_z           = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lmds_vr          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lreadmeta_from_netcdf = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lcheck_inputrecords   = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lequal_azi_alldatasets=logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   labort_if_problems_obsfiles=logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lvoldata_output  = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lfdbk_output     = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lwrite_ready     = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   llookup_interp_mode_dualpol = logbuf(nlbuf)
  !.. Define in public, but only apply behind AUXOUT_OFFLINE (and __COSMO__) compiler directives
  !.. This avoids trouble in case there are accidential leftover switch settings in runscripts
        nlbuf = nlbuf + 1 ;   lmodfield_output = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutnwp          = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   loutvolaux       = logbuf(nlbuf)
  !
        nlbuf = nlbuf + 1 ;   lqc_flag         = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   ldealiase_vr_obs = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   ltestpattern_hydrometeors = logbuf(nlbuf)
        nlbuf = nlbuf + 1 ;   lcomm_nonblocking_output = logbuf(nlbuf)
      END IF

      ! distribute the path for the input radar data:
      CALL distribute_path_radar ( ydirradarin, icomm_radar )

      ! distribute the file name for directory listings:
      CALL distribute_path_radar ( ydirlistfile, icomm_radar )

      ! distribute the paths for the output radar data:
      CALL distribute_path_radar ( ydirradarout, icomm_radar )
      CALL distribute_path_radar ( ysubdirfof  , icomm_radar )
      CALL distribute_path_radar ( ysubdircomp , icomm_radar )
      CALL distribute_path_radar ( composite_file_pattern , icomm_radar )
      CALL distribute_path_radar ( composite_file_pattern_bub , icomm_radar )
      CALL distribute_path_radar ( ready_file_pattern , icomm_radar )

      ! distribute the path for the Mie lookup tables:
      CALL distribute_path_radar ( ydir_mielookup_read, icomm_radar )

      ! distribute the path for the Mie lookup tables:
      CALL distribute_path_radar ( ydir_mielookup_write, icomm_radar )

      ! distribute the path for the ready-files:
      CALL distribute_path_radar ( ydir_ready_write, icomm_radar )
      
      DO i=1, noutstreams_max
        CALL mpi_bcast (voldata_ostream(i), 1, mpi_voldata_ostream_typ, 0, icomm_radar, ierr)
      END DO
      
    END IF

    !----------------------------------------------------------------------------------------
    ! save voldata_ostream, because it has been altered after namelist reading and it will be
    !  read again below from the namelist, overwriting again the above changes:
    voldata_ostream_save = voldata_ostream

    
    IF (lreadmeta_from_netcdf) THEN

      !-------------------------------------------------------------------------------
      ! .. Check ydirradarin: Necessary here and after last call to read_radarnamelist()

      IF (LEN_TRIM(ADJUSTL(ydirradarin)) == 0) THEN
        ydirradarin(:) = ' '
        ydirradarin    = get_model_inputdir ()
      END IF
      ! .. Add slash to the end of the string if not empty:
      IF (LEN_TRIM(ADJUSTL(ydirradarin)) > 0) THEN
        ydirradarin    = TRIM(ADJUSTL(ydirradarin)) // REPEAT(' ',LEN(ydirradarin))
        IF (ydirradarin(LEN_TRIM(ydirradarin):LEN_TRIM(ydirradarin)) /= '/') THEN
          ydirradarin  = TRIM(ydirradarin) // '/'
        END IF
      END IF

      !---------------------------------------------------------------------
      !---------------------------------------------------------------------
      !
      ! .. Read meta data from NetCDF files and collect them
      !    on PE 0.
      !
      ! Goal: redefine nradsta, rs_meta(1:nradsta) and dbz_meta(1:nradsta)
      !       so that they are consistent with the number of stations
      !       and the order of the stations in the (numbered) radar obs files
      !
      !---------------------------------------------------------------------
      !---------------------------------------------------------------------

      !--------------------------------------------------------------
      ! .. Make full meta data lists depending on the country:
      !      rs_meta_dwd(1:nradsta),    dbz_meta_for_country(1:nradsta)
      !    In this process, nradsta will be re-defined on all PEs!!!
      !--------------------------------------------------------------


!!$ Data files: may contain several fields and times, but all scans have to be the same scanstrategy.
!!$  If a station has more than one, each has to be in separate files. Internally, these are
!!$  treated as different radar stations!
!!$
!!$ icountry: find this out for each radar station before reading the data!
!!$
!!$ 1) look in ydirradarin, which files are there and which radars from which country.
!!$    Set up icountry, station_id, scanstr
!!$ 2) background list should contain all known radars, get_meta_network_all()
!!$ 3) for each radar, read the metadata with a subroutine for each country/format from
!!$    the new wrapper routine read_meta_info()
!!$
!!$ Radar station = station_id + scan_id + mean(ele(k))
!!$ Composites: - are voids possible? NO, not at the moment! should enable voids!
!!$             - When are they resetted? Each time that any of the radars has something to do. The
!!$               composite is then generated by all radars that have something to do at the specific time.
!!$               The domain-decomposed output field keeps the composite-field afterwards until the
!!$               next time that any of the radars has something to do.
!!$             - This would enable to choose different radars/scan strategies for different composites.


#ifndef NETCDF
      !--------------------------------------------------------------
      ! .. abort if #ifndef NETCDF and reading of DWD netcdf data = .true.

      errstring(:) = ' '
      WRITE (errstring,'(a)') ' Reading of NetCDF radar data and its meta-data not possible ' &
           //'because preprocessor switch -DNETCDF not set!'
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
           'radar_namelist_read.f90, input_radarnamelist()')
#endif

      SELECT CASE (icountry)

      CASE ( 1 : ncountry_max )

#ifdef NETCDF
        ! .. Get rs_meta_for_country and dbz_meta_for_country, only needed at PE0:
        
        IF (my_radar_id == 0) THEN
          ALLOCATE(rs_meta_for_country(nradsta_all))
          CALL rsm_init_strings_blanks (rs_meta_for_country)
          CALL get_meta_network_all ( dbz_meta_glob, rs_meta_for_country, dbz_meta_for_country(1:nradsta_all) )
        END IF

        ! .. Get rs_meta_ncdf and collect them on my_radar_id == 0:
        !    - rs_meta_ncdf, err_meta , err_miss are only present on PE0
        !    - nradsta is distributed to all PEs (but at this point it represents only the number
        !       (of stations in obs, which might be increased/decreased by NAMELIST below)
        IF (my_radar_id == 0) THEN
          ALLOCATE(rs_meta_ncdf(nradsta_max))
        ELSE
          ALLOCATE(rs_meta_ncdf(1))
        END IF
        CALL rsm_init_strings_blanks (rs_meta_ncdf)
        CALL read_meta_info_all   ( rs_meta_ncdf, nradsta_ncdf, err_meta , err_miss )
#endif

        IF (my_radar_id == 0) THEN

          !-------------------------------------------------------------------------------
          ! .. If no radar files at all are found in the input directory,
          !    abort model run. This abort might be removed in operational mode by
          !    setting switch labort_if_problems_obsfiles=.FALSE. in the namelist.

          IF (labort_if_problems_obsfiles .AND. err_miss(1) > 0) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a)') 'No valid radar data files found!'
            CALL abort_run (my_radar_id, 10089, &
                 'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                 'radar_namelist_read.f90, input_radarnamelist()')
          END IF

          !-------------------------------------------------------------------------------
          ! .. If some other errors occur for one or all stations, either
          !    abort model run for error numbers 1, 2 and 3 and continue otherwise (error 4).
          !    Corresponding WARNINGs and ERRORs have already been issued from within read_meta_info_all().

          DO ista = 1, nradsta_ncdf
            IF (err_meta(ista) /= 0) THEN
              errstring(:) = ' '
              IF (err_meta(ista) == 1) THEN
                WRITE (errstring,'(a)') 'No consistent sets (vr,qv,z,qz) of required radar station files ' &
                     //'could be found! Check previous WARNING message(s)!'
              ELSE IF (err_meta(ista) == 2) THEN
                WRITE (errstring,'(a)') 'Too many radar stations! Increase nradsta_max in radar_data.f90!'
              ELSE IF (err_meta(ista) == 3) THEN
                WRITE (errstring,'(a)') 'Too many obs_times for a radar station! Increase nobstimes_max in radar_data.f90!'
              ELSE IF (err_meta(ista) == 4) THEN
                WRITE (errstring,'(a)') 'Discarded one or more obs_time for a radar station. Check previous WARNING message(s)!'
              ELSE IF (err_meta(ista) == 5) THEN
                WRITE (errstring,'(a)') 'Discarded one or more moments at one or more obs_time for a radar station. ' // &
                     'Check previous WARNING message(s)!'
              END IF
              IF (err_meta(ista) == 2 .OR. err_meta(ista) == 3) THEN
                CALL abort_run (my_radar_id, 10090, &
                     'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                     'radar_namelist_read.f90, '//TRIM(yzroutine))
              END IF
              IF (labort_if_problems_obsfiles) THEN
                IF (err_meta(ista) == 1 .OR. err_meta(ista) == 4 .OR. err_meta(ista) == 5) THEN
                  CALL abort_run (my_radar_id, 10091, &
                       'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                       'radar_namelist_read.f90, '//TRIM(yzroutine))
                END IF
              ELSE
                IF (err_meta(ista) == 1) THEN
                  WRITE (*,'(a,/,a)') 'WARNING: '//TRIM(ADJUSTL(errstring)), &
                       'There will be no radars simulated in this model run!'
                END IF
                IF (err_meta(ista) == 4) THEN
                  WRITE (*,'(a,/,a,i6.6,a,i0)') 'WARNING: '//TRIM(ADJUSTL(errstring)), &
                       'Number of remaining obs_times for ', rs_meta_ncdf(ista)%station_id, ': ',&
                       rs_meta_ncdf(ista)%nobs_times
                END IF
              END IF
            END IF
          END DO

        END IF  ! my_radar_id == 0

      CASE ( i_fakeobs )

        ! TEST SUITE MODE FOR TESTING FEEDBACK FILE OUTPUT WITHOUT HAVING ACTUAL RADAR OBS AT HAND

        IF (my_radar_id == 0) THEN
          ! Do not actually read any meta data from files, but use the rs_meta
          !  from the namelist. To do this, the following variables have to
          !  be specified:
          ALLOCATE(rs_meta_ncdf(nradsta_all))
          ALLOCATE(rs_meta_for_country(nradsta_all))
          nradsta                             = nradsta_namelist
          nradsta_ncdf                        = nradsta_namelist
          rs_meta_ncdf(:)                     = rs_meta(1:nradsta_all)
          rs_meta_for_country(:)              = rs_meta(1:nradsta_all)
          dbz_meta_for_country(1:nradsta_all) = dbz_meta(1:nradsta_all)
        END IF


      CASE default

        errstring(:) = ' '
        WRITE (errstring,'(a)')                                                                                      &
          'Wrong icountry for real data mode! Currently implemented is icountry == 1 (DWD), 2 (MeteoSwiss) and 3 (Italy)!'
        CALL abort_run (my_radar_id, 12077, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')

      END SELECT


      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      ! .. Some more checks and calculations:
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------

      IF (my_radar_id == 0) THEN

        !-----------------------------------------------------------------------------
        ! .. Merge meta information for grid point dbz calculation to the correct
        !    station_id:
        !-----------------------------------------------------------------------------

        DO ista = 1, nradsta_ncdf
          found = .FALSE.
          DO i = 1, nradsta_all
            IF (dbz_meta_for_country(i)%station_id == rs_meta_ncdf(ista)%station_id) THEN
              dbz_meta_ncdf(ista) = dbz_meta_for_country(i)
              found = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT.found) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a,i6.6,a)') 'DBZ metadata not in dbz_meta_for_country() for ' &
                 //'radar station ID ',rs_meta_ncdf(ista)%station_id,' ! Will discard station from the operator!'
            rs_meta_ncdf(ista)%nobs_times     = missval_int
            rs_meta_ncdf(ista)%obs_times      = unused_value
            rs_meta_ncdf(ista)%nobs_times_obs = missval_int
            rs_meta_ncdf(ista)%obs_times_obs  = unused_value
            IF (labort_if_problems_obsfiles) THEN
              CALL abort_run (my_radar_id, 10092, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, '//TRIM(yzroutine))
            ELSE
              WRITE(*,'(a)') 'ERROR '//TRIM(yzroutine)//': '//TRIM(errstring)
            END IF
          END IF
        END DO

        !-----------------------------------------------------------------------------
        ! .. Overwrite radar metadata from obs files with the values from
        !    the namelist:
        !
        !    1) Match the station IDs from NetCDF files with these from the namelist
        !       and sort corresponding entries into the namelist variables:
        !-----------------------------------------------------------------------------

        rs_meta_from_nmlst = rs_meta
        ii = 0
        DO ista = 1, nradsta_ncdf
          DO i = 1, nradsta_namelist
            IF ( rs_meta_ncdf(ista)%station_id      == rs_meta_from_nmlst(i)%station_id .AND. &
                 TRIM(rs_meta_ncdf(ista)%scanname)  == TRIM(rs_meta_from_nmlst(i)%scanname) ) THEN
              ii = ii + 1
              ista_namelist(ii) = i
              ista_ncdf(ii)     = ista
              rs_meta(i)        = rs_meta_ncdf(ista)
              dbz_meta(i)       = dbz_meta_ncdf(ista)
              ! Preparation for re-calc of obs_times after obs file metadata reading and overwriting by namelist:
              IF (rs_meta_from_nmlst(i)%lobstimes_ovwrt_recalc) THEN
                rs_meta(i)%nobs_times = missval_int
              END IF
            END IF
          END DO
        END DO
        nradsta_changed = ii

        IF (nradsta_changed > nradsta_namelist) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a)') 'Strange namelist error nradsta_changed > nradsta_namelist! ' // &
               'Doubled station_id with same scanname in rs_meta_ncdf()?'
          CALL abort_run (my_radar_id, 10381, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF

        !-----------------------------------------------------------------------------
        !    2) For stations which are contained in the namelist but not in the obs files
        !       fill the metadata with their defaults from the background list
        !-----------------------------------------------------------------------------

        DO i = 1, nradsta_namelist
          found = .FALSE.
          DO ista = 1, nradsta_ncdf
            IF ( rs_meta_ncdf(ista)%station_id      == rs_meta(i)%station_id .AND. &
                 TRIM(rs_meta_ncdf(ista)%scanname)  == TRIM(rs_meta(i)%scanname) ) THEN
              found = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. found) THEN
            found = .FALSE.
            DO ista = 1, nradsta_all
              IF ( rs_meta_for_country(ista)%station_id      == rs_meta(i)%station_id .AND. &
                   TRIM(rs_meta_for_country(ista)%scanname)  == TRIM(rs_meta(i)%scanname) ) THEN
                found = .TRUE.
                EXIT
              END IF
            END DO
            IF (found) THEN
              rs_meta(i)  = rs_meta_for_country(ista)
              dbz_meta(i) = dbz_meta_for_country(ista)
            ELSE
              errstring(:) = ' '
              WRITE (errstring,'(a,i3.3,a,i6.6,a)') 'Namelist error: station rs_meta(', i, &
                   ') with station_id ',rs_meta(i)%station_id, ' and scanname '//TRIM(rs_meta(i)%scanname)// &
                   ' is not registered in the internal background list!'
              CALL abort_run (my_radar_id, 10388, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, input_radarnamelist()')
            END IF
          END IF
        END DO

        !-------------------------------------------------------------------------------
        !   3) Set global defaults for ngpsm_h_glob and ngpsm_v_glob in rs_meta(ista)
        !      before reading the namelist again:
        !-------------------------------------------------------------------------------

        IF (ngpsm_h_glob >= 1) THEN
          DO ista = 1, nradsta_max
            rs_meta(ista)%ngpsm_h = ngpsm_h_glob
          END DO
          DO ista = 1, nradsta_ncdf
            rs_meta_ncdf(ista)%ngpsm_h = ngpsm_h_glob
          END DO
        END IF
        IF (ngpsm_v_glob >= 1) THEN
          DO ista = 1, nradsta_max
            rs_meta(ista)%ngpsm_v = ngpsm_v_glob
          END DO
          DO ista = 1, nradsta_ncdf
            rs_meta_ncdf(ista)%ngpsm_v = ngpsm_v_glob
          END DO
        END IF

        !-------------------------------------------------------------------------------
        !    Similarly for eleind_for_comosites_bub_glob and eleindlist_for_composites_glob:
        !-------------------------------------------------------------------------------

        IF (ldo_composite) THEN
          DO i=1, nel_composite
            IF (eleindlist_for_composite_glob(i) >= 1) THEN
              DO ista = 1, nradsta_max
                rs_meta(ista)%eleindlist_for_composite(i) = eleindlist_for_composite_glob(i)
              END DO
              DO ista = 1, nradsta_ncdf
                rs_meta_ncdf(ista)%eleindlist_for_composite(i) = eleindlist_for_composite_glob(i)
              END DO
            END IF
          END DO
        END IF

        IF (ldo_bubbles .AND. eleind_for_composite_bub_glob >= 1) THEN
          DO ista = 1, nradsta_max
            rs_meta(ista)%eleind_for_composite_bub = eleind_for_composite_bub_glob
          END DO
          DO ista = 1, nradsta_ncdf
            rs_meta_ncdf(ista)%eleind_for_composite_bub = eleind_for_composite_bub_glob
          END DO
        END IF

        !-------------------------------------------------------------------------------
        !    Similarly for preparing range coarsening:
        !
        !    Change the range increment rs_meta%ra_inc for input coarsening
        !    to the value of the namelist parameter ra_inc_coarse_glob, if ra_inc_coarse_glob
        !    is larger than 0.0. Later, if rs_meta%ra_inc is larger than
        !    rs_meta%ra_inc_obs, the range coarsening of input data is triggered and
        !    rs_meta%ra_inc and rs_meta%nra will characterize the coarsened range bins:
        !-------------------------------------------------------------------------------

        IF ( ra_inc_coarse_glob > 0.0_dp ) THEN
          DO ista = 1,nradsta_max
            rs_meta(ista)%ra_inc_coarse = ra_inc_coarse_glob
          END DO
          DO ista = 1,nradsta_ncdf
            rs_meta_ncdf(ista)%ra_inc_coarse = ra_inc_coarse_glob
          END DO
        ELSE
          DO ista = 1,nradsta_max
            rs_meta(ista)%ra_inc_coarse = rs_meta(ista)%ra_inc_obs
          END DO
          DO ista = 1,nradsta_ncdf
            rs_meta_ncdf(ista)%ra_inc_coarse = rs_meta_ncdf(ista)%ra_inc_obs
          END DO
        END IF
        

        !***********************************************************************************
        !    3) Fourth pass over the namelist to override rs_meta and dbz_meta
        !       with the values from the namelist:
        !
        !    3a) in order for this to work in case of shortended obs_times list,
        !        reset the list before reading the namelist. Later, this list of obs_times
        !        is matched against the list of obs_times_obs from the obs files in order
        !        to correctly transfer other time-dependent properties to the correct
        !        time index:
        DO ista=1, nradsta_max
          rs_meta(ista)%obs_times(:) = unused_value
        END DO
        !
        !
        CALL read_radarnamelist (izdom)
        !
        !***********************************************************************************

        !-----------------------------------------------------------------------------
        !    4) Debug output for the changed entries of the NetCDF meta data:
        !-----------------------------------------------------------------------------

        IF (ldebug_radsim) THEN
          DO i = 1, nradsta_changed
            WRITE (*,'(a,i6.6,a)') 'Overwriting (some) meta data of radar station ',rs_meta_ncdf(ista_ncdf(i))%station_id, &
                 ' with values from the namelist RADARSIM_PARAMS'
          END DO
        END IF

        !-----------------------------------------------------------------------------
        !   5) Merge rs_meta with rs_meta_ncdf: all stations from obs files which are
        !      not yet in rs_meta will be appended to the end of rs_meta:
        !-----------------------------------------------------------------------------

        ii = 0
        DO ista = 1, nradsta_ncdf
          found = .FALSE.
          DO i = 1, nradsta_namelist
            IF ( rs_meta_ncdf(ista)%station_id      == rs_meta(i)%station_id .AND. &
                 TRIM(rs_meta_ncdf(ista)%scanname)  == TRIM(rs_meta(i)%scanname) ) THEN
              ! rs_meta_ncdf(ista)%station_id is already in rs_meta-list and does not need to be added
              found = .TRUE.
              it = i
              EXIT
            END IF
          END DO
          IF (.NOT.found) THEN
            ! rs_meta_ncdf(ista)%station_id with rs_meta_ncdf(ista)%scanname is not in rs_meta-list, so append it at the end:
            IF (nradsta_namelist+ii < nradsta_max) THEN
              ii = ii + 1
              rs_meta(nradsta_namelist+ii)  = rs_meta_ncdf(ista)
              dbz_meta(nradsta_namelist+ii) = dbz_meta_ncdf(ista)
            ELSE
              errstring(:) = ' '
              WRITE (errstring,'(a,i0,a)') 'Too many radar stations (obs files + namelist)! ' // &
                   'Max allowed nradsta_max = ',nradsta_max,'! Discarding station ', &
                   rs_meta_ncdf(ista)%station_id, ' ', TRIM(rs_meta_ncdf(ista)%scanname)
              IF (labort_if_problems_obsfiles) THEN
                CALL abort_run (my_radar_id, 10876, &
                     'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                     'radar_namelist_read.f90, '//TRIM(yzroutine))
              ELSE
                WRITE(*,'(a)') 'ERROR '//TRIM(yzroutine)//': '//TRIM(errstring)
              END IF
            END IF
          ELSE
            ! rs_meta_ncdf(ista)%station_id with rs_meta_ncdf(ista)%scanname is in rs_meta-list.
            !  Check if obs_times have been explicitly specified after namelist re-reading,
            !  and if not, overtake them from the obs files:
            IF (ALL(rs_meta(it)%obs_times < miss_threshold) .AND. .NOT.rs_meta(it)%lobstimes_ovwrt_recalc) THEN
              rs_meta(it)%obs_times  = rs_meta_ncdf(ista)%obs_times
              rs_meta(it)%nobs_times = rs_meta_ncdf(ista)%nobs_times
            END IF
          END IF
        END DO
        
        ! .. Final number of stations in the merged list rs_meta:
        nradsta = nradsta_namelist + ii

        !----------------------------------------------------------------------------------------
        ! copy back voldata_ostream, because it has been overwritten by the above
        !  namelist reading_
        !----------------------------------------------------------------------------------------

        voldata_ostream = voldata_ostream_save

        !-----------------------------------------------------------------------------
        ! .. Create names for feedback files for each radar station:
        !-----------------------------------------------------------------------------

#ifdef NUDGING
        IF ( lfdbk_output ) THEN

          DO ista = 1, nradsta

            ! .. Name of feedback file incl. reference date:
            !     (here, the reference date is simply the model start time,
            !     but should perhaps be more precisely the verification reference time)
            rs_meta(ista)%fdbkfile(:) = ' '

            WRITE ( rs_meta(ista)%fdbkfile, '(a,i6.6,a)') &
                 'fof_radar_id-', rs_meta(ista)%station_id, '_' // TRIM(rs_meta(ista)%scanname) // &
                 '_' // ydate_ini_mod // '.nc'

            ! .. Initialize existence flag with .false., so that
            !    the file is newly created and an erroneously pre-existing
            !    file of the same name will be deleted (due to NF_CLOBBER mode in mo_t_netcdf_file.f90):
            rs_meta(ista)%lfdbkfile_exist = .FALSE.

          END DO

        END IF
#endif

      END IF   ! my_radar_id == 0

    ELSE       ! lreadmeta_from_netcdf

      nradsta = nradsta_namelist

    END IF     ! .not. lreadmeta_from_netcdf



    IF (my_radar_id == 0) THEN

      !-------------------------------------------------------------------------------
      !-------------------------------------------------------------------------------
      !
      ! .. Some other final (cross-) checks and computations/alterations on
      !    namelist parameters, which are possible only after the
      !    last call to read_radarnamelist()
      !
      !-------------------------------------------------------------------------------
      !-------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      ! .. check rs_meta(ista)%nobstimes:
      !-------------------------------------------------------------------------------

      DO ista = 1, nradsta

        IF (rs_meta(ista)%nobs_times > nobstimes_max) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a,i0,a,i6.6,a)') 'Namelist error: ' &
               //'rs_meta%nobstimes > nobstimes_max for station ',ista, ' (',rs_meta(ista)%station_id, ')!'
          CALL abort_run (my_radar_id, 10076, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF

      END DO

      !-------------------------------------------------------------------------------
      ! .. Check again that every radar has a valid station_id, otherwise problems
      !    might occur later:
      !-------------------------------------------------------------------------------

      DO ista = 1, nradsta
        IF (rs_meta(ista)%station_id == rs_meta_d(1)%station_id .OR. &
             rs_meta(ista)%station_id <= 0) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a,i0,a)') 'No station_id given for station ',ista,'! Please specify a number > 0!'
          rs_meta(ista)%nobs_times = 0
          rs_meta(ista)%obs_times  = unused_value
          IF (labort_if_problems_obsfiles) THEN
            CALL abort_run (my_radar_id, 10875, &
                 'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                 'radar_namelist_read.f90, '//TRIM(yzroutine))
          ELSE
            WRITE(*,'(a)') 'ERROR '//TRIM(yzroutine)//': '//TRIM(errstring)
          END IF
        END IF
      END DO

      !-------------------------------------------------------------------------------
      ! .. check if station IDs are all unique, otherwise naming of output files will
      !    not be unique:
      !-------------------------------------------------------------------------------

      DO ista = 1, nradsta
        DO i = ista+1, nradsta

          IF ( rs_meta(ista)%station_id     == rs_meta(i)%station_id .AND. &
               TRIM(rs_meta(ista)%scanname) == TRIM(rs_meta(i)%scanname) ) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a,i6.6,3a,2(i0,a))') &
                 'Namelist/inputfile error: station_id (',rs_meta(ista)%station_id,&
                 ') & scanname (',TRIM(rs_meta(ista)%scanname),&
                 ') not unique - identical for stations ',ista,' and ',i,' !'  
!!$ UB: should be enough to just cancel one of the two stations:
!!$            rs_meta(ista)%nobs_times = 0
!!$            rs_meta(ista)%obs_times  = unused_value
            rs_meta(i)%nobs_times = 0
            rs_meta(i)%obs_times  = unused_value
            IF (labort_if_problems_obsfiles) THEN
              CALL abort_run (my_radar_id, 10179, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, '//TRIM(yzroutine))
            ELSE
              WRITE (*,'(a)') 'ERROR '//TRIM(yzroutine)//': '//TRIM(errstring)
            END IF
          END IF
          
        END DO
      END DO

      !------------------------------------------------------------------------------
      ! .. check if radar composites from obs and sim are available for bubble trigger. If not, abbort.
      !------------------------------------------------------------------------------

      IF ( ldo_bubbles .AND. .NOT. (lreadmeta_from_netcdf .AND. loutdbz) ) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Namelist error: reflectivity observations not available for warm bubble generator, ' // &
                                'either set lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE., or ldo_bubbles=.FALSE.'
        CALL abort_run (my_radar_id, 10180, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      !------------------------------------------------------------------------------
      ! .. check if radar composites from obs and sim are available for bubble trigger. If not, abbort.
      !------------------------------------------------------------------------------

      IF ( ldo_composite .AND. .NOT. loutdbz ) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Namelist error: reflectivity not available for composite(s), ' // &
                                'either set loutdbz=.TRUE., or ldo_composite=.FALSE.'
        CALL abort_run (my_radar_id, 10180, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      !-------------------------------------------------------------------------------
      !- Section 4b: Calculate some necessary parameters
      !-------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      ! .. Specification of observation times, separately for "normal" output and
      !     for feedback and voldata files
      !
      !    First some preparations:
      !      - Fix the triplets of rs_meta%dt_obs, dt_obs_fdbk and dt_obs_voldata:
      !
      !        make sure that the triplet follows is the new notation <from>,<to>,<incr>
      !        and, if necessary, convert the traditional convention <incr>,<miss>,<miss> to the new
      !
      !      - global background list of obs times in feedback files
      !
      !-------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      !   .. Fix the triplets of rs_meta%dt_obs, dt_obs_fdbk and dt_obs_voldata:
      !        make sure that the triplet follows is the new notation <from>,<to>,<incr>
      !        and, if necessary, convert the traditional convention <incr>,<miss>,<miss> to the new
      !-------------------------------------------------------------------------------

      DO ista = 1, nradsta
        CALL check_and_fix_dtobs (rs_meta(ista)%dt_obs, errstring, ierr)
        IF (ierr /= 0) THEN
          cstation(:) = ' '
          WRITE (cstation, '(a,i0,a,i6.6,a)') ' for station ', ista, ' (',rs_meta(ista)%station_id, ')'
          CALL abort_run (my_radar_id, 10077, &
               'ERROR: problem in input_radarnamelist(): dt_obs '//TRIM(ADJUSTL(errstring))//TRIM(cstation), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
      END DO

      CALL check_and_fix_dtobs (dt_obs_fdbk_glob, errstring, ierr)
      IF (ierr /= 0) THEN
        CALL abort_run (my_radar_id, 10077, &
             'ERROR: problem in input_radarnamelist(): dt_obs_fdbk_glob '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      CALL check_and_fix_dtobs (dt_obs_voldata_glob, errstring, ierr)
      IF (ierr /= 0) THEN
        CALL abort_run (my_radar_id, 10077, &
             'ERROR: problem in input_radarnamelist(): dt_obs_voldata_glob '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      !-------------------------------------------------------------------------------
      !   .. global background list of obs times in feedback files
      !-------------------------------------------------------------------------------

      nobs_times_fdbk_glob = 0
      DO i=1, nobstimes_max
        IF (obs_times_fdbk_glob(i) >= 0.0_dp) THEN
          nobs_times_fdbk_glob = nobs_times_fdbk_glob + 1
          obs_times_fdbk_glob(nobs_times_fdbk_glob) = obs_times_fdbk_glob(i)
        END IF
      END DO

      IF (nobs_times_fdbk_glob == 0 .AND. dt_obs_fdbk_glob(3) > 0.0_dp) THEN

        nobs_times_fdbk_glob = missval_int ! will be set correctly by create_obstimes_from_dtobs
        CALL create_obstimes_from_dtobs (dt_obs_fdbk_glob, nobs_times_fdbk_glob, &
                                         obs_times_fdbk_glob, errstring, ierr)
        IF (ierr /= 0) THEN
          CALL abort_run (my_radar_id, 10079, &
               'ERROR: problem in input_radarnamelist(): dt_obs_fdbk_glob'//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
 
      END IF

      !-------------------------------------------------------------------------------
      !   .. global background list of obs times in volume data files
      !-------------------------------------------------------------------------------

      nobs_times_voldata_glob = 0
      DO i=1, nobstimes_max
        IF (obs_times_voldata_glob(i) >= 0.0_dp) THEN
          nobs_times_voldata_glob = nobs_times_voldata_glob + 1
          obs_times_voldata_glob(nobs_times_voldata_glob) = obs_times_voldata_glob(i)
        END IF
      END DO

      IF (nobs_times_voldata_glob == 0 .AND. dt_obs_voldata_glob(3) > 0.0_dp) THEN

        nobs_times_voldata_glob = missval_int ! will be set correctly by create_obstimes_from_dtobs
        CALL create_obstimes_from_dtobs (dt_obs_voldata_glob, nobs_times_voldata_glob, &
                                         obs_times_voldata_glob, errstring, ierr)
        IF (ierr /= 0) THEN
          CALL abort_run (my_radar_id, 10079, &
               'ERROR: problem in input_radarnamelist(): dt_obs_voldata_glob'//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
 
      END IF

      !-------------------------------------------------------------------------------
      ! .. now actually fill the observation times (and, if needed,
      !       obs_startrec, obs_endrec) for all radar stations
      !-------------------------------------------------------------------------------

      DO ista = 1, nradsta

        ! .. prepare info on radar station in error strings below:
        cstation(:) = ' '
        WRITE (cstation, '(a,i0,a,i6.6,a)') ' for station ', ista, ' (',rs_meta(ista)%station_id, ')'

        IF ( ALL(rs_meta(ista)%obs_times < miss_threshold) .AND. .NOT.lreadmeta_from_netcdf) THEN

          ! .. rs_meta(ista)%dt_obs is the relevant parameter here:
          CALL create_obstimes_from_dtobs (rs_meta(ista)%dt_obs, rs_meta(ista)%nobs_times, &
                                           rs_meta(ista)%obs_times, errstring, ierr)
          IF (ierr /= 0) THEN
            CALL abort_run (my_radar_id, 10078, &
                 'ERROR: problem in input_radarnamelist(): dt_obs'//TRIM(ADJUSTL(errstring))//TRIM(cstation), &
                 'radar_namelist_read.f90, input_radarnamelist()')
          END IF

        ELSE

          IF ( lreadmeta_from_netcdf .AND. ( &
               (rs_meta(ista)%lobstimes_ovwrt_recalc .AND. rs_meta(ista)%nobs_times /= missval_int) .OR. &
               (rs_meta(ista)%nobs_times_obs == missval_int .AND. ALL(rs_meta(ista)%obs_times < miss_threshold)) &
               ) ) THEN

            ! .. Either
            !     this is a station which is in the obs files, has not been eliminated after reading
            !     and has lobstimes_ovwrt_recalc=.TRUE.,
            !    OR
            !     this is a station which is not in the obs files and which has no obs_times from the namelist yet.
            !    In any case, the list of obs_times should be generated from rs_meta(ista)%dt_obs:
            
            ! .. Take care to get rid of previous valid values in the obs_times list:
            rs_meta(ista)%obs_times(:) = unused_value
            ! .. rs_meta(ista)%obs_times are re-calculated from dt_obs:
            CALL create_obstimes_from_dtobs (rs_meta(ista)%dt_obs, rs_meta(ista)%nobs_times, &
                                             rs_meta(ista)%obs_times, errstring, ierr)
            IF (ierr /= 0) THEN
              CALL abort_run (my_radar_id, 10078, &
                   'ERROR: problem in input_radarnamelist(): dt_obs'//TRIM(ADJUSTL(errstring))//TRIM(cstation), &
                   'radar_namelist_read.f90, input_radarnamelist()')
            END IF
          ELSE
            ! .. Set dt_obs to unused_value to indicate in YUSPECIF_RADAR that it has not been used:
            rs_meta(ista)%dt_obs(:) = unused_value
          END IF

          IF ( ANY(rs_meta(ista)%obs_times >= miss_threshold) ) THEN
          
            ! .. Obs times >= 0.0 have been specified in rs_meta(ista)%obs_times
            !    or have been read or re-defined explicitly from actual observation files
            !    or have been created by the specification via dt_obs and nobs_times.
            !    Now, we sort the obs times to the beginning of
            !    the vector and set nobstimes properly:
            !      (for the case of lreadmeta_from_netcdf=.true., we also
            !       reorder obs_cdate and obs_startrec, obs_endrec et al. correspondingly)

            rs_meta_buf(1)%obsfile        = rs_meta(ista)%obsfile
            rs_meta_buf(1)%obsfile_format = rs_meta(ista)%obsfile_format
            rs_meta_buf(1)%naz_ncdf       = rs_meta(ista)%naz_ncdf
            rs_meta_buf(1)%obs_times      = rs_meta(ista)%obs_times
            rs_meta_buf(1)%obs_cdate      = rs_meta(ista)%obs_cdate
            rs_meta_buf(1)%obs_startrec   = rs_meta(ista)%obs_startrec
            rs_meta_buf(1)%obs_endrec     = rs_meta(ista)%obs_endrec
            rs_meta_buf(1)%nobs_times     = rs_meta(ista)%nobs_times
            rs_meta_buf(1)%ext_nyq        = rs_meta(ista)%ext_nyq

            rs_meta(ista)%obsfile(:,:)(:) = ' '
            rs_meta(ista)%obsfile         = obsfile_missingname
            rs_meta(ista)%obsfile_format(:)(:) = ' '
            rs_meta(ista)%naz_ncdf        = missval_int
            rs_meta(ista)%obs_times       = unused_value
            rs_meta(ista)%obs_cdate       = 'YYYYMMDDHHMMSS'
            rs_meta(ista)%obs_startrec    = missval_int
            rs_meta(ista)%obs_endrec      = missval_int
            IF ( lreadmeta_from_netcdf ) THEN
              ! ext_nyq needs to be initialized with the correct default for the station_id/scanname:
              !  (Otherwise it already has the correct default for it's icountry, which is independent of obs_time)
              DO i = 1, nradsta_all
                IF ( rs_meta_for_country(i)%station_id      == rs_meta(ista)%station_id .AND. &
                     TRIM(rs_meta_for_country(i)%scanname)  == TRIM(rs_meta(ista)%scanname) ) THEN
                  rs_meta(ista)%ext_nyq = rs_meta_for_country(i)%ext_nyq
                END IF
              END DO
            END IF
            
            rs_meta(ista)%nobs_times      = 0

            DO i = 1, nobstimes_max
              IF ( check_obstime_within_modelrun(izdom, rs_meta_buf(1)%obs_times(i)) ) THEN
                
                it = -1
                IF (lreadmeta_from_netcdf) THEN
                  DO ii=1, nobstimes_max
                    IF ( rs_meta(ista)%obs_times_obs(ii) >= miss_threshold .AND. &
                         ABS(rs_meta(ista)%obs_times_obs(ii)-rs_meta_buf(1)%obs_times(i)) <= 1e-1_dp ) THEN
                      it = ii
                      EXIT
                    END IF
                  END DO
                END IF

                rs_meta(ista)%nobs_times                          = rs_meta(ista)%nobs_times + 1
                rs_meta(ista)%obs_times(rs_meta(ista)%nobs_times) = rs_meta_buf(1)%obs_times(i)

                IF (it > 0) THEN
                  rs_meta(ista)%obsfile_format(rs_meta(ista)%nobs_times) = rs_meta_buf(1)%obsfile_format(it)
                  rs_meta(ista)%obs_cdate   (rs_meta(ista)%nobs_times  ) = rs_meta_buf(1)%obs_cdate(it)
                  rs_meta(ista)%obs_startrec(rs_meta(ista)%nobs_times,:) = rs_meta_buf(1)%obs_startrec(it,:)
                  rs_meta(ista)%obs_endrec  (rs_meta(ista)%nobs_times,:) = rs_meta_buf(1)%obs_endrec(it,:)
                  rs_meta(ista)%obsfile     (rs_meta(ista)%nobs_times,:) = rs_meta_buf(1)%obsfile(it,:)
                  rs_meta(ista)%naz_ncdf    (rs_meta(ista)%nobs_times,:) = rs_meta_buf(1)%naz_ncdf(it,:)
                  rs_meta(ista)%ext_nyq   (:,rs_meta(ista)%nobs_times  ) = rs_meta_buf(1)%ext_nyq(:,it)
                END IF

              END IF
            END DO

          END IF
        END IF

        ! .. Set obs_cdate, if not already done:
        DO i = 1,rs_meta(ista)%nobs_times
          IF (rs_meta(ista)%obs_cdate(i) == 'YYYYMMDDHHMMSS') THEN
            rs_meta(ista)%obs_cdate(i) = new_datetime(ydate_ini_mod, rs_meta(ista)%obs_times(i))
          END IF
        END DO

        !-----------------------------------------------------------------------------
        ! .. Check if there are obs_times within the time of the model run,
        !    otherwise issue a warning:
        !-----------------------------------------------------------------------------

        IF ( rs_meta(ista)%nobs_times > 0 ) THEN
          IF ( .NOT. ANY(check_obstime_within_modelrun(izdom, rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times))) ) THEN

            WRITE (*,'(a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: no obs_times '// &
                 'within model run for station ',ista, ' (',rs_meta(ista)%station_id,')!'
          END IF
        ELSE
          IF (lreadmeta_from_netcdf) THEN
            WRITE (*,'(a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: no obs_times defined '// &
                 'for station ',ista, ' (',rs_meta(ista)%station_id,') or station discarded after meta data reading!'
          ELSE
            WRITE (*,'(a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: no obs_times defined '// &
                 'for station ',ista, ' (',rs_meta(ista)%station_id,')!'
          END IF
        END IF

        !-----------------------------------------------------------------------------
        ! .. Set obs_times_fdbk for the feedback files from namelist parameters
        !     rs_meta%dt_obs_fdbk or rs_meta%obs_times_fdbk, or if these are not
        !     explicitly specified in the namelist, overtake them from rs_meta%obs_times.
        !    Check against availability of obs for these times:
        !-----------------------------------------------------------------------------

        IF ( lreadmeta_from_netcdf .AND. rs_meta(ista)%nobs_times > 0 ) THEN

          IF ( ALL(rs_meta(ista)%obs_times_fdbk < miss_threshold) .AND. rs_meta(ista)%dt_obs_fdbk(3) > 0.0_dp ) THEN

            ! .. rs_meta(ista)%dt_obs_fdbk is the relevant parameter here, because no valid
            !     rs_meta(ista)%obs_times_fdbk have been specified in the namelist:

            CALL create_obstimes_from_dtobs (rs_meta(ista)%dt_obs_voldata, rs_meta(ista)%nobs_times_voldata, &
                                             rs_meta(ista)%obs_times_voldata, errstring, ierr)
            IF (ierr /= 0) THEN
              CALL abort_run (my_radar_id, 10079, &
                   'ERROR: problem in input_radarnamelist(): dt_obs_voldata'//TRIM(ADJUSTL(errstring))//TRIM(cstation), &
                   'radar_namelist_read.f90, input_radarnamelist()')
            END IF
 
          ELSE IF ( ANY(rs_meta(ista)%obs_times_fdbk >= miss_threshold) ) THEN

            ! .. rs_meta(ista)%obs_times_fdbk have been specified in the namelist, so
            !     we sort them to the beginning of the vector and filter for
            !     obs_times_fdbk during the model run:

            rs_meta_buf(1)%obs_times_fdbk   = rs_meta(ista)%obs_times_fdbk

            rs_meta(ista)%obs_times_fdbk    = unused_value
            rs_meta(ista)%nobs_times_fdbk   = 0
            DO i = 1, nobstimes_max
              IF ( check_obstime_within_modelrun(izdom, rs_meta_buf(1)%obs_times_fdbk(i)) ) THEN
                rs_meta(ista)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk + 1
                rs_meta(ista)%obs_times_fdbk(rs_meta(ista)%nobs_times_fdbk) = &
                     rs_meta_buf(1)%obs_times_fdbk(i)
              END IF
            END DO

            ! .. Set dt_obs_fdbk to unused_value to indicate in YUSPECIF_RADAR that it has not been used:
            rs_meta(ista)%dt_obs_fdbk(:) = unused_value

          ELSE IF (nobs_times_fdbk_glob > 0) THEN

            ! .. global obs_times_fdbk_glob or dt_obs_fdbk_glob have been specified in the namelist, so
            !    apply them now to the actual radar station:

            rs_meta(ista)%obs_times_fdbk    = unused_value
            rs_meta(ista)%nobs_times_fdbk = 0
            DO i = 1, nobs_times_fdbk_glob
              IF ( check_obstime_within_modelrun(izdom, obs_times_fdbk_glob(i)) ) THEN
                rs_meta(ista)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk + 1
                rs_meta(ista)%obs_times_fdbk(rs_meta(ista)%nobs_times_fdbk) = obs_times_fdbk_glob(i)
              END IF
            END DO

          ELSE

            ! .. No extra specification of obs_times_fdbk or dt_obs_fdbk in the namelist,
            !     so we just overtake all the "normal" obs_times to the feedback files:

            rs_meta(ista)%nobs_times_fdbk = rs_meta(ista)%nobs_times
            rs_meta(ista)%obs_times_fdbk  = rs_meta(ista)%obs_times
            rs_meta(ista)%dt_obs_fdbk     = rs_meta(ista)%dt_obs

          END IF

          ! .. Check if the obs_times_fdbk have a corresponding obs_times and eliminate them if not.
          !     Otherwise, there will not be corresponding output in the feedback files:

          rs_meta_buf(1)%obs_times_fdbk  = rs_meta(ista)%obs_times_fdbk
          rs_meta_buf(1)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk
          rs_meta(ista)%obs_times_fdbk   = unused_value
          rs_meta(ista)%nobs_times_fdbk  = 0
          DO i = 1, rs_meta_buf(1)%nobs_times_fdbk
            IF ( ALL( ABS( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) - &
                           rs_meta_buf(1)%obs_times_fdbk(i) ) > 1e-1_dp ) ) THEN
              WRITE (*,'(a,f10.1,a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: obs_times_fdbk = ', &
                   rs_meta_buf(1)%obs_times_fdbk(i), &
                   ' not available for station ',ista, ' (',rs_meta(ista)%station_id,')!'
            ELSE
              rs_meta(ista)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk + 1
              rs_meta(ista)%obs_times_fdbk(rs_meta(ista)%nobs_times_fdbk) = rs_meta_buf(1)%obs_times_fdbk(i)
            END IF
          END DO

          ! .. Check if the remaining obs_times_fdbk really have corresponding obs files and eliminate them if not.
          !     Otherwise, the corresponding feedback files would be unnecessarily large:
          rs_meta_buf(1)%obs_times_fdbk  = rs_meta(ista)%obs_times_fdbk
          rs_meta_buf(1)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk
          rs_meta(ista)%obs_times_fdbk   = unused_value
          rs_meta(ista)%nobs_times_fdbk  = 0
          DO i = 1, rs_meta_buf(1)%nobs_times_fdbk
            search_loop: DO ii = 1, rs_meta(ista)%nobs_times
              IF ( ABS( rs_meta(ista)%obs_times(ii) - rs_meta_buf(1)%obs_times_fdbk(i) ) <= 1e-1_dp ) THEN
                IF ( ALL( rs_meta(ista)%obsfile(ii,:) == obsfile_missingname) ) THEN
                  WRITE (*,'(a,f10.1,a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: no obs file for obs_times_fdbk = ', &
                       rs_meta_buf(1)%obs_times_fdbk(i), &
                       ' for station ',ista, ' (',rs_meta(ista)%station_id,')!'
                ELSE
                  rs_meta(ista)%nobs_times_fdbk = rs_meta(ista)%nobs_times_fdbk + 1
                  rs_meta(ista)%obs_times_fdbk(rs_meta(ista)%nobs_times_fdbk) = rs_meta_buf(1)%obs_times_fdbk(i)
                END IF
                EXIT search_loop
              END IF
            END DO search_loop
          END DO
          
        END IF

        !-----------------------------------------------------------------------------
        ! .. Similarly for obs_times_voldata:
        !-----------------------------------------------------------------------------

        IF ( lvoldata_output .AND. rs_meta(ista)%nobs_times > 0 ) THEN

          IF ( ALL(rs_meta(ista)%obs_times_voldata < miss_threshold) .AND. rs_meta(ista)%dt_obs_voldata(3) > 0.0_dp ) THEN

            ! .. rs_meta(ista)%dt_obs_voldata is the relevant parameter here, because no valid
            !     rs_meta(ista)%obs_times_voldata have been specified in the namelist:

            CALL create_obstimes_from_dtobs (rs_meta(ista)%dt_obs_voldata, rs_meta(ista)%nobs_times_voldata, &
                                             rs_meta(ista)%obs_times_voldata, errstring, ierr)
            IF (ierr /= 0) THEN
              CALL abort_run (my_radar_id, 10079, &
                   'ERROR: problem in input_radarnamelist(): dt_obs_voldata'//TRIM(ADJUSTL(errstring))//TRIM(cstation), &
                   'radar_namelist_read.f90, input_radarnamelist()')
            END IF
 
          ELSE IF ( ANY(rs_meta(ista)%obs_times_voldata >= miss_threshold) ) THEN

            ! .. rs_meta(ista)%obs_times_voldata have been specified in the namelist, so
            !     we sort them to the beginning of the vector and filter for
            !     obs_times_voldata during the model run:

            rs_meta_buf(1)%obs_times_voldata   = rs_meta(ista)%obs_times_voldata

            rs_meta(ista)%obs_times_voldata    = unused_value
            rs_meta(ista)%nobs_times_voldata   = 0
            DO i = 1, nobstimes_max
              IF ( check_obstime_within_modelrun(izdom, rs_meta_buf(1)%obs_times_voldata(i)) ) THEN
                rs_meta(ista)%nobs_times_voldata = rs_meta(ista)%nobs_times_voldata + 1
                rs_meta(ista)%obs_times_voldata(rs_meta(ista)%nobs_times_voldata) = &
                     rs_meta_buf(1)%obs_times_voldata(i)
              END IF
            END DO

            ! .. Set dt_obs_voldata to unused_value to indicate in YUSPECIF_RADAR that it has not been used:
            rs_meta(ista)%dt_obs_voldata(:) = unused_value

          ELSE IF (nobs_times_voldata_glob > 0) THEN

            ! .. global obs_times_voldata_glob or dt_obs_voldata_glob have been specified in the namelist, so
            !    apply them now to the actual radar station:

            rs_meta(ista)%obs_times_voldata    = unused_value
            rs_meta(ista)%nobs_times_voldata = 0
            DO i = 1, nobs_times_voldata_glob
              IF ( check_obstime_within_modelrun(izdom, obs_times_voldata_glob(i)) ) THEN
                rs_meta(ista)%nobs_times_voldata = rs_meta(ista)%nobs_times_voldata + 1
                rs_meta(ista)%obs_times_voldata(rs_meta(ista)%nobs_times_voldata) = obs_times_voldata_glob(i)
              END IF
            END DO

          ELSE

            ! .. No extra specification of obs_times_voldata or dt_obs_voldata in the namelist,
            !     so we just overtake all the "normal" obs_times to the volume data  files:

            rs_meta(ista)%nobs_times_voldata = rs_meta(ista)%nobs_times
            rs_meta(ista)%obs_times_voldata  = rs_meta(ista)%obs_times
            rs_meta(ista)%dt_obs_voldata     = rs_meta(ista)%dt_obs

          END IF

          ! .. Check if the obs_times_voldata have a corresponding obs_times and eliminate them if not.
          !     Otherwise, there will not be corresponding output in the volume data files:

          rs_meta_buf(1)%obs_times_voldata  = rs_meta(ista)%obs_times_voldata
          rs_meta_buf(1)%nobs_times_voldata = rs_meta(ista)%nobs_times_voldata
          rs_meta(ista)%obs_times_voldata   = unused_value
          rs_meta(ista)%nobs_times_voldata  = 0
          DO i = 1, rs_meta_buf(1)%nobs_times_voldata
            IF ( ALL( ABS( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) - &
                           rs_meta_buf(1)%obs_times_voldata(i) ) > 1e-1_dp ) ) THEN
              WRITE (*,'(a,f10.1,a,i0,a,i6.6,a)') '  WARNING input_radarnamelist: obs_times_voldata = ', &
                   rs_meta_buf(1)%obs_times_voldata(i), &
                   ' not available for station ',ista, ' (',rs_meta(ista)%station_id,')!'
            ELSE
              rs_meta(ista)%nobs_times_voldata = rs_meta(ista)%nobs_times_voldata + 1
              rs_meta(ista)%obs_times_voldata(rs_meta(ista)%nobs_times_voldata) = rs_meta_buf(1)%obs_times_voldata(i)
            END IF
          END DO

        END IF

      END DO   ! nradsta


      !-----------------------------------------------------------------------------
      ! .. Overtake the radar wavelength and station_id from rs_meta - structure to
      !    dbz_meta - structure:
      !-----------------------------------------------------------------------------

      DO ista = 1, nradsta

        dbz_meta(ista)%lambda_radar = rs_meta(ista)%lambda
        dbz_meta(ista)%station_id   = rs_meta(ista)%station_id

      END DO

      !-----------------------------------------------------------------------------
      ! .. Compute rs_meta(ista)%alpha3_eff_0 from interpolated table-lookup
      !      (values from Blahak, 2008, JAOTECH)
      !-----------------------------------------------------------------------------

      DO ista = 1, nradsta

        rs_meta(ista)%alpha3_eff_0 = get_alpha3_eff_0 ( &
             rs_meta(ista)%dalpha, rs_meta(ista)%Phi3 )

      END DO

      !-------------------------------------------------------------------------------
      ! .. Set rs_meta(ista)%ind_ele_fdbk and rs_meta(ista)%nel_fdbk correctly,
      !    which are the elevations to write into the feedback file.
      !    Depends on rs_meta(ista)%ind_ele_fdbk, thin_step_ele and ind_ele_fdbk_glob
      !    Added by EB 05/2017 (second test)
      !-------------------------------------------------------------------------------

      ! .. First, prepare ind_ele_fdbk_glob: sort valid entries to the front of the vector and count
      !     the valid entries. Later, the corresponding values in the rs_meta - structure
      !     are overwritten with the _glob-values if nel_fdbk_glob > 0:
      nel_fdbk_glob = 0
      DO i=1, nel_max
        IF (ind_ele_fdbk_glob(i) > 0) THEN
          nel_fdbk_glob = nel_fdbk_glob + 1
          ind_ele_fdbk_glob(nel_fdbk_glob) = ind_ele_fdbk_glob(i)
        END IF
      END DO

      ! .. Then, fill rs_meta(ista)%ind_ele_fdbk with the correct values, depending
      !     on previously readings from the namelist:
      DO ista = 1, nradsta
        IF ( ANY(rs_meta(ista)%ind_ele_fdbk > 0 .AND. rs_meta(ista)%ind_ele_fdbk <= rs_meta(ista)%nel) ) THEN
          ! Valid values for individual station ista have been specified in the namelist,
          !  so overtake these valid values:
          rs_meta(ista)%nel_fdbk = 0
          DO i=1, nel_max
            IF (rs_meta(ista)%ind_ele_fdbk(i) > 0 .AND. rs_meta(ista)%ind_ele_fdbk(i) <= rs_meta(ista)%nel) THEN
              rs_meta(ista)%nel_fdbk                             = rs_meta(ista)%nel_fdbk + 1
              rs_meta(ista)%ind_ele_fdbk(rs_meta(ista)%nel_fdbk) = rs_meta(ista)%ind_ele_fdbk(i)
            END IF
          END DO
        ELSE
          IF (nel_fdbk_glob > 0) THEN
            rs_meta(ista)%nel_fdbk = 0
            DO i=1, nel_fdbk_glob
              IF (ind_ele_fdbk_glob(i) <= rs_meta(ista)%nel) THEN
                rs_meta(ista)%nel_fdbk                             = rs_meta(ista)%nel_fdbk + 1
                rs_meta(ista)%ind_ele_fdbk(rs_meta(ista)%nel_fdbk) = ind_ele_fdbk_glob(i)
              END IF
            END DO
          ELSE
            rs_meta(ista)%nel_fdbk                          = rs_meta(ista)%nel
            rs_meta(ista)%ind_ele_fdbk(1:rs_meta(ista)%nel) = (/ (i, i=1, rs_meta(ista)%nel) /)
          END IF
        END IF
      END DO

      ! .. Last step: elevation thinning
      IF (thin_step_ele > 1) THEN
        DO ista = 1, nradsta
          nel_fdbk_tmp = rs_meta(ista)%nel_fdbk
          rs_meta(ista)%nel_fdbk = 0
          DO i=1, nel_fdbk_tmp, thin_step_ele
            rs_meta(ista)%nel_fdbk                             = rs_meta(ista)%nel_fdbk + 1
            rs_meta(ista)%ind_ele_fdbk(rs_meta(ista)%nel_fdbk) = rs_meta(ista)%ind_ele_fdbk(i)
          END DO
        END DO
      END IF

      !-------------------------------------------------------------------------------
      ! .. Set rs_meta(ista)%ind_ele_voldata and rs_meta(ista)%nel_voldata correctly,
      !    which are the elevations to write into the volume data files.
      !    Depends on rs_meta(ista)%ind_ele_voldata and ind_ele_voldata_glob
      !-------------------------------------------------------------------------------

      ! .. First, prepare ind_ele_voldata_glob: sort valid entries to the front of the vector and count
      !     the valid entries. Later, the corresponding values in the rs_meta - structure
      !     are overwritten with the _glob-values if nel_voldata_glob > 0:
      nel_voldata_glob = 0
      DO i=1, nel_max
        IF (ind_ele_voldata_glob(i) > 0) THEN
          nel_voldata_glob = nel_voldata_glob + 1
          ind_ele_voldata_glob(nel_voldata_glob) = ind_ele_voldata_glob(i)
        END IF
      END DO

      ! .. Then, fill rs_meta(ista)%ind_ele_voldata with the correct values, depending
      !     on previously readings from the namelist:
      DO ista = 1, nradsta
        IF ( ANY(rs_meta(ista)%ind_ele_voldata > 0 .AND. rs_meta(ista)%ind_ele_voldata <= rs_meta(ista)%nel) ) THEN
          ! Valid values for individual station ista have been specified in the namelist,
          !  so overtake these valid values:
          rs_meta(ista)%nel_voldata = 0
          DO i=1, nel_max
            IF (rs_meta(ista)%ind_ele_voldata(i) > 0 .AND. rs_meta(ista)%ind_ele_voldata(i) <= rs_meta(ista)%nel) THEN
              rs_meta(ista)%nel_voldata                                = rs_meta(ista)%nel_voldata + 1
              rs_meta(ista)%ind_ele_voldata(rs_meta(ista)%nel_voldata) = rs_meta(ista)%ind_ele_voldata(i)
            END IF
          END DO
        ELSE
          IF (nel_voldata_glob > 0) THEN
            rs_meta(ista)%nel_voldata = 0
            DO i=1, nel_voldata_glob
              IF (ind_ele_voldata_glob(i) <= rs_meta(ista)%nel) THEN
                rs_meta(ista)%nel_voldata                             = rs_meta(ista)%nel_voldata + 1
                rs_meta(ista)%ind_ele_voldata(rs_meta(ista)%nel_voldata) = ind_ele_voldata_glob(i)
              END IF
            END DO
          ELSE
            rs_meta(ista)%nel_voldata                          = rs_meta(ista)%nel
            rs_meta(ista)%ind_ele_voldata(1:rs_meta(ista)%nel) = (/ (i, i=1, rs_meta(ista)%nel) /)
          END IF
        END IF
      END DO


      !-------------------------------------------------------------------------------
      ! .. Make sure that rs_meta(ista)%eleind_for_composite_bub and
      !    rs_meta(ista)%eleindlist_for_composite are within the allowed range
      !-------------------------------------------------------------------------------

      IF (ldo_composite) THEN
        DO i=1, nel_composite
          DO ista = 1, nradsta
            ! .. If rs_meta(ista)%eleind_for_composite_bub = 98, it will be constructed
            !     from the precipitation scans. If no precip scans are in the data, composite will be empty.
            ! .. If rs_meta(ista)%eleindlist_for_composite(i) = 99, it will be a vertical
            !     maximum composite taken over all elevations.
            !    Otherwise, it will be constructed from one specific elevation:
            IF (rs_meta(ista)%eleindlist_for_composite(i) /= 98 .AND. rs_meta(ista)%eleindlist_for_composite(i) /= 99) THEN
              rs_meta(ista)%eleindlist_for_composite(i) = &
                   MAX( MIN( rs_meta(ista)%eleindlist_for_composite(i), rs_meta(ista)%nel), 1)
            END IF
          END DO
        END DO
      END IF

      IF (ldo_bubbles) THEN
        DO ista = 1, nradsta
          ! .. If rs_meta(ista)%eleind_for_composite_bub = 98, it will be constructed
          !     from the precipitation scans. If no precip scans are in the data, composite will be empty.
          ! .. If rs_meta(ista)%eleind_for_composite_bub = 99, it will be a vertical
          !     maximum composite taken over all elevations.
          !    Otherwise, it will be constructed from one specific elevation:
          IF (rs_meta(ista)%eleind_for_composite_bub /= 98 .AND. rs_meta(ista)%eleind_for_composite_bub /= 99) THEN
            rs_meta(ista)%eleind_for_composite_bub = &
                 MAX( MIN(rs_meta(ista)%eleind_for_composite_bub,  rs_meta(ista)%nel), 1)
          END IF
        END DO
      END IF

      !-------------------------------------------------------------------------------
      ! .. Prepare range coarsening in obs-data mode:
      !
      !    Set rs_meta(i)%ra_inc, rs_meta(i)%nra and rs_meta(i)%n_aggr_ra_obs
      !    for range bin coarsening in obs-data mode:
      !
      !  if (lreadmeta_from_netcdf) then
      !
      !    The range resolution and number of range bins which are actually used in EMVORADO
      !    are the values in rs_meta(i)%ra_inc, rs_meta(i)%nra.
      !    The goal here is to determine a coarsened range resolution
      !    rs_meta(i)%ra_inc as the nearest interger multiple to the user-chosen
      !    rs_meta(i)%ra_inc_coarse  of the observation's range resolution rs_meta(i)%ra_inc_obs,
      !    if rs_meta(i)%ra_inc_coarse > rs_meta(i)%ra_inc_obs.
      !
      !    Otherwise, rs_meta(i)%ra_inc and rs_meta(i)%nra remain unchanged from their
      !    original values, which should be equal at this point to rs_meta(i)%ra_inc_obs
      !    and rs_meta(i)%nra_obs. To ensure this also for cases when actual data
      !    for a specific station are missing at the actual time, we explicitly set
      !    rs_meta(i)%ra_inc and rs_meta(i)%nra to their desired _obs values.
      !
      !  endif
      !
      !-------------------------------------------------------------------------------

      IF (lreadmeta_from_netcdf) THEN
        DO ista = 1, nradsta
          IF ( rs_meta(ista)%ra_inc_coarse > rs_meta(ista)%ra_inc_obs + eps ) THEN
            rs_meta(ista)%n_aggr_ra_obs = MAX(NINT(rs_meta(ista)%ra_inc_coarse/rs_meta(ista)%ra_inc_obs),1)
            rs_meta(ista)%ra_inc = rs_meta(ista)%n_aggr_ra_obs * rs_meta(ista)%ra_inc_obs
            rs_meta(ista)%nra = CEILING(rs_meta(ista)%ra_inc_obs*rs_meta(ista)%nra_obs / rs_meta(ista)%ra_inc)
            ! .. Reduce the coarse rs_meta(ista)%nra if the last coarse range bin would not be entirely filled
            !     by original range bins until its outer boundary:
!!$            rs_meta(ista)%nra = rs_meta(ista)%nra - &
!!$                 CEILING( ( (rs_meta(ista)%nra+0.5_dp)*rs_meta(ista)%ra_inc - &
!!$                            rs_meta(ista)%ra_inc_obs*rs_meta(ista)%nra_obs ) / rs_meta(ista)%ra_inc )
          ELSE
            ! .. for safety in case of missing actual data for a station or for a
            !     wrong explicit setting of ra_inc and nra in the namelist by the user:
            rs_meta(ista)%ra_inc        = rs_meta(ista)%ra_inc_obs
            rs_meta(ista)%ra_inc_coarse = rs_meta(ista)%ra_inc_obs
            rs_meta(ista)%nra           = rs_meta(ista)%nra_obs
            rs_meta(ista)%n_aggr_ra_obs = 1
          END IF
        END DO
      ELSE
        ! .. rs_meta(ista)%ra_inc_coarse may also be used to re-compute ra_inc in case if the user
        !     wants to exactly reproduce what happens to a real radar when using ra_inc_coarse
        !     in the namelist. This might be useful in forecast runs without using obs data,
        !     which are intended to emulate the exact behaviour of forecast runs with using obs data
        !     and range coarsening:
        DO ista = 1, nradsta
          IF ( rs_meta(ista)%ra_inc_coarse > rs_meta(ista)%ra_inc_obs + eps ) THEN
            ! .. ra_inc_coarse has been specified in the namelist to a larger value than ra_inc_obs.
            !    Range-coarsening based on ra_inc_obs and nra_obs, coming either
            !     from the NAMELIST or from the default for this station:
            rs_meta(ista)%n_aggr_ra_obs = MAX(NINT(rs_meta(ista)%ra_inc_coarse/rs_meta(ista)%ra_inc_obs),1)
            rs_meta(ista)%ra_inc = rs_meta(ista)%n_aggr_ra_obs * rs_meta(ista)%ra_inc_obs
            rs_meta(ista)%nra = CEILING(rs_meta(ista)%ra_inc_obs*rs_meta(ista)%nra_obs / rs_meta(ista)%ra_inc)
            ! .. Reduce the coarse rs_meta(ista)%nra if the last coarse range bin would not be entirely filled
            !     by original range bins until its outer boundary:
!!$            rs_meta(ista)%nra = rs_meta(ista)%nra - &
!!$                 CEILING( ( (rs_meta(ista)%nra+0.5_dp)*rs_meta(ista)%ra_inc - &
!!$                            rs_meta(ista)%ra_inc_obs*rs_meta(ista)%nra_obs ) / rs_meta(ista)%ra_inc )
          ELSE
            ! .. ra_inc_coarse has not been specified in the namelist, so ra_inc is used directly and
            !     takes precedence over ra_inc_obs:
            rs_meta(ista)%ra_inc_obs    = rs_meta(ista)%ra_inc
            rs_meta(ista)%ra_inc_coarse = rs_meta(ista)%ra_inc
            rs_meta(ista)%nra_obs       = rs_meta(ista)%nra
            rs_meta(ista)%n_aggr_ra_obs = 1
          END IF
        END DO
      END IF

    END IF  ! my_radar_id == 0

    !-------------------------------------------------------------------------------
    !- Section 5: Part 2:
    !             Distribute station meta data to all other nodes
    !-------------------------------------------------------------------------------

    IF (num_radar > 1) THEN

      ! First, re-distribute nradsta to all other nodes because some stations
      !  could have been eliminated from the list because of:
      !  - errors in their meta data or scan strategies if lreadmeta_from_netcdf=.true.
      !  - station position is outside the model domain

      CALL distribute_values_radar (nradsta, 1, 0, icomm_radar, ierr)

      ! loop over all radar stations
      DO ista = 1, nradsta

        ! distribute the rs_meta structure and the dbz_meta structure:

        CALL mpi_bcast (rs_meta(ista), 1, mpi_radar_meta_typ_alltimes, 0,icomm_radar,ierr)
        CALL mpi_bcast (dbz_meta(ista),1, mpi_dbzcalc_params_typ     , 0,icomm_radar,ierr)

      END DO

      ! Just for completeness: distribute the global ind_ele_XXX_glob:
      CALL distribute_values_radar (ind_ele_fdbk_glob   , nel_max, 0, icomm_radar, ierr)
      CALL distribute_values_radar (ind_ele_voldata_glob, nel_max, 0, icomm_radar, ierr)

    END IF


#ifdef __ICON__
    !-------------------------------------------------------------------------------
    ! .. Set up the environment in radar_interface.f90 to find the nearest
    !    grid cell to a given lon/lat coordinate. This is the basis for all
    !    ICON grid related operations in the following.
    !-------------------------------------------------------------------------------

    IF (lcompute_pe_fwo .AND. nradsta > 0) THEN
      CALL setup_auxgrid_for_cellindex (izdom, nradsta, rs_meta(1:nradsta), ldebug_radsim)
    END IF
#endif

    !-------------------------------------------------------------------------------
    ! .. Eliminate stations outside the model domain from the station list
    !    in case of online beam propagation,
    !    by setting the number of obs_times to 0. During a later check, such stations
    !    will be eliminated from the station list and nradsta reduced
    !    accordingly:
    !-------------------------------------------------------------------------------

    IF (lonline) THEN

      DO ista = 1,nradsta

        ! Check if radar station is outside model domain:

        CALL geo2model_coord_domaincheck (izdom, rs_meta(ista)%lon, rs_meta(ista)%lat, rlon_r, rlat_r, &
                                          is_inside, rlon_min, rlon_max, rlat_min, rlat_max)

        IF ( .NOT. is_inside .AND. my_radar_id == 0 ) THEN
#ifndef __ICON__
          WRITE (*,'(a,i0,a,i6.6,a,2(/,a,f12.5,a,2f14.5,a),/,a)') &
               'WARNING input_radarnamelist(): Radar Nr. ', ista, ' (',rs_meta(ista)%station_id, &
               ') position is outside of model domain,','  lonrot_rad = ', &
               rlon_r, ' ( domain min/max ', rlon_min, rlon_max, ' )','  latrot_rad = ', &
               rlat_r, ' ( domain min/max ', rlat_min, rlat_max, ' ) ', &
               '  Station will be discarded from the forward operator, because lonline=.TRUE.!'
#else
          WRITE (*,'(a,i0,a,i6.6,a,2(/,a,f12.5),/,a)') &
               'WARNING input_radarnamelist(): Radar Nr. ', ista, ' (',rs_meta(ista)%station_id, &
               ') position is outside of model domain! '  , &
               '  lonrot_rad = ', rlon_r, '  latrot_rad = ', rlat_r, &
               '  Station will be discarded from the forward operator, because lonline=.TRUE.!'
#endif
          rs_meta(ista)%nobs_times = 0
          rs_meta(ista)%obs_times  = unused_value

        END IF

      END DO

    END IF

    !---------------------------------------------------------------------------------------
    ! .. Eliminate those stations from the station list where the number of obs_times is 0.
    !    nradsta is reduced accordingly and later re-distributed to all processors:
    !---------------------------------------------------------------------------------------

    nradsta_valid = 0
    DO ista = 1, nradsta
      IF (rs_meta(ista)%nobs_times > 0) THEN
        nradsta_valid = nradsta_valid + 1
        rs_meta(nradsta_valid) = rs_meta(ista)
        dbz_meta(nradsta_valid) = dbz_meta(ista)
      END IF
    END DO
    DO ista = nradsta_valid+1, nradsta
      rs_meta(ista)  = rs_meta_d(1)
      dbz_meta(ista) = dbz_meta_d(1)
    END DO
    nradsta = nradsta_valid


    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !- Section 6: Do some checks which can only be done after the
    !             distribution of the namelist parameters to all nodes
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    ! .. Find processor(s) on whose subdomain(s) the radar station is located,
    !    get orography height at the station and
    !    check station height against orography:
    !-------------------------------------------------------------------------------

    DO ista = 1, nradsta

      ! calculate rotated longitude and latitude of radar position and the
      !  global index for the nearest model grid point of the radar station
      !  for meta information:

      IF (lcompute_pe_fwo) THEN

        CALL geo2model_cellindex (izdom, rs_meta(ista)%lon, rs_meta(ista)%lat, &
                                  rs_meta(ista)%i_nearest_mod, rs_meta(ista)%j_nearest_mod)

        ! Interpolate the orography height at the radar station and broadcast to all PEs:

        CALL interp2D_model2geo_horiz_scalar(izdom, hhl, bottomlevel_stag(), rs_meta(ista)%lon, rs_meta(ista)%lat, &
                                             hsurf_r, is_inside)

        IF ( is_inside ) THEN
          IF ( ldebug_radsim ) THEN
            WRITE (*,'(a,i0,a,i6.6,a,i6)') 'INFO from '//TRIM(yzroutine)// &
                 ': position of radar station ',ista, ' (', rs_meta(ista)%station_id, ') has been found on proc ', my_radar_id
          END IF
        ELSE
          rs_meta(ista)%i_nearest_mod = -9999
          rs_meta(ista)%j_nearest_mod = -9999
        END IF

        ! Distribute hsurf_r by means of a mpi_reduce "MAX" - operation. By using this
        ! trick, the situation is correctly handled that a radar station is
        ! positionend in the boundary zone of a processor region and hsurf_r
        ! has been computed above on more than 1 processor:

      ELSE
        hsurf_r = -9999.99_dp
        is_inside = .FALSE.
        rs_meta(ista)%i_nearest_mod = -9999
        rs_meta(ista)%j_nearest_mod = -9999
      END IF

      IF (num_radar > 1) THEN
        errstring(:) = ' '
        CALL global_values_radar(rs_meta(ista)%i_nearest_mod, "MAX", icomm_radar, -1, errstring, ierr)
        CALL global_values_radar(rs_meta(ista)%j_nearest_mod, "MAX", icomm_radar, -1, errstring, ierr)
        errstring(:) = ' '
        CALL global_values_radar(hsurf_r, "MAX", icomm_radar, -1, errstring, ierr)
      END IF

      IF (my_radar_id == 0) THEN
        IF (hsurf_r < -9000_dp .AND. lonline) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a,i0,a,i6.6,a)') 'ERROR for lonline=.TRUE.: '// &
               'orography height at the station site for station ', ista, ' (', rs_meta(ista)%station_id, &
               ') could not be determined.'
          CALL abort_run (my_radar_id, 10374, &
               'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
               'radar_namelist_read.f90, input_radarnamelist()')
        END IF
      END IF


      !-------------------------------------------------------------------------------
      ! .. Set rs_meta(ista)%alt_msl according to namelist settings,
      !    rs_meta(ista)%alt_msl takes precedence over rs_meta(ista)%alt_agl_mod:

      !   - rs_meta(ista)%alt_agl_mod :  height AGL above the model orography
      !   - rs_meta(ista)%alt_msl     :  height above MSL
      !
      ! IF only rs_meta(ista)%alt_agl_mod or nothing is given in the namelist and
      !  rs_meta(ista)%alt_msl < -9000.0 (default), then the radar height
      !  rs_meta(ista)%alt_msl will be set to hsurf(iradar,jradar) + rs_meta(ista)%alt_agl_mod.
      !  The default for rs_meta(ista)%alt_agl_mod is 50 m.
      !
      ! IF rs_meta(ista)%alt_msl >= -9000.0, this value will be used
      !  directly and takes precedence over the other, because all relevant subroutines
      !  of the radar simulator work with rs_meta(ista)%alt_msl.
      !
      ! IF the radar station is outside the domain and no rs_meta(ista)%alt_msl
      !  is given in the namelist, rs_meta(ista)%alt_msl_true is used
      !  instead. alt_msl_true is available if lreadmeta_from_netcdf=.true. or can be
      !  specified in the namelist.
      !-------------------------------------------------------------------------------

      rs_meta(ista)%msl_mod = hsurf_r

      IF ( rs_meta(ista)%alt_msl < -9000.0_dp ) THEN

        IF (rs_meta(ista)%msl_mod < -9000.0_dp) THEN

          ! in this case the radar station is outside the model domain and
          !  rs_meta(ista)%alt_msl is not given explicitly in the namelist.

          IF (rs_meta(ista)%alt_msl_true > -9000.0_dp) THEN
            ! in this case, lreadmeta_from_netcdf=.true. and the true station
            !  height has been read from the radar file or from the background
            !  meta data list and can be used:
            rs_meta(ista)%alt_msl = rs_meta(ista)%alt_msl_true
          ELSE
            IF (my_radar_id == 0) THEN
              errstring(:) = ' '
              WRITE (errstring,'(a,i0,a,i6.6,a,i3.3,a)') 'Radar Nr. ', ista, ' (', rs_meta(ista)%station_id, &
                   ') outside model domain and namelist parameter'// &
                   ' rs_meta(',ista,')%alt_msl for absolute antenna height MSL is missing in namelist!'
              CALL abort_run (my_radar_id, 10375, &
                   'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
                   'radar_namelist_read.f90, input_radarnamelist()')
            END IF
          END IF
        ELSE
          rs_meta(ista)%alt_msl = rs_meta(ista)%msl_mod + rs_meta(ista)%alt_agl_mod
        END IF

      ELSE

        ! unset rs_meta(ista)%alt_agl_mod to show explicitly that this
        ! value is not used in the simulator:

        rs_meta(ista)%alt_agl_mod = -9999.99

      END IF

      !-------------------------------------------------------------------------------
      ! .. Final check of station height:
      !-------------------------------------------------------------------------------

      IF ( rs_meta(ista)%alt_msl <= hsurf_r ) THEN

        errstring(:) = ' '
        WRITE (errstring,'(a,i0,a,i6.6,a,4(a,es12.5))') 'Radar Nr. ', ista, ' (', rs_meta(ista)%station_id, &
             ') station height is equal or below model orography, ', &
             ' h_rad MSL = ', rs_meta(ista)%alt_msl, &
             ' hsurf     = ', hsurf_r, &
             ' lonrot_rad   = ', rlon_r, ' latrot_rad = ', rlat_r
        CALL abort_run (my_radar_id, 10376, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')

      END IF

      !-------------------------------------------------------------------------------
      ! .. Fill rs_meta(ista)%ista: (At the moment only needed for cdfin-output
      !     in case of the deprecated linput_old_format=.TRUE.)
      !-------------------------------------------------------------------------------

      rs_meta(ista)%ista = ista

      !-------------------------------------------------------------------------------
      ! .. Fill rs_meta(ista)%nel_present and rs_meta(ista)%ind_ele_present(:)
      !     assuming that all elevations in the nominal strategy are present.
      !    In case of obs data reading, this will be adjusted later for each obs time step.
      !-------------------------------------------------------------------------------

      rs_meta(ista)%nel_present = rs_meta(ista)%nel
      rs_meta(ista)%ind_ele_present(:) = missval_int
      rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) = (/ (i, i=1, rs_meta(ista)%nel_present) /)
      
    END DO  ! loop over nradsta


    !-------------------------------------------------------------------------------
    ! .. Check ydirradarin: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    IF (LEN_TRIM(ADJUSTL(ydirradarin)) == 0) THEN
      ydirradarin(:) = ' '
      ydirradarin    = get_model_inputdir ()
    END IF
    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ydirradarin)) > 0) THEN
      ydirradarin    = TRIM(ADJUSTL(ydirradarin)) // REPEAT(' ',LEN(ydirradarin))
      IF (ydirradarin(LEN_TRIM(ydirradarin):LEN_TRIM(ydirradarin)) /= '/') THEN
        ydirradarin  = TRIM(ydirradarin) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ydirradarout: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    IF (LEN_TRIM(ADJUSTL(ydirradarout)) == 0) THEN
      ydirradarout(:) = ' '
      ydirradarout    = get_model_outputdir ()
    END IF
    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ydirradarout)) > 0) THEN
      ydirradarout    = TRIM(ADJUSTL(ydirradarout)) // REPEAT(' ',LEN(ydirradarout))
      IF (ydirradarout(LEN_TRIM(ydirradarout):LEN_TRIM(ydirradarout)) /= '/') THEN
        ydirradarout  = TRIM(ydirradarout) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ysubdirfof: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ysubdirfof)) > 0) THEN
      ysubdirfof    = TRIM(ADJUSTL(ysubdirfof)) // REPEAT(' ',LEN(ysubdirfof))
      IF (ysubdirfof(LEN_TRIM(ysubdirfof):LEN_TRIM(ysubdirfof)) /= '/') THEN
        ysubdirfof  = TRIM(ADJUSTL(ysubdirfof)) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ysubdircomp: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ysubdircomp)) > 0) THEN
      ysubdircomp    = TRIM(ADJUSTL(ysubdircomp)) // REPEAT(' ',LEN(ysubdircomp))
      IF (ysubdircomp(LEN_TRIM(ysubdircomp):LEN_TRIM(ysubdircomp)) /= '/') THEN
        ysubdircomp  = TRIM(ADJUSTL(ysubdircomp)) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ydir_mielookup_write: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ydir_mielookup_write)) > 0) THEN
      ydir_mielookup_write    = TRIM(ADJUSTL(ydir_mielookup_write)) // REPEAT(' ',LEN(ydir_mielookup_write))
      IF (ydir_mielookup_write(LEN_TRIM(ydir_mielookup_write):LEN_TRIM(ydir_mielookup_write)) /= '/') THEN
        ydir_mielookup_write  = TRIM(ADJUSTL(ydir_mielookup_write)) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ydir_mielookup_read: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ydir_mielookup_read)) > 0) THEN
      ydir_mielookup_read    = TRIM(ADJUSTL(ydir_mielookup_read)) // REPEAT(' ',LEN(ydir_mielookup_read))
      IF (ydir_mielookup_read(LEN_TRIM(ydir_mielookup_read):LEN_TRIM(ydir_mielookup_read)) /= '/') THEN
        ydir_mielookup_read  = TRIM(ADJUSTL(ydir_mielookup_read)) // '/'
      END IF
    END IF

    !-------------------------------------------------------------------------------
    ! .. Check ydir_ready_write: only possible after last call to read_radarnamelist()
    ! Default should be ydirradarout to be backward compatible
    !-------------------------------------------------------------------------------

    IF (LEN_TRIM(ADJUSTL(ydir_ready_write)) == 0) THEN
      ydir_ready_write(:) = ' '
      ydir_ready_write    = ydirradarout
    END IF
    ! .. Add slash to the end of the string if not empty:
    IF (LEN_TRIM(ADJUSTL(ydir_ready_write)) > 0) THEN
      ydir_ready_write    = TRIM(ADJUSTL(ydir_ready_write)) // REPEAT(' ',LEN(ydir_ready_write))
      IF (ydir_ready_write(LEN_TRIM(ydir_ready_write):LEN_TRIM(ydir_ready_write)) /= '/') THEN
        ydir_ready_write  = TRIM(ydir_ready_write) // '/'
      END IF
    END IF
    
    !-------------------------------------------------------------------------------
    ! .. Check levelidlist_for_composite_glob: only possible after last call to read_radarnamelist()
    !-------------------------------------------------------------------------------

    IF (ldo_composite) THEN
      DO i=1, nel_composite
        IF (levelidlist_for_composite_glob(i) == -99) THEN
          levelidlist_for_composite_glob(i) = eleindlist_for_composite_glob(i)
        END IF
      END DO
    END IF

    !-----------------------------------------------------------------------------
    ! .. Check existence of feedback files for restart runs and prevent overwrite:
    !      only possible after last call to read_radarnamelist()
    !-----------------------------------------------------------------------------

#ifdef NUDGING
    IF ( lreadmeta_from_netcdf .AND. lfdbk_output .AND. run_is_restart() ) THEN
      DO ista = 1, nradsta
        IF (ysubdirfof(1:1) == '/') THEN
          INQUIRE (file=TRIM(ysubdirfof)//TRIM(rs_meta(ista)%fdbkfile), &
                   exist=rs_meta(ista)%lfdbkfile_exist)
        ELSE
          INQUIRE (file=TRIM(ydirradarout)//TRIM(ysubdirfof)//TRIM(rs_meta(ista)%fdbkfile), &
                   exist=rs_meta(ista)%lfdbkfile_exist)
        END IF
      END DO
    END IF
#endif


    !-------------------------------------------------------------------------------
    ! Initialize the smoothing weigths for each radar (if needed):
    !-------------------------------------------------------------------------------

    IF (lsmooth) THEN
      CALL init_smth_info()
    END IF

    !-------------------------------------------------------------------------------
    ! - Section 7: Write diagnostic output on PE 0
    !-------------------------------------------------------------------------------

    IF (my_radar_id == 0 .AND. ldebug_radsim) THEN

      WRITE(*,*)
      WRITE(*,*) "=================================================================="
      WRITE(*,*) "Information on Radar Simulator"
      WRITE(*,*) "=================================================================="
      WRITE(*,'(a,T35,"=",L8)') '  ltestpattern_hydrometeors', ltestpattern_hydrometeors
      WRITE(*,'(a,T35,"=",L8)') '  lfdbk_output', lfdbk_output
      WRITE(*,'(a,T35,"=",L8)') '  lwrite_ready', lwrite_ready
      WRITE(*,'(a,T35,"=",L8)') '  lvoldata_output', lvoldata_output
#if (defined AUXOUT_OFFLINE && defined __COSMO__)
!#ifdef AUXOUT_OFFLINE
      WRITE(*,'(a,T35,"=",L8)') '  lmodfield_output', lmodfield_output
      WRITE(*,'(a,T35,"=",L8)') '  loutnwp', loutnwp
#else
      WRITE(*,'(a,T35,"=",L8,a)') '  lmodfield_output', lmodfield_output, '    (inactive)'
      WRITE(*,'(a,T35,"=",L8,a)') '  loutnwp', loutnwp, '    (inactive)'
#endif
#ifdef AUXOUT_OFFLINE
      WRITE(*,'(a,T35,"=",L8)') '  loutvolaux', loutvolaux
#else
      WRITE(*,'(a,T35,"=",L8,a)') '  loutvolaux', loutvolaux, '    (inactive)'
#endif
      WRITE(*,'(a,T35,"=",L8)') '  lqc_flag', lqc_flag
      WRITE(*,'(a,T35,"=",L8)') '  lout_geom', lout_geom
      WRITE(*,'(a,T35,"=",L8)') '  loutradwind', loutradwind
      WRITE(*,'(a,T35,"=",L8)') '  ldealiase_vr_obs', ldealiase_vr_obs
      WRITE(*,'(a,T35,"=",L8)') '  loutdbz', loutdbz
      WRITE(*,'(a,T35,"=",L8)') '  loutpolstd', loutpolstd
      WRITE(*,'(a,T35,"=",L8)') '  loutpolall', loutpolall
      WRITE(*,'(a,T35,"=",L8)') '  lextdbz ', lextdbz
      WRITE(*,'(a,T35,"=",L8)') '  lweightdbz ', lweightdbz
      WRITE(*,'(a,T35,"=",L8)') '  lfall ', lfall
      WRITE(*,'(a,T35,"=",L8)') '  lfill_vr_backgroundwind', lfill_vr_backgroundwind
      WRITE(*,'(a,T35,"=",L8)') '  lonline ', lonline
      WRITE(*,'(a,T35,"=",L8)') '  lcomm_nonblocking_online ', lcomm_nonblocking_online
      WRITE(*,'(a,T35,"=",L8)') '  lsode ', lsode
      WRITE(*,'(a,T35,"=",L8)') '  lsmooth', lsmooth
      WRITE(*,'(a,T35,"=",L8)') '  lmds_z', lmds_z
      WRITE(*,'(a,T35,"=",L8)') '  lmds_vr', lmds_vr
      WRITE(*,'(a,T35,"=",L8)') '  lreadmeta_from_netcdf', lreadmeta_from_netcdf
      WRITE(*,'(a,T35,"=",L8)') '  lcheck_inputrecords',   lcheck_inputrecords
      WRITE(*,'(a,T35,"=",L8)') '  ldo_composite',    ldo_composite
      WRITE(*,'(a,T35,"=",L8)') '  lcomposite_output',    lcomposite_output
      WRITE(*,'(a,T35,"=",L8)') '  lcomm_nonblocking_output ', lcomm_nonblocking_output
      WRITE(*,'(a,T45,"=",i8)') '  nel_composite (global no. of elev. for composites) ',  nel_composite
      IF (nel_composite > 0) THEN
        WRITE(*,'(a,T45,"=",i8)') '  first entry in eleindlist_for_composite_glob ', eleindlist_for_composite_glob(1)
        WRITE(*,'(a,T45,"=",i8)') '  last entry in eleindlist_for_composite_glob ', eleindlist_for_composite_glob(nel_composite)
      END IF
      WRITE(*,'(a,T35,"=",L8)')      '  ldo_bubbles',      ldo_bubbles
      WRITE(*,'(a,T45,"=",i8)')      '  eleind_for_composite_bub_glob',   eleind_for_composite_bub_glob
      WRITE(*,'(a,T45,"=",es12.5)')  '  tstart_bubble_search',     tstart_bubble_search
      WRITE(*,'(a,T45,"=",es12.5)')  '  dt_bubble_search',         dt_bubble_search
      WRITE(*,'(a,T45,"=",es12.5)')  '  tend_bubble_search',       tend_bubble_search
      WRITE(*,'(a,T45,"=",es12.5)')  '  t_offset_bubble_trigger_async', t_offset_bubble_trigger_async
      WRITE(*,'(a,T45,"=",es12.5)')  '  prob_bubble',    prob_bubble
      WRITE(*,'(a,T45,"=",L8)')      '  lbub_isolated',  lbub_isolated
      WRITE(*,'(a,T45,"=",es12.5)')  '  maxdim_obs',   maxdim_obs
      WRITE(*,'(a,T45,"=",2es12.5)') '  threshold_obs',threshold_obs
      WRITE(*,'(a,T45,"=",2es12.5)') '  threshold_mod',threshold_mod
      WRITE(*,'(a,T45,"=",2es12.5)') '  areamin_mod',  areamin_mod
      WRITE(*,'(a,T45,"=",2es12.5)') '  areamin_obs',  areamin_obs
      WRITE(*,'(a,T45,"=",es12.5)')  '  mult_dist_obs',     mult_dist_obs
      WRITE(*,'(a,T45,"=",es12.5)')  '  mult_dist_mod',     mult_dist_mod
      WRITE(*,'(a,T45,"=",es12.5)')  '  add_dist_obs',    add_dist_obs
      WRITE(*,'(a,T45,"=",es12.5)')  '  add_dist_mod',     add_dist_mod
      WRITE(*,'(a,T45,"=",es12.5)') '  dt_bubble_advect',           dt_bubble_advect
      WRITE(*,'(a,T45,"=",es12.5)') '  zlow_meanwind_bubble_advect', zlow_meanwind_bubble_advect
      WRITE(*,'(a,T45,"=",es12.5)') '  zup_meanwind_bubble_advect', zup_meanwind_bubble_advect
      WRITE(*,'(a,T45,"=",a)')      '  bubble_type', TRIM(bubble_type)
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_heatingrate', bubble_heatingrate
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_timespan', bubble_timespan
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_dT', bubble_dT
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_centz', bubble_centz
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_radx', bubble_radx
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_rady', bubble_rady
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_radz', bubble_radz
      WRITE(*,'(a,T45,"=",es12.5)') '  bubble_rotangle', bubble_rotangle
      WRITE(*,'(a,T45,"=",L8)')     '  bubble_addnoise_T', bubble_addnoise_T
      WRITE(*,'(a,T45,"=",es12.5)') '    bubble_dT_noise', bubble_dT_noise
      WRITE(*,'(a,T45,"=",L8)')     '  bubble_holdrhconst', bubble_holdrhconst
      WRITE(*,'(a,T35,"=",i8)')     '  itype_supobing', itype_supobing
      WRITE(*,'(a,T35,"=",i8)')     '  supob_nrb', supob_nrb
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_azi_maxsector_vr', supob_azi_maxsector_vr
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_cart_resolution', supob_cart_resolution
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_ave_width', supob_ave_width
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_minrange_vr', supob_minrange_vr
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_minrange_z', supob_minrange_z
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_vrw', supob_vrw
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_rfl', supob_rfl
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_lowthresh_z_obs', supob_lowthresh_z_obs
      WRITE(*,'(a,T35,"=",es12.4)') '  supob_lowthresh_z_sim', supob_lowthresh_z_sim
      WRITE(*,'(a,T35,"=",i8)')     '  itype_obserr_vr', itype_obserr_vr
      WRITE(*,'(a,T35,"=",es12.4)') '    baseval_obserr_vr', baseval_obserr_vr
      WRITE(*,'(a,T35,"=",es12.4)') '    maxval_obserr_vr', maxval_obserr_vr
      WRITE(*,'(a,T35,"=",es12.4)') '    ramp_lowdbz_obserr_vr', ramp_lowdbz_obserr_vr
      WRITE(*,'(a,T35,"=",es12.4)') '    ramp_highdbz_obserr_vr', ramp_highdbz_obserr_vr
      WRITE(*,'(a,T35,"=",i8)')     '  itype_metric_refl_fdbk', itype_metric_refl_fdbk
      WRITE(*,'(a,T35,"=",es12.4)') '    minval_obserr_lwc', minval_obserr_lwc
      WRITE(*,'(a,T35,"=",i8)')     '  thin_step_azi', thin_step_azi
      WRITE(*,'(a,T35,"=",i8)')     '  thin_step_range', thin_step_range
      WRITE(*,'(a,T35,"=",i8)')     '  thin_step_ele', thin_step_ele
      IF (nel_fdbk_glob > 0) THEN
        WRITE(*,'(a,T45,"=",i8)')   '  global no. of elev. for feedback', nel_fdbk_glob
        WRITE(*,'(a,T45,"=",i8)')   '  first entry in ind_ele_fdbk_glob', ind_ele_fdbk_glob(1)
        WRITE(*,'(a,T45,"=",i8)')   '  last entry in ind_ele_fdbk_glob', ind_ele_fdbk_glob(nel_fdbk_glob)
      ELSE
        WRITE(*,'(a)') '  global elevation choice for feedback not active'
        WRITE(*,'(a)') '  ind_ele_fdbk_glob not active'
        WRITE(*,'(a)') '  ind_ele_fdbk_glob not active'
      END IF
      IF (nel_voldata_glob > 0) THEN
        WRITE(*,'(a,T45,"=",i8)') '  global no. of elev. for volume data', nel_voldata_glob
        WRITE(*,'(a,T45,"=",i8)') '  first entry in ind_ele_voldata_glob', ind_ele_voldata_glob(1)
        WRITE(*,'(a,T45,"=",i8)') '  last entry in ind_ele_voldata_glob', ind_ele_voldata_glob(nel_voldata_glob)
      ELSE
        WRITE(*,'(a)') '  global elevation choice for volume data not active'
        WRITE(*,'(a)') '  ind_ele_voldata_glob not active'
        WRITE(*,'(a)') '  ind_ele_voldata_glob not active'
      END IF
      IF (nobs_times_fdbk_glob > 0) THEN
        WRITE(*,'(a,T45,"=",i8)') '  global no. of obs times for feedback', nobs_times_fdbk_glob
        WRITE(*,'(a,T45,"=",f10.1)') '  first entry in obs_times_fdbk_glob', obs_times_fdbk_glob(1)
        WRITE(*,'(a,T45,"=",f10.1)') '  last entry in obs_times_fdbk_glob', obs_times_fdbk_glob(nobs_times_fdbk_glob)
      ELSE
        WRITE(*,'(a)') '  global obs times choice for feedback not active'
        WRITE(*,'(a)') '  obs_times_fdbk_glob not active'
        WRITE(*,'(a)') '  obs_times_fdbk_glob not active'
      END IF
      IF (nobs_times_voldata_glob > 0) THEN
        WRITE(*,'(a,T45,"=",i8)') '  global no. of obs times for volume data', nobs_times_voldata_glob
        WRITE(*,'(a,T45,"=",f10.1)') '  first entry in obs_times_voldata_glob', obs_times_voldata_glob(1)
        WRITE(*,'(a,T45,"=",f10.1)') '  last entry in obs_times_voldata_glob', obs_times_voldata_glob(nobs_times_voldata_glob)
      ELSE
        WRITE(*,'(a)') '  global obs times choice for volume data not active'
        WRITE(*,'(a)') '  obs_times_voldata_glob not active'
        WRITE(*,'(a)') '  obs_times_voldata_glob not active'
      END IF
      WRITE(*,'(a,T35,"=",i8)') '  nradsta', nradsta
      WRITE(*,*)
      WRITE(*,*) "=================================================================="
      WRITE(*,*) "Information on Radar Stations"
      WRITE(*,*) "=================================================================="
      DO ista = 1, nradsta
        WRITE(*,'(a,i6.6,a,i0)') "Radar station: "//TRIM(ADJUSTL(rs_meta(ista)%station_name))//"   ID: ", &
             rs_meta(ista)%station_id,"   ista = ",ista
        WRITE(*,'(" lat: ",F7.3,2X,"lon: ",F7.3,2X,"alt_agl_mod: ",F9.3,'//&
             '2X,"alt_msl: ",F9.3,2X,"msl_mod: ",F9.3)') &
             rs_meta(ista)%lat,rs_meta(ista)%lon,rs_meta(ista)%alt_agl_mod, &
             rs_meta(ista)%alt_msl,rs_meta(ista)%msl_mod
        WRITE(*,'(" range bins (number/increment): ",I4,"/",F8.2,/,'//&
             '" azimuths (number/increment/start): ",I4,"/",F6.2,"/",F6.2)') &
             rs_meta(ista)%nra, rs_meta(ista)%ra_inc, &
             rs_meta(ista)%naz, rs_meta(ista)%az_inc, rs_meta(ista)%az_start
        IF (lsmooth) THEN
          WRITE(*,'("  Phi3: ",f6.2,/,"  Theta3: ",f6.2,/,"  dalpha: ",f6.2,' &
               //'/,"  alpha3_eff_0: ",f7.3)') &
               rs_meta(ista)%Phi3, rs_meta(ista)%Theta3, rs_meta(ista)%dalpha, &
               rs_meta(ista)%alpha3_eff_0
        ENDIF
        IF (lmds_z .OR. lmds_vr) THEN
          WRITE (*,'("  Minimum detectable signal: ",f8.2, " dBZ in ",f8.2," km")') &
               rs_meta(ista)%mds_Z0, rs_meta(ista)%mds_r0*0.001_dp
        END IF
        WRITE(*,*) "elevations: "
        WRITE(*,*) (rs_meta(ista)%el_arr(i),i=1,rs_meta(ista)%nel)
        WRITE(*,*) "observation times:   no. obs. times = ", rs_meta(ista)%nobs_times
        cnobs(:) = ' '
        IF (rs_meta(ista)%nobs_times > 10) THEN
          WRITE(cnobs,'(i10)') 10
          WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |")," ...")') &
               (rs_meta(ista)%obs_times(i),i=1,10)
        ELSE IF (rs_meta(ista)%nobs_times > 0) THEN
          WRITE(cnobs,'(i10)') rs_meta(ista)%nobs_times
          WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |"))') &
               (rs_meta(ista)%obs_times(i),i=1,rs_meta(ista)%nobs_times)
        END IF
        IF (lvoldata_output) THEN
          WRITE(*,*) "observation times in volume data:   no. obs. times = ", rs_meta(ista)%nobs_times_voldata
          cnobs(:) = ' '
          IF (rs_meta(ista)%nobs_times_voldata > 10) THEN
            WRITE(cnobs,'(i10)') 10
            WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |")," ...")') &
                 (rs_meta(ista)%obs_times_voldata(i),i=1,10)
          ELSE IF (rs_meta(ista)%nobs_times_voldata > 0) THEN
            WRITE(cnobs,'(i10)') rs_meta(ista)%nobs_times_voldata
            WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |"))') &
                 (rs_meta(ista)%obs_times_voldata(i),i=1,rs_meta(ista)%nobs_times_voldata)
          END IF
        END IF
        IF (lreadmeta_from_netcdf) THEN
          WRITE(*,*) "observation times in feedback:   no. obs. times = ", rs_meta(ista)%nobs_times_fdbk
          cnobs(:) = ' '
          IF (rs_meta(ista)%nobs_times_fdbk > 10) THEN
            WRITE(cnobs,'(i10)') 10
            WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |")," ...")') &
                 (rs_meta(ista)%obs_times_fdbk(i),i=1,10)
          ELSE IF (rs_meta(ista)%nobs_times_fdbk > 0) THEN
            WRITE(cnobs,'(i10)') rs_meta(ista)%nobs_times_fdbk
            WRITE(*,'('//TRIM(ADJUSTL(cnobs))//'(f8.1," |"))') &
                 (rs_meta(ista)%obs_times_fdbk(i),i=1,rs_meta(ista)%nobs_times_fdbk)
          END IF
        END IF
        WRITE(*,*) "--------------------------------------------------"
      END DO

    END IF

    !--------------------------------------------------------------------------------------
    ! - Section 8: Write namelist parameters and their defaults to file unit ctrlfile_unit:
    !--------------------------------------------------------------------------------------

    IF ( my_radar_id == 0 ) THEN

      ctrlfile(:) = ' '
      IF (LEN_TRIM(nmlctrlfilepath) > 0) THEN
        ctrlfile = TRIM(nmlctrlfilepath) // '/' // TRIM(nmlctrlfile)
      ELSE
        ctrlfile = TRIM(nmlctrlfile)
      ENDIF

      CALL get_free_funit(ctrlfile_unit)
      IF (lfirst_ctrl_output) THEN
        OPEN (ctrlfile_unit, file=TRIM(ctrlfile), form='formatted', status='replace', &
             iostat=ierr)
        lfirst_ctrl_output = .FALSE.
      ELSE
        OPEN (ctrlfile_unit, file=TRIM(ctrlfile), form='formatted', status='old', &
             position='append', iostat=ierr)
      END IF
      IF (ierr /= 0) THEN
        errstring(:) = ' '
        WRITE (errstring,'(a,i0)') 'Error during opening file '//TRIM(ctrlfile)//' for domain ',izdom
        CALL abort_run (my_radar_id, 10085, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')
      END IF

      WRITE (ctrlfile_unit, '(A2)')  '  '
      WRITE (ctrlfile_unit, '(A,I3)') '0    NAMELIST: radarsim_params for domain dom = ', izdom
      WRITE (ctrlfile_unit, '(A)')    '     ----------------------------------------------'
      WRITE (ctrlfile_unit, '(A2)')  '  '
      WRITE (ctrlfile_unit, '(A)')  'Variables for radar simulator'
      WRITE (ctrlfile_unit, TRIM(format_char)) 'Variable', 'Actual Value',   &
           'Default Value', 'Format'
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'dom', dom, missval_int, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'nradsta_namelist', nradsta_namelist, nradsta_namelist_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'nradsta', nradsta, nradsta_namelist_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'icountry', icountry, icountry_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ldebug_radsim', ldebug_radsim, ldebug_radsim_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ltestpattern_hydrometeors', ltestpattern_hydrometeors, ltestpattern_hydrometeors_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lvoldata_output', lvoldata_output, lvoldata_output_d     , ' L '
      DO ii=1, noutstreams_max
        WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T40,A,T55,A,T70,A)')                             &
             'voldata_ostream(',ii,')%format', TRIM(voldata_ostream(ii)%format), TRIM(voldata_ostream_d(ii)%format), 'C*12'
        WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T40,A,T55,A,T70,A)')                             &
             'voldata_ostream(',ii,')%grib2_packingtype', TRIM(voldata_ostream(ii)%grib2_packingtype), &
                                                          TRIM(voldata_ostream_d(ii)%grib2_packingtype), 'C*12'
        IF ( ALL( LEN_TRIM(voldata_ostream(ii)%output_list(:)) == 0 ) ) THEN
          WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T40,A,T55,A,T70,A)')                             &
               'voldata_ostream(',ii,')%output_list(01)', '-', '-', 'C*12'
        ELSE
          DO  i=1, noutput_fields_max
            IF ( LEN_TRIM(voldata_ostream(ii)%output_list(i)) > 0) THEN
              WRITE (ctrlfile_unit, '(T2,A,i2.2,A,i2.2,A,T40,A,T55,A,T70,A)')                             &
                   'voldata_ostream(',ii,')%output_list(',i,')', TRIM(voldata_ostream(ii)%output_list(i)), '-', 'C*12'
            END IF
          END DO
          IF (LEN_TRIM(voldata_ostream(ii)%file_pattern) > 0) THEN
            WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T40,A,/,T40,A,T70,A)')        &
                 'voldata_ostream(',ii,')%file_pattern', TRIM(voldata_ostream(ii)%file_pattern), '-', 'C*300'
          ELSE
            WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T40,A,/,T40,A,T70,A)')       &
                 'voldata_ostream(',ii,')%file_pattern', '-', '-', 'C*300'
          END IF
          WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T50,A,T65,A,T80,A)')                             &
               'voldata_ostream(',ii,')%pat_scantype_volscan', voldata_ostream(ii)%pat_scantype_volscan, &
                                                               voldata_ostream_d(ii)%pat_scantype_volscan, 'C*12'
          WRITE (ctrlfile_unit, '(T2,A,i2.2,A,T50,A,T65,A,T80,A)')                             &
               'voldata_ostream(',ii,')%pat_scantype_precipscan', voldata_ostream(ii)%pat_scantype_precipscan, &
                                                                  voldata_ostream_d(ii)%pat_scantype_precipscan, 'C*12'
          WRITE (ctrlfile_unit, TRIM(format_real_f_vec2))                       &
               'voldata_ostream(',ii,')%content_tref', voldata_ostream(ii)%content_tref, voldata_ostream_d(ii)%content_tref, ' R '
          WRITE (ctrlfile_unit, TRIM(format_real_f_vec2))                       &
               'voldata_ostream(',ii,')%content_dt', voldata_ostream(ii)%content_dt, voldata_ostream_d(ii)%content_dt, ' R '
        END IF
      END DO
#if (defined AUXOUT_OFFLINE && defined __COSMO__)
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lmodfield_output', lmodfield_output, lmodfield_output_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutnwp', loutnwp, loutnwp_d, ' L '
#endif
#ifdef AUXOUT_OFFLINE
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutvolaux', loutvolaux, loutvolaux_d, ' L '
#endif
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lqc_flag', lqc_flag, lqc_flag_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ldealiase_vr_obs', ldealiase_vr_obs, ldealiase_vr_obs_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lfdbk_output', lfdbk_output, lfdbk_output_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lwrite_ready', lwrite_ready, lwrite_ready_d     , ' L '
      IF (LEN_TRIM(ready_file_pattern) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                        &
             'ready_file_pattern', TRIM(ready_file_pattern), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                        &
             'ready_file_pattern', '-', '-', 'C*300'
      END IF
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lout_geom', lout_geom, lout_geom_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutradwind', loutradwind, loutradwind_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutdbz', loutdbz, loutdbz_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcalc_dbz_on_radarbins', lcalc_dbz_on_radarbins, lcalc_dbz_on_radarbins_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutpolstd', loutpolstd, loutpolstd_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'loutpolall', loutpolall, loutpolall_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lextdbz', lextdbz, lextdbz_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'itype_mpipar_lookupgen', itype_mpipar_lookupgen, itype_mpipar_lookupgen_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'pe_start_lookupgen', pe_start_lookupgen, pe_start_lookupgen_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'pe_end_lookupgen', pe_end_lookupgen, pe_end_lookupgen_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'llookup_interp_mode_dualpol', llookup_interp_mode_dualpol, llookup_interp_mode_dualpol_d     , ' L '
      CALL ctrl_output_dbzmeta_nuspec_fwo(ctrlfile_unit, -1, dbz_meta_glob, dbz_meta_glob_d, 'dbz_meta_glob')
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lweightdbz', lweightdbz, lweightdbz_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lfall', lfall, lfall_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lfill_vr_backgroundwind', lfill_vr_backgroundwind, lfill_vr_backgroundwind_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lonline', lonline, lonline_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcomm_nonblocking_online', lcomm_nonblocking_online, lcomm_nonblocking_online_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lsode', lsode, lsode_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lsmooth', lsmooth, lsmooth_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'ngpsm_h_glob', ngpsm_h_glob, ngpsm_h_glob_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'ngpsm_v_glob', ngpsm_v_glob, ngpsm_v_glob_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lmds_z', lmds_z, lmds_z_d        , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lmds_vr', lmds_vr, lmds_vr_d     , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcomm_nonblocking_output', lcomm_nonblocking_output, lcomm_nonblocking_output_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lreadmeta_from_netcdf', lreadmeta_from_netcdf, lreadmeta_from_netcdf_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcheck_inputrecords', lcheck_inputrecords, lcheck_inputrecords_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lequal_azi_alldatasets', lequal_azi_alldatasets, lequal_azi_alldatasets_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'labort_if_problems_obsfiles', labort_if_problems_obsfiles, labort_if_problems_obsfiles_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'labort_if_problems_gribout', labort_if_problems_gribout, labort_if_problems_gribout_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ldo_composite', ldo_composite, ldo_composite_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcomposite_output', lcomposite_output, lcomposite_output_d , ' L '
      IF (LEN_TRIM(composite_file_pattern) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'composite_file_pattern', TRIM(composite_file_pattern), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'composite_file_pattern', '-', '-', 'C*300'
      END IF
      WRITE (ctrlfile_unit, TRIM(format_char))                             &
           'comp_grib2_packingtype', TRIM(comp_grib2_packingtype), TRIM(comp_grib2_packingtype_d), 'C*12'
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,I12,T55,I12,T70,A)')                  &
           'comp_meta', 'ni', comp_meta%ni, comp_meta_d%ni, ' I '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,I12,T55,I12,T70,A)')                  &
           'comp_meta', 'nj', comp_meta%nj, comp_meta_d%nj, ' I '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'pollon', comp_meta%pollon, comp_meta_d%pollon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'pollat', comp_meta%pollat, comp_meta_d%pollat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'polgam', comp_meta%polgam, comp_meta_d%polgam, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'dlon', comp_meta%dlon, comp_meta_d%dlon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'dlat', comp_meta%dlat, comp_meta_d%dlat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'startlon', comp_meta%startlon, comp_meta_d%startlon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta', 'startlat', comp_meta%startlat, comp_meta_d%startlat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.1,T55,f12.1,T70,A)')              &
           'comp_meta', 'r_earth', comp_meta%r_earth, comp_meta_d%r_earth, ' R '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'nel_composite', nel_composite, nel_composite_d, ' I '
      DO i = 1, MAX(nel_composite, 1)
        WRITE (ctrlfile_unit, '(T4,A,i3.3,A,T40,I12,T55,I12,T70,A)')     &
             'eleindlist_for_composite_glob(',i,')', eleindlist_for_composite_glob(i), eleindlist_for_composite_glob_d(i), ' I '
      END DO
      DO i = 1, MAX(nel_composite, 1)
        WRITE (ctrlfile_unit, '(T4,A,i3.3,A,T40,I12,T55,I12,T70,A)')     &
             'levelidlist_for_composite_glob(',i,')', levelidlist_for_composite_glob(i), levelidlist_for_composite_glob_d(i), ' I '
      END DO
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ldo_bubbles', ldo_bubbles, ldo_bubbles_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'ldo_bubbles_manual', ldo_bubbles_manual, ldo_bubbles_manual_d , ' L '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lcomposite_output_bub', lcomposite_output_bub, lcomposite_output_bub_d , ' L '
      IF (LEN_TRIM(composite_file_pattern_bub) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'composite_file_pattern_bub', TRIM(composite_file_pattern_bub), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'composite_file_pattern_bub', '-', '-', 'C*300'
      END IF
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,I12,T55,I12,T70,A)')                  &
           'comp_meta_bub', 'ni', comp_meta_bub%ni, comp_meta_bub_d%ni, ' I '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,I12,T55,I12,T70,A)')                  &
           'comp_meta_bub', 'nj', comp_meta_bub%nj, comp_meta_bub_d%nj, ' I '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'pollon', comp_meta_bub%pollon, comp_meta_bub_d%pollon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'pollat', comp_meta_bub%pollat, comp_meta_bub_d%pollat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'polgam', comp_meta_bub%polgam, comp_meta_bub_d%polgam, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'dlon', comp_meta_bub%dlon, comp_meta_bub_d%dlon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'dlat', comp_meta_bub%dlat, comp_meta_bub_d%dlat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'startlon', comp_meta_bub%startlon, comp_meta_bub_d%startlon, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.4,T55,f12.4,T70,A)')              &
           'comp_meta_bub', 'startlat', comp_meta_bub%startlat, comp_meta_bub_d%startlat, ' R '
      WRITE (ctrlfile_unit, '(T8,A,"%",A,T40,f12.1,T55,f12.1,T70,A)')              &
           'comp_meta_bub', 'r_earth', comp_meta_bub%r_earth, comp_meta_bub_d%r_earth, ' R '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'eleind_for_composite_bub_glob', eleind_for_composite_bub_glob, eleind_for_composite_bub_glob_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'tstart_bubble_search', tstart_bubble_search, tstart_bubble_search_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'dt_bubble_search', dt_bubble_search, dt_bubble_search_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'tend_bubble_search', tend_bubble_search, tend_bubble_search_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           't_offset_bubble_trigger_async', t_offset_bubble_trigger_async, t_offset_bubble_trigger_async_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'prob_bubble', prob_bubble, prob_bubble_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'lbub_isolated', lbub_isolated, lbub_isolated_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'maxdim_obs', maxdim_obs, maxdim_obs_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'threshold_obs(1)', threshold_obs(1), threshold_obs_d(1), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'threshold_obs(2)', threshold_obs(2), threshold_obs_d(2), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'threshold_mod(1)', threshold_mod(1), threshold_mod_d(1), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'threshold_mod(2)',threshold_mod(2) , threshold_mod_d(2), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'areamin_mod(1)', areamin_mod(1), areamin_mod_d(1), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'areamin_mod(2)', areamin_mod(2), areamin_mod_d(2), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'areamin_obs(1)', areamin_obs(1), areamin_obs_d(1), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'areamin_obs(2)', areamin_obs(2), areamin_obs_d(2), ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'mult_dist_obs', mult_dist_obs, mult_dist_obs_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'mult_dist_mod', mult_dist_mod, mult_dist_mod_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'add_dist_obs', add_dist_obs, add_dist_obs_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'add_dist_mod', add_dist_mod, add_dist_mod_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'dt_bubble_advect', dt_bubble_advect, dt_bubble_advect_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'zlow_meanwind_bubble_advect', zlow_meanwind_bubble_advect, zlow_meanwind_bubble_advect_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'zup_meanwind_bubble_advect', zup_meanwind_bubble_advect, zup_meanwind_bubble_advect_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_char))                                 &
           'bubble_type', TRIM(bubble_type), TRIM(bubble_type_d), 'C*12'
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_heatingrate', bubble_heatingrate, bubble_heatingrate_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_timespan', bubble_timespan, bubble_timespan_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_dT', bubble_dT, bubble_dT_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_centz', bubble_centz, bubble_centz_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_radx', bubble_radx, bubble_radx_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_rady', bubble_rady, bubble_rady_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_radz', bubble_radz, bubble_radz_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_rotangle', bubble_rotangle, bubble_rotangle_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'bubble_addnoise_T', bubble_addnoise_T, bubble_addnoise_T_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'bubble_dT_noise', bubble_dT_noise, bubble_dT_noise_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_logical))                       &
           'bubble_holdrhconst', bubble_holdrhconst, bubble_holdrhconst_d, ' L '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'itype_supobing', itype_supobing, itype_supobing_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'supob_nrb', supob_nrb, supob_nrb_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_azi_maxsector_vr', supob_azi_maxsector_vr, supob_azi_maxsector_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_cart_resolution', supob_cart_resolution, supob_cart_resolution_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_ave_width', supob_ave_width, supob_ave_width_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_minrange_vr', supob_minrange_vr, supob_minrange_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_minrange_z', supob_minrange_z, supob_minrange_z_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_vrw', supob_vrw, supob_vrw_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_rfl', supob_rfl, supob_rfl_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_lowthresh_z_obs', supob_lowthresh_z_obs, supob_lowthresh_z_obs_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'supob_lowthresh_z_sim', supob_lowthresh_z_sim, supob_lowthresh_z_sim_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'itype_metric_refl_fdbk', itype_metric_refl_fdbk, itype_metric_refl_fdbk_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'minval_obserr_lwc', minval_obserr_lwc, minval_obserr_lwc_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'itype_obserr_vr', itype_obserr_vr, itype_obserr_vr_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'baseval_obserr_vr', baseval_obserr_vr, baseval_obserr_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'maxval_obserr_vr', maxval_obserr_vr, maxval_obserr_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'ramp_lowdbz_obserr_vr', ramp_lowdbz_obserr_vr, ramp_lowdbz_obserr_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_real))                       &
           'ramp_highdbz_obserr_vr', ramp_highdbz_obserr_vr, ramp_highdbz_obserr_vr_d, ' R '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'thin_step_azi', thin_step_azi, thin_step_azi_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'thin_step_range', thin_step_range, thin_step_range_d, ' I '
      WRITE (ctrlfile_unit, TRIM(format_int))                       &
           'thin_step_ele', thin_step_ele, thin_step_ele_d, ' I '
      DO i = 1, MAX(nel_fdbk_glob, 1)
        WRITE (ctrlfile_unit, TRIM(format_int_vec))     &
             'ind_ele_fdbk_glob(',i,')', ind_ele_fdbk_glob(i), ind_ele_fdbk_glob_d(i), ' I '
      END DO
      DO i = 1, MAX(nel_voldata_glob, 1)
        WRITE (ctrlfile_unit, TRIM(format_int_vec))     &
             'ind_ele_voldata_glob(',i,')', ind_ele_voldata_glob(i), ind_ele_voldata_glob_d(i), ' I '
      END DO
      DO i = 1, 3
        WRITE (ctrlfile_unit, TRIM(format_real_f_vec))                       &
             'dt_obs_fdbk_glob(',i,')', dt_obs_fdbk_glob(i), dt_obs_fdbk_glob_d(i), ' R '
      END DO
      DO i = 1, MAX(nobs_times_fdbk_glob, 1)
        WRITE (ctrlfile_unit, TRIM(format_real_f_vec))     &
             'obs_times_fdbk_glob(',i,')', obs_times_fdbk_glob(i), obs_times_fdbk_glob_d(i), ' R '
      END DO
      DO i = 1, 3
        WRITE (ctrlfile_unit, TRIM(format_real_f_vec))                       &
             'dt_obs_voldata_glob(',i,')', dt_obs_voldata_glob(i), dt_obs_voldata_glob_d(i), ' R '
      END DO
      WRITE (ctrlfile_unit, TRIM(format_real_f))                       &
           'ra_inc_coarse_glob', ra_inc_coarse_glob, ra_inc_coarse_glob_d, ' R '
      DO i = 1, MAX(nobs_times_voldata_glob, 1)
        WRITE (ctrlfile_unit, TRIM(format_real_f_vec))     &
             'obs_times_voldata_glob(',i,')', obs_times_voldata_glob(i), obs_times_voldata_glob_d(i), ' R '
      END DO
      IF (LEN_TRIM(ydirradarin) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirradarin', TRIM(ydirradarin), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirradarin', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ydirlistfile) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirlistfile', TRIM(ydirlistfile), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirlistfile', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ydirradarout) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirradarout', TRIM(ydirradarout), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydirradarout', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ysubdirfof) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ysubdirfof', TRIM(ysubdirfof), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ysubdirfof', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ysubdircomp) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ysubdircomp', TRIM(ysubdircomp), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ysubdircomp', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ydir_mielookup_read) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_mielookup_read', TRIM(ydir_mielookup_read), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_mielookup_read', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ydir_mielookup_write) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_mielookup_write', TRIM(ydir_mielookup_write), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_mielookup_write', '-', '-', 'C*300'
      END IF
      IF (LEN_TRIM(ydir_ready_write) > 0) THEN
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_ready_write', TRIM(ydir_ready_write), '-', 'C*300'
      ELSE
        WRITE (ctrlfile_unit, TRIM(format_path))                             &
             'ydir_ready_write', '-', '-', 'C*300'
      END IF
      DO ista = 1, nradsta
        WRITE (ctrlfile_unit, '(T2,a,i3.3,a)') 'Meta data of radar station ',ista, ' :'
        CALL ctrl_output_ista_nuspec_fwo(ctrlfile_unit, ista, rs_meta(ista), rs_meta_for_country(:), &
             dbz_meta(ista), dbz_meta_d(1))
      END DO

      CLOSE (ctrlfile_unit)

    END IF

    ! Clean up memory:
    IF (ALLOCATED(rs_meta_ncdf)) DEALLOCATE (rs_meta_ncdf)
    IF (ALLOCATED(rs_meta_for_country)) DEALLOCATE (rs_meta_for_country)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    SUBROUTINE read_radarnamelist (iidom)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: iidom

      INTEGER :: iz_err, dom_d, iii

      INQUIRE ( file=TRIM(ADJUSTL(nlradarsimfile)) , exist=nlfile_exists )

      IF ( .NOT.nlfile_exists ) THEN

        WRITE(*,'(a)') 'THE NAMELIST FILE '//TRIM(ADJUSTL(nlradarsimfile))// &
                       ' FOR THE RADAR SIMULATOR DOES NOT EXIST!'
        iz_err = 23456

      ELSE   ! nlfile_exists

        iz_err  = 0
        io_errmsg(:) = ' '
        CALL get_free_funit(nuin)
        OPEN(nuin , FILE=TRIM(ADJUSTL(nlradarsimfile)), FORM='FORMATTED', STATUS='OLD',  &
             IOSTAT=iz_err, IOMSG=io_errmsg)
        IF(iz_err /= 0) THEN
          WRITE(*,'(a)') ' ERROR    *** Error while opening file '// &
               TRIM(ADJUSTL(nlradarsimfile))//'*** '
          WRITE(*,'(a)') ' ERROR from OPEN(FILE='//TRIM(ADJUSTL(nlradarsimfile))//'): '//TRIM(io_errmsg)
          iz_err = -3
        ELSE
          ! Endless loop to find the namelist group for the specific domain iidom:
          !  (There should be one group for each domain in the namelist input file)
          dom_d = missval_int
          DO
            dom = dom_d
            iz_err = 0
            io_errmsg(:) = ' '
            READ (nuin, NML=radarsim_params, IOSTAT=iz_err, IOMSG=io_errmsg)
            IF (iz_err > 0) THEN
              WRITE(*,'(a)') ' ERROR    *** '// &
                             'Error while reading a NAMELIST group /RADARSIM_PARAMS/ in file '// &
                             TRIM(nlradarsimfile)// &
                             ', probably wrong values ***'
              WRITE(*,'(a)') ' ERROR from READ(NML=radarsim_params): '//TRIM(io_errmsg)
              iz_err = -2
            ELSE IF (iz_err < 0) THEN
              WRITE(*,'(a,i0,a)') ' ERROR    *** Premature EOF encountered while trying to read '// &
                                  'NAMELIST /RADARSIM_PARAMS/ for domain ',iidom,' ***'
              WRITE(*,'(a,i0,a)') '          *** '// &
                                  'Probably missing the namelist group for the specific domain dom = ', &
                                  iidom,' ***'
              WRITE(*,'(a)') ' ERROR from READ(NML=radarsim_params): '//TRIM(io_errmsg)
              iz_err = -1
            ENDIF
            IF (iz_err /= 0) EXIT
            IF (dom == dom_d) THEN
              WRITE(*,'(a)') ' ERROR    *** A NAMELIST group /RADARSIM_PARAMS/ in file '//TRIM(nlradarsimfile)// &
                   ' is missing the dom parameter (domain) ***'
               iz_err = -4
            END IF
            ! If correct namelist group for domain iidom has been read, exit the loop:
            IF (dom == iidom) EXIT
          END DO
        ENDIF

        CLOSE(nuin)

      END IF

      IF  (iz_err /= 0) THEN

        errstring(:) = ' '
        WRITE (errstring,'(a)') 'Error during opening/reading file '//TRIM(nlradarsimfile)
        CALL abort_run (my_radar_id, 10077, &
             'ERROR: problem in input_radarnamelist(): '//TRIM(ADJUSTL(errstring)), &
             'radar_namelist_read.f90, input_radarnamelist()')

      END IF

    END SUBROUTINE read_radarnamelist

    SUBROUTINE check_and_fix_dtobs ( dt_obs, errstring, ierr )

      REAL(kind=dp), INTENT(inout)  :: dt_obs(3)
      CHARACTER(len=*), INTENT(out) :: errstring
      INTEGER, INTENT(out)          :: ierr

      ierr = 0
      errstring(:) = ' '

      IF (dt_obs(1) > 0.0_dp .AND. dt_obs(2) < miss_threshold .AND. dt_obs(3) < miss_threshold) THEN
        ! .. "Traditional" notation, only the dt is given in the first element of dt_obs and the
        !    output timesteps span the entire model run and are synced with the time 0.0. Resort this to the new advanced notation:
        dt_obs(3)  = dt_obs(1)
        dt_obs(1)  = get_domain_starttime_in_sec(izdom) - eps_dtobs  ! eps is signal for sync with t = 0.0 s
        dt_obs(2)  = get_domain_endtime_in_sec(izdom)
      ELSE IF (ANY(dt_obs > miss_threshold)) THEN
        ! .. Advanced notation with a triplet: dt_obs(1:3) "t_start":"t_end":"t_incr"
        IF (dt_obs(3) > 0.0_dp) THEN
          IF (dt_obs(1) < miss_threshold) THEN
            dt_obs(1) = get_domain_starttime_in_sec(izdom) - eps_dtobs  ! eps is signal for sync with t = 0.0 s
          ELSE
            dt_obs(1) = MAX( dt_obs(1), get_domain_starttime_in_sec(izdom) )
          END IF
          IF (dt_obs(2) < miss_threshold) THEN
            dt_obs(2) = get_domain_endtime_in_sec(izdom)
          ELSE
            dt_obs(2) = MIN( dt_obs(2), get_domain_endtime_in_sec(izdom) )
          END IF
        ELSE IF (dt_obs(3) > miss_threshold) THEN
          errstring(:) = ' '
          WRITE (errstring,'(a,3(1x,f0.1,","),1x,a)') 'wrong triplet:', dt_obs, &
               'should either be (incr > 0.0, miss ,miss) or ( starttime/miss, endtime/miss, incr > 0.0)'
          ierr = 1
        END IF
      ELSE
        ! .. dt_obs contains only unused_values and is left untouched
        CONTINUE
      END IF
      
    END SUBROUTINE check_and_fix_dtobs
    
    SUBROUTINE create_obstimes_from_dtobs (dt_obs, nobs_times, obs_times, errstring, ierr)

      ! .. dt_obs has to follow the advanced triplet notation (tstart, tend, tinc) and all values of the triplet must be > miss_threshold.
      !    And if dt_obs(1) has not been explicitly given in the namelist (then it has been altered by eps),
      !    the times are synced to t = 0.0 s. Otherwise, they span exactly from dt_obs(1) to dt_obs(2) in steps of dt_obs(3).
      
      REAL(kind=dp), INTENT(in)     :: dt_obs(3)
      INTEGER, INTENT(inout)        :: nobs_times
      REAL(kind=dp), INTENT(inout)  :: obs_times(:)
      CHARACTER(len=*), INTENT(out) :: errstring
      INTEGER, INTENT(out)          :: ierr

      INTEGER       :: i, ii
      REAL(kind=dp) :: tmptime, ztstart, zt_init, zt_end
      
      ierr = 0
      errstring(:) = ' '

      IF (nobs_times < 0) THEN
        nobs_times = num_regular_obstimes ( izdom, dt_obs )
      END IF

      zt_init = get_domain_starttime_in_sec(izdom)
      zt_end  = get_domain_endtime_in_sec(izdom)
      
      IF ( NINT((zt_init-dt_obs(1))/eps_dtobs) == 1 ) THEN
        ! dt_obs(1) not explicitly specified, compute the series of obs_times in such a way that 0.0 s is part of the obs_times:
        ztstart = INT(dt_obs(1)/dt_obs(3)) * dt_obs(3)
      ELSE
        ztstart = dt_obs(1)
      END IF
      ii = 0
      DO i = 1, nobs_times
        tmptime = ztstart + (i-1) * dt_obs(3)
        IF ( check_obstime_within_modelrun(izdom, tmptime) ) THEN
          ii = ii + 1
          IF (ii > SIZE(obs_times)) THEN
            errstring(:) = ' '
            WRITE (errstring,'(a,f0.1,a)') '(3) too small (', dt_obs(3) , ' s) '// &
                 'so that nobs_times > nobstimes_max'
            ierr = 1
            EXIT
          END IF
          obs_times(ii)  = tmptime
        END IF
      END DO
      nobs_times = ii
      
    END SUBROUTINE create_obstimes_from_dtobs

  END SUBROUTINE input_radarnamelist

  !==============================================================================
  !+ Module procedures in radar_src for control output of the radar namelist
  !  parameters in the structures rs_meta(ista) and dbz_meta(ista) to a file
  !  (in COSMO per default to  "YUSPECIF_RADAR")
  !------------------------------------------------------------------------------

  SUBROUTINE ctrl_output_ista_nuspec_fwo(funit, ista, rs_m, rs_m_for_country, dbz_m, dbz_m_d)

    !------------------------------------------------------------------------------
    !
    ! Description: Control output of the radar namelist
    !              parameters in the structures rs_meta(ista) and dbz_meta(ista) to a file
    !              (in COSMO per default to  "YUSPECIF_RADAR").
    !
    ! Method: This is a wrapper routine, calling two other routines for rs_meta(ista)
    !         and dbz_meta(ista)
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    INTEGER,                INTENT(in) :: funit, ista
    TYPE(radar_meta_type),  INTENT(in) :: rs_m
    TYPE(radar_meta_type),  INTENT(in) :: rs_m_for_country(:)
    TYPE(t_dbzcalc_params), INTENT(in) :: dbz_m, dbz_m_d

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    TYPE(radar_meta_type)              :: rs_m_d
    INTEGER                            :: i

    CHARACTER (LEN=*), PARAMETER      :: yzroutine = 'emvorado::ctrl_output_ista_nuspec_fwo()'

    !---------------------------------------------------------------------------------
    ! .. Parameters for struct radar_meta_type:

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    rs_m_d = get_meta_proto( icountry=rs_m%icountry )
    IF (lreadmeta_from_netcdf) THEN
      DO i=1, nradsta_all
        IF (rs_m%station_id == rs_m_for_country(i)%station_id .AND. &
            rs_m%scanname(1:3) == rs_m_for_country(i)%scanname(1:3)) THEN
          rs_m_d = rs_m_for_country(i)
        END IF
      END DO
    END IF
    
    CALL ctrl_output_rsmeta_nuspec_fwo (funit, ista, rs_m , rs_m_d, 'rs_m' )
    CALL ctrl_output_dbzmeta_nuspec_fwo(funit, ista, dbz_m, dbz_m_d, 'dbz_m')

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE ctrl_output_ista_nuspec_fwo


  SUBROUTINE ctrl_output_rsmeta_nuspec_fwo(funit, ista, rs_m, rs_m_d, rsname)

    !------------------------------------------------------------------------------
    !
    ! Description: Control output of the radar namelist
    !              parameters in the structure radar_meta_type to a file
    !              (in COSMO per default to  "YUSPECIF_RADAR")
    !
    ! Method: Wrinting to file with unit "funit". File has to be opened before
    !         calling this subroutine!
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    INTEGER,               INTENT(in) :: funit, ista
    TYPE(radar_meta_type), INTENT(in) :: rs_m, rs_m_d
    CHARACTER (LEN=*)                 :: rsname

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER      :: yzroutine = 'emvorado::ctrl_output_rsmeta_nuspec_fwo()'
    CHARACTER (LEN=80)                :: cname
    INTEGER                           :: i, ii
    CHARACTER(len=10)                 :: cii

    CHARACTER(len=*), PARAMETER       :: format_char             = '(T8,A,"%",A,T46,A,T61,A,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_char_vec         = '(T8,A,"%",A,i3.3,A,T46,A14,T61,A14,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_int              = '(T8,A,"%",A,T46,I12,T61,I12,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_longname_int     = '(T4,A,"%",A,T46,I12,T61,I12,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_int_vec          = '(T8,A,"%",A,I3.3,A,T46,I12,T61,I12,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_longname_int_vec = '(T4,A,"%",A,i3.3,A,T46,I12,T61,I12,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f1          = '(T8,A,"%",A,T46,f12.1,T61,f12.1,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f1_vec      = '(T8,A,"%",A,i3.3,A,T46,f12.1,T61,f12.1,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f3          = '(T8,A,"%",A,T46,f12.3,T61,f12.3,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f3_vec      = '(T8,A,"%",A,i3.3,A,T46,f12.3,T61,f12.3,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_a_real_f3_vec    = '(T8,A,"%",A,i3.3,A,T46,A,T61,f12.3,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f3_a_vec    = '(T8,A,"%",A,i3.3,A,T46,F12.3,T61,A,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f5_vec      = '(T8,A,"%",A,i3.3,A,T46,f12.5,T61,f12.5,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_real_f6          = '(T8,A,"%",A,T46,f12.6,T61,f12.6,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_logical          = '(T8,A,"%",A,T46,L12,T61,L12,T76,A)'
    CHARACTER(len=*), PARAMETER       :: format_logical_vec      = '(T8,A,"%",A,i3.3,A,i2.2,A,T46,L12,T76,A)'
    
    !---------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    cname(:) = ' '
    IF (ista > 0) THEN
      WRITE(cname,'(a,i3.3,a)') TRIM(rsname)//'(',ista,')'
    ELSE
      WRITE(cname,'(a)') TRIM(rsname)
    END IF

    WRITE (funit, TRIM(format_char))                     &
         TRIM(cname), 'station_name',  rs_m%station_name, rs_m_d%station_name, ' C*3 '
    WRITE (funit, TRIM(format_char))                     &
         TRIM(cname), 'scanname',  rs_m%scanname, rs_m_d%scanname, ' C*10 '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'station_id', rs_m%station_id, rs_m_d%station_id, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'icountry', rs_m%icountry, rs_m_d%icountry, ' I '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'lambda', rs_m%lambda, rs_m_d%lambda, ' R '
    WRITE (funit, TRIM(format_real_f6))              &
         TRIM(cname), 'lat', rs_m%lat, rs_m_d%lat, ' R '
    WRITE (funit, TRIM(format_real_f6))              &
         TRIM(cname), 'lon', rs_m%lon, rs_m_d%lon, ' R '
    WRITE (funit, TRIM(format_real_f1))              &
         TRIM(cname), 'alt_agl_mod', rs_m%alt_agl_mod, rs_m_d%alt_agl_mod, ' R '
    WRITE (funit, TRIM(format_real_f1))              &
         TRIM(cname), 'alt_msl', rs_m%alt_msl, rs_m_d%alt_msl, ' R '
    WRITE (funit, TRIM(format_real_f1))              &
         TRIM(cname), 'msl_mod', rs_m%msl_mod, rs_m_d%msl_mod, ' R '
    WRITE (funit, TRIM(format_real_f1))              &
         TRIM(cname), 'alt_msl_true', rs_m%alt_msl_true, rs_m_d%alt_msl_true, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'i_nearest_mod', rs_m%i_nearest_mod, rs_m_d%i_nearest_mod, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'j_nearest_mod', rs_m%j_nearest_mod, rs_m_d%j_nearest_mod, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'nobs_times', rs_m%nobs_times, rs_m_d%nobs_times, ' I '
    DO i = 1, 3
      WRITE (funit, TRIM(format_real_f1_vec))              &
           TRIM(cname), 'dt_obs(',i,')', rs_m%dt_obs(i), rs_m_d%dt_obs(i), ' R '
    END DO
    DO i = 1, MAX(rs_m%nobs_times, 1)
      WRITE (funit, TRIM(format_real_f1_vec))       &
           TRIM(cname), 'obs_times(',i,')', rs_m%obs_times(i), rs_m_d%obs_times(i), ' R '
    END DO
    IF (lreadmeta_from_netcdf) THEN
      WRITE (funit, TRIM(format_logical))              &
           TRIM(cname), 'lobstimes_ovwrt_recalc', rs_m%lobstimes_ovwrt_recalc, rs_m_d%lobstimes_ovwrt_recalc, ' L '
      WRITE (funit, TRIM(format_int))                  &
           TRIM(cname), 'nobs_times_obs', rs_m%nobs_times_obs, rs_m_d%nobs_times_obs, ' I '
      DO i = 1, MAX(rs_m%nobs_times_obs, 1)
        WRITE (funit, TRIM(format_real_f1_vec))       &
             TRIM(cname), 'obs_times_obs(',i,')', rs_m%obs_times_obs(i), rs_m_d%obs_times_obs(i), ' R '
      END DO
      WRITE (funit, TRIM(format_logical))              &
           TRIM(cname), 'lvrad_to_fdbk', rs_m%lvrad_to_fdbk, rs_m_d%lvrad_to_fdbk, ' L '
      WRITE (funit, TRIM(format_real_f1))              &
           TRIM(cname), 'vnyq_min_for_vr_active_fdbk', rs_m%vnyq_min_for_vr_active_fdbk, rs_m_d%vnyq_min_for_vr_active_fdbk, ' R '
      WRITE (funit, TRIM(format_logical))              &
           TRIM(cname), 'ldbzh_to_fdbk', rs_m%ldbzh_to_fdbk, rs_m_d%ldbzh_to_fdbk, ' L '
    END IF
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'nel', rs_m%nel, rs_m_d%nel, ' I '
    DO i = 1, MAX(rs_m%nel,rs_m_d%nel)
      IF (i <= rs_m%nel .AND. i <= rs_m_d%nel) THEN
        WRITE (funit, TRIM(format_real_f3_vec))       &
             TRIM(cname), 'el_arr(',i,')', rs_m%el_arr(i), rs_m_d%el_arr(i), ' R '
      ELSE IF (i <= rs_m_d%nel) THEN
        WRITE (funit, TRIM(format_a_real_f3_vec))       &
             TRIM(cname), 'el_arr(',i,')', '-', rs_m_d%el_arr(i), ' R '
      ELSE
        WRITE (funit, TRIM(format_real_f3_a_vec))       &
             TRIM(cname), 'el_arr(',i,')', rs_m%el_arr(i), '-', ' R '
      END IF
    END DO
    WRITE (funit, TRIM(format_longname_int))                  &
         TRIM(cname), 'eleind_for_composite_bub', rs_m%eleind_for_composite_bub, rs_m_d%eleind_for_composite_bub, ' I '
    DO i = 1, MAX(nel_composite, 1)
      WRITE (funit, TRIM(format_longname_int_vec))     &
           TRIM(cname), 'eleindlist_for_composite(',i,')', rs_m%eleindlist_for_composite(i), &
                                                         rs_m_d%eleindlist_for_composite(i), ' I '
    END DO
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'az_start', rs_m%az_start, rs_m_d%az_start, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'naz', rs_m%naz, rs_m_d%naz, ' I '
    IF (lreadmeta_from_netcdf) THEN
      WRITE (funit, TRIM(format_int))                &
           TRIM(cname), 'naz_ncdf (vr)', rs_m%naz_ncdf(1,1), rs_m_d%naz_ncdf(1,1), ' I '
      WRITE (funit, TRIM(format_int))                &
           TRIM(cname), 'naz_ncdf (qv)', rs_m%naz_ncdf(1,2), rs_m_d%naz_ncdf(1,2), ' I '
      WRITE (funit, TRIM(format_int))                &
           TRIM(cname), 'naz_ncdf ( z)', rs_m%naz_ncdf(1,3), rs_m_d%naz_ncdf(1,3), ' I '
      WRITE (funit, TRIM(format_int))                &
           TRIM(cname), 'naz_ncdf (qz)', rs_m%naz_ncdf(1,4), rs_m_d%naz_ncdf(1,4), ' I '
    END IF
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'az_inc', rs_m%az_inc, rs_m_d%az_inc, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'nra', rs_m%nra, rs_m_d%nra,          ' I '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'ra_inc', rs_m%ra_inc, rs_m_d%ra_inc, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'ra_inc_obs', rs_m%ra_inc_obs, rs_m_d%ra_inc_obs, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'ra_inc_coarse', rs_m%ra_inc_coarse, rs_m_d%ra_inc_coarse, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'n_aggr_ra_obs', rs_m%n_aggr_ra_obs, rs_m_d%n_aggr_ra_obs, ' I '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'mds_Z0', rs_m%mds_Z0, rs_m_d%mds_Z0, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'mds_r0', rs_m%mds_r0, rs_m_d%mds_r0, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Phi3', rs_m%Phi3, rs_m_d%Phi3,       ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Theta3', rs_m%Theta3, rs_m_d%Theta3, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'dalpha', rs_m%dalpha, rs_m_d%dalpha, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'alpha3_eff_0', rs_m%alpha3_eff_0, rs_m_d%alpha3_eff_0, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'smth_interv_fact', rs_m%smth_interv_fact, rs_m_d%smth_interv_fact, ' R '

    ! Print ext_nyq only for the first 3 obs times in a row:
    ii = MIN(MAX(rs_m%nobs_times, 1), 3)
    cii = ' '
    WRITE (cii,'(i0)') ii
    DO i = 1, MAX(rs_m%nel,rs_m_d%nel)
      IF (i <= rs_m%nel .AND. i <= rs_m_d%nel) THEN
        WRITE (funit, '(T8,A,"%",A,i3.3,A,i1.1,A,T46,'//TRIM(cii)//'f12.3," ...",T91,'//TRIM(cii)//'f12.3," ...",T135,A)') &
             TRIM(cname), 'ext_nyq(',i,',1:',ii,')', rs_m%ext_nyq(i,1:ii), rs_m_d%ext_nyq(i,1:ii), ' R '
      ELSE IF (i <= rs_m_d%nel) THEN
        WRITE (funit, '(T8,A,"%",A,i3.3,A,i1.1,A,T46,'//TRIM(cii)//'("-",11x)," ...",T91,'//TRIM(cii)//'f12.3," ...",T135,A)') &
             TRIM(cname), 'ext_nyq(',i,',1:',ii,')', rs_m_d%ext_nyq(i,1:ii), ' R '
      ELSE
        WRITE (funit, '(T8,A,"%",A,i3.3,A,i1.1,A,T46,'//TRIM(cii)//'f12.3," ...",T91,'//TRIM(cii)//'("-",11x)," ...",T135,A)') &
             TRIM(cname), 'ext_nyq(',i,',1:',ii,')', rs_m_d%ext_nyq(i,1:ii), ' R '
        WRITE (funit, TRIM(format_real_f3_a_vec))       &
             TRIM(cname), 'ext_nyq(',i,',1)', rs_m%ext_nyq(i,1), ' R '
      END IF
    END DO

    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'high_nyq(ele1)', rs_m%high_nyq(1), rs_m_d%high_nyq(1), ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'prf(ele1)', rs_m%prf(1), rs_m_d%prf(1),          ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'dualprf_ratio(ele1)', rs_m%dualprf_ratio(1), rs_m_d%dualprf_ratio(1), ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'rngate_len', rs_m%rngate_len, rs_m_d%rngate_len, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'num_gates', rs_m%num_gates, rs_m_d%num_gates, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'num_pulses', rs_m%num_pulses, rs_m_d%num_pulses, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'ngpsm_h', rs_m%ngpsm_h, rs_m_d%ngpsm_h, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'ngpsm_v', rs_m%ngpsm_v, rs_m_d%ngpsm_v, ' I '
    DO i = 1, rs_m%ngpsm_h
      WRITE (funit, TRIM(format_real_f5_vec))       &
           TRIM(cname), 'xabscsm_h(',i,')', rs_m%xabscsm_h(i), rs_m_d%xabscsm_h(i), ' R '
    END DO
    DO i = 1, rs_m%ngpsm_v
      WRITE (funit, TRIM(format_real_f5_vec))       &
           TRIM(cname), 'xabscsm_v(',i,')', rs_m%xabscsm_v(i), rs_m_d%xabscsm_v(i), ' R '
    END DO
    DO i = 1, rs_m%ngpsm_h
      WRITE (funit, TRIM(format_real_f5_vec))       &
           TRIM(cname), 'weigsm_h(',i,')', rs_m%weigsm_h(i), rs_m_d%weigsm_h(i), ' R '
    END DO
    DO i = 1, rs_m%ngpsm_v
      WRITE (funit, TRIM(format_real_f5_vec))       &
           TRIM(cname), 'weigsm_v(',i,')', rs_m%weigsm_v(i), rs_m_d%weigsm_v(i), ' R '
    END DO
    
    IF (lreadmeta_from_netcdf) THEN
      DO i = 1, rs_m%nobs_times
        WRITE (funit, TRIM(format_char_vec))         &
             TRIM(cname), 'obs_cdate(',i,')', rs_m%obs_cdate(i), rs_m_d%obs_cdate(i), ' C*14 '
      END DO
      DO i = 1, rs_m%nobs_times
        WRITE (funit, '(T8,A,"%",A,i3.3,A,T46,'//TRIM(cndatakind)//'I10)') &
             TRIM(cname), 'obs_startrec(',i,',1:'//TRIM(cndatakind)//')', rs_m%obs_startrec(i,1:ndatakind)
        WRITE (funit, '(T8,A,"%",A,i3.3,A,T46,'//TRIM(cndatakind)//'I10)') &
             TRIM(cname), 'obs_endrec  (',i,',1:'//TRIM(cndatakind)//')', rs_m%obs_endrec(i,1:ndatakind)
      END DO
      DO ii=1,ndatakind
        WRITE (funit, '(T8,A,"%",A,"(",i2.2,")",T46,I12,T61,I12,T76,A)')   &
             TRIM(cname), 'nrep_ncdf', ii, rs_m%nrep_ncdf(ii), rs_m_d%nrep_ncdf(ii), ' I '
      END DO
      DO ii = 1, ndatakind
        DO i = 1, rs_m%nobs_times
          WRITE (funit, '(T8,A,"%",A,i3.3,A,i2.2,A,T46,A,T65,A)')                                &
               TRIM(cname), 'obsfile(',i,',',ii,')', TRIM(ADJUSTL(rs_m%obsfile_format(i))), TRIM(ADJUSTL(rs_m%obsfile(i,ii)))
          WRITE (funit, TRIM(format_logical_vec))                        &
               TRIM(cname), 'lobs_avail(',i,',',ii,')', rs_m%lobs_avail(i,ii), ' L '
        END DO
      END DO
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_vrad', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_vrad))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_dbzh', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_dbzh))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_zdr', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_zdr))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_rhv', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_rhv))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_kdp', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_kdp))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_phidp', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_phidp))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_ldr', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_ldr))
      WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
           TRIM(cname), 'obs_hdf5_varname_cflag', TRIM(ADJUSTL(rs_m%obs_hdf5_varname_cflags))

      IF (lfdbk_output) THEN
        WRITE (funit, '(T8,A,"%",A,T46,A)')                                &
             TRIM(cname), 'fdbkfile', TRIM(ADJUSTL(rs_m%fdbkfile))
        WRITE (funit, TRIM(format_int))                  &
             TRIM(cname), 'nobs_times_fdbk', rs_m%nobs_times_fdbk, rs_m_d%nobs_times_fdbk, ' I '
        DO i = 1, 3
          WRITE (funit, TRIM(format_real_f1_vec))              &
               TRIM(cname), 'dt_obs_fdbk(',i,')', rs_m%dt_obs_fdbk(i), rs_m_d%dt_obs_fdbk(i), ' R '
        END DO
        DO i = 1, MAX(rs_m%nobs_times_fdbk, 1)
          WRITE (funit, TRIM(format_real_f1_vec))       &
               TRIM(cname), 'obs_times_fdbk(',i,')', rs_m%obs_times_fdbk(i), rs_m_d%obs_times_fdbk(i), ' R '
        END DO
        DO i = 1, MAX(rs_m%nel_fdbk, 1)
          WRITE (funit, TRIM(format_int_vec))     &
               TRIM(cname), 'ind_ele_fdbk(',i,')', rs_m%ind_ele_fdbk(i), rs_m_d%ind_ele_fdbk(i), ' R '
        END DO
      END IF
    END IF
    
    IF (lvoldata_output) THEN
      WRITE (funit, TRIM(format_int))                  &
           TRIM(cname), 'nobs_times_voldata', rs_m%nobs_times_voldata, rs_m_d%nobs_times_voldata, ' I '
      DO i = 1, 3
        WRITE (funit, TRIM(format_real_f1_vec))              &
             TRIM(cname), 'dt_obs_voldata(',i,')', rs_m%dt_obs_voldata(i), rs_m_d%dt_obs_voldata(i), ' R '
      END DO
      DO i = 1, MAX(rs_m%nobs_times_voldata, 1)
        WRITE (funit, TRIM(format_real_f1_vec))       &
             TRIM(cname), 'obs_times_voldata(',i,')', rs_m%obs_times_voldata(i), rs_m_d%obs_times_voldata(i), ' R '
      END DO
      DO i = 1, MAX(rs_m%nel_voldata, 1)
        WRITE (funit, TRIM(format_int_vec))     &
             TRIM(cname), 'ind_ele_voldata(',i,')', rs_m%ind_ele_voldata(i), rs_m_d%ind_ele_voldata(i), ' R '
      END DO
    END IF


    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE ctrl_output_rsmeta_nuspec_fwo


  SUBROUTINE ctrl_output_dbzmeta_nuspec_fwo(funit, ista, dbz_m, dbz_m_d, dbzname)

    !------------------------------------------------------------------------------
    !
    ! Description: Control output of the namelist
    !              parameters in the structure dbzcalc_params to a file
    !              (in COSMO per default to  "YUSPECIF_RADAR")
    !
    ! Method: Wrinting to file with unit "funit". File has to be opened before
    !         calling this subroutine!
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    IMPLICIT NONE

    INTEGER,                INTENT(in) :: funit, ista
    TYPE(t_dbzcalc_params), INTENT(in) :: dbz_m, dbz_m_d
    CHARACTER (LEN=*)                  :: dbzname

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=*), PARAMETER       :: yzroutine = 'emvorado::ctrl_output_dbzmeta_nuspec_fwo()'
    CHARACTER (LEN=80)                 :: cname

    CHARACTER(len=*), PARAMETER        :: format_char        = '(T8,A,"%",A,T46,A,T61,A,T76,A)'
    CHARACTER(len=*), PARAMETER        :: format_int         = '(T8,A,"%",A,T46,I12,T61,I12,T76,A)'
    CHARACTER(len=*), PARAMETER        :: format_real_f3     = '(T8,A,"%",A,T46,f12.3,T61,f12.3,T76,A)'
    CHARACTER(len=*), PARAMETER        :: format_real_e3     = '(T8,A,"%",A,T46,e12.3,T61,e12.3,T76,A)'
    CHARACTER(len=*), PARAMETER        :: format_logical     = '(T8,A,"%",A,T46,L12,T61,L12,T76,A)'
    
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    cname(:) = ' '
    IF (ista > 0) THEN
      WRITE(cname,'(a,i3.3,a)') TRIM(dbzname)//'(',ista,')'
    ELSE
      WRITE(cname,'(a)') TRIM(dbzname)
    END IF

    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'station_id', dbz_m%station_id, dbz_m_d%station_id, ' I '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'lambda_radar', dbz_m%lambda_radar, dbz_m_d%lambda_radar, ' R '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname), 'itype_refl', dbz_m%itype_refl, dbz_m_d%itype_refl, ' I '
    WRITE (funit, TRIM(format_logical))              &
         TRIM(cname), 'ldynamic_wetgrowth_gh', dbz_m%ldynamic_wetgrowth_gh, dbz_m_d%ldynamic_wetgrowth_gh, ' L '
    WRITE (funit, TRIM(format_logical))              &
         TRIM(cname), 'llookup_mie', dbz_m%llookup_mie, dbz_m_d%llookup_mie, ' L '
    WRITE (funit, '(T8,A,"%",A,T46,6L2,T61,6L2,T76,A)')                  &
         TRIM(cname), 'lhydrom_choice_testing', dbz_m%lhydrom_choice_testing, dbz_m_d%lhydrom_choice_testing, ' L '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname),  'isnow_type', dbz_m%isnow_type, dbz_m_d%isnow_type, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname),  'igraupel_type', dbz_m%igraupel_type, dbz_m_d%igraupel_type, ' I '
    WRITE (funit, TRIM(format_int))                  &
         TRIM(cname),  'itype_Dref_fmelt', dbz_m%itype_Dref_fmelt, dbz_m_d%itype_Dref_fmelt, ' I '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmeltbegin_i', dbz_m%Tmeltbegin_i, dbz_m_d%Tmeltbegin_i, ' R '
    !WRITE (funit, TRIM(format_real_f3))              &
    !     TRIM(cname), 'Tmin_i', dbz_m%Tmin_i, dbz_m_d%Tmin_i, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'meltdegTmin_i', dbz_m%meltdegTmin_i, dbz_m_d%meltdegTmin_i, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_min_i', dbz_m%Tmax_min_i, dbz_m_d%Tmax_min_i, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_max_i', dbz_m%Tmax_max_i, dbz_m_d%Tmax_max_i, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qthresh_i', dbz_m%qthresh_i, dbz_m_d%qthresh_i, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qnthresh_i', dbz_m%qnthresh_i, dbz_m_d%qnthresh_i, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmeltbegin_s', dbz_m%Tmeltbegin_s, dbz_m_d%Tmeltbegin_s, ' R '
    !WRITE (funit, TRIM(format_real_f3))              &
    !     TRIM(cname), 'Tmin_s', dbz_m%Tmin_s, dbz_m_d%Tmin_s, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'meltdegTmin_s', dbz_m%meltdegTmin_s, dbz_m_d%meltdegTmin_s, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_min_s', dbz_m%Tmax_min_s, dbz_m_d%Tmax_min_s, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_max_s', dbz_m%Tmax_max_s, dbz_m_d%Tmax_max_s, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qthresh_s', dbz_m%qthresh_s, dbz_m_d%qthresh_s, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qnthresh_s', dbz_m%qnthresh_s, dbz_m_d%qnthresh_s, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmeltbegin_g', dbz_m%Tmeltbegin_g, dbz_m_d%Tmeltbegin_g, ' R '
    !WRITE (funit, TRIM(format_real_f3))              &
    !     TRIM(cname), 'Tmin_g', dbz_m%Tmin_g, dbz_m_d%Tmin_g, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'meltdegTmin_g', dbz_m%meltdegTmin_g, dbz_m_d%meltdegTmin_g, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_min_g', dbz_m%Tmax_min_g, dbz_m_d%Tmax_min_g, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_max_g', dbz_m%Tmax_max_g, dbz_m_d%Tmax_max_g, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qthresh_g', dbz_m%qthresh_g, dbz_m_d%qthresh_g, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qnthresh_g', dbz_m%qnthresh_g, dbz_m_d%qnthresh_g, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmeltbegin_h', dbz_m%Tmeltbegin_h, dbz_m_d%Tmeltbegin_h, ' R '
    !WRITE (funit, TRIM(format_real_f3))              &
    !     TRIM(cname), 'Tmin_h', dbz_m%Tmin_h, dbz_m_d%Tmin_h, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'meltdegTmin_h', dbz_m%meltdegTmin_h, dbz_m_d%meltdegTmin_h, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_min_h', dbz_m%Tmax_min_h, dbz_m_d%Tmax_min_h, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'Tmax_max_h', dbz_m%Tmax_max_h, dbz_m_d%Tmax_max_h, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qthresh_h', dbz_m%qthresh_h, dbz_m_d%qthresh_h, ' R '
    WRITE (funit, TRIM(format_real_e3))              &
         TRIM(cname), 'qnthresh_h', dbz_m%qnthresh_h, dbz_m_d%qnthresh_h, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'ext_tune_fac_pure', dbz_m%ext_tune_fac_pure, dbz_m_d%ext_tune_fac_pure, ' R '
    WRITE (funit, TRIM(format_real_f3))              &
         TRIM(cname), 'ext_tune_fac_melt', dbz_m%ext_tune_fac_melt, dbz_m_d%ext_tune_fac_melt, ' R '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_drysnow_mie',  dbz_m%ctype_drysnow_mie, dbz_m_d%ctype_drysnow_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wetsnow_mie',  dbz_m%ctype_wetsnow_mie, dbz_m_d%ctype_wetsnow_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_drygraupel_mie',  dbz_m%ctype_drygraupel_mie, dbz_m_d%ctype_drygraupel_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wetgraupel_mie',  dbz_m%ctype_wetgraupel_mie, dbz_m_d%ctype_wetgraupel_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_dryhail_mie',  dbz_m%ctype_dryhail_mie, dbz_m_d%ctype_dryhail_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wethail_mie',  dbz_m%ctype_wethail_mie, dbz_m_d%ctype_wethail_mie, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_drysnow_ray',  dbz_m%ctype_drysnow_ray, dbz_m_d%ctype_drysnow_ray, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wetsnow_ray',  dbz_m%ctype_wetsnow_ray, dbz_m_d%ctype_wetsnow_ray, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_drygraupel_ray', dbz_m%ctype_drygraupel_ray, dbz_m_d%ctype_drygraupel_ray, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wetgraupel_ray', dbz_m%ctype_wetgraupel_ray, dbz_m_d%ctype_wetgraupel_ray, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_dryhail_ray', dbz_m%ctype_dryhail_ray, dbz_m_d%ctype_dryhail_ray, ' C*12 '
    WRITE (funit, TRIM(format_char))                 &
         TRIM(cname), 'ctype_wethail_ray', dbz_m%ctype_wethail_ray, dbz_m_d%ctype_wethail_ray, ' C*12 '

    IF (dbz_m%itype_refl .GT. 4) THEN
      WRITE (funit, TRIM(format_char))                 &
           TRIM(cname), 'polMP_r%ARmodel', dbz_m%polMP_r%ARmodel, dbz_m_d%polMP_r%ARmodel, ' C*5 '
      IF (tolower(TRIM(dbz_m%polMP_r%ARmodel)) .EQ. 'poly') THEN
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_r%c0', dbz_m%polMP_r%c0, dbz_m_d%polMP_r%c0, ' R '
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_r%c1', dbz_m%polMP_r%c1, dbz_m_d%polMP_r%c1, ' R '
      END IF
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_r%ARmin', dbz_m%polMP_r%ARmin, dbz_m_d%polMP_r%ARmin, ' R '
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_r%sig', dbz_m%polMP_r%sig, dbz_m_d%polMP_r%sig, ' R '
      
      WRITE (funit, TRIM(format_char))                 &
           TRIM(cname), 'polMP_i%ARmodel', dbz_m%polMP_i%ARmodel, dbz_m_d%polMP_i%ARmodel, ' C*5 '
      IF (tolower(TRIM(dbz_m%polMP_i%ARmodel)) .EQ. 'poly') THEN
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_i%c0', dbz_m%polMP_i%c0, dbz_m_d%polMP_i%c0, ' R '
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_i%c1', dbz_m%polMP_i%c1, dbz_m_d%polMP_i%c1, ' R '
      END IF
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_i%ARmin', dbz_m%polMP_i%ARmin, dbz_m_d%polMP_i%ARmin, ' R '
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_i%sig', dbz_m%polMP_i%sig, dbz_m_d%polMP_i%sig, ' R '

      WRITE (funit, TRIM(format_char))                 &
           TRIM(cname), 'polMP_s%ARmodel', dbz_m%polMP_s%ARmodel, dbz_m_d%polMP_s%ARmodel, ' C*5 '
      IF (tolower(TRIM(dbz_m%polMP_s%ARmodel)) .EQ. 'poly') THEN
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_s%c0', dbz_m%polMP_s%c0, dbz_m_d%polMP_s%c0, ' R '
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_s%c1', dbz_m%polMP_s%c1, dbz_m_d%polMP_s%c1, ' R '
      END IF
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_s%ARmin', dbz_m%polMP_s%ARmin, dbz_m_d%polMP_s%ARmin, ' R '
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_s%sig', dbz_m%polMP_s%sig, dbz_m_d%polMP_s%sig, ' R '

      WRITE (funit, TRIM(format_char))                 &
           TRIM(cname), 'polMP_g%ARmodel', dbz_m%polMP_g%ARmodel, dbz_m_d%polMP_g%ARmodel, ' C*5 '
      IF (tolower(TRIM(dbz_m%polMP_g%ARmodel)) .EQ. 'poly') THEN
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_g%c0', dbz_m%polMP_g%c0, dbz_m_d%polMP_g%c0, ' R '
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_g%c1', dbz_m%polMP_g%c1, dbz_m_d%polMP_g%c1, ' R '
      END IF
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_g%ARmin', dbz_m%polMP_g%ARmin, dbz_m_d%polMP_g%ARmin, ' R '
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_g%sig', dbz_m%polMP_g%sig, dbz_m_d%polMP_g%sig, ' R '

      WRITE (funit, TRIM(format_char))                 &
           TRIM(cname), 'polMP_h%ARmodel', dbz_m%polMP_h%ARmodel, dbz_m_d%polMP_h%ARmodel, ' C*5 '
      IF (tolower(TRIM(dbz_m%polMP_h%ARmodel)) .EQ. 'poly') THEN
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_h%c0', dbz_m%polMP_h%c0, dbz_m_d%polMP_h%c0, ' R '
        WRITE (funit, TRIM(format_real_e3))              &
             TRIM(cname), 'polMP_h%c1', dbz_m%polMP_h%c1, dbz_m_d%polMP_h%c1, ' R '
      END IF
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_h%ARmin', dbz_m%polMP_h%ARmin, dbz_m_d%polMP_h%ARmin, ' R '
      WRITE (funit, TRIM(format_real_f3))              &
           TRIM(cname), 'polMP_h%sig', dbz_m%polMP_h%sig, dbz_m_d%polMP_h%sig, ' R '
    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE ctrl_output_dbzmeta_nuspec_fwo

!===============================================================================
!===============================================================================

  SUBROUTINE init_smth_info
    !------------------------------------------------------------------------------
    !
    ! Description: Initializes some parameters needed for the
    !              spatial smoothing (short:smth)  with Gauss-Legendre quadrature
    !
    !
    ! Method:
    !
    ! Input: Some namelist parameters from namelist RADARSIM_PARAMS
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'emvorado::init_smth_info()'
    CHARACTER (LEN=80) :: yzerrmsg
    INTEGER            :: i, ista

    ! Integration points: vertical and horizontal, which means elevation
    !                     and azimuth in radar system
    REAL (KIND=dp), ALLOCATABLE    :: xabsc_v(:),xabsc_h(:)
    ! Weights of integration points
    REAL (KIND=dp), ALLOCATABLE    :: weig_v(:) ,weig_h(:)

    CHARACTER(len=cmaxlen) :: errstring



    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    DO ista = 1, nradsta

      rs_meta(ista)%xabscsm_v = 0.0_dp
      rs_meta(ista)%weigsm_v  = 0.0_dp

      ALLOCATE(xabsc_v(rs_meta(ista)%ngpsm_v))
      ALLOCATE(weig_v(rs_meta(ista)%ngpsm_v))

      CALL gauleg(rs_meta(ista)%ngpsm_v,xabsc_v,weig_v)

      rs_meta(ista)%xabscsm_v(1:rs_meta(ista)%ngpsm_v) = xabsc_v
      rs_meta(ista)%weigsm_v(1:rs_meta(ista)%ngpsm_v)  = weig_v

      IF (ldebug_radsim) &
           WRITE(*,*) ista, '  xabsc_v = ', (xabsc_v(i),i=1,rs_meta(ista)%ngpsm_v)


      rs_meta(ista)%xabscsm_h = 0.0_dp
      rs_meta(ista)%weigsm_h  = 0.0_dp

      ALLOCATE(xabsc_h(rs_meta(ista)%ngpsm_h))
      ALLOCATE(weig_h(rs_meta(ista)%ngpsm_h))

      CALL gauleg(rs_meta(ista)%ngpsm_h,xabsc_h,weig_h)

      rs_meta(ista)%xabscsm_h(1:rs_meta(ista)%ngpsm_h) = xabsc_h
      rs_meta(ista)%weigsm_h(1:rs_meta(ista)%ngpsm_h)  = weig_h

      IF (ldebug_radsim) &
           WRITE(*,*) ista, '  xabsc_h = ', (xabsc_h(i),i=1,rs_meta(ista)%ngpsm_h)

      DEALLOCATE(xabsc_v)
      DEALLOCATE(xabsc_h)
      DEALLOCATE(weig_v)
      DEALLOCATE(weig_h)

    END DO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE init_smth_info

  SUBROUTINE  gauleg(ngp, xabsc, weig)

    INTEGER             :: i,j,m
    REAL (KIND=dp)      ::  p1, p2, p3, pp, z, z1
    INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
    REAL (KIND=dp),     INTENT(OUT):: xabsc(ngp), weig(ngp)
    REAL(KIND=dp), PARAMETER :: pi     = 4.0_dp * ATAN(1.0_dp)

    m = (ngp + 1) / 2
    !* Roots are symmetric in the interval - so only need to find half of them  */

    DO i = 1, m                 ! Loop over the desired roots */

      z = COS( pi * (i-0.25_dp) / (ngp+0.5_dp) )
      !*   Starting with the above approximation to the ith root,
      !*          we enter the main loop of refinement by NEWTON'S method   */
      LOOP: DO

        ! 100 p1 = 1.0_dp
        p1 = 1.0_dp
        p2 = 0.0_dp
        !*  Loop up the recurrence relation to get the Legendre
        !*  polynomial evaluated at z                 */

        DO j = 1, ngp
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*j-1.0_dp) * z * p2 - (j-1.0_dp)*p3) / j
        ENDDO

        !* p1 is now the desired Legendre polynomial. We next compute pp,
        !* its derivative, by a standard relation involving also p2, the
        !* polynomial of one lower order.      */
        pp = ngp*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z = z1 - p1/pp                          ! Newton's Method  */

        IF (dabs(z-z1) .GT. 0.000000000000001_dp) CYCLE LOOP
        EXIT LOOP

      END DO LOOP

      xabsc(i) =  - z                           ! Roots will be bewteen -1.0 & 1.0 */
      xabsc(ngp+1-i) =  + z                     ! and symmetric about the origin  */
      weig(i) = 2.0_dp/((1.0_dp-z*z)*pp*pp) ! Compute the weight and its       */
      weig(ngp+1-i) = weig(i)                   ! symmetric counterpart         */

    END DO     ! i loop

  END SUBROUTINE gauleg

!===============================================================================
!===============================================================================


END MODULE radar_namelist_read
