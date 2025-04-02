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


!!! IDEE: replace ndoms_max by actual number of domains ndoms

#if defined TWOMOM_SB_NEW || defined TWOMOM_SB_OLD
#define TWOMOM_SB
#endif

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_organize

!------------------------------------------------------------------------------
!
! Description: Main organizational routine of the Radar forward operator
!              EMVORADO for calculation of simulated radar
!              quantities from a number of radar stations based on the COSMO
!              model variables.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------

!================================================================================
!
! Developer's notes:
! ------------------
!
! TODO:
! ... Correct timing of EMVORADO in COSMO for the new asynchroneous radar IO
! ... Simplify filename convention for radar obs: require the country name as a suffix (create links before running EMVORADO):
!     - DWD: cdfin*.dwd
!     - Swiss: [MP]LA152050000??.001.?.h5.mch
!     - Italy: odim_*_01.h5.arpa
! ... Consolidation of missing_obs and -999.99 values for all variables! replace -999.99 by missing_obs, where appropriate!
!      Only missing_obs should be used outside of data readings, because missing_z, missing_vr, ..., for different
!      stations/countries might be different and needs to be separated from the rest of the code.
! ... Rename global icountry --> icountry_glob in namelist
! ... Loesung fuer 2D Arrays in radar_mie*90: 3D-pointer (:,1,:) auf (:,:) (i,k)?
! ... Bestimmung von Tmax: Modellspezifisches in radar_mie_cosmo_iface.f90

!!$ Changes 09/2018:
!!$ - Eliminated linput_old_format together with the old cdfin_XX_NN input format and related procedures.
!!$ - Asynchroneous radar IO mode implemented, also in COSMO model
!!$ - ireals --> wp
!!$ - renamed some modules of EMVORADO, so that now all module names start with "radar"
!!$ - Reduced rs_meta-structure for single-file meta data reading
!!$ - Register new reduced default scan-strategy no. 9 (0.5, 1.5, 3.5) and 10 for operational radialwind DA
!!$ - Rename namelist filename to INPUT_RADARSIM
!!$ - Add missing registration of ASR Borkum in get_meta_network_dwd

!!$ Changes 30.4.2018:
!!$ - rs_meta%az_start is now hardcoded the center of the sampling interval of the first nominal output azimut,
!!$    and it is no more read from the "azoff" variable in the cdfin-files. It has to be specified
!!$    via the background meta data list or via the namelist.
!!$   The real radar azimuts are mapped to these nominal azimuts via nearest neighbour, without any interpolation.
!!$ - rs_meta%ext_nyq, %high_nyq, %prf, %dualprf_ratio are now a function of elevation index, because of
!!$    the Swiss radars

!!$ New: for the "cdfin" voldata output format and if linput_old_format=.false.,
!!$    it is possible to organize the output files in time batches, the same way as
!!$    it is the case with cdfin-input files if linput_old_format=.false.
!!$   Two new namelist parameters have been added to control this:
!!$
!!$      "cdfin_dt"    :   length of time batches in output "cdfin"-files [s]
!!$      "cdfin_tref"  :   reference time for synchronization of the batches [s since model start]
!!$
!!$   For example, if we want to write 5-min radar data in a half-hour forecast run and set
!!$    cdfin_dt = 900.0 and cdfin_tref=300.0, we get the following files:
!!$
!!$      "cdfin_zrsim_id-006661_200303210000_200303210000_volscan"
!!$      "cdfin_zrsim_id-006661_200303210005_200303210015_volscan"
!!$      "cdfin_zrsim_id-006661_200303210020_200303210030_volscan"

!!$ New: If precip scans are simulated (also alongside normal volume scans), it is
!!$    possible to use them for composites. To do so, one has to specify an
!!$    elevation index of "98" in the composite index list "eleindlist_for_composite_glob"
!!$    or "eleind_for_composite_bub_glob" for the bubble generator.

!!$ New: For composites from volume scan elevation data, it is possible to
!!$    output a vertical column maximum composite over all elevations. To do so, one
!!$    has to specify an elevation index of "99" in the composite index
!!$    list "eleindlist_for_composite_glob" or "eleind_for_composite_bub_glob" for the bubble generator.

!!$ New: Simulation of DWD precip scans implemented, based on the below mentioned data files
!!$
!!$        "cdfin_z_id-010908_201307282200_201307282255_precipscan"
!!$
!!$   Only possible for "linput_old_format=.FALSE.", also in case of lreadmeta_from_netcdf=.FALSE.!
!!$
!!$   The azimut-dependent elevations for each DWD station are stored in the new .incf-file
!!$       "src/radar_elevations_precipscan.incf"
!!$    which is a F90-source file autogenerated by the script "format_precipscan_f90",
!!$    using the text file "elevations_precipscan.txt" as input. The latter has been sent
!!$    by M. Mott in 2017. Script and text file are in the git-repo "cosmo-emvorado-package"
!!$    (git clone git@git.mpimet.mpg.de:cosmo-emvorado-package.git)
!!$
!!$   It is possible to simulate PPI volume scans and precipitation scans in one
!!$    model run. The input files for PPI scan data and precip scan data have to
!!$    be in the same input directory "ydirradarin".
!!$
!!$   Precip- and volume scans of the same station are internally treated as two
!!$    different radar stations. The same would hold if there were different
!!$    volume scan strategies in different data files for the same station.
!!$   To identify such differences for the same station, a "scanname" has been added
!!$    to the rs_meta(ista) meta data structure. It is either "PPIxxxx" or "PRECIP".
!!$    "xxxx" stands for the mean volume scan elevation times 10, e.g.,
!!$    "PPI0080" for the default 10 elevation DWD scan strategy, which has an average
!!$    elevation of 8.0 degrees.
!!$   This scanname is automatically determined from the given scan strategy, where
!!$    a precipitation scan is identified by rs_meta%icountry=1, rs_meta%nel=1, and
!!$    rs_meta%el_arr(1)=0.4, 0.59, 0.8 or 1.3 (for DWD obs data from > 2016, this is always 0.8).
!!$
!!$   It is possible to simulate precip scans without actual data files (noobs- or idealized modes).
!!$   For this, one has to specify a known DWD station ID in an rs_meta(i) - section in the namelist
!!$    and at the same time rs_meta(i)%icountry=1, rs_meta(i)%nel=1, rs_meta(i)%el_arr=0.8.
!!$   Example:
!!$            rs_meta(1)%station_id = 10908
!!$            rs_meta(1)%icountry   = 1
!!$            rs_meta(1)%nel        = 1
!!$            rs_meta(1)%el_arr     = 0.8
!!$
!!$   As a side effect, the output file names for fof's and ASCII-files have been
!!$    extended by the above scanname, which yields to the following
!!$    examples of the new file names:
!!$     "zrsim_id-010950_PPI0080_20130728120000_00001500_polar.dat"
!!$     "zrsim_id-010950_PRECIP_20130728120000_00001500_polar.dat"
!!$     "fof_radar_id-010871_PPI0080_20130728120000.nc"
!!$     "fof_radar_id-010871_PRECIP_20130728120000.nc"
!!$

!!$ New: Switch "linput_old_format" to choose between "classic" cdfin-format (linput_old_format=.TRUE., which is the default)
!!$  and a new format (linput_old_format=.FALSE.) where the "classic" files can be optionally splitted up according to the obs_times.
!!$  and follow a new naming convention. Here, a station can have more than one file per datakind (vr, qv, z, qz).
!!$
!!$  The new naming convention follows the examples
!!$
!!$        "cdfin_z_id-10908_201307282200_201307282255_volscan"     or
!!$        "cdfin_z_id-10908_201307282200_201307282255_precipscan"
!!$
!!$ Pattern:   cdfin_<"vr"|"qv"|"z"|"qz">_<station_id>_<obs_startdate>_<obs_enddate>_<"volscan"|"precipscan">
!!$
!!$ Both the "classic" and new cdfin-format is the result of converting BUFR-data
!!$  from DWD database "sky" to netcdf using the utility "bufrx2netcdf".
!!$
!!$ To generate the new cdfin-format, the option -b <minutes> has been added to the
!!$  script "~/ublahak/BIN/get_radbufr_data.sky" to split up the cdfins into time batches
!!$  of <minutes>. This script now always leads to the above new cdfin filename convention,
!!$  even if -b is not specified.
!!$ It is still possible to read cdfin-files having all the obs_times in one file.

!!$ New: different composites now possible for different elevations
!!$  new namelist parameters:
!!$  - eleind_for_composites_bub_glob
!!$  - nel_composite
!!$  - eleindlist_for_composites_glob(1:nel_composite)

!!$ Set very small absolute values of radial wind (< eps_vr=1e-12) to 0, otherwise results are not reproducible
!!$  for different processor configs! Not really clear, why this happens, but it seems somehow not related
!!$  to the processor grid but to MPI communication truncation error or something like that ...

!!$  New #ifdef WITH_ZLIB in SUBROUTINE output3d_ascii_radar() to enable compression of ascii files
!!$  (if voldata_format = 'ascii-gzip' in namelist)

!!$ New switch lfill_vr_backgroundwind
!!$
!!$ To take into account refl. weighted term. fallspeed at lsmooth=.false. and get vr everywhere:
!!$  lweightdbz = .true., fall = .true., lmds_vr = .true., lfill_vr_backgroundwind = .true.
!!$
!!$ To take into account refl. weighted term. fallspeed at lsmooth=.true. and get vr everywhere:
!!$  lweightdbz = .true., fall = .true., lmds_vr = .true., lfill_vr_backgroundwind = .true.
!!$

!!$ voldata_output_list: possible entries:
!!$
! 'zrsim'       - only if loutdbz=.TRUE.:
!                 simul. radar reflectivity in dBZ  (-999.99=missing value, -99.99=correct 0)
! 'zdrsim'      - only if loutdbz=.TRUE. and loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. radar differential reflectivity in dB  (-999.99=missing value)
! 'rhvsim'      - only if loutdbz=.TRUE. and loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. radar co-polar correlation coefficient in []  (-999.99=missing value)
! 'kdpsim'      - only if loutdbz=.TRUE. and loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. radar specific differential phase in deg/km  (-999.99=missing value)
! 'phidpsim'    - only if loutdbz=.TRUE. and loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. radar differential propagation phase shift in deg  (-999.99=missing value)
! 'vrsim'       - only if loutradwind=.TRUE.:
!                 simul. radial wind in m/s (-999.99=missing value)
! 'zrobs'       - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE.:
!                 obs. radar reflectivity in dBZ
! 'zdrobs'      - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 obs. radar differential reflectivity in dB  (-999.99=missing value)
! 'rhvobs'      - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 obs. radar co-polar correlation coefficient in []  (-999.99=missing value)
! 'kdpobs'      - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                obs. radar specific differential phase in deg/km  (-999.99=missing value)
! 'phidpobs'    - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 obs. radar differential propagation phase shift in deg  (-999.99=missing value)
! 'vrobs'       - only if lreadmeta_from_netcdf=.TRUE. and loutradwind=.TRUE.:
!                 obs. radial wind (dealiasing depends on namelist switch ldealiase_vr_obs!)
! 'qzobs'       - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE.:
!                 quality flags for refl. obs
! 'qzdrobs'     - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 quality flags for diff. refl. obs
! 'qrhvpobs'    - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 quality flags for corr. coeff. obs
! 'qkdpobs'     - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 quality flags for spec. diff. phase obs
! 'qphidpobs'   - only if lreadmeta_from_netcdf=.TRUE. and loutdbz=.TRUE. and loutpolstd=.TRUE.:
!                 quality flags for diff. prop. phase obs
! 'qvobs'       - only if lreadmeta_from_netcdf=.TRUE. and loutradwind=.TRUE.:
!                 quality flags for radial wind obs
! 'losim'       - only if lout_geom=.TRUE.:
!                 simul. geogr. longitude of radar bins degrees
! 'lasim'       - only if lout_geom=.TRUE.:
!                 simul. geogr. latitude of radar bins degrees
! 'hrsim'       - only if lout_geom=.TRUE.:
!                 simul. height of radar bins m MSL
! 'ersim'       - only if lout_geom=.TRUE.:
!                 simul. local beam elevation angle at radar bins degrees
! 'adsim'       - only if lout_geom=.TRUE. and lonline=.TRUE.:
!                 simul. arc distance from radar site (great circle distance)
! 'ahpisim'     - only if lextdbz=.TRUE. and itype_refl=1 or >4:
!                 simul. path integrated attenuation dB
! 'ahsim'       - only if lextdbz=.TRUE. and itype_refl=1 or >4:
!                 simul. specific attenuation (twoway-attenuation coefficient) db/km
! 'adppisim'    - only if lextdbz=.TRUE. and  loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. path integrated differential attenuation in dB
! 'adpsim'      - only if lextdbz=.TRUE. and  loutpolstd=.TRUE. and itype_refl > 4:
!                 simul. specific differential attenuation in dB/km
! 'vrsupsim'    - only if also itype_supobing > 0
! 'vrsupobs'    - only if also itype_supobing > 0
! 'zrsupsim'    - only if also itype_supobing > 0
! 'zrsupobs'    - only if also itype_supobing > 0
! 'zdrsupsim'   - only if also itype_supobing > 0
! 'zdrsupobs'   - only if also itype_supobing > 0
! 'rhvsupsim'   - only if also itype_supobing > 0
! 'rhvsupobs'   - only if also itype_supobing > 0
! 'kdpsupsim'   - only if also itype_supobing > 0
! 'kdpsupobs'   - only if also itype_supobing > 0
! 'losupsim'    - only if also itype_supobing > 0
! 'lasupsim'    - only if also itype_supobing > 0
! 'vasim'       - simul radial wind field with values everywhere, used as a proxy for dealiasing
!                 the observations (in case of ldealiase_vr_obs=.TRUE.)
! 'vrobserr'    - ???
! 'vrsupobserr' - ???
! 'ldrsim'      - only if loutdbz=.TRUE. and loutpolall=.TRUE. and itype_refl > 4:
!                 simul. linear depolarization ratio in []

!!$ renamed lascii_output --> lvoldata_output
!!$ removed lascii_as_fortran_binary from namelist
!!$ added string 'voldata_format' to namelist: 'ascii', 'f90-binary', 'ncdf', 'ascii-gzip'
!!$ added list of output variables 'voldata_output_list' to namelist
!!$ implemented rs_meta%dt_obs_voldata, rs_meta%obs_times_voldata
!!$ implemented dt_obs_voldata_glob, obs_times_voldata_glob
!!$ Implemented elevation choice for voldata along the lines of feedback files
!!$  (ind_ele_voldata_glob (Indiex-list), rs_meta%ind_ele_voldata(Indiex-list) )

!!$ implemented rs_meta%dt_obs_fdbk, rs_meta%obs_times_fdbk
!!$ implemented dt_obs_fdbk_glob, obs_times_fdbk_glob
!!$ Implemented elevation choice for fdbk:
!!$  ind_ele_fdbk_glob (Index-list), rs_meta%ind_ele_fdbk (Indiex-list)

!!$ implemented icountry == 0 for testing of feedbackfiles without having
!!$  actual observations at hand. For this, the DWD radar template is loaded in the
!!$  background, radar station metadata from namelist are taken and the observations
!!$  are faked by copying it from the simulated data.
!!$ You have to explicity specify the obs_times and nobs_times explicitly for each radar
!!$  station in the namelists, although it is unclear why this is necessary.

!!$ New namelist switch "ltestpattern_hydrometeors" for setting an artificial test pattern
!!$  of hydrometeors, useful during development and for testing purposes (testsuite).

!!$ Condition for writing data to fof-files is (obs valid) .AND. (sim valid) , not .OR.

! Meta data of Radar Karlsruhe:
!    rs_meta(1)%station_id          = 10999
!    rs_meta(1)%station_name        = 'KAR'
!    rs_meta(1)%lambda = 0.055
!    rs_meta(1)%lat         = 49.094167
!    rs_meta(1)%lon         =  8.436667
!    rs_meta(1)%alt_agl_mod = 38.0,
!    rs_meta(1)%alt_msl_true  =  148.0,
!    rs_meta(1)%nel         = 14,
!    rs_meta(1)%el_arr      =   0.4, 1.1, 2.0, 3.0, 4.5, 6, 7.5, 9, 11, 13, 16, 20, 24, 30,
!    rs_meta(1)%az_start    = 0.0,
!    rs_meta(1)%naz         = 360,
!    rs_meta(1)%Theta3      = 1.0,
!    rs_meta(1)%Phi3        = 1.0,
!    rs_meta(1)%dalpha      = 1.0,
!    rs_meta(1)%nra         = 240,
!    rs_meta(1)%ra_inc      = 500.0,
!    rs_meta(1)%ext_nyq     = 47.2,
!    rs_meta(1)%prf         = 885.0,
!    rs_meta(1)%dualprf_ratio = 1.33333,
!    rs_meta(1)%dt_obs(1)   = 300.0
!    rs_meta(1)%mds_Z0 = -20.0,         ! Minimum detectable signal [dBZ]
!    rs_meta(1)%mds_r0 = 10000.0,

!================================================================================


  USE radar_kind, ONLY : dp, wpfwo

  USE radar_interface, ONLY : &
       get_model_variables, get_model_hydrometeors, get_model_config_for_radar, &
       get_model_inputdir, get_model_outputdir, set_testpattern_hydrometeors_mg, &
       abort_run,                       &
       get_model_time_sec, get_obstime_ind_of_currtime, get_obstime_ind_of_modtime, it_is_time_for_radar, &
       get_runtime_timings, it_is_time_for_bubblecheck, &
       check_if_currtime_is_obstime, &
       grid_length_model, get_datetime_ini, get_datetime_act, &
       alloc_aux_model_variables, dealloc_aux_model_variables, &
       run_is_restart

#ifdef __ICON__
  USE radar_interface, ONLY :  &
       setup_runtime_timings,  &
#ifdef NUDGING
       set_fdbk_metadata,      &
#endif
       get_domaincenter_global
#endif


  USE radar_interface,     ONLY :                  &
#ifdef __COSMO__
       distribute_composite_cosmo,                 &
       comp_dbzsim,     comp_dbzobs,               &
       comp_dbzsim_bub, comp_dbzobs_bub,           &
#endif
       trigger_warm_bubbles

  USE radar_composites, ONLY : &
       alloc_composite, alloc_composite_bub,       &
       comp_dbzsim_tot,     comp_dbzobs_tot,       &
       comp_dbzsim_bub_tot, comp_dbzobs_bub_tot,   &
       collect_smooth_composite, dealloc_composite, dealloc_composite_bub

  USE radar_bubblegen, ONLY : &
       detect_missing_cells

 !------------------------------------------------------------------------------

   USE radar_utilities, ONLY : get_free_funit, smth_az_horzscan

  !------------------------------------------------------------------------------


  USE radar_parallel_utilities, ONLY :  &
       def_mpi_radar_meta_type,         &
       def_mpi_polmp_type,              &
       def_mpi_dbzcalc_params_type,     &
       def_mpi_compmeta_type,           &
       def_mpi_voldataostream_type

  !------------------------------------------------------------------------------

  USE radar_data,               ONLY :       &
       obsfile_missingname,                  &
       nobstimes_max, nradsta_max, cmaxlen,  &
       switch_to_domain_radar,               &
       rpvect, ipvect,                       &
       mpi_dbzcalc_params_typ,               &
       mpi_polmp_typ, &
       mpi_radar_meta_typ_alltimes,          &
       mpi_radar_meta_typ_onetime,           &
       ydate_ini_mod,                        &
       idom, ndoms_max, ke_fwo,   &
       num_compute_fwo,   & ! number of compute PEs
       num_radar,         & ! number of radar PEs (num_compute + num_radario)
       num_radario,       & ! number of radar-IO PEs
       my_cart_id_fwo,    & ! rank of this PE (=subdomain) in the cartesian communicator
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       my_radario_id,     & ! rank of this PE in the (asynchroneous) radario communicator
       icomm_cart_fwo,    & ! communicator for the virtual cartesian topology
       icomm_radar,       & ! communicator for the group of radar-IO PEs + compute PEs
       icomm_radario,     & ! communicator for the group of radar-IO PEs
       lcompute_pe_fwo,   & ! indicates whether this is a compute PE or not
       lradar_pe,         & ! indicates whether this is a radar PE or not (compute or radar-IO)
       lradario_pe,       & ! indicates whether this is a radar-IO PE or not
       num_radar_dom,     & ! number of radar PEs (num_compute + num_radario_dom) per radar-active model domain
       num_radario_dom,   & ! number of radar-IO PEs per radar-active model domain
       radario_master_dom,& ! root-PEs of the radario group for each active radar domain (in the radar_dom comm., not radario_dom-comm.!!!)
       icomm_radar_dom,   & ! communicator for the group of radar-IO PEs of each domain + compute PEs
       icomm_radario_dom, & ! communicator for the group of radar-IO PEs for each domain
       lradar_pe_dom,     & ! indicates whether this is a radar PE for a certain domain or not (compute or radar-IO)
       lradario_pe_dom,   & ! indicates whether this is a radar-IO PE for a certain domain or not
       my_radar_id_dom,   & ! rank of this PE in the radar communicator (cart+radario_dom)
       my_radario_id_dom, & ! rank of this PE in the (asynchroneous) radario communicator icomm_radario_dom
       i_fwo_prep_compute,& ! Timing flag
       i_fwo_bubbles,     & ! Timing flag
       i_fwo_composites,  & ! Timing flag
       i_fwo_ini,         & ! Timing flag for the initialization of the forward operator
       i_fwo_comm,        & ! Timing flag for MPI-communications
       i_fwo_ongeom,      & ! Timing flag for the ray tracing in online beam propagation
       i_fwo_comppolar,   & ! Timing flag for interpolation of reflectivity and radial wind
                            !  from model grid points to the radar bins/auxiliary azi slice grid
       i_fwo_out,         & ! Timing flag for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                            !  reading obs data, producing feedback files)
       i_fwo_barrier,     & ! Timing flag for barrier waiting in MPI-communications (measure for load imbalance)
       rs_meta, rs_data, rs_grid, cart_data, dbz_meta, i_dbzh, &
       nradsta, itlrad_dyn, itlrad_qx, nbl_az, &
       comp_meta, comp_meta_bub, mpi_comp_meta_typ, mpi_voldata_ostream_typ, &
       i_dwd, i_meteoswiss, i_arpasim

  USE radar_data_namelist, ONLY :  &
       store_domain_radar_nml, switch_to_domain_radar_nml, &
       ldebug_radsim, loutradwind, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       lcheck_inputrecords, lvoldata_output, lfdbk_output, lwrite_ready, ready_file_pattern, &
       ydirradarin, ydirradarout, ysubdircomp, lequal_azi_alldatasets, &
       ydir_mielookup_read, ydir_mielookup_write, ydir_ready_write, icountry, lqc_flag, ldealiase_vr_obs, &
       itype_supobing, supob_nrb, &
       thin_step_azi, thin_step_range, thin_step_ele,  &
       supob_azi_maxsector_vr, supob_cart_resolution, &
       supob_ave_width, supob_minrange_vr, supob_minrange_z, supob_vrw, supob_rfl, &
       lmds_z, lmds_vr, &
       supob_lowthresh_z_obs, supob_lowthresh_z_sim, ltestpattern_hydrometeors, &
       ldo_composite, composite_file_pattern, composite_file_pattern_bub, &
       lfill_vr_backgroundwind, &
       lsmooth_composite_bub_glob,  &
       nsmoothpoints_for_comp_bub_glob, nfilt_for_comp_bub_glob, &
       lsmooth_composite_glob,  &
       nsmoothpoints_for_comp_glob, nfilt_for_comp_glob, nel_composite, levelidlist_for_composite_glob, &
       ldo_bubbles, lcomposite_output, lcomposite_output_bub, &
       dt_bubble_search, prob_bubble, maxdim_obs, lbub_isolated, &
       t_offset_bubble_trigger_async, &
       threshold_obs, threshold_mod, areamin_mod, areamin_obs, &
       mult_dist_obs, add_dist_obs, mult_dist_mod, add_dist_mod, &
       dt_bubble_advect, zlow_meanwind_bubble_advect, zup_meanwind_bubble_advect, &
       tstart_bubble_search, tend_bubble_search, &
       lcomm_nonblocking_online, lcalc_dbz_on_radarbins, &
       labort_if_problems_gribout

  USE radar_obs_meta_list, ONLY : get_meta_proto

  USE radar_namelist_read, ONLY: input_radarnamelist

  !------------------------------------------------------------------------------

  USE radar_model2rays, ONLY : calc_geometry_grid, &
                               calc_geometry_online, &
                               calc_geometry_onlinenew, &
                               calc_geometry_onsmth, &
                               calc_geometry_onsmthnew, &
                               calc_geometry_smth, &
                               calc_geometry_vec, &
                               calc_grd_fallspeed, &
                               calc_grd_reflectivity, &
                               calc_grd_rfridx, &
                               calc_grd_winduvw, &
                               calc_mod_fallspeed, &
                               calc_mod_fallspeed_online, &
                               calc_mod_fallspeed_onsmth, &
                               calc_mod_fallspeed_smth, &
                               calc_mod_radialwind, &
                               calc_mod_radialwind_online, &
                               calc_mod_radialwind_onsmth, &
                               calc_mod_radialwind_smth, &
                               calc_mod_refl_modelgrid, &
                               calc_mod_refl_radarbins, &
                               calc_mod_reflectivity_online, &
                               calc_mod_reflectivity_onsmth, &
                               calc_mod_refl_smth_modelgrid, &
                               calc_mod_refl_smth_radarbins, &
                               calc_sta_rfridx, &
                               distribute_onlineinfos_all_int, &
                               distribute_onlineinfos_all_reals, &
                               distribute_sta_rfridx, &
                               distr_onlinf_all_int_all2allv, &
                               distr_onlinf_all_reals_all2allv

  !------------------------------------------------------------------------------

  USE radar_output_utils,     ONLY : opendiagfiles, closediagfiles
#ifdef GRIBAPI
  USE radar_output_composite, ONLY : write_composite_grib
#endif
  USE radar_output_readyfile, ONLY : write_ready_radar
  
  USE radar_output_driver,  ONLY : init_cart_info, &
                                   output_radar, &
                                   output_radar_smth

  !------------------------------------------------------------------------------

   USE radar_mie_iface_cosmo_driver, ONLY : init_lookup_mie

  !==============================================================================

#ifndef NOMPI
  USE mpi
#endif


!================================================================================
!================================================================================

  IMPLICIT NONE

#ifdef NOMPI
  include "nompi_mpif.h"
#endif

!================================================================================
!================================================================================


#ifdef HAS_MAXRSS
  INTERFACE
    FUNCTION maxrss ()
      INTEGER(kind=KIND(1)) :: maxrss
    END FUNCTION maxrss
  END INTERFACE
#endif

  !==============================================================================

  PRIVATE

  PUBLIC :: any_radar_activity, organize_radar, init_radar

  !==============================================================================

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  SUBROUTINE init_grid_info ( idom_in )

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    INTEGER, INTENT(in)    :: idom_in

    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'emvorado::init_grid_info()'
    INTEGER            :: i, ista, nblaztmp(nradsta_max)


    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! Initialize number of boundlines in auxiliary azimut slice
    ! grid if necessary:
    IF (lsmooth .AND. lonline) THEN
      ! >= 1 line of overlap at the lower and upper azimut index ends for the
      ! auxiliary azimut-slices-grid (must result in more space than needed
      ! by the horizontal averaging stencil)
      ! For simplicity, we take the maximum value needed by any
      ! of the radars.
      ! MIGHT BE CHANGED TO INDIVIDUAL VALUES FOR EACH OF THE RADARS LATERON!
      nblaztmp = 0
      DO ista = 1, nradsta
        nblaztmp(ista) = CEILING( &
             smth_az_horzscan( &
                               rs_meta(ista)%alpha3_eff_0, &
                               rs_meta(ista)%dalpha, &
                               rs_meta(ista)%Phi3, &
                               rs_meta(ista)%smth_interv_fact, &
                               rs_meta(ista)%xabscsm_h(rs_meta(ista)%ngpsm_h), &
                               MAXVAL(rs_meta(ista)%el_arr(1:rs_meta(ista)%nel)) &
                             ) / rs_meta(ista)%az_inc )
        ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.
      END DO
      nbl_az = MAX(MAXVAL(nblaztmp), 1)
!!$ Need modifications for RHI mode!
    ELSE
      nbl_az = 0
    END IF

    IF (ldebug_radsim) THEN
      WRITE (*,*) TRIM(yzroutine)//' : nbl_az = ', nbl_az
    END IF

    IF (lonline) THEN

      DO ista = 1, nradsta

        ! .. The azimut slice grid's horizontal grid distance is
        !    determined as the smaller of the two horizontal model grid
        !    distances:
        rs_grid(ista)%al_inc      = grid_length_model (idom_in)
        rs_grid(ista)%nal         = CEILING(rs_meta(ista)%ra_inc * rs_meta(ista)%nra &
                                            / rs_grid(ista)%al_inc) + 1

        ! needed for index calculation when nbl_az > 0
        ! (extra azimut slices at lower and upper azimut range for lsmooth = .TRUE. .AND. lonline = .true.
        ! in the auxiliary azimuth grid to facilitate interpolation to smoothing points inbetween)
        rs_grid(ista)%naz_nbl     = rs_meta(ista)%naz + 2*nbl_az

        ! Control output:
        IF (ldebug_radsim) THEN
          WRITE (*,'(a,i3,a,f12.2," m")') 'rs_grid(', ista, ')%al_inc  = ',rs_grid(ista)%al_inc
          WRITE (*,'(a,i3,a,i12)')        'rs_grid(', ista, ')%nal     = ',rs_grid(ista)%nal
          WRITE (*,'(a,i3,a,i12)')        'rs_grid(', ista, ')%naz_nbl = ',rs_grid(ista)%naz_nbl
        END IF

      END DO

    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE init_grid_info


  !==============================================================================
  !+ Module procedure in radar_src for determining if there is work to do at the moment:
  !------------------------------------------------------------------------------

  FUNCTION any_radar_activity ( idom_in ) RESULT (there_is_work)

    INTEGER, INTENT(in) :: idom_in
    LOGICAL             :: there_is_work
    INTEGER             :: ista
    REAL(kind=dp)       :: time_mod

    there_is_work = .FALSE.

    IF (lradar_pe) THEN

      time_mod = get_model_time_sec()

      IF (ldebug_radsim .AND. my_radar_id == 0)  WRITE (*,'(a,i3,a,i5,a,f10.0)') &
           'calling any_radar_activity() for domain ',idom_in,' on proc ', my_radar_id, &
           ' for time = ', time_mod


      ! .. switch meta data and data pointers to model domain idom_in and
      !    set internal idom - index of EMVORADO to idom_in (needed in calls to calc_dbz_vec(), calc_fallspeed_vec()):
      CALL switch_to_domain_radar( idom_in )

      ! .. restore the global namelist parameters for this domain from their container storage:
      CALL switch_to_domain_radar_nml ( idom_in )

      IF (loutradwind .OR. loutdbz) THEN

        ! determine if an observation is available for the current time step
        DO ista = 1, nradsta
          there_is_work = (there_is_work .OR. it_is_time_for_radar ( rs_meta(ista)%obs_times ))
        END DO

      END IF

      IF (lradar_pe_dom(idom_in) .AND. num_radario > 0 .AND. lreadmeta_from_netcdf .AND. loutdbz .AND. ldo_bubbles) THEN

        there_is_work = (there_is_work .OR. &
             it_is_time_for_bubblecheck(idom_in, time_mod, &
             &                          tstart_bubble_search, dt_bubble_search, tend_bubble_search, &
             &                          t_offset_bubble_trigger_async))

      END IF

      IF (ldebug_radsim .AND. my_radar_id == 0)  WRITE (*,'(a,i3,a,i5,a,l2,a,f10.0)') &
           'finished any_radar_activity() for domain ',idom_in,' on proc ', my_radar_id, &
           '  activity = ', there_is_work, '  time = ', time_mod

    END IF

  END FUNCTION any_radar_activity

  !==============================================================================
  !+ Module procedure in radar_src for organizing the radar forward operator
  !------------------------------------------------------------------------------

  SUBROUTINE organize_radar (action, itimelevel_dyn, itimelevel_qx, idom_in)

    !------------------------------------------------------------------------------
    !
    ! Description: Organizing routine for the radar forward operator. This routine
    !              is the interface to the radar operator
    !
    ! NOTE: itimelevel_dyn (u,v,w,T,p) and itimelevel_qx (qx,qc,...) have to be consistent with the output
    !       time steps of these fields in the driving model.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    CHARACTER(LEN=*), INTENT(in)  :: action
    INTEGER,          INTENT(in)  :: itimelevel_dyn, itimelevel_qx
    INTEGER,          INTENT(in)  :: idom_in

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER  :: yzroutine = 'emvorado::organize_radar()'
    CHARACTER (LEN=cmaxlen) :: yzerrmsg, zyoutdircomp
    LOGICAL        ::  lcalc(nradsta_max), lcomposite_pe, lobs_composite
    LOGICAL, SAVE  ::  lfirst_cmp   (ndoms_max) = .TRUE.
    LOGICAL, SAVE  ::  lfirst_output(ndoms_max) = .TRUE.
    INTEGER        :: ista, izerror, i, itime
    REAL (KIND=dp) :: time_mod

    TYPE(rpvect) :: pointervec1_real(nradsta_max), pointervec2_real(nradsta_max)
    TYPE(ipvect) :: pointervec1_int(nradsta_max),  pointervec2_int(nradsta_max), pointervec3_int(nradsta_max)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    ! .. switch meta data and data pointers to model domain idom_in and
    !    set internal idom - index of EMVORADO to idom_in (needed in calls to calc_dbz_vec(), calc_fallspeed_vec()):
    CALL switch_to_domain_radar( idom_in )

    ! .. Set the correct time level for the model variables to be used for
    !    the radar simulation: should be nnow for computing steps, consistent with src_output.f90,
    !    and nnew for initialization.
    itlrad_dyn = itimelevel_dyn
    itlrad_qx  = itimelevel_qx

IF (action == 'init') THEN

    ! .. Set some namelist switches to default initial values, because they are used before they
    !     are read from the namelist in init_radar() below:
    ldebug_radsim         = .FALSE.
    lonline               = .FALSE.
    ldo_composite         = .FALSE.
    ldo_bubbles           = .FALSE.
    lreadmeta_from_netcdf = .FALSE.

#ifdef __ICON__
    CALL setup_runtime_timings ()
    CALL get_runtime_timings (i_fwo_ini, lstart=.TRUE.)
#endif

    IF (ldebug_radsim)  WRITE (*,'(a,a,i2,a,i5)') &
         TRIM(yzroutine), ' (initialization stage) for domain ',idom_in,' on proc ', my_radar_id

    ! .. set up model-specific settings like domain sizes, but also internal idom of radar_data.f90
    CALL get_model_config_for_radar ( idom_in )

#ifdef __ICON__
    ! .. determines the center coordinates of each domain and stores into
    !     internal variables in radar_interface.f90 for later use in
    !    CALL get_lonlat_domain_center():
    CALL get_domaincenter_global ( idom_model=idom_in, ldebug=.FALSE. )
#endif

    ! .. Link the needed Model variables:
    IF (lcompute_pe_fwo) THEN
      CALL get_model_hydrometeors (itlrad_qx, idom_in)
      CALL get_model_variables (itlrad_dyn, itlrad_qx, idom_in)

      ! .. allocate model specific data fields for the operator
      CALL alloc_aux_model_variables ( idom_in, .FALSE. )
    END IF

!!$===============================================================
!!$ This has to be called once at the beginning!
!!$===============================================================

    IF (lradar_pe) THEN

      ! .. read namelist for domain idom_in, do cross checks, distribute / initialize all
      !    necessary parameters on all PEs, compute MIE lookup tables if necessary
      !    including former read_meta_info(), read_smth_info() and read_grid_info():
      CALL init_radar ( idom_in )

      ! .. store the global namelist parameters (not those in rs_meta()% or dbz_meta()%)
      !    for the actual domain in a container storage for later re-use, when the operator is again
      !    called for this domain (rs_meta()% or dbz_meta()% have already been linked
      !    to their respective containers above):
      CALL store_domain_radar_nml ( idom_in )

    END IF

#if defined __ICON__ && defined NUDGING
    ! .. Determines some overall metadata for the feedback files
    !     which are the same for all radar stations, distributes
    !     them to all radar PEs and
    !     stores them internally to be retrieved later by
    !    get_fdbk_metadata(idom_model, ldebug) on any of the radar PEs:
    CALL set_fdbk_metadata ( idom_model=idom_in, ldebug=.TRUE. )
#endif

    IF (lcompute_pe_fwo) call dealloc_aux_model_variables ()

#ifdef __COSMO__
    CALL get_runtime_timings (i_fwo_ini)
#endif
#ifdef __ICON__
    CALL get_runtime_timings (-1, lstop=.TRUE.)
#endif


ELSE IF (action == 'compute') THEN   ! full 3D radar simulation

!!$===============================================================
!!$ This has to be called at every timestep!
!!$===============================================================

#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_prep_compute, lstart=.TRUE.)
#endif

    IF (lradar_pe) THEN
      ! .. restore the global namelist parameters for this domain from their container storage:
      CALL switch_to_domain_radar_nml ( idom_in )
    END IF

    lcalc(:) = .FALSE.

    IF (lradar_pe .AND. (loutradwind .OR. loutdbz)) THEN

      ! Time in seconds since model start:
      time_mod = get_model_time_sec ()

      ! determine if an observation is available for the current time step
      DO ista = 1, nradsta
        lcalc(ista) = it_is_time_for_radar ( rs_meta(ista)%obs_times )
      END DO

    END IF

    ! do nothing if nothing has to be done:
    IF (ANY(lcalc)) THEN

      IF (ldebug_radsim) WRITE (*,'(a,a,i2,a,i5)') &
           TRIM(yzroutine), ' for domain ',idom_in,' on proc ', my_radar_id

      ! .. set up model-specific settings like domain sizes, but also internal idom_in of radar_data.f90
      CALL get_model_config_for_radar ( idom_in )

      ! .. Link the needed Model variables:
      IF (lcompute_pe_fwo) THEN
        CALL get_model_hydrometeors (itlrad_qx, idom_in)
        CALL get_model_variables (itlrad_dyn, itlrad_qx, idom_in)
        CALL alloc_aux_model_variables (idom_in, lonline)
        IF (ltestpattern_hydrometeors) THEN
          CALL set_testpattern_hydrometeors_mg ()
        END IF
      END IF

#ifdef HAS_MAXRSS
      WRITE (*,'(a,f10.1,a,i10,a,i5)') &
           'DIAG radar after entry: Memory maximum consumption at t =', time_mod, &
           ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

    END IF ! ANY(lcalc)

    !--------------------------------------------------------------------------------------------
    ! .. Delayed bubble triggering call for asynchroneous IO before the next missing cell detection.
    !    The time delay t_offset_bubble_trigger_async should be smaller or equal to
    !     the bubble search time interval  dt_bubble_search!
    !
    ! THIS SYNCHRONIZES THE WORKER PEs WITH THE ASYNCHR. IO PEs (EXCHANGE OF INFORMATIONS
    !  ON MISSING BUBBLES) AND DESTROYS PART OF THE RUNTIME ADVANTAGE OF ASYNCHR. IO.

    IF (lradar_pe_dom(idom_in) .AND. num_radario > 0 .AND. &
        lreadmeta_from_netcdf .AND. loutdbz .AND. ldo_bubbles) THEN
      IF ( it_is_time_for_bubblecheck(idom_in, time_mod, &
           &                         tstart_bubble_search, dt_bubble_search, tend_bubble_search, &
           &                         t_offset_bubble_trigger_async) ) THEN

        ! .. Find out if there are dbzh observations so that there is valid content in the obs composites:
        lobs_composite = .FALSE.
        DO ista = 1, nradsta
          itime = get_obstime_ind_of_modtime ( time_mod - t_offset_bubble_trigger_async, &
                                               rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) )
          IF (itime > 0) THEN
            IF (rs_meta(ista)%obsfile(itime,i_dbzh) /= obsfile_missingname) THEN
              lobs_composite = .TRUE.
            END IF
          END IF
        END DO

        IF (lobs_composite) THEN
          IF (lcompute_pe_fwo) THEN
            CALL get_model_variables (itlrad_dyn, itlrad_qx, idom_in)
          END IF
          CALL trigger_warm_bubbles (idom_in, time_mod, &
               dt_bubble_advect, zlow_meanwind_bubble_advect, zup_meanwind_bubble_advect)
        END IF
      END IF
    END IF

    !--------------------------------------------------------------------------------------------
    ! .. Main EMVORADO work tree:

    IF (ANY(lcalc)) THEN

      lcomposite_pe = ((num_radario > 0 .AND. lradario_pe_dom(idom_in)) .OR. &
                       (num_radario == 0 .AND. lcompute_pe_fwo))

      IF (ldo_bubbles .AND. lreadmeta_from_netcdf .AND. &
          loutdbz .AND. lcomposite_pe) THEN
        ! allocate comp_dbzobs_bub_tot, comp_dbzsim_bub_tot
        CALL alloc_composite_bub ( ldebug_radsim, comp_meta_bub )
      END IF

      IF (ldo_composite .AND. nel_composite > 0 .AND. &
          loutdbz .AND. lcomposite_pe) THEN
        ! allocate comp_dbzobs_tot, comp_dbzsim_tot
        CALL alloc_composite ( ldebug_radsim, comp_meta, nel_composite, &
                               lreadmeta_from_netcdf )
      END IF

      IF (lonline) THEN

        IF (lfirst_cmp(idom_in)) THEN

#ifdef __ICON__
          CALL get_runtime_timings (i_fwo_ini)
#endif
          IF (lcompute_pe_fwo) THEN

            CALL calc_geometry_grid

#ifdef __COSMO__
            CALL get_runtime_timings ( i_fwo_ini )
#endif
#ifdef __ICON__
            CALL get_runtime_timings (i_fwo_comm)
#endif

            !--------------------- rs_grid%hl_grd -> rs_grid%hl_azgrd ------------
            DO ista = 1,nradsta
              pointervec1_real(ista)%p => rs_grid(ista)%hl_grd
              pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
              pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
              pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
            END DO
            ! should be made for all radars, regardless of lcalc, therefore the ".OR..TRUE." :-)
            IF (lcomm_nonblocking_online) THEN
              CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .OR. .TRUE., &
                   pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            ELSE
              CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .OR. .TRUE., &
                   pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            END IF
            DO ista = 1,nradsta
              IF (ASSOCIATED(rs_grid(ista)%hl_azgrd)) DEALLOCATE(rs_grid(ista)%hl_azgrd)
              rs_grid(ista)%hl_azgrd => pointervec2_real(ista)%p
              NULLIFY(pointervec2_real(ista)%p)
              pointervec2_real(ista)%n = -1
            END DO

            !--------------------- rs_grid%ind_intp(:,2) -> rs_grid%ind_azgrd, rs_grid%nazgrd ----------
            DO ista = 1,nradsta
              pointervec1_int(ista)%p => rs_grid(ista)%ind_intp(:,2)
              pointervec1_int(ista)%n =  rs_grid(ista)%ngrd
              pointervec2_int(ista)%p => rs_grid(ista)%ind_intp(:,2)
              pointervec2_int(ista)%n =  rs_grid(ista)%ngrd
            END DO
            ! should be made for all radars, regardless of lcalc, therefore the ".OR..TRUE." :-)
            IF (lcomm_nonblocking_online) THEN
              CALL distribute_onlineinfos_all_int( lcalc(1:nradsta) .OR. .TRUE., &
                   pointervec1_int(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec3_int(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            ELSE
              CALL distr_onlinf_all_int_all2allv( lcalc(1:nradsta) .OR. .TRUE., &
                   pointervec1_int(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec3_int(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            END IF
            DO ista = 1,nradsta
              IF (ASSOCIATED(rs_grid(ista)%ind_azgrd)) DEALLOCATE(rs_grid(ista)%ind_azgrd)
              rs_grid(ista)%ind_azgrd => pointervec3_int(ista)%p
              rs_grid(ista)%nazgrd    =  pointervec3_int(ista)%n
              NULLIFY(pointervec3_int(ista)%p)
              pointervec3_int(ista)%n = -1
            END DO

          END IF   ! lcompute_pe_fwo

          lfirst_cmp(idom_in) = .FALSE.

#ifdef __COSMO__
          CALL get_runtime_timings (i_fwo_comm)
#endif

        END IF ! lfirst_cmp(idom_in)

#ifdef HAS_MAXRSS
        WRITE (*,'(a,f10.1,a,i10,a,i5)') &
             'DIAG radar after azgrid geometry: Memory maximum consumption at t =', time_mod, &
             ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_comppolar)
#endif
        IF (lcompute_pe_fwo) THEN

          IF (loutradwind .OR. loutdbz) THEN

            CALL calc_sta_rfridx(lcalc)

            CALL calc_grd_rfridx(lcalc)

            IF (loutradwind) THEN
              IF (lfall) THEN
                CALL calc_grd_fallspeed(lcalc)
              END IF
              CALL calc_grd_winduvw(lcalc)        ! defines rs_data(:)%radwind_grd(:)
            END IF

            IF (loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN
              CALL calc_grd_reflectivity(lcalc)      ! defines rs_data(:)%zh_radar_grd(:), rs_data(:)%z_ext_grd(:)
            END IF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after interpolation to azgrid: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

#ifdef __COSMO__
            CALL get_runtime_timings (i_fwo_comppolar)
#endif
#ifdef __ICON__
            CALL get_runtime_timings (i_fwo_comm)
#endif

            IF (loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN

              !--------------------- rs_grid%zh_radar_grd -> rs_grid%zh_radar_azgrd ----------
              DO ista = 1,nradsta
                pointervec1_real(ista)%p => rs_grid(ista)%zh_radar_grd
                pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
              END DO
              IF (lcomm_nonblocking_online) THEN
                CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              ELSE
                CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              END IF
              DO ista = 1,nradsta
                IF (ASSOCIATED(rs_grid(ista)%zh_radar_azgrd)) DEALLOCATE(rs_grid(ista)%zh_radar_azgrd)
                rs_grid(ista)%zh_radar_azgrd => pointervec2_real(ista)%p
                NULLIFY(pointervec2_real(ista)%p)
                pointervec2_real(ista)%n = -1
              END DO

              IF (loutpolstd .OR. loutpolall) THEN

                !--------------------- rs_grid%zv_radar_grd -> rs_grid%zv_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%zv_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%zv_radar_azgrd)) DEALLOCATE(rs_grid(ista)%zv_radar_azgrd)
                  rs_grid(ista)%zv_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

                !--------------------- rs_grid%rrhv_radar_grd -> rs_grid%rrhv_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%rrhv_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%rrhv_radar_azgrd)) DEALLOCATE(rs_grid(ista)%rrhv_radar_azgrd)
                  rs_grid(ista)%rrhv_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

                !--------------------- rs_grid%irhv_radar_grd -> rs_grid%irhv_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%irhv_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%irhv_radar_azgrd)) DEALLOCATE(rs_grid(ista)%irhv_radar_azgrd)
                  rs_grid(ista)%irhv_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

                !--------------------- rs_grid%kdp_radar_grd -> rs_grid%kdp_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%kdp_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%kdp_radar_azgrd)) DEALLOCATE(rs_grid(ista)%kdp_radar_azgrd)
                  rs_grid(ista)%kdp_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

              END IF

              IF (loutpolall) THEN

                !--------------------- rs_grid%zvh_radar_grd -> rs_grid%zvh_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%zvh_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%zvh_radar_azgrd)) DEALLOCATE(rs_grid(ista)%zvh_radar_azgrd)
                  rs_grid(ista)%zvh_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

              END IF

              IF (lextdbz) THEN

                !--------------------- rs_grid%ah_radar_grd -> rs_grid%ah_radar_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%ah_radar_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl==1 .OR. dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                       (dbz_meta(1:nradsta)%itype_refl==1 .OR. dbz_meta(1:nradsta)%itype_refl>4), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%ah_radar_azgrd)) DEALLOCATE(rs_grid(ista)%ah_radar_azgrd)
                  rs_grid(ista)%ah_radar_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

                IF (loutpolstd .OR. loutpolall) THEN

                  !--------------------- rs_grid%adp_radar_grd -> rs_grid%adp_radar_azgrd ----------
                  DO ista = 1,nradsta
                    pointervec1_real(ista)%p => rs_grid(ista)%adp_radar_grd
                    pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                    pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                    pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                  END DO
                  IF (lcomm_nonblocking_online) THEN
                    CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta) .AND. &
                         (dbz_meta(1:nradsta)%itype_refl>4), &
                         pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                         ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                         icomm_cart_fwo,num_compute_fwo)
                  ELSE
                    CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta) .AND. &
                         (dbz_meta(1:nradsta)%itype_refl>4), &
                         pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                         ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                         icomm_cart_fwo,num_compute_fwo)
                  END IF
                  DO ista = 1,nradsta
                    IF (ASSOCIATED(rs_grid(ista)%adp_radar_azgrd)) DEALLOCATE(rs_grid(ista)%adp_radar_azgrd)
                    rs_grid(ista)%adp_radar_azgrd => pointervec2_real(ista)%p
                    NULLIFY(pointervec2_real(ista)%p)
                    pointervec2_real(ista)%n = -1
                  END DO

                END IF

              END IF

            END IF   ! loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))

            IF (loutradwind) THEN

              !--------------------- rs_grid%u_grd -> rs_grid%u_azgrd ----------
              DO ista = 1,nradsta
                pointervec1_real(ista)%p => rs_grid(ista)%u_grd
                pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
              END DO
              IF (lcomm_nonblocking_online) THEN
                CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              ELSE
                CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              END IF
              DO ista = 1,nradsta
                IF (ASSOCIATED(rs_grid(ista)%u_azgrd)) DEALLOCATE(rs_grid(ista)%u_azgrd)
                rs_grid(ista)%u_azgrd => pointervec2_real(ista)%p
                NULLIFY(pointervec2_real(ista)%p)
                pointervec2_real(ista)%n = -1
              END DO

              !--------------------- rs_grid%v_grd -> rs_grid%v_azgrd ----------
              DO ista = 1,nradsta
                pointervec1_real(ista)%p => rs_grid(ista)%v_grd
                pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
              END DO
              IF (lcomm_nonblocking_online) THEN
                CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              ELSE
                CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              END IF
              DO ista = 1,nradsta
                IF (ASSOCIATED(rs_grid(ista)%v_azgrd)) DEALLOCATE(rs_grid(ista)%v_azgrd)
                rs_grid(ista)%v_azgrd => pointervec2_real(ista)%p
                NULLIFY(pointervec2_real(ista)%p)
                pointervec2_real(ista)%n = -1
              END DO

              !--------------------- rs_grid%w_grd -> rs_grid%w_azgrd ----------
              DO ista = 1,nradsta
                pointervec1_real(ista)%p => rs_grid(ista)%w_grd
                pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
              END DO
              IF (lcomm_nonblocking_online) THEN
                CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              ELSE
                CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                     pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                     ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                     icomm_cart_fwo,num_compute_fwo)
              END IF
              DO ista = 1,nradsta
                IF (ASSOCIATED(rs_grid(ista)%w_azgrd)) DEALLOCATE(rs_grid(ista)%w_azgrd)
                rs_grid(ista)%w_azgrd => pointervec2_real(ista)%p
                NULLIFY(pointervec2_real(ista)%p)
                pointervec2_real(ista)%n = -1
              END DO

              IF(lfall) THEN

                !--------------------- rs_grid%vt_grd -> rs_grid%vt_azgrd ----------
                DO ista = 1,nradsta
                  pointervec1_real(ista)%p => rs_grid(ista)%vt_grd
                  pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
                  pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
                  pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
                END DO
                IF (lcomm_nonblocking_online) THEN
                  CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                ELSE
                  CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                       pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                       ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                       icomm_cart_fwo,num_compute_fwo)
                END IF
                DO ista = 1,nradsta
                  IF (ASSOCIATED(rs_grid(ista)%vt_azgrd)) DEALLOCATE(rs_grid(ista)%vt_azgrd)
                  rs_grid(ista)%vt_azgrd => pointervec2_real(ista)%p
                  NULLIFY(pointervec2_real(ista)%p)
                  pointervec2_real(ista)%n = -1
                END DO

              ENDIF

            END IF

            !--------------------- rs_grid%rfridx_grd -> rs_grid%rfridx_azgrd ----------
            DO ista = 1,nradsta
              pointervec1_real(ista)%p => rs_grid(ista)%rfridx_grd
              pointervec1_real(ista)%n =  rs_grid(ista)%ngrd
              pointervec2_int(ista)%p  => rs_grid(ista)%ind_intp(:,2)
              pointervec2_int(ista)%n  =  rs_grid(ista)%ngrd
            END DO
            IF (lcomm_nonblocking_online) THEN
              CALL distribute_onlineinfos_all_reals( lcalc(1:nradsta), &
                   pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            ELSE
              CALL distr_onlinf_all_reals_all2allv( lcalc(1:nradsta), &
                   pointervec1_real(1:nradsta),pointervec2_int(1:nradsta),nradsta, &
                   ones(1:nradsta),rs_meta(1:nradsta)%naz,rs_grid(1:nradsta)%naz_nbl,pointervec2_real(1:nradsta), &
                   icomm_cart_fwo,num_compute_fwo)
            END IF
            DO ista = 1,nradsta
              IF (ASSOCIATED(rs_grid(ista)%rfridx_azgrd)) DEALLOCATE(rs_grid(ista)%rfridx_azgrd)
              rs_grid(ista)%rfridx_azgrd => pointervec2_real(ista)%p
              NULLIFY(pointervec2_real(ista)%p)
              pointervec2_real(ista)%n = -1
            END DO

            CALL distribute_sta_rfridx

#ifdef __COSMO__
            CALL get_runtime_timings (i_fwo_comm)
#endif

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after distribute_onlineinfos: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

#ifdef __ICON__
            CALL get_runtime_timings (i_fwo_ongeom)
#endif

            IF (lsmooth) THEN
              IF (lsode) THEN
                CALL calc_geometry_onsmthnew(lcalc)
              ELSE
                CALL calc_geometry_onsmth(lcalc)
              ENDIF
            ELSE ! lsmooth = .false.
              IF (lsode) THEN
                CALL calc_geometry_onlinenew(lcalc)
              ELSE
                CALL calc_geometry_online(lcalc)
              ENDIF
            END IF ! lsmooth

#ifdef __COSMO__
            CALL get_runtime_timings (i_fwo_ongeom)
#endif

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after online beam propagation: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

#ifdef __ICON__
            CALL get_runtime_timings (i_fwo_comppolar)
#endif
            IF (loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN
              IF (lsmooth) THEN
                CALL calc_mod_reflectivity_onsmth(lcalc)
              ELSE ! lsmooth = .false.
                CALL calc_mod_reflectivity_online(lcalc)
              END IF ! lsmooth
            ENDIF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after online reflectivity: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

            IF (loutradwind) THEN
              IF (lsmooth) THEN
                IF (lfall) THEN
                  CALL calc_mod_fallspeed_onsmth(lcalc)
                ENDIF
                CALL calc_mod_radialwind_onsmth(lcalc)
              ELSE ! lsmooth = .false.
                IF (lfall) THEN
                  CALL calc_mod_fallspeed_online(lcalc)
                ENDIF
                CALL calc_mod_radialwind_online(lcalc)
              END IF ! lsmooth
            ENDIF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after online radial wind: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

          END IF  ! loutradwind or loutdbz

        ELSE      ! lcompute_pe_fwo

          IF (lradario_pe_dom(idom_in)) THEN

            ! Dummy allocation of pointers on the asynchroneous output PEs:
            DO ista = 1, nradsta
              IF (lsmooth) THEN
                rs_data(ista)%nsmth = 0
                rs_data(ista)%nsmth_above_sfc = 0
                ALLOCATE(rs_data(ista)%vt_mod_smth(0))
                ALLOCATE(rs_data(ista)%radwind_mod_smth(0))
                ALLOCATE(rs_data(ista)%zh_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%ah_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%zv_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%rrhv_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%irhv_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%kdp_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%adp_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%zvh_radar_mod_smth(0))
                ALLOCATE(rs_data(ista)%ind_intp_smth(0,2))
              ELSE
                rs_data(ista)%nobs = 0
                rs_data(ista)%nobs_above_sfc = 0
                ALLOCATE(rs_data(ista)%vt_mod(0))
                ALLOCATE(rs_data(ista)%radwind_mod(0))
                ALLOCATE(rs_data(ista)%zh_radar_mod(0))
                ALLOCATE(rs_data(ista)%ah_radar_mod(0))
                ALLOCATE(rs_data(ista)%zv_radar_mod(0))
                ALLOCATE(rs_data(ista)%rrhv_radar_mod(0))
                ALLOCATE(rs_data(ista)%irhv_radar_mod(0))
                ALLOCATE(rs_data(ista)%kdp_radar_mod(0))
                ALLOCATE(rs_data(ista)%adp_radar_mod(0))
                ALLOCATE(rs_data(ista)%zvh_radar_mod(0))
                ALLOCATE(rs_data(ista)%ind_intp(0,2))
              END IF
              ALLOCATE(rs_data(ista)%hl_loc(0))
              ALLOCATE(rs_data(ista)%el_loc(0))
              ALLOCATE(rs_data(ista)%s_loc(0))
            END DO

          END IF

        END IF    ! lcompute_pe_fwo

#ifdef __COSMO__
        CALL get_runtime_timings (i_fwo_comppolar)
#endif
#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_out)
#endif

        IF (lradar_pe_dom(idom_in)) THEN   ! lcompute_pe_fwo or lradario_pe_dom(idom_in)

          IF  (loutradwind .OR. loutdbz) THEN

            IF (ldebug_radsim) THEN
              IF (lfirst_output(idom_in)) THEN
                IF (run_is_restart()) THEN
                  ! In case of restart runs, the diagfiles might or might not exist.
                  ! If existing, they should be opened and appended, if not, they should be created.
                  ! This is achieved  by status "unknown", position "append".
                  CALL opendiagfiles ( "unknown", "append" )
                ELSE
                  CALL opendiagfiles ( "replace", "rewind" )
                END IF
                lfirst_output(idom_in) = .FALSE.
              ELSE
                CALL opendiagfiles ( "old", "append" )
              END IF
            END IF

            IF (lsmooth) THEN
              CALL output_radar_smth(lcalc,time_mod)
            ELSE ! lsmooth = .false.
              CALL output_radar(lcalc,time_mod)
            END IF ! lsmooth

            IF (ldebug_radsim) THEN
              CALL closediagfiles ()
            END IF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after online output: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

            ! clean up memory:
            DO ista = 1,nradsta
              IF (ASSOCIATED(rs_grid(ista)%vt_grd))              DEALLOCATE(rs_grid(ista)%vt_grd)
              IF (ASSOCIATED(rs_grid(ista)%u_grd))               DEALLOCATE(rs_grid(ista)%u_grd)
              IF (ASSOCIATED(rs_grid(ista)%v_grd))               DEALLOCATE(rs_grid(ista)%v_grd)
              IF (ASSOCIATED(rs_grid(ista)%w_grd))               DEALLOCATE(rs_grid(ista)%w_grd)
              IF (ASSOCIATED(rs_grid(ista)%zh_radar_grd))        DEALLOCATE(rs_grid(ista)%zh_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%ah_radar_grd))        DEALLOCATE(rs_grid(ista)%ah_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%zv_radar_grd))        DEALLOCATE(rs_grid(ista)%zv_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%rrhv_radar_grd))      DEALLOCATE(rs_grid(ista)%rrhv_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%irhv_radar_grd))      DEALLOCATE(rs_grid(ista)%irhv_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%kdp_radar_grd))       DEALLOCATE(rs_grid(ista)%kdp_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%adp_radar_grd))       DEALLOCATE(rs_grid(ista)%adp_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%zvh_radar_grd))       DEALLOCATE(rs_grid(ista)%zvh_radar_grd)
              IF (ASSOCIATED(rs_grid(ista)%rfridx_grd))          DEALLOCATE(rs_grid(ista)%rfridx_grd)
              IF (ASSOCIATED(rs_grid(ista)%vt_azgrd))            DEALLOCATE(rs_grid(ista)%vt_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%u_azgrd))             DEALLOCATE(rs_grid(ista)%u_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%v_azgrd))             DEALLOCATE(rs_grid(ista)%v_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%w_azgrd))             DEALLOCATE(rs_grid(ista)%w_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%zh_radar_azgrd))      DEALLOCATE(rs_grid(ista)%zh_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%ah_radar_azgrd))      DEALLOCATE(rs_grid(ista)%ah_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%zv_radar_azgrd))      DEALLOCATE(rs_grid(ista)%zv_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%rrhv_radar_azgrd))    DEALLOCATE(rs_grid(ista)%rrhv_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%irhv_radar_azgrd))    DEALLOCATE(rs_grid(ista)%irhv_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%kdp_radar_azgrd))     DEALLOCATE(rs_grid(ista)%kdp_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%adp_radar_azgrd))     DEALLOCATE(rs_grid(ista)%adp_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%zvh_radar_azgrd))     DEALLOCATE(rs_grid(ista)%zvh_radar_azgrd)
              IF (ASSOCIATED(rs_grid(ista)%rfridx_azgrd))        DEALLOCATE(rs_grid(ista)%rfridx_azgrd)
              IF (ASSOCIATED(rs_data(ista)%vt_mod))              DEALLOCATE(rs_data(ista)%vt_mod)
              IF (ASSOCIATED(rs_data(ista)%vt_mod_smth))         DEALLOCATE(rs_data(ista)%vt_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%radwind_mod))         DEALLOCATE(rs_data(ista)%radwind_mod)
              IF (ASSOCIATED(rs_data(ista)%radwind_mod_smth))    DEALLOCATE(rs_data(ista)%radwind_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zh_radar_mod))        DEALLOCATE(rs_data(ista)%zh_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%zh_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zh_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%ah_radar_mod))        DEALLOCATE(rs_data(ista)%ah_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%ah_radar_mod_smth))   DEALLOCATE(rs_data(ista)%ah_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zv_radar_mod))        DEALLOCATE(rs_data(ista)%zv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%zv_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod))      DEALLOCATE(rs_data(ista)%rrhv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod))      DEALLOCATE(rs_data(ista)%irhv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%irhv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod))       DEALLOCATE(rs_data(ista)%kdp_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%kdp_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%adp_radar_mod))       DEALLOCATE(rs_data(ista)%adp_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%adp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%adp_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod))       DEALLOCATE(rs_data(ista)%zvh_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod_smth))  DEALLOCATE(rs_data(ista)%zvh_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%w_intp))              DEALLOCATE(rs_data(ista)%w_intp)
              IF (ASSOCIATED(rs_data(ista)%w_intp_smth))         DEALLOCATE(rs_data(ista)%w_intp_smth)
              IF (ASSOCIATED(rs_data(ista)%ind_intp))            DEALLOCATE(rs_data(ista)%ind_intp)
              IF (ASSOCIATED(rs_data(ista)%ind_intp_smth))       DEALLOCATE(rs_data(ista)%ind_intp_smth)
              IF (ASSOCIATED(rs_data(ista)%hl_loc))              DEALLOCATE(rs_data(ista)%hl_loc)
              IF (ASSOCIATED(rs_data(ista)%el_loc))              DEALLOCATE(rs_data(ista)%el_loc)
              IF (ASSOCIATED(rs_data(ista)%s_loc))               DEALLOCATE(rs_data(ista)%s_loc)
              IF (ASSOCIATED(cart_data(ista)%dis))               DEALLOCATE(cart_data(ista)%dis)
              IF (ASSOCIATED(cart_data(ista)%calt))              DEALLOCATE(cart_data(ista)%calt)
              IF (ASSOCIATED(cart_data(ista)%rdata_ind))         DEALLOCATE(cart_data(ista)%rdata_ind)
              IF (ASSOCIATED(cart_data(ista)%cartdata_ind))      DEALLOCATE(cart_data(ista)%cartdata_ind)
              IF (ASSOCIATED(cart_data(ista)%aw))                DEALLOCATE(cart_data(ista)%aw)
            END DO

          ENDIF  ! loutradwind .or. loutdbz

        END IF   ! lradar_pe_dom(idom_in)


      ELSE ! lonline=.false.

        IF (lfirst_cmp(idom_in)) THEN

#ifdef __ICON__
          CALL get_runtime_timings (i_fwo_ini)
#endif

          IF (lcompute_pe_fwo) THEN

            IF (loutradwind .OR. loutdbz) THEN

              IF (lsmooth) THEN
                CALL calc_geometry_smth                  ! defines rs_data(:)%w_intp(:,:)
              ELSE
                CALL calc_geometry_vec                   ! defines rs_data(:)%w_intp(:,:)
              END IF
              !         rs_data(:)%ind_intp(:,:)
              !         rs_data(:)%nobs
            END IF

          END IF

          lfirst_cmp(idom_in) = .FALSE.

#ifdef __COSMO__
          CALL get_runtime_timings (i_fwo_ini)
#endif

        END IF  ! lfirst_cmp(idom_in)

#ifdef HAS_MAXRSS
        IF (ldebug_radsim) WRITE (*,'(a,f10.1,a,i10,a,i5)') &
             'DIAG radar after geometry: Memory maximum consumption at t =', time_mod, &
             ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_comppolar)
#endif
        IF (loutradwind .OR. loutdbz) THEN

          IF (lcompute_pe_fwo) THEN

            IF (loutradwind) THEN
              If (lsmooth) THEN
                 IF (lfall) THEN
                    CALL calc_mod_fallspeed_smth(lcalc)
                 ENDIF
                 CALL calc_mod_radialwind_smth(lcalc)        ! defines rs_data(:)%radwind_mod(:)
              ELSE
                 IF (lfall) THEN
                    CALL calc_mod_fallspeed(lcalc)
                 ENDIF
                 CALL calc_mod_radialwind(lcalc)             ! defines rs_data(:)%radwind_mod(:)
              END IF

            END IF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after radialwind: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

            IF (loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN
              IF (lsmooth) THEN
                IF (lcalc_dbz_on_radarbins) THEN
                  CALL calc_mod_refl_smth_radarbins(lcalc)      ! defines rs_data(:)%zh_radar_mod(:)
                ELSE
                  CALL calc_mod_refl_smth_modelgrid(lcalc)      ! defines rs_data(:)%zh_radar_mod(:)
                END IF
              ELSE
                IF (lcalc_dbz_on_radarbins) THEN
                  CALL calc_mod_refl_radarbins(lcalc)           ! defines rs_data(:)%zh_radar_mod(:)
                ELSE
                  CALL calc_mod_refl_modelgrid(lcalc)           ! defines rs_data(:)%zh_radar_mod(:)
                END IF
              END IF
            END IF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after reflectivity: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

          ELSE      ! lradario_pe

            IF (lradario_pe_dom(idom_in)) THEN

              ! Dummy allocation of data pointers on the asynchroneous output PEs:
              DO ista = 1,nradsta
                IF (lsmooth) THEN
                  rs_data(ista)%nsmth = 0
                  rs_data(ista)%nsmth_above_sfc = 0
                  ALLOCATE(rs_data(ista)%vt_mod_smth(0))
                  ALLOCATE(rs_data(ista)%radwind_mod_smth(0))
                  ALLOCATE(rs_data(ista)%zh_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%ah_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%zv_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%rrhv_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%irhv_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%kdp_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%adp_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%zvh_radar_mod_smth(0))
                  ALLOCATE(rs_data(ista)%ind_intp_smth(0,2))
                ELSE
                  rs_data(ista)%nobs = 0
                  rs_data(ista)%nobs_above_sfc = 0
                  ALLOCATE(rs_data(ista)%vt_mod(0))
                  ALLOCATE(rs_data(ista)%radwind_mod(0))
                  ALLOCATE(rs_data(ista)%zh_radar_mod(0))
                  ALLOCATE(rs_data(ista)%ah_radar_mod(0))
                  ALLOCATE(rs_data(ista)%zv_radar_mod(0))
                  ALLOCATE(rs_data(ista)%rrhv_radar_mod(0))
                  ALLOCATE(rs_data(ista)%irhv_radar_mod(0))
                  ALLOCATE(rs_data(ista)%kdp_radar_mod(0))
                  ALLOCATE(rs_data(ista)%adp_radar_mod(0))
                  ALLOCATE(rs_data(ista)%zvh_radar_mod(0))
                  ALLOCATE(rs_data(ista)%ind_intp(0,2))
                END IF
                IF (.NOT.ASSOCIATED(rs_data(ista)%hl_loc)) ALLOCATE(rs_data(ista)%hl_loc(0))
                IF (.NOT.ASSOCIATED(rs_data(ista)%el_loc)) ALLOCATE(rs_data(ista)%el_loc(0))
                IF (.NOT.ASSOCIATED(rs_data(ista)%s_loc))  ALLOCATE(rs_data(ista)%s_loc(0))
              END DO

            END IF

          END IF    ! lcompute_pe_fwo

#ifdef __COSMO__
          CALL get_runtime_timings (i_fwo_comppolar)
#endif
#ifdef __ICON__
          CALL get_runtime_timings (i_fwo_out)
#endif

          IF (lradar_pe_dom(idom_in)) THEN   ! lcompute_pe_fwo .or. lradario_pe_dom(idom_in)

            IF (ldebug_radsim) THEN
              IF (lfirst_output(idom_in)) THEN
                IF (run_is_restart()) THEN
                  ! In case of restart runs, the diagfiles might or might not exist.
                  ! If existing, they should be opened and appended, if not, they should be created.
                  ! This is achieved  by status "unknown", position "append".
                  CALL opendiagfiles ( "unknown", "append" )
                ELSE
                  CALL opendiagfiles ( "replace", "rewind" )
                END IF
                lfirst_output(idom_in) = .FALSE.
              ELSE
                CALL opendiagfiles ( "old", "append" )
              END IF
            END IF

            IF (lsmooth) THEN
              CALL output_radar_smth(lcalc,time_mod)
            ELSE
              CALL output_radar(lcalc,time_mod)
            END IF

            IF (ldebug_radsim) THEN
              CALL closediagfiles ()
            END IF

#ifdef HAS_MAXRSS
            WRITE (*,'(a,f10.1,a,i10,a,i5)') &
                 'DIAG radar after output: Memory maximum consumption at t =', time_mod, &
                 ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

            ! clean up memory:
            DO ista = 1,nradsta
              IF (ASSOCIATED(rs_data(ista)%vt_mod_smth))         DEALLOCATE(rs_data(ista)%vt_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%vt_mod))              DEALLOCATE(rs_data(ista)%vt_mod)
              IF (ASSOCIATED(rs_data(ista)%radwind_mod_smth))    DEALLOCATE(rs_data(ista)%radwind_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%radwind_mod))         DEALLOCATE(rs_data(ista)%radwind_mod)
              IF (ASSOCIATED(rs_data(ista)%zh_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zh_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zh_radar_mod))        DEALLOCATE(rs_data(ista)%zh_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%ah_radar_mod_smth))   DEALLOCATE(rs_data(ista)%ah_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%ah_radar_mod))        DEALLOCATE(rs_data(ista)%ah_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%zv_radar_mod_smth))   DEALLOCATE(rs_data(ista)%zv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zv_radar_mod))        DEALLOCATE(rs_data(ista)%zv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod))      DEALLOCATE(rs_data(ista)%rrhv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%rrhv_radar_mod_smth)) DEALLOCATE(rs_data(ista)%rrhv_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%irhv_radar_mod))      DEALLOCATE(rs_data(ista)%irhv_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%kdp_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%kdp_radar_mod))       DEALLOCATE(rs_data(ista)%kdp_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%adp_radar_mod_smth))  DEALLOCATE(rs_data(ista)%adp_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%adp_radar_mod))       DEALLOCATE(rs_data(ista)%adp_radar_mod)
              IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod_smth))  DEALLOCATE(rs_data(ista)%zvh_radar_mod_smth)
              IF (ASSOCIATED(rs_data(ista)%zvh_radar_mod))       DEALLOCATE(rs_data(ista)%zvh_radar_mod)
              IF (ASSOCIATED(cart_data(ista)%dis))               DEALLOCATE(cart_data(ista)%dis)
              IF (ASSOCIATED(cart_data(ista)%calt))              DEALLOCATE(cart_data(ista)%calt)
              IF (ASSOCIATED(cart_data(ista)%rdata_ind))         DEALLOCATE(cart_data(ista)%rdata_ind)
              IF (ASSOCIATED(cart_data(ista)%cartdata_ind))      DEALLOCATE(cart_data(ista)%cartdata_ind)
              IF (ASSOCIATED(cart_data(ista)%aw))                DEALLOCATE(cart_data(ista)%aw)
            END DO

          END IF   ! lradar_pe_dom(idom_in)

        END IF  ! loutradwind .OR. loutdbz

      END IF ! lonline

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_bubbles)
#endif

      ! .. Find out if there are dbzh observations so that there is valid content in the obs composites:
      lobs_composite = .FALSE.
      IF (lreadmeta_from_netcdf .AND. loutdbz) THEN
        DO ista = 1, nradsta
          itime = get_obstime_ind_of_currtime ( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) )
          IF (itime > 0) THEN
            IF (rs_meta(ista)%obsfile(itime,i_dbzh) /= obsfile_missingname) THEN
              lobs_composite = .TRUE.
            END IF
          END IF
        END DO
      END IF

      IF (ldo_bubbles .AND. lreadmeta_from_netcdf .AND. loutdbz) THEN
        IF (lcomposite_pe) THEN
          CALL collect_smooth_composite (idom_model=idom_in, compdata_tot=comp_dbzsim_bub_tot, &
                                         comp_meta=comp_meta_bub, &
                                         lsmooth_composite=lsmooth_composite_bub_glob, &
                                         nsmoothpoints_for_comp=nsmoothpoints_for_comp_bub_glob, &
                                         nfilt_for_comp=nfilt_for_comp_bub_glob, &
                                         ldebug=ldebug_radsim)
#ifdef __COSMO__
          CALL distribute_composite_cosmo (idom_model=idom_in, compdata_tot=comp_dbzsim_bub_tot, &
                                           compdata=comp_dbzsim_bub, &
                                           comp_meta=comp_meta_bub)
#endif
          IF (lobs_composite) THEN
            CALL collect_smooth_composite (idom_model=idom_in, compdata_tot=comp_dbzobs_bub_tot, &
                                           comp_meta=comp_meta_bub, &
                                           lsmooth_composite=lsmooth_composite_bub_glob, &
                                           nsmoothpoints_for_comp=nsmoothpoints_for_comp_bub_glob, &
                                           nfilt_for_comp=nfilt_for_comp_bub_glob, &
                                           ldebug=ldebug_radsim)
#ifdef __COSMO__
            CALL distribute_composite_cosmo (idom_model=idom_in, compdata_tot=comp_dbzobs_bub_tot, &
                                             compdata=comp_dbzobs_bub, &
                                             comp_meta=comp_meta_bub)
#endif
          END IF
        END IF

        ! Grib-output of composite:
        IF (my_radar_id_dom(idom_in) == radario_master_dom(idom_in) .AND. lcomposite_output_bub) THEN
#ifdef GRIBAPI
          zyoutdircomp(:) = ' '
          IF (ysubdircomp(1:1) == '/') THEN
            zyoutdircomp = TRIM(ysubdircomp)
          ELSE
            zyoutdircomp = TRIM(ydirradarout)//TRIM(ysubdircomp)
          END IF
          CALL write_composite_grib(outdir=TRIM(zyoutdircomp), file_pattern=TRIM(composite_file_pattern_bub), &
                                    cmp_meta=comp_meta_bub, comptyp='sim_bub', ilevel=0, &
                                    comp2d_tot=comp_dbzsim_bub_tot, error=izerror, errmsg=yzerrmsg)
          IF (izerror /= 0) THEN
            IF (labort_if_problems_gribout) THEN
              CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
            ELSE
              WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//'::write_composite_grib: '//TRIM(yzerrmsg)
            END IF
          END IF
          IF (lobs_composite) THEN
            CALL write_composite_grib(outdir=TRIM(zyoutdircomp), file_pattern=TRIM(composite_file_pattern_bub), &
                                      cmp_meta=comp_meta_bub, comptyp='obs_bub', ilevel=0, &
                                      comp2d_tot=comp_dbzobs_bub_tot, error=izerror, errmsg=yzerrmsg)
            IF (izerror /= 0) THEN
              IF (labort_if_problems_gribout) THEN
                CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
              ELSE
                WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//'::write_composite_grib: '//TRIM(yzerrmsg)
              END IF
            END IF
          END IF
#else
          yzerrmsg(:) = ' '
          yzerrmsg = 'ERROR: composite output in grib2-format chosen (lcomposite_output_bub=.TRUE.)'// &
               ' but code is not compiled with -DGRIBAPI!'
          CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
#endif
        END IF

        ! Find warm bubbles (detect_missing_cells) and communicate them to the compute PEs (trigger_warm_bubbles)
        IF (lradar_pe_dom(idom_in) .AND. lobs_composite) THEN

          IF (.NOT.lcomposite_pe) THEN
            ! dummy allocation only:
            ALLOCATE (comp_dbzobs_bub_tot(1,1), comp_dbzsim_bub_tot(1,1))
          END IF

          IF (my_radar_id_dom(idom_in) == radario_master_dom(idom_in)) THEN
            ! Detect missing cells in the global composite fields and store them
            ! in the "bubble_list"-TYPE on the radario_master_dom PE:
            IF ( it_is_time_for_bubblecheck(idom=idom_in, time_mod=time_mod, &
                 &                          tstart_bubblecheck=tstart_bubble_search, &
                 &                          dt_bubblecheck=dt_bubble_search, &
                 &                          tend_bubblecheck=tend_bubble_search, &
                 &                          t_offset_bubblecheck=0.0_dp) ) THEN
              CALL detect_missing_cells(idom_in, time_mod, comp_dbzobs_bub_tot,comp_dbzsim_bub_tot,comp_meta_bub, &
                                        dt_bubble_search, prob_bubble, maxdim_obs, lbub_isolated, &
                                        threshold_obs, threshold_mod, areamin_mod, areamin_obs, &
                                        mult_dist_obs, add_dist_obs, mult_dist_mod, add_dist_mod, ldebug_radsim )
            END IF
          END IF

          IF (num_radario == 0) THEN
            ! This routine triggers the bubbles present in the "bubble_list" on the radario_master_dom PE and resets
            ! the bubble list afterwards on all PEs:
            IF ( it_is_time_for_bubblecheck(idom=idom_in, time_mod=time_mod, &
                 &                          tstart_bubblecheck=tstart_bubble_search, &
                 &                          dt_bubblecheck=dt_bubble_search, &
                 &                          tend_bubblecheck=tend_bubble_search, &
                 &                          t_offset_bubblecheck=0.0_dp) ) THEN
              CALL trigger_warm_bubbles (idom_in, time_mod, &
                   dt_bubble_advect, zlow_meanwind_bubble_advect, zup_meanwind_bubble_advect)
            END IF
          END IF

          IF (.NOT.lcomposite_pe) THEN
            DEALLOCATE (comp_dbzobs_bub_tot, comp_dbzsim_bub_tot)
          END IF

        END IF
        IF (lcomposite_pe) THEN
          CALL dealloc_composite_bub ()  ! deallocate comp_dbzobs_bub_tot, comp_dbzsim_bub_tot
        END IF
      END IF

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_composites)
#endif

      IF (ldo_composite .AND. nel_composite > 0 .AND. loutdbz) THEN
        IF (lcomposite_pe) THEN
          DO i= 1, nel_composite
            CALL collect_smooth_composite (idom_model=idom_in, compdata_tot=comp_dbzsim_tot(:,:,i), &
                                           comp_meta=comp_meta, &
                                           lsmooth_composite=lsmooth_composite_glob, &
                                           nsmoothpoints_for_comp=nsmoothpoints_for_comp_glob, &
                                           nfilt_for_comp=nfilt_for_comp_glob, &
                                           ldebug=ldebug_radsim)
#ifdef __COSMO__
            CALL distribute_composite_cosmo (idom_model=idom_in, compdata_tot=comp_dbzsim_tot(:,:,i), &
                                             compdata=comp_dbzsim(:,:,i), &
                                             comp_meta=comp_meta)
#endif
            IF (lreadmeta_from_netcdf .AND. lobs_composite) THEN
              CALL collect_smooth_composite (idom_model=idom_in, compdata_tot=comp_dbzobs_tot(:,:,i), &
                                             comp_meta=comp_meta, &
                                             lsmooth_composite=lsmooth_composite_glob, &
                                             nsmoothpoints_for_comp=nsmoothpoints_for_comp_glob, &
                                             nfilt_for_comp=nfilt_for_comp_glob, &
                                             ldebug=ldebug_radsim)
#ifdef __COSMO__
              CALL distribute_composite_cosmo (idom_model=idom_in, compdata_tot=comp_dbzobs_tot(:,:,i), &
                                               compdata=comp_dbzobs(:,:,i), &
                                               comp_meta=comp_meta)
#endif
            END IF

            ! Grib-output of composite:
            IF (my_radar_id_dom(idom_in) == radario_master_dom(idom_in) .AND. lcomposite_output) THEN
#ifdef GRIBAPI
              zyoutdircomp(:) = ' '
              IF (ysubdircomp(1:1) == '/') THEN
                zyoutdircomp = TRIM(ysubdircomp)
              ELSE
                zyoutdircomp = TRIM(ydirradarout)//TRIM(ysubdircomp)
              END IF
              CALL write_composite_grib(outdir=TRIM(zyoutdircomp), file_pattern=TRIM(composite_file_pattern), &
                                        cmp_meta=comp_meta, comptyp='sim', &
                                        ilevel=levelidlist_for_composite_glob(i), &
                                        comp2d_tot=comp_dbzsim_tot(:,:,i), error=izerror, errmsg=yzerrmsg)
              IF (izerror /= 0) THEN
                IF (labort_if_problems_gribout) THEN
                  CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
                ELSE
                  WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//'::write_composite_grib: '//TRIM(yzerrmsg)
                END IF
              END IF
              IF (lreadmeta_from_netcdf .AND. lobs_composite) THEN
                CALL write_composite_grib(outdir=TRIM(zyoutdircomp), file_pattern=TRIM(composite_file_pattern), &
                                          cmp_meta=comp_meta, comptyp='obs', &
                                          ilevel=levelidlist_for_composite_glob(i), &
                                          comp2d_tot=comp_dbzobs_tot(:,:,i), error=izerror, errmsg=yzerrmsg)
                IF (izerror /= 0) THEN
                  IF (labort_if_problems_gribout) THEN
                    CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
                  ELSE
                    WRITE (*,'(a)') 'WARNING '//TRIM(yzroutine)//'::write_composite_grib: '//TRIM(yzerrmsg)
                  END IF
                END IF
              END IF
#else
              yzerrmsg(:) = ' '
              yzerrmsg = 'ERROR: composite output in grib2-format chosen (lcomposite_output=.TRUE.)'// &
                   ' but code is not compiled with -DGRIBAPI!'
              CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
#endif
            END IF

          END DO

          IF (lcomposite_pe) THEN
            CALL dealloc_composite ()  ! deallocate comp_dbzobs_tot, comp_dbzsim_tot
          END IF

        END IF     ! lcomposite_pe
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_out)
#endif

      IF ( (loutradwind .OR. loutdbz) .AND. lradar_pe_dom(idom_in) ) THEN

#ifdef __ICON__
        CALL get_runtime_timings (i_fwo_barrier)
#endif

        IF (num_radar > 1) THEN

          IF (ldebug_radsim) THEN
            WRITE (*,*) 'Reached MPI_BARRIER at end of organize_radar on proc ', &
                 my_radar_id, ' obs_time (model) = ', time_mod
          END IF

          IF (num_radario == 0) THEN

            CALL MPI_BARRIER(icomm_radar, izerror)

            IF (lwrite_ready .AND. my_radar_id_dom(idom_in) == radario_master_dom(idom_in)) THEN
              CALL write_ready_radar ( file_pattern = ready_file_pattern, &
                   &                   path     = ydir_ready_write, &
                   &                   content  = 'Finished emvorado for time '//TRIM(get_datetime_act()),         &
                   &                   ierror   = izerror )
              IF (izerror /= 0) THEN
                yzerrmsg(:) = ' '
                yzerrmsg = 'ERROR creating ready file for EMVORADO, '//&
                     'maybe problem with filename or it already exists!'
                CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
              END IF
            END IF

          ELSE

            IF (lradario_pe_dom(idom_in)) THEN

              CALL MPI_BARRIER(icomm_radario_dom(idom_in), izerror)

              IF (lwrite_ready .AND. my_radar_id_dom(idom_in) == radario_master_dom(idom_in)) THEN
                CALL write_ready_radar ( file_pattern = ready_file_pattern, &
                     &                   path     = ydir_ready_write, &
                     &                   content  = 'Finished emvorado for time '//TRIM(get_datetime_act()),         &
                     &                   ierror   = izerror )
                IF (izerror /= 0) THEN
                  yzerrmsg(:) = ' '
                  yzerrmsg = 'ERROR creating ready file for EMVORADO, '//&
                       'maybe problem with filename or it already exists!'
                  CALL abort_run (my_radar_id, izerror, TRIM(yzerrmsg), 'radar_organize.f90, '//TRIM(yzroutine))
                END IF
              END IF

            END IF

          END IF

        END IF

#ifdef __COSMO__
        CALL get_runtime_timings (i_fwo_barrier)
#endif

      END IF   ! lradar_pe_dom(idom_in)

    END IF ! lcalc

    IF (lcompute_pe_fwo) CALL dealloc_aux_model_variables ()


#ifdef HAS_MAXRSS
    WRITE (*,'(a,f10.1,a,i10,a,i5)') &
         'DIAG radar end: Memory maximum consumption at t =', time_mod, &
         ' s: ', maxrss(), ' MB on PE ', my_radar_id
#endif

    IF (ldebug_radsim) WRITE (*,'(a,i2,a,i5)')  &
         'done with '//TRIM(yzroutine)//' for domain ',idom_in,' on proc ', my_radar_id

#ifdef __ICON__
    ! -999 is only a dummy here. This call stops the last started timer:
    CALL get_runtime_timings (-1, lstop=.TRUE.)
#endif

  ENDIF   ! action == 'compute'

  END SUBROUTINE organize_radar


  !==============================================================================
  !+ Module procedure in radar_src for the initialization of the radar operator
  !------------------------------------------------------------------------------

  SUBROUTINE init_radar (idom_in)

    !------------------------------------------------------------------------------
    !
    ! Description: Reads namelist RADARSIM_PARAMS, gets metadata from NetCDF,
    !              initializes MODULE variables, opens control output files
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    INTEGER, INTENT(in)    :: idom_in

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'emvorado::init_radar()'
    CHARACTER (LEN=80) :: yzerrmsg
    CHARACTER (LEN= 9) :: yinput       ! Namelist INPUT file
    INTEGER            :: ista, ipe_rad

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE init_radar
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! Initial date and time of the forecast run in format "YYYYMMDDhhmmss"
    ydate_ini_mod(:) = '0'
    ydate_ini_mod = get_datetime_ini ()

    IF (num_radar > 1) THEN
      ! .. Define the mpi types for the structure rs_meta:
      CALL def_mpi_radar_meta_type     (nobstimes_max, mpi_radar_meta_typ_alltimes  )
      CALL def_mpi_radar_meta_type     (1            , mpi_radar_meta_typ_onetime   )
      ! .. Define the mpi type for the structure polMP (used in dbz_meta):
      CALL def_mpi_polmp_type          (mpi_polmp_typ)
      ! .. Define the mpi type for the structure dbz_meta:
      CALL def_mpi_dbzcalc_params_type (mpi_polmp_typ, mpi_dbzcalc_params_typ)
      ! .. Define the mpi type for the structure comp_meta:
      CALL def_mpi_compmeta_type       (mpi_comp_meta_typ)
      CALL def_mpi_voldataostream_type (mpi_voldata_ostream_typ)
    END IF

    ! Input the namelist variables for the radar simulator:
    !  and, if lreadmeta_from_netcdf = .true., get radar station metadata
    !  from radar netcdf-files. These metadata take precedence over the corresponding
    !  namelist parameters.
    !  (comprises all parameters from the former read_meta_info() and
    !   read_smth_info() )
    CALL input_radarnamelist (idom_in)

    ! Initialize Mie lookup tables in parallel runs (compute or read from file). This
    !  is achieved by calling calc_dbz_vec_modelgrid() for each of the different
    !  Mie-reflectivity configurations (type dbzcalc_params) on different processors.
    !  (Nothing is done for serial runs, because there the tables are generated
    !  at the first regular call(s) to calc_dbz_vec_modelgrid().)
    IF ((loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) .AND. lcompute_pe_fwo) THEN
      CALL init_lookup_mie (nradsta, dbz_meta(1:nradsta), idom_in, &
           TRIM(ydir_mielookup_read), TRIM(ydir_mielookup_write), ldebug_radsim, TRIM(yzroutine))
    END IF

    ! Initialize the parameters for the azimut slice grid for each radar (if needed):
    !  (this needs the parameter nbl_az from above)
    IF (lcompute_pe_fwo) THEN
      CALL init_grid_info( idom_in )
    END IF

    IF (itype_supobing /= 0) THEN
      CALL init_cart_info( idom_in )
    END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE init_radar



END  MODULE radar_organize
