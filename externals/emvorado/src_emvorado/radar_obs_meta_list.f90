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


MODULE radar_obs_meta_list

!------------------------------------------------------------------------------
!
! Description:
!   This module provides methods to define the default/background meta data of radar stations
!   for the radar forward operator EMVORADO, depending on country and WMO radar
!   station identifier.
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
  
  USE radar_dbzcalc_params_type, ONLY : t_dbzcalc_params

  USE radar_data, ONLY : &
       unused_value, missval_int, &
       obsfile_missingname, fdbkfile_missingname, &
       cmaxlen, radar_meta_type, &
       radar_meta_type_onetime, rsm_init_strings_blanks, &
       my_radar_id,       &
       idom, &
       nel_composite_max, nel_max, nobstimes_max, &
       i_fakeobs, i_dwd, i_meteoswiss, i_arpasim, i_belgium, i_denmark, i_france, &
       i_poland, i_czech, i_netherlands, i_slovakia, &
       nradsta_all,     &
       nradsta_dwd,     &
       nradsta_swiss,   &
       nradsta_france,   &
       nradsta_italy,   &
       nradsta_belgium, &
       nradsta_netherlands, &
       nradsta_denmark, &
       nradsta_poland, &
       nradsta_czech, &
       nradsta_slovakia

  USE radar_data_namelist, ONLY :  ldebug_radsim

  USE radar_interface, ONLY: get_lonlat_domain_center, abort_run

  !================================================================================
  !================================================================================

  IMPLICIT NONE

  PRIVATE

  !==============================================================================
  ! Public Subroutines:

  PUBLIC :: set_scanname, get_meta_proto, get_meta_network_belgium, &
       &    get_meta_network_dwd, get_meta_network_italy, get_meta_network_swiss,         &
       &    get_meta_network_all, get_elarr_precipscan

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  INTERFACE set_scanname
    MODULE PROCEDURE     &
         set_scanname_m, &
         set_scanname_o
  END INTERFACE set_scanname

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for DWD radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_dwd ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_dwd'
    
    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for DWD radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 50.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value /)   ! triplet for increment of obs times in seconds, to construct m%obs_times list (from, to, incr)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 10               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'DWD radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0/)  ! nominal elevations (from new DWD radars...)
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 18               ! if metadata are read from obs radar files
    m%nel_default(2) = 18
    m%nel_default(3) = 10
    m%nel_default(4) = 10
    m%nel_default(5) =  1
    m%nel_default(6) =  1
    m%nel_default(7) =  1
    m%nel_default(8) =  1
    m%nel_default(9) =  3
    m%nel_default(10) =  3
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'DWD radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.7, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, &
         9.5,11.0,13.0,15.0,17.0,19.0,23.0,29.0,37.0/)  ! erroneous nominal elevations from old DWD radars, which was in effect at some ASRs.
    m%el_arr_default(1:m%nel_default(2),2)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, &
         9.5,11.0,13.0,15.0,17.0,19.0,23.0,29.0,37.0/)  ! nominal elevations (from old DWD radars...)
    m%el_arr_default(1:m%nel_default(3),3)  = &
         (/0.8, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0/)  ! erroneous nominal elevations (from new DWD radars...), which was in effect sometimes
    m%el_arr_default(1:m%nel_default(4),4)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0/)  ! nominal elevations (from new DWD radars...)
    m%el_arr_default(1:m%nel_default(5),5)  = (/ 0.4 /)  !  precipitation scan is indicated by 0.4° for some stations
    m%el_arr_default(1:m%nel_default(6),6)  = (/ 0.59/)  !  precipitation scan is indicated by 0.59° for some stations
    m%el_arr_default(1:m%nel_default(7),7)  = (/ 0.8 /)  !  precipitation scan is indicated by 0.8° for some stations
    m%el_arr_default(1:m%nel_default(8),8)  = (/ 1.3 /)  !  precipitation scan is indicated by 1.3° for some stations
    m%el_arr_default(1:m%nel_default(9),9)  = (/ 0.5, 1.5, 3.5 /)  ! reduced elevations for operational radial wind DA
    m%el_arr_default(1:m%nel_default(10),10)  = (/ 0.8, 1.5, 3.5 /)  ! reduced elevations for operational radial wind DA

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 124          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 1000.0_dp    ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 32.5_dp                ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 32.5_dp ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1000.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.    ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.    ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRADH'   ! for ras7 and OPERA; ras11 would have 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR' 
    m%obs_hdf5_varname_rhv     = 'RHOHV' 
    m%obs_hdf5_varname_kdp     = 'KDP'  
    m%obs_hdf5_varname_phidp   = 'PHIDP'
    m%obs_hdf5_varname_ldr     = 'LDR' 
    m%obs_hdf5_varname_cflags  = 'CFLAGS'

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_dwd

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for MeteoSwiss radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_swiss ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)      :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_swiss'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Swiss radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 50.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true     = -9999.99_dp  ! True station height MSL (m), taken from obs data
    m%i_nearest_mod    = -9999            ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod    = -9999            ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)     ! triplet for increment of obs times in seconds, to construct m%obs_times list (from, to, incr)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int            ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int            ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)       = unused_value
    m%nel         = 20               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Swiss radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/-0.2, 0.4, 1.0, 1.6, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, &
            8.5, 9.5,11.0,13.0,16.0,20.0,25.0,30.0,35.0,40.0/)  ! nominal elevations (from Swiss radars...)
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 20               ! if metadata are read from obs radar files
    m%nel_default(2) =  5               ! reduced elevs for OPERA data hub
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Swiss radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
          (/-0.2, 0.4, 1.0, 1.6, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, &
             8.5, 9.5,11.0,13.0,16.0,20.0,25.0,30.0,35.0,40.0/)  ! original nominal elevations from Swiss radars
    m%el_arr_default(1:m%nel_default(2),2)  = &
          (/-0.2, 0.4, 1.0, 1.6, 2.5/)  ! elevations which are distributet to OPERA

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360            ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 492          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 500.0_dp     ! Actually used ra_inc for the radar simulations [m]
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -14.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp  ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp  ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp  ! azimutal averaging interval for calculation of one averaged pulse
                                           ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp  ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                   ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                   ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp  ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max)  = 11.5_dp                ! Nyquist velocity of the first obs time, dummy value, will be read from obsfiles
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 11.5_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)         = 600.0_dp  ! Low PRF
    m%dualprf_ratio(1:nel_max) = 1.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 83.3333333_dp          ! Range gate length
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 12
    m%num_pulses  = 32
    m%obsfile(:,:)      = obsfile_missingname  ! names input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file for vr
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .FALSE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.    ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'  ! For both native and OPERA hdf5
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_swiss

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Italian radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_italy ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)      :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_italy'
    
    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Italian radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
!!$ FIXME: does this really have to be such exact? If other countries in the same band (C-, X-, S-) are present in the same run, it is more efficient to make them the same wavelength!
    m%lambda      = 0.05357143_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)     ! triplet for increment of obs times in seconds, to construct m%obs_times list (from, to, incr)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int            ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value    ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int            ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 11               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Italian radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0, 7.0, 9.5, &
         13.0, 18.0, 25.0/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 11               ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Italian radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0, 7.0, 9.5, &
         13.0, 18.0, 25.0/)  ! nominal elevations

    m%az_start    = 0.45_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 400
    m%naz_ncdf(:,:) = 400            ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 0.9_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 212          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 1000.0_dp    ! Actually used ra_inc for the radar simulations [m]
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 0.9_dp  ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 0.9_dp  ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 0.9_dp  ! azimutal averaging interval for calculation of one averaged pulse
                                           ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp  ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                   ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                   ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp  ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 20.1_dp                ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 60.3_dp ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)         = 1500.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)      = obsfile_missingname  ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file for vr
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.    ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.    ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_italy

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Belgian radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_belgium ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_belgium'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Belgian radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/) ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 9                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Belgian radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.5, 1.2, 2.1, 3.4, 4.8, 6.5, 9.0, 13.0, 25.0/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 9                ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Belgian radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.5, 1.2, 2.1, 3.4, 4.8, 6.5, 9.0, 13.0, 25.0/)  ! nominal elevations

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 480          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 500.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 53.5_dp                ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 53.5_dp ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1000.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.    ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.    ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_belgium

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Danish radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_denmark ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_denmark'
    
    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Denmark radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)      ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 11               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Denmark radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    ! Danish radars have alternating strategies from 5' to 5', and this is the superset, but where we
    !  later treat the sometimes appearing 8.4 elevation as 8.5:
    m%el_arr(1:m%nel)  = &
         (/0.5, 0.7, 1.0, 1.5, 2.4, 4.5, 4.8, 8.5, 10.0, 13.0, 15.0/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 11               ! if metadata are read from Denmark netcdf radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Denmark radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.5, 0.7, 1.0, 1.5, 2.4, 4.5, 4.8, 8.5, 10.0, 13.0, 15.0/)
    
    m%az_start    = 0.0_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 480          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 500.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 47.2_dp                ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 47.2_dp ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 800.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'
    m%obs_hdf5_varname_rhv     = 'RHOHV'
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'
    m%obs_hdf5_varname_ldr     = 'LDR'
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_denmark

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Danish radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_france ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_france'
    
    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for French radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)      ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 11               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'French radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    ! French radars have different scan strategies for each station, so we define one
    !  prototype strategy and refine later for each single station:
    m%el_arr(1:m%nel)  = &
         (/0.8, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 11               ! if metadata are read from French netcdf radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'French radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.8, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5/)
    
    m%az_start    = 0.0_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 267          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 960.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 58.89_dp                ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 58.89_dp ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 800.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRADH'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'
    m%obs_hdf5_varname_rhv     = 'RHOHV'
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'
    m%obs_hdf5_varname_ldr     = 'LDR'
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_france

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Polish radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_poland ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_poland'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Polish radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)  ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 10               ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Polish radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.5, 1.4, 2.4, 3.4, 5.3, 7.7, 10.6, 14.1, 18.5, 23.8/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 10               ! if metadata are read from obs radar files
    m%nel_default(2) = 10               ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Polish radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.5, 1.4, 2.4, 3.4, 5.3, 7.7, 10.6, 14.1, 18.5, 23.8/)  ! nominal elevations
    m%el_arr_default(1:m%nel_default(2),2)  = &
         (/0.5, 1.4, 2.1, 3.4, 5.5, 7.7, 10.6, 14.3, 18.6, 23.8/)  ! nominal elevations

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 250          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 1000.0_dp    ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 7.7_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 7.7_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1000.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .FALSE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.    ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRADH'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_poland

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Czech radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_czech ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_czech'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Czech radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)  ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 9                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Czech radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.1, 0.5, 0.9, 1.3, 1.7, 2.2, 3.2, 4.5, 6.3/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 9                ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Czech radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.1, 0.5, 0.9, 1.3, 1.7, 2.2, 3.2, 4.5, 6.3/)  ! default set 1 for nominal elevations, also used as a whitelist to eliminate all other elevations which are not contained

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 650          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 400.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 1.0_dp       ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 1.0_dp       ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 7.25_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 7.25_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1000.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .FALSE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_czech

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Dutch radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_netherlands ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_netherlands'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Netherlands radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod    = 20.0_dp
    m%alt_msl        = -9999.99_dp  ! will be actually used for all computations;
                                    ! < 9000.0 means that it will be determined from
                                    ! m%alt_agl_mod + model orography height at the
                                    ! radar station coordinates
    m%msl_mod        = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true   = -9999.99_dp  ! True station height MSL (m), taken from obs data
    m%i_nearest_mod  = -9999        ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod  = -9999        ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times     = missval_int  ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)  ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs = missval_int      ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk= missval_int      ! similar for the feedback file output
    m%dt_obs_fdbk(:)    = unused_value     ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int ! similar for the volume data file output
    m%dt_obs_voldata(:) = unused_value     ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 5                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Netherlands radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.3, 0.8, 1.2, 2.0, 2.8/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 5                ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Netherlands radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.3, 0.8, 1.2, 2.0, 2.8/)  ! default set 1 for nominal elevations, also used as a whitelist to eliminate all other elevations which are not contained

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 838          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 223.5_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 0.92_dp      ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 0.92_dp      ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                 !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                 !   factor to the effective 3-dB-oneway beamwidths.
                                 !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 32.0_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 32.0_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1000.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 4.0_dp / 3.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRADH'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_netherlands

  !==============================================================================
  !+ Module procedure in radar_src for defining a prototype meta data structure
  !  for Slovak radars as a basis for later refinement for each station.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto_slovakia ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_slovakia'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Slovak radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)   ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 12                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Slovak radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.0, 0.5, 1.0, 1.5, 2.0, 2.7, 3.4, 4.4, 7.0, 11.4, 18.3, 26.7/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 12               ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'Slovak radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.0, 0.5, 1.0, 1.5, 2.0, 2.7, 3.4, 4.4, 7.0, 11.4, 18.3, 26.7/)  ! default set 1 for nominal elevations, also used as a whitelist to eliminate all other elevations which are not contained

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 960          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 250.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 0.95_dp      ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 0.95_dp      ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 40.0_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 40.0_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 600.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 5.0_dp / 4.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 25.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_slovakia

  ! no own icountry, is used for icountry = i_dwd
  FUNCTION get_meta_proto_kitcband ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_kitcband'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Slovak radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)   ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 14                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'KIT C-band radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.4, 1.0, 2.0, 3.0, 4.5, 6.0, 7.5, 9.0, 11.0, 13.0, 16.0, 20.0, 24.0, 30.0/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 14               ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'KIT C-band radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.4, 1.1, 2.0, 3.0, 4.5, 6.0, 7.5, 9.0, 11.0, 13.0, 16.0, 20.0, 24.0, 30.0/)  ! default set 1 for nominal elevations, also used as a whitelist to eliminate all other elevations which are not contained

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 240          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 500.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 0.88_dp      ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 0.88_dp      ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 47.2413_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 47.2413_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 1180.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 3.0_dp / 2.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 15.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'Vh'
    m%obs_hdf5_varname_dbzh    = 'Zh'
    m%obs_hdf5_varname_zdr     = 'ZDR'
    m%obs_hdf5_varname_rhv     = 'RHOHV'
    m%obs_hdf5_varname_kdp     = 'KDP'
    m%obs_hdf5_varname_phidp   = 'PHIDP'
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    ! Varlist from hdf5-file: "'Zh:1 Zv:1 UZh:1 UZv:1 ZDR:1 Vh:1 Wh:1 PHIDP:1 KDP:1 RHOHV:1 SQIh:1 CCORh:1 SNRh:1'"

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_kitcband

  ! no own icountry, is used for icountry = i_dwd
  FUNCTION get_meta_proto_kitcube ( icountry ) RESULT(m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring
    INTEGER                 :: i

    CHARACTER(len=*), PARAMETER :: yzroutine = 'get_meta_proto_kitcube'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ! .. Initialize all strings in the type with blanks:
    CALL rsm_init_strings_blanks (m)
    
    ! .. Prototype parameters in TYPE(radar_meta_type) for Slovak radar station parameters:
    !
    m%station_name  = "XXX"                    ! Short name of station
    m%station_id    = 999999
    m%ista          = 999                      ! internal station index; will be automatically set after namelist reading
    m%lambda      = 0.055_dp
    m%icountry    = icountry

    CALL get_lonlat_domain_center (idom, m%lon, m%lat)

    ! Height AGL. Will not be used for actual computations,
    ! but takes effect if m%alt_msl < -9000.0, in that it
    ! will be added to the model orography height at the
    ! radar station and written to m%alt_msl.
    m%alt_agl_mod = 20.0_dp
    m%alt_msl     = -9999.99_dp  ! will be actually used for all computations;
                                     ! < 9000.0 means that it will be determined from
                                     ! m%alt_agl_mod + model orography height at the
                                     ! radar station coordinates
    m%msl_mod     = -9999.99_dp  ! oro height MSL at station (m). Dummy value, will be determined automatically
    m%alt_msl_true  = -9999.99_dp   ! True station height MSL (m), taken from obs data
    m%i_nearest_mod = -9999             ! i-Index of nearest neighbour grid point to radar station (also dummy)
    m%j_nearest_mod = -9999             ! j-Index of nearest neighbour grid point to radar station (also dummy)
    m%lobstimes_ovwrt_recalc = .FALSE.  ! Flag to enable re-calculation of obs_times from dt_obs and nobs_times after obs data file reading
    m%nobs_times  = missval_int         ! number of obs times
    m%dt_obs(:)   = (/ 300.0_dp, unused_value, unused_value/)   ! triplet for increment of obs times in seconds,
                                        ! to construct m%obs_times list: either <incr>,<miss>,<miss> or <from>,<to>,<incr>)
    m%obs_times(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start, this is the relevant list
    m%nobs_times_obs  = missval_int     ! number of obs times in obs files
    m%obs_times_obs(1:nobstimes_max) = unused_value ! actual obs times in seconds since model start from the obs files
    m%nobs_times_fdbk  = missval_int        ! similar for the feedback file output
    m%dt_obs_fdbk(:)      = unused_value   ! similar for the feedback file output
    m%obs_times_fdbk(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%nobs_times_voldata  = missval_int         ! similar for the volume data file output
    m%dt_obs_voldata(:)      = unused_value    ! similar for the volume data file output
    m%obs_times_voldata(1:nobstimes_max) = unused_value    ! similar for the feedback file output
    m%el_arr   (1:nel_max)  = unused_value
    m%nel            = 10                ! actually used nel
    IF (m%nel > nel_max) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'KIT Cube radar metadata', &
           ' have more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10073, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr(1:m%nel)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0/)  ! nominal elevations
    CALL set_scanname ( m )  ! Short name for scan strategy

    ! default sets of nominal elevations; come into play if metadata are read from radar files:
    m%nel_default(:) = -9999            ! nel of some default scanstrategies; come into play
    m%nel_default(1) = 10               ! if metadata are read from obs radar files
    IF (ANY(m%nel_default(:) > nel_max)) THEN
      errstring(:) = ' '
      WRITE (errstring,'(a,a,i6)') 'KIT Cube radar metadata default scan strategy list ', &
           ' has more elevations than allowed by nel_max = ', nel_max
      CALL abort_run (my_radar_id, 10074, &
           'ERROR: problem in '//TRIM(yzroutine)//'(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, '//TRIM(yzroutine)//'()')
    END IF
    m%el_arr_default(:,:) = unused_value
    m%el_arr_default(1:m%nel_default(1),1)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 8.0, 12.0, 17.0, 25.0/)  ! default set 1 for nominal elevations, also used as a whitelist to eliminate all other elevations which are not contained

    m%az_start    = 0.5_dp       ! center of "regular" first azimut bin [deg]
    m%naz         = 360
    m%naz_ncdf(:,:) = 360        ! dummy obs file azi dimension, will be read from obs files
    m%az_inc      = 1.0_dp
!!$    m%ra_start    = 1000.0_dp
    m%nra         = 200          ! max. number of range bins occuring in a volume scan, actually used for radar simulations
    m%nra_obs     = m%nra        ! Original max. number of range bins from observation files
    m%ra_inc      = 500.0_dp     ! Actually used ra_inc for the radar simulations [m].
    m%ra_inc_obs  = m%ra_inc     ! Original range increment of the input observation data [m]. May be smaller than the actually used ra_inc in obs-data mode for input coarsening option.
    m%ra_inc_coarse=unused_value ! Approximate range increment for range coarsening [m]. ra_inc will be set to the nearest value which is an integer multiple of ra_inc_obs.
    m%n_aggr_ra_obs = 1          ! Number of range bins to aggregate when obs are read (= NINT(ra_inc/ra_inc_obs)); will be automatically determined after meta data reading
    m%mds_Z0      = -20.0_dp     ! Minimum detectable signal at reference range [dBZ]
    m%mds_r0      = 10000.0_dp   ! Reference range for minimum detectable signal [m]
    m%Theta3      = 0.95_dp      ! vertical 3-dB-oneway beam width (deg)
    m%Phi3        = 0.95_dp      ! horizontal 3-dB-oneway  beam width (deg)
    m%dalpha      = 1.0_dp       ! azimutal averaging interval for calculation of one averaged pulse
                                 ! (averaging over the number of statistically independent pulses)
    m%alpha3_eff_0 = 1.461_dp    ! dummy value for the effective horizontal 3-dB-oneway beam width at elevation = 0.0.
                                 ! Depends on the ratio (dalpha/phi3) - has to be exaclty determined later!
                                 ! The provisional value 1.461 is valid for dalpha/phi3 = 1.0
    m%smth_interv_fact = 1.29_dp ! Factor to determine the azimutal and elevational integration range
                                      !   for the smoothing over the beam function. The ranges are computed by multiplying this
                                      !   factor to the effective 3-dB-oneway beamwidths.
                                      !   ( a value of 1.29 leads to the 90-%-weight-range of the beam function)
    m%ngpsm_v     = 1            ! number of vertical smoothing points for Gauss-Legendre quadrature
    m%ngpsm_h     = 1            ! number of horizontal smoothing points for Gauss-Legendre quadrature
    m%xabscsm_v   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%xabscsm_h   = 0.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_v    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%weigsm_h    = 1.0_dp   ! initial dummy value, will be overwritten later and should be omitted in namelist!
    m%high_nyq(1:nel_max) = 16.0_dp                 ! Nyquist velocity of the first obs time
    m%ext_nyq(1:nel_max,1:nobstimes_max) = 16.0_dp  ! Nyquist velocity incl. techniques like Dual-PRF
    m%prf(1:nel_max)      = 667.0_dp  ! Low PRF (???)
    m%dualprf_ratio(1:nel_max) = 3.0_dp / 2.0_dp  ! Dual PRF Ratio
    m%rngate_len  = 125.0_dp               ! Range gate length of observation input data [m]
    m%obs_cdate(1:nobstimes_max) = "YYYYMMDDHHMMSS"       ! array of observation times as string
    m%obs_startrec(1:nobstimes_max,:) = missval_int              ! array of observation start records
    m%obs_endrec(1:nobstimes_max,:)   = missval_int              ! array of observation end records
    m%nrep_ncdf   = -1
    m%num_gates   = 0
    m%num_pulses  = 0
    m%obsfile(:,:)     = obsfile_missingname   ! names of input files for observations
    m%obsfile_format(:)(:) = ' '               ! format of input files for observations
    m%lobs_avail(:,:)   = .FALSE.              ! Flags for success or failure of data reading from NetCDF
    m%fdbkfile          = fdbkfile_missingname ! name of NetCDF feedback file
    m%lfdbkfile_exist   = .FALSE.              ! Flag for existence of feedback file

    m%ind_ele_present(:) = missval_int
    m%nel_present        = m%nel            ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    m%ind_ele_present(1:m%nel_present) = (/ (i, i=1, m%nel_present) /) ! has to be re-set correctly after namelist reading and after obs files reading for the actual timestep
    
    m%ind_ele_fdbk(:)    = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_fdbk           = 0                ! has to be set correctly after namelist reading
    m%ind_ele_voldata(:) = missval_int      ! has to be set correctly after namelist reading, because it is a vector and there is a ind_ele_fdbk_glob
    m%nel_voldata        = 0                ! has to be set correctly after namelist reading

    m%eleind_for_composite_bub = 1         ! elevation index for the construction of the warm bubble generator composite
    m%eleindlist_for_composite(1:nel_composite_max) = (/ (i, i=1, nel_composite_max) /) ! elevation index list for the other composites (can be more than one!)

    ! .. Determine which fields of this station should be written to feedback files:
    m%lvrad_to_fdbk = .TRUE.   ! if true, write radial wind of this stationt to feedback files
    m%vnyq_min_for_vr_active_fdbk = 15.0  !  if lvrad_to_fdbk = .TRUE., only set VRAD active if the actual v_nyq is larger than this value
    m%ldbzh_to_fdbk = .TRUE.   ! if true, write horizontal reflectivity of this station to feedback files

    ! .. Define the exact HDF5 shortnames of the radar moments, which should be used from the obs data:
    m%obs_hdf5_varname_vrad    = 'VRAD'
    m%obs_hdf5_varname_dbzh    = 'DBZH'
    m%obs_hdf5_varname_zdr     = 'ZDR'     ! up to now dummy
    m%obs_hdf5_varname_rhv     = 'RHOHV'   ! up to now dummy
    m%obs_hdf5_varname_kdp     = 'KDP'     ! up to now dummy
    m%obs_hdf5_varname_phidp   = 'PHIDP'   ! up to now dummy
    m%obs_hdf5_varname_ldr     = 'LDR'     ! up to now dummy
    m%obs_hdf5_varname_cflags  = 'CFLAGS'  ! up to now dummy

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END FUNCTION get_meta_proto_kitcube

  !==============================================================================
  !+ Wrapper procedure in radar_src for initializing the meta data prototypes
  !  for the different countries.
  !------------------------------------------------------------------------------

  FUNCTION get_meta_proto ( icountry ) RESULT (m)
    IMPLICIT NONE
    INTEGER                 :: icountry
    TYPE(radar_meta_type)   :: m
    CHARACTER(len=cmaxlen)  :: errstring


    SELECT CASE (icountry)
    CASE (i_dwd)
      m = get_meta_proto_dwd ( icountry )
    CASE (i_meteoswiss)
      m = get_meta_proto_swiss ( icountry )
    CASE (i_arpasim)
      m = get_meta_proto_italy ( icountry )
    CASE (i_belgium)
      m = get_meta_proto_belgium ( icountry )
    CASE (i_denmark)
      m = get_meta_proto_denmark ( icountry )
    CASE (i_france)
      m = get_meta_proto_france ( icountry )
    CASE (i_poland)
      m = get_meta_proto_poland ( icountry )
    CASE (i_czech)
      m = get_meta_proto_czech ( icountry )
    CASE (i_netherlands)
      m = get_meta_proto_netherlands ( icountry )
    CASE (i_slovakia)
      m = get_meta_proto_slovakia ( icountry )
    CASE (i_fakeobs)
      m = get_meta_proto_dwd ( i_dwd )
      m%icountry = i_fakeobs
    CASE default
      errstring(:) = ' '
      WRITE (errstring,'(a)') 'Wrong icountry! Currently implemented are 1 (DWD), 2 (Meteoswiss)' // &
           ', 3 (Italy), 4 (Belgium), 5 (Denmark), 6 (France)!'
      CALL abort_run (my_radar_id, 12077, &
           'ERROR: problem in get_meta_proto(): '//TRIM(ADJUSTL(errstring)), &
           'radar_obs_meta_list.f90, get_meta_proto()')
    END SELECT

  END FUNCTION get_meta_proto
  
  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Italian (Emilia-Romagna) network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_italy (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_italy)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_italy)

    TYPE(radar_meta_type) :: mproto


    ! Scan strategies (from THOMAS):
    ! These are set individually for each station below.
    ! Note: it is not necessary to overtake really all of the different flavours for a certain station.
    !       If one strategy is a subset of another, it can be left out, because it is allowed that
    !       some elevations are missing on input. Their obs data will be filled by miss_values.
!!$    m%el_arr_default(:,:) = unused_value
!!$    m%el_arr_default(1:m%nel_default(1),1)  = &                         ! Gattatico by DPC
!!$         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0, 7.0, 9.5, 13.0, 18.0/)
!!$    m%el_arr_default(1:m%nel_default(2),2)  = &                         ! Gattatico and  S. P.
!!$         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0/)                               ! Capofiume (good weather
!!$                                                                        ! conditions) by DPC
!!$    m%el_arr_default(1:m%nel_default(3),3)  = &                         ! Bric
!!$         (/-0.1,0.5,1.2,2.0,3.0,4.4,5.8,7.4,10.0,15.0/)
!!$    m%el_arr_default(1:m%nel_default(4),4)  = &                         ! Settepani
!!$         (/-0.3,0.7,2.1, 4.0,6.4,9.7,15.0/)
!!$    m%el_arr_default(1:m%nel_default(5),5)  = &                         ! Il Monte, Lauro, Crocione
!!$         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0,11.0,13.5,16.0/)      ! (rarely)
!!$    m%el_arr_default(1:m%nel_default(6),6)  = &                         ! Serano, Zoufplan, Armidda
!!$         (/-0.2, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0,9.0,11.0,13.5,16.0/)
!!$    m%el_arr_default(1:m%nel_default(7),7)  = &                         ! Crocione
!!$         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.1,11.0,13.5,16.0/)
!!$    m%el_arr_default(1:m%nel_default(8),8)  = &                         ! Pettinascura
!!$         (/-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.1,9.1,11.0,13.5,16.0/)
!!$    m%el_arr_default(1:m%nel_default(9),9)  = &                         ! Pettinascura (rarely)
!!$         (/-0.1, 0.6, 1.5, 2.5, 3.5, 4.6, 5.6, 7.1,9.1,11.0,13.5,16.0/)
!!$    m%el_arr_default(1:m%nel_default(10),10)  = &                       ! Grande
!!$         (/0.8, 1.4, 2.4, 3.4, 4.4, 5.9, 7.9, 9.9, 14.9/)
!!$    m%el_arr_default(1:m%nel_default(11),11)  = &                       ! Sagittaria
!!$         (/1.0, 1.4, 2.4, 3.4, 4.4, 5.9, 7.9, 9.9, 14.9/)
!!$    m%el_arr_default(1:m%nel_default(12),12)  = &                       ! Rasu
!!$         (/0.0, 0.7, 1.5, 2.5, 3.6, 5.0, 6.5, 8.4, 10.6, 13.3, 16.3, 20.0/)
!!$    m%el_arr_default(1:m%nel_default(13),13)  = &                       ! Midia
!!$         (/0.5, 1.5, 2.5, 3.5/)
!!$    m%el_arr_default(1:m%nel_default(14),14)  = &                       ! Macaion 
!!$         (/0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 14.0/)
!!$    m%el_arr_default(1:m%nel_default(15),15)  = &                       ! S. P. Capofiume by Arpae
!!$         (/0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 25.0/)                         ! (good weather conditions)
!!$    m%el_arr_default(1:m%nel_default(16),16)  = &                       ! S. P. Capofiume by Arpae
!!$         (/0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 8.0, 11.0, 18.0, 25.0/)        ! (bad  weather conditions)
!!$    m%el_arr_default(1:m%nel_default(17),17)  = &                       ! S. P. Capofiume by DPC
!!$         (/0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 8.0, 11.0, 18.0/)              ! (bad  weather conditions)


    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_arpasim )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): San Pietro Capofiume
    !
    m(1)%station_name           = 'SPC'
    m(1)%station_id             = 16144
    m(1)%lambda                 = 0.055
    m(1)%lat                    = 44.655
    m(1)%lon                    = 11.624
    m(1)%alt_agl_mod            = 20.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl                = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true           = 11.0
    m(1)%mds_Z0                 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0                 = 10000.0
    m(1)%dt_obs(1)              = 300.0         ! 5 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(1)%nel_default(1)         = 10
    m(1)%el_arr_default(:,1)    = unused_value
    m(1)%el_arr_default(1:m(1)%nel_default(1),1)  = &              ! S. P. Capofiume by Arpae (bad weather conditions)
         (/0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 8.0, 11.0, 18.0, 25.0/)
    m(1)%nel_default(2)         = 7
    m(1)%el_arr_default(:,2)    = unused_value
    m(1)%el_arr_default(1:m(1)%nel_default(2),2)  = &              ! S. P. Capofiume by Arpae (good weather conditions)  
         (/0.5, 1.4, 2.3, 3.2, 4.1, 5.0, 25.0/)                    ! (contains the elevations distributed by DPC)
    !
    m(1)%nel                    = m(1)%nel_default(1)    ! Default elevation set is the first of the default elevations
    m(1)%el_arr(:)              = unused_value
    m(1)%el_arr(1:m(1)%nel)     = m(1)%el_arr_default(1:m(1)%nel_default(1),1)
    !
    dbzparams(1)%station_id     = 16144
    !
    CALL set_scanname ( m(1) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Gattatico
    !
    m(2)%station_name           = 'GAT'
    m(2)%station_id             = 16199
    m(2)%lambda                 = 0.055
    m(2)%lat                    = 44.7914
    m(2)%lon                    = 10.4992
    m(2)%alt_agl_mod            = 20.0
    m(2)%alt_msl                = -9999.99
    m(2)%alt_msl_true = 40.0
    m(2)%mds_Z0                 = -20.0
    m(2)%mds_r0 = 10000.0
    m(2)%dt_obs(1)              = 300.0         ! 5 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(2)%nel_default(1)         = 10
    m(2)%el_arr_default(:,1)    = unused_value
    m(2)%el_arr_default(1:m(2)%nel_default(1),1)  = &             ! Gattatico by DPC
         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0, 7.0, 9.5, 13.0, 18.0/)
    m(2)%nel_default(2)         = 6
    m(2)%el_arr_default(:,2)    = unused_value
    m(2)%el_arr_default(1:m(2)%nel_default(2),2)  = &             ! Gattatico good weather conditions) by DPC
         (/0.5, 1.4, 2.3, 3.2, 4.2, 5.0/)
    !
    m(2)%nel                    = m(2)%nel_default(1)    ! Default elevation set is the first of the default elevations
    m(2)%el_arr(:)              = unused_value
    m(2)%el_arr(1:m(2)%nel)     = m(2)%el_arr_default(1:m(2)%nel_default(1),1)
    !
    dbzparams(2)%station_id     = 16199
    !
    CALL set_scanname ( m(2) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): Settepani
    !
    m(3)%station_name           = 'SET'
    m(3)%station_id             = 16999
    m(3)%lambda                 = 0.055
    m(3)%lat                    = 44.245897
    m(3)%lon                    = 8.197409
    m(3)%alt_agl_mod            = 11.0
    m(3)%alt_msl                = 1400.0
    m(3)%alt_msl_true           = 1400.0
    m(3)%mds_Z0                 = -20.0
    m(3)%mds_r0                 = 10000.0
    m(3)%dt_obs(1)              = 300.0         ! 5 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(3)%az_start               = 0.0
    m(3)%naz                    = 360
    m(3)%az_inc                 = 1.0
    m(3)%nra                    = 172
    m(3)%ra_inc                 = 1000.0
    !
    m(3)%nel_default(1)         = 7
    m(3)%el_arr_default(:,1)    = unused_value
    m(3)%el_arr_default(1:m(3)%nel_default(1),1)  = &
         (/-0.3, 0.7, 2.1, 4.0, 6.4, 9.7, 15.0/)
    !
    m(3)%nel                    = m(3)%nel_default(1)
    m(3)%el_arr(:)              = unused_value
    m(3)%el_arr(1:m(3)%nel)     = m(3)%el_arr_default(1:m(3)%nel_default(1),1)
    !
    dbzparams(3)%station_id     = 16999
    !
    CALL set_scanname ( m(3) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 4: m(4) and dbzparams(4): Bric della Croce
    !
    m(4)%station_name           = 'BRC'
    m(4)%station_id             = 16998
    m(4)%lambda                 = 0.055
    m(4)%lat                    = 45.034
    m(4)%lon                    = 7.733
    m(4)%alt_msl                = 773.0
    m(4)%alt_msl_true           = 773.0
    m(4)%mds_Z0                 = -20.0
    m(4)%mds_r0                 = 10000.0
    m(4)%dt_obs(1)              = 300.0         ! 5 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(4)%az_start               = 0.0
    m(4)%naz                    = 360
    m(4)%az_inc                 = 1.0
    m(4)%nra                    = 170
    m(4)%ra_inc                 = 1000.0
    !
    m(4)%nel_default(1)         = 10
    m(4)%el_arr_default(:,1)    = unused_value
    m(4)%el_arr_default(1:m(4)%nel_default(1),1)  = &
         (/-0.1, 0.5, 1.2, 2.0, 3.0, 4.4, 5.8, 7.4, 10.0, 15.0/)
    !
    m(4)%nel                    = m(4)%nel_default(1)
    m(4)%el_arr(:)              = unused_value
    m(4)%el_arr(1:m(4)%nel)     = m(4)%el_arr_default(1:m(4)%nel_default(1),1)
    !
    dbzparams(4)%station_id     = 16998
    !
    CALL set_scanname ( m(4) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 5: m(5) and dbzparams(5): Armidda
    !
    m(5)%station_name           = 'ARM'
    m(5)%station_id             = 16112
    m(5)%lambda                 = 0.055
    m(5)%lat                    = 39.8822
    m(5)%lon                    = 9.49377
    m(5)%alt_agl_mod            = 12.0
    m(5)%alt_msl                = 1261.0
    m(5)%alt_msl_true           = 1261.0
    m(5)%mds_Z0                 = -20.0
    m(5)%mds_r0                 = 10000.0
    m(5)%dt_obs(1)              = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(5)%az_start               = 0.0
    m(5)%naz                    = 360
    m(5)%az_inc                 = 1.0
    m(5)%nra                    = 200
    m(5)%ra_inc                 = 1000.0
    !
    m(5)%nel_default(1)         = 12
    m(5)%el_arr_default(:,1)    = unused_value
    m(5)%el_arr_default(1:m(5)%nel_default(1),1)  = &
         (/-0.2, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.5, 16.0/)
    !
    m(5)%nel                    = m(5)%nel_default(1)
    m(5)%el_arr(:)              = unused_value
    m(5)%el_arr(1:m(5)%nel)     = m(5)%el_arr_default(1:m(5)%nel_default(1),1)
    !
    dbzparams(5)%station_id     = 16112
    !
    CALL set_scanname ( m(5) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 6: m(6) and dbzparams(6): Crocione
    !
    m(6)%station_name           = 'CRO'
    m(6)%station_id             = 16101
    m(6)%lambda                 = 0.055
    m(6)%lat                    = 43.9615
    m(6)%lon                    = 10.6103
    m(6)%alt_agl_mod            = 12.0
    m(6)%alt_msl                = 1044.0
    m(6)%alt_msl_true           = 1044.0
    m(6)%mds_Z0                 = -20.0
    m(6)%mds_r0                 = 10000.0
    m(6)%dt_obs(1)              = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(6)%az_start               = 0.0
    m(6)%naz                    = 360
    m(6)%az_inc                 = 1.0
    m(6)%nra                    = 200
    m(6)%ra_inc                 = 1000.0
    !
    m(6)%nel_default(1)         = 11
    m(6)%el_arr_default(:,1)    = unused_value
    m(6)%el_arr_default(1:m(6)%nel_default(1),1)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.5, 16.0/)  ! Sometimes 9.1 instead of 9.0, but this is within the tolerance of 0.2 deg
    m(6)%nel_default(2)         = 11
    m(6)%el_arr_default(:,2)    = unused_value
    m(6)%el_arr_default(1:m(6)%nel_default(2),2)  = &   ! rarely
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0 ,15.0, 16.0/)
    !
    m(6)%nel                    = m(6)%nel_default(1)    ! Default elevation set is the first of the default elevations
    m(6)%el_arr(:)              = unused_value
    m(6)%el_arr(1:m(6)%nel)     = m(6)%el_arr_default(1:m(6)%nel_default(1),1)
    !
    dbzparams(6)%station_id     = 16101
    !
    CALL set_scanname ( m(6) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 7: m(7) and dbzparams(7): Il Monte
    !
    m(7)%station_name           = 'MON'
    m(7)%station_id             = 16102
    m(7)%lambda                 = 0.055
    m(7)%lat                    = 41.9394
    m(7)%lon                    = 14.6208
    m(7)%alt_agl_mod            = 12.0
    m(7)%alt_msl                = 710.0
    m(7)%alt_msl_true           = 710.0
    m(7)%mds_Z0                 = -20.0
    m(7)%mds_r0                 = 10000.0
    m(7)%dt_obs(1)              = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(7)%az_start               = 0.0
    m(7)%naz                    = 360
    m(7)%az_inc                 = 1.0
    m(7)%nra                    = 200
    m(7)%ra_inc                 = 1000.0
    !
    m(7)%nel_default(1)         = 11
    m(7)%el_arr_default(:,1)    = unused_value
    m(7)%el_arr_default(1:m(7)%nel_default(1),1)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0 ,13.5, 16.0/)
    m(7)%nel_default(2)         = 11
    m(7)%el_arr_default(:,2)    = unused_value
    m(7)%el_arr_default(1:m(7)%nel_default(2),2)  = &   ! rarely
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0 ,15.0, 16.0/)
    !
    m(7)%nel                    = m(7)%nel_default(1)
    m(7)%el_arr(:)              = unused_value
    m(7)%el_arr(1:m(7)%nel)     = m(7)%el_arr_default(1:m(7)%nel_default(1),1)
    !
    dbzparams(7)%station_id     = 16102
    !
    CALL set_scanname ( m(7) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 8: m(8) and dbzparams(8): Lauro
    !
    m(8)%station_name           = 'LAU'
    m(8)%station_id             = 16107
    m(8)%lambda                 = 0.055
    m(8)%lat                    = 37.1126
    m(8)%lon                    = 14.8357
    m(8)%alt_agl_mod            = 12.0
    m(8)%alt_msl                = 980.0
    m(8)%alt_msl_true           = 980.0
    m(8)%mds_Z0                 = -20.0
    m(8)%mds_r0                 = 10000.0
    m(8)%dt_obs(1)              = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(8)%az_start               = 0.0
    m(8)%naz                    = 360
    m(8)%az_inc                 = 1.0
    m(8)%nra                    = 200
    m(8)%ra_inc                 = 1000.0
    !
    m(8)%nel_default(1)         = 11
    m(8)%el_arr_default(:,1)    = unused_value
    m(8)%el_arr_default(1:m(8)%nel_default(1),1)  = &
         (/0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0 ,13.5, 16.0/)
    !
    m(8)%nel                    = m(8)%nel_default(1)
    m(8)%el_arr(:)              = unused_value
    m(8)%el_arr(1:m(8)%nel)     = m(8)%el_arr_default(1:m(8)%nel_default(1),1)
    !
    dbzparams(8)%station_id     = 16107
    !
    CALL set_scanname ( m(8) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 9: m(9) and dbzparams(9): Pettinascura
    !
    m(9)%station_name           = 'PET'
    m(9)%station_id             = 16105
    m(9)%lambda                 = 0.055
    m(9)%lat                    = 39.36981
    m(9)%lon                    = 16.61833
    m(9)%alt_agl_mod            = 12.0
    m(9)%alt_msl                = 1725.0
    m(9)%alt_msl_true           = 1725.0
    m(9)%mds_Z0                 = -20.0
    m(9)%mds_r0                 = 10000.0
    m(9)%dt_obs(1)              = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(9)%az_start               = 0.0
    m(9)%naz                    = 360
    m(9)%az_inc                 = 1.0
    m(9)%nra                    = 200
    m(9)%ra_inc                 = 1000.0
    !
    m(9)%nel_default(1)         = 12
    m(9)%el_arr_default(:,1)    = unused_value
    m(9)%el_arr_default(1:m(9)%nel_default(1),1)  = &
         (/-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.1, 9.1, 11.0, 13.5, 16.0/)
    m(9)%nel_default(2)         = 12
    m(9)%el_arr_default(:,2)    = unused_value
    m(9)%el_arr_default(1:m(9)%nel_default(2),2)  = &                     ! (rarely)
         (/-0.1, 0.6, 1.5, 2.5, 3.5, 4.6, 5.6, 7.1, 9.1, 11.0, 13.5, 16.0/)
    !
    m(9)%nel                    = m(9)%nel_default(1)
    m(9)%el_arr(:)              = unused_value
    m(9)%el_arr(1:m(9)%nel)     = m(9)%el_arr_default(1:m(9)%nel_default(1),1)
    !
    dbzparams(9)%station_id     = 16105
    !
    CALL set_scanname ( m(9) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 10: m(10) and dbzparams(10): Serano
    !
    m(10)%station_name          = 'SER'
    m(10)%station_id            = 16103
    m(10)%lambda                = 0.055
    m(10)%lat                   = 42.8659
    m(10)%lon                   = 12.8002
    m(10)%alt_agl_mod           = 12.0
    m(10)%alt_msl               = 1400.0
    m(10)%alt_msl_true          = 1400.0
    m(10)%mds_Z0                = -20.0
    m(10)%mds_r0                = 10000.0
    m(10)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(10)%az_start              = 0.0
    m(10)%naz                   = 360
    m(10)%az_inc                = 1.0
    m(10)%nra                   = 200
    m(10)%ra_inc                = 1000.0
    !
    m(10)%nel_default(1)        = 12
    m(10)%el_arr_default(:,1)   = unused_value
    m(10)%el_arr_default(1:m(10)%nel_default(1),1)  = &
         (/-0.2, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.5, 16.0/)
    !
    m(10)%nel                   = m(10)%nel_default(1)
    m(10)%el_arr(:)             = unused_value
    m(10)%el_arr(1:m(10)%nel)   = m(10)%el_arr_default(1:m(10)%nel_default(1),1)
    !
    dbzparams(10)%station_id    = 16103
    !
    CALL set_scanname ( m(10) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 11: m(11) and dbzparams(11): Zoufplan
    !
    m(11)%station_name          = 'ZOU'
    m(11)%station_id            = 16106
    m(11)%lambda                = 0.055
    m(11)%lat                   = 46.5625
    m(11)%lon                   = 12.9703
    m(11)%alt_agl_mod           = 12.0
    m(11)%alt_msl               = 2001.0
    m(11)%alt_msl_true          = 2001.0
    m(11)%mds_Z0                = -20.0
    m(11)%mds_r0                = 10000.0
    m(11)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(11)%az_start              = 0.0
    m(11)%naz                   = 360
    m(11)%az_inc                = 1.0
    m(11)%nra                   = 200
    m(11)%ra_inc                = 1000.0
    !
    m(11)%nel_default(1)        = 12
    m(11)%el_arr_default(:,1)   = unused_value
    m(11)%el_arr_default(1:m(11)%nel_default(1),1)  = &
         (/-0.2, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.5, 16.0/)
    !
    m(11)%nel                   = m(11)%nel_default(1)
    m(11)%el_arr(:)             = unused_value
    m(11)%el_arr(1:m(11)%nel)   = m(11)%el_arr_default(1:m(11)%nel_default(1),1)
    !
    dbzparams(11)%station_id    = 16106
    !
    CALL set_scanname ( m(11) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 12: m(12) and dbzparams(12): Grande
    !
    m(12)%station_name          = 'GRA'
    m(12)%station_id            = 16903
    m(12)%lambda                = 0.055
    m(12)%lat                   = 45.3619
    m(12)%lon                   = 11.6728
    m(12)%alt_agl_mod           = 12.0
    m(12)%alt_msl               = 480.0
    m(12)%alt_msl_true          = 480.0
    m(12)%mds_Z0                = -20.0
    m(12)%mds_r0                = 10000.0
    m(12)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(12)%az_start              = 0.0
    m(12)%naz                   = 360
    m(12)%az_inc                = 1.0
    m(12)%nra                   = 120
    m(12)%ra_inc                = 1000.0
    !
    m(12)%nel_default(1)        = 9
    m(12)%el_arr_default(:,1)   = unused_value
    m(12)%el_arr_default(1:m(12)%nel_default(1),1)  = &
         (/0.8, 1.4, 2.4, 3.4, 4.4, 5.9, 7.9, 9.9, 14.9/)
    !
    m(12)%nel                   = m(12)%nel_default(1)
    m(12)%el_arr(:)             = unused_value
    m(12)%el_arr(1:m(12)%nel)   = m(12)%el_arr_default(1:m(12)%nel_default(1),1)
    !
    dbzparams(12)%station_id    = 16903
    !
    CALL set_scanname ( m(12) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 13: m(13) and dbzparams(13): Sagittaria
    !
    m(13)%station_name          = 'SAG'
    m(13)%station_id            = 16905
    m(13)%lambda                = 0.055
    m(13)%lat                   = 45.6942
    m(13)%lon                   = 12.7859
    m(13)%alt_agl_mod           = 12.0
    m(13)%alt_msl               = 15.0
    m(13)%alt_msl_true          = 15.0
    m(13)%mds_Z0                = -20.0
    m(13)%mds_r0                = 10000.0
    m(13)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(13)%az_start              = 0.0
    m(13)%naz                   = 360
    m(13)%az_inc                = 1.0
    m(13)%nra                   = 120
    m(13)%ra_inc                = 1000.0
    !
    m(13)%nel_default(1)        = 9
    m(13)%el_arr_default(:,1)   = unused_value
    m(13)%el_arr_default(1:m(13)%nel_default(1),1)  = &
         (/1.0, 1.4, 2.4, 3.4, 4.4, 5.9, 7.9, 9.9, 14.9/)
    !
    m(13)%nel                   = m(13)%nel_default(1)
    m(13)%el_arr(:)             = unused_value
    m(13)%el_arr(1:m(13)%nel)   = m(13)%el_arr_default(1:m(13)%nel_default(1),1)
    !
    dbzparams(13)%station_id    = 16905
    !
    CALL set_scanname ( m(13) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 14: m(14) and dbzparams(14): Rasu
    !
    m(14)%station_name          = 'RAS'
    m(14)%station_id            = 16907
    m(14)%lambda                = 0.055
    m(14)%lat                   = 40.4217
    m(14)%lon                   = 9.00463
    m(14)%alt_agl_mod           = 12.0
    m(14)%alt_msl               = 1250.0
    m(14)%alt_msl_true          = 1250.0
    m(14)%mds_Z0                = -20.0
    m(14)%mds_r0                = 10000.0
    m(14)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(14)%az_start              = 0.0
    m(14)%naz                   = 360
    m(14)%az_inc                = 1.0
    m(14)%nra                   = 200
    m(14)%ra_inc                = 1000.0
    !
    m(14)%nel_default(1)        = 12
    m(14)%el_arr_default(:,1)   = unused_value
    m(14)%el_arr_default(1:m(14)%nel_default(1),1)  = &
         (/0.0, 0.7, 1.5, 2.5, 3.6, 5.0, 6.5, 8.4, 10.6, 13.3, 16.3, 20.0/)
    !
    m(14)%nel                   = m(14)%nel_default(1)
    m(14)%el_arr(:)             = unused_value
    m(14)%el_arr(1:m(14)%nel)   = m(14)%el_arr_default(1:m(14)%nel_default(1),1)
    !
    dbzparams(14)%station_id    = 16907
    !
    CALL set_scanname ( m(14) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 15: m(15) and dbzparams(15): Midia
    !
    m(15)%station_name          = 'MID'
    m(15)%station_id            = 16996
    m(15)%lambda                = 0.055
    m(15)%lat                   = 42.0578
    m(15)%lon                   = 13.1772
    m(15)%alt_agl_mod           = 12.0
    m(15)%alt_msl               = 1727.0
    m(15)%alt_msl_true          = 1727.0
    m(15)%mds_Z0                = -20.0
    m(15)%mds_r0                = 10000.0
    m(15)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(15)%az_start              = 0.0
    m(15)%naz                   = 360
    m(15)%az_inc                = 1.0
    m(15)%nra                   = 198
    m(15)%ra_inc                = 1000.0
    !
    m(15)%nel_default(1)        = 4
    m(15)%el_arr_default(:,1)   = unused_value
    m(15)%el_arr_default(1:m(15)%nel_default(1),1)  = &
         (/0.5, 1.5, 2.5, 3.5/)
    !
    m(15)%nel                   = m(15)%nel_default(1)
    m(15)%el_arr(:)             = unused_value
    m(15)%el_arr(1:m(15)%nel)   = m(15)%el_arr_default(1:m(15)%nel_default(1),1)
    !
    dbzparams(15)%station_id    = 16996
    !
    CALL set_scanname ( m(15) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 16: m(16) and dbzparams(16): Macaion
    !
    m(16)%station_name          = 'MAC'
    m(16)%station_id            = 16997
    m(16)%lambda                = 0.055
    m(16)%lat                   = 46.495
    m(16)%lon                   = 11.21
    m(16)%alt_agl_mod           = 12.0
    m(16)%alt_msl               = 1876.0
    m(16)%alt_msl_true          = 1876.0
    m(16)%mds_Z0                = -20.0
    m(16)%mds_r0                = 10000.0
    m(16)%dt_obs(1)             = 600.0         ! 10 minutes, traditional notation <incr>,<miss>,<miss>
    !
    m(16)%az_start              = 0.0
    m(16)%naz                   = 360
    m(16)%az_inc                = 1.0
    m(16)%nra                   = 120
    m(16)%ra_inc                = 1000.0
    !
    m(16)%nel_default(1)        = 7
    m(16)%el_arr_default(:,1)   = unused_value
    m(16)%el_arr_default(1:m(16)%nel_default(1),1)  = &
         (/0.0, 1.0, 2.0, 3.0, 4.0, 8.0, 14.0/)
    !
    m(16)%nel                   = m(16)%nel_default(1)
    m(16)%el_arr(:)             = unused_value
    m(16)%el_arr(1:m(16)%nel)   = m(16)%el_arr_default(1:m(16)%nel_default(1),1)
    !
    dbzparams(16)%station_id    = 16997
    !
    CALL set_scanname ( m(16) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================

  END SUBROUTINE get_meta_network_italy

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the German (DWD) network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_dwd (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_dwd)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_dwd)

    TYPE(radar_meta_type) :: mproto
    INTEGER               :: i

    !========================================================
    ! .. start from the protoype to generate volscans:
    !========================================================

!!$ UB: this failed on the CRAY, only the first element m(1) was properly set:
!!$    m(:)         = get_meta_proto_dwd ()
    mproto       = get_meta_proto ( icountry=i_dwd )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !=============================================================
    ! Individual settings for each station, default scan strategy:
    !    scanname = 'PPI0080'
    !=============================================================
    !
    ! STATION 1: m(1) and dbzparams(1): BOO Boosted
    !
    m(1)%station_name        = 'BOO'
    m(1)%station_id          = 10132
    m(1)%lambda = 0.055
    m(1)%lat = 54.00438
    m(1)%lon = 10.04687
    m(1)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 125.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 10132
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): DRS Dresden
    !
    m(2)%station_name        = 'DRS'
    m(2)%station_id          = 10488
    m(2)%lambda = 0.055
    m(2)%lat = 51.12452
    m(2)%lon = 13.76867
    m(2)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 262.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 10488
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): EIS Eisberg
    !
    m(3)%station_name        = 'EIS'
    m(3)%station_id          = 10780
    m(3)%lambda = 0.055
    m(3)%lat = 49.54066
    m(3)%lon = 12.40278
    m(3)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true     = 799.0
    m(3)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    !
    dbzparams(3)%station_id          = 10780
    !
    !========================================================
    !
    ! STATION 4: m(4) and dbzparams(4): EMD Emden
    !
    m(4)%station_name        = 'EMD'
    m(4)%station_id          = 10204
    m(4)%lambda = 0.055
    m(4)%lat = 53.33872
    m(4)%lon = 7.02377
    m(4)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(4)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(4)%alt_msl_true     = 58.0
    m(4)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(4)%mds_r0 = 10000.0
    !
    dbzparams(4)%station_id          = 10204

    !========================================================
    !
    ! STATION 5: m(5) and dbzparams(5): ESS Essen (new site, slightly moved tower)
    !
    m(5)%station_name        = 'ESS'
    m(5)%station_id          = 10410
    m(5)%lambda = 0.055
    m(5)%lat = 51.40563
    m(5)%lon = 6.96712
    m(5)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(5)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(5)%alt_msl_true     = 185.0
    m(5)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(5)%mds_r0 = 10000.0
    !
    dbzparams(5)%station_id          = 10410
    !
    !========================================================
    !
    ! STATION 6: m(6) and dbzparams(6): FBG Feldberg
    !
    m(6)%station_name        = 'FBG'
    m(6)%station_id          = 10908
    m(6)%lambda = 0.055
    m(6)%lat = 47.87361
    m(6)%lon = 8.00361
    m(6)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(6)%alt_msl     = 1516.0 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(6)%alt_msl_true     = 1516.0
    m(6)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(6)%mds_r0 = 10000.0
    !
    dbzparams(6)%station_id          = 10908
    !
    !========================================================
    !
    ! STATION 7: m(7) and dbzparams(7): FLD Flechtdorf
    !
    m(7)%station_name        = 'FLD'
    m(7)%station_id          = 10440
    m(7)%lambda = 0.055
    m(7)%lat = 51.31119
    m(7)%lon = 8.80206
    m(7)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(7)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(7)%alt_msl_true     = 623.0
    m(7)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(7)%mds_r0 = 10000.0
    !
    dbzparams(7)%station_id          = 10440
    !
    !========================================================
    !
    ! STATION 8: m(8) and dbzparams(8): HNR Hannover (new)
    !
    m(8)%station_name        = 'HNR'
    m(8)%station_id          = 10339
    m(8)%lambda = 0.055
    m(8)%lat = 52.46008
    m(8)%lon = 9.69452
    m(8)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(8)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(8)%alt_msl_true     = 98.0
    m(8)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(8)%mds_r0 = 10000.0
    !
    dbzparams(8)%station_id          = 10339
    !
    !========================================================
    !
    ! STATION 9: m(9) and dbzparams(9): NEU Neuhaus
    !
    m(9)%station_name        = 'NEU'
    m(9)%station_id          = 10557
    m(9)%lambda = 0.055
    m(9)%lat = 50.50012
    m(9)%lon = 11.13504
    m(9)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(9)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(9)%alt_msl_true     = 878.0
    m(9)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(9)%mds_r0 = 10000.0
    !
    dbzparams(9)%station_id          = 10557
    !
    !========================================================
    !
    ! STATION 10: m(10) and dbzparams(10): NHB Neuheilenbach
    !
    m(10)%station_name        = 'NHB'
    m(10)%station_id          = 10605
    m(10)%lambda = 0.055
    m(10)%lat = 50.10973
    m(10)%lon = 6.54853
    m(10)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(10)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(10)%alt_msl_true     = 585.0
    m(10)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(10)%mds_r0 = 10000.0
    !
    dbzparams(10)%station_id          = 10605
    !
    !========================================================
    !
    ! STATION 11: m(11) and dbzparams(11): OFT Offenthal
    !
    m(11)%station_name        = 'OFT'
    m(11)%station_id          = 10629
    m(11)%lambda = 0.055
    m(11)%lat = 49.98478
    m(11)%lon = 8.71293
    m(11)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(11)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(11)%alt_msl_true     = 246.0
    m(11)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(11)%mds_r0 = 10000.0
    !
    dbzparams(11)%station_id          = 10629
    !
    !========================================================
    !
    ! STATION 12: m(12) and dbzparams(12): PRO Proetzel
    !
    m(12)%station_name        = 'PRO'
    m(12)%station_id          = 10392
    m(12)%lambda = 0.055
    m(12)%lat = 52.64867
    m(12)%lon = 13.85821
    m(12)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(12)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(12)%alt_msl_true     = 194.0
    m(12)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(12)%mds_r0 = 10000.0
    !
    dbzparams(12)%station_id          = 10392
    !
    !========================================================
    !
    ! STATION 13: m(13) and dbzparams(13): MEM Memmingen
    !
    m(13)%station_name        = 'MEM'
    m(13)%station_id          = 10950
    m(13)%lambda = 0.055
    m(13)%lat = 48.04214
    m(13)%lon = 10.21924
    m(13)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(13)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(13)%alt_msl_true     = 724.0
    m(13)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(13)%mds_r0 = 10000.0
    !
    dbzparams(13)%station_id          = 10950
    !
    !========================================================
    !
    ! STATION 14: m(14) and dbzparams(14): ROS Rostock
    !
    m(14)%station_name        = 'ROS'
    m(14)%station_id          = 10169
    m(14)%lambda = 0.055
    m(14)%lat = 54.17566
    m(14)%lon = 12.05808
    m(14)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(14)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(14)%alt_msl_true     = 37.0
    m(14)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(14)%mds_r0 = 10000.0
    !
    dbzparams(14)%station_id          = 10169
    !
    !========================================================
    !
    ! STATION 15: m(15) and dbzparams(15): ISN Isen (old name: Schnaupping)
    !
    m(15)%station_name        = 'ISN'
    m(15)%station_id          = 10873
    m(15)%lambda = 0.055
    m(15)%lat = 48.17470
    m(15)%lon = 12.10179
    m(15)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(15)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(15)%alt_msl_true     = 678.0
    m(15)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(15)%mds_r0 = 10000.0
    !
    dbzparams(15)%station_id          = 10873
    !
    !========================================================
    !
    ! STATION 16: m(16) and dbzparams(16): TUR Tuerkheim
    !
    m(16)%station_name        = 'TUR'
    m(16)%station_id          = 10832
    m(16)%lambda = 0.055
    m(16)%lat = 48.58528
    m(16)%lon = 9.78278
    m(16)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(16)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(16)%alt_msl_true     = 768.0
    m(16)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(16)%mds_r0 = 10000.0
    !
    dbzparams(16)%station_id          = 10832
    !
    !========================================================
    !
    ! STATION 17: m(17) and dbzparams(17): UMD Ummendorf
    !
    m(17)%station_name        = 'UMD'
    m(17)%station_id          = 10356
    m(17)%lambda = 0.055
    m(17)%lat = 52.16009
    m(17)%lon = 11.17609
    m(17)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(17)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(17)%alt_msl_true     = 183.0
    m(17)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(17)%mds_r0 = 10000.0
    !
    dbzparams(17)%station_id          = 10356
    !
    !========================================================
    !
    ! STATION 18: m(18) and dbzparams(18): BLN Berlin Tempelhof
    !
    m(18)%station_name        = 'BLN'
    m(18)%station_id          = 10384
    m(18)%lambda = 0.055
    m(18)%lat = 52.47787
    m(18)%lon = 13.38697
    m(18)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(18)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(18)%alt_msl_true     = 80.0
    m(18)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(18)%mds_r0 = 10000.0
    !
    dbzparams(18)%station_id          = 10384
    !
    !========================================================
    !
    ! STATION 19: m(19) and dbzparams(19): ASE Essen (Fail Safe Radar system)
    !
    m(19)%station_name        = 'ASE'
    m(19)%station_id          = 10412
    m(19)%lambda = 0.055
    m(19)%lat = 51.40514
    m(19)%lon = 6.96384
    m(19)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(19)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(19)%alt_msl_true     = 188.0
    m(19)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(19)%mds_r0 = 10000.0
    !
    dbzparams(19)%station_id          = 10412
    !
    !========================================================
    !
    ! STATION 20: m(20) and dbzparams(20): HAM Hamburg
    !
    m(20)%station_name        = 'HAM'
    m(20)%station_id          = 10147
    m(20)%lambda = 0.055
    m(20)%lat = 53.62123
    m(20)%lon = 9.99662
    m(20)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(20)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(20)%alt_msl_true     = 46.0
    m(20)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(20)%mds_r0 = 10000.0
    !
    dbzparams(20)%station_id          = 10147
    !
    !========================================================
    !
    ! STATION 21: m(21) and dbzparams(21): FRI Frankfurt Flughafen
    !
    m(21)%station_name        = 'FRI'
    m(21)%station_id          = 10630
    m(21)%lambda = 0.055
    m(21)%lat = 50.02244
    m(21)%lon = 8.55854
    m(21)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(21)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(21)%alt_msl_true     = 145.0
    m(21)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(21)%mds_r0 = 10000.0
    !
    dbzparams(21)%station_id          = 10630
    !
    !========================================================
    !
    ! STATION 22: m(22) and dbzparams(22): MUC Muenchen Fuerholzen
    !
    m(22)%station_name        = 'MUC'
    m(22)%station_id          = 10871
    m(22)%lambda = 0.055
    m(22)%lat = 48.33638
    m(22)%lon = 11.61170
    m(22)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(22)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(22)%alt_msl_true     = 511.0
    m(22)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(22)%mds_r0 = 10000.0
    !
    dbzparams(22)%station_id          = 10871
    !
    !========================================================
    !
    ! STATION 23: m(23) and dbzparams(23): ASW Rostock (Fail Safe Radar system)
    !
    m(23)%station_name        = 'ASW'
    m(23)%station_id          = 10089
    m(23)%lambda = 0.055
    m(23)%lat = 54.17312
    m(23)%lon = 12.10703
    m(23)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(23)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(23)%alt_msl_true     = 35.0
    m(23)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(23)%mds_r0 = 10000.0
    !
    dbzparams(23)%station_id          = 10089
    !
    !========================================================
    !
    ! STATION 24: m(24) and dbzparams(24): ASF Feldberg (Fail Safe Radar system)
    !
    m(24)%station_name        = 'ASF'
    m(24)%station_id          = 10907
    m(24)%lambda = 0.055
    m(24)%lat = 47.87260
    m(24)%lon = 8.00683
    m(24)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(24)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(24)%alt_msl_true     = 1502.0
    m(24)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(24)%mds_r0 = 10000.0
    !
    dbzparams(24)%station_id          = 10907
    !
    !========================================================
    !
    ! STATION 25: m(25) and dbzparams(25): HAN Hannover (old)
    !
    m(25)%station_name        = 'HAN'
    m(25)%station_id          = 10338
    m(25)%lambda = 0.055
    m(25)%lat = 52.46306
    m(25)%lon = 9.69832
    m(25)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(25)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(25)%alt_msl_true     = 81.0
    m(25)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(25)%mds_r0 = 10000.0
    !
    dbzparams(25)%station_id          = 10338
    !
    !========================================================
    !
    ! STATION 26: m(26) and dbzparams(26): MOHP Hohenpeissenberg
    !
    m(26)%station_name        = 'MHP'
    m(26)%station_id          = 10962
    m(26)%lambda = 0.055
    m(26)%lat = 47.80151
    m(26)%lon = 11.00929
    m(26)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(26)%alt_msl     = 1006.0 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(26)%alt_msl_true     = 1006.0
    m(26)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(26)%mds_r0 = 10000.0
    !
    dbzparams(26)%station_id          = 10962
    !
    !========================================================
    !
    ! STATION 27: m(27) and dbzparams(27): ASB Ausfallsicherungsradar Borkum
    !
    m(27)%station_name        = 'ASB'
    m(27)%station_id          = 10103
    m(27)%lambda = 0.055
    m(27)%lat = 53.56401
    m(27)%lon = 6.74829
    m(27)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(27)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(27)%alt_msl_true     = 36.0
    m(27)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(27)%mds_r0 = 10000.0
    !
    dbzparams(27)%station_id          = 10103
    !
    !========================================================
    !
    ! STATION 28: m(28) and dbzparams(28): XXX Fake radar for OSSEs
    !
    m(28)%station_name        = 'YY1'
    m(28)%station_id          = 10991
    m(28)%lambda = 0.055
    m(28)%lat = 0.0
    m(28)%lon = 0.0
    m(28)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(28)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(28)%alt_msl_true     = 50.0
    m(28)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(28)%mds_r0 = 10000.0
    !
    dbzparams(28)%station_id          = 10991
    !
    !========================================================
    !
    ! STATION 29: m(29) and dbzparams(29): XXX Fake radar for OSSEs
    !
    m(29)%station_name        = 'YY2'
    m(29)%station_id          = 10992
    m(29)%lambda = 0.055
    m(29)%lat = 0.0
    m(29)%lon = 0.0
    m(29)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(29)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(29)%alt_msl_true     = 50.0
    m(29)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(29)%mds_r0 = 10000.0
    !
    dbzparams(29)%station_id          = 10992
    !
    !========================================================
    !
    ! STATION 30: m(30) and dbzparams(30): XXX Fake radar for OSSEs
    !
    m(30)%station_name        = 'YY3'
    m(30)%station_id          = 10993
    m(30)%lambda = 0.055
    m(30)%lat = 0.0
    m(30)%lon = 0.0
    m(30)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(30)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(30)%alt_msl_true     = 50.0
    m(30)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(30)%mds_r0 = 10000.0
    !
    dbzparams(30)%station_id          = 10993
    !
    !========================================================
    !
    ! STATION 30: m(30) and dbzparams(30): XXX Fake radar for OSSEs
    !
    m(31)%station_name        = 'YY4'
    m(31)%station_id          = 10994
    m(31)%lambda = 0.055
    m(31)%lat = 0.0
    m(31)%lon = 0.0
    m(31)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(31)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(31)%alt_msl_true     = 50.0
    m(31)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(31)%mds_r0 = 10000.0
    !
    dbzparams(31)%station_id          = 10994
    !
    !========================================================

    !========================================================
    ! .. Add the precip scans for the 31 stations:
    !========================================================

    DO i= 32, 62
      ! clone the first 31 stations
      m(i)                    = m(i-31)
      dbzparams(i)%station_id = m(i)%station_id
      ! Reset scan strategy to precip scans:
      m(i)%nel = 1
      m(i)%el_arr(:) = unused_value
      m(i)%el_arr(1) = 0.8           !  precipitation scan is indicated by 0.8°
      CALL set_scanname ( m(i) )     ! scanname = 'PRECIP'
    END DO

    !========================================================
    !
    ! STATION 63: m(63) and dbzparams(63): KIT Karlsruhe C-Band
    !
    !    NOTE: files are mmms, but do not conform to ODIM-hdf5 attribute/dataset naming conventions!
    !    Filename convention:
    !
    !      scan-sidpol-120km-14_20001_20230712124001_00.h5
    !
    m(63) = get_meta_proto_kitcband ( icountry=i_dwd )
    
    m(63)%station_name        = 'KAR'
    m(63)%station_id          = 20001
    m(63)%lambda = 0.055
    m(63)%lat = 49.094167
    m(63)%lon = 8.436667
    m(63)%alt_agl_mod = 28.0     ! only effective if alt_msl < -9000
    m(63)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(63)%alt_msl_true     = 148.0
    m(63)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(63)%mds_r0 = 10000.0
    !
    dbzparams(63)%station_id          = 20001
    !
    !========================================================
    !
    ! STATION 64: m(64) and dbzparams(64): KIT-Cube X-Band at Bondorf for Swabian Moses Campaign
    !
    !       NOTE: files are smms in ODIM-hdf5, rename original files to OPERA T_PA convention:
    !
    !       2023071213050900dBZ.vol.h5  -> T_PAGZ99_C_KITC_20230712130509.h5
    !       2023071213050900V.vol.h5    -> T_PAHZ99_C_KITC_20230712130509.h5
    !
    !
    m(64) = get_meta_proto_kitcube ( icountry=i_dwd )
    
    m(64)%station_name        = 'BON'
    m(64)%station_id          = 10000
    m(64)%lambda = 0.032
    m(64)%lat = 47.8252
    m(64)%lon = 8.34774
    m(64)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(64)%alt_msl     = 901.0 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(64)%alt_msl_true     = 901.0
    m(64)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(64)%mds_r0 = 10000.0
    !
    dbzparams(64)%station_id          = 10000
    !

  END SUBROUTINE get_meta_network_dwd

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Swiss (MCH) network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_swiss (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_swiss)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_swiss)

    TYPE(radar_meta_type) :: mproto
    INTEGER               :: i

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_meteoswiss )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !============================================================
    ! Individual setting for each station, default scan strategy:
    !    scanname = 'PPI0120'
    !============================================================
    !
    ! STATION 1: m(1) and dbzparams(1): ALB Albis
    !
    m(1)%station_name        = 'ALB'
    m(1)%station_id          = 6661
    m(1)%lambda = 0.055
    m(1)%lat = 47.284333
    m(1)%lon = 8.512
    m(1)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 938.0
    m(1)%mds_Z0 = -14.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 6661
    !
    !========================================================
    !
    ! STATION : m(2) and dbzparams(2): DOL La Dole
    !
    m(2)%station_name        = 'DOL'
    m(2)%station_id          = 6699
    m(2)%lambda = 0.055
    m(2)%lat = 46.425113
    m(2)%lon = 6.099415
    m(2)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 1682.0
    m(2)%mds_Z0 = -14.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 6699
    !
    !========================================================
    !
    ! STATION : m(3) and dbzparams(3): MLE Monte Lema
    !
    m(3)%station_name        = 'MLE'
    m(3)%station_id          = 6768
    m(3)%lambda = 0.055
    m(3)%lat = 46.040761
    m(3)%lon = 8.833217
    m(3)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true = 1626.0
    m(3)%mds_Z0 = -14.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    !
    dbzparams(3)%station_id          = 6768
    !
    !========================================================
    !
    ! STATION : m(4) and dbzparams(4): PMO Plaine Morte
    !
    m(4)%station_name        = 'PMO'
    m(4)%station_id          = 6726
    m(4)%lambda = 0.055
    m(4)%lat = 46.370646
    m(4)%lon = 7.486552
    m(4)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(4)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(4)%alt_msl_true = 2937.0
    m(4)%mds_Z0 = -14.0         ! Minimum detectable signal [dBZ]
    m(4)%mds_r0 = 10000.0
    !
    dbzparams(4)%station_id          = 6726
    !
    !========================================================
    !
    ! STATION : m(5) and dbzparams(5): WEI Weissfluhgipfel
    !
    m(5)%station_name        = 'WEI'
    m(5)%station_id          = 6776
    m(5)%lambda = 0.055
    m(5)%lat = 46.834974
    m(5)%lon = 9.794458
    m(5)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(5)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(5)%alt_msl_true = 2840.0
    m(5)%mds_Z0 = -14.0         ! Minimum detectable signal [dBZ]
    m(5)%mds_r0 = 10000.0
    !
    dbzparams(5)%station_id          = 6776
    !
    !========================================================

    !========================================================
    ! .. Add the reduced scan strategy (5 lowest elev) for data
    !    from the OPERA data hub:
    !========================================================

    DO i = 6, 10
      ! clone the first 5 stations
      m(i)                    = m(i-5)
      dbzparams(i)%station_id = m(i)%station_id
      ! Reset scan strategy to precip scans:
      m(i)%nel                = m(i)%nel_default(2)
      m(i)%el_arr(:)          = unused_value
      m(i)%el_arr(1:m(i)%nel) = m(i)%el_arr_default(1:m(i)%nel_default(2),2)
      CALL set_scanname ( m(i) )  ! scanname = 'PPI0011'
    END DO


  END SUBROUTINE get_meta_network_swiss

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Italian (Emilia-Romagna) network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_belgium (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_belgium)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_belgium)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_belgium )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Wideumont
    !
    m(1)%station_name        = 'WID'
    m(1)%station_id          = 6477
    m(1)%lambda = 0.055
    m(1)%lat = 49.9135
    m(1)%lon = 5.5044
    m(1)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 590.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 6477
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Jabekke
    !
    m(2)%station_name        = 'JAB'
    m(2)%station_id          = 6410
    m(2)%lambda = 0.055
    m(2)%lat = 51.1919
    m(2)%lon = 3.0641
    m(2)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 4.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 6410
    !
    !========================================================

  END SUBROUTINE get_meta_network_belgium

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Danish network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_denmark (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_denmark)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_denmark)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_denmark )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Sindal
    !
    m(1)%station_name        = 'SIN'
    m(1)%station_id          = 6034
    m(1)%lambda = 0.055
    m(1)%lat = 57.4893
    m(1)%lon = 10.1365
    m(1)%alt_agl_mod = 10.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 109.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 6034
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Römö
    !
    m(2)%station_name        = 'ROM'
    m(2)%station_id          = 6096
    m(2)%lambda = 0.055
    m(2)%lat = 55.1731
    m(2)%lon = 8.552
    m(2)%alt_agl_mod = 10.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 15.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 6096
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): Virring
    !
    m(3)%station_name        = 'VIR'
    m(3)%station_id          = 6103
    m(3)%lambda = 0.055
    m(3)%lat = 56.024
    m(3)%lon = 10.0246
    m(3)%alt_agl_mod = 10.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true = 142.0
    m(3)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    !
    dbzparams(3)%station_id          = 6103
    !
    !========================================================
    !
    ! STATION 4: m(4) and dbzparams(4): Stevns
    !
    m(4)%station_name        = 'STE'
    m(4)%station_id          = 6173
    m(4)%lambda = 0.055
    m(4)%lat = 55.3262
    m(4)%lon = 12.4493
    m(4)%alt_agl_mod = 10.0     ! only effective if alt_msl < -9000
    m(4)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(4)%alt_msl_true = 53.0
    m(4)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(4)%mds_r0 = 10000.0
    !
    dbzparams(4)%station_id          = 6173
    !
    !========================================================
    !
    ! STATION 5: m(5) and dbzparams(5): Bornholm
    !
    m(5)%station_name        = 'BOR'
    m(5)%station_id          = 6194
    m(5)%lambda = 0.055
    m(5)%lat = 55.1127
    m(5)%lon = 14.8875
    m(5)%alt_agl_mod = 10.0     ! only effective if alt_msl < -9000
    m(5)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(5)%alt_msl_true = 171.0
    m(5)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(5)%mds_r0 = 10000.0
    !
    dbzparams(5)%station_id          = 6194
    !
    !
    !========================================================

  END SUBROUTINE get_meta_network_denmark

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Danish network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_france (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_france)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_france)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_france )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto


    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): AJA
    !
    m(1)%station_name        = 'AJA'
    m(1)%station_id          = 7760
    m(1)%lambda = 0.055
    m(1)%lat = 41.95310
    m(1)%lon = 8.70054
    m(1)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 784.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    m(1)%nel_default(1)      = 7
    m(1)%el_arr_default(:,1) = unused_value
    m(1)%el_arr_default(1:m(1)%nel_default(1),1)  = (/0.4, 0.8, 1.4, 2.2, 3.6, 4.8, 6.0/)
    m(1)%nel                 = m(1)%nel_default(1)
    m(1)%el_arr(:)           = unused_value
    m(1)%el_arr(1:m(1)%nel)  = m(1)%el_arr_default(1:m(1)%nel_default(1),1)
    !
    dbzparams(1)%station_id          = 7760
    !
    CALL set_scanname ( m(1) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): ABB
    !
    m(2)%station_name        = 'ABB'
    m(2)%station_id          = 7005
    m(2)%lambda = 0.055
    m(2)%lat = 50.13600
    m(2)%lon = 1.83470
    m(2)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 83.6
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    m(2)%nel_default(1)      = 11
    m(2)%el_arr_default(:,1) = unused_value
    m(2)%el_arr_default(1:m(2)%nel_default(1),1)  = (/0.4, 0.9, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5/)
    m(2)%nel                 = m(2)%nel_default(1)
    m(2)%el_arr(:)           = unused_value
    m(2)%el_arr(1:m(2)%nel)  = m(2)%el_arr_default(1:m(2)%nel_default(1),1)
    !
    dbzparams(2)%station_id          = 7005
    !
    CALL set_scanname ( m(2) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): BOR
    !
    m(3)%station_name        = 'BOR'
    m(3)%station_id          = 7510
    m(3)%lambda = 0.055
    m(3)%lat = 44.83150
    m(3)%lon = -0.69188
    m(3)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true = 70.0
    m(3)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    m(3)%nel_default(1)      = 9
    m(3)%el_arr_default(:,1) = unused_value
    m(3)%el_arr_default(1:m(3)%nel_default(1),1)  = (/0.4, 0.8, 1.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0/)
    m(3)%nel                 = m(3)%nel_default(1)
    m(3)%el_arr(:)           = unused_value
    m(3)%el_arr(1:m(3)%nel)  = m(3)%el_arr_default(1:m(3)%nel_default(1),1)
    !
    dbzparams(3)%station_id          = 7510
    !
    CALL set_scanname ( m(3) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 4: m(4) and dbzparams(4): BOU
    !
    m(4)%station_name        = 'BOU'
    m(4)%station_id          = 7255
    m(4)%lambda = 0.055
    m(4)%lat = 47.05860
    m(4)%lon = 2.35955
    m(4)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(4)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(4)%alt_msl_true = 173.5
    m(4)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(4)%mds_r0 = 10000.0
    m(4)%nel_default(1)      = 11
    m(4)%el_arr_default(:,1) = unused_value
    m(4)%el_arr_default(1:m(4)%nel_default(1),1)  = (/0.7, 1.3, 2.2, 3.2, 4.2, 5.2, 6.0, 7.0, 8.0, 9.0, 10.0/)
    m(4)%nel                 = m(4)%nel_default(1)
    m(4)%el_arr(:)           = unused_value
    m(4)%el_arr(1:m(4)%nel)  = m(4)%el_arr_default(1:m(4)%nel_default(1),1)
    !
    dbzparams(4)%station_id          = 7255
    !
    CALL set_scanname ( m(4) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 5: m(5) and dbzparams(5): GRE
    !
    m(5)%station_name        = 'GRE'
    m(5)%station_id          = 7436
    m(5)%lambda = 0.055
    m(5)%lat = 45.10440
    m(5)%lon = 1.36967
    m(5)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(5)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(5)%alt_msl_true = 360.3
    m(5)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(5)%mds_r0 = 10000.0
    m(5)%nel_default(1)      = 9
    m(5)%el_arr_default(:,1) = unused_value
    m(5)%el_arr_default(1:m(5)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6, 7.6/)
    m(5)%nel                 = m(5)%nel_default(1)
    m(5)%el_arr(:)           = unused_value
    m(5)%el_arr(1:m(5)%nel)  = m(5)%el_arr_default(1:m(5)%nel_default(1),1)
    !
    dbzparams(5)%station_id          = 7436
    !
    CALL set_scanname ( m(5) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 6: m(6) and dbzparams(6): CAE
    !
    m(6)%station_name        = 'CAE'
    m(6)%station_id          = 7027
    m(6)%lambda = 0.055
    m(6)%lat = 48.92720
    m(6)%lon = -0.14955
    m(6)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(6)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(6)%alt_msl_true = 167.4
    m(6)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(6)%mds_r0 = 10000.0
    m(6)%nel_default(1)      = 11
    m(6)%el_arr_default(:,1) = unused_value
    m(6)%el_arr_default(1:m(6)%nel_default(1),1)  = (/0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/)
    m(6)%nel                 = m(6)%nel_default(1)
    m(6)%el_arr(:)           = unused_value
    m(6)%el_arr(1:m(6)%nel)  = m(6)%el_arr_default(1:m(6)%nel_default(1),1)
    !
    dbzparams(6)%station_id          = 7027
    !
    CALL set_scanname ( m(6) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 7: m(7) and dbzparams(7): NAN
    !
    m(7)%station_name        = 'NAN'
    m(7)%station_id          = 7180
    m(7)%lambda = 0.055
    m(7)%lat = 48.71580
    m(7)%lon = 6.58156
    m(7)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(7)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(7)%alt_msl_true = 296.3
    m(7)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(7)%mds_r0 = 10000.0
    m(7)%nel_default(1)      = 11
    m(7)%el_arr_default(:,1) = unused_value
    m(7)%el_arr_default(1:m(7)%nel_default(1),1)  = (/0.7, 1.3, 1.9, 2.6, 3.6, 4.6, 5.6, 6.6, 7.6, 8.6, 9.6/)
    m(7)%nel                 = m(7)%nel_default(1)
    m(7)%el_arr(:)           = unused_value
    m(7)%el_arr(1:m(7)%nel)  = m(7)%el_arr_default(1:m(7)%nel_default(1),1)
    !
    dbzparams(7)%station_id          = 7180
    !
    CALL set_scanname ( m(7) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 8: m(8) and dbzparams(8): PLA
    !
    m(8)%station_name        = 'PLA'
    m(8)%station_id          = 7108
    m(8)%lambda = 0.055
    m(8)%lat = 48.46090
    m(8)%lon = -4.42983
    m(8)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(8)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(8)%alt_msl_true = 110.8
    m(8)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(8)%mds_r0 = 10000.0
    m(8)%nel_default(1)      = 9
    m(8)%el_arr_default(:,1) = unused_value
    m(8)%el_arr_default(1:m(8)%nel_default(1),1)  = (/0.4, 0.8, 1.8, 2.8, 4.0, 5.0, 6.0, 9.0, 12.0/)
    m(8)%nel                 = m(8)%nel_default(1)
    m(8)%el_arr(:)           = unused_value
    m(8)%el_arr(1:m(8)%nel)  = m(8)%el_arr_default(1:m(8)%nel_default(1),1)
    !
    dbzparams(8)%station_id          = 7108
    !
    CALL set_scanname ( m(8) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 9: m(9) and dbzparams(9): TOU
    !
    m(9)%station_name        = 'TOU'
    m(9)%station_id          = 7629
    m(9)%lambda = 0.055
    m(9)%lat = 43.57430
    m(9)%lon = 1.37630
    m(9)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(9)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(9)%alt_msl_true = 187.1
    m(9)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(9)%mds_r0 = 10000.0
    m(9)%nel_default(1)      = 11
    m(9)%el_arr_default(:,1) = unused_value
    m(9)%el_arr_default(1:m(9)%nel_default(1),1)  = (/0.8, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5/)
    m(9)%nel                 = m(9)%nel_default(1)
    m(9)%el_arr(:)           = unused_value
    m(9)%el_arr(1:m(9)%nel)  = m(9)%el_arr_default(1:m(9)%nel_default(1),1)
    !
    dbzparams(9)%station_id          = 7629
    !
    CALL set_scanname ( m(9) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 10: m(10) and dbzparams(10): TRA
    !
    m(10)%station_name        = 'TRA'
    m(10)%station_id          = 7145
    m(10)%lambda = 0.055
    m(10)%lat = 48.77460
    m(10)%lon = 2.00831
    m(10)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(10)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(10)%alt_msl_true = 190.2
    m(10)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(10)%mds_r0 = 10000.0
    m(10)%nel_default(1)      = 11
    m(10)%el_arr_default(:,1) = unused_value
    m(10)%el_arr_default(1:m(10)%nel_default(1),1)  = (/0.4, 0.8, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5/)
    m(10)%nel                 = m(10)%nel_default(1)
    m(10)%el_arr(:)           = unused_value
    m(10)%el_arr(1:m(10)%nel) = m(10)%el_arr_default(1:m(10)%nel_default(1),1)
    !
    dbzparams(10)%station_id          = 7145
    !
    CALL set_scanname ( m(10) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 11: m(11) and dbzparams(11): TRO
    !
    m(11)%station_name        = 'TRO'
    m(11)%station_id          = 7168
    m(11)%lambda = 0.055
    m(11)%lat = 48.46210
    m(11)%lon = 4.30932
    m(11)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(11)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(11)%alt_msl_true = 164.8
    m(11)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(11)%mds_r0 = 10000.0
    m(11)%nel_default(1)      = 11
    m(11)%el_arr_default(:,1) = unused_value
    m(11)%el_arr_default(1:m(11)%nel_default(1),1)  = (/0.4, 1.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0/)
    m(11)%nel                 = m(11)%nel_default(1)
    m(11)%el_arr(:)           = unused_value
    m(11)%el_arr(1:m(11)%nel) = m(11)%el_arr_default(1:m(11)%nel_default(1),1)
    !
    dbzparams(11)%station_id          = 7168
    !
    CALL set_scanname ( m(11) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 12: m(12) and dbzparams(12): LEP
    !
    m(12)%station_name        = 'LEP'
    m(12)%station_id          = 7471
    m(12)%lambda = 0.055
    m(12)%lat = 45.28920
    m(12)%lon = 3.70948
    m(12)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(12)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(12)%alt_msl_true = 1143.4
    m(12)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(12)%mds_r0 = 10000.0
    m(12)%nel_default(1)      = 11
    m(12)%el_arr_default(:,1) = unused_value
    m(12)%el_arr_default(1:m(12)%nel_default(1),1)  = (/0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/)
    m(12)%nel                 = m(12)%nel_default(1)
    m(12)%el_arr(:)           = unused_value
    m(12)%el_arr(1:m(12)%nel) = m(12)%el_arr_default(1:m(12)%nel_default(1),1)
    !
    dbzparams(12)%station_id          = 7471
    !
    CALL set_scanname ( m(12) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 13: m(13) and dbzparams(13): TRE
    !
    m(13)%station_name        = 'TRE'
    m(13)%station_id          = 7223
    m(13)%lambda = 0.055
    m(13)%lat = 47.33740
    m(13)%lon = -1.65632
    m(13)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(13)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(13)%alt_msl_true = 81.0
    m(13)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(13)%mds_r0 = 10000.0
    m(13)%nel_default(1)      = 11
    m(13)%el_arr_default(:,1) = unused_value
    m(13)%el_arr_default(1:m(13)%nel_default(1),1)  = (/0.4, 0.8, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0/)
    m(13)%nel                 = m(13)%nel_default(1)
    m(13)%el_arr(:)           = unused_value
    m(13)%el_arr(1:m(13)%nel) = m(13)%el_arr_default(1:m(13)%nel_default(1),1)
    !
    dbzparams(13)%station_id          = 7223
    !
    CALL set_scanname ( m(13) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 14: m(14) and dbzparams(14): NIZ
    !
    m(14)%station_name        = 'NIZ'
    m(14)%station_id          = 7381
    m(14)%lambda = 0.055
    m(14)%lat = 46.06780
    m(14)%lon = 4.44535
    m(14)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(14)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(14)%alt_msl_true = 919.8
    m(14)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(14)%mds_r0 = 10000.0
    m(14)%nel_default(1)      = 11
    m(14)%el_arr_default(:,1) = unused_value
    m(14)%el_arr_default(1:m(14)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6, 7.6, 8.6, 9.6/)
    m(14)%nel                 = m(14)%nel_default(1)
    m(14)%el_arr(:)           = unused_value
    m(14)%el_arr(1:m(14)%nel) = m(14)%el_arr_default(1:m(14)%nel_default(1),1)
    !
    dbzparams(14)%station_id          = 7381
    !
    CALL set_scanname ( m(14) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 15: m(15) and dbzparams(15): MCL
    !
    m(15)%station_name        = 'MCL'
    m(15)%station_id          = 7637
    m(15)%lambda = 0.055
    m(15)%lat = 43.99050
    m(15)%lon = 2.60962
    m(15)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(15)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(15)%alt_msl_true = 678.7
    m(15)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(15)%mds_r0 = 10000.0
    m(15)%nel_default(1)      = 9
    m(15)%el_arr_default(:,1) = unused_value
    m(15)%el_arr_default(1:m(15)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5/)
    m(15)%nel                 = m(15)%nel_default(1)
    m(15)%el_arr(:)           = unused_value
    m(15)%el_arr(1:m(15)%nel) = m(15)%el_arr_default(1:m(15)%nel_default(1),1)
    !
    dbzparams(15)%station_id          = 7637
    !
    CALL set_scanname ( m(15) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 16: m(16) and dbzparams(16): AVE
    !
    m(16)%station_name        = 'AVE'
    m(16)%station_id          = 7083
    m(16)%lambda = 0.055
    m(16)%lat = 50.12830
    m(16)%lon = 3.81181
    m(16)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(16)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(16)%alt_msl_true = 208.8
    m(16)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(16)%mds_r0 = 10000.0
    m(16)%nel_default(1)      = 8
    m(16)%el_arr_default(:,1) = unused_value
    m(16)%el_arr_default(1:m(16)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.6, 3.6, 4.6, 6.0, 8.0/)
    m(16)%nel                 = m(16)%nel_default(1)
    m(16)%el_arr(:)           = unused_value
    m(16)%el_arr(1:m(16)%nel) = m(16)%el_arr_default(1:m(16)%nel_default(1),1)
    !
    dbzparams(16)%station_id          = 7083
    !
    CALL set_scanname ( m(16) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 17: m(17) and dbzparams(17): CHE
    !
    m(17)%station_name        = 'CHE'
    m(17)%station_id          = 7336
    m(17)%lambda = 0.055
    m(17)%lat = 46.69860
    m(17)%lon = 0.06555
    m(17)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(17)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(17)%alt_msl_true = 174.4
    m(17)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(17)%mds_r0 = 10000.0
    m(17)%nel_default(1)      = 8
    m(17)%el_arr_default(:,1) = unused_value
    m(17)%el_arr_default(1:m(17)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.6, 3.6, 4.6, 6.0, 8.0/)
    m(17)%nel                 = m(17)%nel_default(1)
    m(17)%el_arr(:)           = unused_value
    m(17)%el_arr(1:m(17)%nel) = m(17)%el_arr_default(1:m(17)%nel_default(1),1)
    !
    dbzparams(17)%station_id          = 7336
    !
    CALL set_scanname ( m(17) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 18: m(18) and dbzparams(18): BLA
    !
    m(18)%station_name        = 'BLA'
    m(18)%station_id          = 7274
    m(18)%lambda = 0.055
    m(18)%lat = 47.35520
    m(18)%lon = 4.77594
    m(18)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(18)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(18)%alt_msl_true = 607.2
    m(18)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(18)%mds_r0 = 10000.0
    m(18)%nel_default(1)      = 8
    m(18)%el_arr_default(:,1) = unused_value
    m(18)%el_arr_default(1:m(18)%nel_default(1),1)  = (/0.5, 1.0, 1.6, 2.6, 3.6, 4.6, 6.0, 8.0/)
    m(18)%nel                 = m(18)%nel_default(1)
    m(18)%el_arr(:)           = unused_value
    m(18)%el_arr(1:m(18)%nel) = m(18)%el_arr_default(1:m(18)%nel_default(1),1)
    !
    dbzparams(18)%station_id          = 7274
    !
    CALL set_scanname ( m(18) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 19: m(19) and dbzparams(19): MOM
    !
    m(19)%station_name        = 'MOM'
    m(19)%station_id          = 7606
    m(19)%lambda = 0.055
    m(19)%lat = 43.62450
    m(19)%lon = -0.60940
    m(19)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(19)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(19)%alt_msl_true = 145.7
    m(19)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(19)%mds_r0 = 10000.0
    m(19)%nel_default(1)      = 9
    m(19)%el_arr_default(:,1) = unused_value
    m(19)%el_arr_default(1:m(19)%nel_default(1),1)  = (/0.4, 1.0, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6, 7.6/)
    m(19)%nel                 = m(19)%nel_default(1)
    m(19)%el_arr(:)           = unused_value
    m(19)%el_arr(1:m(19)%nel) = m(19)%el_arr_default(1:m(19)%nel_default(1),1)
    !
    dbzparams(19)%station_id          = 7606
    !
    CALL set_scanname ( m(19) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================
    !
    ! STATION 20: m(20) and dbzparams(20): MTC
    !
    m(20)%station_name        = 'MTC'
    m(20)%station_id          = 7291
    m(20)%lambda = 0.055
    m(20)%lat = 47.36860
    m(20)%lon = 7.01897
    m(20)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(20)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(20)%alt_msl_true = 925.7
    m(20)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(20)%mds_r0 = 10000.0
    m(20)%nel_default(1)      = 9
    m(20)%el_arr_default(:,1) = unused_value
    m(20)%el_arr_default(1:m(20)%nel_default(1),1)  = (/0.4, 0.7, 1.2, 2.2, 3.0, 3.9, 5.0, 6.0, 7.0/)
    m(20)%nel                 = m(20)%nel_default(1)
    m(20)%el_arr(:)           = unused_value
    m(20)%el_arr(1:m(20)%nel) = m(20)%el_arr_default(1:m(20)%nel_default(1),1)
    !
    dbzparams(20)%station_id          = 7291
    !
    CALL set_scanname ( m(20) )  ! Scan strategy has been adapted, so adapt also the scanname
    !
    !========================================================

  END SUBROUTINE get_meta_network_france

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Polish network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_poland (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_poland)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_poland)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_poland )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Legionowo
    !
    m(1)%station_name        = 'LEG'
    m(1)%station_id          = 12374
    m(1)%lambda = 0.055
    m(1)%lat = 52.405
    m(1)%lon = 20.9609
    m(1)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 119.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 12374
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Ramza
    !
    m(2)%station_name        = 'RAM'
    m(2)%station_id          = 12514
    m(2)%lambda = 0.055
    m(2)%lat = 50.1517
    m(2)%lon = 18.7267
    m(2)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 358.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 12514
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): Pastewnik
    !
    m(3)%station_name        = 'PAS'
    m(3)%station_id          = 12544
    m(3)%lambda = 0.055
    m(3)%lat = 50.892
    m(3)%lon = 16.0395
    m(3)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true = 688.0
    m(3)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    !
    dbzparams(3)%station_id          = 12544
    !
    !========================================================
    !
    ! STATION 4: m(4) and dbzparams(4): Rzeszow
    !
    m(4)%station_name        = 'RZE'
    m(4)%station_id          = 12579
    m(4)%lambda = 0.055
    m(4)%lat = 50.1141
    m(4)%lon = 22.037
    m(4)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(4)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(4)%alt_msl_true = 241.0
    m(4)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(4)%mds_r0 = 10000.0
    !
    dbzparams(4)%station_id          = 12579
    !
    !========================================================
    !
    ! STATION 5: m(5) and dbzparams(5): Poznan
    !
    m(5)%station_name        = 'POZ'
    m(5)%station_id          = 12331
    m(5)%lambda = 0.055
    m(5)%lat = 52.4133
    m(5)%lon = 16.7971
    m(5)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(5)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(5)%alt_msl_true = 130.0
    m(5)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(5)%mds_r0 = 10000.0
    !
    dbzparams(5)%station_id          = 12331
    !
    !========================================================
    !
    ! STATION 6: m(6) and dbzparams(6): Swidwin
    !
    m(6)%station_name        = 'SWI'
    m(6)%station_id          = 12220
    m(6)%lambda = 0.055
    m(6)%lat = 53.7903
    m(6)%lon = 15.8311
    m(6)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(6)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(6)%alt_msl_true = 146.0
    m(6)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(6)%mds_r0 = 10000.0
    !
    dbzparams(6)%station_id          = 12220
    !
    !========================================================
    !
    ! STATION 7: m(7) and dbzparams(7): Gdansk
    !
    m(7)%station_name        = 'GDA'
    m(7)%station_id          = 12151
    m(7)%lambda = 0.055
    m(7)%lat = 54.3843
    m(7)%lon = 18.4563
    m(7)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(7)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(7)%alt_msl_true = 158.0
    m(7)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(7)%mds_r0 = 10000.0
    !
    dbzparams(7)%station_id          = 12151
    !
    !========================================================
    !
    ! STATION 8: m(8) and dbzparams(8): Brzuchania
    !
    m(8)%station_name        = 'BRZ'
    m(8)%station_id          = 12568
    m(8)%lambda = 0.055
    m(8)%lat = 50.3942
    m(8)%lon = 20.0797
    m(8)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(8)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(8)%alt_msl_true = 453.0
    m(8)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(8)%mds_r0 = 10000.0
    !
    dbzparams(8)%station_id          = 12568
    !
    !========================================================

  END SUBROUTINE get_meta_network_poland

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Czech network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_czech (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_czech)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_czech)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_czech )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Skalky
    !
    m(1)%station_name        = 'SKA'
    m(1)%station_id          = 11718
    m(1)%lambda = 0.055
    m(1)%lat = 49.5011
    m(1)%lon = 16.7885
    m(1)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 767.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 11718
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Brdy
    !
    m(2)%station_name        = 'BRD'
    m(2)%station_id          = 11480
    m(2)%lambda = 0.055
    m(2)%lat = 49.6583
    m(2)%lon = 13.8178
    m(2)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 916.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 11480
    !
    !========================================================

  END SUBROUTINE get_meta_network_czech

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Dutch network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_netherlands (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_netherlands)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_netherlands)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_netherlands )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Den Helder
    !
    m(1)%station_name        = 'DHL'
    m(1)%station_id          = 6234
    m(1)%lambda = 0.055
    m(1)%lat = 52.9528
    m(1)%lon = 4.79061
    m(1)%alt_agl_mod = 50.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 55.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 6234
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Herwijnen
    !
    m(2)%station_name        = 'HRW'
    m(2)%station_id          = 6356
    m(2)%lambda = 0.055
    m(2)%lat = 51.8369
    m(2)%lon = 5.1381
    m(2)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 25.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 6356
    !
    !========================================================

  END SUBROUTINE get_meta_network_netherlands

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for each radar of the Slovak network.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_slovakia (dbzparams_proto, m , dbzparams)
    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_slovakia)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_slovakia)

    TYPE(radar_meta_type) :: mproto

    !========================================================
    ! .. start from the protoype:
    !========================================================

    mproto       = get_meta_proto ( icountry=i_slovakia )
    m(:)         = mproto
    dbzparams(:) = dbzparams_proto

    !========================================================
    ! INDIVIDUAL SETTINGS FOR EACH RADAR STATION:
    !========================================================
    !
    ! STATION 1: m(1) and dbzparams(1): Maly Javornik
    !
    m(1)%station_name        = 'JAV'
    m(1)%station_id          = 11812
    m(1)%lambda = 0.055
    m(1)%lat = 48.25610
    m(1)%lon = 17.15310
    m(1)%alt_agl_mod = 15.0     ! only effective if alt_msl < -9000
    m(1)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(1)%alt_msl_true = 600.0
    m(1)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(1)%mds_r0 = 10000.0
    !
    dbzparams(1)%station_id          = 11812
    !
    !========================================================
    !
    ! STATION 2: m(2) and dbzparams(2): Kojsovska Hola
    !
    m(2)%station_name        = 'KOJ'
    m(2)%station_id          = 11958
    m(2)%lambda = 0.055
    m(2)%lat = 48.78290
    m(2)%lon = 20.98730
    m(2)%alt_agl_mod = 16.0     ! only effective if alt_msl < -9000
    m(2)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(2)%alt_msl_true = 1256.0
    m(2)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(2)%mds_r0 = 10000.0
    !
    dbzparams(2)%station_id          = 11958
    !
    !========================================================
    !
    ! STATION 3: m(3) and dbzparams(3): Kubinska Hola
    !
    m(3)%station_name        = 'KUB'
    m(3)%station_id          = 11887
    m(3)%lambda = 0.055
    m(3)%lat = 49.27170
    m(3)%lon = 19.24930
    m(3)%alt_agl_mod = 20.0     ! only effective if alt_msl < -9000
    m(3)%alt_msl     = -9999.99 ! value < -9000 indicates that alt_agl_mod should be used instead
    m(3)%alt_msl_true = 1425.0
    m(3)%mds_Z0 = -20.0         ! Minimum detectable signal [dBZ]
    m(3)%mds_r0 = 10000.0
    !
    dbzparams(3)%station_id          = 11887
    !
    !========================================================

  END SUBROUTINE get_meta_network_slovakia

  !==============================================================================
  !+ Module procedure in radar_src for defining the meta data structures
  !  for all radars known from all background lists.
  !------------------------------------------------------------------------------

  SUBROUTINE get_meta_network_all (dbzparams_proto, m, dbzparams)

    IMPLICIT NONE
    TYPE(t_dbzcalc_params), INTENT(in)  :: dbzparams_proto
    TYPE(radar_meta_type),  INTENT(out) :: m(nradsta_all)
    TYPE(t_dbzcalc_params), INTENT(out) :: dbzparams(nradsta_all)

    INTEGER                    :: iu, io
    CHARACTER (LEN=32)         :: yzroutine

    !========================================================
    ! .. combine from the single countries:
    !========================================================

    yzroutine(:) = ' '
    yzroutine = 'get_meta_network_all'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    iu = 1
    io = nradsta_dwd
    CALL get_meta_network_dwd (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_swiss - 1
    CALL get_meta_network_swiss (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_italy - 1
    CALL get_meta_network_italy (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_belgium - 1
    CALL get_meta_network_belgium (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_denmark - 1
    CALL get_meta_network_denmark (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_france - 1
    CALL get_meta_network_france (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_poland - 1
    CALL get_meta_network_poland (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_czech - 1
    CALL get_meta_network_czech (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_netherlands - 1
    CALL get_meta_network_netherlands (dbzparams_proto, m(iu:io), dbzparams(iu:io))
    iu = io + 1
    io = iu + nradsta_slovakia - 1
    CALL get_meta_network_slovakia (dbzparams_proto, m(iu:io), dbzparams(iu:io))

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE get_meta_network_all

  !==============================================================================
  !+ Module procedure in radar_src for retrieving the azimut-dependend
  !  elevation angles for the DWD precip scans.
  !------------------------------------------------------------------------------

  SUBROUTINE get_elarr_precipscan (rsm, elarr)

    TYPE(radar_meta_type), INTENT(in)              :: rsm
    REAL (KIND=dp), INTENT(inout), ALLOCATABLE :: elarr(:)

    INTEGER                 :: ierr
    CHARACTER(len=cmaxlen)      :: errstring

    ierr = 0
    errstring(:) = ' '

    IF (.NOT.ALLOCATED(elarr)) ALLOCATE (elarr(rsm%naz))
    elarr = 0.0_dp

    IF (rsm%nel /= 1) THEN
      ierr = 1
      WRITE (errstring, '(a,i6.6,a)') 'number of nominal elevations in precip scan is not 1 for station ',rsm%station_id, '!'
    ELSE
      SELECT CASE (rsm%station_id)
        ! include the file "elevations_precipscan.incf", which was automatically created from
        !  "elevations_precipscan.txt" by the script "format_precipscan_f90". The file contains
        !  CASE () clauses for each DWD radar station and fills the vector elarr(:):
        INCLUDE "radar_elevations_precipscan.incf"
      CASE default
        ierr = 3
        WRITE (errstring, '(a,i6.6,a)') 'No nominal elevations for precipscan of station ',rsm%station_id, &
             ' defined! Add a CASE clause!'
      END SELECT

    END IF

    IF (ierr /= 0) THEN
      CALL abort_run (my_radar_id, 13075+ierr, &
           'ERROR: problem in get_elarr_precipscan(): '//TRIM(errstring), &
           'radar_obs_meta_list.f90, get_elarr_precipscan()')
    END IF

  END SUBROUTINE get_elarr_precipscan

  SUBROUTINE set_scanname_m ( rsm )

    TYPE(radar_meta_type),   INTENT(inout) :: rsm

    rsm%scanname(:) = ' '
    IF (rsm%nel == 1 .AND. rsm%icountry == i_dwd .and. &
         ANY( ABS( rsm%el_arr(1) - (/ 0.4_dp, 0.59_dp, 0.8_dp, 1.3_dp /) ) < 1e-6_dp ) ) THEN
      WRITE (rsm%scanname, '(a)') 'PRECIP'
    ELSE
      WRITE (rsm%scanname, '(a,i4.4)') 'PPI', NINT(SUM(rsm%el_arr(1:rsm%nel))/rsm%nel*10.0_dp)
    END IF

  END SUBROUTINE set_scanname_m

  SUBROUTINE set_scanname_o ( rsm )

    TYPE(radar_meta_type_onetime),   INTENT(inout) :: rsm

    rsm%scanname(:) = ' '
    IF (rsm%nel == 1 .AND. rsm%icountry == i_dwd .and. &
         ANY( ABS( rsm%el_arr(1) - (/ 0.4_dp, 0.59_dp, 0.8_dp, 1.3_dp /) ) < 1e-6_dp ) ) THEN
      WRITE (rsm%scanname, '(a)') 'PRECIP'
    ELSE
      WRITE (rsm%scanname, '(a,i4.4)') 'PPI', NINT(SUM(rsm%el_arr(1:rsm%nel))/rsm%nel*10.0_dp)
    END IF

  END SUBROUTINE set_scanname_o


END MODULE radar_obs_meta_list
