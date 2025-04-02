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

MODULE radar_output_country_obs

!------------------------------------------------------------------------------
!
! Description: Modules of the radar forward operator EMVORADO for processing
!              of the various output methods, data and formats:
!              volume data output, feedback file output and reflectivity composite
!              generation and output.
!              This module contains methods to collect/MPI-copy simulated radar data
!              from the compute PEs to the output PEs, for reading obs data files,
  !              for producing superobservations. It uses procedures from other modules
!              for writing feedback files, for writing volume data files (ASCII, NETCDF, or BIN-format)
!              and for writing grib2 reflectivity composites on rotated lat/lon grids.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

!!$ TODO:
!!$ - height above radar station in feedback instead of height MSL? However, station height is also given in the feedback file header ...

  USE radar_kind, ONLY : dp
  
  USE radar_data, ONLY :     &
       miss_threshold, miss_value, miss_thresh_rhv, miss_value_rhv, &
       zero_value, reject_value, shield_value, &
       Z_crit_radar, dBZ_crit_radar, &
       shield_low_threshold, shield_up_threshold, &
       shield_low_thresh_rhv, shield_up_thresh_rhv, &
       missthr_int, missval_int, &
       i_vrad, i_qualvrad, i_dbzh, i_qualdbzh, &
       i_zdr, i_kdp, i_phidp, i_rhv, &
       ncountry_max, i_dwd, i_meteoswiss, i_arpasim, i_fakeobs, &
       ydate_ini_mod,                           &
       idom,                          &
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       raddeg,              &
       i_fwo_composites, i_fwo_out, &
       rs_meta, rs_data, cart_data, dbz_meta, &
       missing_obs, &
       fdbk_meta_type, &
       comp_meta, comp_meta_bub

  USE radar_composites, ONLY : &
       composite2D_dbz_maxmethod_ista,  &
       comp_dbzsim_tot,     comp_dbzobs_tot, &
       comp_dbzsim_bub_tot, comp_dbzobs_bub_tot

  USE radar_interface, ONLY : &
       abort_run,             &
       get_runtime_timings, &
       check_if_currtime_is_obstime, &
       it_is_time_for_radar, &
#ifdef NUDGING
       get_fdbk_metadata, &
#endif
       get_obstime_ind_of_modtime

  !------------------------------------------------------------------------------

  USE radar_utilities, ONLY :   &
                               ind2sub3D, sub2ind3D,       &
                               phirot2phi,                 &
                               rlarot2rla,                 &
                               phi2phirot,                 &
                               rla2rlarot,                 &
                               init_vari,                  &
                               set_missing_and_correct0

  !------------------------------------------------------------------------------

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, lequal_azi_alldatasets, loutradwind, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       lfdbk_output, &
       lqc_flag, ldealiase_vr_obs, &
       itype_supobing, supob_nrb, &
       supob_vrw, supob_rfl, &
       itype_obserr_vr, baseval_obserr_vr, maxval_obserr_vr, ramp_lowdbz_obserr_vr, ramp_highdbz_obserr_vr, &
       baseval_obserr_dbz, &
       itype_metric_refl_fdbk, minval_obserr_lwc, lmds_z, lmds_vr, &
       supob_lowthresh_z_obs, supob_lowthresh_z_sim, &
       ldo_composite, lcomposite_output, &
       nel_composite, ldo_bubbles, &
       htop      ! Maximum height MSL for radar computations ( <= model domain top height )

#ifdef NETCDF
  USE radar_obs_data_read, ONLY : read_obs_rad
#endif

  USE radar_model2rays, ONLY : get_azvec

  USE radar_data_io, ONLY : radgeomoutputunit, radwindoutputunit, &
       &                    radwindobsoutputunit, radrefloutputunit, &
       &                    radreflobsoutputunit, & !extrefloutputunit, &
       &                    zdroutputunit, zdrobsoutputunit, &
       &                    rhvoutputunit, rhvobsoutputunit, &
       &                    kdpoutputunit, kdpobsoutputunit, &
       &                    ahoutputunit, adpoutputunit, &
       &                    ldroutputunit, ldrobsoutputunit

  USE radar_output_utils, ONLY : get_fileprefix_ascii_output, control_output

  USE radar_output3d_ascii, ONLY : output3d_ascii_radar

#if defined(NUDGING) && defined(NETCDF)
  USE radar_output_fdbk, ONLY : obs_write_cdf_feedobs
#endif


  !------------------------------------------------------------------------------

!================================================================================
!================================================================================

  IMPLICIT NONE

!================================================================================
!================================================================================

  PRIVATE

  PUBLIC :: output_radar_country_obs

  !==============================================================================

  !.. Various small epsilon-thresholds:
  REAL (KIND=dp), PARAMETER :: eps_vr    = 1e-12_dp  ! |radial winds| < eps_vr are set to 0.0_dp

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !==============================================================================
  !+ Module procedure in radar_src for:
  !  - output of obs data of one radar station
  !  - output of fdbk files (sim and obs) including superobing
  !  on a specific PE. Is called from output_my_ista() or output_my_ista_smth()
  !  within the parallelization loops.
  !------------------------------------------------------------------------------

  SUBROUTINE output_radar_country_obs ( ista, time_mod, itime, vr_mod, &
                                        zr_mod, zdr_mod, kdp_mod, phidp_mod, rhv_mod, &
                                        hr_mod, lat_mod, lon_mod, &
                                        vr_mod_for_dealiasing )

    !------------------------------------------------------------------------------
    !
    ! Description: Reads the radar obs data of station ista from NetCDF files
    !              (one file per data type - vr, qv, z, qz - for the entire time frame)
    !              and writes the data to disc. This is done in
    !              ascii files, netcdf files or bin files,
    !              as well as NetCDF feedback files.
    !
    !              Before writing, the radar data are sorted into a 3D array
    !              according to range, azimut and elevation indices,
    !              data(range,azi,ele),
    !              range runs from 1 ... Nrange ( r = range * ra_inc )
    !              azi   runs from 1 ... Naz    ( az = az_start + (azi-1)*az_inc )
    !              ele   runs from 1 ... Nele   ( as defined in el_arr(ele) meta data )
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    REAL(KIND=dp), INTENT(IN) :: time_mod     ! seconds since model start
    INTEGER,       INTENT(IN) :: ista         ! index of radar station
    INTEGER,       INTENT(in) :: itime
    
    ! simulated radar fields for feedback files with TARGET attribute
    REAL(KIND=dp), INTENT(IN), TARGET, &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
              vr_mod, lat_mod, lon_mod, &
              zr_mod, zdr_mod, kdp_mod, phidp_mod, rhv_mod

    ! simulated radar fields for feedback files without TARGET attribute:
    REAL(KIND=dp), INTENT(IN), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
              hr_mod
    REAL(KIND=dp), OPTIONAL, DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         vr_mod_for_dealiasing                          ! radial wind for dealising purposes,
                                                        ! which means that we need a full volume coverage
                                                        ! without "empty" regions, caused by, e.g., min. detect. signal

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'output_radar_country_obs'

    INTEGER    :: i, iobs, nobs_obs, tmppos_obs, j, k, tot_num_observables_fdbk

    INTEGER, ALLOCATABLE :: &
         m_obs(:),  & ! azimuthal indices for all radar points of a radar station in observations
         n_obs(:),  & ! radial    indices for all radar points of a radar station in observations
         o_obs(:)     ! elevation indices for all radar points of a radar station in observations

    REAL(KIND=dp), ALLOCATABLE, TARGET :: vrpolar_obs(:,:,:)    , zrpolar_obs(:,:,:)    , &
                                          vrpolar_obs_q(:,:,:)  , zrpolar_obs_q(:,:,:)  , &
                                          vrpolar_obs_err(:,:,:), zrpolar_obs_err(:,:,:), &
                                          vrobs_supob(:,:,:)    , vrmod_supob(:,:,:)    , &
                                          zrobs_supob(:,:,:)    , zrmod_supob(:,:,:)    , zrmod_tmp(:,:,:), &
                                          vrobserr_supob(:,:,:) , zrobserr_supob(:,:,:) , &
                                          lat_supob(:,:,:)      , lon_supob(:,:,:)

    REAL(KIND=dp), ALLOCATABLE, TARGET :: zdrpolar_obs(:,:,:), zdrobs_supob(:,:,:), zdrmod_supob(:,:,:), &
                                          kdppolar_obs(:,:,:), kdpobs_supob(:,:,:), kdpmod_supob(:,:,:), &
                                          phidppolar_obs(:,:,:), phidpobs_supob(:,:,:), phidpmod_supob(:,:,:), &
                                          rhvpolar_obs(:,:,:), rhvobs_supob(:,:,:), rhvmod_supob(:,:,:)

    REAL(KIND=dp), POINTER :: p_polar_obs(:,:,:), p_polar_obs_q(:,:,:), &
                              p_polar_obserr(:,:,:), &
                              p_polar_obs_supob(:,:,:), &
                              p_polar_mod(:,:,:), p_polar_mod_q(:,:,:), &
                              p_lat_mod(:,:,:)  , p_lon_mod(:,:,:)

    CHARACTER(len=100) :: fileprefix

    LOGICAL            :: zlradwind_ok, zlzradar_ok, zlzradar_avail, zlrefend, zlastpr, zlwrite_veridata
    LOGICAL            :: it_is_time_for_feedback, zlread_dbz_obs
    LOGICAL            :: zlzdrradar_ok, zlkdpradar_ok, zlphidpradar_ok, zlrhvradar_ok

#ifdef NUDGING
    CHARACTER (LEN=300)        :: yerrstring
    INTEGER                    :: fdbk_err, nexce_rep, nexce_bdy
    REAL(KIND=dp), ALLOCATABLE ::   &
         rpress(:,:,:), zazivec(:), zrvec(:)
    REAL(KIND=dp)              :: t11000, p11000, exp_poly
    TYPE(fdbk_meta_type)       :: fdbk_meta_data
#endif

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE output_radar_country_obs
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    !.. Read obs radar data from station number "ista" at time "time_mod":

    SELECT CASE (rs_meta(ista)%icountry)
    CASE ( 1 : ncountry_max )

#ifndef NETCDF
      IF (rs_meta(ista)%icountry == i_dwd) THEN
        CALL abort_run (my_radar_id, 15075, &
             'ERROR: problem in '//TRIM(yzroutine)//': Reading of NetCDF radar data not possible, missing flag -DNETCDF!', &
             'radar_process_output.f90, '//TRIM(yzroutine)//'()')
      END IF
#endif
#ifndef HDF5_RADAR_INPUT
      IF (rs_meta(ista)%icountry == i_dwd .OR. rs_meta(ista)%icountry == i_meteoswiss .OR. rs_meta(ista)%icountry == i_arpasim) THEN
        CALL abort_run (my_radar_id, 15075, &
             'ERROR: problem in '//TRIM(yzroutine)//': Reading of HDF5 radar data not possible, missing flag -DHDF5_RADAR_INPUT!', &
             'radar_process_output.f90, '//TRIM(yzroutine)//'()')
      END IF
#endif

#ifdef NETCDF
      CALL read_obs_rad(ista, time_mod, 'read')
#endif

    CASE ( i_fakeobs )

      ! For testsuite purposes: do not actually read any data,
      ! but fake obs data below. To enable this, some parameters
      ! have to be set properly here:
      rs_data(ista)%nobs_obs(:)      = 0
      rs_meta(ista)%lobs_avail(:,:)   = .TRUE.

    CASE default
      CALL abort_run (my_radar_id, 15077, &
           'ERROR: problem in '//TRIM(yzroutine)//                                                      &
           ': Wrong icountry! the only implemented options are 1 (DWD), 2 (MeteoSwiss), and 3 (Itlay - Emilia Romana)!', &
           'radar_process_output.f90, '//TRIM(yzroutine)//'()')
    END SELECT

    ! .. nobs_obs from read_XXX_rad:
    nobs_obs = MAXVAL(rs_data(ista)%nobs_obs)

    ! .. check which observations (and quality flags) are available for output:
    zlradwind_ok   = ( rs_meta(ista)%lobs_avail(itime,i_vrad) .AND. rs_meta(ista)%lobs_avail(itime,i_qualvrad) ) .AND. loutradwind
    zlzradar_avail = ( rs_meta(ista)%lobs_avail(itime,i_dbzh) .AND. rs_meta(ista)%lobs_avail(itime,i_qualdbzh) )
    zlzradar_ok    = zlzradar_avail .AND. loutdbz
!!$ UB: no quality flags used for ZDR yet. Could use refl. quality flag in the future, but needs to be discussed!
    zlzdrradar_ok = rs_meta(ista)%lobs_avail(itime,i_zdr) .AND. (loutpolstd .OR. loutpolall)
    zlkdpradar_ok = rs_meta(ista)%lobs_avail(itime,i_kdp) .AND. (loutpolstd .OR. loutpolall)
    zlphidpradar_ok = rs_meta(ista)%lobs_avail(itime,i_phidp) .AND. (loutpolstd .OR. loutpolall)
    zlrhvradar_ok = rs_meta(ista)%lobs_avail(itime,i_rhv) .AND. (loutpolstd .OR. loutpolall)

#ifdef NUDGING
    ! .. check if it is time for feedback output:
    it_is_time_for_feedback = it_is_time_for_radar ( rs_meta(ista)%obs_times_fdbk )
    ! .. check if reflectivity obs have been read and are needed:
    zlread_dbz_obs = ( (zlzradar_ok) .OR. &
         (it_is_time_for_feedback .AND. zlradwind_ok .AND. itype_obserr_vr > 0) )
#else
    it_is_time_for_feedback = .TRUE.   ! .TRUE. will trigger an appropriate error message below
    zlread_dbz_obs = zlzradar_ok
#endif

    ! .. determine the total number of observables that will go to the feedback files:
    tot_num_observables_fdbk = 0
    IF (lfdbk_output .AND. it_is_time_for_feedback .AND. loutradwind .AND. rs_meta(ista)%lvrad_to_fdbk) THEN
      ! Check if there are any fdbk output times with valid radial wind obs. If yes, increment counter for number of observables:
      DO i=1, rs_meta(ista)%nobs_times_fdbk
        k = get_obstime_ind_of_modtime ( rs_meta(ista)%obs_times_fdbk(i) , &
                                         rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) )
        IF (k > missthr_int) THEN
          IF (rs_meta(ista)%obs_startrec(k,i_vrad) > missthr_int .AND. &
               ( .NOT.lqc_flag .OR. ( lqc_flag .AND. rs_meta(ista)%obs_startrec(k,i_qualvrad) > missthr_int ) ) ) THEN
            tot_num_observables_fdbk = tot_num_observables_fdbk + 1
            EXIT
          END IF
        END IF
      END DO
    END IF
    IF (lfdbk_output .AND. it_is_time_for_feedback .AND. loutdbz .AND. rs_meta(ista)%ldbzh_to_fdbk) THEN
      ! Check if there are any fdbk output times with valid reflectivity obs. If yes, increment counter for number of observables:
      DO i=1, rs_meta(ista)%nobs_times_fdbk
        k = get_obstime_ind_of_modtime ( rs_meta(ista)%obs_times_fdbk(i) , &
                                         rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) )
        IF (k > missthr_int) THEN
          IF (rs_meta(ista)%obs_startrec(k,i_dbzh) > missthr_int .AND. &
               ( .NOT.lqc_flag .OR. ( lqc_flag .AND. rs_meta(ista)%obs_startrec(k,i_qualdbzh) > missthr_int ) ) ) THEN
            tot_num_observables_fdbk = tot_num_observables_fdbk + 1
            EXIT
          END IF
        END IF
      END DO
    END IF

    ! .. If radwind obs or reflectivity obs are available, allocate local arrays
    !     and get 3D indices for sorting:
    !----------------------------------------------------------------------------

    IF ( zlradwind_ok .OR. zlzradar_ok) THEN

      ALLOCATE( vrpolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
                vrpolar_obs_q(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
                zrpolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
                zrpolar_obs_q(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
      !.. Initialize with standard missing value:
      CALL init_vari(vrpolar_obs  , miss_value)
      CALL init_vari(vrpolar_obs_q, miss_value)
      CALL init_vari(zrpolar_obs  , miss_value)
      CALL init_vari(zrpolar_obs_q, miss_value)

      IF (itype_supobing > 0) THEN
        ALLOCATE( vrobs_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  vrobserr_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  vrmod_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  zrobs_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  zrobserr_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  zrmod_supob(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  lat_supob  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  lon_supob  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel)   )
        !.. Initialize with standard missing value:
        CALL init_vari(vrobs_supob, miss_value)
        CALL init_vari(vrobserr_supob, baseval_obserr_vr )
        CALL init_vari(vrmod_supob, miss_value)
        CALL init_vari(zrobs_supob, miss_value)
        CALL init_vari(zrobserr_supob, baseval_obserr_dbz )
        CALL init_vari(zrmod_supob, miss_value)
        CALL init_vari(lat_supob  , miss_value)
        CALL init_vari(lon_supob  , miss_value)
      END IF

!!$ fixed      IF ((lfdbk_output .AND. it_is_time_for_feedback) .OR. itype_obserr_vr > 0) THEN
      IF ((lfdbk_output .AND. it_is_time_for_feedback) .OR. itype_supobing > 0) THEN
        ALLOCATE( zrpolar_obs_err(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel),  &
                  vrpolar_obs_err(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
        CALL init_vari(zrpolar_obs_err, baseval_obserr_dbz)  ! Standard value
        CALL init_vari(vrpolar_obs_err, baseval_obserr_vr )  ! Standard value
      END IF

    END IF

    IF ( zlzdrradar_ok ) THEN
      ALLOCATE( zdrpolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
      !.. Initialize with standard missing value:
      CALL init_vari(zdrpolar_obs, miss_value)
    END IF
    IF ( zlkdpradar_ok ) THEN
      ALLOCATE( kdppolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
      !.. Initialize with standard missing value:
      CALL init_vari(kdppolar_obs, miss_value)
    END IF
    IF ( zlphidpradar_ok ) THEN
      ALLOCATE( phidppolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
      !.. Initialize with standard missing value:
      CALL init_vari(phidppolar_obs, miss_value)
    END IF
    IF ( zlrhvradar_ok ) THEN
      ALLOCATE( rhvpolar_obs  (rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
      !.. Initialize with standard missing value:
      CALL init_vari(rhvpolar_obs, miss_value)
    END IF

    IF ( zlradwind_ok .OR. zlzradar_ok .OR. zlzdrradar_ok) THEN

      ALLOCATE( m_obs(nobs_obs), n_obs(nobs_obs), o_obs(nobs_obs) )

      IF ( lequal_azi_alldatasets ) THEN
        CALL init_vari(m_obs, -1)
        CALL init_vari(n_obs, -1)
        CALL init_vari(o_obs, -1)

!$omp parallel do private(iobs, tmppos_obs)
        DO iobs = 1, nobs_obs
          tmppos_obs = rs_data(ista)%ind_intp_obs(iobs)
          IF (tmppos_obs > missthr_int) THEN
            ! .. existing datum in NetCDF file (the data itself might still be missing values):
            CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                 m_obs(iobs), n_obs(iobs), o_obs(iobs))
          ENDIF
        END DO
!$omp end parallel do
      END IF

      !----------------------------------------------------------------------------------
      ! .. preparations for feedback file output:
      !     coordinate vectors and temporary storage fiels in order (range, azi, ele):

#if defined(NUDGING) && defined(NETCDF)
      IF ( lfdbk_output .AND. it_is_time_for_feedback) THEN

        zlrefend = .FALSE.

        ! .. zlastpr has to be .true. if we are at the last output time step:
        zlastpr  = check_if_currtime_is_obstime ( rs_meta(ista)%obs_times_fdbk(rs_meta(ista)%nobs_times_fdbk) )

        zlwrite_veridata = .TRUE.   ! write also the simulated data to feedback file


        CALL get_azvec (rs_meta(ista), zazivec)

        ALLOCATE (zrvec(rs_meta(ista)%nra))
        DO i=1, rs_meta(ista)%nra
          zrvec(i) = i * rs_meta(ista)%ra_inc
        END DO

        ! .. pressure at radar bin heights according to ICAO standard atmosphere:
 !!$ Caution: data fields in order (azi, range, ele), not (range, azi, ele)!
        ALLOCATE ( rpress(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel) )
        CALL init_vari(rpress, miss_value)

        exp_poly = 9.816_dp / (287.05_dp*0.0065_dp)
        t11000   = 293.16_dp - 0.0065_dp * 11000.0_dp
        p11000   = 1013.25e2_dp * EXP(exp_poly*LOG(t11000/293.16_dp))
!$omp parallel do private(i,j,k)
        DO k=1, rs_meta(ista)%nel
          DO j=1, rs_meta(ista)%nra
            DO i=1, rs_meta(ista)%naz
              IF (hr_mod(i,j,k) <= 11000.0_dp .AND. hr_mod (i,j,k)>= miss_threshold) THEN
                rpress(i,j,k) = 1013.25e2_dp * EXP(exp_poly*LOG(1.0_dp - 0.0065_dp*hr_mod(i,j,k)/293.16_dp))
              ELSE IF (hr_mod(i,j,k) > 11000.0_dp .AND. hr_mod(i,j,k) <= htop) THEN
                rpress(i,j,k) =  p11000 * EXP(-9.816_dp * (hr_mod(i,j,k)-11000.0_dp) / (287.05_dp*t11000))
              END IF
            END DO
          END DO
        END DO
!$omp end parallel do
!!$        WHERE (hr_mod <= 11000.0_dp .AND. hr_mod >= miss_threshold)
!!$          rpress = 1013.25e2_dp * (1.0_dp - 0.0065_dp*hr_mod/293.16_dp) ** exp_poly
!!$        ELSEWHERE (hr_mod > 11000.0_dp .AND. hr_mod <= htop)
!!$          rpress =  p11000 * EXP(-9.816_dp * (hr_mod-11000.0_dp) / (287.05_dp*t11000))
!!$        END WHERE

      END IF
#endif

    END IF

    !------------------------------------------------------------------------------------------------
    ! .. If reflectivity obs are available and needed, read them in, sort them into regular 3D arrays
    !     and apply quality flags:
    !------------------------------------------------------------------------------------------------

    IF ( zlread_dbz_obs ) THEN

      IF ( .NOT.lequal_azi_alldatasets ) THEN

        CALL init_vari(m_obs, -1)
        CALL init_vari(n_obs, -1)
        CALL init_vari(o_obs, -1)
!$omp parallel do private(iobs, tmppos_obs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          tmppos_obs = rs_data(ista)%ind_intp_obs_z(iobs)
          IF (tmppos_obs > missthr_int) THEN
            ! .. existing datum in NetCDF file (the data itself might still be missing values):
            CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                 m_obs(iobs), n_obs(iobs), o_obs(iobs))
          ENDIF
        END DO
!$omp end parallel do
      END IF

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!$omp parallel do private(iobs)
      DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
        IF ( rs_data(ista)%ind_intp_obs_z(iobs) > missthr_int .AND. &
             rs_data(ista)%zh_radar_obs(iobs) /= missing_obs ) THEN
          ! .. Observed reflectivity is already in logarithmic dBZ units:
          zrpolar_obs( m_obs(iobs),&
                       n_obs(iobs),&
                       o_obs(iobs) ) = rs_data(ista)%zh_radar_obs(iobs)
        END IF
      END DO
!$omp end parallel do

      IF (lqc_flag) THEN

        IF ( .NOT.lequal_azi_alldatasets ) THEN

          CALL init_vari(m_obs, -1)
          CALL init_vari(n_obs, -1)
          CALL init_vari(o_obs, -1)
!$omp parallel do private(iobs, tmppos_obs)
          DO iobs = 1, rs_data(ista)%nobs_obs(i_qualdbzh)
            tmppos_obs = rs_data(ista)%ind_intp_obs_qz(iobs)
            IF (tmppos_obs > missthr_int) THEN
              ! .. existing datum in NetCDF file (the data itself might still be missing values):
              CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                   m_obs(iobs), n_obs(iobs), o_obs(iobs))
            ENDIF
          END DO
!$omp end parallel do

        END IF

        !.. Sort the collected radar points into the 3D field according
        !   to sorted azimut, range and elevation:
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_qualdbzh)
          IF ( rs_data(ista)%ind_intp_obs_qz(iobs) > missthr_int .AND. &
               rs_data(ista)%zh_radar_obs_q(iobs) /= missing_obs ) THEN
            ! .. Observed reflectivity is already in logarithmic dBZ units:
            zrpolar_obs_q( m_obs(iobs),&
                           n_obs(iobs),&
                           o_obs(iobs) ) = rs_data(ista)%zh_radar_obs_q(iobs)
          END IF
        END DO
!$omp end parallel do

      END IF  !  lqc_flag

!!$ Special faking of obs data for testing purposes or for testsuite mode,
!!$  if the actual obs data reader does not yet exist:
      IF (rs_meta(ista)%icountry == i_fakeobs) THEN
        WRITE (*,*) '*** WARNING: Experimental run with faked reflectivity obs! ***'
        zrpolar_obs  (:,:,:) = zr_mod(:,:,:)
        zrpolar_obs_q(:,:,:) = 0
      END IF
!!$ End of special faking!

      ! ..  Eliminate values above model top and for blocked ray parts:
!$omp parallel do private(i,j,k) collapse(3)
      DO k=1, rs_meta(ista)%nel
        DO j=1, rs_meta(ista)%nra
          DO i=1, rs_meta(ista)%naz
            IF (hr_mod(i,j,k) > htop .OR. hr_mod(i,j,k) < miss_threshold) THEN
              zrpolar_obs(i,j,k)   = miss_value
              zrpolar_obs_q(i,j,k) = miss_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do

      ! Set correct 0-value (zero_value) for zrpolar (missing values have been set above (miss_value)):
      call set_missing_and_correct0(zrpolar_obs)

      ! quality control flag evaluation for the z observations:
      IF (lqc_flag) THEN
        CALL qcrad_flag_eval(ista, 'z', zrpolar_obs_q, zrpolar_obs)
      END IF

      IF (ldebug_radsim) THEN
        CALL control_output(radreflobsoutputunit(ista), time_mod, rs_meta(ista), "Z [dBZ]", &
           zrpolar_obs, (zrpolar_obs > miss_threshold), miss_value, nobs_obs)
      END IF

    END IF  ! zlread_dbz_obs

    !---------------------------------------------------------------------
    ! .. If radwind obs are available, sort them into regular 3D data set:
    !---------------------------------------------------------------------

    IF ( zlradwind_ok ) THEN

      IF ( .NOT.lequal_azi_alldatasets ) THEN

        CALL init_vari(m_obs, -1)
        CALL init_vari(n_obs, -1)
        CALL init_vari(o_obs, -1)
!$omp parallel do private(iobs, tmppos_obs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_vrad)
          tmppos_obs = rs_data(ista)%ind_intp_obs_vr(iobs)
          IF (tmppos_obs > missthr_int) THEN
            ! .. existing datum in NetCDF file (the data itself might still be missing values):
            CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                 m_obs(iobs), n_obs(iobs), o_obs(iobs))
          ENDIF
        END DO
!$omp end parallel do
      END IF

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!$omp parallel do private(iobs)
      DO iobs = 1, rs_data(ista)%nobs_obs(i_vrad)
        IF ( rs_data(ista)%ind_intp_obs_vr(iobs) > missthr_int.AND. &
             rs_data(ista)%radwind_obs(iobs) /= missing_obs ) THEN
          vrpolar_obs( m_obs(iobs),&
                       n_obs(iobs),&
                       o_obs(iobs) ) = rs_data(ista)%radwind_obs(iobs)
        END IF
      END DO
!$omp end parallel do

      IF (lqc_flag) THEN

        IF ( .NOT.lequal_azi_alldatasets ) THEN

          CALL init_vari(m_obs, -1)
          CALL init_vari(n_obs, -1)
          CALL init_vari(o_obs, -1)
!$omp parallel do private(iobs, tmppos_obs)
          DO iobs = 1, rs_data(ista)%nobs_obs(i_qualvrad)
            tmppos_obs = rs_data(ista)%ind_intp_obs_qv(iobs)
            IF (tmppos_obs > missthr_int) THEN
            ! .. existing datum in NetCDF file (the data itself might still be missing values):
              CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                   m_obs(iobs), n_obs(iobs), o_obs(iobs))
            ENDIF
          END DO
!$omp end parallel do
        END IF

        !.. Sort the collected radar points into the 3D field according
        !   to sorted azimut, range and elevation:
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_qualvrad)
          IF ( rs_data(ista)%ind_intp_obs_qv(iobs) > missthr_int .AND. &
               rs_data(ista)%radwind_obs_q(iobs) /= missing_obs ) THEN
            vrpolar_obs_q( m_obs(iobs),&
                           n_obs(iobs),&
                           o_obs(iobs) ) = rs_data(ista)%radwind_obs_q(iobs)
          END IF
        END DO
!$omp end parallel do

      END IF  ! lqc_flag

!!$ Special faking of obs data for testing purpuses, when the actual obs
!!$  data reader does not yet exist:
      IF (rs_meta(ista)%icountry == i_fakeobs) THEN
        WRITE (*,*) '*** WARNING: Experimental run with faked radial wind obs! ***'
        vrpolar_obs  (:,:,:) = vr_mod(:,:,:)
        vrpolar_obs_q(:,:,:) = 0
      END IF
!!$ End of special faking!

      ! ..  Eliminate values above model top:
!!$      WHERE (hr_mod > htop)
!!$        vrpolar_obs   = miss_value
!!$        vrpolar_obs_q = miss_value
!!$      END WHERE
!$omp parallel do private(i,j,k) collapse(3)
      DO k=1, rs_meta(ista)%nel
        DO j=1, rs_meta(ista)%nra
          DO i=1, rs_meta(ista)%naz
            IF (hr_mod(i,j,k) > htop .OR. hr_mod(i,j,k) < miss_threshold) THEN
              vrpolar_obs(i,j,k)   = miss_value
              vrpolar_obs_q(i,j,k) = miss_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do

      ! quality control flag evaluation for the vr observations:
      IF (lqc_flag) THEN
        CALL qcrad_flag_eval(ista, 'vr', vrpolar_obs_q, vrpolar_obs)
      END IF

      ! dealiasing radial wind for the vr observations:
      IF (ldealiase_vr_obs) THEN
        IF (PRESENT(vr_mod_for_dealiasing)) THEN
          CALL dealiasing_vr(ista, itime, vr_mod_for_dealiasing, vrpolar_obs)
        ELSE
          CALL abort_run (my_radar_id, 15089, &
               'ERROR in '//TRIM(yzroutine)//': OPTIONAL argument vr_mod_for_dealiasing not present, '//&
               'but is needed for dealiasing_vr()', &
               'radar_process_output.f90, '//TRIM(yzroutine)//'()')
        END IF
      END IF

    END IF


    !---------------------------------------------------------------------
    ! .. If radwind obs are available, do different kinds of output:
    !---------------------------------------------------------------------

    IF ( zlradwind_ok ) THEN

      !.. Output the sorted polar data set into a standard file format:
      !   - binary on the SX9, convert to simple ASCII with program "bin2ascii_convrates3d" from Ulrich Blahak
      !   - simple ASCII on all other systems
      CALL get_fileprefix_ascii_output ('vrobs', rs_meta(ista), fileprefix)
      CALL output3d_ascii_radar(itime, &
           rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
           vrpolar_obs, TRIM(ADJUSTL(fileprefix)), &
           'Observed radial velocity, no aliasing', 'm/s', 'polar', &
            rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.NOT.ldealiase_vr_obs)

      IF (lqc_flag) THEN
        CALL get_fileprefix_ascii_output ('qvobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrpolar_obs_q, TRIM(ADJUSTL(fileprefix)), &
             'Quality of observed radial velocity, no aliasing', '-', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.NOT.ldealiase_vr_obs)
      END IF

      IF (ldebug_radsim) THEN
        CALL control_output(radwindobsoutputunit(ista), time_mod, rs_meta(ista), "V_r [m/s]", &
             vrpolar_obs, (vrpolar_obs > miss_threshold), miss_value, nobs_obs)
      END IF
      
      SELECT CASE (itype_supobing)
      CASE (0)
      CASE (1)
        IF (itype_obserr_vr == 2) THEN
          ! Compute vrpolar_obs_err and vrobs_supob in any case for diagnostic output of the superob'ed quantities:
          CALL set_obserr_radwind( ista=ista, zrdata=zrpolar_obs, missingthresh=miss_threshold, &
                                   dBZ_0=ramp_lowdbz_obserr_vr, e_o_0=maxval_obserr_vr, &
                                   dBZ_1=ramp_highdbz_obserr_vr, e_o_1=baseval_obserr_vr, &
                                   rdata_obserr=vrpolar_obs_err )
        ENDIF
        CALL averagesuperobing_radar   (ista, 'vr', lat_mod, lon_mod, vrpolar_obs, vr_mod, vrpolar_obs_err, vrobs_supob, &
                                        vrmod_supob, vrobserr_supob, lat_supob, lon_supob)
        IF (itype_obserr_vr == 1 .AND. lfdbk_output .AND. it_is_time_for_feedback .AND. &
             zlzradar_avail .AND. rs_meta(ista)%lvrad_to_fdbk) THEN
          ! This is only needed for actual output of feedback files:
          ! Extra superobing of un-clipped reflectivity. Only zrobs_supob will be used
          !  below for vr observation error specification. zrmod_supob is a dummy.
          CALL averagesuperobing_radar(ista, 'z', lat_mod, lon_mod, zrpolar_obs, zr_mod, zrpolar_obs_err, zrobs_supob, &
                                       zrmod_supob, zrobserr_supob, lat_supob, lon_supob )
        END IF
      CASE (2)
        IF (itype_obserr_vr == 2) THEN
          ! Compute vrpolar_obs_err and vrobs_supob in any case for diagnostic output of the superob'ed quantities:
          CALL set_obserr_radwind( ista=ista, zrdata=zrpolar_obs, missingthresh=miss_threshold, &
                                   dBZ_0=ramp_lowdbz_obserr_vr, e_o_0=maxval_obserr_vr, &
                                   dBZ_1=ramp_highdbz_obserr_vr, e_o_1=baseval_obserr_vr, &
                                   rdata_obserr=vrpolar_obs_err )
        ENDIF
        CALL mediansuperobing_radar_vec(ista, 'vr', lat_mod, lon_mod, vrpolar_obs, vr_mod, vrpolar_obs_err, vrobs_supob, &
                                        vrmod_supob, lat_supob, lon_supob)
        IF (itype_obserr_vr == 1 .AND. lfdbk_output .AND. it_is_time_for_feedback .AND. &
             zlzradar_avail .AND. rs_meta(ista)%lvrad_to_fdbk) THEN
          ! This is only needed for actual output of feedback files:
          ! Extra superobing of un-clipped reflectivity. Only zrobs_supob will be used
          !  below for vr observation error specification. zrmod_supob is a dummy.
          CALL mediansuperobing_radar_vec(ista, 'z', lat_mod, lon_mod, zrpolar_obs, zr_mod, zrpolar_obs_err, zrobs_supob, &
                                          zrmod_supob, lat_supob, lon_supob )
        END IF
      CASE default
        CALL abort_run (my_radar_id, 15079, &
             'ERROR: problem in '//TRIM(yzroutine)//': Wrong itype_supobing! Possible values are 0, 1,or 2', &
             'radar_process_output.f90, '//TRIM(yzroutine)//'()')
      END SELECT

      IF (itype_supobing > 0) THEN

        CALL get_fileprefix_ascii_output ('vrsupobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrobs_supob, TRIM(ADJUSTL(fileprefix)), &
             'Observed radial velocity, no aliasing', 'm/s', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.NOT.ldealiase_vr_obs)
          
        CALL get_fileprefix_ascii_output ('vrsupsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrmod_supob, TRIM(ADJUSTL(fileprefix)), &
             'Simulated radial velocity, no aliasing', 'm/s', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          
        CALL get_fileprefix_ascii_output ('losupsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             lon_supob, TRIM(ADJUSTL(fileprefix)), &
             'Longitude of superobservation grid points', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          
        CALL get_fileprefix_ascii_output ('lasupsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             lat_supob, TRIM(ADJUSTL(fileprefix)), &
             'Latitude of superobservation grid points', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

      ENDIF
    
      !----------------------------------------------------------------------------
      ! .. Output of NetCDF feedback file (vr and z are written to the same file!):
      !----------------------------------------------------------------------------

      IF ( lfdbk_output .AND. it_is_time_for_feedback .AND. rs_meta(ista)%lvrad_to_fdbk) THEN

#if defined(NUDGING) && defined(NETCDF)

        yerrstring(:) = ' '
        fdbk_err = 0

        ! .. computation of radial wind obserr:
        IF (ABS(baseval_obserr_vr-1.0_dp) > 1e-6_dp) THEN
          WRITE (*,'(a,f0.1,a)') 'WARNING '//TRIM(yzroutine)// &
               ': baseval_obserr_vr for radial wind has been set by namelist to ', &
               baseval_obserr_vr, ', not to its desired relative value 1.0!'
        END IF

        SELECT CASE (itype_obserr_vr)
        CASE (0)
          CALL init_vari(vrpolar_obs_err, baseval_obserr_vr)
        CASE (1)
          IF (zlzradar_avail) THEN
            IF (itype_supobing == 0) THEN
              p_polar_obs_supob => zrpolar_obs
            ELSE
              p_polar_obs_supob => zrobs_supob
            END IF
            ! if baseval_obserr_vr = 1.0, vrpolar_obs_err is the relative error increase which can be superimposed on some
            !  other obs err settings in DACE:
            CALL set_obserr_radwind( ista=ista, zrdata=p_polar_obs_supob, missingthresh=miss_threshold, &
                                     dBZ_0=ramp_lowdbz_obserr_vr, e_o_0=maxval_obserr_vr, &
                                     dBZ_1=ramp_highdbz_obserr_vr, e_o_1=baseval_obserr_vr, &
                                     rdata_obserr=vrpolar_obs_err )
          ELSE
            WRITE (*,'(a,/,a,i6.6,a,/,a,f0.1,a)') 'WARNING '//TRIM(yzroutine)// &
                 ' for reflectivity dependent obs error for radial wind in fof-files:', &
                 '   MISSING reflectivity obs data for station=', rs_meta(ista)%station_id, &
                 ', scanname='//TRIM(rs_meta(ista)%scanname)//' and date/time='// &
                 TRIM(rs_meta(ista)%obs_cdate(itime)), &
                 '   Re-setting radial wind obs error to constant max value ',maxval_obserr_vr,'!'
            CALL init_vari(vrpolar_obs_err, maxval_obserr_vr)
          END IF
        CASE (2)
          IF (zlzradar_avail) THEN
            SELECT CASE (itype_supobing)
            CASE (0)
              CALL set_obserr_radwind( ista=ista, zrdata=zrpolar_obs, missingthresh=miss_threshold, &
                                       dBZ_0=ramp_lowdbz_obserr_vr, e_o_0=maxval_obserr_vr, &
                                       dBZ_1=ramp_highdbz_obserr_vr, e_o_1=baseval_obserr_vr, &
                                       rdata_obserr=vrpolar_obs_err )
            CASE (1,2)
              ! vrpolar_obs_err or vrobserr_supob has already been set in above call averagesuperobing_radar()
              CONTINUE
            END SELECT
          ELSE
            WRITE (*,'(a,/,a,i6.6,a,/,a,f0.1,a)') 'WARNING '//TRIM(yzroutine)// &
                 ' for reflectivity dependent obs error for radial wind in fof-files:', &
                 '   MISSING reflectivity obs data for station=', rs_meta(ista)%station_id, &
                 ', scanname='//TRIM(rs_meta(ista)%scanname)//' and date/time='// &
                 TRIM(rs_meta(ista)%obs_cdate(itime)), &
                 '   Re-setting radial wind obs error to constant max value ',maxval_obserr_vr,'!'
            SELECT CASE (itype_supobing)
            CASE (0,1)
              CALL init_vari(vrpolar_obs_err, maxval_obserr_vr)
            CASE (2)
              CALL init_vari(vrobserr_supob, maxval_obserr_vr)
            END SELECT
          END IF
        CASE default
          CALL abort_run (my_radar_id, 15057, &
               'ERROR: problem in '//TRIM(yzroutine)//': Wrong itype_obserr_vr! Possible values are 0, 1 or 2', &
               'radar_process_output.f90, '//TRIM(yzroutine)//'()')
        END SELECT

        ! Link correct input fields for obs_write_cdf_feedobs depending on superobing:
        SELECT CASE (itype_supobing)
        CASE (0)
          p_polar_obs => vrpolar_obs
          p_polar_obserr => vrpolar_obs_err
          p_polar_mod => vr_mod
          p_lat_mod   => lat_mod
          p_lon_mod   => lon_mod
        CASE (1,2)
          p_polar_obs => vrobs_supob
          SELECT CASE (itype_obserr_vr)
          CASE(0,1)
            p_polar_obserr => vrpolar_obs_err
          CASE (2)
            p_polar_obserr => vrobserr_supob
          END SELECT
          p_polar_mod => vrmod_supob
          p_lat_mod   => lat_supob
          p_lon_mod   => lon_supob
        END SELECT


        CALL get_fileprefix_ascii_output ('vrobserr', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrpolar_obs_err, TRIM(ADJUSTL(fileprefix)), &
             'Observation error for radial wind', 'm/s', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (itype_supobing > 0) THEN
          CALL get_fileprefix_ascii_output ('vrsupobserr', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               p_polar_obserr, TRIM(ADJUSTL(fileprefix)), &
               'Superobe''d Observation error for radial wind', 'm/s', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF

        CALL get_fdbk_metadata ( idom, fdbk_meta_data )

        CALL obs_write_cdf_feedobs ( ista, itime, tot_num_observables_fdbk, 'vr', time_mod/60.0, &
             p_polar_obs, vrpolar_obs_q, p_polar_obserr, &
             p_polar_mod, p_lat_mod, p_lon_mod, hr_mod, rpress, &
             zazivec, zrvec, rs_meta(ista)%el_arr(1:rs_meta(ista)%nel), &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel,   &
             htop, &
             (zlastpr .AND. .NOT.loutdbz), zlrefend, zlwrite_veridata,  &
             ydate_ini_mod,                 &
             fdbk_meta_data,                                   &
             nexce_rep, nexce_bdy, fdbk_err, yerrstring,       &
             obs_is_aliased = .NOT.ldealiase_vr_obs, mod_is_aliased = .FALSE. )

        IF (fdbk_err /= 0) THEN
          WRITE (*,'(a)') 'WARNING: problems in '//TRIM(yzroutine)//' writing vr: '//TRIM(yerrstring)
        END IF

#else

        CALL abort_run (my_radar_id, 15078, &
             'ERROR: problem in '//TRIM(yzroutine)//                                                      &
             ': Writing of NetCDF feedback files not possible, missing flags -DNETCDF and/or -DNUDGING!'//&
             ' In ICON, you have to enable submodule dace_icon!', &
             'radar_process_output.f90, '//TRIM(yzroutine)//'()')

#endif

      END IF  ! lfdbk_output

    END IF  ! zlradwind_ok


    !---------------------------------------------------------------------
    ! .. If reflectivity obs are available, do different kinds of output:
    !--------------------------------------------------------------------

    IF ( zlzradar_ok ) THEN

      !.. Output the sorted polar data set into a standard file format:
      !   - binary on the SX9, convert to simple ASCII with program "bin2ascii_convrates3d" from Ulrich Blahak
      !   - simple ASCII on all other systems
      CALL get_fileprefix_ascii_output ('zrobs', rs_meta(ista), fileprefix)
      CALL output3d_ascii_radar(itime, &
           rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
           zrpolar_obs, TRIM(ADJUSTL(fileprefix)), &
           'Observed radar reflectivity', 'dBZ', 'polar', &
           rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)

      IF (lqc_flag) THEN
        CALL get_fileprefix_ascii_output ('qzobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zrpolar_obs_q, TRIM(ADJUSTL(fileprefix)), &
             'Observed radar reflectivity', 'dBZ', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)
      END IF

      !.. Add the contribution of this radar station to the global sim- and obs- composites
      !    on the model grid:

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_composites)
#endif
      IF (ldo_bubbles .AND. zlzradar_ok) THEN
        IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
          IF (rs_meta(ista)%eleind_for_composite_bub == 98) THEN
            ! .. Construct the composite from the precip scans (not volume scan):
            !     The elevation for lat/lon computations will not be the true elevations
            !     but the nominal one (either 0.4, 0.59, 0.8, or 1.3)
            CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                 comp_meta=comp_meta_bub, &
                 eleind_for_composite=1, compdata_tot=comp_dbzobs_bub_tot, &
                 ldebug=ldebug_radsim)
          END IF
        ELSE
          IF (rs_meta(ista)%eleind_for_composite_bub == 99) THEN
            ! .. If eleind_for_composite_bub = 99, then make a vertical maximum composite of all elevations:
            DO i=1, rs_meta(ista)%nel
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                   comp_meta=comp_meta_bub, &
                   eleind_for_composite=i, compdata_tot=comp_dbzobs_bub_tot, &
                   ldebug=ldebug_radsim)
            END DO
          ELSE  IF (rs_meta(ista)%eleind_for_composite_bub /= 98) THEN
            ! .. Else, construct the composite from one single elevation:
            CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                 comp_meta=comp_meta_bub, &
                 eleind_for_composite=rs_meta(ista)%ind_ele_present( &
                                        MIN( rs_meta(ista)%eleind_for_composite_bub, rs_meta(ista)%nel_present) ), &
                 compdata_tot=comp_dbzobs_bub_tot, &
                 ldebug=ldebug_radsim)
          END IF
        END IF
      END IF

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_composites)
#endif
      IF (ldo_composite .AND. zlzradar_ok) THEN
        DO i = 1, nel_composite
          IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
            IF (rs_meta(ista)%eleindlist_for_composite(i) == 98) THEN
              ! .. Construct the composite from the precip scans (not volume scan):
              !     The elevation for lat/lon computations will not be the true elevations
              !     but the nominal one (either 0.4, 0.59, 0.8, or 1.3)
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                   comp_meta=comp_meta, &
                   eleind_for_composite=1, compdata_tot=comp_dbzobs_tot(:,:,i), &
                   ldebug=ldebug_radsim)
            END IF
          ELSE

            IF (rs_meta(ista)%eleindlist_for_composite(i) == 99) THEN
              ! .. If eleindlist_for_composite(i) = 99, then make a vertical maximum composite of all elevations:
              DO j=1, rs_meta(ista)%nel
                CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                     comp_meta=comp_meta, &
                     eleind_for_composite=j, compdata_tot=comp_dbzobs_tot(:,:,i), &
                     ldebug=ldebug_radsim)
              END DO
            ELSE IF (rs_meta(ista)%eleindlist_for_composite(i) /= 98) THEN
              ! .. Else, construct the composite from one single elevation:
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar_obs, rsm=rs_meta(ista), &
                   comp_meta=comp_meta, &
                   eleind_for_composite=rs_meta(ista)%ind_ele_present( &
                                          MIN( rs_meta(ista)%eleindlist_for_composite(i), rs_meta(ista)%nel_present) ), &
                   compdata_tot=comp_dbzobs_tot(:,:,i), &
                   ldebug=ldebug_radsim)
            END IF
          END IF
        END DO
      END IF
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_out)
#endif


      IF (itype_supobing > 0) THEN
        ! .. Optionally, apply a lower threshold to the data before superobing. This
        !    threshold will be the value for valid no-rain data and has to
        !    be consistent to the similar threshold in the LETKF-software:
        IF (supob_lowthresh_z_obs > miss_threshold) THEN
          CALL apply_lower_thresh(ista, supob_lowthresh_z_obs, miss_threshold, zrpolar_obs)
        END IF
        ! zrmod_tmp is needed because zr_mod has intent(in) and cannot be changed by apply_lower_thresh():
        ALLOCATE (zrmod_tmp(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        zrmod_tmp = zr_mod
        IF (supob_lowthresh_z_sim > miss_threshold) THEN
          CALL apply_lower_thresh(ista, supob_lowthresh_z_sim, miss_threshold, zrmod_tmp )
        END IF
      END IF

      SELECT CASE (itype_supobing)
      CASE (0)
      CASE (1)
        CALL averagesuperobing_radar(ista, 'z', lat_mod, lon_mod, zrpolar_obs, zrmod_tmp, zrpolar_obs_err, zrobs_supob, &
                                                                      zrmod_supob, zrobserr_supob, lat_supob, lon_supob )
      CASE (2)
        CALL mediansuperobing_radar_vec(ista, 'z', lat_mod, lon_mod, zrpolar_obs, zrmod_tmp, zrpolar_obs_err, zrobs_supob, &
                                                                                         zrmod_supob, lat_supob, lon_supob )
      CASE default
          CALL abort_run (my_radar_id, 15079, &
               'ERROR: problem in '//TRIM(yzroutine)//': Wrong itype_supobing! Possible values are 0, 1,or 2', &
               'radar_process_output.f90, '//TRIM(yzroutine)//'()')
      END SELECT

      IF (ALLOCATED(zrmod_tmp))  DEALLOCATE(zrmod_tmp)


      IF (itype_supobing > 0) THEN

        CALL get_fileprefix_ascii_output ('zrsupobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zrobs_supob, TRIM(ADJUSTL(fileprefix)), &
             'Observed radar reflectivity', 'dBZ', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)

        CALL get_fileprefix_ascii_output ('zrsupsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zrmod_supob, TRIM(ADJUSTL(fileprefix)), &
             'Simulated radar reflectivity', 'dBZ', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)

        IF (.NOT.zlradwind_ok) THEN
          ! ... otherwise losupsim and lasupsim have already been written to file ...
          CALL get_fileprefix_ascii_output ('losupsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               lon_supob, TRIM(ADJUSTL(fileprefix)), &
               'Longitude of superobservation grid points', 'deg', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

          CALL get_fileprefix_ascii_output ('lasupsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               lat_supob, TRIM(ADJUSTL(fileprefix)), &
               'Latitude of superobservation grid points', 'deg', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF
        
      ENDIF



      !----------------------------------------------------------------------------
      ! .. Output of NetCDF feedback file (vr and z are written to the same file!):
      !----------------------------------------------------------------------------

      IF ( lfdbk_output .AND. it_is_time_for_feedback .AND. rs_meta(ista)%ldbzh_to_fdbk) THEN

#if defined(NUDGING) && defined(NETCDF)

        yerrstring(:) = ' '
        fdbk_err = 0

        ! Link correct input fields for obs_write_cdf_feedobs depending on superobing:
        SELECT CASE (itype_supobing)
        CASE (0)
          p_polar_obs => zrpolar_obs
          p_polar_obserr => zrpolar_obs_err
          p_polar_mod => zr_mod
          p_lat_mod   => lat_mod
          p_lon_mod   => lon_mod
        CASE (1,2)
          p_polar_obs => zrobs_supob
          p_polar_obserr => zrobserr_supob
          p_polar_mod => zrmod_supob
          p_lat_mod   => lat_supob
          p_lon_mod   => lon_supob
        END SELECT

        ! Check namelist setting of itype_metric_refl_fdbk, before calling obs_write_cdf_feedobs(),
        !  to ensure that it has a valid value within this routine:
        IF (itype_metric_refl_fdbk < 1 .OR. itype_metric_refl_fdbk > 3) THEN
          WRITE (yerrstring(:), '(a,i0,a)') 'Wrong itype_metric_refl_fdbk (', itype_metric_refl_fdbk, &
                                            ')! Possible values are 1, 2 or 3'
          CALL abort_run (my_radar_id, 35079, &
               'ERROR: problem in '//TRIM(yzroutine)//': '//TRIM(yerrstring), &
               'radar_process_output.f90, '//TRIM(yzroutine)//'()')
        END IF

        zrpolar_obs_err = baseval_obserr_dbz

        CALL get_fdbk_metadata ( idom, fdbk_meta_data )

        CALL obs_write_cdf_feedobs ( ista, itime, tot_num_observables_fdbk, 'z', time_mod/60.0, &
             p_polar_obs, zrpolar_obs_q, p_polar_obserr, &
             p_polar_mod, p_lat_mod, p_lon_mod, hr_mod, rpress,         &
             zazivec, zrvec, rs_meta(ista)%el_arr(1:rs_meta(ista)%nel), &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel,   &
             htop, &
             zlastpr, zlrefend, zlwrite_veridata,                &
             ydate_ini_mod,                   &
             fdbk_meta_data,                                     &
             nexce_rep, nexce_bdy, fdbk_err, yerrstring )

        IF (fdbk_err /= 0) THEN
          WRITE (*,'(a)') 'WARNING: problems in '//TRIM(yzroutine)//' writing z: '//TRIM(yerrstring)
        END IF

#else

        CALL abort_run (my_radar_id, 15078, &
             'ERROR: problem in '//TRIM(yzroutine)//                                                      &
             ': Writing of NetCDF feedback files not possible, missing flags -DNETCDF and/or -DNUDGING!'//&
             ' In ICON, you have to enable submodule dace_icon!', &
             'radar_process_output.f90, '//TRIM(yzroutine)//'()')

#endif

      END IF  ! lfdbk_output

    END IF  ! zlzradar_ok


    !---------------------------------------------------------------------
    ! .. If ZDR obs are available, do different kinds of output:
    !--------------------------------------------------------------------

    IF ( zlzdrradar_ok .OR. zlkdpradar_ok .OR. zlphidpradar_ok .OR. zlrhvradar_ok) THEN

      IF ( .NOT.lequal_azi_alldatasets ) THEN

        CALL init_vari(m_obs, -1)
        CALL init_vari(n_obs, -1)
        CALL init_vari(o_obs, -1)
!$omp parallel do private(iobs, tmppos_obs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          tmppos_obs = rs_data(ista)%ind_intp_obs_z(iobs)
          IF (tmppos_obs > missthr_int) THEN
            ! .. existing datum in NetCDF file (the data itself might still be missing values):
            CALL ind2sub3D(tmppos_obs, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                 m_obs(iobs), n_obs(iobs), o_obs(iobs))
          ENDIF
        END DO
!$omp end parallel do
      END IF

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
      IF (zlzdrradar_ok) THEN
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          IF ( rs_data(ista)%ind_intp_obs_z(iobs) > missthr_int .AND. rs_data(ista)%zdr_radar_obs(iobs) /= missing_obs) THEN
            ! .. Observed differential reflectivity is already in logarithmic dBZ units:
            zdrpolar_obs( m_obs(iobs),n_obs(iobs),o_obs(iobs) ) = rs_data(ista)%zdr_radar_obs(iobs)
          END IF
        END DO
!$omp end parallel do
      END IF

          
      IF (zlkdpradar_ok) THEN
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          IF ( rs_data(ista)%ind_intp_obs_z(iobs) > missthr_int .AND. rs_data(ista)%kdp_radar_obs(iobs) /= missing_obs) THEN
            ! .. Observed kdp is in deg/km:
            kdppolar_obs( m_obs(iobs),n_obs(iobs),o_obs(iobs) ) = rs_data(ista)%kdp_radar_obs(iobs)
          END IF
        END DO
!$omp end parallel do
      END IF

      IF (zlphidpradar_ok) THEN
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          IF ( rs_data(ista)%ind_intp_obs_z(iobs) > missthr_int .AND. rs_data(ista)%phidp_radar_obs(iobs) /= missing_obs) THEN
            ! .. Observed phidp is in deg:
            phidppolar_obs( m_obs(iobs),n_obs(iobs),o_obs(iobs) ) = rs_data(ista)%phidp_radar_obs(iobs)
          END IF
        END DO
!$omp end parallel do
      END IF

      IF (zlrhvradar_ok) THEN
!$omp parallel do private(iobs)
        DO iobs = 1, rs_data(ista)%nobs_obs(i_dbzh)
          IF ( rs_data(ista)%ind_intp_obs_z(iobs) > missthr_int .AND. rs_data(ista)%rhv_radar_obs(iobs) /= missing_obs) THEN
            rhvpolar_obs( m_obs(iobs),n_obs(iobs),o_obs(iobs) ) = rs_data(ista)%rhv_radar_obs(iobs)
          END IF
        END DO
!$omp end parallel do
      END IF

!!$ Special faking of obs data for testing purposes or for testsuite mode,
!!$  if the actual obs data reader does not yet exist or no obs are used:
      IF (rs_meta(ista)%icountry == i_fakeobs) THEN
        WRITE (*,*) '*** WARNING: Experimental run with faked differential reflectivity obs! ***'
        IF (zlzdrradar_ok) zdrpolar_obs  (:,:,:) = zdr_mod(:,:,:)
        IF (zlkdpradar_ok) kdppolar_obs  (:,:,:) = kdp_mod(:,:,:)
        IF (zlphidpradar_ok) phidppolar_obs  (:,:,:) = phidp_mod(:,:,:)
        IF (zlrhvradar_ok) rhvpolar_obs  (:,:,:) = rhv_mod(:,:,:)
      END IF
!!$ End of special faking!

      ! ..  Eliminate values above model top:
!$omp parallel do private(i,j,k) collapse(3)
      DO k=1, rs_meta(ista)%nel
        DO j=1, rs_meta(ista)%nra
          DO i=1, rs_meta(ista)%naz
            IF (hr_mod(i,j,k) > htop .OR. hr_mod(i,j,k) < miss_threshold) THEN
              IF (zlzdrradar_ok) zdrpolar_obs(i,j,k)   = miss_value
              IF (zlkdpradar_ok) kdppolar_obs(i,j,k)   = miss_value
              IF (zlphidpradar_ok) phidppolar_obs(i,j,k)   = miss_value
              IF (zlrhvradar_ok) rhvpolar_obs(i,j,k)   = miss_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do

      IF (zlzdrradar_ok) THEN
        CALL get_fileprefix_ascii_output ('zdrobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zdrpolar_obs, TRIM(ADJUSTL(fileprefix)), 'Observed differential radar reflectivity', 'dB', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)
        IF (ldebug_radsim) THEN
          CALL control_output(zdrobsoutputunit(ista), time_mod, rs_meta(ista), "ZDR [dB]", &
               zdrpolar_obs, (zdrpolar_obs > miss_threshold), miss_value, nobs_obs)
        END IF
      END IF
      IF (zlkdpradar_ok) THEN
        CALL get_fileprefix_ascii_output ('kdpobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             kdppolar_obs, TRIM(ADJUSTL(fileprefix)), 'Observed spec. differential phase shift', 'deg/km', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)
        IF (ldebug_radsim) THEN
          CALL control_output(kdpobsoutputunit(ista), time_mod, rs_meta(ista), "KDP [deg/km]", &
               kdppolar_obs, (kdppolar_obs > miss_threshold), miss_value, nobs_obs)
        END IF
      END IF
      IF (zlphidpradar_ok) THEN
        CALL get_fileprefix_ascii_output ('phidpobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             phidppolar_obs, TRIM(ADJUSTL(fileprefix)), 'Observed differential propagation phase', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)
      END IF
      IF (zlrhvradar_ok) THEN
        CALL get_fileprefix_ascii_output ('rhvobs', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             rhvpolar_obs, TRIM(ADJUSTL(fileprefix)), 'Observed H/V correllation coeff.', '-', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.TRUE.)
        IF (ldebug_radsim) THEN
          CALL control_output(rhvobsoutputunit(ista), time_mod, rs_meta(ista), "RHOHV [-]", &
               rhvpolar_obs, (rhvpolar_obs > miss_threshold), miss_value, nobs_obs)
        END IF
      END IF

    END IF


    !---------------------------------------------------------------------
    ! .. Cleanup:
    !---------------------------------------------------------------------

    IF ( zlradwind_ok .OR. zlzradar_ok ) THEN
      DEALLOCATE (vrpolar_obs, vrpolar_obs_q, zrpolar_obs, zrpolar_obs_q)
      DEALLOCATE (m_obs, n_obs, o_obs)
      IF (itype_supobing > 0) THEN
        DEALLOCATE (vrobs_supob, vrobserr_supob, vrmod_supob, zrobs_supob, zrobserr_supob, zrmod_supob, lon_supob, lat_supob)
      END IF
      IF (ALLOCATED(vrpolar_obs_err)) DEALLOCATE (vrpolar_obs_err, zrpolar_obs_err)
    END IF

    IF (ALLOCATED(zdrpolar_obs))   DEALLOCATE (zdrpolar_obs)
    IF (ALLOCATED(kdppolar_obs))   DEALLOCATE (kdppolar_obs)
    IF (ALLOCATED(phidppolar_obs)) DEALLOCATE (phidppolar_obs)
    IF (ALLOCATED(rhvpolar_obs))   DEALLOCATE (rhvpolar_obs)

#ifdef NETCDF
    IF (.NOT.(rs_meta(ista)%icountry == i_fakeobs)) THEN
      CALL read_obs_rad(ista, time_mod, 'cleanup')
    END IF
#endif
#if defined(NUDGING) && defined(NETCDF)
    IF (ALLOCATED(rpress)) DEALLOCATE (rpress)
#endif

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE output_radar_country_obs

  !==============================================================================
  !==============================================================================

  !==============================================================================
  !
  ! Service routines for output_radar_country_obs()
  !
  !==============================================================================


  SUBROUTINE qcrad_flag_eval(ista, cvarname, rdata_obs_q, rdata_obs)

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'qcrad_flag_eval'

    CHARACTER (LEN= * )      , INTENT (IN)  :: &
         cvarname               ! Tag for type of output field: 'vr' or 'z'

    INTEGER, INTENT(IN) :: ista         ! index of radar station

    REAL(KIND=dp), INTENT(IN), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obs_q

    REAL(KIND=dp), INTENT(INOUT), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obs

    INTEGER   :: i, j, k
    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE qcrad_flag_eval
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (TRIM(cvarname) == 'vr') THEN

!$omp parallel do private(i,j,k) collapse(3)
      DO k=1, rs_meta(ista)%nel
        DO j=1, rs_meta(ista)%nra
          DO i=1, rs_meta(ista)%naz
            IF (rdata_obs_q(i,j,k) /= 0.0_dp) THEN
              rdata_obs(i,j,k) = miss_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do

    ELSEIF (TRIM(cvarname) == 'z') THEN

!$omp parallel do private(i,j,k) collapse(3)
      DO k=1, rs_meta(ista)%nel
        DO j=1, rs_meta(ista)%nra
          DO i=1, rs_meta(ista)%naz
            IF (rdata_obs_q(i,j,k) /= 0.0_dp) THEN
              rdata_obs(i,j,k) = miss_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do

    ELSE

      CALL abort_run (my_radar_id, 80071, &
           'ERROR '//TRIM(yzroutine)//': quality control for cvarname='//TRIM(cvarname)//' not implemented!', &
           'radar_process_output.f90, '//TRIM(yzroutine))

    ENDIF

  END SUBROUTINE qcrad_flag_eval


  SUBROUTINE dealiasing_vr(ista, itime, rdata_mod, rdata_obs)

    !==============================================================================
    ! Dealiasing of observed radial wind rdata_obs with the help of
    !  simulated radial wind field rdata_mod as a first guess. It is assumed
    !  that the model data give valid radial wind values everywhere, i.e.,
    !  they should not have been filtered for minimum detectable signal, etc.
    !
    ! Method: Simple algorithm based on the assumption that the first guess
    !         from the model is "sufficiently" similar to the "true"
    !         observed value. Usually works well if the Nyquist velocity
    !         is large enough so that multiple foldings to not occur.
    !
    !==============================================================================

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'dealiasing_vr'

    INTEGER, INTENT(IN) :: ista         ! index of radar station
    INTEGER, INTENT(IN) :: itime        ! output time index

    REAL(KIND=dp), INTENT(IN), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_mod

    REAL(KIND=dp), INTENT(INOUT), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obs

    INTEGER    ::  iel, ira, iaz, K
    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE dealiasing_vr
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

!$omp parallel do private(iel,ira,iaz,K)
    DO iel = 1,  rs_meta(ista)%nel
      DO ira = 1,  rs_meta(ista)%nra
        DO iaz = 1,  rs_meta(ista)%naz

          IF (rdata_obs(iaz,ira,iel) > miss_threshold .AND. rdata_mod(iaz,ira,iel) > miss_threshold) THEN

            K = NINT((rdata_obs(iaz,ira,iel) - rdata_mod(iaz,ira,iel))/(2.0_dp*rs_meta(ista)%ext_nyq(iel,itime)))

            ! rs_meta(ista)%high_nyq: Nyquist-velocity corresponding to the higher of the two PRF
            ! rs_meta(ista)%ext_nyq : Effective extended Nyquist velocity after Dual PRF processing
            rdata_obs(iaz,ira,iel) =  rdata_obs(iaz,ira,iel) - 2.0_dp*K*rs_meta(ista)%ext_nyq(iel,itime)

          ELSE

            ! Enforce missing value for obs data at bins where no dealiasing was possible:
            rdata_obs(iaz,ira,iel) =  miss_value

          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE dealiasing_vr

  SUBROUTINE apply_lower_thresh(ista, lowthresh, missingthresh, rdata)

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'apply_lower_thresh'

    INTEGER, INTENT(IN) :: ista         ! index of radar station

    REAL(KIND=dp), INTENT(IN)       :: lowthresh, missingthresh

    REAL(KIND=dp), INTENT(INOUT), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata

    INTEGER :: i, j, k

!$omp parallel do private(i,j,k)
    DO k=1, rs_meta(ista)%nel
      DO j=1, rs_meta(ista)%nra
        DO i=1, rs_meta(ista)%naz
          IF (missingthresh < rdata(i,j,k) .AND. rdata(i,j,k) < lowthresh) THEN
            rdata(i,j,k) = lowthresh
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

  END SUBROUTINE apply_lower_thresh

  SUBROUTINE set_obserr_radwind(ista, zrdata, missingthresh, dBZ_0, e_o_0, dBZ_1, e_o_1, rdata_obserr)

    ! Sketch of observation error for radial wind (e_o) as function of dBZ:
    ! =====================================================================
    !
    ! Linear ramp function:   e_o(dBZ) = e_o_0 + (e_o_1 - e_o_0) * min( max( (dBZ - dBZ_0) / (dBZ_1 - dBZ_0) , 0.0) , 1.0)
    !
    !   e_o ^
    !       |               |             |
    ! e_o_0 |---------------\             |
    !       |               | \           |
    !       |               |   \         |
    !       |               |     \       |
    !       |               |       \     |
    !       |               |         \   |
    !       |               |           \ |
    ! e_o_1 | - - - - - - - | - - - - - - \------------------------------
    !       |               |             |
    !     -------------------------------------------------------------------->
    !       |             dBZ_0          dBZ_1                             dBZ

    INTEGER, INTENT(IN) :: ista         ! index of radar station
    REAL(kind=dp), INTENT(IN)       :: missingthresh, dBZ_0, e_o_0, dBZ_1, e_o_1
    REAL(kind=dp), INTENT(IN), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         zrdata
    REAL(kind=dp), INTENT(INOUT), &
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obserr

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'set_obserr_radwind'

    INTEGER :: i, j, k

    IF (ldebug_radsim) THEN
      WRITE (*,'(2a,i6.6)') TRIM(yzroutine), ' for station ', rs_meta(ista)%station_id
    END IF

!$omp parallel do private(i,j,k)
    DO k=1, rs_meta(ista)%nel
      DO j=1, rs_meta(ista)%nra
        DO i=1, rs_meta(ista)%naz
          IF (missingthresh < zrdata(i,j,k)) THEN
            rdata_obserr(i,j,k) = e_o_0 + (e_o_1 - e_o_0) * MIN( MAX( (zrdata(i,j,k)-dBZ_0)/(dBZ_1-dBZ_0), 0.0_dp), 1.0_dp)
          ELSE
            rdata_obserr(i,j,k) = MAX(e_o_0, e_o_1)
          END IF
        END DO
      END DO
    END DO
!$omp end parallel do

    IF (ldebug_radsim) THEN
      WRITE (*,'(3a,i6.6,4x,f0.1)') 'Done with ', TRIM(yzroutine), ' for station ', &
           rs_meta(ista)%station_id, MAXVAL(rdata_obserr)
    END IF

  END SUBROUTINE set_obserr_radwind

  SUBROUTINE averagesuperobing_radar(ista, cvarname, lat_mod, lon_mod, rdata_obs, rdata_mod, rdata_obserr, &
                                     rdata_obs_supob, rdata_mod_supob, rdata_obserr_supob, lat_supob, lon_supob)

    CHARACTER (LEN= * )      , INTENT (IN)  :: &
         cvarname               ! Tag for type of output field: 'vr' or 'z'

    !REAL(KIND=dp), INTENT(IN)       :: time_mod     ! seconds since model start
    INTEGER, INTENT(IN) :: ista         ! index of radar station
    REAL(KIND=dp), INTENT(IN), &                    ! simulated radar fields on original polar grid
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         lat_mod, lon_mod, rdata_obs, rdata_mod, rdata_obserr

    REAL(KIND=dp), INTENT(INOUT), &                 ! simulated geogr. coordinates, to be overwritten
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         lat_supob, lon_supob                           !  with the "quasi-cartesian" geographic coordinates for locations
                                                        !  which are stored in the superob'ed output fields

    REAL(KIND=dp), INTENT(OUT), &                   ! superobing radar fields for feedback files
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obs_supob, rdata_mod_supob, rdata_obserr_supob


    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'averagesuperobing_radar'

    INTEGER    :: i, j, iel, ira, iaz, ira_sub, iaz_sub, i_cart, icp, ip, np, &
         naz_average, iaz_average, ira_average, iaz_range, nobsmax, ncp, ncploc, aw_loc

    REAL(KIND=dp)                       :: w_tmp, dis_tmp, clon_r, clat_r, z_tmp, rlat_r, rlon_r, &
         dhori, w_obserr   ! dhori: horizontal radius of cressman weighting domain;
    ! dvert: vertical radius of cressman weighting domain


    INTEGER, ALLOCATABLE, SAVE :: count_obs(:), count_mod(:), indloc(:)
!$omp threadprivate(count_obs, count_mod, indloc)

    REAL    (KIND=dp), ALLOCATABLE, SAVE :: &
         w_obs(:)            ,&
         w_mod(:)            ,&
         rdata_obs_wsum(:)   ,&
         rdata_obserr_wsum(:),&
         rdata_mod_wsum(:)   ,&
         rdata_obs_sum(:)    ,&
         rdata_mod_sum(:)    ,&
         rdata_obs_varsum(:) ,&
         rdata_mod_varsum(:)
!$omp threadprivate(w_obs,w_mod,rdata_obs_wsum,rdata_obserr_wsum,rdata_mod_wsum,rdata_obs_sum,rdata_mod_sum,rdata_obs_varsum,rdata_mod_varsum)

    REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE superobing_radar
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    SELECT CASE (TRIM(cvarname))
    CASE ('vr','z')
      CONTINUE
    CASE default
      CALL abort_run (my_radar_id, 80073, &
          'ERROR '//TRIM(yzroutine)//': superobing for cvarname='//TRIM(cvarname)//' not implemented!', &
          'radar_process_output.f90, '//TRIM(yzroutine))
    END SELECT

    ! allocate and initialize local aux arrays:
    !  number of radar points:
    nobsmax = rs_meta(ista)%nra * rs_meta(ista)%naz * rs_meta(ista)%nel
    !  total number of cartesian points:
    ncp     =  cart_data(ista)%ni_tot *  cart_data(ista)%nj_tot
    !  maximum number of cartesian points within radar range enclosing rectangle:
    ncploc  =  cart_data(ista)%ni     *  cart_data(ista)%nj

    ALLOCATE(cart_data(ista)%dis(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%calt(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%rdata_ind(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%cartdata_ind(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%aw(rs_meta(ista)%nel,ncploc))

    ! initialize output variables with missing values:
    CALL init_vari(rdata_obs_supob, miss_value)
    CALL init_vari(rdata_mod_supob, miss_value)
    IF (TRIM(cvarname) == 'vr') THEN
      CALL init_vari(rdata_obserr_supob, baseval_obserr_vr)
    ELSE IF (TRIM(cvarname) == 'z') THEN
      CALL init_vari(rdata_obserr_supob, baseval_obserr_dbz)
    END IF
    CALL init_vari(cart_data(ista)%dis          , HUGE(0.0_dp) )
    CALL init_vari(cart_data(ista)%calt         , miss_value   )
    CALL init_vari(cart_data(ista)%rdata_ind    , missval_int  )
    CALL init_vari(cart_data(ista)%cartdata_ind , missval_int  )
    CALL init_vari(cart_data(ista)%aw           , 0            )

    ! Specify the cartesian coordinats and search for the
    ! closest radar point of each cartesian grid point. Search its corresponding distance and azimuthal averaging interval

!$omp parallel
    ! allocate work arrays:
    ALLOCATE(indloc(ncp))

    ALLOCATE(count_mod(ncploc))
    ALLOCATE(count_obs(ncploc))
    ALLOCATE(rdata_mod_sum(ncploc))
    ALLOCATE(rdata_obs_sum(ncploc))
    ALLOCATE(rdata_mod_varsum(ncploc))
    ALLOCATE(rdata_obs_varsum(ncploc))
    ALLOCATE(w_obs(ncploc))
    ALLOCATE(w_mod(ncploc))
    ALLOCATE(rdata_obs_wsum(ncploc))
    ALLOCATE(rdata_obserr_wsum(ncploc))
    ALLOCATE(rdata_mod_wsum(ncploc))

!$omp do private(iel,ira,iaz,np,rlon_r,rlat_r,i,j,ip,icp,clon_r,clat_r,dis_tmp), &
!$omp& private(naz_average,iaz_average,ira_average,aw_loc,w_tmp,w_obserr,ira_sub,iaz_sub,z_tmp)
    DO iel = 1,  rs_meta(ista)%nel

      np = 0
      indloc = -999999

!CDIR NODEP
      DO ira = 1,  rs_meta(ista)%nra
!CDIR NODEP
        DO iaz = 1,  rs_meta(ista)%naz

          rlon_r = rla2rlarot(lat_mod(iaz,ira,iel), lon_mod(iaz,ira,iel), &
               cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)
          rlat_r = phi2phirot(lat_mod(iaz,ira,iel), lon_mod(iaz,ira,iel), &
               cart_data(ista)%pollat, cart_data(ista)%pollon)

          ! Determine the indices of the closest cartesian point:
          i = NINT((rlon_r - cart_data(ista)%startlon-cart_data(ista)%width_deg) / cart_data(ista)%dlon) + 1
          j = NINT((rlat_r - cart_data(ista)%startlat-cart_data(ista)%width_deg) / cart_data(ista)%dlat) + 1

          IF (i >= 1 .AND. i <= cart_data(ista)%ni_tot .AND. j >= 1 .AND. j <= cart_data(ista)%nj_tot) THEN

            ! Continuous combined cartesian index to be stored in the local index array indloc:
            icp =  i + (j - 1)*cart_data(ista)%ni_tot
            IF (indloc(icp) == -999999) THEN
              ! This is a new point. Generate a new local index and store it in indloc:
              np = np + 1
              indloc(icp) = np
              ip = np
            ELSE
              ! This is a known point. Retrieve its local index from indloc:
              ip = indloc(icp)
            END IF

            ! Rotated coordinates of the nearest cartesian point:
            clon_r = cart_data(ista)%startlon + cart_data(ista)%width_deg + (i-1)*cart_data(ista)%dlon
            clat_r = cart_data(ista)%startlat + cart_data(ista)%width_deg + (j-1)*cart_data(ista)%dlat

            ! "Fake" Euklidean distance to the radar point:
            dis_tmp = SQRT((rlat_r-clat_r)**2 + (rlon_r-clon_r)**2)

            ! If this radar point is nearer to the cartesian point than any previous radar point, store it:
            IF (dis_tmp < cart_data(ista)%dis(iel,ip)) THEN
              cart_data(ista)%dis(iel,ip)          = dis_tmp
              cart_data(ista)%calt(iel,ip)         = miss_value  ! not further used so far
              cart_data(ista)%rdata_ind(iel,ip)    = iaz + (ira-1)*rs_meta(ista)%naz
              cart_data(ista)%cartdata_ind(iel,ip) = icp
              ! half azimuthal averaging interval
              cart_data(ista)%aw(iel,ip) = CEILING ( &
                   ATAN(cart_data(ista)%width/(ira*rs_meta(ista)%ra_inc))*raddeg / &
                   rs_meta(ista)%az_inc )
            END IF

          END IF

        END DO  ! end of nra
      END DO  ! end of naz

      count_mod(:)     = 0
      count_obs(:)     = 0
      w_obs(:)            = 0.0_dp
      w_mod(:)            = 0.0_dp
      rdata_obs_wsum(:)   = 0.0_dp
      rdata_obserr_wsum(:)   = 0.0_dp
      rdata_mod_wsum(:)   = 0.0_dp
      rdata_mod_sum(:) = 0.0_dp
      rdata_obs_sum(:) = 0.0_dp
      rdata_mod_varsum(:) = 0.0_dp
      rdata_obs_varsum(:) = 0.0_dp

      IF (TRIM(cvarname) == 'vr') THEN

        naz_average = MIN( MAXVAL(cart_data(ista)%aw(iel,:)), cart_data(ista)%aw_max)

        DO iaz_average = -naz_average, naz_average, 1
          DO ira_average = -cart_data(ista)%nra_average, cart_data(ista)%nra_average, 1
            DO ip = 1, np

                ! half azimuthal averaging interval, not bigger than aw_max
              aw_loc = MIN(cart_data(ista)%aw_max, cart_data(ista)%aw(iel,ip) )

              iaz = MODULO((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1

              ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

              ira_sub = ira + ira_average

              IF  (ABS(iaz_average) <= aw_loc          .AND. &
                   ira_sub >= 1                        .AND. &
                   ira_sub <=rs_meta(ista)%nra) THEN

                iaz_sub = MODULO(iaz + iaz_average - 1, rs_meta(ista)%naz) + 1


                ! increased radius influence by 5 % to avoid 0 weight for the outermost points:
                w_tmp =  (REAL((1.05*aw_loc)**2 - iaz_average**2)   &
                        / REAL((1.05*aw_loc)**2 + iaz_average**2))  &
                        *(REAL((1.05*cart_data(ista)%nra_average)**2 - ira_average**2) &
                        / REAL((1.05*cart_data(ista)%nra_average)**2 + ira_average**2))

                w_obserr = 1.0 / rdata_obserr(iaz_sub,ira_sub,iel)

                ! superobing the model data with Cressman weighting function
                IF (rdata_mod(iaz_sub,ira_sub,iel) > miss_threshold                 &
                     .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_vr) THEN

                  count_mod(ip) = count_mod(ip) + 1
                  w_mod(ip) = w_mod(ip) + (w_tmp * w_obserr)

                  rdata_mod_sum(ip)    = rdata_mod_sum(ip)    + rdata_mod(iaz_sub,ira_sub,iel)
                  rdata_mod_varsum(ip) = rdata_mod_varsum(ip) + rdata_mod(iaz_sub,ira_sub,iel)*rdata_mod(iaz_sub,ira_sub,iel)
                  rdata_mod_wsum(ip)   = rdata_mod_wsum(ip)   + w_tmp*w_obserr*rdata_mod(iaz_sub,ira_sub,iel)

                END IF

                ! superobing the obs data with Cressman weighting function
                IF (rdata_obs(iaz_sub,ira_sub,iel) > miss_threshold       &
                     .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_vr) THEN

                  count_obs(ip) = count_obs(ip) + 1
                  w_obs(ip) = w_obs(ip) + (w_tmp * w_obserr)

                  rdata_obs_sum(ip)     = rdata_obs_sum(ip)     + rdata_obs(iaz_sub,ira_sub,iel)
                  rdata_obs_varsum(ip)  = rdata_obs_varsum(ip)  + rdata_obs(iaz_sub,ira_sub,iel)*rdata_obs(iaz_sub,ira_sub,iel)
                  rdata_obs_wsum(ip)    = rdata_obs_wsum(ip)    + w_tmp*w_obserr*rdata_obs(iaz_sub,ira_sub,iel)
                  rdata_obserr_wsum(ip) = rdata_obserr_wsum(ip) + w_tmp*w_obserr*rdata_obserr(iaz_sub,ira_sub,iel)

                END IF

              END IF  ! radar point within averaging region

            END DO  ! np
          END DO ! end of nra_average
        END DO ! end of naz_average

        ! Compute variance of data within the superobing regions:
        DO ip = 1, np
          IF (count_mod(ip) >= supob_nrb) THEN
            rdata_mod_varsum(ip) =  (rdata_mod_varsum(ip) - rdata_mod_sum(ip)*rdata_mod_sum(ip)/count_mod(ip)) / count_mod(ip)
          END IF
          IF (count_obs(ip) >= supob_nrb) THEN
            rdata_obs_varsum(ip) =  (rdata_obs_varsum(ip) - rdata_obs_sum(ip)*rdata_obs_sum(ip)/count_obs(ip)) / count_obs(ip)
          END IF
        END DO

        DO ip = 1, np

          iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
          ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

          IF (count_mod(ip) >= supob_nrb .AND. w_mod(ip) > 0.0_dp .AND. &
               rdata_mod_varsum(ip) <= supob_vrw*supob_vrw) THEN
            rdata_mod_supob(iaz,ira,iel) = rdata_mod_wsum(ip)/w_mod(ip)
            IF (ABS(rdata_mod_supob(iaz,ira,iel)) < eps_vr) THEN
              rdata_mod_supob(iaz,ira,iel) = 0.0_dp
            END IF
          END IF

          IF (count_obs(ip) >= supob_nrb .AND. w_obs(ip) > 0.0_dp .AND. &
               rdata_obs_varsum(ip) <= supob_vrw*supob_vrw) THEN
            rdata_obs_supob(iaz,ira,iel)    = rdata_obs_wsum(ip)    / w_obs(ip)
            rdata_obserr_supob(iaz,ira,iel) = rdata_obserr_wsum(ip) / w_obs(ip)
          ELSE
            rdata_obs_supob(iaz,ira,iel)    = reject_value
            ! rdata_obserr_supob stays at its initial value (baseval_obserr_vr)
          END IF

        END DO

      ELSE IF (TRIM(cvarname) == 'z') THEN

        naz_average = MAXVAL(cart_data(ista)%aw(iel,:))

        DO iaz_average = -naz_average, naz_average, 1
          DO ira_average = -cart_data(ista)%nra_average, cart_data(ista)%nra_average, 1
            DO ip = 1, np

                ! half azimuthal averaging interval
!!$              aw_loc = MIN(cart_data(ista)%aw_max, cart_data(ista)%aw(iel,ip) )
              aw_loc = cart_data(ista)%aw(iel,ip)

              iaz = MODULO((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
              ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1
              ira_sub = ira + ira_average

              IF  (ABS(iaz_average) <= aw_loc          .AND. &
                   ira_sub >= 1                        .AND. &
                   ira_sub <=rs_meta(ista)%nra) THEN

                iaz_sub = MODULO(iaz + iaz_average - 1, rs_meta(ista)%naz) + 1

!!$ UB: increased radius influence by 5 % to avoid 0 weight for the outermost points:
                w_tmp = ( REAL((1.05*aw_loc)**2 - iaz_average**2)   &
                        / REAL((1.05*aw_loc)**2 + iaz_average**2))  &
                        *(REAL((1.05*cart_data(ista)%nra_average)**2 - ira_average**2) &
                        / REAL((1.05*cart_data(ista)%nra_average)**2 + ira_average**2))

                ! superobing the model data with Cressman weighting function, taking into account valid 0's:
                IF (rdata_mod(iaz_sub,ira_sub,iel) > miss_threshold &
                     .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_z) THEN

                  count_mod(ip) = count_mod(ip) + 1
                  w_mod(ip) = w_mod(ip) + w_tmp

                  IF (rdata_mod(iaz_sub,ira_sub,iel) >= dBZ_crit_radar) THEN
!!$                    z_tmp = 10.0_dp**(0.1_dp* rdata_mod(iaz_sub,ira_sub,iel))
                    z_tmp = EXP(rdata_mod(iaz_sub,ira_sub,iel)*ln10_o10)
                    rdata_mod_sum(ip)    = rdata_mod_sum(ip)    + z_tmp
                    rdata_mod_varsum(ip) = rdata_mod_varsum(ip) + z_tmp*z_tmp
                    rdata_mod_wsum(ip)   = rdata_mod_wsum(ip)   + w_tmp*z_tmp
                  END IF
                END IF

                ! superobing the obs data with Cressman weighting function, taking into account valid 0's:
                IF (rdata_obs(iaz_sub,ira_sub,iel) > miss_threshold       &
                     .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_z) THEN

                  count_obs(ip) = count_obs(ip) + 1
                  w_obs(ip) = w_obs(ip) + w_tmp

                  IF (rdata_obs(iaz_sub,ira_sub,iel) >= dBZ_crit_radar) THEN
!!$                    z_tmp = 10.0_dp**(0.1_dp* rdata_obs(iaz_sub,ira_sub,iel))
                    z_tmp = EXP(rdata_obs(iaz_sub,ira_sub,iel)*ln10_o10)
                    rdata_obs_sum(ip)     = rdata_obs_sum(ip)     + z_tmp
                    rdata_obs_varsum(ip)  = rdata_obs_varsum(ip)  + z_tmp*z_tmp
                    rdata_obs_wsum(ip)    = rdata_obs_wsum(ip)    + w_tmp*z_tmp
                    rdata_obserr_wsum(ip) = rdata_obserr_wsum(ip) + w_tmp
                  END IF

                END IF

              END IF ! radar point within averaging region


            END DO  ! np
          END DO ! end of nra_average
        END DO ! end of naz_average

        ! Compute variance of data within the superobing regions:
        DO ip = 1, np
          IF (count_mod(ip) >= supob_nrb) THEN
            rdata_mod_varsum(ip) =  (rdata_mod_varsum(ip) - rdata_mod_sum(ip)*rdata_mod_sum(ip)/count_mod(ip)) / count_mod(ip)
          END IF
          IF (count_obs(ip) >= supob_nrb) THEN
            rdata_obs_varsum(ip) =  (rdata_obs_varsum(ip) - rdata_obs_sum(ip)*rdata_obs_sum(ip)/count_obs(ip)) / count_obs(ip)
          END IF
        END DO

        DO ip = 1, np

          iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
          ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

          ! .. look if we have enough valid data:
          IF (count_mod(ip) >= supob_nrb .AND. rdata_mod_varsum(ip) <= supob_rfl*supob_rfl) THEN
            ! .. is the mean value > 0.0?
            IF (rdata_mod_wsum(ip) >= Z_crit_radar*w_mod(ip) .AND. w_mod(ip) > 0.0_dp) THEN
              rdata_mod_supob(iaz,ira,iel) = 10.0_dp*LOG10(rdata_mod_wsum(ip)/w_mod(ip))
            ELSE
              ! .. this is a valid mean value of 0 (zero_value in log space):
              rdata_mod_supob(iaz,ira,iel) = zero_value
            END IF
          END IF

          ! .. look if we have enough valid data:
          IF (count_obs(ip) >= supob_nrb .AND. rdata_obs_varsum(ip) <= supob_rfl*supob_rfl) THEN
            ! .. is the mean value > 0.0?
            IF (rdata_obs_wsum(ip) >= Z_crit_radar*w_obs(ip) .AND. w_obs(ip) > 0.0_dp) THEN
              rdata_obs_supob(iaz,ira,iel)    = 10.0_dp*LOG10(rdata_obs_wsum(ip) / w_obs(ip))
              rdata_obserr_supob(iaz,ira,iel) = (rdata_obserr_wsum(ip)           / w_obs(ip))
            ELSE
              ! .. this is a valid mean value of 0 (zero_value in log space):
              rdata_obs_supob(iaz,ira,iel)    = zero_value
              rdata_obserr_supob(iaz,ira,iel) = zero_value
            END IF
          ELSE
            rdata_obs_supob(iaz,ira,iel)    = reject_value
            ! rdata_obserr_supob stays at its initial value (baseval_obserr_dbz)
          END IF


        END DO

      ENDIF   !  cvarname = 'vr' or 'z'

      ! .. write the geographic (quasi-cartesian) coordinates to the INOUT fields lat_supob and lon_supob
      !     at the suberob'ed locations:
      DO ip = 1, np

        i = MOD((cart_data(ista)%cartdata_ind(iel,ip)-1),cart_data(ista)%ni_tot) + 1
        j = (cart_data(ista)%cartdata_ind(iel,ip)-1)/cart_data(ista)%ni_tot  + 1

        iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
        ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

        clon_r = cart_data(ista)%startlon + cart_data(ista)%width_deg + (i-1)*cart_data(ista)%dlon
        clat_r = cart_data(ista)%startlat + cart_data(ista)%width_deg + (j-1)*cart_data(ista)%dlat

        lat_supob(iaz,ira,iel) = phirot2phi(clat_r, clon_r, &
             cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)
        lon_supob(iaz,ira,iel) = rlarot2rla(clat_r, clon_r, &
             cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)

      END DO

    END DO ! end of nel
!$omp end do

    DEALLOCATE(indloc)
    DEALLOCATE(count_obs)
    DEALLOCATE(count_mod)
    DEALLOCATE(rdata_mod_sum)
    DEALLOCATE(rdata_obs_sum)
    DEALLOCATE(rdata_mod_varsum)
    DEALLOCATE(rdata_obs_varsum)
    DEALLOCATE(w_obs)
    DEALLOCATE(w_mod)
    DEALLOCATE(rdata_obs_wsum)
    DEALLOCATE(rdata_obserr_wsum)
    DEALLOCATE(rdata_mod_wsum)
!$omp end parallel

    IF (ldebug_radsim) WRITE (*,*) 'Done with ', TRIM(yzroutine), ' on proc ', my_radar_id

  END SUBROUTINE averagesuperobing_radar


  SUBROUTINE mediansuperobing_radar_vec(ista, cvarname, lat_mod, lon_mod, rdata_obs, rdata_mod, rdata_obserr, &
                                        rdata_obs_supob, rdata_mod_supob, lat_supob, lon_supob)

!!$ rdata_obserr just a dummy input at the moment. No superobing of obserr implemented yet!

    CHARACTER (LEN= * )      , INTENT (IN)  :: &
         cvarname               ! Tag for type of output field: 'vr' or 'z'

    INTEGER, INTENT(IN) :: ista         ! index of radar station
    REAL(KIND=dp), INTENT(IN), &                    ! simulated radar fields for feedback files
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         lat_mod, lon_mod, rdata_obs, rdata_mod, rdata_obserr

    REAL(KIND=dp), INTENT(INOUT), &                 ! simulated geogr. coordinates, to be overwritten
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         lat_supob, lon_supob                           !  with the "quasi-cartesian" geographic coordinates for locations
                                                        !  which are stored in the superob'ed output fields
    REAL(KIND=dp), INTENT(OUT), &                    ! superobing radar fields for feedback files
         DIMENSION( rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel ) :: &
         rdata_obs_supob, rdata_mod_supob

    CHARACTER (LEN=*), PARAMETER :: yzroutine = 'mediansuperobing_radar_vec'

    INTEGER    :: i, j, iel, ira, iaz, ira_sub, iaz_sub, icp, ip, np, ict, jct, nct, &
         naz_average, iaz_average, ira_average, nobsmax, ncp, ncploc, aw_loc

    REAL(KIND=dp)                       :: dis_tmp, w_tmp, clon_r, clat_r, z_tmp, rlat_r, rlon_r, &
         mod_tmp, obs_tmp, dhori ! dhori: horizontal radius of cressman weighting domain;
    ! dvert: vertical radius of cressman weighting domain

    INTEGER  , ALLOCATABLE :: count_mod(:), count_obs(:), indloc(:)

    REAL    (KIND=dp), ALLOCATABLE :: &
         rdata_obs_tmp(:,:,:)  ,&
         rdata_mod_tmp(:,:,:)  ,&
         rdata_mod_sum(:)      ,&
         rdata_obs_sum(:)      ,&
         rdata_mod_ave(:)      ,&
         rdata_obs_ave(:)      ,&
         rdata_mod_varsum(:)   ,&
         rdata_obs_varsum(:)

    REAL(kind=dp), PARAMETER :: ln10_o10 = LOG(10.0_dp)*0.1_dp


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE mediansuperobing_radar_vec
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

!!! THIS ROUTINE IS NOT FULLY FUNCTIONAL AT THE MOMENT!!!
        CALL abort_run (my_radar_id, 15080, &
               'ERROR: '//TRIM(yzroutine)//' IS CURRENTLY NOT FUNCTIONAL! itype_supobing = 2 is not possible!', &
               'radar_process_output.f90, '//TRIM(yzroutine)//'()')


    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    SELECT CASE (TRIM(cvarname))
    CASE ('vr','z')
      CONTINUE
    CASE default
      CALL abort_run (my_radar_id, 80074, &
          'ERROR '//TRIM(yzroutine)//': superobing for cvarname='//TRIM(cvarname)//' not implemented!', &
          'radar_process_output.f90, '//TRIM(yzroutine))
    END SELECT

    ! allocate and initialize local aux arrays:
    !  number of radar points:
    nobsmax = rs_meta(ista)%nra * rs_meta(ista)%naz * rs_meta(ista)%nel
    !  total number of cartesian points:
    ncp     =  cart_data(ista)%ni_tot *  cart_data(ista)%nj_tot
    !  maximum number of cartesian points within radar range enclosing rectangle:
    ncploc  =  cart_data(ista)%ni     *  cart_data(ista)%nj

    ! allocate work arrays:
    ALLOCATE(indloc(ncp))

    ALLOCATE(count_mod(ncploc))
    ALLOCATE(count_obs(ncploc))
    ALLOCATE(rdata_mod_sum(ncploc))
    ALLOCATE(rdata_obs_sum(ncploc))
    ALLOCATE(rdata_mod_ave(ncploc))
    ALLOCATE(rdata_obs_ave(ncploc))
    ALLOCATE(rdata_mod_varsum(ncploc))
    ALLOCATE(rdata_obs_varsum(ncploc))

    ALLOCATE(cart_data(ista)%dis(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%calt(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%rdata_ind(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%cartdata_ind(rs_meta(ista)%nel,ncploc))
    ALLOCATE(cart_data(ista)%aw(rs_meta(ista)%nel,ncploc))

!CDIR BEGIN COLLAPSE
    rdata_obs_supob = miss_value
    rdata_mod_supob = miss_value

    cart_data(ista)%dis        = HUGE(0.0_dp)
    cart_data(ista)%calt       = miss_value
    cart_data(ista)%rdata_ind  = missval_int
    cart_data(ista)%cartdata_ind  = missval_int
    cart_data(ista)%aw         = 0
!CDIR END

    ! Specify the cartesian coordinats and search for the
    ! closest radar point of each cartesian grid point. Search its corresponding distance and azimuthal averaging interval

    DO iel = 1,  rs_meta(ista)%nel

      np = 0
      indloc = -999999

!CDIR NODEP
      DO ira = 1,  rs_meta(ista)%nra
!CDIR NODEP
        DO iaz = 1,  rs_meta(ista)%naz

          rlon_r = rla2rlarot(lat_mod(iaz,ira,iel), lon_mod(iaz,ira,iel), &
               cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)
          rlat_r = phi2phirot(lat_mod(iaz,ira,iel), lon_mod(iaz,ira,iel), &
               cart_data(ista)%pollat, cart_data(ista)%pollon)

          ! Determine the indices of the closest cartesian point:
          i = NINT((rlon_r - cart_data(ista)%startlon-cart_data(ista)%width_deg) / cart_data(ista)%dlon) + 1
          j = NINT((rlat_r - cart_data(ista)%startlat-cart_data(ista)%width_deg) / cart_data(ista)%dlat) + 1

          IF (i < 1 .OR. i > cart_data(ista)%ni_tot .OR. j < 1 .OR. j > cart_data(ista)%nj_tot) THEN
           CYCLE
          END IF

          ! Continuous combined cartesian index to be stored in the local index array indloc:
          icp =  i + (j - 1)*cart_data(ista)%ni_tot
          IF (indloc(icp) == -999999) THEN
            ! This is a new point. Generate a new local index and store it in indloc:
            np = np + 1
            indloc(icp) = np
            ip = np
          ELSE
            ! This is a known point. Retrieve its local index from indloc:
            ip = indloc(icp)
          END IF

          ! Rotated coordinates of the nearest cartesian point:
          clon_r = cart_data(ista)%startlon + cart_data(ista)%width_deg + (i-1)*cart_data(ista)%dlon
          clat_r = cart_data(ista)%startlat + cart_data(ista)%width_deg + (j-1)*cart_data(ista)%dlat

          ! "Fake" Euklidean distance to the radar point:
          dis_tmp = sqrt((rlat_r-clat_r)**2 + (rlon_r-clon_r)**2)

          ! If this radar point is nearer to the cartesian point than any previous radar point, store it:
          IF (dis_tmp < cart_data(ista)%dis(iel,ip)) THEN

            cart_data(ista)%dis(iel,ip)          = dis_tmp
            cart_data(ista)%calt(iel,ip)         = miss_value  ! not further used so far
            cart_data(ista)%rdata_ind(iel,ip)    = iaz + (ira-1)*rs_meta(ista)%naz
            cart_data(ista)%cartdata_ind(iel,ip) = icp

            ! half azimuthal averaging interval
            cart_data(ista)%aw(iel,ip) = CEILING ( &
                 ATAN(cart_data(ista)%width/(ira*rs_meta(ista)%ra_inc))*raddeg / &
                 rs_meta(ista)%az_inc )

          END IF

        END DO  ! end of nra
      END DO  ! end of naz


      count_mod = 0
      count_obs = 0

      rdata_mod_sum = 0.0_dp
      rdata_obs_sum = 0.0_dp

      rdata_mod_ave = 0.0_dp
      rdata_obs_ave = 0.0_dp

      rdata_mod_varsum = 0.0_dp
      rdata_obs_varsum = 0.0_dp

      IF (TRIM(cvarname) == 'vr') THEN

        naz_average = MIN( MAXVAL(cart_data(ista)%aw(iel,:)), cart_data(ista)%aw_max)
        !  maximum number of raw radar bins for one superobing point:
        nct = (2*cart_data(ista)%aw_max+1)*(2*cart_data(ista)%nra_average+1)

        ALLOCATE(rdata_obs_tmp(rs_meta(ista)%nel,ncploc,nct))
        ALLOCATE(rdata_mod_tmp(rs_meta(ista)%nel,ncploc,nct))
!CDIR BEGIN COLLAPSE
        rdata_obs_tmp   = miss_value
        rdata_mod_tmp   = miss_value
!CDIR END

        DO iaz_average = -naz_average, naz_average, 1
          DO ira_average = -cart_data(ista)%nra_average, cart_data(ista)%nra_average, 1
            DO ip = 1, np

              IF (cart_data(ista)%rdata_ind(iel,ip) > 0) THEN

                ! half azimuthal averaging interval, not bigger than aw_max
                aw_loc = MIN(cart_data(ista)%aw_max, cart_data(ista)%aw(iel,ip) )

                iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
                ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1
                ira_sub = ira + ira_average

                IF  (ABS(iaz_average) <= aw_loc .AND. &
                     ira_sub >= 1                                        .AND. &
                     ira_sub <=rs_meta(ista)%nra)                             THEN

                  iaz_sub = MODULO(iaz + iaz_average - 1, rs_meta(ista)%naz) + 1

                  IF (rdata_mod(iaz_sub,ira_sub,iel) > miss_threshold &
                       .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_vr) THEN
                    count_mod(ip) =  count_mod(ip) + 1
                    rdata_mod_sum(ip) = rdata_mod_sum(ip) +  rdata_mod(iaz_sub,ira_sub,iel)
                    rdata_mod_tmp(iel,ip,count_mod(ip)) = rdata_mod(iaz_sub,ira_sub,iel)
                  END IF

                  IF (rdata_obs(iaz_sub,ira_sub,iel) > miss_threshold &
                       .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_vr) THEN
                    count_obs(ip) =  count_obs(ip) + 1
                    rdata_obs_sum(ip) = rdata_obs_sum(ip) +  rdata_obs(iaz_sub,ira_sub,iel)
                    rdata_obs_tmp(iel,ip,count_obs(ip)) = rdata_obs(iaz_sub,ira_sub,iel)
                  END IF

                END IF ! for ip with valid vaules

              END IF

            END DO  ! ncp
          END DO ! end of nra_average
        END DO ! end of naz_average

        DO ip = 1, np
          IF (count_mod(ip) >= supob_nrb) THEN
            rdata_mod_ave(ip) = rdata_mod_sum(ip)/count_mod(ip)
          END IF
          IF (count_obs(ip) >= supob_nrb) THEN
            rdata_obs_ave(ip) = rdata_obs_sum(ip)/count_obs(ip)
          END IF
        END DO


        DO ict = 1, nct
          DO jct = nct, ict+1, -1
            DO ip = 1, np

              IF (count_mod(ip) >= supob_nrb) THEN
                IF (rdata_mod_tmp(iel,ip,jct) < rdata_mod_tmp(iel,ip,jct-1) &
                     .AND. rdata_mod_tmp(iel,ip,jct) > miss_threshold) THEN

                  mod_tmp =rdata_mod_tmp(iel,ip,jct-1)
                  rdata_mod_tmp(iel,ip,jct-1) = rdata_mod_tmp(iel,ip,jct)
                  rdata_mod_tmp(iel,ip,jct) = mod_tmp

                END IF
              END IF

              IF (count_obs(ip) >= supob_nrb) THEN
                IF (rdata_obs_tmp(iel,ip,jct) < rdata_obs_tmp(iel,ip,jct-1)  &
                     .AND. rdata_obs_tmp(iel,ip,jct) > miss_threshold) THEN

                  obs_tmp =rdata_obs_tmp(iel,ip,jct-1)
                  rdata_obs_tmp(iel,ip,jct-1) = rdata_obs_tmp(iel,ip,jct)
                  rdata_obs_tmp(iel,ip,jct) = obs_tmp

                END IF
              END IF

            END DO
          END DO
        END DO

!!$ UB: this can be made more efficient, like in averagesuperobing_vec()!
        DO ict = 1, nct
          DO ip = 1, np

            IF (count_mod(ip) >= supob_nrb .AND. ict <= count_mod(ip)) THEN
              rdata_mod_varsum(ip) =  rdata_mod_varsum(ip) + (rdata_mod_tmp(iel,ip,ict) - rdata_mod_ave(ip))**2
            END IF

            IF (count_obs(ip) >= supob_nrb .AND. ict <= count_obs(ip)) THEN
              rdata_obs_varsum(ip) =  rdata_obs_varsum(ip) + (rdata_obs_tmp(iel,ip,ict) - rdata_obs_ave(ip))**2
            END IF

          END DO
        END DO

        DO ip = 1, np

          IF (cart_data(ista)%rdata_ind(iel,ip) > 0) THEN

            iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
            ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

            IF (count_mod(ip) >= supob_nrb .AND. rdata_mod_varsum(ip)/count_mod(ip) <= supob_vrw*supob_vrw) THEN
              rdata_mod_supob(iaz,ira,iel) = rdata_mod_tmp(iel,ip,FLOOR(count_mod(ip)/2.0)+1)
              IF (ABS(rdata_mod_supob(iaz,ira,iel)) < eps_vr) THEN
                rdata_mod_supob(iaz,ira,iel) = 0.0_dp
              END IF
            END IF

            IF (count_obs(ip) >= supob_nrb .AND. rdata_obs_varsum(ip)/count_obs(ip) <= supob_vrw*supob_vrw) THEN
              rdata_obs_supob(iaz,ira,iel) = rdata_obs_tmp(iel,ip,FLOOR(count_obs(ip)/2.0)+1)
            END IF

          END IF

        END DO

      ELSE IF (TRIM(cvarname) == 'z') THEN

        naz_average = MAXVAL(cart_data(ista)%aw(iel,:))
        !  maximum number of raw radar bins for one superobing point:
        nct = (2*naz_average+1)*(2*cart_data(ista)%nra_average+1)

        ALLOCATE(rdata_obs_tmp(rs_meta(ista)%nel,ncploc,nct))
        ALLOCATE(rdata_mod_tmp(rs_meta(ista)%nel,ncploc,nct))
!CDIR BEGIN COLLAPSE
        rdata_obs_tmp   = miss_value
        rdata_mod_tmp   = miss_value
!CDIR END

        DO iaz_average = -naz_average, naz_average, 1
          DO ira_average = -cart_data(ista)%nra_average, cart_data(ista)%nra_average, 1
            DO ip = 1, np

              IF (cart_data(ista)%rdata_ind(iel,ip) > 0) THEN

!!$ UB:
!!$                aw_loc = MIN(cart_data(ista)%aw_max, cart_data(ista)%aw(iel,ip) )
                aw_loc = cart_data(ista)%aw(iel,ip)

                iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
                ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1
                ira_sub = ira + ira_average

                IF  (ABS(iaz_average) <= aw_loc .AND. &
                     ira_sub >= 1               .AND. &
                     ira_sub <=rs_meta(ista)%nra) THEN

                  iaz_sub = MODULO(iaz + iaz_average - 1, rs_meta(ista)%naz) + 1

                  IF (rdata_mod(iaz_sub,ira_sub,iel) > miss_threshold &
                       .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_z) THEN
                    count_mod(ip) =  count_mod(ip) + 1
                    IF (rdata_mod(iaz_sub,ira_sub,iel) >= dBZ_crit_radar) THEN
                      rdata_mod_sum(ip) = rdata_mod_sum(ip) + EXP(ln10_o10*rdata_mod(iaz_sub,ira_sub,iel))
                      rdata_mod_tmp(iel,ip,count_mod(ip)) = EXP(ln10_o10*rdata_mod(iaz_sub,ira_sub,iel))
                    END IF
                  END IF

                  IF (rdata_obs(iaz_sub,ira_sub,iel) > miss_threshold &
                       .AND. ira*rs_meta(ista)%ra_inc >= cart_data(ista)%minra_z) THEN
                    count_obs(ip) =  count_obs(ip) + 1
                    IF (rdata_obs(iaz_sub,ira_sub,iel) >= dBZ_crit_radar) THEN
                      rdata_obs_sum(ip) = rdata_obs_sum(ip) + EXP(ln10_o10*rdata_obs(iaz_sub,ira_sub,iel))
                      rdata_obs_tmp(iel,ip,count_obs(ip)) = EXP(ln10_o10*rdata_obs(iaz_sub,ira_sub,iel))
                    END IF
                  END IF

                END IF ! for ip with valid vaules

              END IF

            END DO  ! np
          END DO ! end of nra_average
        END DO ! end of naz_average

        DO ip = 1, np
          IF (count_mod(ip) >= supob_nrb) THEN
            rdata_mod_ave(ip) = rdata_mod_sum(ip)/count_mod(ip)
          END IF
          IF (count_obs(ip) >= supob_nrb) THEN
            rdata_obs_ave(ip) = rdata_obs_sum(ip)/count_obs(ip)
          END IF
        END DO

        DO ict = 1, nct
          DO jct = nct, ict+1, -1
            DO ip = 1, np

              IF (count_mod(ip) >= supob_nrb) THEN

                IF (rdata_mod_tmp(iel,ip,jct) < rdata_mod_tmp(iel,ip,jct-1) &
                     .AND. rdata_mod_tmp(iel,ip,jct) > miss_threshold) THEN

                  mod_tmp =rdata_mod_tmp(iel,ip,jct-1)
                  rdata_mod_tmp(iel,ip,jct-1) = rdata_mod_tmp(iel,ip,jct)
                  rdata_mod_tmp(iel,ip,jct) = mod_tmp

                END IF
              END IF

              IF (count_obs(ip) > supob_nrb) THEN

                IF (rdata_obs_tmp(iel,ip,jct) < rdata_obs_tmp(iel,ip,jct-1)  &
                     .AND. rdata_obs_tmp(iel,ip,jct) > miss_threshold) THEN

                  obs_tmp =rdata_obs_tmp(iel,ip,jct-1)
                  rdata_obs_tmp(iel,ip,jct-1) = rdata_obs_tmp(iel,ip,jct)
                  rdata_obs_tmp(iel,ip,jct) = obs_tmp

                END IF
              END IF

            END DO
          END DO
        END DO

!!$ UB: this can be made more efficient, like in averagesuperobing_vec()!
        DO ict = 1, nct
          DO ip = 1, np

            IF (count_mod(ip) >= supob_nrb .AND. ict <= count_mod(ip)) THEN
              rdata_mod_varsum(ip) =  rdata_mod_varsum(ip) + (rdata_mod_tmp(iel,ip,ict) - rdata_mod_ave(ip))**2
            END IF

            IF (count_obs(ip) >= supob_nrb .AND. ict <= count_obs(ip)) THEN
              rdata_obs_varsum(ip) =  rdata_obs_varsum(ip) + (rdata_obs_tmp(iel,ip,ict) - rdata_obs_ave(ip))**2
            END IF

          END DO
        END DO

        DO ip = 1, np

          IF (cart_data(ista)%rdata_ind(iel,ip) > 0) THEN

            iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
            ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

            IF (count_mod(ip) >= supob_nrb) THEN
              IF (rdata_mod_varsum(ip)/count_mod(ip) <= supob_rfl*supob_rfl) THEN
                IF (rdata_mod_tmp(iel,ip,FLOOR(count_mod(ip)/2.0)+1) > Z_crit_radar) THEN
                  rdata_mod_supob(iaz,ira,iel) = 10.0_dp*LOG10(rdata_mod_tmp(iel,ip,FLOOR(count_mod(ip)/2.0)+1))
                ELSE
                  rdata_mod_supob(iaz,ira,iel) = zero_value
                END IF
              END IF
            END IF

            IF (count_obs(ip) >= supob_nrb) THEN
              IF (rdata_obs_varsum(ip)/count_obs(ip) <= supob_rfl*supob_rfl) THEN
                IF (rdata_obs_tmp(iel,ip,FLOOR(count_obs(ip)/2.0)+1) > Z_crit_radar) THEN
                  rdata_obs_supob(iaz,ira,iel) = 10.0_dp*LOG10(rdata_obs_tmp(iel,ip,FLOOR(count_obs(ip)/2.0)+1))
                ELSE
                  rdata_obs_supob(iaz,ira,iel) = zero_value
                END IF
              END IF
            END IF

          END IF

        END DO

      ENDIF

      ! .. write the geographic (quasi-cartesian) coordinates to the INOUT fields lat_supob and lon_supob
      !     at the suberob'ed locations:
      DO ip = 1, np

        i = MOD((cart_data(ista)%cartdata_ind(iel,ip)-1),cart_data(ista)%ni_tot) + 1
        j = (cart_data(ista)%cartdata_ind(iel,ip)-1)/cart_data(ista)%ni_tot  + 1

        iaz = MOD((cart_data(ista)%rdata_ind(iel,ip)-1),rs_meta(ista)%naz) + 1
        ira = (cart_data(ista)%rdata_ind(iel,ip)-1)/rs_meta(ista)%naz + 1

        clon_r = cart_data(ista)%startlon + cart_data(ista)%width_deg + (i-1)*cart_data(ista)%dlon
        clat_r = cart_data(ista)%startlat + cart_data(ista)%width_deg + (j-1)*cart_data(ista)%dlat

        lat_supob(iaz,ira,iel) = phirot2phi(clat_r, clon_r, &
             cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)
        lon_supob(iaz,ira,iel) = rlarot2rla(clat_r, clon_r, &
             cart_data(ista)%pollat, cart_data(ista)%pollon, cart_data(ista)%polgam)

      END DO

    END DO ! end of nel

    DEALLOCATE(indloc)

    DEALLOCATE(count_obs)
    DEALLOCATE(count_mod)

    DEALLOCATE(rdata_mod_sum)
    DEALLOCATE(rdata_obs_sum)

    DEALLOCATE(rdata_mod_ave)
    DEALLOCATE(rdata_obs_ave)

    DEALLOCATE(rdata_mod_varsum)
    DEALLOCATE(rdata_obs_varsum)

    DEALLOCATE(rdata_obs_tmp)
    DEALLOCATE(rdata_mod_tmp)

    IF (ldebug_radsim) WRITE (*,*) 'Done with ', TRIM(yzroutine), ' on proc ', my_radar_id

  END SUBROUTINE mediansuperobing_radar_vec

END MODULE radar_output_country_obs
