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

MODULE radar_output_station

!------------------------------------------------------------------------------
!
! Description: Modules of the radar forward operator EMVORADO for processing
!              of the various output methods, data and formats:
!              volume data output, feedback file output and reflectivity composite
!              generation and output.
!              This module contains methods to collect/MPI-copy simulated radar data
!              from the compute PEs to the output PEs, for reading obs data files,
!              for producing superobservations. It uses methods from another module radar_output_methods.f90
!              (soon to be splitted from the present module!!!!)
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
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       i_fwo_bubbles,     & ! Timing flag
       i_fwo_composites,  & ! Timing flag
       i_fwo_out,         & ! Timing flag for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                            !  reading obs data, producing feedback files)
       degrad,            &
       rs_meta, cart_data, dbz_meta, &
       comp_meta, comp_meta_bub

  USE radar_composites, ONLY : &
       composite2D_dbz_maxmethod_ista,  &
       comp_dbzsim_tot,      &
       comp_dbzsim_bub_tot

  USE radar_interface, ONLY : &
       abort_run,             &
       get_runtime_timings, &
       get_obstime_ind_of_currtime

 !------------------------------------------------------------------------------

   USE radar_utilities, ONLY :   &
                               ind2sub3D, sub2ind3D,       &
                               init_vari,                  &
                               set_missing_and_correct0

  !------------------------------------------------------------------------------

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, loutradwind, &
       loutdbz, loutpolstd, loutpolall, lextdbz, lout_geom, &
       lweightdbz, lfall, lonline, lsode, lsmooth, lreadmeta_from_netcdf, &
       ldealiase_vr_obs, &
       lmds_z, lmds_vr, &
       lfill_vr_backgroundwind, &
       ldo_composite, lcomposite_output, &
       nel_composite, ldo_bubbles

#ifdef AUXOUT_OFFLINE
  USE radar_data_namelist, only: &
       loutvolaux
#endif

  USE radar_model2rays, ONLY : online2geo, rad2geo_const

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
  
  USE radar_output_country_obs, ONLY : output_radar_country_obs
  
!================================================================================
!================================================================================

  IMPLICIT NONE

!================================================================================
!================================================================================

  PRIVATE

  PUBLIC ::  output_my_ista

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
  !+ Module procedure in radar_src for output of one radar station on
  !  a specific PE. Is called from output_radar() within the parallelization loops.
  !------------------------------------------------------------------------------

  SUBROUTINE output_my_ista (time_mod, ista, &
                             nobs, radpos_all, &
                             vt_mod_all, radwind_mod_all,&
                             zh_radar_mod_all, ah_radar_mod_all, &
                             zv_radar_mod_all, rrhv_radar_mod_all, irhv_radar_mod_all, &
                             kdp_radar_mod_all, adp_radar_mod_all, zvh_radar_mod_all, &
                             hl_loc_all, el_loc_all, s_loc_all, &
                             geom_written)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    REAL (KIND=dp), INTENT(IN)     :: time_mod          ! seconds since model start
    INTEGER, INTENT(in)  :: ista, nobs
    INTEGER, POINTER  :: radpos_all(:) ! radar points indices of a radar station
    REAL    (KIND=dp), POINTER   :: vt_mod_all(:), radwind_mod_all(:), &
                                    zh_radar_mod_all(:), ah_radar_mod_all(:), &
                                    zv_radar_mod_all(:), &
                                    rrhv_radar_mod_all(:), irhv_radar_mod_all(:), &
                                    kdp_radar_mod_all(:), adp_radar_mod_all(:), &
                                    zvh_radar_mod_all(:), &
                                    hl_loc_all(:), el_loc_all(:), s_loc_all(:)
    LOGICAL, INTENT(inout)  :: geom_written

    !------------------------------------------------------------------------------
    !
    ! Local variables:

    CHARACTER(LEN=*), PARAMETER :: yzroutine = 'output_my_ista'

    INTEGER       :: m,n,o,iobs,i,j,k,izerror,itime
    INTEGER, ALLOCATABLE :: &
         m_all(:),       & ! azimuthal indices for all radar points of a radar station in simul
         n_all(:),       & ! radial    indices for all radar points of a radar station in simul
         o_all(:)          ! elevation indices for all radar points of a radar station in simul

    REAL(KIND=dp) :: ext_coeff !, maxext, minext

    REAL(KIND=dp), PARAMETER :: cext = 20.0/LOG(10.0)

    REAL(KIND=dp), ALLOCATABLE :: &
         vrpolar(:,:,:), zrpolar(:,:,:), &
#ifdef AUXOUT_OFFLINE
         zvpolar(:,:,:), zvhpolar(:,:,:), &
#endif
         zdrpolar(:,:,:), &
         rhvpolar(:,:,:), &
         kdppolar(:,:,:), phidppolar(:,:,:), &
         ldrpolar(:,:,:), &
         zepolar(:,:,:), zetpolar(:,:,:), &
         zedpolar(:,:,:), zedtpolar(:,:,:), &
         hrpolar(:,:,:), erpolar(:,:,:), srpolar(:,:,:), &
         lonpolar(:,:,:), latpolar(:,:,:), &
         mds_dbz(:), hdummy(:,:), &
         vrpolar_for_dealiasing(:,:,:)

    CHARACTER(len=100) :: fileprefix


    IF (ldebug_radsim) THEN
      WRITE (*,'(a,i0,a,i6.6,a,i0)') TRIM(yzroutine)//': output of radar station ', &
           ista, ' (ID: ', rs_meta(ista)%station_id, ' on proc ', my_radar_id
    END IF

    itime = get_obstime_ind_of_currtime ( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) ) ! missval_int if not found

    ! .. Reset the list of present elevations to the list of nominal elevations:
    rs_meta(ista)%nel_present = rs_meta(ista)%nel
    rs_meta(ista)%ind_ele_present(:) = missval_int
    rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) = (/ (i, i=1, rs_meta(ista)%nel_present) /)
    
    ! .. Allocate fields for 3D radar data:
    !     (allocate in any case because these are arguments to a subroutine)
    ALLOCATE(vrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
             zrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
    ! .. Initialize with standard missing value (= point outside model domain):
    CALL init_vari(vrpolar, miss_value)
    CALL init_vari(zrpolar, miss_value)

    IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
      ALLOCATE(zdrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               rhvpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               kdppolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               phidppolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
      CALL init_vari(zdrpolar, miss_value)
      CALL init_vari(rhvpolar, miss_value)
      CALL init_vari(kdppolar, miss_value)
      CALL init_vari(phidppolar, 0.0_dp)
#ifdef AUXOUT_OFFLINE
      IF (loutvolaux) THEN
        ALLOCATE(zvpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(zvpolar, miss_value)
      END IF
#endif
      IF (loutpolall) THEN
        ALLOCATE(ldrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(ldrpolar, miss_value)
#ifdef AUXOUT_OFFLINE
        IF (loutvolaux) THEN
          ALLOCATE(zvhpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
          CALL init_vari(zvhpolar, miss_value)
        END IF
#endif
      END IF
    END IF

    IF (lextdbz .AND. &
        (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
      ALLOCATE(zepolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               zetpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
      CALL init_vari(zepolar, miss_value)
      CALL init_vari(zetpolar, 0.0_dp)
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(zedpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
                 zedtpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(zedpolar, miss_value)
        CALL init_vari(zedtpolar, 0.0_dp)
      END IF
    END IF

    ! .. Allocate index fields for re-copying radar data vectors to 3D data sets:
    ALLOCATE(m_all(nobs), n_all(nobs), o_all(nobs))
    CALL init_vari(m_all, -1)
    CALL init_vari(n_all, -1)
    CALL init_vari(o_all, -1)

!$omp parallel do private(iobs)
    DO iobs = 1, nobs
      CALL ind2sub3D(radpos_all(iobs), rs_meta(ista)%naz, rs_meta(ista)%nra, &
           m_all(iobs), n_all(iobs), o_all(iobs))
    END DO
!$omp end parallel do

    IF ( ANY( n_all > rs_meta(ista)%nra .OR. n_all < 1 .OR. &
              m_all > rs_meta(ista)%naz .OR. m_all < 1 .OR. &
              o_all > rs_meta(ista)%nel .OR. o_all < 1 ) ) THEN
      WRITE (*,*) TRIM(yzroutine)//' ERROR: SO NE SCH***!!!'
      WRITE (*,'(a,T70,6a12)') ' ERROR debug list radpos_all', 'nra', 'naz', 'nel', 'n_all', 'm_all', 'o_all'
      DO iobs=1, nobs
        IF (  n_all(iobs) > rs_meta(ista)%nra .OR. n_all(iobs) < 1 .OR. &
              m_all(iobs) > rs_meta(ista)%naz .OR. m_all(iobs) < 1 .OR. &
              o_all(iobs) > rs_meta(ista)%nel .OR. o_all(iobs) < 1 ) THEN
          WRITE (*,'(a,i5,a,i6.6,a,i8,a,i12,a,6i12)') ' ERROR proc=', my_radar_id, ' station=', rs_meta(ista)%station_id, &
                     ' radpos_all(', iobs,') = ', radpos_all(iobs), ':', &
                     rs_meta(ista)%nra, rs_meta(ista)%naz, rs_meta(ista)%nel, &
                     n_all(iobs), m_all(iobs), o_all(iobs)
        END IF
      END DO
    END IF


    IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written)) .OR. lreadmeta_from_netcdf ) THEN

      ! .. Store height and local elevation on polar fields:

      ALLOCATE(hrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               erpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
      ! .. Initialize with standard missing value (= point outside model domain):

      CALL init_vari(hrpolar, miss_value)
      CALL init_vari(erpolar, miss_value)

!$omp parallel do private(iobs)
      DO iobs = 1, nobs
        hrpolar( m_all(iobs),&
                 n_all(iobs),&
                 o_all(iobs) ) = hl_loc_all(iobs)
        erpolar( m_all(iobs),&
                 n_all(iobs),&
                 o_all(iobs) ) = el_loc_all(iobs)
      END DO
!$omp end parallel do

      ! .. Compute lon/lat of radar points in a way that there are no missing values
      !     in the resulting lat/lon polar arrays (i.e., all coordinates are defined,
      !     even if the data are missing values):

      ALLOCATE(lonpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               latpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
      CALL init_vari(lonpolar, miss_value)
      CALL init_vari(latpolar, miss_value)

      IF (lonline) THEN

        ALLOCATE(srpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(srpolar, miss_value)

!$omp parallel do private(iobs)
        DO iobs = 1, nobs
          srpolar( m_all(iobs),&
                   n_all(iobs),&
                   o_all(iobs) ) = s_loc_all(iobs)
        END DO
!$omp end parallel do

        ! .. lat/lon are computed from rs_data(ista)%s_loc, and for blocked/missing ray
        !    parts, the 4/3 earth radius model is used for extrapolation (srpolar, hrpolar are filled):
        CALL online2geo(rs_meta(ista), srpolar, hrpolar, erpolar, latpolar, lonpolar)

      ELSE

        IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
          ALLOCATE(hdummy(rs_meta(ista)%naz,rs_meta(ista)%nra))
          CALL rad2geo_const (rs_meta(ista), latpolar, lonpolar, hdummy)
          ! .. Fill missing values in hrpolar (above model top) with values from hdummy
          !     (should be continuous along rays, because for both the same 4/3 earth radius model has been applied):
!$omp parallel do collapse(3) private(i,j,k)
          DO i = 1, rs_meta(ista)%nel
            DO j = 1, rs_meta(ista)%nra
              DO k = 1, rs_meta(ista)%naz
                IF (hrpolar(k,j,i) < miss_threshold) THEN
                  hrpolar(k,j,i) = hdummy(k,j)
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do
        ELSE
          ALLOCATE(hdummy(rs_meta(ista)%nra,rs_meta(ista)%nel))
          CALL rad2geo_const (rs_meta(ista), latpolar, lonpolar, hdummy)
          ! .. Fill missing values in hrpolar (above model top) with values from hdummy
          !     (should be continuous along rays, because for both the same 4/3 earth radius model has been applied):
!$omp parallel do collapse(3) private(i,j,k)
          DO i = 1, rs_meta(ista)%nel
            DO j = 1, rs_meta(ista)%nra
              DO k = 1, rs_meta(ista)%naz
                IF (hrpolar(k,j,i) < miss_threshold) THEN
                  hrpolar(k,j,i) = hdummy(j,i)
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do
        END IF

        DEALLOCATE(hdummy)

      END IF


      IF ( lout_geom .AND. (lonline .OR. .NOT.geom_written) ) THEN

        IF (.NOT. lonline) THEN
          ! Zero out shielded pixels in case of 4/3 earth model. For online ray propagation, this is not necessary
          !   because the ray tracing ends on impact on orography!
          CALL elim_shielded(erpolar, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, miss_value)
        END IF

        !.. Output the sorted polar data set into a standard file format:
        !   - binary on the SX9, convert to simple ASCII with program "bin2ascii_convrates3d" from Ulrich Blahak
        !   - simple ASCII on all other systems

        IF ( lonline ) THEN
          CALL get_fileprefix_ascii_output ('adsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               srpolar, TRIM(ADJUSTL(fileprefix)), &
               'Simul. arc distance at MSL from radar station', 'm', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

          IF (ldebug_radsim) THEN
            CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "arc_dist [m]", &
                                srpolar, (srpolar > miss_threshold), miss_value, nobs)
          END IF
        END IF

        CALL get_fileprefix_ascii_output ('losim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             lonpolar, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radar bin geographic longitude', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "lon_polar [deg]", &
                              lonpolar, (lonpolar > miss_threshold), miss_value, nobs)
        END IF

        CALL get_fileprefix_ascii_output ('lasim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             latpolar, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radar bin geographic latitude', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "lat_polar [deg]", &
                              latpolar, (latpolar > miss_threshold), miss_value, nobs)
        END IF

        CALL get_fileprefix_ascii_output ('hrsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             hrpolar, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radar bin height above MSL', 'm', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "hl_polar [m]", &
                              hrpolar, (hrpolar > miss_threshold), miss_value, nobs)
        END IF

        CALL get_fileprefix_ascii_output ('ersim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             erpolar, TRIM(ADJUSTL(fileprefix)), &
             'Simul. local elevation angle', 'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "el_polar [deg]", &
                              erpolar, (erpolar > miss_threshold), miss_value, nobs)
        END IF

        ! Set flag that geometry for this radar has been written to a file:
        geom_written = .TRUE.

      END IF

    END IF

    IF (lmds_z .OR. lmds_vr) THEN
      ! Minimum detectable signal as function of range in dBZ:
      ALLOCATE( mds_dbz(rs_meta(ista)%nra) )
!$omp parallel do private(iobs)
      DO iobs = 1, rs_meta(ista)%nra
        mds_dbz(iobs) = rs_meta(ista)%mds_Z0 + 20.0_dp * LOG10(iobs*rs_meta(ista)%ra_inc/rs_meta(ista)%mds_r0)
      END DO
!$omp end parallel do
    END IF

    IF (loutdbz .OR. (loutradwind .AND. lmds_vr)) THEN
      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:

!$omp parallel do private(iobs)
      DO iobs = 1, nobs
        ! Convert reflectivity from linear to dBZ values:
        IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
          zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 10.0_dp * LOG10(zh_radar_mod_all(iobs))
        ELSEIF (zh_radar_mod_all(iobs) >= miss_threshold) THEN
          zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
        ELSEIF (zh_radar_mod_all(iobs) >= shield_low_threshold) THEN
          ! Points below the surface:
          zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zh_radar_mod_all(iobs)
        END IF
        ! in other cases miss_value or shield_value (below model orography) ...
      END DO
!$omp end parallel do

      IF (.NOT. lonline) THEN
        ! Zero out shielded pixels in case of 4/3 earth model.
        !   For online ray propagation, this is not necessary
        !   because the ray tracing ends on impact on orography!
        CALL elim_shielded(zrpolar, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             miss_value)
      END IF

#ifdef AUXOUT_OFFLINE
      IF (loutvolaux) THEN
        ! Output unattenuated ZH
        IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          CALL get_fileprefix_ascii_output ('zruasim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               zrpolar, TRIM(ADJUSTL(fileprefix)), &
               'Simul. unattenuated radar reflectivity', &
               'dBZ', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF
      END IF
#endif

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

        ! Note: Output to file for all of them only after lextdbz loop (where
        !       also zrpolar is output), since some are affected by attenuation!

#ifdef AUXOUT_OFFLINE
        IF (loutvolaux) THEN
          ! Calculate ZV explicitly
!$omp parallel do private(iobs)
          DO iobs = 1, nobs
            IF (zv_radar_mod_all(iobs) >= Z_crit_radar) THEN
              zvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 10.0_dp * LOG10(zv_radar_mod_all(iobs))
            ELSEIF (zv_radar_mod_all(iobs) >= miss_threshold) THEN
              zvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
            ELSEIF (zh_radar_mod_all(iobs) >= shield_low_threshold) THEN
              zvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zv_radar_mod_all(iobs)
            END IF
          END DO
!$omp end parallel do

          IF (.NOT. lonline) THEN
            CALL elim_shielded(zvpolar, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 miss_value)
          END IF

          ! Output unattenuated ZV
          IF (lextdbz .AND. &
              (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
            CALL get_fileprefix_ascii_output ('zvuasim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 zvpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Simul. unattenuated vertically polarized radar reflectivity', &
                 'dBZ', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          END IF
        END IF
#endif

        !..ZDR
!$omp parallel do private(iobs)
        DO iobs = 1, nobs
          ! Calculate zdr=zh/zv and convert to dB values
          IF (zv_radar_mod_all(iobs) >= Z_crit_radar) THEN
            IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
              zdrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
                   zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) - &
                   10.0_dp * LOG10(zv_radar_mod_all(iobs))
            ELSEIF (zh_radar_mod_all(iobs) >= miss_threshold) THEN
              zdrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
            ELSEIF (zh_radar_mod_all(iobs) >= shield_low_threshold) THEN
              ! Points below the surface:
              zdrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zh_radar_mod_all(iobs)
            END IF
          ! ZV too small - No need to reset, already initialzed like that
          !ELSE
          !  zdrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = miss_value
          END IF
          ! in other cases miss_value or shield_value (below model orography) ...
        END DO
!$omp end parallel do

        IF (.NOT. lonline) THEN
          CALL elim_shielded(zdrpolar, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               miss_value)
        END IF

#ifdef AUXOUT_OFFLINE
        IF (loutvolaux) THEN
          ! Output unattenuated ZDR
          IF (lextdbz .AND. &
              (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
            CALL get_fileprefix_ascii_output ('zdruasim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 zdrpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Simul. unattenuated Differential Reflectivity', &
                 'dB', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          END IF
        END IF
#endif

        !..RHV
!$omp parallel do private(iobs)
        DO iobs = 1, nobs
          ! Calculate rhv = SQRT( (rrhv^2+irhv^2) / (zh*zv) )
          IF ((zv_radar_mod_all(iobs) >= Z_crit_radar) .AND. &
              (zh_radar_mod_all(iobs) >= Z_crit_radar)) THEN
            IF ((rrhv_radar_mod_all(iobs) >= miss_thresh_rhv) .AND. &
                (irhv_radar_mod_all(iobs) >= miss_thresh_rhv)) THEN
              rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
                   rrhv_radar_mod_all(iobs)**2 + irhv_radar_mod_all(iobs)**2
              IF (rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) >= Z_crit_radar**2) THEN
                rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
                     SQRT( rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) / &
                           (zv_radar_mod_all(iobs) * zh_radar_mod_all(iobs)) )
              ELSE
                rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
              END IF
            ELSEIF ((rrhv_radar_mod_all(iobs) >= shield_low_thresh_rhv) .OR. &
                    (irhv_radar_mod_all(iobs) >= shield_low_thresh_rhv)) THEN
!              rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
!                   MIN(rrhv_radar_mod_all(iobs),irhv_radar_mod_all(iobs))
              rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = shield_value
            END IF
          ! ZH or ZV too small - No need to reset, already initialzed like that
          !ELSE
          !  rhvpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = miss_value
          END IF
        END DO
!$omp end parallel do

        IF (.NOT. lonline) THEN
          CALL elim_shielded(rhvpolar, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               miss_value)
        END IF

        !..KDP
!$omp parallel do private(iobs)
        DO iobs = 1, nobs
          kdppolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = kdp_radar_mod_all(iobs)
        END DO
!$omp end parallel do

        IF (.NOT. lonline) THEN
          CALL elim_shielded(kdppolar, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               miss_value)
        END IF

!$omp parallel do private(i,j,k)
        DO i = 1, rs_meta(ista)%nel
          DO j = 1, rs_meta(ista)%nra
            DO k = 1, rs_meta(ista)%naz
              ! Clip missing values miss_value resp. shield_value to 0.0
              ! Can't use the MAX(field,0) as for ah, since kdp can be negative
!!$ UB: the previous formulation would have been wrong if a miss_value appears along a ray, because then the j-chain would be broken.
!!$     normally this does not really happen for simulated data, but we correct it anyways.
              IF (kdppolar(k,j,i) >= miss_threshold) THEN
                ext_coeff = 2.0_dp * kdppolar(k,j,i)
                ! Convert KDP from deg/m to deg/km
                kdppolar(k,j,i) = 1000.0_dp * kdppolar(k,j,i)
              ELSE
                ext_coeff = 0.0_dp
              END IF
              ! Sum up total phase shift:
              phidppolar(k,j,i) = phidppolar(k,MAX(j-1,1),i) + ext_coeff * rs_meta(ista)%ra_inc
            END DO
          END DO
        END DO
!$omp end parallel do

        ! Eliminate phidppolar where kdppolar is missing or shielded:
        WHERE (kdppolar < miss_threshold) phidppolar = miss_value 
        
        ! FIXME:
        ! Should a phase wrapping (to [0..360] or [-180..180]deg) be applied?
        ! In order to compare to obs that seems reasonable.
        ! But so far, we don't.
        !
        ! YES it should: the obs of DWD are from -180 to 180

        IF (loutpolall) THEN
          !..LDR
!$omp parallel do private(iobs)
          DO iobs = 1, nobs
            !! Calculate ldr=zvh/zh and convert to dB values
            !IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
            !  IF (zvh_radar_mod_all(iobs) >= Z_crit_radar) THEN
            !    ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
            !         10.0_dp * LOG10(zvh_radar_mod_all(iobs)) - &
            !         zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) )
            !  ELSEIF (zvh_radar_mod_all(iobs) >= miss_threshold) THEN
            !    ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
            !  ELSEIF (zvh_radar_mod_all(iobs) >= shield_low_threshold) THEN
            !    ! Points below the surface:
            !    ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zvh_radar_mod_all(iobs)
            !  END IF
            !ELSE ! ZH too small
            !  ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = miss_value
            !END IF
            !! in other cases miss_value or shield_value (below model orography) ...

            ! JM201025:
            ! Thought to output LDR in dBZ to be able to see small values, but
            ! confuses me, particularly when checking qMie results (provides
            ! values -30..-90, which is likely correct considering that lin Ldr
            ! should be 0). Therefore changing (back & for now) to linear
            ! output.
            IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
              IF (zvh_radar_mod_all(iobs) >= Z_crit_radar) THEN
                ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
                     zvh_radar_mod_all(iobs) / zh_radar_mod_all(iobs)
              ELSEIF (zvh_radar_mod_all(iobs) >= miss_threshold) THEN
                ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 0.0_dp
              ELSEIF (zvh_radar_mod_all(iobs) >= shield_low_threshold) THEN
                ! Points below the surface:
                ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zvh_radar_mod_all(iobs)
              END IF
            ELSE ! ZH too small
              ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = miss_value
            END IF
          END DO
!$omp end parallel do

          IF (.NOT. lonline) THEN
            CALL elim_shielded(ldrpolar, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 miss_value)
          END IF

#ifdef AUXOUT_OFFLINE
          IF (loutvolaux) THEN
            ! Output unattenuated LDR
            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              CALL get_fileprefix_ascii_output ('ldruasim', rs_meta(ista), fileprefix)
              CALL output3d_ascii_radar(itime, &
                   rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                   ldrpolar, TRIM(ADJUSTL(fileprefix)), &
                   'Simul. unattenuated Linear Depolarization Ratio', &
                   !'[dB]', 'polar', &
                   '[-]', 'polar', &
                   rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
            END IF

            ! Calculate ZVH explicitly
!$omp parallel do private(iobs)
            DO iobs = 1, nobs
              IF (zvh_radar_mod_all(iobs) >= Z_crit_radar) THEN
                zvhpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 10.0_dp * LOG10(zvh_radar_mod_all(iobs))
              ELSEIF (zv_radar_mod_all(iobs) >= miss_threshold) THEN
                zvhpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
              ELSEIF (zh_radar_mod_all(iobs) >= shield_low_threshold) THEN
                zvhpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zvh_radar_mod_all(iobs)
              END IF
            END DO
!$omp end parallel do

            IF (.NOT. lonline) THEN
              CALL elim_shielded(zvhpolar, &
                   rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                   miss_value)
            END IF

            ! Output unattenuated ZVH
            IF (lextdbz .AND. &
                (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              CALL get_fileprefix_ascii_output ('zvhuasim', rs_meta(ista), fileprefix)
              CALL output3d_ascii_radar(itime, &
                   rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                   zvhpolar, TRIM(ADJUSTL(fileprefix)), &
                   'Simul. unattenuated H-turned-V polarized radar reflectivity', &
                   'dBZ', 'polar', &
                   rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
            END IF
          END IF
#endif

        END IF ! loutpolall

      END IF ! loutpolstd

      !..Attenuation and attenuation effects on other parameters
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN

!$omp parallel do private(iobs)
        DO iobs = 1, nobs
          zepolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = ah_radar_mod_all(iobs)
        END DO
!$omp end parallel do

        IF (.NOT. lonline) THEN
          ! Zero out shielded pixels in case of 4/3 earth model.
          !   For online ray propagation, this is not necessary
          !   because the ray tracing ends on impact on orography!
          CALL elim_shielded(zepolar, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               miss_value)
        END IF

!!$ no collapse(3) because of j-1 dependence!
!$omp parallel do private(i,j,k,ext_coeff)
        DO i = 1, rs_meta(ista)%nel
          DO j = 1, rs_meta(ista)%nra
            DO k = 1, rs_meta(ista)%naz
              ! Clip missing values miss_value resp. shield_value to 0.0:
              ext_coeff = cext * MAX(zepolar(k,j,i), 0.0_dp)  ! two-way [dB/m]
              ! Sum up total extinction:
              zetpolar(k,j,i) = zetpolar(k,MAX(j-1,1),i) + ext_coeff * rs_meta(ista)%ra_inc
              IF (zrpolar(k,j,i) >= miss_threshold) &
                zrpolar(k,j,i) = zrpolar(k,j,i) - zetpolar(k,j,i)
            END DO
          END DO
        END DO
!$omp end parallel do

        ! Set correct missing values for zepolar and convert units:
!!$        WHERE (zepolar >= miss_threshold)
!!$          zepolar = 1000.0_dp * cext * zepolar(k,j,i)  ! unit is now dB/km
!!$        ELSEWHERE
!!$          zepolar = miss_value
!!$        END WHERE
!$omp parallel do private(i,j,k) collapse(3)
        DO i = 1, rs_meta(ista)%nel
          DO j = 1, rs_meta(ista)%nra
            DO k = 1, rs_meta(ista)%naz
              IF (zepolar(k,j,i) >= miss_threshold) THEN
                zepolar(k,j,i) = 1000.0_dp * cext * zepolar(k,j,i)  ! unit is now dB/km
              ELSE
                zepolar(k,j,i) = miss_value
              END IF
            END DO
          END DO
        END DO
!$omp end parallel do

        ! Output to file(s):
        IF (loutdbz) THEN
          CALL get_fileprefix_ascii_output ('ahsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               zepolar, TRIM(ADJUSTL(fileprefix)), &
               'Twoway attenuation coefficient for reflectivity', &
               'dB/km', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF

        IF (ldebug_radsim) THEN
          CALL control_output(&
               ahoutputunit(ista), time_mod, rs_meta(ista), &
               "k_2 [dB/km]", &
               zepolar, (zepolar > miss_threshold), miss_value, nobs)
        END IF

        IF (loutdbz) THEN
          CALL get_fileprefix_ascii_output ('ahpisim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               zetpolar, TRIM(ADJUSTL(fileprefix)), &
               'Path integrated attenuation for reflectivity', &
               'dB', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

          !..ADP
!$omp parallel do private(iobs)
          DO iobs = 1, nobs
            zedpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = adp_radar_mod_all(iobs)
          END DO
!$omp end parallel do


          IF (.NOT. lonline) THEN
            CALL elim_shielded(zedpolar, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 miss_value)
          END IF

!$omp parallel do private(i,j,k,ext_coeff)
          DO i = 1, rs_meta(ista)%nel
            DO j = 1, rs_meta(ista)%nra
              DO k = 1, rs_meta(ista)%naz
                ! Clip missing values miss_value resp. shield_value to 0.0
                ! Can't use the MAX(field,0) as for ah, since adp can be negative
                IF (zedpolar(k,j,i) >= miss_threshold) THEN
                  ext_coeff = cext * zedpolar(k,j,i)  ! two-way [dB/m]
                ELSE
                  ext_coeff = 0.0_dp
                END IF
                ! Sum up total extinction:
                zedtpolar(k,j,i) = zedtpolar(k,MAX(j-1,1),i) + ext_coeff * rs_meta(ista)%ra_inc
                !..attenuated ZDR
                IF (zdrpolar(k,j,i) >= miss_threshold) THEN
                  zdrpolar(k,j,i) = zdrpolar(k,j,i) - zedtpolar(k,j,i)
                END IF
#ifdef AUXOUT_OFFLINE
                IF (loutvolaux) THEN
                  !..attenuated ZV
                  IF (zvpolar(k,j,i) >= miss_threshold) THEN
                    zvpolar(k,j,i) = zvpolar(k,j,i) - (zetpolar(k,j,i)-zedtpolar(k,j,i))
                  END IF
                  IF (loutpolall) THEN
                    !..attenuated ZVH
                    IF (zvhpolar(k,j,i) >= miss_threshold) THEN
                      zvhpolar(k,j,i) = zvhpolar(k,j,i) - (zetpolar(k,j,i)-0.5*zedtpolar(k,j,i))
                    END IF
                  END IF
                END IF
#endif
              END DO
            END DO
          END DO
!$omp end parallel do

          ! Eliminate zedtpolar where zedpolar is missing or shielded:
          WHERE (zedpolar < miss_threshold) zedtpolar = miss_value

        ! Set correct missing values for zedpolar and convert units:
!$omp parallel do private(i,j,k) collapse(3)
          DO i = 1, rs_meta(ista)%nel
            DO j = 1, rs_meta(ista)%nra
              DO k = 1, rs_meta(ista)%naz
                IF (zedpolar(k,j,i) >= miss_threshold) THEN
                  zedpolar(k,j,i) = 1000.0_dp * cext * zedpolar(k,j,i)  ! unit is now dB/km
                ELSE
                  zedpolar(k,j,i) = miss_value
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do

          ! Output to file(s):
          IF (loutdbz) THEN
            CALL get_fileprefix_ascii_output ('adpsim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 zedpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Twoway differential attenuation coefficient for reflectivity', &
                 'dB/km', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          END IF

          IF (ldebug_radsim) THEN
            CALL control_output(&
                 adpoutputunit(ista), time_mod, rs_meta(ista), &
                 "dk_2 [dB/km]", &
                 zedpolar, (zedpolar > miss_threshold), miss_value, nobs)
          END IF

          IF (loutdbz) THEN
            CALL get_fileprefix_ascii_output ('adppisim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 zedtpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Path integrated differential attenuation for reflectivity', &
                 'dB', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          END IF

          IF (loutpolall) THEN
            !..attenuated LDR
!$omp parallel do private(i,j,k,ext_coeff)
            DO i = 1, rs_meta(ista)%nel
              DO j = 1, rs_meta(ista)%nra
                DO k = 1, rs_meta(ista)%naz
                  IF (ldrpolar(k,j,i) >= miss_threshold) THEN
                    !ldrpolar(k,j,i) = ldrpolar(k,j,i) + 0.5_dp * zedtpolar(k,j,i)  ! log-space LDR
                    ldrpolar(k,j,i) = ldrpolar(k,j,i) * EXP(zedtpolar(k,j,i)/cext) ! for lin-space Ldr
                  END IF
                END DO
              END DO
            END DO
!$omp end parallel do
          END IF

        END IF   ! loutpol

      END IF   ! lextdbz


      ! Set correct 0-value (zero_value) and missing value (miss_value) for zrpolar:
!!$      WHERE (zrpolar >= miss_threshold .AND. zrpolar < dBZ_crit_radar)
!!$        zrpolar = zero_value
!!$      END WHERE
!!$      WHERE (zrpolar < miss_threshold)
!!$        zrpolar = miss_value
!!$      END WHERE
      CALL set_missing_and_correct0(zrpolar)

      IF (lmds_z) THEN
        ! Take into account the minimum detectable signal:
        ! (for simplicity assumed equal to that of the reflectivity)
!$omp parallel do private(o,n,m)
        DO o=1, rs_meta(ista)%nel
          DO n=1, rs_meta(ista)%nra
            DO m=1, rs_meta(ista)%naz
              IF (zrpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                ! Z is below the MDS, so set zrpolar to a "correct 0":
                zrpolar(m,n,o) = zero_value
              END IF
            END DO
          END DO
        END DO
!$omp end parallel do
      END IF

      IF (ldebug_radsim) THEN
        CALL control_output(&
             radrefloutputunit(ista), time_mod, rs_meta(ista), &
             "Z [dBZ]", &
             zrpolar, (zrpolar > miss_threshold), miss_value, nobs)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

#ifdef AUXOUT_OFFLINE
        IF (loutvolaux) THEN
          !..Output attenuated ZV
          CALL set_missing_and_correct0(zvpolar)
          IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
            DO o=1, rs_meta(ista)%nel
              DO n=1, rs_meta(ista)%nra
                DO m=1, rs_meta(ista)%naz
                  ! Minimum detectable signal masking should be consistent over
                  ! all radar moments, so use zrpolar instead of zvpolar:
                  IF (zvpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                    zvpolar(m,n,o) = zero_value
                  END IF
                END DO
              END DO
            END DO
!$omp end parallel do
          END IF
          CALL get_fileprefix_ascii_output ('zvsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               zvpolar, TRIM(ADJUSTL(fileprefix)), &
               'Simul. vertically polarized radar reflectivity', &
               'dBZ', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
        END IF
#endif

        ! FIXME:
        ! Should those get a minimum signal adjustments like zrpolar?
        ! Correct-0, at least, should not apply here as lin-ZDR = 0 for Zh = 0,
        ! which should be as invalid as Zv=0, hence both cases (actually any
        ! Z < Z_crit) makes lin-ZDR = missing_value

        !..ZDR
        !WHERE (zdrpolar < miss_threshold)
        !  zdrpolar = miss_value
        !END WHERE
        IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
          DO o=1, rs_meta(ista)%nel
            DO n=1, rs_meta(ista)%nra
              DO m=1, rs_meta(ista)%naz
                ! Minimum detectable signal masking should be consistent over
                ! all radar moments, so use zrpolar:
                IF (zdrpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                  zdrpolar(m,n,o) = miss_value
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do
        END IF
        CALL get_fileprefix_ascii_output ('zdrsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zdrpolar, TRIM(ADJUSTL(fileprefix)), &
             'Differential Reflectivity', &
             'dB', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(&
               zdroutputunit(ista), time_mod, rs_meta(ista), &
               "ZDR [dB]", &
               zdrpolar, (zdrpolar > miss_threshold), miss_value, nobs)
        END IF

        !..RHV
        IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
          DO o=1, rs_meta(ista)%nel
            DO n=1, rs_meta(ista)%nra
              DO m=1, rs_meta(ista)%naz
                ! Minimum detectable signal masking should be consistent over
                ! all radar moments, so use zrpolar:
                IF (rhvpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                  rhvpolar(m,n,o) = miss_value
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do
        END IF
        CALL get_fileprefix_ascii_output ('rhvsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             rhvpolar, TRIM(ADJUSTL(fileprefix)), &
             'Co-Polar Correlation Coefficient', &
             '[]', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(&
               rhvoutputunit(ista), time_mod, rs_meta(ista), &
               "RhoHV []", &
               rhvpolar, (rhvpolar > miss_threshold), miss_value, nobs)
        END IF

        !..KDP
        IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
          DO o=1, rs_meta(ista)%nel
            DO n=1, rs_meta(ista)%nra
              DO m=1, rs_meta(ista)%naz
                ! Minimum detectable signal masking should be consistent over:
                ! all radar moments, so use zrpolar:
                IF (kdppolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                  kdppolar(m,n,o) = miss_value
                END IF
              END DO
            END DO
          END DO
!$omp end parallel do
        END IF
        CALL get_fileprefix_ascii_output ('kdpsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             kdppolar, TRIM(ADJUSTL(fileprefix)), &
             'Specific differential phase', &
             'deg/km', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(&
               kdpoutputunit(ista), time_mod, rs_meta(ista), &
               "KDP [deg/km]", &
               kdppolar, (kdppolar > miss_threshold), miss_value, nobs)
        END IF

        !..PHIDP (as path integrated quantity, without lmds_z application)
        CALL get_fileprefix_ascii_output ('phidpsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             phidppolar, TRIM(ADJUSTL(fileprefix)), &
             'Differential propagation phase', &
             'deg', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (loutpolall) THEN
          !..LDR
          IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
            DO o=1, rs_meta(ista)%nel
              DO n=1, rs_meta(ista)%nra
                DO m=1, rs_meta(ista)%naz
                  ! Minimum detectable signal masking should be consistent over
                  ! all radar moments, so use zrpolar:
                  IF (ldrpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                    ldrpolar(m,n,o) = miss_value
                  END IF
                END DO
              END DO
            END DO
!$omp end parallel do
          END IF
          CALL get_fileprefix_ascii_output ('ldrsim', rs_meta(ista), fileprefix)
          CALL output3d_ascii_radar(itime, &
               rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
               ldrpolar, TRIM(ADJUSTL(fileprefix)), &
               'Linear Depolarization Ratio', &
!!$               '[dB]', 'polar', &
               '[-]', 'polar', &
               rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

          IF (ldebug_radsim) THEN
            CALL control_output(&
                 ldroutputunit(ista), time_mod, rs_meta(ista), &
                 "LDR [-]", &
                 ldrpolar, (ldrpolar > miss_threshold), miss_value, nobs)
          END IF

#ifdef AUXOUT_OFFLINE
          IF (loutvolaux) THEN
            !..Output attenuated ZVH
            CALL set_missing_and_correct0(zvhpolar)
            IF (lmds_z) THEN
!$omp parallel do private(o,n,m)
              DO o=1, rs_meta(ista)%nel
                DO n=1, rs_meta(ista)%nra
                  DO m=1, rs_meta(ista)%naz
                    ! Minimum detectable signal masking should be consistent over
                    ! all radar moments, so use zrpolar:
                    IF (zvhpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n)) THEN
                      zvhpolar(m,n,o) = zero_value
                    END IF
                  END DO
                END DO
              END DO
!$omp end parallel do
            END IF
            CALL get_fileprefix_ascii_output ('zvhsim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 zvhpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Simul. unattenuated H-turned-V polarized radar reflectivity', &
                 'dBZ', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
          END IF
#endif
        END IF   ! loutpolall

      END IF   ! loutpol
    END IF   ! loutdbz .OR. (loutradwind .AND. lmds_vr)

    IF (loutdbz) THEN
      
      !.. Output the sorted polar data set into a standard file format:
      !   - binary on the SX9, convert to simple ASCII with program
      !     "bin2ascii_convrates3d" from Ulrich Blahak
      !   - simple ASCII on all other systems
      CALL get_fileprefix_ascii_output ('zrsim', rs_meta(ista), fileprefix)
      CALL output3d_ascii_radar(itime, &
           rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
           zrpolar, TRIM(ADJUSTL(fileprefix)), &
           'Simul. radar reflectivity', &
           'dBZ', 'polar', &
           rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

    END IF
    
    IF (loutradwind) THEN

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!CDIR NODEP,VOVERTAKE,VOB
!$omp parallel do private(iobs)
      DO iobs = 1, nobs

        IF (ABS(radwind_mod_all(iobs)) < eps_vr) THEN
          radwind_mod_all(iobs) = 0.0_dp
        END IF

        IF (lfall .AND. radwind_mod_all(iobs) >= miss_threshold) THEN
          ! only impose terminal fall velocity if valid value and not below orography:
          vrpolar( m_all(iobs),&
                    n_all(iobs),&
                    o_all(iobs) ) = radwind_mod_all(iobs) - SIN(el_loc_all(iobs)*degrad)*vt_mod_all(iobs)
        ELSE
          ! otherwise, just overtake all the valid and missing (miss_value) and below-orography (shield_value) values:
          vrpolar( m_all(iobs),&
                   n_all(iobs),&
                   o_all(iobs) ) = radwind_mod_all(iobs)
        END IF

      END DO
!$omp end parallel do

      IF ((ldealiase_vr_obs .AND. lreadmeta_from_netcdf) .OR. lfill_vr_backgroundwind) THEN
        ! save a radial wind field as a reference for dealiasing. This field is not affected
        !  by lmds_vr and elim_shielded below, so that we have a radial wind at all locations
        !  which are not below the orography (lonline=.false.) or at least
        !  are not blocked by the orography (lonline=.true.).
        ALLOCATE(vrpolar_for_dealiasing(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        vrpolar_for_dealiasing = vrpolar
        ! possibility to print it to file, if desired:
        CALL get_fileprefix_ascii_output ('vasim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrpolar_for_dealiasing, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radial velocity, internally used for obs dealiasing', 'm/s', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radwindoutputunit(ista), time_mod, rs_meta(ista), "V_r (for dealiasing) [m/s]", &
                              vrpolar_for_dealiasing, (vrpolar_for_dealiasing > miss_threshold), miss_value, nobs)
        END IF
      END IF

      IF (.NOT. lonline) THEN
        ! Zero out shielded pixels in case of 4/3 earth model. For online ray propagation, this is not necessary
        !   because the ray tracing ends on impact on orography!
        ! IF WE DO THIS AFTER STORING OF vrpolar_for_dealiasing, SHIELDED AREAS WILL GET
        !  A MODEL WIND ESTIMATE FOR vrpolar WHEN lfill_vr_backgroundwind=.TRUE.
        ! IF SHIELDED AREAS SHOULD INSTEAD BE MADE EMPTY, THIS CALL TO elim_shielded() HAS
        !  TO BE MOVED BEFORE THE ABOVE BLOCK WITH "vrpolar_for_dealiasing = vrpolar"!
        CALL elim_shielded(vrpolar, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, miss_value)
      END IF

      IF (lmds_vr) THEN
        ! Take into account the minimum detectable signal (for simplicity assumed equal to that of the reflectivity):
!$omp parallel do private(o,n,m)
        DO o=1, rs_meta(ista)%nel
          DO n=1, rs_meta(ista)%nra
            DO m=1, rs_meta(ista)%naz
              IF (vrpolar(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < MAX(mds_dbz(n),dBZ_crit_radar) ) THEN
                ! Z is below the MDS, so set vrpolar to a MISSING value:
                vrpolar(m,n,o) = miss_value
              END IF
            END DO
          END DO
        END DO
!$omp end parallel do
      END IF

      IF (lfill_vr_backgroundwind) THEN
        ! Fill values with too low reflectivity for a "physically meaningful" radar simulation
        !  with the background radial wind, to provide full data coverage (useful for data assimilation purposes):
!$omp parallel do private(o,n,m)
        DO o=1, rs_meta(ista)%nel
          DO n=1, rs_meta(ista)%nra
            DO m=1, rs_meta(ista)%naz
              IF (vrpolar(m,n,o) < miss_threshold .AND. vrpolar_for_dealiasing(m,n,o) >= miss_threshold) THEN
                vrpolar(m,n,o) = vrpolar_for_dealiasing(m,n,o)
              END IF
            END DO
          END DO
        END DO
!$omp end parallel do
      END IF



      !.. Output the sorted polar data set into a standard file format:
      !   - binary on the SX9, convert to simple ASCII with program "bin2ascii_convrates3d" from Ulrich Blahak
      !   - simple ASCII on all other systems
      CALL get_fileprefix_ascii_output ('vrsim', rs_meta(ista), fileprefix)
      CALL output3d_ascii_radar(itime, &
           rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
           vrpolar, TRIM(ADJUSTL(fileprefix)), &
           'Simul. radial velocity, no aliasing', 'm/s', 'polar', &
           rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

      IF (ldebug_radsim) THEN
        CALL control_output(radwindoutputunit(ista), time_mod, rs_meta(ista), "V_r [m/s]", &
                            vrpolar, (vrpolar > miss_threshold), miss_value, nobs)
      END IF

    END IF ! loutradwind


    IF (lreadmeta_from_netcdf) THEN

      ! Dummy allocation of polarization parameters to avoid calling
      !  the next subroutine with a non-allocated argument. Some compilers don't like this:
      IF (.NOT. ALLOCATED(zdrpolar)) ALLOCATE(zdrpolar(0,0,0))
      IF (.NOT. ALLOCATED(kdppolar)) ALLOCATE(kdppolar(0,0,0))
      IF (.NOT. ALLOCATED(phidppolar)) ALLOCATE(phidppolar(0,0,0))
      IF (.NOT. ALLOCATED(rhvpolar)) ALLOCATE(rhvpolar(0,0,0))

      IF (ldealiase_vr_obs) THEN
        CALL output_radar_country_obs (ista, time_mod, itime, vrpolar, zrpolar, &
                                       zdrpolar, kdppolar, phidppolar, rhvpolar, &
                                       hrpolar, latpolar, lonpolar, &
                                       vr_mod_for_dealiasing=vrpolar_for_dealiasing)
      ELSE
        CALL output_radar_country_obs (ista, time_mod, itime, vrpolar, zrpolar, &
                                       zdrpolar, kdppolar, phidppolar, rhvpolar, &
                                       hrpolar, latpolar, lonpolar)
      END IF

    END IF

    ! .. Output simulated reflectivity and composits. The list of present elevations in obs files has been set when reading the obs data or has been set to the nominal set if no obs have been found:
    IF (loutdbz) THEN

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_bubbles)
#endif
      IF (ldo_bubbles .AND. lreadmeta_from_netcdf) THEN

        IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
          IF (rs_meta(ista)%eleind_for_composite_bub == 98) THEN
            ! .. Construct the composite from the precip scans (not volume scan):
            !     The elevation for lat/lon computations will not be the true elevations
            !     but the nominal one (either 0.4, 0.59, 0.8, or 1.3)
            CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                 comp_meta=comp_meta_bub, &
                                                 eleind_for_composite=1, compdata_tot=comp_dbzsim_bub_tot, &
                                                 ldebug=ldebug_radsim)
          END IF
        ELSE
          IF (rs_meta(ista)%eleind_for_composite_bub == 99) THEN
            ! .. If eleind_for_composite_bub = 99, then make a vertical maximum composite of all elevations:
            DO i=1, rs_meta(ista)%nel
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                   comp_meta=comp_meta_bub, &
                                                   eleind_for_composite=i, compdata_tot=comp_dbzsim_bub_tot, &
                                                   ldebug=ldebug_radsim)
            END DO
          ELSE IF (rs_meta(ista)%eleind_for_composite_bub /= 98) THEN
            ! .. Else, construct the composite from one single elevation of the volume scans, not precip scan
            !      (index eleind_for_composite_bub is relative to the list of actually present elevations):
            CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                 comp_meta=comp_meta_bub, &
                                                 eleind_for_composite=rs_meta(ista)%ind_ele_present( &
                                      MIN( rs_meta(ista)%eleind_for_composite_bub, rs_meta(ista)%nel_present) ), &
                                                 compdata_tot=comp_dbzsim_bub_tot, &
                                                 ldebug=ldebug_radsim)
          END IF
        END IF
      END IF

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_composites)
#endif
      IF (ldo_composite) THEN
        DO i=1, nel_composite
          IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
            IF (rs_meta(ista)%eleindlist_for_composite(i) == 98) THEN
              ! .. Construct the composite from the precip scans (not volume scan):
              !     The elevation for lat/lon computations will not be the true elevations
              !     but the nominal one (either 0.4, 0.59, 0.8, or 1.3)
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                   comp_meta=comp_meta, &
                                                   eleind_for_composite=1, compdata_tot=comp_dbzsim_tot(:,:,i), &
                                                   ldebug=ldebug_radsim)
            END IF
          ELSE
            IF (rs_meta(ista)%eleindlist_for_composite(i) == 99) THEN
              ! .. If eleindlist_for_composite(i) = 99, then make a vertical maximum composite of all elevations:
              DO j=1, rs_meta(ista)%nel
                CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                     comp_meta=comp_meta, &
                                                     eleind_for_composite=j, compdata_tot=comp_dbzsim_tot(:,:,i), &
                                                     ldebug=ldebug_radsim)
              END DO
            ELSE IF (rs_meta(ista)%eleindlist_for_composite(i) /= 98) THEN
              ! .. Else, construct the composite from one single elevation
              !      (index eleindlist_for_composite(i) is relative to the list of actually present elevations):
              CALL composite2D_dbz_maxmethod_ista (zpolar=zrpolar, rsm=rs_meta(ista), &
                                                   comp_meta=comp_meta, &
                                                   eleind_for_composite=rs_meta(ista)%ind_ele_present( &
                                        MIN( rs_meta(ista)%eleindlist_for_composite(i), rs_meta(ista)%nel_present) ), &
                                                   compdata_tot=comp_dbzsim_tot(:,:,i), &
                                                   ldebug=ldebug_radsim)
            END IF
          END IF
        END DO
      END IF
#ifdef __ICON__
    CALL get_runtime_timings (i_fwo_out)
#endif

    END IF ! loutdbz
    
    ! .. Clean up memory:  (UB: This could be made more elegant by if allocated(...) for all fields):
    DEALLOCATE(vrpolar, zrpolar, m_all, n_all, o_all)
    IF (ALLOCATED(zdrpolar))   DEALLOCATE(zdrpolar)
    IF (ALLOCATED(rhvpolar))   DEALLOCATE(rhvpolar)
    IF (ALLOCATED(kdppolar))   DEALLOCATE(kdppolar)
    IF (ALLOCATED(phidppolar)) DEALLOCATE(phidppolar)
    IF (ALLOCATED(ldrpolar))   DEALLOCATE(ldrpolar)
    IF (ALLOCATED(zepolar))    DEALLOCATE(zepolar)
    IF (ALLOCATED(zetpolar))   DEALLOCATE(zetpolar)
    IF (ALLOCATED(zedpolar))   DEALLOCATE(zedpolar)
    IF (ALLOCATED(zedtpolar))  DEALLOCATE(zedtpolar)

#ifdef AUXOUT_OFFLINE
    IF (ALLOCATED(zvpolar))  DEALLOCATE(zvpolar)
    IF (ALLOCATED(zvhpolar)) DEALLOCATE(zvhpolar)
#endif

    IF (ALLOCATED(hrpolar))  DEALLOCATE(hrpolar)
    IF (ALLOCATED(erpolar))  DEALLOCATE(erpolar)
    IF (ALLOCATED(srpolar))  DEALLOCATE(srpolar)
    IF (ALLOCATED(lonpolar)) DEALLOCATE(lonpolar, latpolar)
    IF (ALLOCATED(mds_dbz))  DEALLOCATE(mds_dbz)
    IF (ALLOCATED(vrpolar_for_dealiasing)) DEALLOCATE(vrpolar_for_dealiasing)


    IF (ldebug_radsim) THEN
      WRITE (*,'(a,i0,a,i6.6,a,i0)') TRIM(yzroutine)//': FINISHED output of radar station #', &
           ista, ' (ID: ', rs_meta(ista)%station_id, ') on proc ', my_radar_id
    END IF

  END SUBROUTINE output_my_ista

  !================================================================================

  SUBROUTINE elim_shielded(f, naz, nra, nel, shieldedval)

    !------------------------------------------------------------------------------
    !
    ! Description: Set values in field f to a value of "shieldedval"
    !              outwards in range, beginning from first
    !              shielded value (=988.88) encountered along the ray.
    !
    ! The field f (az_ind, ra_ind, el_ind) is modified in place.
    !
    !------------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(inout), DIMENSION(:,:,:) :: f
    INTEGER, INTENT(in)             :: naz, nra, nel
    REAL(KIND=dp), INTENT(in)       :: shieldedval   ! missing value for shielded rangebins

    INTEGER  :: i, j, k
    LOGICAL  :: flag(naz,nel)

!CDIR COLLAPSE
    flag = .FALSE.

!$omp parallel private(j,k,i)
    DO j=1, nra
!$omp do collapse(2)
      DO k=1, nel
        DO i=1, naz
          IF (.NOT. flag(i,k)) THEN
            IF (f(i,j,k) < shield_up_threshold .AND. f(i,j,k) > shield_low_threshold) THEN
              f(i,j,k) = shieldedval
              flag(i,k) = .TRUE.
            END IF
          ELSE
            f(i,j,k) = shieldedval
          END IF
        END DO
      END DO
!$omp end do
    END DO
!$omp end parallel

  END SUBROUTINE elim_shielded

END MODULE radar_output_station
