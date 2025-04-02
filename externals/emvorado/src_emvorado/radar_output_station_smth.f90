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

MODULE radar_output_station_smth

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
       nradsta_max,                             &
       idom, ndoms_max,                         &
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       i_fwo_bubbles,     & ! Timing flag
       i_fwo_composites,  & ! Timing flag
       i_fwo_out,         & ! Timing flag for output (collecting simulated data on one PE per station, sorting, ASCII-output,
                            !  reading obs data, producing feedback files)
       degrad,             &
       rs_meta, rs_data, dbz_meta, &
       nradsta, &
       missing_obs, &
       comp_meta, comp_meta_bub

  USE radar_composites, ONLY : &
       composite2D_dbz_maxmethod_ista,  &
       comp_dbzsim_tot, &
       comp_dbzsim_bub_tot

  USE radar_interface, ONLY : &
       abort_run,             &
       get_runtime_timings, &
       get_obstime_ind_of_currtime

 !------------------------------------------------------------------------------

   USE radar_utilities, ONLY :   &
                               ind2sub3D, sub2ind3D,       &
                               ind2sub5D, sub2ind5D,       &
                               f4_eff_horzscan,            &
                               smth_el_horzscan,           &
                               smth_az_horzscan,           &
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

  USE radar_obs_meta_list, ONLY : get_elarr_precipscan

  USE radar_model2rays, ONLY : get_azvec, online2geo, rad2geo_const

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

  PUBLIC ::  output_my_ista_smth

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

  SUBROUTINE output_my_ista_smth (time_mod, ista, nsmth, &
                                  radpos_all_smth, vt_mod_all_smth, radwind_mod_all_smth,&
                                  zh_radar_mod_all_smth, ah_radar_mod_all_smth, &
                                  zv_radar_mod_all_smth, &
                                  rrhv_radar_mod_all_smth, irhv_radar_mod_all_smth, &
                                  kdp_radar_mod_all_smth, adp_radar_mod_all_smth, &
                                  zvh_radar_mod_all_smth, &
                                  hl_loc_all_smth, el_loc_all_smth, s_loc_all_smth, &
                                  geom_written)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    REAL(KIND=dp), INTENT(IN) :: time_mod          ! seconds since model start
    INTEGER,       INTENT(IN) :: ista, nsmth
    INTEGER,       POINTER    :: radpos_all_smth(:) ! radar points indices of a radar station
    REAL(KIND=dp), POINTER    :: vt_mod_all_smth(:), radwind_mod_all_smth(:),    &
                                 zh_radar_mod_all_smth(:), ah_radar_mod_all_smth(:), &
                                 zv_radar_mod_all_smth(:), &
                                 rrhv_radar_mod_all_smth(:), irhv_radar_mod_all_smth(:), &
                                 kdp_radar_mod_all_smth(:), adp_radar_mod_all_smth(:), &
                                 zvh_radar_mod_all_smth(:), &
                                 hl_loc_all_smth(:), el_loc_all_smth(:), s_loc_all_smth(:)

    LOGICAL, INTENT(inout)  :: geom_written

    !------------------------------------------------------------------------------
    !
    ! Local variables:

    CHARACTER(LEN=*), PARAMETER:: yzroutine = 'output_my_ista_smth'

    INTEGER                    :: iobs,nobs,i,j,k,m,n,o,iaz,iel,nrp,ismth,irp,iv,ih,iaerr,izs,itime

    REAL(KIND=dp)              :: el,az,el_tmp,az_tmp,ext_coeff

    REAL(KIND=dp), POINTER     :: pchk_sh(:)

    REAL(KIND=dp), ALLOCATABLE :: hrpolar(:,:,:),hl_loc_smth(:,:,:),hl_loc(:),hl_all(:),intgrl_hl(:)
    REAL(KIND=dp), ALLOCATABLE :: lonpolar(:,:,:),latpolar(:,:,:), hdummy(:,:)
    REAL(KIND=dp), ALLOCATABLE :: erpolar(:,:,:), el_loc_smth(:,:,:), el_loc(:), el_all(:), intgrl_el(:)
    REAL(KIND=dp), ALLOCATABLE :: srpolar(:,:,:), s_loc_smth(:,:,:), s_loc(:), s_all(:), intgrl_s(:)
    REAL(KIND=dp), ALLOCATABLE :: vrpolar(:,:,:), radwind_mod_all(:),      &
                                  radwind_mod_smth(:,:,:), intgrl_wind(:), &
                                  vt_mod_smth(:,:,:), intgrl_vt(:),        &
                                  vrpolar_for_dealiasing(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: zrpolar(:,:,:), &
                                  zepolar(:,:,:),zetpolar(:,:,:), &
                                  zdrpolar(:,:,:), &
                                  rhvpolar(:,:,:), &
                                  kdppolar(:,:,:), &
                                  phidppolar(:,:,:), &
                                  ldrpolar(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: zh_radar_mod_all(:), zv_radar_mod_all(:), zvh_radar_mod_all(:), &
                                  kdp_radar_mod_all(:)
    REAL(KIND=dp), ALLOCATABLE :: zh_radar_mod_smth(:,:,:), zv_radar_mod_smth(:,:,:), zvh_radar_mod_smth(:,:,:), &
                                  kdp_radar_mod_smth(:,:,:), rrhv_radar_mod_smth(:,:,:), irhv_radar_mod_smth(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: zh_radar_mod_nmovh(:,:,:,:,:), zv_radar_mod_nmovh(:,:,:,:,:), &
                                  zvh_radar_mod_nmovh(:,:,:,:,:), &
                                  rrhv_radar_mod_nmovh(:,:,:,:,:), irhv_radar_mod_nmovh(:,:,:,:,:), &
                                  ah_radar_mod_nmovh(:,:,:,:,:),aht_radar_mod_nmovh(:,:,:,:,:), &
                                  adp_radar_mod_nmovh(:,:,:,:,:),avt_radar_mod_nmovh(:,:,:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: intgrl_z(:), intgrl_zv(:), intgrl_pol(:), intgrl_pf(:), intgrl_pf_nonblocked(:), &
                                  intgrl_rrhv(:), intgrl_irhv(:), &
                                  pf_value(:,:,:), cos_value(:,:,:), &
                                  mds_dbz(:), pfcos_value(:,:,:)
    REAL(KIND=dp), ALLOCATABLE :: azarr(:), el_sta(:,:), el_precip(:)

    INTEGER      , ALLOCATABLE :: iv_all(:), ih_all(:), m_all(:), n_all(:), o_all(:), &
                                  flag_smth(:,:,:), counter(:), radpos(:), radpos_all(:)

    !!$ zrdpolar et al not yet implemented, just allocated and deallocated as dummies:
    REAL(KIND=dp), ALLOCATABLE :: zdr_radar_mod_all(:), zdr_radar_mod_smth(:,:,:)

    CHARACTER(len=100)               :: fileprefix

    ! For checking the radar points below sfc in the first call:
    LOGICAL, SAVE :: firstcall_offline(nradsta_max,ndoms_max) = .TRUE.


    iaerr = 0
    izs = 0


    IF (ldebug_radsim) THEN
      WRITE (*,'(a,i0,a,i6.6,a,i0)') TRIM(yzroutine)//': output of radar station ', &
           ista, ' (ID: ', rs_meta(ista)%station_id, ' on proc ', my_radar_id
    END IF

    itime = get_obstime_ind_of_currtime ( rs_meta(ista)%obs_times(1:rs_meta(ista)%nobs_times) ) ! missval_int if not found

    ! .. Reset the lsit of present elevations to the list of nominal elevations:
    rs_meta(ista)%nel_present = rs_meta(ista)%nel
    rs_meta(ista)%ind_ele_present(:) = missval_int
    rs_meta(ista)%ind_ele_present(1:rs_meta(ista)%nel_present) = (/ (i, i=1, rs_meta(ista)%nel_present) /)


    ALLOCATE(radpos(nsmth), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(iv_all(nsmth), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(ih_all(nsmth), stat=izs); iaerr = iaerr + ABS(izs)

    CALL init_vari(radpos, 0)
    CALL init_vari(iv_all, 0)
    CALL init_vari(ih_all, 0)

!$omp parallel do private (ismth,irp,m,n,o)
    DO ismth = 1, nsmth
      irp = radpos_all_smth(ismth)
      CALL ind2sub5D( irp, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
                      m, n, o, iv_all(ismth), ih_all(ismth) )
      CALL sub2ind3D( m, n, o, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                      radpos(ismth))
    ENDDO
!$omp end parallel do
    
    nrp  = MAXVAL(radpos)

    ! Create a flag matrix indicating if the radar point
    ! is within the model domain, i.e., > miss_value:
    ALLOCATE(flag_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
             stat=izs); iaerr = iaerr + ABS(izs)
    CALL init_vari(flag_smth, 0)

    IF (.NOT. lonline .AND. firstcall_offline(ista,idom)) THEN

      ! Prepare checking field for shielded pixels depending on configuration:
      ! NOTE: either loutradwind or loutdbz is .true. otherwise output_radar_smth() is not called!
      NULLIFY (pchk_sh)
      IF (loutradwind) THEN
        pchk_sh => radwind_mod_all_smth
      ELSEIF (loutdbz) THEN
        pchk_sh => zh_radar_mod_all_smth
      END IF

      ALLOCATE(zh_radar_mod_nmovh(rs_meta(ista)%naz,rs_meta(ista)%nra, &
               rs_meta(ista)%nel,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
      CALL init_vari(zh_radar_mod_nmovh, miss_value)

      ! n_sh(:,:,:,:) = Range index for beginning of beam shielding:
      IF  (.NOT.ASSOCIATED(rs_data(ista)%n_sh)) THEN
        ALLOCATE(rs_data(ista)%n_sh(rs_meta(ista)%naz, &
                 rs_meta(ista)%nel,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
      END IF
      CALL init_vari(rs_data(ista)%n_sh, rs_meta(ista)%nra+1)

      ! Determine first occurence of beam shielding along each radar ray
      ! and store range index in n_sh(:,:,:,:):
!$omp parallel do private (m,n,o)
      DO ismth = 1, nsmth
        CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                        m, n, o )
        zh_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) = &
             pchk_sh(ismth)
      END DO
!$omp end parallel do

!$omp parallel do private (ih,iv,i,j,k)
      DO ih = 1, rs_meta(ista)%ngpsm_h
        DO iv = 1, rs_meta(ista)%ngpsm_v
          DO i = 1, rs_meta(ista)%nel
            DO j = 1, rs_meta(ista)%nra
              DO k = 1, rs_meta(ista)%naz

                IF ( zh_radar_mod_nmovh(k,j,i,iv,ih) > shield_low_threshold .AND. &
                     zh_radar_mod_nmovh(k,j,i,iv,ih) < shield_up_threshold .AND. &
                     j < rs_data(ista)%n_sh(k,i,iv,ih) ) THEN

                  rs_data(ista)%n_sh(k,i,iv,ih) = j

                END IF

              END DO
            END DO
          END DO
        END DO
      END DO
!$omp end parallel do

      NULLIFY (pchk_sh)
      DEALLOCATE (zh_radar_mod_nmovh)

      firstcall_offline(ista,idom) = .FALSE.

    END IF

    IF (lonline) THEN

      ! Shielded pixels are not contained in the data vectors because
      !  the ray tracing ends on impact on orography!
!$omp parallel do
      DO ismth = 1, nsmth
        flag_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
             flag_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) + 1
      ENDDO
!$omp end parallel do

    ELSE

      ! Flag out shielded pixels in case of 4/3 earth model. For online ray propagation, this is not necessary
      !  because the ray tracing ends on impact on orography!

!$omp parallel do private (ismth,m,n,o)
      DO ismth = 1, nsmth
        CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                        m, n, o )
        IF (n < rs_data(ista)%n_sh(m,o,iv_all(ismth),ih_all(ismth))) THEN
!!! first index is sorted according to azimut, range, and elevation!
          flag_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
               flag_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) + 1
        END IF
      ENDDO
!$omp end parallel do
    END IF

    IF (ANY(flag_smth > 1)) THEN
      WRITE (*,*) 'WARNING: More than one data value for some smth-points!'
    END IF

    ALLOCATE(radpos_all(nrp), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(counter(nrp),    stat=izs); iaerr = iaerr + ABS(izs)

    CALL init_vari(radpos_all, 0)
    CALL init_vari(counter, 0)

    DO ih = 1, rs_meta(ista)%ngpsm_h
      DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do private (irp)
        DO irp = 1, nrp
          IF(flag_smth(irp,iv,ih) /= 0) THEN
            counter(irp) = counter(irp) + 1
          ENDIF
        ENDDO
!$omp end parallel do
      ENDDO
    ENDDO

    nobs = 0
!!! This loop cannot easily be OMP parallelized ...
    DO irp = 1, nrp
      IF (counter(irp) > 0) THEN
        nobs = nobs + 1
!!! nobs are sorted according to azimut, range, and elevation!
        radpos_all(nobs) = irp
      ENDIF
    ENDDO

    ALLOCATE(m_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(n_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(o_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)

    CALL init_vari(m_all, -1)
    CALL init_vari(n_all, -1)
    CALL init_vari(o_all, -1)

!$omp parallel do private (iobs,irp)
    DO iobs = 1, nobs
      irp = radpos_all(iobs)
      CALL ind2sub3D( irp, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                       m_all(iobs), n_all(iobs), o_all(iobs) )
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


!!$==============================================================================================
!!$==============================================================================================
!!$
!!$ !!!!!!!!!!!!!!!!!!!! Compute coordinates and output heights and local elevations !!!!!!!!!!!!


    IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written)) .OR. lreadmeta_from_netcdf ) THEN

      ! .. Compute bin-averaged height of auxiliary rays and store on hrpolar:

      ALLOCATE(hl_loc_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), stat=izs); iaerr = iaerr + ABS(izs)
      ALLOCATE(hl_loc(nrp), stat=izs); iaerr = iaerr + ABS(izs)
      ALLOCATE(intgrl_hl(nrp), stat=izs); iaerr = iaerr + ABS(izs)

      CALL init_vari(hl_loc_smth, miss_value)
      CALL init_vari(hl_loc, 0.0_dp)
      CALL init_vari(intgrl_hl, 0.0_dp)

!$omp parallel do
      DO ismth = 1, nsmth
        hl_loc_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) =  hl_loc_all_smth(ismth)
      ENDDO
!$omp end parallel do

      DO ih = 1, rs_meta(ista)%ngpsm_h
        DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do private (irp)
          DO irp = 1, nrp
            IF (flag_smth(irp,iv,ih) > 0 .AND. hl_loc_smth(irp,iv,ih) > miss_threshold) THEN
              hl_loc(irp) = hl_loc(irp) +  &
                   rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih) * &
                   hl_loc_smth(irp,iv,ih)
              intgrl_hl(irp) = intgrl_hl(irp) + &
                   rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih)
            END IF
          ENDDO
!$omp end parallel do
        ENDDO
      ENDDO

      DEALLOCATE(hl_loc_smth)

      ALLOCATE(hl_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)

!$omp parallel do private (iobs,irp)
      DO iobs = 1, nobs
        irp = radpos_all(iobs)
        IF (intgrl_hl(irp) >= 1e-20_dp) THEN
          hl_all(iobs) = hl_loc(irp) / intgrl_hl(irp)
        ELSE
          hl_all(iobs) = miss_value
        END IF
      ENDDO
!$omp end parallel do

      ALLOCATE(hrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               stat=izs); iaerr = iaerr + ABS(izs)
      CALL init_vari(hrpolar, miss_value)

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!$omp parallel do
      DO iobs = 1, nobs
        hrpolar( m_all(iobs),&
                 n_all(iobs),&
                 o_all(iobs) ) = hl_all(iobs)
      END DO
!$omp end parallel do

      !.. Clean up some memory:
      IF(ALLOCATED(hl_all))           DEALLOCATE(hl_all)
      IF(ALLOCATED(hl_loc))           DEALLOCATE(hl_loc)
      IF(ALLOCATED(intgrl_hl))        DEALLOCATE(intgrl_hl)

    END IF  ! calculate hl_all, output hrpolar



    IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written)) .OR. (loutradwind .AND. lfall) &
                                                            .OR. lreadmeta_from_netcdf) THEN

      ! .. Compute bin-averaged local elevation of auxiliary rays and store on erpolar:

      ALLOCATE(el_loc_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), stat=izs); iaerr = iaerr + ABS(izs)
      ALLOCATE(el_loc(nrp), stat=izs); iaerr = iaerr + ABS(izs)
      ALLOCATE(intgrl_el(nrp), stat=izs); iaerr = iaerr + ABS(izs)

      CALL init_vari(el_loc_smth, miss_value)
      CALL init_vari(el_loc, 0.0_dp)
      CALL init_vari(intgrl_el, 0.0_dp)

!$omp parallel do
      DO ismth = 1, nsmth
        el_loc_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) =  el_loc_all_smth(ismth)
      ENDDO
!$omp end parallel do

      DO ih = 1, rs_meta(ista)%ngpsm_h
        DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do private (irp)
          DO irp = 1, nrp
            IF (flag_smth(irp,iv,ih) > 0 .AND. el_loc_smth(irp,iv,ih) > miss_threshold) THEN
              el_loc(irp) = el_loc(irp) + &
                   rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih) * &
                   el_loc_smth(irp,iv,ih)
              intgrl_el(irp) = intgrl_el(irp) + &
                   rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih)
            ENDIF
          ENDDO
!$omp end parallel do
        ENDDO
      ENDDO

      ALLOCATE(el_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)

!$omp parallel do private (iobs,irp)
      DO iobs = 1, nobs
        irp = radpos_all(iobs)
        IF (counter(irp) > 0 .AND. intgrl_el(irp) >= 1e-20_dp) THEN
          el_all(iobs) = el_loc(irp) / intgrl_el(irp)
        ELSE
          el_all(iobs) = miss_value
        END IF
      ENDDO
!$omp end parallel do


      ALLOCATE(erpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               stat=izs); iaerr = iaerr + ABS(izs)
      CALL init_vari(erpolar, miss_value)

      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!$omp parallel do private (iobs)
      DO iobs = 1, nobs
        erpolar( m_all(iobs),&
                 n_all(iobs),&
                 o_all(iobs) ) = el_all(iobs)
      END DO
!$omp end parallel do

      !.. Clean up some memory:
      IF(ALLOCATED(el_all))           DEALLOCATE(el_all)
      IF(ALLOCATED(el_loc))           DEALLOCATE(el_loc)
      IF(ALLOCATED(intgrl_el))        DEALLOCATE(intgrl_el)

    END IF  ! calculate el_all


    IF ( (lout_geom .AND. (lonline .OR. .NOT.geom_written)) .OR. lreadmeta_from_netcdf ) THEN

      ! .. Compute lon/lat of radar points in a way that there are no missing values
      !     in the resulting lat/lon polar arrays (i.e., all coordinates are defined,
      !     even if the data are missing values):

      ALLOCATE(lonpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               latpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
      CALL init_vari(lonpolar, miss_value)
      CALL init_vari(latpolar, miss_value)

      IF (lonline) THEN

        ! .. Compute bin-averaged arc distance of auxiliary rays and store on erpolar:

        ALLOCATE(s_loc_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h))
        ALLOCATE(s_loc(nrp))
        ALLOCATE(intgrl_s(nrp))

        CALL init_vari(s_loc_smth, miss_value)
        CALL init_vari(s_loc, 0.0_dp)
        CALL init_vari(intgrl_s, 0.0_dp)

!$omp parallel do
        DO ismth = 1, nsmth
          s_loc_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) =  s_loc_all_smth(ismth)
        ENDDO
!$omp end parallel do

        DO ih = 1, rs_meta(ista)%ngpsm_h
          DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do
            DO irp = 1, nrp
              IF (flag_smth(irp,iv,ih) > 0 .AND. s_loc_smth(irp,iv,ih) > miss_threshold) THEN
                s_loc(irp) = s_loc(irp) + &
                     rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih) * &
                     s_loc_smth(irp,iv,ih)
                intgrl_s(irp) = intgrl_s(irp) + &
                     rs_meta(ista)%weigsm_v(iv) * rs_meta(ista)%weigsm_h(ih)
              ENDIF
            ENDDO
!$omp end parallel do
          ENDDO
        ENDDO

        ALLOCATE(s_all(nobs))

!$omp parallel do private (iobs,irp)
        DO iobs = 1, nobs
          irp = radpos_all(iobs)
          IF (counter(irp) > 0 .AND. intgrl_s(irp) >= 1e-20_dp) THEN
            s_all(iobs) = s_loc(irp) / intgrl_s(irp)
          ELSE
            s_all(iobs) = miss_value
          END IF
        ENDDO
!$omp end parallel do


        ALLOCATE(srpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(srpolar, miss_value)
!$omp parallel do private (iobs)
        DO iobs = 1, nobs
          srpolar( m_all(iobs),&
                   n_all(iobs),&
                   o_all(iobs) ) = s_all(iobs)
        END DO
!$omp end parallel do

        !.. Clean up some memory:
        IF(ALLOCATED(s_all))           DEALLOCATE(s_all)
        IF(ALLOCATED(s_loc))           DEALLOCATE(s_loc)
        IF(ALLOCATED(intgrl_s))        DEALLOCATE(intgrl_s)
        IF(ALLOCATED(s_loc_smth))      DEALLOCATE(s_loc_smth)

        ! .. lat/lon are computed from rs_data(ista)%s_loc, and for blocked/missing ray
        !    parts, the 4/3 earth radius model is used for extrapolation (srpolar, hrpolar are filled):
        CALL online2geo(rs_meta(ista), srpolar, hrpolar, erpolar, &
                        latpolar, lonpolar)

      ELSE

        IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
          ALLOCATE(hdummy(rs_meta(ista)%naz,rs_meta(ista)%nra))
          CALL rad2geo_const (rs_meta(ista), latpolar, lonpolar, hdummy)
          ! .. Fill missing values in hrpolar (above model top) with values from hdummy
          !     (should be continuous along rays, because for both the same 4/3 earth radius model has been applied):
!$omp parallel do private(i,j,k) collapse (3)
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
!$omp parallel do private(i,j,k) collapse (3)
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

    END IF

    IF ( lout_geom .AND. (lonline .OR. .NOT.geom_written) ) THEN

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
        CALL control_output(radgeomoutputunit(ista), time_mod, rs_meta(ista), "el_polar [m]", &
             erpolar, (erpolar > miss_threshold), miss_value, nobs)
      END IF

      ! Set flag that geometry for this radar has been written to a file:
      geom_written = .TRUE.

    END IF



!!$ !!!!!!!!!!!!!!!!! End of computations and output of heights and local elevations !!!!!!!!!!!!
!!$
!!$==============================================================================================
!!$==============================================================================================

    ! Allocate fields for sum of integration weights, beam function and functional determinant
    ! of the coordinate trafo kartesian -> spherical:
    ALLOCATE(intgrl_pf(nrp), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(intgrl_pf_nonblocked(nrp), stat=izs); iaerr = iaerr + ABS(izs)
    ALLOCATE(pfcos_value(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), stat=izs); iaerr = iaerr + ABS(izs)

    CALL init_vari(intgrl_pf, 0.0_dp)
    CALL init_vari(intgrl_pf_nonblocked, 0.0_dp)
    CALL init_vari(pfcos_value, 0.0_dp)

    ! allocate and fill up azimuth array
    CALL get_azvec (rs_meta(ista), azarr)

    ! fill up antenna elevation array:
    ALLOCATE(el_sta(rs_meta(ista)%naz,rs_meta(ista)%nel))
    IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
      CALL get_elarr_precipscan (rs_meta(ista), el_precip)
      el_sta(:,1) = el_precip
    ELSE
      DO iel=1, rs_meta(ista)%nel
        el_sta(:,iel) = rs_meta(ista)%el_arr(iel)
      END DO
    END IF

!$omp parallel do private (ismth,irp,az,el,el_tmp,az_tmp,m,n,o)
    DO ismth = 1, nsmth

      irp = radpos(ismth)

      IF (counter(irp) /= 0) THEN

        CALL ind2sub3D( irp, rs_meta(ista)%naz, rs_meta(ista)%nra, &
                        m, n, o )

        az  = azarr( m )

        el  = rs_meta(ista)%el_arr( o )

        el_tmp =  el_sta(m,o) + smth_el_horzscan( rs_meta(ista)%Theta3, &
                                                  rs_meta(ista)%smth_interv_fact, &
                                                  rs_meta(ista)%xabscsm_v(iv_all(ismth)) )

        az_tmp =  az + smth_az_horzscan( &
                                         rs_meta(ista)%alpha3_eff_0, &
                                         rs_meta(ista)%dalpha, &
                                         rs_meta(ista)%Phi3, &
                                         rs_meta(ista)%smth_interv_fact, &
                                         rs_meta(ista)%xabscsm_h(ih_all(ismth)), &
                                         el &
                                       )
        ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.

        pfcos_value(irp,iv_all(ismth),ih_all(ismth)) = &
             f4_eff_horzscan( rs_meta(ista)%alpha3_eff_0, &
                              rs_meta(ista)%dalpha, &
                              rs_meta(ista)%Phi3, &
                              rs_meta(ista)%Theta3, &
                              az,el,az_tmp,el_tmp ) * COS(el_tmp*degrad)

      ENDIF

    ENDDO
!$omp end parallel do

    DEALLOCATE (azarr, el_sta)
    IF (ALLOCATED(el_precip)) DEALLOCATE (el_precip)

    DO ih = 1, rs_meta(ista)%ngpsm_h
      DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do private (irp)
        DO irp = 1, nrp

          ! Integral over beam weighting function including missing points.
          ! This makes sense for reflectivity calculations, where this integral
          ! is the denominator. In this way, missing values (e.g., shielded values)
          ! enter the beam-function averaged value as zeros with full weight.

          ! Sum up all values when at least one of them is not a missing value:
          IF (counter(irp) > 0) THEN

            intgrl_pf(irp) = intgrl_pf(irp) + &
                 rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih) * &
                 pfcos_value(irp,iv,ih)

          END IF

          ! Sum up all values which are not blocked:
          IF (flag_smth(irp,iv,ih) > 0) THEN

            intgrl_pf_nonblocked(irp) = intgrl_pf(irp) + &
                 rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih) * &
                 pfcos_value(irp,iv,ih)

          END IF

        ENDDO
!$omp end parallel do
      ENDDO
    ENDDO

    IF (lmds_z .OR. lmds_vr) THEN
      ! Minimum detectable signal as function of range in dBZ:
      ALLOCATE( mds_dbz(rs_meta(ista)%nra) )
!$omp parallel do private (iobs)
      DO iobs = 1, rs_meta(ista)%nra
        mds_dbz(iobs) = rs_meta(ista)%mds_Z0 + 20.0_dp * LOG10(iobs*rs_meta(ista)%ra_inc/rs_meta(ista)%mds_r0)
      END DO
!$omp end parallel do
    END IF

    IF (loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))) THEN

      ! zh_radar_mod_smth should hold the linear value of reflectivity with index order
      ! (radpos, phi_points, theta_points):
      ALLOCATE(zh_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), stat=izs); iaerr = iaerr + ABS(izs)
      CALL init_vari(zh_radar_mod_smth, miss_value)


      IF (loutdbz .AND.(loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

        ALLOCATE(zv_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
        CALL init_vari(zv_radar_mod_smth, miss_value)

        ALLOCATE(rrhv_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
        CALL init_vari(rrhv_radar_mod_smth, miss_value)
        ALLOCATE(irhv_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
        CALL init_vari(irhv_radar_mod_smth, miss_value)
        
        ALLOCATE(kdp_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
        CALL init_vari(kdp_radar_mod_smth, miss_value)
        
        IF (loutpolall) THEN
          ALLOCATE(zvh_radar_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(zvh_radar_mod_smth, miss_value)
        END IF

      END IF

      ! Consider attenuation effects
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN

        ! Non-polarimetric parameters
        ALLOCATE(zh_radar_mod_nmovh(&
             rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
             rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
             stat=izs); iaerr = iaerr + ABS(izs)
        ALLOCATE(ah_radar_mod_nmovh(&
             rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
             rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
             stat=izs); iaerr = iaerr + ABS(izs)
        ! UB: Initialize with standard missing value (= point outside model domain):
        CALL init_vari(zh_radar_mod_nmovh, miss_value)
        CALL init_vari(ah_radar_mod_nmovh, miss_value)

        ! JM201005:
        ! aht_radar_mod_nmovh is allocated as (naz,nra) only, ie without nel and
        ! without ngpsm_v/h dimensions. This is as aht_radar_mod_nmovh is only
        ! used as a temporary container.
        ! NOTE: Neither zepolar nor zetpolar are output to when smoothed (because
        ! their correct calc is a bit more tricky and their usefulness is
        ! somewhat limited.
        ALLOCATE(aht_radar_mod_nmovh(&
             rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
             rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
             stat=izs); iaerr = iaerr + ABS(izs)

!$omp parallel do private(ismth,m,n,o)
        DO ismth = 1, nsmth
          CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                          m, n, o )
          zh_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) = &
               zh_radar_mod_all_smth(ismth)
          ah_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
               ah_radar_mod_all_smth(ismth)
        END DO
!$omp end parallel do

        aht_radar_mod_nmovh(:,1,:,:,:) = 0.0_dp
!$omp parallel do private(ih,iv,i,j,k,ext_coeff) collapse(2)
        DO ih = 1, rs_meta(ista)%ngpsm_h
          DO iv = 1, rs_meta(ista)%ngpsm_v
            DO i = 1, rs_meta(ista)%nel
              DO j = 1, rs_meta(ista)%nra
                DO k = 1, rs_meta(ista)%naz
                  ext_coeff = MAX(ah_radar_mod_nmovh(k,j,i,iv,ih), 0.0_dp)  ! one-way [1/m]
                  aht_radar_mod_nmovh(k,j,i,iv,ih) = aht_radar_mod_nmovh(k,MAX(j-1,1),i,iv,ih) + &
                       ext_coeff*rs_meta(ista)%ra_inc
                  ! Note: z_radar_mod_nmovh is in linear units, so correct for
                  ! extinction in linear space:
                  IF (zh_radar_mod_nmovh(k,j,i,iv,ih) >= 0.0_dp) THEN
                    ! preserving missing values < 0 !
                    zh_radar_mod_nmovh(k,j,i,iv,ih) = zh_radar_mod_nmovh(k,j,i,iv,ih) * &
                         EXP(-2.0_dp * aht_radar_mod_nmovh(k,j,i,iv,ih))
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO
!$omp end parallel do

!$omp parallel do private (ismth,m,n,o)
        DO ismth = 1, nsmth
          CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                          m, n, o )
          ! still linear reflectivity value:
          zh_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
               zh_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth))
        END DO
!$omp end parallel do

        DEALLOCATE(zh_radar_mod_nmovh)

        !========================================================================
        !
        ! Polarimetric parameters (only needed if loutdbz=.true.)
        ! JM201005:
        ! Here, the extinction effects are calculated and imprinted onto the
        ! extinction affected parameters. Hence, we only treat those here (ie
        ! ZDR, LDR).
        ! Beam integration is done further below.
        ! Range integration as required for PhiDP is not done at all yet.
        !

        IF (loutdbz .AND. (loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN

          ALLOCATE(zv_radar_mod_nmovh(&
               rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
               rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          ALLOCATE(rrhv_radar_mod_nmovh(&
               rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
               rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          ALLOCATE(irhv_radar_mod_nmovh(&
               rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
               rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          ALLOCATE(adp_radar_mod_nmovh(&
               rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
               rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          ALLOCATE(avt_radar_mod_nmovh(&
               rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
               rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(zv_radar_mod_nmovh, miss_value)
          CALL init_vari(rrhv_radar_mod_nmovh, miss_value)
          CALL init_vari(irhv_radar_mod_nmovh, miss_value)
          CALL init_vari(adp_radar_mod_nmovh, miss_value)
          CALL init_vari(avt_radar_mod_nmovh, 0.0_dp)

          IF (loutpolall) THEN
            ALLOCATE(zvh_radar_mod_nmovh(&
                     rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel,&
                     rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                     stat=izs); iaerr = iaerr + ABS(izs)
            CALL init_vari(zvh_radar_mod_nmovh, miss_value)
          END IF

!$omp parallel do private(ismth,m,n,o)
          DO ismth = 1, nsmth
            CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                            m, n, o )
            adp_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
                 adp_radar_mod_all_smth(ismth)
            zv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
                 zv_radar_mod_all_smth(ismth)
            rrhv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
                 rrhv_radar_mod_all_smth(ismth)
            irhv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
                 irhv_radar_mod_all_smth(ismth)
            IF (loutpolall) THEN
              zvh_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth)) =  &
                   zvh_radar_mod_all_smth(ismth)
            END IF
          END DO
!$omp end parallel do

          aht_radar_mod_nmovh(:,1,:,:,:) = 0.0_dp
          avt_radar_mod_nmovh(:,1,:,:,:) = 0.0_dp
!$omp parallel
!$omp do private(ih,iv,i,j,k,ext_coeff) collapse(2)
          DO ih = 1, rs_meta(ista)%ngpsm_h
            DO iv = 1, rs_meta(ista)%ngpsm_v
              DO i = 1, rs_meta(ista)%nel
                DO j = 1, rs_meta(ista)%nra
                  DO k = 1, rs_meta(ista)%naz
                    ! adp can be < 0, hence MAX(val,0) is not applicable here
                    IF (adp_radar_mod_nmovh(k,j,i,iv,ih) >= miss_threshold) THEN
                      ! ext_coeff V in linear units [1/m], therefore no cext factor:
                      ext_coeff = MAX(ah_radar_mod_nmovh(k,j,i,iv,ih),0.0_dp) - adp_radar_mod_nmovh(k,j,i,iv,ih)  ! one-way [1/m]
                    ELSE
                      ext_coeff = MAX(ah_radar_mod_nmovh(k,j,i,iv,ih),0.0_dp)
                    END IF
                    avt_radar_mod_nmovh(k,j,i,iv,ih) = avt_radar_mod_nmovh(k,MAX(j-1,1),i,iv,ih) + &
                         ext_coeff*rs_meta(ista)%ra_inc
                    IF (zv_radar_mod_nmovh(k,j,i,iv,ih) >= miss_threshold) THEN
                      ! zv_radar_mod_nmovh is in linear units, so correct for extinction
                      !  in linear units:
                      zv_radar_mod_nmovh(k,j,i,iv,ih) = zv_radar_mod_nmovh(k,j,i,iv,ih) * &
                           EXP(-2.0_dp * avt_radar_mod_nmovh(k,j,i,iv,ih))
                    END IF
                    IF (rrhv_radar_mod_nmovh(k,j,i,iv,ih) >= miss_threshold) THEN
                      rrhv_radar_mod_nmovh(k,j,i,iv,ih) = rrhv_radar_mod_nmovh(k,j,i,iv,ih) * &
                           EXP(-aht_radar_mod_nmovh(k,j,i,iv,ih)-avt_radar_mod_nmovh(k,j,i,iv,ih))
                      irhv_radar_mod_nmovh(k,j,i,iv,ih) = irhv_radar_mod_nmovh(k,j,i,iv,ih) * &
                           EXP(-aht_radar_mod_nmovh(k,j,i,iv,ih)-avt_radar_mod_nmovh(k,j,i,iv,ih))
                    END IF
                    IF (loutpolall .AND. zvh_radar_mod_nmovh(k,j,i,iv,ih) >= miss_threshold) THEN
                      zvh_radar_mod_nmovh(k,j,i,iv,ih) = zvh_radar_mod_nmovh(k,j,i,iv,ih) * &
                           EXP(-aht_radar_mod_nmovh(k,j,i,iv,ih)-avt_radar_mod_nmovh(k,j,i,iv,ih))
                    END IF
                  END DO
                END DO
              END DO
            END DO
          END DO
!$omp end do

!$omp do private (ismth,m,n,o)
          DO ismth = 1, nsmth
            CALL ind2sub3D( radpos(ismth), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                            m, n, o )
            zv_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
                 zv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth))
            rrhv_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
                 rrhv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth))
            irhv_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
                 irhv_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth))
            IF (loutpolall) THEN
              zvh_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
                   zvh_radar_mod_nmovh(m,n,o,iv_all(ismth),ih_all(ismth))
            END IF
          END DO
!$omp end do
!$omp end parallel
          
        END IF

        IF (ALLOCATED(adp_radar_mod_nmovh)) DEALLOCATE(adp_radar_mod_nmovh)
        IF (ALLOCATED(avt_radar_mod_nmovh)) DEALLOCATE(avt_radar_mod_nmovh)
        IF (ALLOCATED(zvh_radar_mod_nmovh)) DEALLOCATE(zvh_radar_mod_nmovh)
        IF (ALLOCATED(zv_radar_mod_nmovh))  DEALLOCATE(zv_radar_mod_nmovh)
        IF (ALLOCATED(rrhv_radar_mod_nmovh))  DEALLOCATE(rrhv_radar_mod_nmovh)
        IF (ALLOCATED(irhv_radar_mod_nmovh))  DEALLOCATE(irhv_radar_mod_nmovh)

      ELSE  ! .not. lextdbz

        CALL sort_ismth_to_radpos_iv_ih (zh_radar_mod_all_smth, zh_radar_mod_smth)

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          CALL sort_ismth_to_radpos_iv_ih (rrhv_radar_mod_all_smth, rrhv_radar_mod_smth)
          CALL sort_ismth_to_radpos_iv_ih (irhv_radar_mod_all_smth, irhv_radar_mod_smth)
          CALL sort_ismth_to_radpos_iv_ih (zv_radar_mod_all_smth, zv_radar_mod_smth)
          IF (loutpolall) THEN
            CALL sort_ismth_to_radpos_iv_ih (zvh_radar_mod_all_smth, zvh_radar_mod_smth)
          END IF
        END IF

      END IF  ! lextdbz

      IF (ALLOCATED(ah_radar_mod_nmovh))        DEALLOCATE(ah_radar_mod_nmovh)
      IF (ALLOCATED(aht_radar_mod_nmovh))       DEALLOCATE(aht_radar_mod_nmovh)
      
      ! .. 2D-integral over (ze_radar/l_ext^2)*f4*cos(theta) dtheta dphi:
      !   ( NOTE: the integral weights for Gauss-Legendre are already
      !     pre-multiplied with the total norm. integration range [-1;1],
      !     but the transformation factor (b-a)/2 to the "true" ranges,
      !     which are, e.g.,
      !     rs_meta%smth_interv_fact*[-theta3/2;+theta3/2]
      !     i.e.,
      !     rs_meta%smth_interv_fact * theta3, are missing here,
      !     since the integrals only appear later as quotients,
      !     so that a constant factor cancels out. )
      !
      ALLOCATE(intgrl_z(nrp), stat=izs); iaerr = iaerr + ABS(izs)

      ! Sum up rs_meta(ista)%weigsm_v*rs_meta(ista)%weigsm_h*zh_radar_mod_smth*pfcos_value:
      CALL weightedsum_beamfunction (zh_radar_mod_smth, 0.0_dp, intgrl_z)

      ALLOCATE(zrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
               stat=izs); iaerr = iaerr + ABS(izs)
      CALL init_vari(zrpolar, miss_value)

      IF (loutdbz .OR. (loutradwind .AND. lmds_vr)) THEN

        ALLOCATE(zh_radar_mod_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)

        ! Compute zh_radar_mod_all = intgrl_z / intgrl_pf for all points where counter(irp) > 0, set miss_value otherwise:
        CALL divide_intpol_by_intpf (intgrl_z, intgrl_pf, zh_radar_mod_all)

        !.. Sort the collected radar points into the 3D field according
        !   to sorted azimut, range and elevation:
        DO iobs = 1, nobs
          IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
            zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 10.0_dp * LOG10(zh_radar_mod_all(iobs))
          ELSEIF (zh_radar_mod_all(iobs) >= miss_threshold) THEN
            zrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zero_value
          END IF
        END DO

        ! Set correct 0-value (zero_value) and missing value (miss_value) for zrpolar:
        CALL set_missing_and_correct0(zrpolar)

        IF (lmds_z) THEN
          CALL mask_mindetectsignal ( zrpolar, zero_value )
        END IF

      END IF ! loutdbz .or. (loutradwind .and. lmds_vr)


      IF (loutdbz) THEN

        !.. Output the sorted polar data set into a standard file format:
        !   - binary on the SX9, convert to simple ASCII with program "bin2ascii_convrates3d" from Ulrich Blahak

        IF (ldebug_radsim) THEN
          CALL control_output(radrefloutputunit(ista), time_mod, rs_meta(ista), "Z [dBZ]", &
               zrpolar, (zrpolar > miss_threshold), miss_value, nobs)
        END IF

        !   - simple ASCII on all other systems
        CALL get_fileprefix_ascii_output ('zrsim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             zrpolar, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radar reflectivity', 'dBZ', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)


        IF ((loutpolall .OR. loutpolstd) .AND.  dbz_meta(ista)%itype_refl > 4) THEN

          !=========================================================================
          ! .. Polarization parameters:
          !=========================================================================


          ! .. ZV and ZDR:
          !=========================================================================
          
          ALLOCATE(intgrl_zv(nrp), stat=izs); iaerr = iaerr + ABS(izs)
          
          ! Sum up rs_meta(ista)%weigsm_v*rs_meta(ista)%weigsm_h*zv_radar_mod_smth*pfcos_value:
          CALL weightedsum_beamfunction (zv_radar_mod_smth, 0.0_dp, intgrl_zv)

          ALLOCATE(zv_radar_mod_all(nobs),  stat=izs); iaerr = iaerr + ABS(izs)
          ! Compute zvh_radar_mod_all = intgrl_zv / intgrl_pf for all points where counter(irp) > 0, set miss_value otherwise:
          CALL  divide_intpol_by_intpf (intgrl_zv, intgrl_pf, zv_radar_mod_all)

          ALLOCATE(zdrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(zdrpolar, miss_value)

          !.. Sort the collected radar points into the 3D field according
          !   to sorted azimut, range and elevation and scale zdr to log for output:
          DO iobs = 1, nobs
            IF (zh_radar_mod_all(iobs) >= Z_crit_radar .AND. zv_radar_mod_all(iobs) >= Z_crit_radar) THEN
              zdrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = 10.0_dp * LOG10(zh_radar_mod_all(iobs)/zv_radar_mod_all(iobs))
            END IF
          END DO
          
          ! Impose minimum detectable signal:
          IF (lmds_z) THEN
            CALL mask_mindetectsignal ( zdrpolar, miss_value )
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

          ! .. KPD, PHIDP:
          !=========================================================================
          
          ALLOCATE(intgrl_pol(nrp), stat=izs); iaerr = iaerr + ABS(izs)

          CALL sort_ismth_to_radpos_iv_ih (kdp_radar_mod_all_smth, kdp_radar_mod_smth)

          ! Sum up rs_meta(ista)%weigsm_v*rs_meta(ista)%weigsm_h*kdp_radar_mod_smth*pfcos_value:
          CALL weightedsum_beamfunction (kdp_radar_mod_smth, miss_threshold, intgrl_pol)

          ALLOCATE(kdp_radar_mod_all(nobs),  stat=izs); iaerr = iaerr + ABS(izs)
          ! Compute kdp_radar_mod_all = intgrl_pol / intgrl_pf_nonblocked for all points where counter(irp) > 0,
          !  set miss_value otherwise. The normalization with the "non-blocked" part of the beam function
          !  is done because we are dealing with a phase shift, not an energy:
          CALL divide_intpol_by_intpf (intgrl_pol, intgrl_pf_nonblocked, kdp_radar_mod_all)

          ALLOCATE(kdppolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(kdppolar, miss_value)
          CALL sort_into_polarcoords (kdp_radar_mod_all, miss_threshold, kdppolar)  ! [deg/m]
          DEALLOCATE (kdp_radar_mod_all)
          
          ALLOCATE(phidppolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(phidppolar, 0.0_dp)
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
!$omp parallel
          WHERE (kdppolar < miss_threshold) phidppolar = miss_value 
!$omp end parallel

          ! Impose minimum detectable signal:
          IF (lmds_z) THEN
            CALL mask_mindetectsignal ( kdppolar, miss_value )
          END IF

        ! FIXME:
        ! Should a phase wrapping (to [0..360] or [-180..180]deg) be applied?
        ! In order to compare to obs that seems reasonable.
        ! But so far, we don't.
        !
        ! YES it should: the obs of DWD are from -180 to 180

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

          ! .. RHOHV:
          !=========================================================================

          ALLOCATE(intgrl_rrhv(nrp), stat=izs); iaerr = iaerr + ABS(izs)
          ALLOCATE(intgrl_irhv(nrp), stat=izs); iaerr = iaerr + ABS(izs)

          ! Attenuated rrhv and irhv integrated over beam weighting function (not averaged!):
          CALL weightedsum_beamfunction (rrhv_radar_mod_smth, miss_thresh_rhv, intgrl_rrhv)
          CALL weightedsum_beamfunction (irhv_radar_mod_smth, miss_thresh_rhv, intgrl_irhv)

          ALLOCATE(rhvpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), stat=izs); iaerr = iaerr + ABS(izs)
          CALL init_vari(rhvpolar, miss_value)
          DO iobs = 1, nobs
            irp = radpos_all(iobs)
            ! Simplified criterion: if beam averaged zh is > Z_crit_radar, intgrl_z*intgrl_zv should also be > 0.
            IF (counter(irp) > 0 .AND. zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
              rhvpolar(m_all(iobs),n_all(iobs),o_all(iobs)) = &
                   SQRT( ( intgrl_rrhv(irp)**2 + intgrl_irhv(irp)**2 ) / (intgrl_z(irp)*intgrl_zv(irp)) )
            ELSE
              rhvpolar(m_all(iobs),n_all(iobs),o_all(iobs)) = miss_value
            END IF
          END DO

          DEALLOCATE(intgrl_zv)
          DEALLOCATE(intgrl_rrhv, intgrl_irhv)
          
          ! Impose minimum detectable signal:
          IF (lmds_z) THEN
            CALL mask_mindetectsignal ( rhvpolar, miss_value )
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

          ! .. LDR:
          !=========================================================================

          IF (loutpolall) THEN
             
            ! Sum up rs_meta(ista)%weigsm_v*rs_meta(ista)%weigsm_h*zvh_radar_mod_smth*pfcos_value:
            CALL weightedsum_beamfunction (zvh_radar_mod_smth, 0.0_dp, intgrl_pol)

            ALLOCATE(zvh_radar_mod_all(nobs),  stat=izs); iaerr = iaerr + ABS(izs)
            ! Compute zvh_radar_mod_all = intgrl_pol / intgrl_pf for all points where counter(irp) > 0, set miss_value otherwise:
            CALL divide_intpol_by_intpf (intgrl_pol, intgrl_pf, zvh_radar_mod_all)
            
            ALLOCATE(ldrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), stat=izs); iaerr = iaerr + ABS(izs)
            CALL init_vari(ldrpolar, miss_value)

            !.. Sort the collected radar points into the 3D field according
            !   to sorted azimut, range and elevation:
            DO iobs = 1, nobs
              IF (zh_radar_mod_all(iobs) >= Z_crit_radar) THEN
                ldrpolar( m_all(iobs),n_all(iobs),o_all(iobs) ) = zvh_radar_mod_all(iobs) / zh_radar_mod_all(iobs)
              END IF
            END DO

            ! Impose minimum detectable signal mask:
            IF (lmds_z) THEN
              CALL mask_mindetectsignal ( ldrpolar, miss_value )
            END IF
          
            CALL get_fileprefix_ascii_output ('ldrsim', rs_meta(ista), fileprefix)
            CALL output3d_ascii_radar(itime, &
                 rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                 ldrpolar, TRIM(ADJUSTL(fileprefix)), &
                 'Linear Depolarization Ratio', &
!!$                 '[dB]', 'polar', &
                 '[-]', 'polar', &
                 rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)
            
            IF (ldebug_radsim) THEN
              CALL control_output(&
                   ldroutputunit(ista), time_mod, rs_meta(ista), &
                   "LDR [-]", &
                   ldrpolar, (ldrpolar > miss_threshold), miss_value, nobs)
            END IF

          END IF
          
        END IF  ! (loutpolall .or. loutpolstd) .and. itype_refl > 4

      END IF ! loutdbz

      IF (ALLOCATED(intgrl_pol))        DEALLOCATE(intgrl_pol)
      IF (ALLOCATED(zh_radar_mod_all))  DEALLOCATE(zh_radar_mod_all)
      IF (ALLOCATED(zv_radar_mod_all))  DEALLOCATE(zv_radar_mod_all)
      IF (ALLOCATED(zvh_radar_mod_all)) DEALLOCATE(zvh_radar_mod_all)    
      IF (ALLOCATED(zv_radar_mod_all))  DEALLOCATE(zv_radar_mod_all)
      IF (ALLOCATED(zv_radar_mod_smth)) DEALLOCATE(zv_radar_mod_smth)

    ENDIF ! loutdbz .OR. (loutradwind .AND. (lweightdbz .OR. lmds_vr))


    ALLOCATE(vrpolar(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel), &
             stat=izs); iaerr = iaerr + ABS(izs)
    CALL init_vari(vrpolar, miss_value)

    IF (loutradwind) THEN

      ALLOCATE(intgrl_wind(nrp)); iaerr = iaerr + ABS(izs)
      ALLOCATE(radwind_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
               stat=izs); iaerr = iaerr + ABS(izs)

      CALL init_vari(intgrl_wind, 0.0_dp)
      CALL init_vari(radwind_mod_smth, miss_value)

      ! Re-initialize intgrl_pf for re-calculation for radwind, disregarding missing values::
      CALL init_vari(intgrl_pf, 0.0_dp)

      ALLOCATE(intgrl_vt(nrp)); iaerr = iaerr + ABS(izs)
      CALL init_vari(intgrl_vt, 0.0_dp)

      IF (lfall) THEN
        ALLOCATE(vt_mod_smth(nrp,rs_meta(ista)%ngpsm_v,rs_meta(ista)%ngpsm_h), &
                 stat=izs); iaerr = iaerr + ABS(izs)
        CALL init_vari(vt_mod_smth, 0.0_dp)
      END IF

!$omp parallel do
      DO ismth = 1, nsmth
        radwind_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = radwind_mod_all_smth(ismth)
        IF (ABS(radwind_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth))) < eps_vr) THEN
          radwind_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = 0.0_dp
        END IF
        IF (lfall) THEN
          vt_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = vt_mod_all_smth(ismth)
        END IF
      ENDDO
!$omp end parallel do

      IF ((ldealiase_vr_obs .AND. lreadmeta_from_netcdf) .OR. lfill_vr_backgroundwind) THEN
        ! save a radial wind field as a reference for dealiasing. This field is not affected
        !  by lweightdbz and lmds_vr below, so that we have a radial wind at those locations
        !  which do not belong to radar bins below the orography (lonline=.false.) or
        !  are at least not blocked by the orography (lonline=.true.).
        DO ih = 1, rs_meta(ista)%ngpsm_h
          DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do
            DO irp = 1, nrp
              IF (radwind_mod_smth(irp,iv,ih) >= miss_threshold) THEN
                intgrl_wind(irp) = intgrl_wind(irp) + &
                     rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*radwind_mod_smth(irp,iv,ih)* &
                     pfcos_value(irp,iv,ih)
                intgrl_pf(irp) = intgrl_pf(irp) + &
                     rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*pfcos_value(irp,iv,ih)
              END IF
            END DO
!$omp end parallel do
          END DO
        END DO

        ALLOCATE(vrpolar_for_dealiasing(rs_meta(ista)%naz,rs_meta(ista)%nra,rs_meta(ista)%nel))
        CALL init_vari(vrpolar_for_dealiasing, miss_value)
        DO iobs = 1, nobs
          irp = radpos_all(iobs)
          IF( intgrl_pf(irp) >= 1e-20_dp ) THEN
            vrpolar_for_dealiasing ( m_all(iobs),n_all(iobs),o_all(iobs) ) = &
                 intgrl_wind(irp) / intgrl_pf(irp)
          END IF
          IF (ABS(vrpolar_for_dealiasing ( m_all(iobs),n_all(iobs),o_all(iobs) )) < eps_vr)  THEN
            vrpolar_for_dealiasing ( m_all(iobs),n_all(iobs),o_all(iobs) ) = 0.0_dp
          END IF
        END DO
        ! re-initialize the used summation variables, they are also used below:
        CALL init_vari(intgrl_wind, 0.0_dp)
        CALL init_vari(intgrl_pf, 0.0_dp)

        ! and print to file if desired:
        CALL get_fileprefix_ascii_output ('vasim', rs_meta(ista), fileprefix)
        CALL output3d_ascii_radar(itime, &
             rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
             vrpolar_for_dealiasing, TRIM(ADJUSTL(fileprefix)), &
             'Simul. radial velocity, internally used for obs dealiasing', 'm/s', 'polar', &
             rmeta=rs_meta(ista), dbzmeta=dbz_meta(ista), is_aliased=.FALSE.)

        IF (ldebug_radsim) THEN
          CALL control_output(radwindoutputunit(ista), time_mod, rs_meta(ista), "V_r for dealiasing [m/s]", &
               vrpolar_for_dealiasing, (vrpolar_for_dealiasing > miss_threshold), miss_value, nobs)
        END IF

      END IF

      IF (lweightdbz) THEN

        DO ih = 1, rs_meta(ista)%ngpsm_h
          DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do
            DO irp = 1, nrp
              IF (counter(irp) > 0) THEN
                IF ( flag_smth(irp,iv,ih) > 0 .AND. zh_radar_mod_smth(irp,iv,ih) > 0.0_dp .AND. &
                     radwind_mod_smth(irp,iv,ih) > miss_threshold ) THEN
                  intgrl_wind(irp) = intgrl_wind(irp) + &
                       rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*radwind_mod_smth(irp,iv,ih)* &
                       pfcos_value(irp,iv,ih)*zh_radar_mod_smth(irp,iv,ih)
                  IF (lfall) THEN
                    intgrl_vt(irp) = intgrl_vt(irp) + &
                         rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*      &
                         SIN(el_loc_smth(irp,iv,ih)*degrad)* vt_mod_smth(irp,iv,ih)* &
                         pfcos_value(irp,iv,ih)*zh_radar_mod_smth(irp,iv,ih)
                  END IF
                END IF
              END IF
            ENDDO
!$omp end parallel do
          ENDDO
        ENDDO

      ELSE   ! .not. lweightdbz

        DO ih = 1, rs_meta(ista)%ngpsm_h
          DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do
            DO irp = 1, nrp
              IF (counter(irp) > 0) THEN
                IF ( flag_smth(irp,iv,ih) > 0 .AND. radwind_mod_smth(irp,iv,ih) > miss_threshold ) THEN
                  intgrl_wind(irp) = intgrl_wind(irp) + &
                       rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*radwind_mod_smth(irp,iv,ih)* &
                       pfcos_value(irp,iv,ih)
                  IF(lfall) THEN
                    intgrl_vt(irp) = intgrl_vt(irp) + &
                         rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*     &
                         SIN(el_loc_smth(irp,iv,ih)*degrad)*vt_mod_smth(irp,iv,ih)* &
                         pfcos_value(irp,iv,ih)
                  END IF
                  ! Re-calculate intgrl_pf. This time, weights for missing values are not counted,
                  ! which makes sense for radial wind:
                  intgrl_pf(irp) = intgrl_pf(irp) + &
                       rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*pfcos_value(irp,iv,ih)
                END IF
              END IF
            ENDDO
!$omp end parallel do
          ENDDO
        ENDDO

      END IF

      ALLOCATE(radwind_mod_all(nobs), stat=izs); iaerr = iaerr + ABS(izs)

      IF (lweightdbz) THEN

        DO iobs = 1, nobs

          irp = radpos_all(iobs)
          IF( intgrl_z(irp) < Z_crit_radar .OR. counter(irp) == 0 ) THEN
            radwind_mod_all(iobs) = miss_value
          ELSE
            ! In case of lfall=false, intgrl_vt = 0.0, so this equation covers both cases
            radwind_mod_all(iobs) = (intgrl_wind(irp) - intgrl_vt(irp))/intgrl_z(irp)
          ENDIF

          IF (ABS(radwind_mod_all(iobs)) < eps_vr) radwind_mod_all(iobs) = 0.0_dp

        ENDDO

      ELSE ! .not. lweightdbz

        DO iobs = 1, nobs

          irp = radpos_all(iobs)
          IF( counter(irp) > 0 .AND. intgrl_pf(irp) >= 1e-20_dp) THEN
            ! In case of lfall=false, intgrl_vt = 0.0, so this equation covers both cases
            radwind_mod_all(iobs) = (intgrl_wind(irp) -  intgrl_vt(irp))/intgrl_pf(irp)
          ELSE
            radwind_mod_all(iobs) = miss_value
          END IF

          IF (ABS(radwind_mod_all(iobs)) < eps_vr) radwind_mod_all(iobs) = 0.0_dp

        ENDDO

      ENDIF  ! lweightdbz


      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
      DO iobs = 1, nobs
        vrpolar( m_all(iobs),&
                 n_all(iobs),&
                 o_all(iobs) ) = radwind_mod_all(iobs)
      END DO

      IF (lmds_vr) THEN
        ! Take into account the minimum detectable signal (for simplicity assumed equal to that of the reflectivity):
!$omp parallel do private (o,n,m)
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
!$omp parallel do private (o,n,m)
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


    !.. Trap allocation errors:
    IF (iaerr /= 0) THEN
      CALL abort_run (my_radar_id, 39999, &
           'ERROR: allocation error in ' //TRIM(yzroutine)// &
           ', probably not enough memory! Decrease number of smoothing points! ', &
           TRIM(yzroutine)//', distribution of refractive index at radar stations')
    END IF

    !.. First memory cleanup:
    DEALLOCATE( radpos,                                            &
                iv_all, ih_all, flag_smth, intgrl_pf, intgrl_pf_nonblocked, pfcos_value, &
                counter, radpos_all,          &
                m_all, n_all, o_all )

    IF (ALLOCATED(mds_dbz))          DEALLOCATE(mds_dbz)

    IF(ALLOCATED(intgrl_z))           DEALLOCATE(intgrl_z)
    IF(ALLOCATED(zh_radar_mod_smth))  DEALLOCATE(zh_radar_mod_smth)

    IF(ALLOCATED(intgrl_wind))       DEALLOCATE(intgrl_wind)
    IF(ALLOCATED(radwind_mod_smth))  DEALLOCATE(radwind_mod_smth)
    IF(ALLOCATED(radwind_mod_all))   DEALLOCATE(radwind_mod_all)

    IF(ALLOCATED(intgrl_vt))         DEALLOCATE(intgrl_vt)
    IF(ALLOCATED(vt_mod_smth))       DEALLOCATE(vt_mod_smth)

    IF(ALLOCATED(erpolar))           DEALLOCATE(erpolar)
    IF(ALLOCATED(el_loc_smth))       DEALLOCATE(el_loc_smth)

    IF (lreadmeta_from_netcdf) THEN

      ! Dummy allocation of polarization parameters to avoid calling
      !  the next subroutine with a non-allocated argument. Some compilers don't like this:
      IF (.NOT. ALLOCATED(zdrpolar))   ALLOCATE(zdrpolar(0,0,0))
      IF (.NOT. ALLOCATED(kdppolar))   ALLOCATE(kdppolar(0,0,0))
      IF (.NOT. ALLOCATED(phidppolar)) ALLOCATE(phidppolar(0,0,0))
      IF (.NOT. ALLOCATED(rhvpolar))   ALLOCATE(rhvpolar(0,0,0))

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
            ! .. Else, construct the composite from one single elevation of the volume scans (not precip scan):
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
              ! .. Else, construct the composite from one single elevation:
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

    END IF
    
    !.. Second memory cleanup:
    IF (ALLOCATED(vrpolar))    DEALLOCATE(vrpolar)
    IF (ALLOCATED(zrpolar))    DEALLOCATE(zrpolar)
    IF (ALLOCATED(zdrpolar))   DEALLOCATE(zdrpolar)
    IF (ALLOCATED(kdppolar))   DEALLOCATE(kdppolar)
    IF (ALLOCATED(phidppolar)) DEALLOCATE(phidppolar)
    IF (ALLOCATED(ldrpolar))   DEALLOCATE(ldrpolar)
    IF (ALLOCATED(rhvpolar))   DEALLOCATE(rhvpolar)
    IF (ALLOCATED(srpolar))    DEALLOCATE(srpolar)
    IF (ALLOCATED(hrpolar))    DEALLOCATE(hrpolar)
    IF (ALLOCATED(lonpolar))   DEALLOCATE(lonpolar, latpolar)
    IF (ALLOCATED(vrpolar_for_dealiasing)) DEALLOCATE(vrpolar_for_dealiasing)


    IF (ldebug_radsim) THEN
      WRITE (*,'(a,i0,a,i6.6,a,i0)') TRIM(yzroutine)//': FINISHED output of radar station #', &
           ista, ' (ID: ', rs_meta(ista)%station_id, ') on proc ', my_radar_id
    END IF

  CONTAINS

    SUBROUTINE sort_ismth_to_radpos_iv_ih (field_radar_mod_all_smth, field_radar_mod_smth)
      REAL(kind=dp), INTENT(in)    :: field_radar_mod_all_smth(:)
      REAL(kind=dp), INTENT(inout) :: field_radar_mod_smth(:,:,:)
!$omp parallel do
      DO ismth = 1, nsmth
        field_radar_mod_smth(radpos(ismth),iv_all(ismth),ih_all(ismth)) = &
             field_radar_mod_all_smth(ismth)
      ENDDO
!$omp end parallel do
    END SUBROUTINE sort_ismth_to_radpos_iv_ih
    
    SUBROUTINE weightedsum_beamfunction (field_radar_mod_smth, thresh, intgrl_pol)
      REAL(kind=dp), INTENT(in)    :: field_radar_mod_smth(:,:,:)
      REAL(kind=dp), INTENT(in)    :: thresh
      REAL(kind=Dp), INTENT(inout) :: intgrl_pol(:)

      CALL init_vari(intgrl_pol, 0.0_dp)
      DO ih = 1, rs_meta(ista)%ngpsm_h
        DO iv = 1, rs_meta(ista)%ngpsm_v
!$omp parallel do private (irp)
          DO irp = 1, nrp
            IF (counter(irp) > 0) THEN
              IF (flag_smth(irp,iv,ih) > 0 .AND. field_radar_mod_smth(irp,iv,ih) > thresh) THEN
                intgrl_pol(irp) = intgrl_pol(irp) + &
                     rs_meta(ista)%weigsm_v(iv)*rs_meta(ista)%weigsm_h(ih)*&
                     field_radar_mod_smth(irp,iv,ih)*pfcos_value(irp,iv,ih)
              END IF
            END IF
          ENDDO
!$omp end parallel do
        ENDDO
      ENDDO
      
    END SUBROUTINE weightedsum_beamfunction

    SUBROUTINE divide_intpol_by_intpf (intgrl_pol, intgrl_pf, field_radar_mod_all)

      REAL(kind=dp), INTENT(in)    :: intgrl_pol(:), intgrl_pf(:)
      REAL(kind=dp), INTENT(inout) :: field_radar_mod_all(:)

!$omp parallel do private (iobs,irp)
      DO iobs = 1, nobs
        irp = radpos_all(iobs)
        ! Assume linear field values:
        IF (counter(irp) > 0) THEN
          field_radar_mod_all(iobs) = intgrl_pol(irp) / intgrl_pf(irp)
        ELSE
          field_radar_mod_all(iobs) = miss_value
        END IF
      ENDDO
!$omp end parallel do
    END SUBROUTINE divide_intpol_by_intpf

    SUBROUTINE sort_into_polarcoords (field_radar_mod_all, thresh, field_polar)
      REAL(kind=dp), INTENT(in)    :: field_radar_mod_all(:)
      REAL(kind=dp), INTENT(in)    :: thresh
      REAL(kind=dp), INTENT(inout) :: field_polar(:,:,:) 
      
      !.. Sort the collected radar points into the 3D field according
      !   to sorted azimut, range and elevation:
!$omp parallel do private (iobs)
      DO iobs = 1, nobs
        IF (field_radar_mod_all(iobs) >= thresh) THEN
          field_polar( m_all(iobs),n_all(iobs),o_all(iobs) ) = field_radar_mod_all(iobs)
        END IF
      END DO
!$omp end parallel do
      
    END SUBROUTINE sort_into_polarcoords
    
    ! Take into account the minimum detectable signal:
    ! ------------------------------------------------
    !  If the horizontal reflectivity zrpolar(:,:,:) is below a range-dependend threshold mds_dbz(range),
    !  then zfield will be assumed not estimable and set to the fill_value.
    ! NOTE: zrpolar(:,:,:) and mds_dbz(:) must be allocated and their values properly defined!
    SUBROUTINE mask_mindetectsignal (zfield, fill_value)
      REAL(kind=dp), INTENT(inout) :: zfield (:,:,:)
      REAL(kind=dp), INTENT(in)    :: fill_value
      INTEGER :: o, n, m
!$omp parallel do private(o,n,m)
      DO o=1, rs_meta(ista)%nel
        DO n=1, rs_meta(ista)%nra
          DO m=1, rs_meta(ista)%naz
            IF (zfield(m,n,o) >= miss_threshold .AND. zrpolar(m,n,o) < mds_dbz(n) ) THEN
              ! reflectivity is below the MDS, so set zfield to a fill_value:
              zfield(m,n,o) = fill_value
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do
    END SUBROUTINE mask_mindetectsignal
    
  END SUBROUTINE output_my_ista_smth

END MODULE radar_output_station_smth
