! Source module for the radar forward operator EMVORADO
!
! ---------------------------------------------------------------
! Copyright (C) 2017-2024, DWD
! Contact information: ulrich.blahak (at) dwd.de 
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE radar_composites

!------------------------------------------------------------------------------
!
! Description:
!   This module provides some routines which are necessary for generating
!   radar composites on a rotated lat/lon grid in the radar forward operator EMVORADO.
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
  
  USE radar_data, ONLY :  &
       miss_threshold, miss_value, zero_value, Z_crit_radar, dBZ_crit_radar, &
       radar_meta_type, &
       num_radar, icomm_radario_dom, icomm_radar, radar_master, &
       num_radario, lradario_pe_dom, lcompute_pe_fwo, &
       my_radar_id, my_radar_id_dom, my_radario_id_dom, &
       composite_meta_type, degrad

  USE radar_utilities, ONLY : &
       phirot2phi,            &
       rlarot2rla,            &
       phi2phirot,            &
       rla2rlarot,            &
       init_vari,             &
       geo_heading, geo_dist, &
       get_range2, smoother

  USE radar_interface, ONLY : &
       abort_run

 !==============================================================================

#ifndef NOMPI
  USE mpi
#endif

!==============================================================================

  IMPLICIT NONE

!==============================================================================

! default private
  PRIVATE

!==============================================================================

#ifdef NOMPI
  INCLUDE "nompi_mpif.h"
#endif

!==============================================================================

  ! These global composite fields are needed to build
  !  the composites in the radar forward operator.
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: &
       comp_dbzsim_tot, comp_dbzobs_tot
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:), TARGET :: &
       comp_dbzsim_bub_tot, comp_dbzobs_bub_tot

PUBLIC comp_dbzsim_tot, comp_dbzobs_tot, comp_dbzsim_bub_tot, comp_dbzobs_bub_tot

!==============================================================================

PUBLIC comp2geo_coord, geo2comp_coord, comp_cellindex2geo, geo2comp_cellindex, &
       alloc_composite, dealloc_composite, &
       alloc_composite_bub, dealloc_composite_bub, &
       composite2D_dbz_maxmethod_ista, collect_smooth_composite


!==============================================================================
! Interface Blocks for overloaded procedures:
!==============================================================================

INTERFACE comp_cellindex2geo
  MODULE PROCEDURE            &
       comp_cellindex2geo_int,           &
       comp_cellindex2geo_real
END INTERFACE comp_cellindex2geo

INTERFACE geo2comp_cellindex
  MODULE PROCEDURE            &
       geo2comp_cellindex_int,           &
       geo2comp_cellindex_real
END INTERFACE geo2comp_cellindex


!==============================================================================
!==============================================================================

CONTAINS

!==============================================================================
!==============================================================================


  !============================================================================
  !
  ! Subroutines related to composite generation on a rotated lat/lon grid
  !
  !============================================================================

  SUBROUTINE comp2geo_coord (x_model, y_model, comp_meta, lon_geo, lat_geo)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)          :: x_model, y_model    ! x and y-coord of the rotated lon/lat grid
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    REAL(KIND=dp), INTENT(out)         :: lon_geo, lat_geo    ! geogr. lon/lat

    lon_geo = rlarot2rla(y_model,x_model,comp_meta%pollat,comp_meta%pollon,comp_meta%polgam)
    lat_geo = phirot2phi(y_model,x_model,comp_meta%pollat,comp_meta%pollon,comp_meta%polgam)


  END SUBROUTINE comp2geo_coord

  SUBROUTINE geo2comp_coord (lon_geo, lat_geo, comp_meta, x_model, y_model)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)         :: lon_geo, lat_geo    ! geogr. lon/lat
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    REAL(KIND=dp), INTENT(out)        :: x_model, y_model    ! x and y-coord of the rotated lon/lat grid

    x_model = rla2rlarot(lat_geo,lon_geo,comp_meta%pollat,comp_meta%pollon,comp_meta%polgam) ! rotated longitude
    y_model = phi2phirot(lat_geo,lon_geo,comp_meta%pollat,comp_meta%pollon)               ! rotated latitude


  END SUBROUTINE geo2comp_coord

  SUBROUTINE comp_cellindex2geo_int (i_model, j_model, comp_meta, lon_geo, lat_geo)

    IMPLICIT NONE

    INTEGER, INTENT(in)   :: i_model, j_model  ! horizontal index of rotated lon/lat grid cell
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    REAL(KIND=dp), INTENT(out)        :: lon_geo, lat_geo  ! geogr. lon/lat of cell center


    REAL(KIND=dp) :: x_model, y_model

    x_model = comp_meta%startlon + (i_model-1) * comp_meta%dlon
    y_model = comp_meta%startlat + (j_model-1) * comp_meta%dlat
    CALL comp2geo_coord (x_model, y_model, comp_meta, lon_geo, lat_geo)



  END SUBROUTINE comp_cellindex2geo_int

  SUBROUTINE comp_cellindex2geo_real (i_model, j_model, comp_meta, lon_geo, lat_geo)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)   :: i_model, j_model  ! horizontal index of rotated lon/lat grid cell
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    REAL(KIND=dp), INTENT(out)        :: lon_geo, lat_geo  ! geogr. lon/lat of cell center


    REAL(KIND=dp) :: x_model, y_model

    x_model = comp_meta%startlon + (i_model-1.0_dp) * comp_meta%dlon
    y_model = comp_meta%startlat + (j_model-1.0_dp) * comp_meta%dlat
    CALL comp2geo_coord (x_model, y_model, comp_meta, lon_geo, lat_geo)



  END SUBROUTINE comp_cellindex2geo_real


  SUBROUTINE geo2comp_cellindex_int (lon_geo, lat_geo, comp_meta, i_model, j_model)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)         :: lon_geo, lat_geo  ! geogr. lon/lat
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    INTEGER, INTENT(out)  :: i_model, j_model  ! horizontal index of rotated lat/lon grid cell
                                                               !  which contains lon_geo/lat_geo

    REAL(KIND=dp) :: x_model, y_model

    CALL geo2comp_coord (lon_geo, lat_geo, comp_meta, x_model, y_model)
    i_model = NINT((x_model-comp_meta%startlon) / comp_meta%dlon) + 1           ! For ICON: "running" horizontal index
    j_model = NINT((y_model-comp_meta%startlat) / comp_meta%dlat) + 1           ! For ICON: perhaps just set to 1


  END SUBROUTINE geo2comp_cellindex_int

   SUBROUTINE geo2comp_cellindex_real (lon_geo, lat_geo, comp_meta, i_model, j_model)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)         :: lon_geo, lat_geo  ! geogr. lon/lat
    TYPE(composite_meta_type), INTENT(in) :: comp_meta
    REAL(KIND=dp), INTENT(out)  :: i_model, j_model  ! horizontal index of rotated lat/lon grid cell
                                                               !  which contains lon_geo/lat_geo

    REAL(KIND=dp) :: x_model, y_model

    CALL geo2comp_coord (lon_geo, lat_geo, comp_meta, x_model, y_model)
    i_model = ((x_model-comp_meta%startlon) / comp_meta%dlon) + 1.0_dp          ! For ICON: "running" horizontal index
    j_model = ((y_model-comp_meta%startlat) / comp_meta%dlat) + 1.0_dp          ! For ICON: perhaps just set to 1


  END SUBROUTINE geo2comp_cellindex_real


  SUBROUTINE alloc_composite ( ldebug, comp_meta, nel, lalloc_obs )

    ! Allocates the total-domain (ie_tot,je_tot) composites (obs and sim) which
    !  are needed on the radar output-PEs and the root node to compute the composites.
    !  These are then splitted into the local PE-domains and distributed to all compute PEs.
    !
    ! call from beginning of 'compute' of organize_radar on all PEs
    !  (compute- and output-PEs, IF async)

    IMPLICIT NONE

    LOGICAL, INTENT(in)  :: ldebug, lalloc_obs
    TYPE(composite_meta_type), INTENT(in)  :: comp_meta
    INTEGER, INTENT(in)  :: nel

    IF (ldebug) WRITE(*,*) 'alloc_composite on proc ', my_radar_id

    ALLOCATE ( comp_dbzsim_tot(comp_meta%ni,comp_meta%nj,nel) )
    ! initialize composites with missing value:
    CALL init_vari(comp_dbzsim_tot(:,:,1:nel) , miss_value)

    IF (lalloc_obs) THEN
      ALLOCATE ( comp_dbzobs_tot(comp_meta%ni,comp_meta%nj,nel) )
      CALL init_vari(comp_dbzobs_tot(:,:,1:nel) , miss_value)
    END IF

    ! will be deallocated in deallocate_composit() at the end of each radar timestep

    IF (ldebug) WRITE(*,*) 'Done with alloc_composite on proc ', my_radar_id

  END SUBROUTINE alloc_composite

  SUBROUTINE alloc_composite_bub ( ldebug, comp_meta )

    ! Allocates the total-domain (ie_tot,je_tot) composites (obs and sim) which
    !  are needed on the radar output-PEs and the root node to compute the composite for the bubble generator.
    !  These are then splitted into the local PE-domains and distributed to all compute PEs.
    !
    ! call from beginning of 'compute' of organize_radar on all PEs
    !  (compute- and output-PEs, IF async)

    IMPLICIT NONE

    LOGICAL, INTENT(in)  :: ldebug
    TYPE(composite_meta_type), INTENT(in)  :: comp_meta

    IF (ldebug) WRITE(*,*) 'alloc_composite on proc ', my_radar_id

    ALLOCATE ( comp_dbzsim_bub_tot(comp_meta%ni,comp_meta%nj) , &
               comp_dbzobs_bub_tot(comp_meta%ni,comp_meta%nj) )

    ! initialize composites with missing value:
    CALL init_vari(comp_dbzsim_bub_tot(:,:) , miss_value)
    CALL init_vari(comp_dbzobs_bub_tot(:,:) , miss_value)

    ! will be deallocated in deallocate_composit() at the end of each radar timestep

    IF (ldebug) WRITE(*,*) 'Done with alloc_composite on proc ', my_radar_id

  END SUBROUTINE alloc_composite_bub

  ! Make a composite of polar dbz-data set "zpolar" using data of the elevation with index
  !  "eleind_for_composite". The radar meta data of station ista have to be provided in the struct rsm.
  !  The composite is made on a rotated lat/lon grid and overlayed (by a maximum operation) on
  !  to the INOUT 2D field compdata, which has to be a total-domain 2D model grid field
  !  (normally one of comp_dbzsim or comp_dbzobs).
  ! The metadata of the rotated lat/lon grid have to be provided on input by the
  !  type(composite_meta_type) comp_meta.

  SUBROUTINE composite2D_dbz_maxmethod_ista (zpolar, rsm, eleind_for_composite, compdata_tot, ldebug, comp_meta)

    ! call from output_radar_country_obs() for each single radar station ista on the output PEs

    IMPLICIT NONE

    ! Input variables:
    REAL(KIND=dp)      , INTENT(in)    :: zpolar(:,:,:)
    TYPE(radar_meta_type)  , INTENT(in)    :: rsm
    INTEGER, INTENT(in)    :: eleind_for_composite
    REAL(kind=dp)       , INTENT(inout)    :: compdata_tot(:,:)  ! either comp_dbzsim_tot or comp_dbzobs_tot
    LOGICAL, INTENT(in)                    :: ldebug
    TYPE(composite_meta_type), INTENT(in)  :: comp_meta

    ! Local variables:
    REAL(kind=dp), ALLOCATABLE       :: compsum(:,:)
    REAL(kind=dp)                    :: idispl, jdispl, latrot, lonrot, lat, lon, &
                                        azi, s, ele, range, span_i, span_j, tmpdat, ln10_o10
    INTEGER, ALLOCATABLE :: compnr(:,:)
    INTEGER          :: i, j, iover, jover, nover, ni, nj, &
                                        i_start, i_end, j_start, j_end, iaz_near, ira_near, &
                                        i_nearest_sta, j_nearest_sta

    CHARACTER(LEN=50)       :: yzroutine

    yzroutine(:) = ' '
    yzroutine    = 'composite2D_dbz_maxmethod_ista'

    IF (ldebug) WRITE(*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    IF (eleind_for_composite < 1 .OR. (eleind_for_composite > rsm%nel .AND. eleind_for_composite /= 99) ) THEN
      WRITE(*,'(a,i3,a,i6.6,a)') 'INFO '//TRIM(yzroutine)//': elevation index ',  eleind_for_composite, &
           ' for composite out of range and discarded for station ', rsm%station_id, ' with scanname '//TRIM(rsm%scanname)
      RETURN
    END IF

    ln10_o10 = LOG(10._dp)*0.1_dp

    ni = SIZE(compdata_tot, dim=1)  ! comp_meta%ni
    nj = SIZE(compdata_tot, dim=2)  ! comp_meta%nj

    ! nearest lat/lon grid point to the radar station:
    CALL geo2comp_cellindex (rsm%lon, rsm%lat, comp_meta, &
                             i_nearest_sta, j_nearest_sta)

    ! helper fields to compute oversampled average values for each composite point:
    ALLOCATE(compsum(ni, nj), compnr(ni, nj))
    CALL init_vari(compsum , 0.0_dp)
    CALL init_vari(compnr  , 0     )

    ! Oversampling points of the lat/lon grid:
    !  Set number of points in such a way that at least every 0.5*range_gate there is an
    !  oversampling point:
    nover = CEILING( comp_meta%r_earth * degrad * MAX(comp_meta%dlon, comp_meta%dlat) / (0.5_dp * rsm%ra_inc) )
    ! make sure that nover is an odd number:
    nover = nover + MOD(nover+1, 2)
    ! max. nover is 7 (to limit comp. burden):
    nover = MIN(nover, 7)

    ! limit the oversampled range of the lat/lon domain to the (almost) smallest rectangle
    !  which surrounds the radar station and PPI scan:
    latrot = comp_meta%startlat + (j_nearest_sta-1) * comp_meta%dlat  ! approx. rotated radar latitude
    span_i  = 1.2 * rsm%nra * rsm%ra_inc / (comp_meta%r_earth*degrad*COS(latrot*degrad)*comp_meta%dlon)
    span_j  = 1.2 * rsm%nra * rsm%ra_inc / (comp_meta%r_earth*degrad*comp_meta%dlat)
    i_start = i_nearest_sta - CEILING(span_i) - 1
    i_end   = i_nearest_sta + CEILING(span_i) + 1
    j_start = j_nearest_sta - CEILING(span_j) - 1
    j_end   = j_nearest_sta + CEILING(span_j) + 1

    ! limit the range to the actual lat/lon domain, excluding the outermost row to
    !  save an if-statement in the loops below:
    i_start = MAX(i_start, 2       )
    i_end   = MIN(i_end,   comp_meta%ni-1)
    j_start = MAX(j_start, 2       )
    j_end   = MIN(j_end,   comp_meta%nj-1)

    ! Oversampling by
    DO jover = 0, nover-1
      jdispl = -0.5_dp + (jover+0.5_dp) / nover
      DO iover = 0, nover-1
        idispl = -0.5_dp + (iover+0.5_dp) / nover

!$omp parallel do private(i,j,lonrot,latrot,lon,lat,azi,s,ele,range,ira_near,iaz_near,tmpdat)
        DO j = j_start, j_end
          DO i = i_start, i_end

            ! rotatet lon/lat of oversampled lat/lon grid point:
            lonrot = comp_meta%startlon + (i-1+idispl) * comp_meta%dlon
            latrot = comp_meta%startlat + (j-1+jdispl) * comp_meta%dlat

            ! convert rotated to geographic coords:
            CALL comp2geo_coord (lonrot, latrot, comp_meta, lon, lat)
!!$ The same as:
!!$            lon = rlarot2rla(latrot,lonrot,comp_meta%pollat,comp_meta%pollon,comp_meta%polgam)
!!$            lat = phirot2phi(latrot,lonrot,comp_meta%pollat,comp_meta%pollon,comp_meta%polgam)

            ! compute azimut, arc length s and radar range from radar station
            !  to model point:
            azi   = geo_heading(rsm%lon, rsm%lat, lon, lat)
            s     = geo_dist   (rsm%lon, rsm%lat, lon, lat, comp_meta%r_earth, rsm%alt_msl)
            ele   = rsm%el_arr(eleind_for_composite)
            range = get_range2 (s, ele, comp_meta%r_earth, rsm%alt_msl)

            ! nearest neighbour interpolation from zpolar to point azi, range
            ! average over all oversampled points,
            ! maximum of pre-existing value in composite and new value:
            ira_near = NINT(range / rsm%ra_inc)
            iaz_near = MODULO( NINT( (azi-rsm%az_start) / rsm%az_inc ), rsm%naz) + 1

            IF ( 1 <= ira_near .AND. ira_near <= rsm%nra .AND. &
                 1 <= iaz_near .AND. iaz_near <= rsm%naz ) THEN
              tmpdat = zpolar(iaz_near,ira_near,eleind_for_composite)
              IF (tmpdat >= dBZ_crit_radar) THEN
                compsum(i,j) = compsum(i,j) + EXP(ln10_o10*tmpdat)
              END IF
              IF (tmpdat >= miss_threshold) THEN
                compnr (i,j)  = compnr(i,j) + 1
              END IF
            END IF

          END DO
        END DO
!$omp end parallel do

      END DO
    END DO

    ! Maximum of actual mean and pre-existing value:
!$omp parallel do private(i,j)
    DO j = j_start, j_end
      DO i = i_start, i_end
        IF (compnr(i,j) > 0) THEN
          IF (compsum(i,j) >= Z_crit_radar) THEN
            compdata_tot(i,j) = MAX(compdata_tot(i,j), 10.0_dp*LOG10(compsum(i,j)/compnr(i,j)) )
          ELSE
            compdata_tot(i,j) = MAX(compdata_tot(i,j), zero_value)
          END IF
        END IF
      END DO
    END DO
!$omp end parallel do

    DEALLOCATE(compsum, compnr)

    IF (ldebug) WRITE(*,*) 'Done with '//TRIM(yzroutine)//' on proc ', my_radar_id, MAXVAL(compdata_tot)

  END SUBROUTINE composite2D_dbz_maxmethod_ista

  SUBROUTINE collect_smooth_composite ( idom_model, compdata_tot, comp_meta, lsmooth_composite, &
                                           nsmoothpoints_for_comp, nfilt_for_comp, ldebug )

    ! call from end of organize_radar
    ! call on all PEs (compute and output)

    ! collect composite from all output PEs on , filter/smooth it and distribute it to all compute PEs:
    !  - either use MPI_ALLREDUCE(MAX) on all PEs (non-output-PEs were initialized with miss_value)
    !  - or create point-to-point communication using MPI_SEND to PE 0 and a 3D-field (ie_tot, je_tot, anz_output_PEs)
    !    to compute the MAX over all output PEs. For this, you can use a flag which has been set on the
    !     output PEs and not on the other PEs. Then, distribute to all PEs using MPI_BCAST().

    IMPLICIT NONE

    REAL(kind=dp), INTENT(inout)  :: &
         compdata_tot(:,:)    ! either comp_dbzsim_tot or comp_dbzobs_tot
    INTEGER, INTENT(in):: &
         idom_model, &              ! domain identifier in the hosting model [1-ndoms]
         nsmoothpoints_for_comp, &  ! width of symetric 2D binomial smoother in grid points
         nfilt_for_comp             ! number of consecutive applications of the smoother
    TYPE(composite_meta_type), INTENT(in)  :: comp_meta
    LOGICAL      , INTENT(in)     :: ldebug, lsmooth_composite

    INTEGER         :: ni, nj, &
         icomm_for_mpi,          &  ! communicator on which the composites have been built (icomm_radar or icomm_radario)
         ipe, irecv_for_mpi
    REAL(kind=dp), ALLOCATABLE  :: sendbuf(:,:), recvbuf(:,:)
    REAL(kind=dp)               :: tmp, ln10_o10

    CHARACTER(LEN=50)     :: yzroutine
    CHARACTER(LEN=255)    :: yerrmsg
    INTEGER               :: ierror, i, j


    yzroutine(:) = ' '
    yzroutine    = 'distribute_smooth_composite'
    yerrmsg(:)   = ' '
    ierror       = 0

    IF (ldebug) WRITE(*,*) TRIM(yzroutine)//' on proc ', my_radar_id

    ni = SIZE(compdata_tot, dim=1)  ! comp_meta%ni
    nj = SIZE(compdata_tot, dim=2)  ! comp_meta%nj

    ln10_o10 = LOG(10._dp)*0.1_dp

    IF (num_radar > 1) THEN

      IF (num_radario > 0) THEN
        icomm_for_mpi = icomm_radario_dom(idom_model)
        irecv_for_mpi = 0                        ! target PE in icomm_radario where to collect the composites
        ipe           = my_radario_id_dom(idom_model)  ! ID of irecv_for_mpi in in comm. icomm_radario_dom
      ELSE
        icomm_for_mpi = icomm_radar
        irecv_for_mpi = radar_master           ! = 0
        ipe           = my_radar_id_dom(idom_model)  ! = 0
      END IF

      IF ((num_radario > 0 .AND. lradario_pe_dom(idom_model)) .OR. (num_radario == 0 .AND. lcompute_pe_fwo)) THEN

        ALLOCATE(sendbuf(ni,nj))
!$omp parallel do private(i,j)
        DO j=1, nj
          DO i=1, ni
            IF (compdata_tot(i,j) > dBZ_crit_radar) THEN
              sendbuf(i,j) = EXP(ln10_o10*compdata_tot(i,j))
            ELSE IF (compdata_tot(i,j) > miss_threshold) THEN
              sendbuf(i,j) = 0.0_dp
            ELSE
              sendbuf(i,j) = compdata_tot(i,j)    ! some missing value
            END IF
          END DO
        END DO
!$omp end parallel do

        IF (irecv_for_mpi == ipe) THEN
          ALLOCATE(recvbuf(ni,nj))
        ELSE
          ! just to enable compilation with bounds checking:
          ALLOCATE(recvbuf(0:0,0:0))
        END IF

!!$ Alberto: MPI_ALLREDUCE wäre schneller, braucht lediglich die Allokierung von recvbuf(ni,nj) auf allen PEs.
        CALL MPI_REDUCE(sendbuf, recvbuf, ni*nj, MPI_DOUBLE_PRECISION, MPI_MAX, irecv_for_mpi, icomm_for_mpi, ierror)
        IF ( ierror /= 0 ) THEN
          WRITE (yerrmsg,'(a)') 'Error in MPI_REDUCE() of composites'
          CALL abort_run(my_radar_id, ierror, yerrmsg, yzroutine)
        ENDIF
        DEALLOCATE(sendbuf)

        IF (irecv_for_mpi == ipe) THEN

          !.. Some spatial filter operations on compdata_tot (apply binomial filter nfilt times):
          IF (lsmooth_composite) THEN
            CALL smoother(f_2d_field=recvbuf, ie=comp_meta%ni, je=comp_meta%nj, &
                 nlength=INT(nsmoothpoints_for_comp, KIND(1)), nfilt=INT(nfilt_for_comp, KIND(1)) )
          END IF

!$omp parallel do private(i,j)
        DO j=1, nj
          DO i=1, ni
            IF (recvbuf(i,j) >= Z_crit_radar) THEN
              compdata_tot(i,j) = 10.0_dp*LOG10(recvbuf(i,j))
            ELSE IF (recvbuf(i,j) >= miss_threshold) THEN
              compdata_tot(i,j) = zero_value   ! correct 0
            END IF
          END DO
        END DO
!$omp end parallel do
        END IF
        DEALLOCATE(recvbuf)

        ! The composite lives now on PE 0 of the radario-group on a 2D total domain field.

      END IF

    ELSE   ! num_compute_fwo == num_radar == 1

!$omp parallel do private(i,j)
      DO j=1, nj
        DO i=1, ni
          IF (compdata_tot(i,j) > dBZ_crit_radar) THEN
            compdata_tot(i,j) = EXP(ln10_o10*compdata_tot(i,j))
          ELSEIF (compdata_tot(i,j) > miss_threshold) THEN
            compdata_tot(i,j) = 0.0_dp
          END IF
        END DO
      END DO
!$omp end parallel do

      IF (lsmooth_composite) THEN
        CALL smoother(f_2d_field=compdata_tot, ie=comp_meta%ni, je=comp_meta%nj, &
                      nlength=INT(nsmoothpoints_for_comp, KIND(1)), nfilt=INT(nfilt_for_comp, KIND(1)) )
      END IF

      ! The composite lives now on PE 0 on a 2D total domain field.

      DO j=1, nj
        DO i=1, ni
          IF (compdata_tot(i,j) >= Z_crit_radar) THEN
            compdata_tot(i,j) = 10.0_dp*LOG10(compdata_tot(i,j))
          ELSEIF (compdata_tot(i,j) >= miss_threshold) THEN
            compdata_tot(i,j) = zero_value   ! correct 0
          END IF
        END DO
      END DO

    END IF   ! num_radar == 1

    IF (ldebug) WRITE(*,*) 'Done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE collect_smooth_composite

  SUBROUTINE dealloc_composite

    IF (ALLOCATED(comp_dbzobs_tot)) DEALLOCATE(comp_dbzobs_tot)
    IF (ALLOCATED(comp_dbzsim_tot)) DEALLOCATE(comp_dbzsim_tot)

  END SUBROUTINE dealloc_composite

  SUBROUTINE dealloc_composite_bub

    IF (ALLOCATED(comp_dbzobs_bub_tot)) DEALLOCATE(comp_dbzobs_bub_tot)
    IF (ALLOCATED(comp_dbzsim_bub_tot)) DEALLOCATE(comp_dbzsim_bub_tot)

  END SUBROUTINE dealloc_composite_bub

!================================================================================

END MODULE radar_composites
