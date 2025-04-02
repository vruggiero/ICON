!!! Seems to be not necessary, therefore commented out:
!!!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

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

MODULE radar_model2rays

!------------------------------------------------------------------------------
!
! Description: Module in radar forward operator EMVORADO to provide the functionalities
!              for interpolation/aggregation
!              of model variables/radar observables to the rays/bins of simulated
!              radar stations. This module contains higher-level routines for
!              for this task considering different configuration options from
!              EMVORADO namelist and uses low-level procedures from
!              other modules.
!              The procedures provided herein are called mainly from radar_organize.f90.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind, ONLY : dp, wp

  USE radar_data, ONLY : &
       miss_threshold, miss_value, miss_thresh_rhv, miss_value_rhv, &
       missval_int, shield_value, shield_value_rhv

  USE radar_interface, ONLY : &    ! This makes available the prognostic model variables and config!
       u     ,   & ! U
       v     ,   & ! V
       w     ,   & ! W
       t     ,   & ! T
       p     ,   & ! P
       rho   ,   & ! RHO
       qv    ,   & ! QV
       vapor_pres, & ! rho*qv*r_v*t
       qc    ,   & ! QC
       qi    ,   & ! QI
       qr    ,   & ! QR
       qs    ,   & ! QS
       qg    ,   & ! QG
       qh    ,   & ! QH
       qnc   ,   & ! NCCLOUD
       qni   ,   & ! NCICE
       qnr   ,   & ! NCRRAIN
       qns   ,   & ! NCSNOW
       qng   ,   & ! NCGRAUPEL
       qnh   ,   & ! NCHAIL
       qgl   ,   & ! QGL
       qhl   ,   & ! QHL
       qnc_s ,   & ! QNC_S
       zh_radar, ah_radar, &
       zv_radar, rrhv_radar, irhv_radar, kdp_radar, adp_radar, zvh_radar, &
       vt_radar, hfl


  USE radar_interface, ONLY : &
       abort_run, nlevels,                      &
       bottomlevel, bottomlevel_stag, toplevel, levelincr, &
       one_level_up, one_level_down, &
       get_runtime_timings, &
#ifdef TWOMOM_SB
       initialize_tmax_2mom_vec_par,                   &
#endif
       initialize_tmax_1mom_vec_par, finalize_tmax

  USE radar_interface,     ONLY :                                                  &
       setup_model2radarbins_vec, setup_model2azislices_vec,                       &
       interp_model2radarbins_scalar, interp2d_model2radarbins_scalar,             &
       interp_model2radarbins_vr,                                                  &
       interp_model2azislices_scalar, interp_model2azislices_vr,                   &
       interp2D_model2geo_horiz_scalar, interp3D_model2geo_scalar

   USE radar_mie_iface_cosmo_driver, ONLY : &
        calc_dbz_vec_modelgrid, calc_fallspeed_vec_modelgrid, init_lookup_mie,     &
        calc_dbz_vec_generic

 !------------------------------------------------------------------------------

   USE radar_utilities, ONLY : el_loc_43, refr_index_air,  &
                               polar2geo, polar2geo_old,   &
                               ind2sub2D, sub2ind2D,       &
                               ind2sub3D, sub2ind3D,       &
                               ind2sub4D, sub2ind4D,       &
                               ind2sub5D, sub2ind5D,       &
                               phirot2phi,                 &
                               rlarot2rla,                 &
                               phi2phirot,                 &
                               rla2rlarot,                 &
                               smth_el_horzscan,           &
                               smth_az_horzscan,           &
                               init_vari,                  &
                               set_missing_and_correct0,   &
                               dbz_to_linear

  !------------------------------------------------------------------------------


  USE radar_parallel_utilities, ONLY :  distribute_values_radar

  !------------------------------------------------------------------------------

  USE radar_data,               ONLY :                                                         &
       nradsta_max,                                                    &
       radar_meta_type,              &
       rpvect, ipvect,                                                                         &
       ydate_ini_mod,                                                                          &
       idom,   &
       imp_fwo_double,        & ! determines the correct REAL type used in the model for MPI
       imp_fwo_integers,      & ! determines the correct INTEGER type used in the model for MPI
       num_compute_fwo,   & ! number of compute PEs
       my_cart_id_fwo,    & ! rank of this PE (=subdomain) in the cartesian communicator
       my_radar_id,       & ! rank of this PE in the radar communicator (cart+radario)
       icomm_cart_fwo,    & ! communicator for the virtual cartesian topology
       i_fwo_compgrid,    & ! Timing flag for computations on the model grid
       i_fwo_comppolar,   & ! Timing flag for interpolation of reflectivity and radial wind
                            !  from model grid points to the radar bins/auxiliary azi slice grid
       degrad,             &
       r_earth_dp,       & ! mean radius of the earth
       rs_meta, rs_data, rs_grid, dbz_meta, &
       nradsta, itlrad_dyn, itlrad_qx, nbl_az

  USE radar_data_namelist, ONLY :  &
       ldebug_radsim, lsode, &
       loutpolstd, loutpolall, lextdbz, lweightdbz, &
       ydir_mielookup_read, ydir_mielookup_write


  USE radar_data_mie, ONLY : &
       Tmax_i_modelgrid, Tmax_s_modelgrid, Tmax_g_modelgrid, Tmax_h_modelgrid, &
       Tmin_g_modelgrid, Tmin_h_modelgrid

  USE radar_obs_meta_list, ONLY : get_elarr_precipscan

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

  PUBLIC :: calc_geometry_grid, &
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
            distr_onlinf_all_reals_all2allv, &
            get_azvec, &
            online2geo, &
            rad2geo_const

  !==============================================================================

  !==============================================================================
  ! Interface blocks for overloaded procedures:

  !==============================================================================
  ! Public and Private Subroutines

  !==============================================================================
  ! Module procedures
  !==============================================================================

CONTAINS

  !==============================================================================
  !+ Module procedure in radar_src for the computation of radar beam geometry
  !------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_vec

    !------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar beam geometry and the interpolation weights
    !              and interpolation indices for the interpolation of the
    !              model data to the radar points. These interpolation
    !              weights and the interpolation indices are constant in time since
    !              the radar geometry is also constant in time and this subroutine
    !              is only called in the first call of the operator.
    !
    ! Method:
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER :: ista, iaz, ira, iel, nobsmax, irp, n, m, o

    REAL    (KIND=dp), ALLOCATABLE :: &
         raarr(:)        ,& ! array of radial bins
         el(:,:)         ,& ! array of elevations
         el_precip(:)    ,& ! array of elevations
         lat_r(:)        ,& ! array of geographical latitudes for each radar point (nrp)
         lon_r(:)        ,& ! array of geographical longitudes for each radar point (nrp)
         alt_r(:)           ! array of altitudes above sea level for each radar point (nrp)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_vec
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_vec'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! allocate and initialize local aux arrays
      nobsmax = rs_meta(ista)%nra * rs_meta(ista)%naz * rs_meta(ista)%nel

      ALLOCATE(lat_r(nobsmax))
      ALLOCATE(lon_r(nobsmax))
      ALLOCATE(alt_r(nobsmax))

!CDIR BEGIN COLLAPSE
      lat_r    = 0.0_dp
      lon_r    = 0.0_dp
      alt_r    = -1.0_dp
!CDIR END

      ! calculate geographical coordinates (lat,lon,alt) of radar points
      CALL rad2geo_const_vec( rs_meta(ista), lat_r, lon_r, alt_r )

      ! Find nearest upper left model grid point to each radar bin and store its 1D-index as well as the interpolation
      !  weights in rs_data(ista)-structure. The data are stored in vectors with points above the surface
      !  sorted to the beginning of the vectors and points below the surface thereafter. The vectors have a length
      !  equal to all radar points with locations inside the local processor domain only!
      ! The continuous 1D-index of radar points is stored in order to be able to retrieve azi, range, ele from
      !  this index:
      CALL setup_model2radarbins_vec (idom, rs_meta(ista), nobsmax, lon_r, lat_r, alt_r, &
           rs_data(ista)%nobs, rs_data(ista)%nobs_above_sfc, rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%hl_loc )

      ! As a result, the weights and indices and radar bin height are now stored here:
      !  rs_data(ista)%w_intp(1:rs_data(ista)%nobs,:)   = interpolation weights for model points relative to an upper left point
      !  rs_data(ista)%ind_intp(1:rs_data(ista)%nobs,1) = continuous 1D-index of upper left model point of the interpol. neighbourhood
      !  rs_data(ista)%ind_intp(1:rs_data(ista)%nobs,2) = continuous 1D-index of radar points (azi, range, ele)
      !  rs_data(ista)%hl_loc(1:rs_data(ista)%nobs)     = Height of radar point above MSL

      DEALLOCATE(lat_r)
      DEALLOCATE(lon_r)
      DEALLOCATE(alt_r)

      ! Compute local elevation and store in rs_data(ista)%el_loc:

      IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
      ALLOCATE(rs_data(ista)%el_loc(rs_data(ista)%nobs))

      !  a) points above the surface:
      CALL get_rangevec (rs_meta(ista), raarr)

      ALLOCATE(el(rs_meta(ista)%naz,rs_meta(ista)%nel))
      IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
        CALL get_elarr_precipscan (rs_meta(ista), el_precip)
        el(:,1) = el_precip
      ELSE
        DO iel=1, rs_meta(ista)%nel
          el(:,iel) = rs_meta(ista)%el_arr(iel)
        END DO
      END IF

!$omp parallel do private (m, n, o)
      DO irp = 1, rs_data(ista)%nobs_above_sfc
        ! Retrieve the indices of azi, range, ele from the continuous 1D-index stored in
        !  rs_data(ista)%ind_intp(:,2):
        CALL ind2sub3D( rs_data(ista)%ind_intp(irp,2), rs_meta(ista)%naz, rs_meta(ista)%nra, &
                        m, n, o )
        rs_data(ista)%el_loc(irp) = &
             el_loc_43 (raarr(n), el(m,o), rs_meta(ista)%alt_msl, r_earth_dp)
      END DO
!$omp end parallel do

      !  b) points below the surface get a certain missing value:
      DO irp = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
        rs_data(ista)%el_loc(irp) = shield_value
      END DO

      ! .. clean up:
      DEALLOCATE(raarr,el)
      IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_vec


  !==============================================================================
  !+ Module procedure in radar_src for the computation of radar subbeams geometry
  !------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_smth

    !------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar main- and subbeams geometry and the
    !              interpolation weights and interpolation indices for the
    !              interpolation to basis points of Gauss-Legendre Quadrature.
    !              These interpolation weights and the interpolation indices
    !              are constant in time since the radar geometry is also constant
    !              in time and this subroutine is only called in the first call of
    !              the operator.
    !
    !
    ! Method:
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER :: ista,iaz,ira,iel,irp,nobsmax,m,n,o,iv,ih

    REAL    (KIND=dp)    ::                   &
         el_sta              ! elevation angle of sub-beam at radar station

    REAL    (KIND=dp),      ALLOCATABLE    :: &
         raarr(:)         ,& ! array of radial bins
         el(:,:)          ,& ! array of elevations
         el_precip(:)     ,& ! array of elevations
         lat_smth(:)      ,& ! array of geographical latitudes for each radar point (nrp)
         lon_smth(:)      ,& ! array of geographical longitudes for each radar point (nrp)
         alt_smth(:)         ! array of altitudes above sea level for each radar point (nrp)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_smth
    !-----------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_smth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! allocate and initialize local aux arrays
      nobsmax = rs_meta(ista)%nra * rs_meta(ista)%naz * rs_meta(ista)%nel * rs_meta(ista)%ngpsm_v * rs_meta(ista)%ngpsm_h

      ALLOCATE(lat_smth(nobsmax))
      ALLOCATE(lon_smth(nobsmax))
      ALLOCATE(alt_smth(nobsmax))

!CDIR BEGIN COLLAPSE
      lat_smth   = 0.0_dp
      lon_smth   = 0.0_dp
      alt_smth   = -1.0_dp
!CDIR END


      ! calculate geographical coordinates (lat,lon,alt) of center- and auxiliary radar points
      CALL rad2geo_const_smth(rs_meta(ista),lat_smth,lon_smth,alt_smth)


      ! Find nearest upper left model grid point to each radar bin and store its 1D-index as well as the interpolation
      !  weights in rs_data(ista)-structure. The data are stored in vectors with points above the surface
      !  sorted to the beginning of the vectors and points below the surface thereafter. The vectors have a length
      !  equal to all radar points with locations inside the local processor domain only!
      ! The continuous 1D-index of radar points is stored in order to be able to retrieve azi, range, ele from
      !  this index:
      CALL setup_model2radarbins_vec (idom, rs_meta(ista), nobsmax, lon_smth, lat_smth, alt_smth, &
           rs_data(ista)%nsmth, rs_data(ista)%nsmth_above_sfc, rs_data(ista)%ind_intp_smth, &
           rs_data(ista)%w_intp_smth, rs_data(ista)%hl_loc )

      ! As a result, the weights and indices and radar bin height are now stored here:
      !  rs_data(ista)%w_intp_smth(1:rs_data(ista)%nobs,:)   = interpolation weights for model points relative to an upper left point
      !  rs_data(ista)%ind_intp_smth(1:rs_data(ista)%nobs,1) = continuous 1D-index of upper left model point of the interpol. neighbourhood
      !  rs_data(ista)%ind_intp_smth(1:rs_data(ista)%nobs,2) = continuous 1D-index of radar points (azi, range, ele)
      !  rs_data(ista)%hl_loc(1:rs_data(ista)%nobs)     = Height of radar point above MSL

      DEALLOCATE(lat_smth)
      DEALLOCATE(lon_smth)
      DEALLOCATE(alt_smth)

      ! Compute local elevation and store in rs_data(ista)%el_loc:

      IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
      ALLOCATE(rs_data(ista)%el_loc(rs_data(ista)%nsmth))

      !  a) points above the surface:
      CALL get_rangevec (rs_meta(ista), raarr)

      ALLOCATE(el(rs_meta(ista)%naz,rs_meta(ista)%nel))
      IF (TRIM(rs_meta(ista)%scanname) == 'PRECIP') THEN
        CALL get_elarr_precipscan (rs_meta(ista), el_precip)
        el(:,1) = el_precip
      ELSE
        DO iel=1, rs_meta(ista)%nel
          el(:,iel) = rs_meta(ista)%el_arr(iel)
        END DO
      END IF

!$omp parallel do private (m, n, o, iv, ih, el_sta)
      DO irp = 1, rs_data(ista)%nsmth_above_sfc
        ! Retrieve the indices of azi, range, ele from the continuous 1D-index stored in
        !  rs_data(ista)%ind_intp_smth(:,2):
        CALL ind2sub5D( rs_data(ista)%ind_intp_smth(irp,2), rs_meta(ista)%naz, rs_meta(ista)%nra, &
             rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
             m, n, o, iv, ih )

        el_sta = el(m,o) + smth_el_horzscan( rs_meta(ista)%Theta3, &
                                             rs_meta(ista)%smth_interv_fact, &
                                             rs_meta(ista)%xabscsm_v(iv) &
                                            )

        rs_data(ista)%el_loc(irp) = &
             el_loc_43 (raarr(n), el_sta, rs_meta(ista)%alt_msl, r_earth_dp)

      END DO
!$omp end parallel do

      !  b) points below the surface get a certain missing value:
      DO irp = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
        rs_data(ista)%el_loc(irp) = shield_value
      END DO

      ! .. clean up:
      DEALLOCATE(raarr, el)
      IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_smth


  !=======================================================================================
  !+ Module procedure in radar_src for the computation of the geometry for auxiliary grids
  !---------------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_grid

    !------------------------------------------------------------------------------
    !
    ! Description: Calculates the geometry for auxiliary grids and the interpolation weights
    !              and interpolation indices for the interpolation of the model
    !              data to the auxiliary grids. These interpolation
    !              weights and the interpolation indices are constant in time since
    !              the radar geometry is also constant in time and this subroutine
    !              is only called in the first call of the operator.
    !              A halo region of nbl_az at the azimutal field margins is
    !              taken into account (important for lsmooth=.true.).
    !
    ! Method:
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) :: yzroutine
    CHARACTER (LEN=80) :: yzerrmsg

    INTEGER            :: ista,i,j,k,iaz,ial

    REAL    (KIND=dp),            ALLOCATABLE    :: &
                                        azarr(:)         ,& ! array of azimuths
                                        alarr(:)         ,& ! array of arc length
                                        lat_g(:,:)       ,& ! array of geographical latitudes for each auxiliary grid (az,al)
                                        lon_g(:,:)          ! array of geographical longitudes for each auxiliary grid (az,al)

    !- End of header


    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_grid
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_grid'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ALLOCATE(lat_g(rs_grid(ista)%naz_nbl,rs_grid(ista)%nal+1))
      ALLOCATE(lon_g(rs_grid(ista)%naz_nbl,rs_grid(ista)%nal+1))
!CDIR BEGIN COLLAPSE
      lat_g    = 0.0_dp
      lon_g    = 0.0_dp
!CDIR END

      ! fill up azimuth array
      ALLOCATE(azarr(rs_grid(ista)%naz_nbl))
      CALL get_azvec_values ( rs_grid(ista)%naz_nbl, rs_meta(ista)%az_inc, &
                              rs_meta(ista)%az_start, nbl_az, azarr)

      ! fill up arc length array
      ALLOCATE(alarr(rs_grid(ista)%nal+1))
      DO ial = 1,rs_grid(ista)%nal+1
        alarr(ial) = (ial-1)*rs_grid(ista)%al_inc
      END DO

      ! calculate geographical coordinates (lat,lon) of grid points
      call grd2geo(rs_meta(ista)%lat,       &
           rs_meta(ista)%lon,       &
           rs_grid(ista)%naz_nbl,   &
           rs_grid(ista)%nal+1,     &
           azarr,                   &
           alarr,                   &
           lat_g,lon_g)

      ! Find nearest model grid points for interpolation and store their index and corresponding
      !  interpolation weights:
      CALL setup_model2azislices_vec (idom, rs_grid(ista), lon_g, lat_g)

      ! As a result, the weights and indices are now stored here:
      !  rs_grid(ista)%w_intp(1:rs_grid(ista)%ngrd,:)   = interpolation weights for model points relative to an upper left point
      !  rs_grid(ista)%ind_intp(1:rs_grid(ista)%ngrd,1) = continuous 1D-index of upper left model point of the interpol. neighbourhood
      !  rs_grid(ista)%ind_intp(1:rs_grid(ista)%ngrd,2) = continuous 1D-index of azi slice grid points (azi, arc distance, height)
      !  rs_grid(ista)%hl_grd(1:rs_grid(ista)%ngrd)     = height of aux grid point above MSL

      ! Deallocate auxiliary fields
      DEALLOCATE(azarr)
      DEALLOCATE(alarr)
      DEALLOCATE(lat_g)
      DEALLOCATE(lon_g)

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_grid


  !==============================================================================
  !+ Module procedure in radar_src for the computation of online radar beam geometry
  !------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_online(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar beam geometry and the interpolation weights
    !              and interpolation indices for the interpolation of the
    !              model data to the radar points.This routine will be called every time
    !              when the operator is called.
    !
    ! Method:      A method called TORE, described in
    !              Zeng et al. (2014): Radar beam tracing methods based on refractive index,
    !              Jour. Atmos. Ocean. Tech.
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                    :: ista,iaz,igrd,nae,nobsmax,nobs,m,n,k,&
         ira,iae,iel,irp,nrp,nrp2,naz,mmax,mmin,ko,no,ke_
    INTEGER                    :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL    (KIND=dp)          :: d1,d2,elsign, el_inc,  &
         rfridx11, rfridx21, rfridx12, rfridx22,              &
         wl,wk


    REAL    (KIND=dp), ALLOCATABLE    ::         &
         hl_azgrd(:,:,:)     ,&
         rfridx_azgrd(:,:,:) ,&
         hl(:,:)             ,&
         al(:,:)             ,&
         el(:,:)             ,&
         el_sta(:,:)         ,&
         el_precip(:)        ,&
         el_loc(:)           ,&
         rfridx(:,:)         ,&
         w_intp(:,:)         ,& ! array of horizontal interpolation weights for each observation (first dimension) and each spatial direction i,j,k (second dimension)
         w_intptmp(:,:)      ,& ! array of horizontal interpolation weights for each temporary radar point
         w_intptmp2(:,:)     ,&
         hl_tmp(:)           ,&
         hl_tmp2(:)          ,&
         hl_loc(:)           ,&
         s_loc(:)


    INTEGER, ALLOCATABLE ::         &
         flag(:,:)           ,&
         azarr_idx(:)        ,&
         ind_intptmp(:,:)    ,&
         ind_intptmp2(:,:)   ,&
         ind_intp(:,:)  ! array of grid indices for each observation (first dimension)
    ! the second dimension contains of
    ! 1:     the continuous number of the model grid cell associated with the observation
    ! 2: the observation indices in azimuthal, radial and elevation direction (irp)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1
    REAL(kind=dp), PARAMETER :: pi = 4.0d0 * atan(1.0d0)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_online
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_online'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (nbl_az /= 0) THEN
      CALL abort_run (my_radar_id, 20574, &
           'ERROR: problem in ' //TRIM(yzroutine)// &
           ': nbl_az > 0, which is not applicable here!', &
           'radar_model2rays.f90, '//TRIM(yzroutine))
    END IF

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      IF (.NOT.lcalc(ista)) CYCLE

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        naz     = iend(my_cart_id_fwo,ista) - istart(my_cart_id_fwo,ista) + 1

        ! fill up azimuth array
        ALLOCATE(azarr_idx(naz))
        DO iaz = 1,naz
          azarr_idx(iaz) = istart(my_cart_id_fwo,ista) + (iaz-1)
        END DO

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

        ! allocate and initialize local aux arrays

        nobsmax = rs_meta(ista)%nra * naz * rs_meta(ista)%nel
        nae     = naz * rs_meta(ista)%nel

        nobs = 0

        ALLOCATE(hl_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(rfridx_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(hl(    nae, 0:rs_meta(ista)%nra))  ! includes hl(0) = alt_msl to start beam propagation calculation at the antenna
        ALLOCATE(al(    nae, 0:rs_meta(ista)%nra))  ! includes al(0) = 0.0 to start beam propagation calculation at the antenna
        ALLOCATE(el(    nae,-2:rs_meta(ista)%nra))  ! includes el(-2:0) to fill with el_sta for the TORE elevation sign change criterion
        ALLOCATE(rfridx(nae, 0:rs_meta(ista)%nra))  ! includes rfridx(0) = rfridx_sta to start beam propagation calculation at the antenna
        ALLOCATE(w_intp(nobsmax,2))
        ALLOCATE(ind_intp(nobsmax,2))
        ALLOCATE(hl_loc(nobsmax))
        ALLOCATE(el_loc(nobsmax))
        ALLOCATE(s_loc(nobsmax))
        ALLOCATE(flag(nae,0:rs_meta(ista)%nra))

!CDIR BEGIN COLLAPSE
        hl_azgrd = miss_value
        rfridx_azgrd = miss_value
        hl       = miss_value
        al       = 0.0_dp
        el       = -99.0_dp
        el_loc   = -99.0_dp
        hl_loc   = miss_value
        s_loc    = miss_value
        rfridx   = miss_value
        w_intp   = -1.0_dp
        ind_intp = -1
        flag(:,0)= 1
        flag(:,1:rs_meta(ista)%nra) = 0
!CDIR END

!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          hl_azgrd(m,n,k) = rs_grid(ista)%hl_azgrd(igrd)

          rfridx_azgrd(m,n,k) = rs_grid(ista)%rfridx_azgrd(igrd)

        ENDDO


        ALLOCATE(w_intptmp(nae,2))
        ALLOCATE(ind_intptmp(nae,4))
        ALLOCATE(w_intptmp2(nae,2))
        ALLOCATE(ind_intptmp2(nae,4))
        ALLOCATE(hl_tmp(nae))
        ALLOCATE(hl_tmp2(nae))

        DO ira = 1, rs_meta(ista)%nra

!CDIR BEGIN COLLAPSE
          w_intptmp = -1.0_dp
          ind_intptmp = -1
          w_intptmp2 = -1.0_dp
          ind_intptmp2 = -1
          hl_tmp = miss_value
          hl_tmp2 = miss_value
          nrp = 0
          nrp2= 0
!CDIR END

          IF (ira == 1) THEN

            DO iae = 1, nae

              CALL ind2sub2D(iae, naz, iaz, iel)

              el(iae,ira-1)     = el_sta( azarr_idx(iaz), iel )
              el(iae,ira-2)     = el(iae,ira-1)
              el(iae,ira-3)     = el(iae,ira-1)
              rfridx(iae,ira-1) = rs_grid(ista)%rfridx_sta
              hl(iae,ira-1)     = rs_meta(ista)%alt_msl

            END DO

          END IF

          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0) THEN

              hl(iae,ira) = SQRT((r_earth_dp + hl(iae,ira-1))**2 +  rs_meta(ista)%ra_inc**2 - &
                   2*rs_meta(ista)%ra_inc*(r_earth_dp + hl(iae,ira-1)) * COS(el(iae,ira-1)*degrad + 0.5*pi)) - &
                   (r_earth_dp)

              d1 = SIN(el(iae,ira-1)*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))

              d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)

              al(iae,ira) = al(iae,ira-1) + r_earth_dp*ASIN(d1)

            ENDIF   ! flag(iae,ira-1) > 0

          ENDDO ! Loop over nae

          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0 .AND. &
                 FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 <= rs_grid(ista)%nal) THEN

              nrp = nrp + 1

              CALL ind2sub2D(iae, naz, iaz, iel)

              CALL sub2ind3D( azarr_idx(iaz), ira, iel, rs_meta(ista)%naz, rs_meta(ista)%nra, ind_intptmp(nrp,1))

              ! save indices for az, el, ra
              ind_intptmp(nrp,2) = iaz

              ind_intptmp(nrp,3) = FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 ! save indices of al

              w_intptmp(nrp,1) =  al(iae,ira)/rs_grid(ista)%al_inc - FLOOR(al(iae,ira)/rs_grid(ista)%al_inc)

              hl_tmp(nrp) =  hl(iae,ira)

            ENDIF

          ENDDO

          IF (nrp > 0) THEN

            CALL calc_vert_weight_online(nrp,hl_tmp(1:nrp), &
                 hl_azgrd,                                  &  ! hl_azgrd(istart_az:iend_az,:,:)
                 azarr_idx, naz,                       &  ! index of iaz in the aux. grid
                 istart(my_cart_id_fwo,ista), iend(my_cart_id_fwo,ista),      &
                 rs_grid(ista)%nal + 1,                     &
                 ke_,                                       &
                 ind_intptmp(1:nrp,2),                      &
                 ind_intptmp(1:nrp,3),                      &
                 w_intptmp(1:nrp,1),                        &
                 ind_intptmp(1:nrp,4),                      &
                 w_intptmp(1:nrp,2))

            ! This Loop needed for a better vecorization of the next loop
            DO irp = 1, nrp

              IF (ind_intptmp(irp,4) > 0) THEN

                nrp2 = nrp2 + 1

                ind_intptmp2(nrp2,1) = ind_intptmp(irp,1)
                ind_intptmp2(nrp2,2) = ind_intptmp(irp,2)
                ind_intptmp2(nrp2,3) = ind_intptmp(irp,3)
                ind_intptmp2(nrp2,4) = ind_intptmp(irp,4)

                w_intptmp2(nrp2,1)   = w_intptmp(irp,1)
                w_intptmp2(nrp2,2)   = w_intptmp(irp,2)

                hl_tmp2(nrp2)        = hl_tmp(irp)

              ENDIF

            ENDDO

!NEC$ ivdep
            DO irp = 1, nrp2

              nobs = nobs + 1

              w_intp(nobs,1) = w_intptmp2(irp,1)

              w_intp(nobs,2) = w_intptmp2(irp,2)

              ! save position of arounding grid points: iaz,ial,ihl
              CALL sub2ind3D(azarr_idx(ind_intptmp2(irp,2)), ind_intptmp2(irp,3), ind_intptmp2(irp,4), &
                   rs_meta(ista)%naz, rs_grid(ista)%nal+1, ind_intp(nobs,1))

              ! save position of radar point: iaz,ira,iel
              ind_intp(nobs,2) =  ind_intptmp2(irp,1)

              hl_loc(nobs) = hl_tmp2(irp)

              iaz = ind_intptmp2(irp,2)
              iel = (ind_intp(nobs,2)-1)/(rs_meta(ista)%naz*rs_meta(ista)%nra) + 1
              CALL sub2ind2D(iaz, iel, naz, iae)

              flag(iae,ira) = 1

              m = azarr_idx(ind_intptmp2(irp,2))
              n = ind_intptmp2(irp,3)
              no= MIN(n+1,rs_grid(ista)%nal+1)
              k = ind_intptmp2(irp,4)
              ko= one_level_down(k)

              rfridx11 = rfridx_azgrd(m,n,k)
              rfridx21 = rfridx_azgrd(m,no,k)
              rfridx12 = rfridx_azgrd(m,n,ko)
              rfridx22 = rfridx_azgrd(m,no,ko)

              wl = w_intp(nobs,1)
              wk = w_intp(nobs,2)


!!$ Missing value treatment at the 4 neighbour points:
!!$ 1) All 4 points are missing value => rfridx(iae,ira) = missing value
!!$ 2) In any other case:
!!$    Try to fill the missing value by surrounding points
!!$    to replace any missing value
!!$
              IF (rfridx11 < miss_threshold .AND. rfridx21 >= miss_threshold) THEN
                rfridx11 = rfridx21
              ELSEIF (rfridx11 >= miss_threshold .AND. rfridx21 < miss_threshold) THEN
                rfridx21 = rfridx11
              ELSEIF (rfridx11 < miss_threshold .AND. rfridx21 < miss_threshold) THEN
                wk = 1.0
              ELSE
                CONTINUE
              ENDIF


              IF (rfridx12 < miss_threshold .AND. rfridx22 >= miss_threshold) THEN
                rfridx12 = rfridx22
              ELSEIF (rfridx12 >= miss_threshold .AND. rfridx22 < miss_threshold) THEN
                rfridx22 = rfridx12
              ELSEIF (rfridx12 < miss_threshold .AND. rfridx22 < miss_threshold) THEN
                wk = 0.0
              ELSE
                CONTINUE
              ENDIF

              rfridx(iae,ira) = (rfridx11*(1.0_dp-wl)+rfridx21*wl)*(1.0_dp - wk) &
                              + (rfridx12*(1.0_dp-wl)+rfridx22*wl)*wk


              elsign = SIGN(1.0_dp,el(iae,ira-1))

              d2 = (r_earth_dp + hl(iae,ira-1))/ &
                   (r_earth_dp + hl(iae,ira))*rfridx(iae,ira-1)/rfridx(iae,ira)*COS(el(iae,ira-1)*degrad)

              IF ( d2 > 1.0 ) THEN
                ! .. This criterion triggers sign changes under ducting conditions, including
                !    a geometric correction for the reflection angle taking the earth curvature
                !    along a finite ray segment into account
                el_inc      =  rs_meta(ista)%ra_inc*COS(el(iae,ira-1))/(r_earth_dp + hl(iae,ira-1))
                el(iae,ira) = -(el(iae,ira-1) + el_inc)
              ELSEIF( el(iae,ira-1) < 0.0 .AND. (el(iae,ira-1) + (el(iae,ira-1) - el(iae,ira-3))) > 0.0 ) THEN
                ! .. This criterion triggers sign changes under normal conditions and negative elevations
                el(iae,ira) = -elsign*ACOS(MIN(MAX(d2,-1.0_dp),1.0_dp))/degrad
              ELSE
                ! .. Normal conditions and no sign change expected
                el(iae,ira) = elsign*ACOS(MIN(MAX(d2,-1.0_dp),1.0_dp))/degrad
              ENDIF

              el_loc(nobs) =  el(iae,ira)
              s_loc (nobs) =  al(iae,ira)

            ENDDO ! loop over nrp2

          END IF  ! nrp > 0

        ENDDO ! loop over nra


        DEALLOCATE(w_intptmp)
        DEALLOCATE(ind_intptmp)
        DEALLOCATE(hl_tmp)
        DEALLOCATE(w_intptmp2)
        DEALLOCATE(ind_intptmp2)
        DEALLOCATE(hl_tmp2)

        rs_data(ista)%nobs     = nobs

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp))   DEALLOCATE (rs_data(ista)%w_intp)
        IF (ASSOCIATED(rs_data(ista)%ind_intp)) DEALLOCATE (rs_data(ista)%ind_intp)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))   DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))    DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp(nobs,2))
        ALLOCATE(rs_data(ista)%ind_intp(nobs,2))
        ALLOCATE(rs_data(ista)%hl_loc(nobs))
        ALLOCATE(rs_data(ista)%el_loc(nobs))
        ALLOCATE(rs_data(ista)%s_loc(nobs))

        IF (nobs > 0) THEN

          ! copy values from auxiliary variables
          rs_data(ista)%w_intp   = w_intp(1:nobs,:)
          rs_data(ista)%ind_intp = ind_intp(1:nobs,:)
          rs_data(ista)%hl_loc   = hl_loc(1:nobs)
          rs_data(ista)%el_loc   = el_loc(1:nobs)
          rs_data(ista)%s_loc    = s_loc(1:nobs)

        END IF

        ! deallocate auxiliary fields
        DEALLOCATE(hl_azgrd)
        DEALLOCATE(rfridx_azgrd)
        DEALLOCATE(hl)
        DEALLOCATE(al)
        DEALLOCATE(el)
        IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)
        DEALLOCATE(el_sta)
        DEALLOCATE(rfridx)
        DEALLOCATE(w_intp)
        DEALLOCATE(ind_intp)
        DEALLOCATE(hl_loc)
        DEALLOCATE(el_loc)
        DEALLOCATE(s_loc)
        DEALLOCATE(flag)
        DEALLOCATE(azarr_idx)

      ELSE    ! iend(my_cart_id_fwo) == -1

        nobs = 0

        rs_data(ista)%nobs     = nobs

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp))   DEALLOCATE (rs_data(ista)%w_intp)
        IF (ASSOCIATED(rs_data(ista)%ind_intp)) DEALLOCATE (rs_data(ista)%ind_intp)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))   DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))    DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp(nobs,2))
        ALLOCATE(rs_data(ista)%ind_intp(nobs,2))
        ALLOCATE(rs_data(ista)%hl_loc(nobs))
        ALLOCATE(rs_data(ista)%el_loc(nobs))
        ALLOCATE(rs_data(ista)%s_loc(nobs))

      END IF   ! iend(my_cart_id_fwo) > -1

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_online


  SUBROUTINE calc_geometry_onlinenew(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar beam geometry and the interpolation weights
    !              and interpolation indices for the interpolation of the
    !              model data to the radar points.This routine will be called every time
    !              when the operator is called.
    !
    ! Method:      A method called SODE, described in
    !              Zeng et al. (2014): Radar beam tracing methods based on refractive index,
    !              Jour. Atmos. Ocean. Tech.
    !
    ! Input files:
    !
    ! Output files:
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER        :: ista,iaz,igrd,nstep,nae,nobsmax,nobs,m,n,k,&
                      ira,iae,iel,irp,nrp,nrp2,naz,mmax,mmin,ko,no,ke_
    INTEGER        :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL (KIND=dp) :: d1,d2,elsign,  &
                      rfridx11, rfridx21, rfridx12, rfridx22, &
                      wl,wk


    REAL    (KIND=dp), ALLOCATABLE    ::         &
         hl_azgrd(:,:,:)     ,&
         rfridx_azgrd(:,:,:) ,&
         hl(:,:)             ,&
         al(:,:)             ,&
         el(:,:)             ,&
         el_sta(:,:)         ,&
         el_precip(:)        ,&
         dhdr(:,:)           ,&
         el_loc(:)           ,&
         rfridx(:,:)         ,&
         w_intp(:,:)         ,& ! array of horizontal interpolation weights for each observation (first dimension) and each spatial direction i,j,k (second dimension)
         w_intptmp(:,:)      ,& ! array of horizontal interpolation weights for each temporary radar point
         w_intptmp2(:,:)     ,&
         hl_tmp(:)           ,&
         hl_tmp2(:)          ,&
         hl_loc(:)           ,&
         s_loc(:)            ,&
         alt_stavec(:)       ,&
         dhdr_stavec(:)      ,&
         rfridx_stavec(:)

    INTEGER, ALLOCATABLE ::   &
         flag(:,:)           ,&
         azarr_idx(:)        ,&
         ind_intptmp(:,:)    ,&
         ind_intptmp2(:,:)   ,&
         ind_intp(:,:)  ! array of grid indices for each observation (first dimension)
                        ! the second dimension contains of
                        ! 1:     the continuous number of the model grid cell associated with the observation
                        ! 2: the observation indices in azimuthal, radial and elevation direction (irp)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1
    REAL(kind=dp), PARAMETER :: pi = 4.0d0 * atan(1.0d0)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_onlinenew
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_onlinenew'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    IF (nbl_az /= 0) THEN
      CALL abort_run (my_radar_id, 20574, &
           'ERROR: problem in ' //TRIM(yzroutine)// &
           ': nbl_az > 0, which is not applicable here!', &
           'radar_model2rays.f90, '//TRIM(yzroutine))
    END IF

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      IF (.NOT.lcalc(ista)) CYCLE

      ! HERE NBL_AZ=0, SINCE NO HORIZONTAL INTERPOLATION WITHIN THE AZIMUT-SLICE-DATA!

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        naz     = iend(my_cart_id_fwo,ista) - istart(my_cart_id_fwo,ista) + 1

        ! fill up azimuth array
        ALLOCATE(azarr_idx(naz))
        DO iaz = 1,naz
          azarr_idx(iaz) = istart(my_cart_id_fwo,ista) + (iaz-1)
        END DO

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

        ! allocate and initialize local aux arrays
        nobsmax = rs_meta(ista)%nra * naz * rs_meta(ista)%nel
        nae     = naz * rs_meta(ista)%nel
        nstep = 1

        ALLOCATE(hl_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(rfridx_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(hl(nae,rs_meta(ista)%nra))
        ALLOCATE(al(nae,rs_meta(ista)%nra))
        ALLOCATE(el(nae,rs_meta(ista)%nra))
        ALLOCATE(dhdr(nae,rs_meta(ista)%nra))
        ALLOCATE(rfridx(nae,rs_meta(ista)%nra))
        ALLOCATE(w_intp(nobsmax,2))
        ALLOCATE(ind_intp(nobsmax,2))
        ALLOCATE(hl_loc(nobsmax))
        ALLOCATE(el_loc(nobsmax))
        ALLOCATE(s_loc(nobsmax))
        ALLOCATE(flag(nae,0:rs_meta(ista)%nra))
        ALLOCATE(alt_stavec(nae))
        ALLOCATE(dhdr_stavec(nae))
        ALLOCATE(rfridx_stavec(nae))

        hl_azgrd     = miss_value
        rfridx_azgrd = miss_value
        hl           = miss_value
        al           = 0.0_dp
        el           = miss_value
        dhdr         = miss_value
        el_loc       = miss_value
        s_loc        = miss_value
        rfridx       = miss_value
        w_intp       = -1.0_dp
        ind_intp     = -1
        flag(:,0)    = 1
        flag(:,1:rs_meta(ista)%nra) = 0

        alt_stavec        = rs_meta(ista)%alt_msl

        DO iae = 1, nae
          CALL ind2sub2D(iae, naz, iaz, iel)
          dhdr_stavec(iae) =  SIN( el_sta(azarr_idx(iaz),iel)*degrad )
        ENDDO

        rfridx_stavec     = rs_grid(ista)%rfridx_sta

        nobs = 0

        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          hl_azgrd(m,n,k) = rs_grid(ista)%hl_azgrd(igrd)

          rfridx_azgrd(m,n,k) = rs_grid(ista)%rfridx_azgrd(igrd)

        ENDDO

        ALLOCATE(w_intptmp(nae,2))
        ALLOCATE(ind_intptmp(nae,4))
        ALLOCATE(w_intptmp2(nae,2))
        ALLOCATE(ind_intptmp2(nae,4))
        ALLOCATE(hl_tmp(nae))
        ALLOCATE(hl_tmp2(nae))

        DO ira = 1, rs_meta(ista)%nra

          w_intptmp = -1.0_dp
          ind_intptmp = -1
          w_intptmp2 = -1.0_dp
          ind_intptmp2 = -1
          hl_tmp = miss_value
          hl_tmp2 = miss_value
          nrp = 0
          nrp2= 0

          IF (ira == 1) THEN

            CALL calc_height_rk4(nstep,nae,rs_meta(ista)%nra,alt_stavec,dhdr_stavec,rfridx_stavec, &
                 rs_meta(ista)%ra_inc,flag(:,0),hl(:,ira),dhdr(:,ira))

          ELSE

            CALL calc_height_rk4(nstep,nae,rs_meta(ista)%nra,hl(:,ira-1),dhdr(:,ira-1),rfridx(:,ira-1),&
                 rs_meta(ista)%ra_inc,flag(:,ira-1),hl(:,ira),dhdr(:,ira))

          ENDIF

          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0) THEN

              IF (ira == 1) THEN

                CALL ind2sub2D(iae, naz, iaz, iel)

                el(iae,ira) = ASIN(dhdr(iae,ira))/degrad
                d1 = SIN(el_sta(azarr_idx(iaz),iel)*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))
                d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)
                al(iae,ira) = r_earth_dp*ASIN(d1)

              ELSE

                el(iae,ira-1) = ASIN(dhdr(iae,ira-1))/degrad
                el(iae,ira) = ASIN(dhdr(iae,ira))/degrad
                d1 = SIN(el(iae,ira-1)*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))
                d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)
                al(iae,ira) = al(iae,ira-1) + r_earth_dp*ASIN(d1)

              ENDIF

            ENDIF   ! flag(iae,ira-1) > 0

          ENDDO ! Loop over nae

          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0 .AND. &
                 FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 <= rs_grid(ista)%nal) THEN

              nrp = nrp + 1

              CALL ind2sub2D(iae, naz, iaz, iel)

              ! save indices for az, el, ra
              CALL sub2ind3D( azarr_idx(iaz), ira, iel, rs_meta(ista)%naz, rs_meta(ista)%nra, ind_intptmp(nrp,1))

              ! save local indices of az
              ind_intptmp(nrp,2) = iaz

              ind_intptmp(nrp,3) = FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 ! save indices of al

              w_intptmp(nrp,1) =  al(iae,ira)/rs_grid(ista)%al_inc - FLOOR(al(iae,ira)/rs_grid(ista)%al_inc)

              hl_tmp(nrp) =  hl(iae,ira)

            ENDIF

          ENDDO

          IF (nrp > 0) THEN

            CALL calc_vert_weight_online(nrp,hl_tmp(1:nrp), &
                 hl_azgrd,                                  &  ! hl_azgrd(istart_az:iend_az,:,:)
                 azarr_idx, naz,                       &  ! index of iaz in the aux. grid
                 istart(my_cart_id_fwo,ista), iend(my_cart_id_fwo,ista),      &
                 rs_grid(ista)%nal + 1,                     &
                 ke_,                                       &
                 ind_intptmp(1:nrp,2),                      &
                 ind_intptmp(1:nrp,3),                      &
                 w_intptmp(1:nrp,1),                        &
                 ind_intptmp(1:nrp,4),                      &
                 w_intptmp(1:nrp,2))

            ! This Loop needed for a better vecorization of the next loop
            DO irp = 1, nrp

              IF (ind_intptmp(irp,4) > 0) THEN

                nrp2 = nrp2 + 1

                ind_intptmp2(nrp2,1) = ind_intptmp(irp,1)
                ind_intptmp2(nrp2,2) = ind_intptmp(irp,2)
                ind_intptmp2(nrp2,3) = ind_intptmp(irp,3)
                ind_intptmp2(nrp2,4) = ind_intptmp(irp,4)

                w_intptmp2(nrp2,1)   = w_intptmp(irp,1)
                w_intptmp2(nrp2,2)   = w_intptmp(irp,2)

                hl_tmp2(nrp2)        = hl_tmp(irp)

              ENDIF

            ENDDO


            IF (nrp2 > 0) THEN

              DO irp = 1, nrp2

                iaz = ind_intptmp2(irp,2)
                iel = (ind_intptmp2(irp,1)-1)/(rs_meta(ista)%naz*rs_meta(ista)%nra) + 1
                CALL sub2ind2D(iaz, iel, naz, iae)

                IF (flag(iae,ira-1) > 0) THEN

                  flag(iae,ira) = 1

                  nobs = nobs + 1

                  w_intp(nobs,1) = w_intptmp2(irp,1)

                  w_intp(nobs,2) = w_intptmp2(irp,2)

                  ! save position of arounding grid points: iaz,ial,ihl
                  CALL sub2ind3D(azarr_idx(ind_intptmp2(irp,2)), ind_intptmp2(irp,3), ind_intptmp2(irp,4), &
                       rs_meta(ista)%naz, rs_grid(ista)%nal+1, ind_intp(nobs,1))

                  ! save position of radar point: iaz,ira,iel
                  ind_intp(nobs,2) =  ind_intptmp2(irp,1)

                  hl_loc(nobs) = hl_tmp2(irp)

                  m = azarr_idx(ind_intptmp2(irp,2))
                  n = ind_intptmp2(irp,3)
                  no= MIN(n+1,rs_grid(ista)%nal+1)
                  k = ind_intptmp2(irp,4)
                  ko= one_level_down(k)

                  rfridx11 = rfridx_azgrd(m,n,k)
                  rfridx21 = rfridx_azgrd(m,no,k)
                  rfridx12 = rfridx_azgrd(m,n,ko)
                  rfridx22 = rfridx_azgrd(m,no,ko)

                  wl = w_intp(nobs,1)
                  wk = w_intp(nobs,2)


!!$ Missing value treatment at the 4 neighbour points:
!!$ 1) All 4 points are missing value => rfridx(iae,ira) = missing value
!!$ 2) In any other case:
!!$    Try to fill the missing value by surrounding points
!!$    to replace any missing value
!!$
                  IF (rfridx11 < miss_threshold .AND. rfridx21 >= miss_threshold) THEN
                    rfridx11 = rfridx21
                  ELSEIF (rfridx11 >= miss_threshold .AND. rfridx21 < miss_threshold) THEN
                    rfridx21 = rfridx11
                  ELSEIF (rfridx11 < miss_threshold .AND. rfridx21 < miss_threshold) THEN
                    wk = 1.0
                  ELSE
                    CONTINUE
                  ENDIF


                  IF (rfridx12 < miss_threshold .AND. rfridx22 >= miss_threshold) THEN
                    rfridx12 = rfridx22
                  ELSEIF (rfridx12 >= miss_threshold .AND. rfridx22 < miss_threshold) THEN
                    rfridx22 = rfridx12
                  ELSEIF (rfridx12 < miss_threshold .AND. rfridx22 < miss_threshold) THEN
                    wk = 0.0
                  ELSE
                    CONTINUE
                  ENDIF

                  rfridx(iae,ira) = (rfridx11*(1.0_dp-wl)+rfridx21*wl)*(1.0_dp - wk) &
                                  + (rfridx12*(1.0_dp-wl)+rfridx22*wl)*wk

                  el_loc(nobs) = el(iae,ira)
                  s_loc (nobs) = al(iae,ira)

                END IF

              ENDDO ! loop over nrp2

            ENDIF

          ENDIF

        ENDDO ! loop over nra

        DEALLOCATE(w_intptmp)
        DEALLOCATE(ind_intptmp)
        DEALLOCATE(hl_tmp)
        DEALLOCATE(w_intptmp2)
        DEALLOCATE(ind_intptmp2)
        DEALLOCATE(hl_tmp2)

        rs_data(ista)%nobs     = nobs

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp))   DEALLOCATE (rs_data(ista)%w_intp)
        IF (ASSOCIATED(rs_data(ista)%ind_intp)) DEALLOCATE (rs_data(ista)%ind_intp)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))   DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))    DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp(nobs,2))
        ALLOCATE(rs_data(ista)%ind_intp(nobs,2))
        ALLOCATE(rs_data(ista)%hl_loc(nobs))
        ALLOCATE(rs_data(ista)%el_loc(nobs))
        ALLOCATE(rs_data(ista)%s_loc(nobs))

        IF ( nobs > 0 ) THEN

          ! copy values from auxiliary variables
          rs_data(ista)%w_intp   = w_intp(1:nobs,:)
          rs_data(ista)%ind_intp = ind_intp(1:nobs,:)
          rs_data(ista)%hl_loc   = hl_loc(1:nobs)
          rs_data(ista)%el_loc   = el_loc(1:nobs)
          rs_data(ista)%s_loc    = s_loc(1:nobs)

        END IF

        ! deallocate auxiliary fields
        DEALLOCATE(hl_azgrd)
        DEALLOCATE(rfridx_azgrd)
        DEALLOCATE(hl)
        DEALLOCATE(al)
        DEALLOCATE(el)
        IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)
        DEALLOCATE(el_sta)
        DEALLOCATE(dhdr)
        DEALLOCATE(rfridx)
        DEALLOCATE(w_intp)
        DEALLOCATE(ind_intp)
        DEALLOCATE(hl_loc)
        DEALLOCATE(el_loc)
        DEALLOCATE(s_loc)
        DEALLOCATE(flag)
        DEALLOCATE(azarr_idx)
        DEALLOCATE(alt_stavec)
        DEALLOCATE(dhdr_stavec)
        DEALLOCATE (rfridx_stavec)

      ELSE     ! iend(my_cart_id_fwo) == -1

        nobs = 0

        rs_data(ista)%nobs     = nobs

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp))   DEALLOCATE (rs_data(ista)%w_intp)
        IF (ASSOCIATED(rs_data(ista)%ind_intp)) DEALLOCATE (rs_data(ista)%ind_intp)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))   DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))   DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))    DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp(nobs,2))
        ALLOCATE(rs_data(ista)%ind_intp(nobs,2))
        ALLOCATE(rs_data(ista)%hl_loc(nobs))
        ALLOCATE(rs_data(ista)%el_loc(nobs))
        ALLOCATE(rs_data(ista)%s_loc(nobs))

      END IF   !iend(my_cart_id_fwo) > -1

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_onlinenew

  !=======================================================================================
  !+ Module procedure in radar_src for the computation of radar main- and subbeam geometry
  !---------------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_onsmth (lcalc)

    !------------------------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar main- and subbeam geometry and the interpolation weights
    !              and interpolation indices for the interpolation of the auxiliary grids to the
    !              smoothing points.This routine will be called every time when the operator is called
    !
    ! Method:      A method called TORE, described in
    !              Zeng et al. (2014): Radar beam tracing methods based on refractive index,
    !              Jour. Atmos. Ocean. Tech.
    !
    ! Input files:
    !
    ! Output files:
    !
    !-------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !-------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER        :: ista,iaz,igrd,nae,nobsmax,nsmth,m,n,k,ira,iae,iel, &
                      irp,nrp,nrp2,naz,mmax,mmin,mo,ko,no,iv,ih,idummy,ke_
    INTEGER        :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL (KIND=dp) :: azstart,azend,d1,d2,elsign,el_inc,  &
         rfridx111,rfridx211,rfridx121,rfridx221,rfridx112,rfridx212,rfridx122,rfridx222, &
         wa,wl,wk


    REAL    (KIND=dp), ALLOCATABLE :: &
         hl_azgrd(:,:,:)     ,&
         rfridx_azgrd(:,:,:) ,&
         hl(:,:)             ,&
         al(:,:)             ,&
         el(:,:)             ,&
         el_sta(:,:)         ,&
         el_precip(:)        ,&
         el_loc(:)           ,&
         rfridx(:,:)         ,&
         w_intp(:,:)         ,& ! array of horizontal interpolation weights for each observation (first dimension) and each spatial direction i,j,k (second dimension)
         w_intptmp(:,:)      ,& ! array of horizontal interpolation weights for each temporary radar point
         w_intptmp2(:,:)     ,&
         hl_tmp(:)           ,&
         hl_tmp2(:)          ,&
         hl_loc(:)           ,&
         s_loc(:)            ,&
         azimuth_vec(:)

    INTEGER, ALLOCATABLE :: &
         flag(:,:)          ,&
         azarr_idx(:)       ,&
         ind_intptmp(:,:)   ,&
         ind_intptmp2(:,:)  ,&
         ind_intp(:,:)  ! array of grid indices for each observation (first dimension)
    ! the second dimension contains of
    ! 1:     the continuous number of the model grid cell associated with the observation
    ! 2: the observation indices in azimuthal, radial and elevation direction (irp)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1
    REAL(kind=dp), PARAMETER :: pi = 4.0d0 * atan(1.0d0)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_onsmth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_onsmth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      IF (.NOT.lcalc(ista)) CYCLE

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        naz       = iend(my_cart_id_fwo,ista) - istart(my_cart_id_fwo,ista) + 1 - 2*nbl_az

        ! fill up azimuth array
        ALLOCATE(azarr_idx(naz))
        DO iaz = 1,naz
          azarr_idx(iaz) = istart(my_cart_id_fwo,ista) + (iaz-1)
        END DO

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

        ! allocate and initialize local aux arrays
        nobsmax = rs_meta(ista)%nra * naz * rs_meta(ista)%nel * rs_meta(ista)%ngpsm_v * rs_meta(ista)%ngpsm_h
        nae     = naz * rs_meta(ista)%nel * rs_meta(ista)%ngpsm_v * rs_meta(ista)%ngpsm_h

!!$ UB>> Von "vollem" Azimut auf "Quasi-" Index (real) umgestellt, deswegen umformuliert:
        azstart = azarr_idx(1)   - nbl_az   ! die "+- nbl_az" sollen die Zwischenbereiche zwischen den PEs abdecken
        azend   = azarr_idx(naz) + nbl_az

        ALLOCATE(hl_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(rfridx_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(    hl(nae, 0:rs_meta(ista)%nra))  ! includes hl(0) = alt_msl to start beam propagation calculation at the antenna
        ALLOCATE(    al(nae, 0:rs_meta(ista)%nra))  ! includes al(0) = 0.0 to start beam propagation calculation at the antenna
        ALLOCATE(    el(nae,-2:rs_meta(ista)%nra))  ! includes el(-2:0) to fill with el_sta for the TORE elevation sign change criterion
        ALLOCATE(rfridx(nae, 0:rs_meta(ista)%nra))  ! includes rfridx(0) = rfridx_sta to start beam propagation calculation at the antenna
        ALLOCATE(w_intp(nobsmax,3))
        ALLOCATE(ind_intp(nobsmax,2))
        ALLOCATE(hl_loc(nobsmax))
        ALLOCATE(el_loc(nobsmax))
        ALLOCATE(s_loc(nobsmax))
        ALLOCATE(flag(nae,0:rs_meta(ista)%nra))

        hl_azgrd = miss_value
        rfridx_azgrd = miss_value
        hl       = miss_value
        al       = 0.0_dp
        el       = -99.0_dp
        el_loc   = -99.0_dp
        s_loc        = miss_value
        rfridx   = miss_value
        w_intp   = -1.0_dp
        ind_intp = -1
        flag(:,0)= 1
        flag(:,1:rs_meta(ista)%nra) = 0  ! see above

        ! .. Helper vector for the azimuth in degrees:
        ALLOCATE(azimuth_vec(nae))
        DO iae = 1, nae

          CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

          ! elevation of the beam center at the radar station:

!!$ inlined code from smth_az_horzscan(ista,ih,el_sta):
!!$        phi3_eff = rs_meta(ista)%alpha3_eff_0 + &
!!$             ( COS(el_sta*degrad) - 1.0 ) * rs_meta(ista)%dalpha * &
!!$             ( 1.0 - EXP(-1.5*rs_meta(ista)%dalpha/rs_meta(ista)%Phi3 ) )
!!$
!!$        smthpoint = 0.5 * phi3_eff * rs_meta(ista)%smth_interv_fact * &
!!$                   rs_meta(ista)%xabscsm_h(ih)
!!$
!!$        azimuth_vec(iae)   = azarr_idx(iaz) + smthpoint / rs_meta(ista)%az_inc


          ! azimut of the smoothing point, depending on el_sta:
          azimuth_vec(iae)   = azarr_idx(iaz) + &
               smth_az_horzscan( &
                                 rs_meta(ista)%alpha3_eff_0, &
                                 rs_meta(ista)%dalpha, &
                                 rs_meta(ista)%Phi3, &
                                 rs_meta(ista)%smth_interv_fact, &
                                 rs_meta(ista)%xabscsm_h(ih), &
                                 rs_meta(ista)%el_arr(iel) &
                               ) / rs_meta(ista)%az_inc
          ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.

        END DO

        nsmth = 0

        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          hl_azgrd(m,n,k) = rs_grid(ista)%hl_azgrd(igrd)

          rfridx_azgrd(m,n,k) = rs_grid(ista)%rfridx_azgrd(igrd)

        ENDDO


        ALLOCATE(w_intptmp(nae,3))
        ALLOCATE(ind_intptmp(nae,4))
        ALLOCATE(w_intptmp2(nae,3))
        ALLOCATE(ind_intptmp2(nae,4))
        ALLOCATE(hl_tmp(nae))
        ALLOCATE(hl_tmp2(nae))

        DO ira = 1, rs_meta(ista)%nra

          w_intptmp = -1.0_dp
          ind_intptmp = -1
          w_intptmp2 = -1.0_dp
          ind_intptmp2 = -1
          hl_tmp = miss_value
          hl_tmp2 = miss_value
          nrp = 0
          nrp2= 0

          IF (ira == 1) THEN

            DO iae = 1, nae

              CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

!!$ inlined code from smth_el_horzscan(ista,iv):
!!$              theta3_eff = rs_meta(ista)%Theta3
!!$
!!$              smthpoint = 0.5 * theta3_eff * rs_meta(ista)%smth_interv_fact * &
!!$                   rs_meta(ista)%xabscsm_v(iv)
!!$
!!$              el(iae,ira-1) = el_sta(azarr_idx(iaz),iel) + &
!!$                     smthpoint

              el(iae,ira-1)     = el_sta(azarr_idx(iaz),iel) + &
                                     smth_el_horzscan( rs_meta(ista)%Theta3, &
                                                       rs_meta(ista)%smth_interv_fact, &
                                                       rs_meta(ista)%xabscsm_v(iv)  &
                                                     )
              el(iae,ira-2)     = el(iae,ira-1)
              el(iae,ira-3)     = el(iae,ira-1)
              rfridx(iae,ira-1) = rs_grid(ista)%rfridx_sta
              hl(iae,ira-1)     = rs_meta(ista)%alt_msl

            END DO

          END IF

          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0) THEN

              hl(iae,ira) = SQRT((r_earth_dp + hl(iae,ira-1))**2 +  rs_meta(ista)%ra_inc**2       - &
                   2*rs_meta(ista)%ra_inc*(r_earth_dp+hl(iae,ira-1))*COS(el(iae,ira-1)*degrad + 0.5*pi)) - &
                   (r_earth_dp)

              d1 = SIN(el(iae,ira-1)*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))

              d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)

              al(iae,ira) = al(iae,ira-1) + (r_earth_dp)*ASIN(d1)

            ENDIF    ! flag(iae,ira-1) > 0

          ENDDO ! Loop over nae

          DO iae = 1, nae

            IF ( flag(iae,ira-1) > 0 .AND. &
                 FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 <= rs_grid(ista)%nal .AND. &
                 azstart <= azimuth_vec(iae) .AND. azimuth_vec(iae) <= azend ) THEN

              nrp = nrp + 1

              ! Compute iaz, iel, iv, ih from local hash index iae:
              CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

              ! Store 5D hash index, take azimut index without nbl_az:
              CALL sub2ind5D(azarr_idx(iaz), ira, iel, iv, ih, &
                   rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
                   ind_intptmp(nrp,1))

              ! save index of az of the left neighbouring
              ! azimut slice in the auxiliary grid, including nbl_az
              ind_intptmp(nrp,2) = FLOOR(azimuth_vec(iae)) + nbl_az

              ind_intptmp(nrp,3) = FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 ! save indices of al

              w_intptmp(nrp,1) = azimuth_vec(iae) - FLOOR(azimuth_vec(iae))

              w_intptmp(nrp,2) = al(iae,ira)/rs_grid(ista)%al_inc - FLOOR(al(iae,ira)/rs_grid(ista)%al_inc)

              hl_tmp(nrp) =  hl(iae,ira)

            ENDIF

          ENDDO

          IF (nrp > 0) THEN

            CALL calc_vert_weight_onsmth(nrp,hl_tmp(1:nrp), &
                 hl_azgrd,                                  &
                 istart(my_cart_id_fwo,ista),               &
                 iend(my_cart_id_fwo,ista),                 &
                 rs_grid(ista)%nal+1,                       &
                 ke_,                                       &
                 ind_intptmp(1:nrp,2),                      &
                 ind_intptmp(1:nrp,3),                      &
                 w_intptmp(1:nrp,1),                        &
                 w_intptmp(1:nrp,2),                        &
                 ind_intptmp(1:nrp,4),                      &
                 w_intptmp(1:nrp,3))


            ! This Loop needed for a better vecorization of the next loop
            DO irp = 1, nrp

              IF (ind_intptmp(irp,4) > 0) THEN

                nrp2 = nrp2 + 1

                ind_intptmp2(nrp2,1) = ind_intptmp(irp,1)
                ind_intptmp2(nrp2,2) = ind_intptmp(irp,2)
                ind_intptmp2(nrp2,3) = ind_intptmp(irp,3)
                ind_intptmp2(nrp2,4) = ind_intptmp(irp,4)

                w_intptmp2(nrp2,1)   = w_intptmp(irp,1)
                w_intptmp2(nrp2,2)   = w_intptmp(irp,2)
                w_intptmp2(nrp2,3)   = w_intptmp(irp,3)

                hl_tmp2(nrp2)        = hl_tmp(irp)

              ENDIF

            ENDDO

!NEC$ ivdep
            DO irp = 1, nrp2

              nsmth = nsmth + 1

              w_intp(nsmth,1) = w_intptmp2(irp,1)

              w_intp(nsmth,2) = w_intptmp2(irp,2)

              w_intp(nsmth,3) = w_intptmp2(irp,3)

              ! save postion of arounding grid points: iaz,ial,ihl
              CALL sub2ind3D( ind_intptmp2(irp,2), ind_intptmp2(irp,3), ind_intptmp2(irp,4), &
                   rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, ind_intp(nsmth,1))

              ! save position of radar point: iaz,ira,iel
              ind_intp(nsmth,2) =  ind_intptmp2(irp,1)

              hl_loc(nsmth) = hl_tmp2(irp)

              ! Helper indices to calculate iae and setting of flag:
              CALL ind2sub5D(ind_intp(nsmth,2), rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
                   iaz, idummy, iel, iv, ih)
              iaz = iaz - azarr_idx(1) + 1

              ! hash index iae:
              CALL sub2ind4D(iaz, iel, iv, ih, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iae)

              flag(iae,ira) = 1

              ! determine indices for interpolation in az, al, hl direction for current radar point
              ! (Indices relative to azimute-slice-grid, i.e., shifted by nbl_az)
              m = ind_intptmp2(irp,2)
              mo= MIN(m+1,iend(my_cart_id_fwo,ista))
              n = ind_intptmp2(irp,3)
              no= MIN(n+1,rs_grid(ista)%nal+1)
              k = ind_intptmp2(irp,4)
              ko= one_level_down(k)

              ! iae = MOD(ind_intp(nsmth,1)-1,rs_meta(ista)%naz*rs_meta(ista)%nel) + 1  ! campared with line 2694
              rfridx111 = rfridx_azgrd(m,n,k)
              rfridx211 = rfridx_azgrd(mo,n,k)
              rfridx121 = rfridx_azgrd(m,no,k)
              rfridx221 = rfridx_azgrd(mo,no,k)
              rfridx112 = rfridx_azgrd(m,n,ko)
              rfridx212 = rfridx_azgrd(mo,n,ko)
              rfridx122 = rfridx_azgrd(m,no,ko)
              rfridx222 = rfridx_azgrd(mo,no,ko)

              wa = w_intp(nsmth,1)
              wl = w_intp(nsmth,2)
              wk = w_intp(nsmth,3)

              IF (rfridx111 < miss_threshold .AND. rfridx211 >= miss_threshold) THEN
                rfridx111 = rfridx211
              ELSEIF (rfridx111 >= miss_threshold .AND. rfridx211 < miss_threshold) THEN
                rfridx211 = rfridx111
              ELSEIF (rfridx111 < miss_threshold .AND. rfridx211 < miss_threshold) THEN
                wl = 1.0
              ELSE
                CONTINUE
              ENDIF

              IF (rfridx121 < miss_threshold .AND. rfridx221 >= miss_threshold) THEN
                rfridx121 = rfridx221
              ELSEIF (rfridx121 >= miss_threshold .AND. rfridx221 < miss_threshold) THEN
                rfridx221 = rfridx121
              ELSEIF (rfridx121 < miss_threshold .AND. rfridx221 < miss_threshold) THEN
                wl = 0.0
              ELSE
                CONTINUE
              ENDIF

              IF (rfridx112 < miss_threshold .AND. rfridx212 >= miss_threshold) THEN
                rfridx112 = rfridx212
              ELSEIF (rfridx112 >= miss_threshold .AND. rfridx212 < miss_threshold) THEN
                rfridx212 = rfridx112
              ELSEIF (rfridx112 < miss_threshold .AND. rfridx212 < miss_threshold) THEN
                wl = 1.0
              ELSE
                CONTINUE
              ENDIF

              IF (rfridx122 < miss_threshold .AND. rfridx222 >= miss_threshold) THEN
                rfridx122 = rfridx222
              ELSEIF (rfridx122 >= miss_threshold .AND. rfridx222 < miss_threshold) THEN
                rfridx222 = rfridx122
              ELSEIF (rfridx122 < miss_threshold .AND. rfridx222 < miss_threshold) THEN
                wl = 0.0
              ELSE
                CONTINUE
              ENDIF

              rfridx(iae,ira) = ((rfridx111*(1.0_dp-wa)+rfridx211*wa)*(1.0_dp - wl)   + &
                   (rfridx121*(1.0_dp-wa)+rfridx221*wa)*wl)*(1.0_dp-wk) + &
                   ((rfridx112*(1.0_dp-wa)+rfridx212*wa)*(1.0_dp - wl)   + &
                   (rfridx122*(1.0_dp-wa)+rfridx222*wa)*wl)*wk

              elsign = SIGN(1.0_dp,el(iae,ira-1))

              d2 = (r_earth_dp + hl(iae,ira-1))/ &
                   (r_earth_dp + hl(iae,ira))*rfridx(iae,ira-1)/rfridx(iae,ira)*COS(el(iae,ira-1)*degrad)

              IF ( d2 > 1.0 ) THEN
                ! .. This criterion triggers sign changes under ducting conditions, including
                !    a geometric correction for the reflection angle taking the earth curvature
                !    along a finite ray segment into account
                el_inc      =  rs_meta(ista)%ra_inc*COS(el(iae,ira-1))/(r_earth_dp + hl(iae,ira-1))
                el(iae,ira) = -(el(iae,ira-1) + el_inc)
              ELSEIF( el(iae,ira-1) < 0.0 .AND. (el(iae,ira-1) + (el(iae,ira-1) - el(iae,ira-3))) > 0.0 ) THEN
                ! .. This criterion triggers sign changes under normal conditions and negative elevations
                el(iae,ira) = -elsign*ACOS(MIN(MAX(d2,-1.0_dp),1.0_dp))/degrad
              ELSE
                ! .. Normal conditions and no sign change expected
                el(iae,ira) = elsign*ACOS(MIN(MAX(d2,-1.0_dp),1.0_dp))/degrad
              ENDIF

              el_loc(nsmth) =  el(iae,ira)
              s_loc (nsmth) =  al(iae,ira)

            ENDDO ! loop over nrp2

          END IF  ! if nrp > 0


        ENDDO ! loop over nra

        DEALLOCATE(w_intptmp)
        DEALLOCATE(ind_intptmp)
        DEALLOCATE(hl_tmp)
        DEALLOCATE(w_intptmp2)
        DEALLOCATE(ind_intptmp2)
        DEALLOCATE(hl_tmp2)


        rs_data(ista)%nsmth     = nsmth

!!$        WRITE(*,*)"nsmth",nsmth
!!$        WRITE(*,*)"max hl",MAXVAL(hl)
!!$        WRITE(*,*)"min hl",MINVAL(hl)
!!$        WRITE(*,*)"max el",MAXVAL(el)
!!$        WRITE(*,*)"min el",MINVAL(el)

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp_smth))   DEALLOCATE (rs_data(ista)%w_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%ind_intp_smth)) DEALLOCATE (rs_data(ista)%ind_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))        DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))        DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))         DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp_smth(nsmth,3))
        ALLOCATE(rs_data(ista)%ind_intp_smth(nsmth,2))
        ALLOCATE(rs_data(ista)%hl_loc(nsmth))
        ALLOCATE(rs_data(ista)%el_loc(nsmth))
        ALLOCATE(rs_data(ista)%s_loc(nsmth))

        IF (nsmth > 0) THEN

          ! copy values from auxiliary variables
          rs_data(ista)%w_intp_smth   = w_intp(1:nsmth,:)
          rs_data(ista)%ind_intp_smth = ind_intp(1:nsmth,:)
          rs_data(ista)%hl_loc        = hl_loc(1:nsmth)
          rs_data(ista)%el_loc        = el_loc(1:nsmth)
          rs_data(ista)%s_loc         = s_loc(1:nsmth)

        END IF

        ! deallocate auxiliary fields
        DEALLOCATE(hl_azgrd)
        DEALLOCATE(rfridx_azgrd)
        DEALLOCATE(hl)
        DEALLOCATE(al)
        DEALLOCATE(el)
        IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)
        DEALLOCATE(el_sta)
        DEALLOCATE(w_intp)
        DEALLOCATE(ind_intp)
        DEALLOCATE(hl_loc)
        DEALLOCATE(el_loc)
        DEALLOCATE(s_loc)
        DEALLOCATE(flag)
        DEALLOCATE(rfridx)
        DEALLOCATE(azarr_idx)
        DEALLOCATE(azimuth_vec)

      ELSE     ! iend(my_cart_id_fwo) == -1

        nsmth = 0

        rs_data(ista)%nsmth     = nsmth

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp_smth))   DEALLOCATE (rs_data(ista)%w_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%ind_intp_smth)) DEALLOCATE (rs_data(ista)%ind_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))        DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))        DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))         DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp_smth(nsmth,3))
        ALLOCATE(rs_data(ista)%ind_intp_smth(nsmth,2))
        ALLOCATE(rs_data(ista)%hl_loc(nsmth))
        ALLOCATE(rs_data(ista)%el_loc(nsmth))
        ALLOCATE(rs_data(ista)%s_loc(nsmth))

      END IF   !iend(my_cart_id_fwo) > -1

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_onsmth


  !=======================================================================================
  !+ Module procedure in radar_src for the computation of radar main- and subbeam geometry
  !---------------------------------------------------------------------------------------

  SUBROUTINE calc_geometry_onsmthnew (lcalc)

    !------------------------------------------------------------------------------------------------
    !
    ! Description: Calculates the radar main- and subbeam geometry and the interpolation weights
    !              and interpolation indices for the interpolation of the auxiliary grids to the
    !              smoothing points.This routine will be called every time when the operator is called
    !
    ! Method:      A method called SODE, described in
    !              Zeng et al. (2014): Radar beam tracing methods based on refractive index,
    !              Jour. Atmos. Ocean. Tech.
    ! Method:
    !
    ! Input files:
    !
    ! Output files:
    !
    !-------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !-------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER        :: ista,iaz,igrd,nstep,nae,nobsmax,nsmth,m,n,k,ira,iae,iel, &
                      irp,nrp,nrp2,naz,mmax,mmin,mo,ko,no,iv,ih,idummy,ke_
    INTEGER        :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL (KIND=dp) :: azstart,azend,el_tmp,d1,d2,elsign,  &
         rfridx111,rfridx211,rfridx121,rfridx221,rfridx112,rfridx212,rfridx122,rfridx222, &
         wa,wl,wk


    REAL    (KIND=dp), ALLOCATABLE :: &
         hl_azgrd(:,:,:)     ,&
         rfridx_azgrd(:,:,:) ,&
         hl(:,:)             ,&
         al(:,:)             ,&
         el(:,:)             ,&
         el_sta(:,:)         ,&
         el_precip(:)        ,&
         dhdr(:,:)           ,&
         el_loc(:)           ,&
         rfridx(:,:)         ,&
         w_intp(:,:)         ,& ! array of horizontal interpolation weights for each observation (first dimension) and each spatial direction i,j,k (second dimension)
         w_intptmp(:,:)      ,& ! array of horizontal interpolation weights for each temporary radar point
         w_intptmp2(:,:)     ,&
         hl_tmp(:)           ,&
         hl_tmp2(:)          ,&
         hl_loc(:)           ,&
         s_loc(:)            ,&
         alt_stavec(:)       ,&
         dhdr_stavec(:)      ,&
         rfridx_stavec(:)    ,&
         azimuth_vec(:)


    INTEGER, ALLOCATABLE :: &
         flag(:,:)          ,&
         azarr_idx(:)       ,&
         ind_intptmp(:,:)   ,&
         ind_intptmp2(:,:)  ,&
         ind_intp(:,:)  ! array of grid indices for each observation (first dimension)
                        ! the second dimension contains of
                        ! 1:     the continuous number of the model grid cell associated with the observation
                        ! 2: the observation indices in azimuthal, radial and elevation direction (irp)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1
    REAL(kind=dp), PARAMETER :: pi = 4.0d0 * atan(1.0d0)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_geometry_onsmthnew
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_geometry_onsmthnew'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      IF (.NOT.lcalc(ista)) CYCLE

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        naz       = iend(my_cart_id_fwo,ista) - istart(my_cart_id_fwo,ista) + 1 - 2*nbl_az

        ! fill up azimuth array
        ALLOCATE(azarr_idx(naz))
        DO iaz = 1,naz
          azarr_idx(iaz) = istart(my_cart_id_fwo,ista) + (iaz-1)
        END DO

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

        ! allocate and initialize local aux arrays
        nobsmax = rs_meta(ista)%nra * naz * rs_meta(ista)%nel * rs_meta(ista)%ngpsm_v * rs_meta(ista)%ngpsm_h
        nae     = naz * rs_meta(ista)%nel * rs_meta(ista)%ngpsm_v * rs_meta(ista)%ngpsm_h
        nstep   = 1

!!$ UB>> Von "vollem" Azimut auf "Quasi-" Index (real) umgestellt, deswegen umformuliert:
        azstart = azarr_idx(1)   - nbl_az   ! die "+- nbl_az" sollen die Zwischenbereiche zwischen den PEs abdecken
        azend   = azarr_idx(naz) + nbl_az

        ALLOCATE(hl_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(rfridx_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(hl(nae,rs_meta(ista)%nra))
        ALLOCATE(al(nae,rs_meta(ista)%nra))
        ALLOCATE(el(nae,rs_meta(ista)%nra))
        ALLOCATE(dhdr(nae,rs_meta(ista)%nra))
        ALLOCATE(rfridx(nae,rs_meta(ista)%nra))
        ALLOCATE(w_intp(nobsmax,3))
        ALLOCATE(ind_intp(nobsmax,2))
        ALLOCATE(hl_loc(nobsmax))
        ALLOCATE(el_loc(nobsmax))
        ALLOCATE(s_loc(nobsmax))
        ALLOCATE(flag(nae,0:rs_meta(ista)%nra))
        ALLOCATE(alt_stavec(nae))
        ALLOCATE(dhdr_stavec(nae))
        ALLOCATE(rfridx_stavec(nae))

        hl_azgrd = miss_value
        rfridx_azgrd = miss_value
        hl       = miss_value
        al       = 0.0_dp
        el       = miss_value
        dhdr     = miss_value
        el_loc   = miss_value
        rfridx   = miss_value
        w_intp   = -1.0_dp
        ind_intp = -1
        flag(:,0)= 1
        flag(:,1:rs_meta(ista)%nra) = 0  ! see above

        alt_stavec        = rs_meta(ista)%alt_msl

        ! .. Helper vector for the azimuth in degrees:
        ALLOCATE(azimuth_vec(nae))
        DO iae = 1, nae

          ! Compute iaz, iel, iv, ih from local hash index iae:
          CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

!!$ code inlined from smth_el_horzscan:
!!$        theta3_eff = rs_meta(ista)%Theta3
!!$
!!$        smthpoint = 0.5 * theta3_eff * rs_meta(ista)%smth_interv_fact * &
!!$             rs_meta(ista)%xabscsm_v(iv)
!!$
!!$        el_tmp = el_sta(azarr_idx(iaz),iel) + smthpoint

          ! elevation of the beam center at the radar station:
          el_tmp = el_sta(azarr_idx(iaz),iel) + &
                       smth_el_horzscan( rs_meta(ista)%Theta3, &
                                         rs_meta(ista)%smth_interv_fact, &
                                         rs_meta(ista)%xabscsm_v(iv)  &
                                       )
          dhdr_stavec(iae) = SIN(el_tmp*degrad)

!!$ code inlined from smth_az_horzscan;
!!$        phi3_eff = rs_meta(ista)%alpha3_eff_0 + &
!!$             ( COS(el_sta*degrad) - 1.0 ) * rs_meta(ista)%dalpha * &
!!$             ( 1.0 - EXP(-1.5*rs_meta(ista)%dalpha/rs_meta(ista)%Phi3 ) )
!!$
!!$        smthpoint = 0.5 * phi3_eff * rs_meta(ista)%smth_interv_fact * &
!!$             rs_meta(ista)%xabscsm_h(ih)
!!$
!!$        azimuth_vec(iae)   = azarr_idx(iaz) + smthpoint / rs_meta(ista)%az_inc

          ! azimut of the smoothing point, depending on el_sta:
          azimuth_vec(iae)   = azarr_idx(iaz) + &
               smth_az_horzscan( &
                                 rs_meta(ista)%alpha3_eff_0, &
                                 rs_meta(ista)%dalpha, &
                                 rs_meta(ista)%Phi3, &
                                 rs_meta(ista)%smth_interv_fact, &
                                 rs_meta(ista)%xabscsm_h(ih), &
                                 rs_meta(ista)%el_arr(iel) &
                               ) / rs_meta(ista)%az_inc
          ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.

        END DO


        rfridx_stavec     = rs_grid(ista)%rfridx_sta

        nsmth = 0

        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          hl_azgrd(m,n,k) = rs_grid(ista)%hl_azgrd(igrd)

          rfridx_azgrd(m,n,k) = rs_grid(ista)%rfridx_azgrd(igrd)

        ENDDO


        ALLOCATE(w_intptmp(nae,3))
        ALLOCATE(ind_intptmp(nae,4))
        ALLOCATE(w_intptmp2(nae,3))
        ALLOCATE(ind_intptmp2(nae,4))
        ALLOCATE(hl_tmp(nae))
        ALLOCATE(hl_tmp2(nae))

        DO ira = 1, rs_meta(ista)%nra


          w_intptmp = -1.0_dp
          ind_intptmp = -1
          w_intptmp2 = -1.0_dp
          ind_intptmp2 = -1
          hl_tmp = miss_value
          hl_tmp2 = miss_value
          nrp = 0
          nrp2= 0

          IF (ira == 1) THEN

            CALL calc_height_rk4(nstep,nae,rs_meta(ista)%nra,alt_stavec,dhdr_stavec,rfridx_stavec,&
                 rs_meta(ista)%ra_inc,flag(:,0),hl(:,ira),dhdr(:,ira))

          ELSE

            CALL calc_height_rk4(nstep,nae,rs_meta(ista)%nra,hl(:,ira-1),dhdr(:,ira-1),rfridx(:,ira-1),&
                 rs_meta(ista)%ra_inc,flag(:,ira-1),hl(:,ira),dhdr(:,ira))

          ENDIF


          DO iae = 1, nae

            IF (flag(iae,ira-1) > 0) THEN

              IF (ira == 1) THEN

                CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

                el_tmp = el_sta(azarr_idx(iaz),iel) + &
                     smth_el_horzscan( rs_meta(ista)%Theta3, &
                                       rs_meta(ista)%smth_interv_fact, &
                                       rs_meta(ista)%xabscsm_v(iv))
                el(iae,ira) = ASIN(dhdr(iae,ira))/degrad
                d1 = SIN(el_tmp*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))
                d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)
                al(iae,ira) = (r_earth_dp)*ASIN(d1)

              ELSE

                el(iae,ira-1) = ASIN(dhdr(iae,ira-1))/degrad
                el(iae,ira) = ASIN(dhdr(iae,ira))/degrad
                d1 = SIN(el(iae,ira-1)*degrad + 0.5*pi) * rs_meta(ista)%ra_inc/(r_earth_dp + hl(iae,ira))
                d1 = MIN(MAX(d1,-1.0_dp),1.0_dp)
                al(iae,ira) = al(iae,ira-1) + (r_earth_dp)*ASIN(d1)

              ENDIF

            ENDIF    ! flag(iae,ira-1) > 0

          ENDDO ! Loop over nae

          DO iae = 1, nae

            IF ( flag(iae,ira-1) > 0 .AND. &
                 FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 <= rs_grid(ista)%nal .AND. &
                 REAL(azstart,dp) <= azimuth_vec(iae) .AND. azimuth_vec(iae) <= REAL(azend,dp)) THEN

              nrp = nrp + 1

              ! Compute iaz, iel, iv, ih from local hash index iae:
              CALL ind2sub4D(iae, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iaz, iel, iv, ih)

              ! Store 5D hash index, take azimut index without nbl_az:
              CALL sub2ind5D(azarr_idx(iaz), ira, iel, iv, ih, &
                   rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
                   ind_intptmp(nrp,1))

              ! save index of az of the left neighbouring
              ! azimut slice in the auxiliary grid, including nbl_az
              ind_intptmp(nrp,2) = FLOOR(azimuth_vec(iae)) + nbl_az

              ind_intptmp(nrp,3) = FLOOR(al(iae,ira)/rs_grid(ista)%al_inc) + 1 ! save indices of al

              w_intptmp(nrp,1) = azimuth_vec(iae) - FLOOR(azimuth_vec(iae))

              w_intptmp(nrp,2) = al(iae,ira)/rs_grid(ista)%al_inc - FLOOR(al(iae,ira)/rs_grid(ista)%al_inc)

              hl_tmp(nrp) =  hl(iae,ira)

!!$              IF (hl_tmp(nrp) < 0.0_dp) WRITE(*,*)'hl_tmp is wrong'

            ENDIF

          ENDDO

          IF (nrp > 0) THEN

            CALL calc_vert_weight_onsmth(nrp,hl_tmp(1:nrp), &
                 hl_azgrd,                                  &
                 istart(my_cart_id_fwo,ista),                   &
                 iend(my_cart_id_fwo,ista),                     &
                 rs_grid(ista)%nal+1,                       &
                 ke_,                                       &
                 ind_intptmp(1:nrp,2),                      &
                 ind_intptmp(1:nrp,3),                      &
                 w_intptmp(1:nrp,1),                        &
                 w_intptmp(1:nrp,2),                        &
                 ind_intptmp(1:nrp,4),                      &
                 w_intptmp(1:nrp,3))

            ! This Loop needed for a better vecorization of the next loop
            DO irp = 1, nrp

              IF (ind_intptmp(irp,4) > 0) THEN

                nrp2 = nrp2 + 1

                ind_intptmp2(nrp2,1) = ind_intptmp(irp,1)
                ind_intptmp2(nrp2,2) = ind_intptmp(irp,2)
                ind_intptmp2(nrp2,3) = ind_intptmp(irp,3)
                ind_intptmp2(nrp2,4) = ind_intptmp(irp,4)

                w_intptmp2(nrp2,1)   = w_intptmp(irp,1)
                w_intptmp2(nrp2,2)   = w_intptmp(irp,2)
                w_intptmp2(nrp2,3)   = w_intptmp(irp,3)

                hl_tmp2(nrp2)        = hl_tmp(irp)

              ENDIF

            ENDDO


            IF (nrp2 > 0) THEN

              DO irp = 1, nrp2

                ! Helper indices to calculate iae and setting of flag:
                CALL ind2sub5D(ind_intptmp2(irp,1), rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, &
                     rs_meta(ista)%ngpsm_v, iaz, idummy, iel, iv, ih)
                iaz = iaz  - azarr_idx(1) + 1

                ! Hash index iae:
                CALL sub2ind4D(iaz, iel, iv, ih, naz, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, iae)

                IF (flag(iae,ira-1) > 0) THEN

                  flag(iae,ira) = 1

                  nsmth = nsmth + 1

                  w_intp(nsmth,1) = w_intptmp2(irp,1)

                  w_intp(nsmth,2) = w_intptmp2(irp,2)

                  w_intp(nsmth,3) = w_intptmp2(irp,3)

                  ! save postion of arounding grid points: iaz,ial,ihl
                  CALL sub2ind3D( ind_intptmp2(irp,2), ind_intptmp2(irp,3), ind_intptmp2(irp,4), &
                       rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, ind_intp(nsmth,1))

                  ! save position of radar point: iaz,ira,iel
                  ind_intp(nsmth,2) =  ind_intptmp2(irp,1)

                  hl_loc(nsmth) = hl_tmp2(irp)

                  ! determine indices for interpolation in az, al, hl direction for current radar point
                  ! (Indices relative to azimute-slice-grid, i.e., shifted by nbl_az)
                  m = ind_intptmp2(irp,2)
                  mo= MIN(m+1,iend(my_cart_id_fwo,ista))
                  n = ind_intptmp2(irp,3)
                  no= MIN(n+1,rs_grid(ista)%nal+1)
                  k = ind_intptmp2(irp,4)
                  ko= one_level_down(k)


                  rfridx111 = rfridx_azgrd(m,n,k)
                  rfridx211 = rfridx_azgrd(mo,n,k)
                  rfridx121 = rfridx_azgrd(m,no,k)
                  rfridx221 = rfridx_azgrd(mo,no,k)
                  rfridx112 = rfridx_azgrd(m,n,ko)
                  rfridx212 = rfridx_azgrd(mo,n,ko)
                  rfridx122 = rfridx_azgrd(m,no,ko)
                  rfridx222 = rfridx_azgrd(mo,no,ko)

                  wa = w_intp(nsmth,1)
                  wl = w_intp(nsmth,2)
                  wk = w_intp(nsmth,3)

                  IF (rfridx111 < miss_threshold .AND. rfridx211 >= miss_threshold) THEN
                    rfridx111 = rfridx211
                  ELSEIF (rfridx111 >= miss_threshold .AND. rfridx211 < miss_threshold) THEN
                    rfridx211 = rfridx111
                  ELSEIF (rfridx111 < miss_threshold .AND. rfridx211 < miss_threshold) THEN
                    wl = 1.0
                  ELSE
                    CONTINUE
                  ENDIF

                  IF (rfridx121 < miss_threshold .AND. rfridx221 >= miss_threshold) THEN
                    rfridx121 = rfridx221
                  ELSEIF (rfridx121 >= miss_threshold .AND. rfridx221 < miss_threshold) THEN
                    rfridx221 = rfridx121
                  ELSEIF (rfridx121 < miss_threshold .AND. rfridx221 < miss_threshold) THEN
                    wl = 0.0
                  ELSE
                    CONTINUE
                  ENDIF

                  IF (rfridx112 < miss_threshold .AND. rfridx212 >= miss_threshold) THEN
                    rfridx112 = rfridx212
                  ELSEIF (rfridx112 >= miss_threshold .AND. rfridx212 < miss_threshold) THEN
                    rfridx212 = rfridx112
                  ELSEIF (rfridx112 < miss_threshold .AND. rfridx212 < miss_threshold) THEN
                    wl = 1.0
                  ELSE
                    CONTINUE
                  ENDIF

                  IF (rfridx122 < miss_threshold .AND. rfridx222 >= miss_threshold) THEN
                    rfridx122 = rfridx222
                  ELSEIF (rfridx122 >= miss_threshold .AND. rfridx222 < miss_threshold) THEN
                    rfridx222 = rfridx122
                  ELSEIF (rfridx122 < miss_threshold .AND. rfridx222 < miss_threshold) THEN
                    wl = 0.0
                  ELSE
                    CONTINUE
                  ENDIF

                  rfridx(iae,ira) = ((rfridx111*(1.0_dp-wa)+rfridx211*wa)*(1.0_dp - wl)   + &
                                     (rfridx121*(1.0_dp-wa)+rfridx221*wa)*wl)*(1.0_dp-wk) + &
                                    ((rfridx112*(1.0_dp-wa)+rfridx212*wa)*(1.0_dp - wl)   + &
                                     (rfridx122*(1.0_dp-wa)+rfridx222*wa)*wl)*wk


                  el_loc(nsmth) =  el(iae,ira)
                  s_loc (nsmth) =  al(iae,ira)

                END IF

              ENDDO ! loop over nrp2

            END IF  ! if nrp2 > 0

          END IF  ! if nrp > 0

        ENDDO ! loop over nra

        DEALLOCATE(w_intptmp)
        DEALLOCATE(ind_intptmp)
        DEALLOCATE(hl_tmp)
        DEALLOCATE(w_intptmp2)
        DEALLOCATE(ind_intptmp2)
        DEALLOCATE(hl_tmp2)


        rs_data(ista)%nsmth     = nsmth

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp_smth))   DEALLOCATE (rs_data(ista)%w_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%ind_intp_smth)) DEALLOCATE (rs_data(ista)%ind_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))        DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))        DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))         DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp_smth(nsmth,3))
        ALLOCATE(rs_data(ista)%ind_intp_smth(nsmth,2))
        ALLOCATE(rs_data(ista)%hl_loc(nsmth))
        ALLOCATE(rs_data(ista)%el_loc(nsmth))
        ALLOCATE(rs_data(ista)%s_loc(nsmth))

        IF (nsmth > 0) THEN

          ! copy values from auxiliary variables
          rs_data(ista)%w_intp_smth   = w_intp(1:nsmth,:)
          rs_data(ista)%ind_intp_smth = ind_intp(1:nsmth,:)
          rs_data(ista)%hl_loc        = hl_loc(1:nsmth)
          rs_data(ista)%el_loc        = el_loc(1:nsmth)
          rs_data(ista)%s_loc         = s_loc(1:nsmth)

        END IF

        ! deallocate auxiliary fields
        DEALLOCATE(hl_azgrd)
        DEALLOCATE(rfridx_azgrd)
        DEALLOCATE(hl)
        DEALLOCATE(al)
        DEALLOCATE(el)
        IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)
        DEALLOCATE(el_sta)
        DEALLOCATE(dhdr)
        DEALLOCATE(w_intp)
        DEALLOCATE(ind_intp)
        DEALLOCATE(hl_loc)
        DEALLOCATE(el_loc)
        DEALLOCATE(s_loc)
        DEALLOCATE(flag)
        DEALLOCATE(rfridx)
        DEALLOCATE(azarr_idx)
        DEALLOCATE(azimuth_vec)
        DEALLOCATE(alt_stavec)
        DEALLOCATE(dhdr_stavec)
        DEALLOCATE(rfridx_stavec)

      ELSE     ! iend(my_cart_id_fwo) == -1

        nsmth = 0

        rs_data(ista)%nsmth     = nsmth

        ! allocate arrays of interpolation coefficients and indices of current radar station
        !  (for safety: deallocate first, if associated)
        IF (ASSOCIATED(rs_data(ista)%w_intp_smth))   DEALLOCATE (rs_data(ista)%w_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%ind_intp_smth)) DEALLOCATE (rs_data(ista)%ind_intp_smth)
        IF (ASSOCIATED(rs_data(ista)%hl_loc))        DEALLOCATE (rs_data(ista)%hl_loc)
        IF (ASSOCIATED(rs_data(ista)%el_loc))        DEALLOCATE (rs_data(ista)%el_loc)
        IF (ASSOCIATED(rs_data(ista)%s_loc))         DEALLOCATE (rs_data(ista)%s_loc)

        ALLOCATE(rs_data(ista)%w_intp_smth(nsmth,3))
        ALLOCATE(rs_data(ista)%ind_intp_smth(nsmth,2))
        ALLOCATE(rs_data(ista)%hl_loc(nsmth))
        ALLOCATE(rs_data(ista)%el_loc(nsmth))
        ALLOCATE(rs_data(ista)%s_loc(nsmth))

      END IF   !iend(my_cart_id_fwo) > -1

    END DO ! loop over stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_geometry_onsmthnew

  !======================================================================================
  !+ Module procedure in radar_src for the computation of the model hydrometeor fallspeed
  !--------------------------------------------------------------------------------------

  SUBROUTINE calc_mod_fallspeed(lcalc)
   !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the fallspeed for each
    !              radar station.
    !
    ! Method:
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista, iobs

    !- End of header
    !==============================================================================


    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_fallspeed
    !-----------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_mod_fallspeed'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate fallspeed on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'
      vt_radar = 0.0_dp
      CALL calc_fallspeed_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), lweightdbz, ldebug_radsim, vt_radar)
#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of model reflectivities:
      ALLOCATE(rs_data(ista)%vt_mod(rs_data(ista)%nobs))

      ! Set data points below surface to shield_value:
      DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
        rs_data(ista)%vt_mod(iobs) = shield_value
      END DO

      CALL interp_model2radarbins_scalar ( vt_radar, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%vt_mod )

      IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
        WRITE(*,*)'max rs_data(ista)%vt_mod = ', MAXVAL(rs_data(ista)%vt_mod), '  nobs = ', rs_data(ista)%nobs
        WRITE(*,*)'min rs_data(ista)%vt_mod = ', MINVAL(rs_data(ista)%vt_mod), '  nobs = ', rs_data(ista)%nobs
      END IF

#ifdef __COSMO__
     CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_fallspeed

  !======================================================================================
  !+ Module procedure in radar_src for the computation of the model hydrometeor fallspeed
  !--------------------------------------------------------------------------------------

  SUBROUTINE calc_mod_fallspeed_smth(lcalc)

    !--------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the fall speed at smoothing points for each
    !              radar station.
    !
    ! Method:
    !
    !--------------------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg


    INTEGER         :: ista,i,j,k,m,n,o,offset_i,offset_j,izerror,&
                       nk,np,ii,jj,kk,iu,io,ju,jo,ku,ko,iv,ih,ismth
    REAL (KIND=dp)  :: wi_smth,                        & !
                       wj_smth,                        & !
                       wk_smth                           !
    REAL (KIND=dp)  ::                                 &
         vt111, vt211, vt121, vt221,     &
         vt112, vt212, vt122, vt222

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_fallspeed_smth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine = 'calc_mod_fallspeed_smth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id


    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate reflectivity on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'
      vt_radar = 0.0_dp
      CALL calc_fallspeed_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), lweightdbz, ldebug_radsim, vt_radar)
#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of model reflectivities:
      ALLOCATE(rs_data(ista)%vt_mod_smth(rs_data(ista)%nsmth))

      ! Set data points below surface to shield_value (back side of the vectors):
      DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
        rs_data(ista)%vt_mod_smth(ismth) = shield_value
      END DO

      CALL interp_model2radarbins_scalar ( vt_radar, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%vt_mod_smth )

      IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
        WRITE(*,*)'max rs_data(ista)%vt_mod_smth = ', MAXVAL(rs_data(ista)%vt_mod_smth), '  nsmth = ', rs_data(ista)%nsmth
        WRITE(*,*)'min rs_data(ista)%vt_mod_smth = ', MINVAL(rs_data(ista)%vt_mod_smth), '  nsmth = ', rs_data(ista)%nsmth
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with calc_mod_vt_smth() on proc ', my_radar_id

  END SUBROUTINE calc_mod_fallspeed_smth


  !===================================================================================
  !+ Module procedure in radar_src for the computation of the model fall speed
  !  under consideration of dynamical propagation path
  !-----------------------------------------------------------------------------------
  SUBROUTINE calc_mod_fallspeed_online(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model fall speed for each
    !              radar station
    !
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                :: ista,igrd,nobsmax,nobs,i,j,k,m,n,o,offset_i,offset_j,iobs,nk,np,no,ko,ke_
    INTEGER                :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL (KIND=dp)         ::                    &
                              wl                ,& ! interpolation weight in i-direction
                              wk                ,& ! interpolation weight in j-direction
                              vt11, vt21, vt12, vt22
    REAL    (KIND=dp), ALLOCATABLE      ::  vt_azgrd(:,:,:)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_fallspeed_online
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_fallspeed_online'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model terminal fallspeed:
      ALLOCATE(rs_data(ista)%vt_mod(rs_data(ista)%nobs))

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(vt_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
!!$ UB>> Initialize with standard missing value (= point outside model domain):
        vt_azgrd = miss_value

        !sort the collected smoothing points data into 3d field according to sorted azimuth, arc length, height.
!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          vt_azgrd(m,n,k) = rs_grid(ista)%vt_azgrd(igrd)

        ENDDO

        rs_data(ista)%vt_mod = miss_value

        ! loop over observation points
!CDIR NODEP,VOVERTAKE,VOB
        DO iobs = 1, rs_data(ista)%nobs

          ! for each radar point determine the continuous number nk of the grid cube
          ! surrounded by the 8 nearest model grid points
          nk = rs_data(ista)%ind_intp(iobs,1)

          ! for each radar point determine the continuous number np of radar points
          np = rs_data(ista)%ind_intp(iobs,2)

          ! determine indices for interpolation grids in atimuthal, arc length, vertical direction for current radar point
          CALL ind2sub3D(nk, rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, o)

          ! determine interpolation weights
          wl = rs_data(ista)%w_intp(iobs,1)
          wk = rs_data(ista)%w_intp(iobs,2)

          no= MIN(n+1,rs_grid(ista)%nal+1)
          ko= one_level_down(k)

          vt11 = vt_azgrd(m,n,k)
          vt21 = vt_azgrd(m,no,k)
          vt12 = vt_azgrd(m,n,ko)
          vt22 = vt_azgrd(m,no,ko)

          ! If the value on a grid point is missing, it will be replaced by value from the horizontal symmetric grid point;
          ! If both value are missing, neglect these two points while intepolation
          IF (vt11 < miss_threshold .AND. vt21 >= miss_threshold)     THEN
            vt11 = vt21
          ELSEIF (vt11 >= miss_threshold .AND. vt21 < miss_threshold) THEN
            vt21 = vt11
          ELSEIF (vt11 < miss_threshold .AND. vt21 < miss_threshold) THEN
            wk = 1.0
          ELSE
            CONTINUE
          ENDIF


          IF (vt12 < miss_threshold .AND. vt22 >= miss_threshold)     THEN
            vt12 = vt22
          ELSEIF (vt12 >= miss_threshold .AND. vt22 < miss_threshold) THEN
            vt22 = vt12
          ELSEIF (vt12 < miss_threshold .AND. vt22 < miss_threshold) THEN
            wk = 0.0
          ELSE
            CONTINUE
          ENDIF

          rs_data(ista)%vt_mod(iobs) = (vt11*(1.0_dp-wl)+vt21*wl)*(1.0_dp - wk) + (vt12*(1.0_dp-wl)+vt22*wl)*wk

        END DO    ! loop over radar points

        IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
          WRITE(*,*) 'max rs_data(ista)%vt_mod = ', &
                     MAXVAL(rs_data(ista)%vt_mod), '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%vt_mod = ', &
                     MINVAL(rs_data(ista)%vt_mod), '  nobs = ', rs_data(ista)%nobs
        END IF

        DEALLOCATE(vt_azgrd)

      END IF   ! iend(my_cart_id_fwo) > -1

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_fallspeed_online


  !=======================================================================================================
  !+ Module procedure in radar_src for the computation of the model fall speed at smoothing points
  !  under consideration of dynamical propagation path
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mod_fallspeed_onsmth(lcalc)

    !---------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model fall speed  for each
    !              radar station
    !
    !---------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)

    INTEGER                :: ista,igrd,nobsmax,ismth,i,j,k,m,n,o,offset_i,offset_j,iobs,nk,mo,no,ko,ke_
    REAL (KIND=dp)         ::      &
         wa  ,& ! interpolation weight in i-direction
         wl  ,& ! interpolation weight in j-direction
         wk  ,&
         vt111, vt211, vt121, vt221, vt112, vt212, vt122, vt222

    REAL    (KIND=dp), ALLOCATABLE      :: vt_azgrd(:,:,:)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_fallspeed_onsmth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_fallspeed_onsmth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model fall speeds:
      ALLOCATE(rs_data(ista)%vt_mod_smth(rs_data(ista)%nsmth))

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(vt_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ! UB>> Initialize with standard missing value (= point outside model domain):
        vt_azgrd = miss_value

        !sort the collected smoothing points data into 3d filed according to sorted azimuth, arc length, heigh
!NEC$ ivdep
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          vt_azgrd(m,n,k) = rs_grid(ista)%vt_azgrd(igrd)

        ENDDO

        rs_data(ista)%vt_mod_smth = miss_value

        ! loop over observation points
!NEC$ ivdep
        DO ismth = 1, rs_data(ista)%nsmth

          ! for each radar point determine the continuous number nk of the grid cube
          ! surrounded by the 8 nearest model grid points
          nk = rs_data(ista)%ind_intp_smth(ismth,1)

          CALL ind2sub3D(nk, rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, o)

          mo= MIN(m+1,iend(my_cart_id_fwo,ista))
          no= MIN(n+1,rs_grid(ista)%nal+1)
          ko= one_level_down(k)

          ! If the value on a grid point is missing (=miss_value) (i.e., is outside the model domain),
          ! it will be replaced by VALUE from the horizontal symmetric grid point;
          ! If both values are missing, neglect these two points while interpolation
          vt111 = vt_azgrd(m,n,k)
          vt211 = vt_azgrd(mo,n,k)
          vt121 = vt_azgrd(m,no,k)
          vt221 = vt_azgrd(mo,no,k)
          vt112 = vt_azgrd(m,n,ko)
          vt212 = vt_azgrd(mo,n,ko)
          vt122 = vt_azgrd(m,no,ko)
          vt222 = vt_azgrd(mo,no,ko)

          wa = rs_data(ista)%w_intp_smth(ismth,1)
          wl = rs_data(ista)%w_intp_smth(ismth,2)
          wk = rs_data(ista)%w_intp_smth(ismth,3)

          IF (vt111 < miss_threshold .AND. vt211 >= miss_threshold)     THEN
            vt111 = vt211
          ELSEIF (vt111 >= miss_threshold .AND. vt211 < miss_threshold) THEN
            vt211 = vt111
          ELSEIF (vt111 < miss_threshold .AND. vt211 < miss_threshold) THEN
            wl = 1.0
          ELSE
            CONTINUE
          ENDIF

          IF (vt121 < miss_threshold .AND. vt221 >= miss_threshold)     THEN
            vt121 = vt221
          ELSEIF (vt121 >= miss_threshold .AND. vt221 < miss_threshold) THEN
            vt221 = vt121
          ELSEIF (vt121 < miss_threshold .AND. vt221 < miss_threshold) THEN
            wl = 0.0
          ELSE
            CONTINUE
          ENDIF

          IF (vt112 < miss_threshold .AND. vt212 >= miss_threshold)     THEN
            vt112 = vt212
          ELSEIF (vt112 >= miss_threshold .AND. vt212 < miss_threshold) THEN
            vt212 = vt112
          ELSEIF (vt112 < miss_threshold .AND. vt212 < miss_threshold) THEN
            wl = 1.0
          ELSE
            CONTINUE
          ENDIF

          IF (vt122 < miss_threshold .AND. vt222 >= miss_threshold)     THEN
            vt122 = vt222
          ELSEIF (vt122 >= miss_threshold .AND. vt222 < miss_threshold) THEN
            vt222 = vt122
          ELSEIF (vt122 < miss_threshold .AND. vt222 < miss_threshold) THEN
            wl = 0.0
          ELSE
            CONTINUE
          ENDIF

          IF (vt111 < miss_threshold .OR. vt211 <  miss_threshold .OR. &
               vt121 < miss_threshold .OR. vt221 <  miss_threshold .OR. &
               vt112 < miss_threshold .OR. vt212 <  miss_threshold .OR. &
               vt122 < miss_threshold .OR. vt222 <  miss_threshold) THEN

            rs_data(ista)%vt_mod_smth(ismth) =   miss_value

          ELSE

            rs_data(ista)%vt_mod_smth(ismth) = ((vt111*(1.0_dp-wa)+vt211*wa)*(1.0_dp - wl)   + &
                                                (vt121*(1.0_dp-wa)+vt221*wa)*wl)*(1.0_dp-wk) + &
                                               ((vt112*(1.0_dp-wa)+vt212*wa)*(1.0_dp - wl)   + &
                                                (vt122*(1.0_dp-wa)+vt222*wa)*wl)*wk

          END IF

        END DO    ! loop over radar points

!!$        IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
!!$          WRITE(*,*) 'max rs_data(ista)%vt_mod_smth = ', &
!!$                     MAXVAL(rs_data(ista)%vt_mod_smth), '  nsmth = ', rs_data(ista)%nsmth
!!$          WRITE(*,*) 'min rs_data(ista)%vt_mod_smth = ', &
!!$                     MINVAL(rs_data(ista)%vt_mod_smth), '  nsmth = ', rs_data(ista)%nsmth
!!$        END IF

        DEALLOCATE(vt_azgrd)

      END IF

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_fallspeed_onsmth



  !================================================================================================
  !+ Module procedure in radar_src for the computation of the radar reflectivity at auxiliary grids
  !------------------------------------------------------------------------------------------------
  SUBROUTINE calc_grd_fallspeed(lcalc)

    !-------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the fall speed at auxiliary grids for each
    !              radar station. Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method:
    !
    !-------------------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_grd_fallspeed
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_grd_fallspeed'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate reflectivity on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'

      vt_radar = 0.0_dp
      CALL calc_fallspeed_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), lweightdbz, ldebug_radsim, vt_radar)
#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_fallspeed_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of reflectivities at auxiliary grids:
      ALLOCATE(rs_grid(ista)%vt_grd(rs_grid(ista)%ngrd))

      CALL interp_model2azislices_scalar (REAL(vt_radar, kind=wp), rs_grid(ista)%ngrd, &
           rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%vt_grd)

      IF (ldebug_radsim .AND. rs_grid(ista)%ngrd > 0) THEN
        WRITE(*,*) 'max rs_grid(ista)%vt_grd = ', &
                   MAXVAL(rs_grid(ista)%vt_grd), ' ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%vt_grd = ', &
                   MINVAL(rs_grid(ista)%vt_grd), ' ngrd = ', rs_grid(ista)%ngrd
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_grd_fallspeed


  !==============================================================================
  !+ Module procedure in radar_src for the computation of the model radial winds
  !------------------------------------------------------------------------------

  SUBROUTINE calc_mod_radialwind(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radial winds for each
    !              radar station by
    !              - bilinear interpolation of u,v and w to radar points
    !              - projection onto radar beam
    !
    ! Method:      Following Jaervinen et al., 2009, Tellus
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:

    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                            :: ista, iobs, np, m, n, o
    REAL    (KIND=dp)                  :: range, azimuth, alpha
    REAL    (KIND=dp), ALLOCATABLE     ::  u_rp(:) ,v_rp(:) ,w_rp(:)   ! work arrays for the vel. components

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_radialwind
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_radialwind'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model radial winds
      ALLOCATE(rs_data(ista)%radwind_mod(rs_data(ista)%nobs))

      ! Set data points below surface to shield_value (back side of the vectors):
      DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
        rs_data(ista)%radwind_mod(iobs) = shield_value
      END DO

      ALLOCATE(u_rp(rs_data(ista)%nobs_above_sfc))
      ALLOCATE(v_rp(rs_data(ista)%nobs_above_sfc))
      ALLOCATE(w_rp(rs_data(ista)%nobs_above_sfc))

      CALL interp_model2radarbins_vr (u, v, w, &
           rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, &
           u_rp, v_rp, w_rp)

!$omp parallel do private(np, m, n, o, range, azimuth, alpha)
      DO iobs = 1, rs_data(ista)%nobs_above_sfc

        ! for each radar point determine the continuous number np of radar points
        np = rs_data(ista)%ind_intp(iobs,2)

        ! determine indices in azimuthal radial and elevation direction for current radar point
        CALL ind2sub3D(np, rs_meta(ista)%naz, rs_meta(ista)%nra, m, n ,o)

        ! calculate range, azimuth and elevation of radar point
        range     = n*rs_meta(ista)%ra_inc
        azimuth   = (m-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start

        ! Formerly: correction angle due to curvature of earth (see Jaervinen et al., 2009).
        ! Now: has been replaced by rs_data(ista)%el_loc, the exact local elevation angle:
        alpha =  rs_data(ista)%el_loc(iobs)*degrad

        rs_data(ista)%radwind_mod(iobs) = (u_rp(iobs)*SIN(azimuth*degrad) + v_rp(iobs)*COS(azimuth*degrad)) * &
                                           COS(alpha) + w_rp(iobs)*SIN(alpha)

      END DO    ! loop over radar points
!$omp end parallel do

      DEALLOCATE(u_rp, v_rp, w_rp)

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_radialwind


  !============================================================================================
  !+ Module procedure in radar_src for the computation of the model radial winds with smoothing
  !--------------------------------------------------------------------------------------------

  SUBROUTINE calc_mod_radialwind_smth(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radial winds for each
    !              radar station by
    !              - trilinear interpolation of u,v and w to smoothing points
    !              - projection onto radar main- or subbeam
    !              - smoothing with Gauss-Legendre Quadraure
    ! Method:      Following Jaervinen et al., 2009, Tellus
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                            :: ista,m,n,o,np,iv,ih,ismth
    REAL    (KIND=dp)                  :: range,azimuth,alpha
    REAL    (KIND=dp), ALLOCATABLE     ::  u_rp(:) ,v_rp(:) ,w_rp(:)   ! work arrays for the vel. components

    !- End of header
    !==============================================================================


    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_radialwind_smth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine = 'calc_mod_radialwind_smth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ALLOCATE(rs_data(ista)%radwind_mod_smth(rs_data(ista)%nsmth))

      ! Set data points below surface to shield_value (back side of the vectors):
      DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
        rs_data(ista)%radwind_mod_smth(ismth) = shield_value
      END DO

      ALLOCATE(u_rp(rs_data(ista)%nsmth_above_sfc))
      ALLOCATE(v_rp(rs_data(ista)%nsmth_above_sfc))
      ALLOCATE(w_rp(rs_data(ista)%nsmth_above_sfc))

      CALL interp_model2radarbins_vr (u, v, w, &
           rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, &
           u_rp, v_rp, w_rp)

!$omp parallel do private(np, m, n, o, iv, ih, range, azimuth, alpha)
      DO ismth = 1, rs_data(ista)%nsmth_above_sfc

        ! for each radar point determine the continuous number np of radar points
        np = rs_data(ista)%ind_intp_smth(ismth,2)

        ! determine m, n, o, iv, ih:
        CALL ind2sub5D(np, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, m, n, o, iv, ih)

        ! calculate range, azimuth and elevation of smoothing point
        range     = n*rs_meta(ista)%ra_inc

!!$ inlined code from smth_az_horzscan:
!!$        phi3_eff = rs_meta(ista)%alpha3_eff_0 + &
!!$             ( COS(rs_meta(ista)%el_arr(o)*degrad) - 1.0 ) * rs_meta(ista)%dalpha * &
!!$             ( 1.0 - EXP(-1.5*rs_meta(ista)%dalpha/rs_meta(ista)%Phi3 ) )
!!$
!!$        smthpoint = 0.5 * phi3_eff * rs_meta(ista)%smth_interv_fact * &
!!$                   rs_meta(ista)%xabscsm_h(ih)
!!$
!!$        azimuth   = (m-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start + smthpoint

        azimuth   = (m-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start + &
                    smth_az_horzscan( &
                                      rs_meta(ista)%alpha3_eff_0, &
                                      rs_meta(ista)%dalpha, &
                                      rs_meta(ista)%Phi3, &
                                      rs_meta(ista)%smth_interv_fact, &
                                      rs_meta(ista)%xabscsm_h(ih), &
                                      rs_meta(ista)%el_arr(o) &
                                    )
        ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.


        ! Formerly: correction angle due to curvature of earth (see Jaervinen et al., 2009).
        ! Now: has been replaced by rs_data(ista)%el_loc, the exact local elevation angle:
        alpha =  rs_data(ista)%el_loc(ismth)*degrad


        ! calculate model radial wind at smoothing point
        rs_data(ista)%radwind_mod_smth(ismth) = (u_rp(ismth)*SIN(azimuth*degrad) + v_rp(ismth)*COS(azimuth*degrad)) * &
                                                                      COS(alpha) + w_rp(ismth)*SIN(alpha)


      END DO    ! loop over radar smoothing points
!$omp end parallel do

      DEALLOCATE(u_rp, v_rp, w_rp)

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_radialwind_smth


  !==============================================================================
  !+ Module procedure in radar_src for the computation of the model radial winds
  !  under consideration of dynamical propagation path
  !------------------------------------------------------------------------------

  SUBROUTINE calc_mod_radialwind_online(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radial winds for each
    !              radar station by
    !              - bilinear interpolation of u,v and w to radar points
    !              - projection onto radar beam
    !
    ! Method:      Following Jaervinen et al., 2009, Tellus
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                        :: ista,igrd,iaz,ira,iel,nobsmax,nobs,i,j,k,m,n,o,offset_i,offset_j,iobs,nk,np,no,ko,ke_
    INTEGER                        :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL    (KIND=dp)                   ::                     &
         wl                 ,& ! interpolation weight in i-direction
         wk                 ,& ! interpolation weight in j-direction
         u11, u21, u12, u22 ,&
         v11, v21, v12, v22 ,&
         w11, w21, w12, w22 ,&
         u,v,w,azimuth,alpha

    REAL    (KIND=dp), ALLOCATABLE      :: u_azgrd(:,:,:),v_azgrd(:,:,:),w_azgrd(:,:,:)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_radialwind_online
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_mod_radialwind_online'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model radial winds:
      ALLOCATE(rs_data(ista)%radwind_mod(rs_data(ista)%nobs))

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(u_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(v_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(w_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))

!!$ UB>> Initialize with standard missing value (= point outside model domain):
        u_azgrd = miss_value
        v_azgrd = miss_value
        w_azgrd = miss_value

        !sort the collected smoothing points data into 3d field according to sorted azimuth, arc length, height.
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          u_azgrd(m,n,k) = rs_grid(ista)%u_azgrd(igrd)

          v_azgrd(m,n,k) = rs_grid(ista)%v_azgrd(igrd)

          w_azgrd(m,n,k) = rs_grid(ista)%w_azgrd(igrd)

        ENDDO

        rs_data(ista)%radwind_mod = miss_value

        ! loop over observation points
!CDIR NODEP,VOVERTAKE,VOB
        DO iobs = 1, rs_data(ista)%nobs

          ! for each radar point determine the continuous number nk of the grid cube
          ! surrounded by the 4 nearest auxiliary grid points
          nk = rs_data(ista)%ind_intp(iobs,1)

          ! for each radar point determine the continuous number np of radar points
          np = rs_data(ista)%ind_intp(iobs,2)

          ! determine indices in azimuthal radial and elevation direction for current radar point
          CALL ind2sub3D(np, rs_meta(ista)%naz, rs_meta(ista)%nra, iaz, ira, iel)

          ! determine indices for interpolation grids in atimuthal, arc length, vertical direction for current radar point
          CALL ind2sub3D(nk, rs_meta(ista)%naz, rs_grid(ista)%nal+1, m, n, k)
          m = m + nbl_az

          ! determine interpolation weights
          wl = rs_data(ista)%w_intp(iobs,1)
          wk = rs_data(ista)%w_intp(iobs,2)

          no= MIN(n+1,rs_grid(ista)%nal+1)
          ko= one_level_down(k)

          u11 = u_azgrd(m,n,k)
          u21 = u_azgrd(m,no,k)
          u12 = u_azgrd(m,n,ko)
          u22 = u_azgrd(m,no,ko)

          v11 = v_azgrd(m,n,k)
          v21 = v_azgrd(m,no,k)
          v12 = v_azgrd(m,n,ko)
          v22 = v_azgrd(m,no,ko)

          w11 = w_azgrd(m,n,k)
          w21 = w_azgrd(m,no,k)
          w12 = w_azgrd(m,n,ko)
          w22 = w_azgrd(m,no,ko)

          ! If the value on a grid point is missing, it will be replaced by value from the horizontal symmetric grid point;
          ! If both value are missing, neglect these two points while intepolation
          IF (u11 < miss_threshold .AND. u21 >= miss_threshold) THEN
            u11 = u21
            v11 = v21
            w11 = w21
          ELSEIF (u11 >= miss_threshold .AND. u21 < miss_threshold) THEN
            u21 = u11
            v21 = v11
            w21 = w11
          ELSEIF (u11 < miss_threshold .AND. u21 < miss_threshold) THEN
            wk = 1.0
          ELSE
            CONTINUE
          ENDIF

          IF (u12 < miss_threshold .AND. u22 >= miss_threshold) THEN
            u12 = u22
            v12 = v22
            w12 = w22
          ELSEIF (u12 >= miss_threshold .AND. u22 < miss_threshold) THEN
            u22 = u12
            v22 = v12
            w22 = w12
          ELSEIF (u12 < miss_threshold .AND. u22 < miss_threshold) THEN
            wk = 0.0
          ELSE
            CONTINUE
          ENDIF

          u = (u11*(1.0_dp-wl)+u21*wl)*(1.0_dp - wk) + (u12*(1.0_dp-wl)+u22*wl)*wk

          v = (v11*(1.0_dp-wl)+v21*wl)*(1.0_dp - wk) + (v12*(1.0_dp-wl)+v22*wl)*wk

          w = (w11*(1.0_dp-wl)+w21*wl)*(1.0_dp - wk) + (w12*(1.0_dp-wl)+w22*wl)*wk


          ! UB<<==========================================================================

          ! calculate azimuth of radar point
          azimuth   = (iaz-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start
          azimuth   = MODULO(azimuth,360.0_dp)

          ! Formerly: correction angle due to curvature of earth (see Jaervinen et al., 2009).
          ! Now: has been replaced by the true local elevation angle:
          alpha = rs_data(ista)%el_loc(iobs)*degrad

          ! calculate model radial wind at radar point
          rs_data(ista)%radwind_mod(iobs) = &
               (u*COS(alpha)*SIN(azimuth*degrad) + v*COS(azimuth*degrad)) * COS(alpha) + w*SIN(alpha)

        END DO    ! loop over radar points

        IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
          WRITE(*,*) 'max rs_data(ista)%radarwind_mod = ', &
                     MAXVAL(rs_data(ista)%radwind_mod), '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%radarwind_mod = ', &
                     MINVAL(rs_data(ista)%radwind_mod), '  nobs = ', rs_data(ista)%nobs
        END IF

        DEALLOCATE(u_azgrd)
        DEALLOCATE(v_azgrd)
        DEALLOCATE(w_azgrd)

      END IF

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_radialwind_online



  !=================================================================================================
  !+ Module procedure in radar_src for the computation of the model radial winds at smoothing points
  !  under consideration of dynamical propagation path
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mod_radialwind_onsmth(lcalc)

    !-----------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radial winds on smoothing points for each
    !              radar station by
    !              - trilinear interpolation of u,v and w to smoothing points
    !              - projection onto radar main- or subbeam
    !
    ! Method:      Following Jaervinen et al., 2009, Tellus
    !
    !-----------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER              :: ista,igrd,iaz,ira,iel, &
                            nobsmax,k,m,n,offset_i,offset_j,ismth, &
                            nk,np,mo,no,ko,ih,iv,ke_
    INTEGER              :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL (KIND=dp)       ::                                                &
                            wa                                           , &
                            wl                                           , & ! interpolation weight in i-direction
                            wk                                           , & ! interpolation weight in j-direction
                            u111, u211, u121, u221,u112, u212, u122, u222, &
                            v111, v211, v121, v221,v112, v212, v122, v222, &
                            w111, w211, w121, w221,w112, w212, w122, w222, &
                            u,v,w,range,azimuth,elevation,alpha
    REAL    (KIND=dp), ALLOCATABLE     ::u_azgrd(:,:,:),v_azgrd(:,:,:),w_azgrd(:,:,:)

    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_radialwind_onsmth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_radialwind_onsmth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model radial winds:
      ALLOCATE(rs_data(ista)%radwind_mod_smth(rs_data(ista)%nsmth))

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(u_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(v_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))
        ALLOCATE(w_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),rs_grid(ista)%nal+1,ke_))

!!$ UB>> Initialize with standard missing value (= point outside model domain):
        u_azgrd = miss_value
        v_azgrd = miss_value
        w_azgrd = miss_value

        !sort the collected smoothing points data into 3d filed according to sorted azimuth, arc length, height.
!NEC$ ivdep
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current smoothing point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          u_azgrd(m,n,k) = rs_grid(ista)%u_azgrd(igrd)

          v_azgrd(m,n,k) = rs_grid(ista)%v_azgrd(igrd)

          w_azgrd(m,n,k) = rs_grid(ista)%w_azgrd(igrd)

        ENDDO

        rs_data(ista)%radwind_mod_smth = miss_value

        ! loop over observation points
!NEC$ ivdep
        DO ismth = 1, rs_data(ista)%nsmth

          ! for each smoothing point determine the continuous number nk of the grid cube
          ! surrounded by the 8 nearest auxiliary grid points
          nk = rs_data(ista)%ind_intp_smth(ismth,1)

          ! for each smoothing point determine the continuous number np of smoothing points
          np = rs_data(ista)%ind_intp_smth(ismth,2)

          ! determine indices for interpolation grids in atimuthal, arc length, vertical direction for current smoothing point
          CALL ind2sub3D(nk, rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          ! determine interpolation weights
          wa = rs_data(ista)%w_intp_smth(ismth,1)
          wl = rs_data(ista)%w_intp_smth(ismth,2)
          wk = rs_data(ista)%w_intp_smth(ismth,3)

          mo = MIN(m+1,iend(my_cart_id_fwo,ista))
          no = MIN(n+1,rs_grid(ista)%nal+1)
          ko = one_level_down(k)

          u111 = u_azgrd(m,n,k)
          u211 = u_azgrd(mo,n,k)
          u121 = u_azgrd(m,no,k)
          u221 = u_azgrd(mo,no,k)
          u112 = u_azgrd(m,n,ko)
          u212 = u_azgrd(mo,n,ko)
          u122 = u_azgrd(m,no,ko)
          u222 = u_azgrd(mo,no,ko)

          v111 = v_azgrd(m,n,k)
          v211 = v_azgrd(mo,n,k)
          v121 = v_azgrd(m,no,k)
          v221 = v_azgrd(mo,no,k)
          v112 = v_azgrd(m,n,ko)
          v212 = v_azgrd(mo,n,ko)
          v122 = v_azgrd(m,no,ko)
          v222 = v_azgrd(mo,no,ko)

          w111 = w_azgrd(m,n,k)
          w211 = w_azgrd(mo,n,k)
          w121 = w_azgrd(m,no,k)
          w221 = w_azgrd(mo,no,k)
          w112 = w_azgrd(m,n,ko)
          w212 = w_azgrd(mo,n,ko)
          w122 = w_azgrd(m,no,ko)
          w222 = w_azgrd(mo,no,ko)

          ! If the value on a grid point is missing, it will be replaced by value from the horizontal symmetric grid point;
          ! If both value are missing, neglect these two points while interpolation
          IF (u111 < miss_threshold .AND. u211 >= miss_threshold) THEN
            u111 = u211
            v111 = v211
            w111 = w211
          ELSEIF (u111 >= miss_threshold .AND. u211 < miss_threshold) THEN
            u211 = u111
            v211 = v111
            w211 = w111
          ELSEIF (u111 < miss_threshold .AND. u211 < miss_threshold) THEN
            wl = 1.0
          ELSE
            CONTINUE
          ENDIF

          IF (u121 < miss_threshold .AND. u221 >= miss_threshold) THEN
            u121 = u221
            v121 = v221
            w121 = w221
          ELSEIF (u121 >= miss_threshold .AND. u221 < miss_threshold) THEN
            u221 = u121
            v221 = v121
            w221 = w121
          ELSEIF (u121 < miss_threshold .AND. u221 < miss_threshold) THEN
            wl = 0.0
          ELSE
            CONTINUE
          ENDIF

          IF (u112 < miss_threshold .AND. u212 >= miss_threshold) THEN
            u112 = u212
            v112 = v212
            w112 = w212
          ELSEIF (u112 >= miss_threshold .AND. u212 < miss_threshold) THEN
            u212 = u112
            v212 = v112
            w212 = w112
          ELSEif (u112 < miss_threshold .AND. u212 < miss_threshold) then
            wl = 1.0
          ELSE
            CONTINUE
          ENDIF

          IF (u122 < miss_threshold .AND. u222 >= miss_threshold) THEN
            u122 = u222
            v122 = v222
            w122 = w222
          ELSEIF (u122 >= miss_threshold .AND. u222 < miss_threshold) THEN
            u222 = u122
            v222 = v122
            w222 = w122
          ELSEif (u122 < miss_threshold .AND. u222 < miss_threshold) then
            wl = 0.0
          ELSE
            CONTINUE
          ENDIF

!!$ do not interpolate if one of the model grid points is invalid:
!!$     (could completely replace the above IF-clauses)
          IF (u111 < miss_threshold .OR. u211 <  miss_threshold .OR. &
              u121 < miss_threshold .OR. u221 <  miss_threshold .OR. &
              u112 < miss_threshold .OR. u212 <  miss_threshold .OR. &
              u122 < miss_threshold .OR. u222 <  miss_threshold .OR. &
              v111 < miss_threshold .OR. v211 <  miss_threshold .OR. &
              v121 < miss_threshold .OR. v221 <  miss_threshold .OR. &
              v112 < miss_threshold .OR. v212 <  miss_threshold .OR. &
              v122 < miss_threshold .OR. v222 <  miss_threshold .OR. &
              w111 < miss_threshold .OR. w211 <  miss_threshold .OR. &
              w121 < miss_threshold .OR. w221 <  miss_threshold .OR. &
              w112 < miss_threshold .OR. w212 <  miss_threshold .OR. &
              w122 < miss_threshold .OR. w222 <  miss_threshold ) THEN

            rs_data(ista)%radwind_mod_smth(ismth) =   miss_value

          ELSE

            u = ((u111*(1.0_dp-wa)+u211*wa)*(1.0_dp - wl)   + &
                 (u121*(1.0_dp-wa)+u221*wa)*wl)*(1.0_dp-wk) + &
                 ((u112*(1.0_dp-wa)+u212*wa)*(1.0_dp - wl)  + &
                 (u122*(1.0_dp-wa)+u222*wa)*wl)*wk


            v = ((v111*(1.0_dp-wa)+v211*wa)*(1.0_dp - wl)   + &
                 (v121*(1.0_dp-wa)+v221*wa)*wl)*(1.0_dp-wk) + &
                 ((v112*(1.0_dp-wa)+v212*wa)*(1.0_dp - wl)  + &
                 (v122*(1.0_dp-wa)+v222*wa)*wl)*wk

            w = ((w111*(1.0_dp-wa)+w211*wa)*(1.0_dp - wl)   + &
                 (w121*(1.0_dp-wa)+w221*wa)*wl)*(1.0_dp-wk) + &
                 ((w112*(1.0_dp-wa)+w212*wa)*(1.0_dp - wl)  + &
                 (w122*(1.0_dp-wa)+w222*wa)*wl)*wk

            alpha = rs_data(ista)%el_loc(ismth)*degrad

            CALL ind2sub5D(np, rs_meta(ista)%naz, rs_meta(ista)%nra, rs_meta(ista)%nel, rs_meta(ista)%ngpsm_v, &
                           iaz, ira, iel, iv, ih)

!!$ inlined code from smth_az_horzscan:
!!$        phi3_eff = rs_meta(ista)%alpha3_eff_0 + &
!!$             ( COS(rs_meta(ista)%el_arr(iel)*degrad) - 1.0 ) * rs_meta(ista)%dalpha * &
!!$             ( 1.0 - EXP(-1.5*rs_meta(ista)%dalpha/rs_meta(ista)%Phi3 ) )
!!$
!!$        smthpoint = 0.5 * phi3_eff * rs_meta(ista)%smth_interv_fact * &
!!$                   rs_meta(ista)%xabscsm_h(ih)
!!$
!!$        azimuth  = (iaz-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start + smthpoint
!!$

            azimuth  = (iaz-1)*rs_meta(ista)%az_inc + rs_meta(ista)%az_start + &
                        smth_az_horzscan( &
                                         rs_meta(ista)%alpha3_eff_0, &
                                         rs_meta(ista)%dalpha, &
                                         rs_meta(ista)%Phi3, &
                                         rs_meta(ista)%smth_interv_fact, &
                                         rs_meta(ista)%xabscsm_h(ih), &
                                         rs_meta(ista)%el_arr(iel) &
                                        )
            ! Note: for DWD precip scan, rs_meta(ista)%el_arr(1) is a suitable approximation.

            azimuth  = MODULO(azimuth, 360.0_dp)

            rs_data(ista)%radwind_mod_smth(ismth) = &
                 (u*COS(alpha)*SIN(azimuth*degrad)+v*COS(azimuth*degrad))*COS(alpha) + w*SIN(alpha)

          END IF

        END DO    ! loop over radar points

        IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
          WRITE(*,*) 'max rs_data(ista)%radwind_mod_smth = ', &
                     MAXVAL(rs_data(ista)%radwind_mod_smth), ' nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%radwind_mod_smth = ', &
                     MINVAL(rs_data(ista)%radwind_mod_smth), ' nsmth = ', rs_data(ista)%nsmth
        END IF

        DEALLOCATE(u_azgrd)
        DEALLOCATE(v_azgrd)
        DEALLOCATE(w_azgrd)

      END IF

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_radialwind_onsmth


  !================================================================================================
  !+ Module procedure in radar_src for the computation of the model u,v,w winds at auxiliary grids
  !------------------------------------------------------------------------------------------------

  SUBROUTINE calc_grd_winduvw(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model u,v,w winds for each
    !              radar station by
    !              - bilinear interpolation of u,v and w to grid points
    !
    !
    !
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                :: ista

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_grd_winduvw
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_grd_winduvw'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      IF (lcalc(ista) .AND. rs_grid(ista)%ngrd >= 0) THEN

        ! allocate array of model radial winds
        ALLOCATE(rs_grid(ista)%u_grd(rs_grid(ista)%ngrd))
        ALLOCATE(rs_grid(ista)%v_grd(rs_grid(ista)%ngrd))
        ALLOCATE(rs_grid(ista)%w_grd(rs_grid(ista)%ngrd))

        CALL interp_model2azislices_vr (u, v, w, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, &
             rs_grid(ista)%u_grd, rs_grid(ista)%v_grd, rs_grid(ista)%w_grd)

      END IF

      IF (ldebug_radsim .AND. rs_grid(ista)%ngrd > 0) THEN
        WRITE(*,*) 'max rs_grid(ista)%u_grd = ', &
                   MAXVAL(rs_grid(ista)%u_grd), ' ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%u_grd = ', &
                   MINVAL(rs_grid(ista)%u_grd), ' ngrd = ', rs_grid(ista)%ngrd

        WRITE(*,*) 'max vs_grid(ista)%_grd = ', &
                   MAXVAL(rs_grid(ista)%v_grd), ' ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%v_grd = ', &
                   MINVAL(rs_grid(ista)%v_grd), ' ngrd = ', rs_grid(ista)%ngrd

        WRITE(*,*) 'max rs_grid(ista)%w_grd = ', &
                   MAXVAL(rs_grid(ista)%w_grd), ' ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%w_grd = ', &
                   MINVAL(rs_grid(ista)%w_grd), ' ngrd = ', rs_grid(ista)%ngrd
      END IF


    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_grd_winduvw


  !==============================================================================
  !+ Module procedure in radar_src for the computation of the radar reflectivity
  !------------------------------------------------------------------------------
  SUBROUTINE calc_mod_refl_modelgrid(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the radar reflectivity and, if 
    !              applicable, polarimetric radar quantities for each radar
    !              station.
    !              Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method:
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista, iobs, jb

    !- End of header
    !==============================================================================


    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_reflectivity
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_mod_refl_modelgrid'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate reflectivity on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'
      CALL init_vari(zh_radar, 0.0_dp)
      CALL init_vari(ah_radar, 0.0_dp)

      CALL init_vari(zv_radar, 0.0_dp)
      CALL init_vari(rrhv_radar, 0.0_dp)
      CALL init_vari(irhv_radar, 0.0_dp)
      CALL init_vari(kdp_radar, 0.0_dp)
      CALL init_vari(adp_radar, 0.0_dp)
      CALL init_vari(zvh_radar, 0.0_dp)

      CALL calc_dbz_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), &
                                  .TRUE., ldebug_radsim, &
                                  TRIM(ydir_mielookup_read), TRIM(ydir_mielookup_write), &
                                  zh_radar=zh_radar,ah_radar=ah_radar,&
                                  zv_radar=zv_radar,rrhv_radar=rrhv_radar,irhv_radar=irhv_radar,&
                                  kdp_radar=kdp_radar,adp_radar=adp_radar,zvh_radar=zvh_radar)
      ! convert dBZ to linear unit mm^6/m^3:
!!$      CALL dbz_to_linear(zh_radar)

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of model reflectivities:
      ALLOCATE(rs_data(ista)%zh_radar_mod(rs_data(ista)%nobs))
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        ALLOCATE(rs_data(ista)%ah_radar_mod(rs_data(ista)%nobs))
      END IF
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%irhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%kdp_radar_mod(rs_data(ista)%nobs))
        IF (loutpolall) THEN
          ALLOCATE(rs_data(ista)%zvh_radar_mod(rs_data(ista)%nobs))
        END IF
        IF (lextdbz) THEN
          ALLOCATE(rs_data(ista)%adp_radar_mod(rs_data(ista)%nobs))
        END IF
      END IF

      ! Set data points below surface to shield_value:
!$omp parallel do private(iobs)
      DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
        rs_data(ista)%zh_radar_mod(iobs) = shield_value
      END DO
!$omp end parallel do

      CALL interp_model2radarbins_scalar (zh_radar, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%zh_radar_mod)

      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
!$omp parallel do private(iobs)
        DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
          rs_data(ista)%ah_radar_mod(iobs)   = shield_value
        END DO
!$omp end parallel do
        CALL interp_model2radarbins_scalar (ah_radar, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%ah_radar_mod)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
!$omp parallel do private(iobs)
        DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
          rs_data(ista)%zv_radar_mod(iobs)    = shield_value
          rs_data(ista)%rrhv_radar_mod(iobs)  = shield_value_rhv
          rs_data(ista)%irhv_radar_mod(iobs)  = shield_value_rhv
          rs_data(ista)%kdp_radar_mod(iobs)   = shield_value
        END DO
!$omp end parallel do
        CALL interp_model2radarbins_scalar (zv_radar, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%zv_radar_mod)
        CALL interp_model2radarbins_scalar (rrhv_radar, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%rrhv_radar_mod)
        CALL interp_model2radarbins_scalar (irhv_radar, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%irhv_radar_mod)
        CALL interp_model2radarbins_scalar (kdp_radar, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%kdp_radar_mod)

        IF (loutpolall) THEN
!$omp parallel do private(iobs)
          DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
            rs_data(ista)%zvh_radar_mod(iobs)   = shield_value
          END DO
!$omp end parallel do
          CALL interp_model2radarbins_scalar (zvh_radar, rs_data(ista)%nobs_above_sfc, &
               rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%zvh_radar_mod)
        END IF

        IF (lextdbz) THEN
!$omp parallel do private(iobs)
          DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
            rs_data(ista)%adp_radar_mod(iobs)   = shield_value
          END DO
!$omp end parallel do
          CALL interp_model2radarbins_scalar (adp_radar, rs_data(ista)%nobs_above_sfc, &
               rs_data(ista)%ind_intp, rs_data(ista)%w_intp, rs_data(ista)%adp_radar_mod)
        END IF
      END IF

      IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
        WRITE(*,*) 'max rs_data(ista)%zh_radar_mod = ', &
                   MAXVAL(rs_data(ista)%zh_radar_mod), &
                   '  nobs = ', rs_data(ista)%nobs
        WRITE(*,*) 'min rs_data(ista)%zh_radar_mod = ', &
                   MINVAL(rs_data(ista)%zh_radar_mod), &
                   '  nobs = ', rs_data(ista)%nobs

        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          WRITE(*,*) 'max rs_data(ista)%ah_radar_mod = ', &
                     MAXVAL(rs_data(ista)%ah_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%ah_radar_mod = ', &
                     MINVAL(rs_data(ista)%ah_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
        END IF

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          WRITE(*,*) 'max rs_data(ista)%zv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%zv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%zv_radar_mod = ', &
                     MINVAL(rs_data(ista)%zv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%rrhv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%rrhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%rrhv_radar_mod = ', &
                     MINVAL(rs_data(ista)%rrhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%irhv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%irhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%irhv_radar_mod = ', &
                     MINVAL(rs_data(ista)%irhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%kdp_radar_mod = ', &
                     MAXVAL(rs_data(ista)%kdp_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%kdp_radar_mod = ', &
                     MINVAL(rs_data(ista)%kdp_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          IF (loutpolall) THEN
            WRITE(*,*) 'max rs_data(ista)%zvh_radar_mod = ', &
                       MAXVAL(rs_data(ista)%zvh_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%zvh_radar_mod = ', &
                       MINVAL(rs_data(ista)%zvh_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
          END IF
          IF (lextdbz) THEN
            WRITE(*,*) 'max rs_data(ista)%adp_radar_mod = ', &
                       MAXVAL(rs_data(ista)%adp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%adp_radar_mod = ', &
                       MINVAL(rs_data(ista)%adp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_refl_modelgrid


  SUBROUTINE calc_mod_refl_radarbins(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the radar reflectivity for each
    !              radar station. Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method: In contrast to calc_mod_refl_modelgrid(), the model hydrometeor fields
    !         are first interpolated to the radar bins before computing dBZ!
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista, iobs

    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)     :: tmp_loc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:)   :: &
         Tmax_i_loc, Tmax_s_loc, Tmax_g_loc, Tmax_h_loc, qnc_s_loc, &
         Tmin_g_loc, Tmin_h_loc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
         rho_loc, t_loc, qc_loc, qr_loc, qi_loc, qs_loc, qg_loc, qh_loc, &
         qnc_loc, qnr_loc, qni_loc, qns_loc, qng_loc, qnh_loc, qgl_loc, qhl_loc
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: &
         zh_radar_loc, ah_radar_loc,   &
         zv_radar_loc, rrhv_radar_loc, irhv_radar_loc, &
         kdp_radar_loc, adp_radar_loc, &
         zvh_radar_loc

    ! Actual value for neigh_tmax_melt (computation neighbourhood for the Tmax
    !  in the degree-of-melting parameterization of radar_mie_meltdegree.f90)
#ifdef __ICON__
    REAL(kind=dp), PARAMETER           :: neigh_tmax_melt = miss_value ! m
#else
    REAL(kind=dp), PARAMETER           :: neigh_tmax_melt = 5000.0 ! m
#endif

    !- End of header
    !==============================================================================


    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_refl_radarbins
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_mod_refl_radarbins'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      IF (.NOT.ASSOCIATED(qng)) THEN
        ! Set the Tmax_XX_modelgrid parameter for the parameterization of the degree of melting
        !  of each hydrometeor category as 2D fields on the model grid:
        CALL initialize_tmax_1mom_vec_par(neigh_tmax_melt, dbz_meta(ista))
      ELSE
#ifdef TWOMOM_SB
        CALL initialize_tmax_2mom_vec_par(neigh_tmax_melt, dbz_meta(ista))
#endif
      END IF

      ALLOCATE(rs_data(ista)%zh_radar_mod(rs_data(ista)%nobs))
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        ALLOCATE(rs_data(ista)%ah_radar_mod(rs_data(ista)%nobs))
      END IF
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%irhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%kdp_radar_mod(rs_data(ista)%nobs))
        IF (loutpolall) &
          ALLOCATE(rs_data(ista)%zvh_radar_mod(rs_data(ista)%nobs))
        IF (lextdbz) &
          ALLOCATE(rs_data(ista)%adp_radar_mod(rs_data(ista)%nobs))
      END IF

      ! Set data points below surface to shield_value:
      DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
        rs_data(ista)%zh_radar_mod(iobs) = shield_value
      END DO

      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
          rs_data(ista)%ah_radar_mod(iobs)   = shield_value
        END DO
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
          rs_data(ista)%zv_radar_mod(iobs)    = shield_value
          rs_data(ista)%rrhv_radar_mod(iobs)   = shield_value_rhv
          rs_data(ista)%irhv_radar_mod(iobs)   = shield_value_rhv
          rs_data(ista)%kdp_radar_mod(iobs)   = shield_value
        END DO

        IF (loutpolall) THEN
          DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
            rs_data(ista)%zvh_radar_mod(iobs)   = shield_value
          END DO
        END IF

        IF (lextdbz) THEN
          DO iobs = rs_data(ista)%nobs_above_sfc+1, rs_data(ista)%nobs
            rs_data(ista)%adp_radar_mod(iobs)   = shield_value
          END DO
        END IF
      END IF

      IF (rs_data(ista)%nobs_above_sfc == 0) THEN
        CALL finalize_tmax
        CYCLE
      END IF

      ALLOCATE(tmp_loc(rs_data(ista)%nobs_above_sfc))

      ALLOCATE(rho_loc(rs_data(ista)%nobs_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (rho, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
      rho_loc(:,1,1) = tmp_loc
      ALLOCATE(t_loc(rs_data(ista)%nobs_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (t, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
      t_loc(:,1,1) = tmp_loc
      ALLOCATE(qc_loc(rs_data(ista)%nobs_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (qc, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
      qc_loc(:,1,1) = tmp_loc
      ALLOCATE(qnc_s_loc(rs_data(ista)%nobs_above_sfc,1))
      CALL interp2d_model2radarbins_scalar (qnc_s, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
      qnc_s_loc(:,1) = tmp_loc
      ALLOCATE(qr_loc(rs_data(ista)%nobs_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (qr, rs_data(ista)%nobs_above_sfc, &
           rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
      qr_loc(:,1,1) = tmp_loc

      IF (ASSOCIATED(qi)) THEN
        ALLOCATE(qi_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qi, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qi_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_i_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_i_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmax_i_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qs)) THEN
        ALLOCATE(qs_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qs, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qs_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_s_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_s_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmax_s_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qg)) THEN
        ALLOCATE(qg_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qg, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qg_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_g_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_g_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmax_g_loc(:,1) = tmp_loc
        ALLOCATE(Tmin_g_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmin_g_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmin_g_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qh)) THEN
        ALLOCATE(qh_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qh, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qh_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_h_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_h_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmax_h_loc(:,1) = tmp_loc
        ALLOCATE(Tmin_h_loc(rs_data(ista)%nobs_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmin_h_modelgrid, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        Tmin_h_loc(:,1) = tmp_loc
      END IF

      IF (ASSOCIATED(qnc)) THEN
        ALLOCATE(qnc_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnc, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qnc_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qnr)) THEN
        ALLOCATE(qnr_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnr, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qnr_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qni)) THEN
        ALLOCATE(qni_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qni, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qni_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qns)) THEN
        ALLOCATE(qns_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qns, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qns_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qng)) THEN
        ALLOCATE(qng_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qng, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qng_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qnh)) THEN
        ALLOCATE(qnh_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnh, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qnh_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qgl)) THEN
        ALLOCATE(qgl_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qgl, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qgl_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qhl)) THEN
        ALLOCATE(qhl_loc(rs_data(ista)%nobs_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qhl, rs_data(ista)%nobs_above_sfc, &
             rs_data(ista)%ind_intp, rs_data(ista)%w_intp, tmp_loc)
        qhl_loc(:,1,1) = tmp_loc
      END IF

      ALLOCATE(zh_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               ah_radar_loc(rs_data(ista)%nobs_above_sfc,1,1))

      ALLOCATE(zv_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               rrhv_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               irhv_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               kdp_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               adp_radar_loc(rs_data(ista)%nobs_above_sfc,1,1), &
               zvh_radar_loc(rs_data(ista)%nobs_above_sfc,1,1))

      CALL calc_dbz_vec_generic (namlist=dbz_meta(ista), &
           ldebug=ldebug_radsim, &
           ydir_lookup_read=TRIM(ydir_mielookup_read), &
           ydir_lookup_write=TRIM(ydir_mielookup_write), &
           rho=rho_loc, &
           t=t_loc, &
           qc=qc_loc, &
           qr=qr_loc, &
           qi=qi_loc, &
           qs=qs_loc, &
           qg=qg_loc, &
           qh=qh_loc, &
           qnc=qnc_loc, &
           qnr=qnr_loc, &
           qni=qni_loc, &
           qns=qns_loc, &
           qng=qng_loc, &
           qnh=qnh_loc, &
           qgl=qgl_loc, &
           qhl=qhl_loc, &
           qnc_s=qnc_s_loc, &
           Tmax_i=Tmax_i_loc, &
           Tmax_s=Tmax_s_loc, &
           Tmax_g=Tmax_g_loc, &
           Tmax_h=Tmax_h_loc, &
           Tmin_g=Tmin_g_loc, &
           Tmin_h=Tmin_h_loc, &
           zh_radar=zh_radar_loc, &
           ah_radar=ah_radar_loc, &
           zv_radar=zv_radar_loc, &
           rrhv_radar=rrhv_radar_loc, &
           irhv_radar=irhv_radar_loc, &
           kdp_radar=kdp_radar_loc, &
           adp_radar=adp_radar_loc, &
           zvh_radar=zvh_radar_loc)

      ! convert dBZ to linear unit mm^6/m^3:
!!$      CALL dbz_to_linear(zh_radar_loc)

      rs_data(ista)%zh_radar_mod(1:rs_data(ista)%nobs_above_sfc) = zh_radar_loc(:,1,1)
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        rs_data(ista)%ah_radar_mod(1:rs_data(ista)%nobs_above_sfc) = ah_radar_loc(:,1,1)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        rs_data(ista)%zv_radar_mod(1:rs_data(ista)%nobs_above_sfc) = zv_radar_loc(:,1,1)
        rs_data(ista)%rrhv_radar_mod(1:rs_data(ista)%nobs_above_sfc) = rrhv_radar_loc(:,1,1)
        rs_data(ista)%irhv_radar_mod(1:rs_data(ista)%nobs_above_sfc) = irhv_radar_loc(:,1,1)
        rs_data(ista)%kdp_radar_mod(1:rs_data(ista)%nobs_above_sfc) = kdp_radar_loc(:,1,1)
        IF (loutpolall) THEN
          rs_data(ista)%zvh_radar_mod(1:rs_data(ista)%nobs_above_sfc) = zvh_radar_loc(:,1,1)
        END IF
        IF (lextdbz) THEN
          rs_data(ista)%adp_radar_mod(1:rs_data(ista)%nobs_above_sfc) = adp_radar_loc(:,1,1)
        END IF
      END IF

      DEALLOCATE(tmp_loc, rho_loc, t_loc, qc_loc, qnc_s_loc, qr_loc)
      IF (ALLOCATED(qi_loc)) DEALLOCATE(qi_loc)
      IF (ALLOCATED(qs_loc)) DEALLOCATE(qs_loc)
      IF (ALLOCATED(qg_loc)) DEALLOCATE(qg_loc)
      IF (ALLOCATED(qh_loc)) DEALLOCATE(qh_loc)
      IF (ALLOCATED(qnc_loc)) DEALLOCATE(qnc_loc)
      IF (ALLOCATED(qnr_loc)) DEALLOCATE(qnr_loc)
      IF (ALLOCATED(qni_loc)) DEALLOCATE(qni_loc)
      IF (ALLOCATED(qns_loc)) DEALLOCATE(qns_loc)
      IF (ALLOCATED(qng_loc)) DEALLOCATE(qng_loc)
      IF (ALLOCATED(qnh_loc)) DEALLOCATE(qnh_loc)
      IF (ALLOCATED(qgl_loc)) DEALLOCATE(qgl_loc)
      IF (ALLOCATED(qhl_loc)) DEALLOCATE(qhl_loc)
      IF (ALLOCATED(Tmax_i_loc)) DEALLOCATE(Tmax_i_loc)
      IF (ALLOCATED(Tmax_s_loc)) DEALLOCATE(Tmax_s_loc)
      IF (ALLOCATED(Tmax_g_loc)) DEALLOCATE(Tmax_g_loc)
      IF (ALLOCATED(Tmax_h_loc)) DEALLOCATE(Tmax_h_loc)
      IF (ALLOCATED(Tmin_g_loc)) DEALLOCATE(Tmin_g_loc)
      IF (ALLOCATED(Tmin_h_loc)) DEALLOCATE(Tmin_h_loc)

      DEALLOCATE(zh_radar_loc, ah_radar_loc)
      DEALLOCATE(zv_radar_loc, rrhv_radar_loc, irhv_radar_loc, &
                 kdp_radar_loc, adp_radar_loc, zvh_radar_loc )

      IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
        WRITE(*,*) 'max rs_data(ista)%zh_radar_mod = ', &
                   MAXVAL(rs_data(ista)%zh_radar_mod), &
                   '  nobs = ', rs_data(ista)%nobs
        WRITE(*,*) 'min rs_data(ista)%zh_radar_mod = ', &
                   MINVAL(rs_data(ista)%zh_radar_mod), &
                   '  nobs = ', rs_data(ista)%nobs

        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          WRITE(*,*) 'max rs_data(ista)%ah_radar_mod = ', &
                     MAXVAL(rs_data(ista)%ah_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%ah_radar_mod = ', &
                     MINVAL(rs_data(ista)%ah_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
        END IF

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          WRITE(*,*) 'max rs_data(ista)%zv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%zv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%zv_radar_mod = ', &
                     MINVAL(rs_data(ista)%zv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%rrhv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%rrhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%rrhv_radar_mod = ', &
                     MINVAL(rs_data(ista)%rrhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%irhv_radar_mod = ', &
                     MAXVAL(rs_data(ista)%irhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%irhv_radar_mod = ', &
                     MINVAL(rs_data(ista)%irhv_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'max rs_data(ista)%kdp_radar_mod = ', &
                     MAXVAL(rs_data(ista)%kdp_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%kdp_radar_mod = ', &
                     MINVAL(rs_data(ista)%kdp_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          IF (loutpolall) THEN
            WRITE(*,*) 'max rs_data(ista)%zvh_radar_mod = ', &
                       MAXVAL(rs_data(ista)%zvh_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%zvh_radar_mod = ', &
                       MINVAL(rs_data(ista)%zvh_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
          END IF
          IF (lextdbz) THEN
            WRITE(*,*) 'max rs_data(ista)%adp_radar_mod = ', &
                       MAXVAL(rs_data(ista)%adp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%adp_radar_mod = ', &
                       MINVAL(rs_data(ista)%adp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

      CALL finalize_tmax

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_refl_radarbins


  !==============================================================================
  !+ Module procedure in radar_src for the computation of the radar reflectivity
  !  on auxiliary rays for the beam smoothing option
  !------------------------------------------------------------------------------
  SUBROUTINE calc_mod_refl_smth_modelgrid(lcalc)

    !--------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the radar reflectivity at smoothing points for each
    !              radar station. Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method:
    !
    !--------------------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg


    INTEGER    :: ista, ismth

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_refl_smth_modelgrid
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine = 'calc_mod_refl_smth_modelgrid'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id


    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate reflectivity on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'
      CALL init_vari(zh_radar, 0.0_dp)
      CALL init_vari(ah_radar, 0.0_dp)

      CALL init_vari(zv_radar, 0.0_dp)
      CALL init_vari(rrhv_radar, 0.0_dp)
      CALL init_vari(irhv_radar, 0.0_dp)
      CALL init_vari(kdp_radar, 0.0_dp)
      CALL init_vari(adp_radar, 0.0_dp)
      CALL init_vari(zvh_radar, 0.0_dp)

      CALL calc_dbz_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), &
                                  .TRUE., ldebug_radsim, &
                                  TRIM(ydir_mielookup_read), TRIM(ydir_mielookup_write), &
                                  zh_radar=zh_radar,ah_radar=ah_radar,&
                                  zv_radar=zv_radar,rrhv_radar=rrhv_radar,irhv_radar=irhv_radar,&
                                  kdp_radar=kdp_radar,adp_radar=adp_radar,zvh_radar=zvh_radar)
      ! convert dBZ to linear unit mm^6/m^3:
!!$      CALL dbz_to_linear(zh_radar)

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of model reflectivities:
      ALLOCATE(rs_data(ista)%zh_radar_mod_smth(rs_data(ista)%nsmth))
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
        ALLOCATE(rs_data(ista)%ah_radar_mod_smth(rs_data(ista)%nsmth))
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%irhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%kdp_radar_mod_smth(rs_data(ista)%nsmth))
        IF (loutpolall) THEN
          ALLOCATE(rs_data(ista)%zvh_radar_mod_smth(rs_data(ista)%nsmth))
        END IF
        IF (lextdbz) THEN
          ALLOCATE(rs_data(ista)%adp_radar_mod_smth(rs_data(ista)%nsmth))
        END IF
      END IF

      ! Set data points below surface to shield_value (back side of the vectors):
!$omp parallel do private(ismth)
      DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
        rs_data(ista)%zh_radar_mod_smth(ismth) = shield_value
      END DO
!$omp end parallel do

      CALL interp_model2radarbins_scalar (zh_radar, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, &
           rs_data(ista)%zh_radar_mod_smth)

      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
!$omp parallel do private(ismth)
        DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
          rs_data(ista)%ah_radar_mod_smth(ismth)   = shield_value
        END DO
!$omp end parallel do
        CALL interp_model2radarbins_scalar (ah_radar, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%ah_radar_mod_smth)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
!$omp parallel do private(ismth)
        DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
          rs_data(ista)%zv_radar_mod_smth(ismth)   = shield_value
          rs_data(ista)%rrhv_radar_mod_smth(ismth) = shield_value_rhv
          rs_data(ista)%irhv_radar_mod_smth(ismth) = shield_value_rhv
          rs_data(ista)%kdp_radar_mod_smth(ismth)  = shield_value
        END DO
!$omp end parallel do
        CALL interp_model2radarbins_scalar (zv_radar, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%zv_radar_mod_smth)
        CALL interp_model2radarbins_scalar (rrhv_radar, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%rrhv_radar_mod_smth)
        CALL interp_model2radarbins_scalar (irhv_radar, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%irhv_radar_mod_smth)
        CALL interp_model2radarbins_scalar (kdp_radar, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%kdp_radar_mod_smth)

        IF (loutpolall) THEN
!$omp parallel do private(ismth)
          DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
            rs_data(ista)%zvh_radar_mod_smth(ismth) = shield_value
          END DO
!$omp end parallel do
          CALL interp_model2radarbins_scalar (zvh_radar, rs_data(ista)%nsmth_above_sfc, &
               rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%zvh_radar_mod_smth)
        END IF

        IF (lextdbz) THEN
!$omp parallel do private(ismth)
          DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
            rs_data(ista)%adp_radar_mod_smth(ismth) = shield_value
          END DO
!$omp end parallel do
          CALL interp_model2radarbins_scalar (adp_radar, rs_data(ista)%nsmth_above_sfc, &
               rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, rs_data(ista)%adp_radar_mod_smth)
        END IF
      END IF

      IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
        WRITE(*,*) 'max rs_data(ista)%zh_radar_mod_smth = ', &
                   MAXVAL(rs_data(ista)%zh_radar_mod_smth), &
                   '  nsmth = ', rs_data(ista)%nsmth
        WRITE(*,*) 'min rs_data(ista)%zh_radar_mod_smth = ', &
                   MINVAL(rs_data(ista)%zh_radar_mod_smth), &
                   '  nsmth = ', rs_data(ista)%nsmth
        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          WRITE(*,*) 'max rs_data(ista)%ah_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%ah_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%ah_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%ah_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
        END IF

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          WRITE(*,*) 'max rs_data(ista)%zv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%zv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%zv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%zv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%rrhv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%rrhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%rrhv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%rrhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%irhv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%irhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%irhv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%irhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%kdp_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%kdp_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%kdp_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%kdp_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          IF (loutpolall) THEN
            WRITE(*,*) 'max rs_data(ista)%zvh_radar_mod_smth = ', &
                       MAXVAL(rs_data(ista)%zvh_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
            WRITE(*,*) 'min rs_data(ista)%zvh_radar_mod_smth = ', &
                       MINVAL(rs_data(ista)%zvh_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
          END IF
          IF (lextdbz) THEN
            WRITE(*,*) 'max rs_data(ista)%adp_radar_mod_smth = ', &
                       MAXVAL(rs_data(ista)%adp_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
            WRITE(*,*) 'min rs_data(ista)%adp_radar_mod_smth = ', &
                       MINVAL(rs_data(ista)%adp_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_refl_smth_modelgrid

  
  SUBROUTINE calc_mod_refl_smth_radarbins(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the radar reflectivity for each
    !              radar station. Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method: In contrast to calc_mod_refl_smth_modelgrid(), the model hydrometeor fields
    !         are first interpolated to the radar bins before computing dBZ!
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista, ismth

    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)     :: tmp_loc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:)   :: &
         Tmax_i_loc, Tmax_s_loc, Tmax_g_loc, Tmax_h_loc, qnc_s_loc, &
         Tmin_g_loc, Tmin_h_loc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
         rho_loc, t_loc, qc_loc, qr_loc, qi_loc, qs_loc, qg_loc, qh_loc, &
         qnc_loc, qnr_loc, qni_loc, qns_loc, qng_loc, qnh_loc, qgl_loc, qhl_loc
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: &
         zh_radar_loc, ah_radar_loc,   &
         zv_radar_loc, rrhv_radar_loc, irhv_radar_loc, &
         kdp_radar_loc, adp_radar_loc, &
         zvh_radar_loc

    ! Actual value for neigh_tmax_melt (computation neighbourhood for the Tmax
    !  in the degree-of-melting parameterization of radar_mie_meltdegree.f90)
#ifdef __ICON__
    REAL(kind=dp), PARAMETER           :: neigh_tmax_melt = miss_value ! m
#else
    REAL(kind=dp), PARAMETER           :: neigh_tmax_melt = 5000.0 ! m
#endif

    !- End of header
    !==============================================================================


    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_refl_smth_radarbins
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_mod_refl_smth_radarbins'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      IF (.NOT.ASSOCIATED(qng)) THEN
        ! Set the Tmax_XX_modelgrid parameter for the parameterization of the degree of melting
        !  of each hydrometeor category as 2D fields on the model grid:
        CALL initialize_tmax_1mom_vec_par(neigh_tmax_melt, dbz_meta(ista))
      ELSE
#ifdef TWOMOM_SB
        CALL initialize_tmax_2mom_vec_par(neigh_tmax_melt, dbz_meta(ista))
#endif
      END IF

      ALLOCATE(rs_data(ista)%zh_radar_mod_smth(rs_data(ista)%nsmth))
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        ALLOCATE(rs_data(ista)%ah_radar_mod_smth(rs_data(ista)%nsmth))
      END IF
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%irhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%kdp_radar_mod_smth(rs_data(ista)%nsmth))
        IF (loutpolall) &
          ALLOCATE(rs_data(ista)%zvh_radar_mod_smth(rs_data(ista)%nsmth))
        IF (lextdbz) &
          ALLOCATE(rs_data(ista)%adp_radar_mod_smth(rs_data(ista)%nsmth))
      END IF

      ! Set data points below surface to shield_value:
      DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
        rs_data(ista)%zh_radar_mod_smth(ismth) = shield_value
      END DO

      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
          rs_data(ista)%ah_radar_mod_smth(ismth)   = shield_value
        END DO
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
          rs_data(ista)%zv_radar_mod_smth(ismth)    = shield_value
          rs_data(ista)%rrhv_radar_mod_smth(ismth)   = shield_value_rhv
          rs_data(ista)%irhv_radar_mod_smth(ismth)   = shield_value_rhv
          rs_data(ista)%kdp_radar_mod_smth(ismth)   = shield_value
        END DO

        IF (loutpolall) THEN
          DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
            rs_data(ista)%zvh_radar_mod_smth(ismth)   = shield_value
          END DO
        END IF

        IF (lextdbz) THEN
          DO ismth = rs_data(ista)%nsmth_above_sfc+1, rs_data(ista)%nsmth
            rs_data(ista)%adp_radar_mod_smth(ismth)   = shield_value
          END DO
        END IF
      END IF

      IF (rs_data(ista)%nsmth_above_sfc == 0) THEN
        CALL finalize_tmax
        CYCLE
      END IF

      ALLOCATE(tmp_loc(rs_data(ista)%nsmth_above_sfc))

      ALLOCATE(rho_loc(rs_data(ista)%nsmth_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (rho, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
      rho_loc(:,1,1) = tmp_loc
      ALLOCATE(t_loc(rs_data(ista)%nsmth_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (t, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
      t_loc(:,1,1) = tmp_loc
      ALLOCATE(qc_loc(rs_data(ista)%nsmth_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (qc, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
      qc_loc(:,1,1) = tmp_loc
      ALLOCATE(qnc_s_loc(rs_data(ista)%nsmth_above_sfc,1))
      CALL interp2d_model2radarbins_scalar (qnc_s, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
      qnc_s_loc(:,1) = tmp_loc
      ALLOCATE(qr_loc(rs_data(ista)%nsmth_above_sfc,1,1))
      CALL interp_model2radarbins_scalar (qr, rs_data(ista)%nsmth_above_sfc, &
           rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
      qr_loc(:,1,1) = tmp_loc

      IF (ASSOCIATED(qi)) THEN
        ALLOCATE(qi_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qi, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qi_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_i_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_i_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmax_i_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qs)) THEN
        ALLOCATE(qs_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qs, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qs_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_s_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_s_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmax_s_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qg)) THEN
        ALLOCATE(qg_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qg, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qg_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_g_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_g_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmax_g_loc(:,1) = tmp_loc
        ALLOCATE(Tmin_g_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmin_g_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmin_g_loc(:,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qh)) THEN
        ALLOCATE(qh_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qh, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qh_loc(:,1,1) = tmp_loc
        ALLOCATE(Tmax_h_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmax_h_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmax_h_loc(:,1) = tmp_loc
        ALLOCATE(Tmin_h_loc(rs_data(ista)%nsmth_above_sfc,1))
        CALL interp2d_model2radarbins_scalar (Tmin_h_modelgrid, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        Tmin_h_loc(:,1) = tmp_loc
      END IF

      IF (ASSOCIATED(qnc)) THEN
        ALLOCATE(qnc_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnc, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qnc_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qnr)) THEN
        ALLOCATE(qnr_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnr, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qnr_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qni)) THEN
        ALLOCATE(qni_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qni, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qni_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qns)) THEN
        ALLOCATE(qns_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qns, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qns_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qng)) THEN
        ALLOCATE(qng_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qng, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qng_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qnh)) THEN
        ALLOCATE(qnh_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qnh, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qnh_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qgl)) THEN
        ALLOCATE(qgl_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qgl, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qgl_loc(:,1,1) = tmp_loc
      END IF
      IF (ASSOCIATED(qhl)) THEN
        ALLOCATE(qhl_loc(rs_data(ista)%nsmth_above_sfc,1,1))
        CALL interp_model2radarbins_scalar (qhl, rs_data(ista)%nsmth_above_sfc, &
             rs_data(ista)%ind_intp_smth, rs_data(ista)%w_intp_smth, tmp_loc)
        qhl_loc(:,1,1) = tmp_loc
      END IF

      ALLOCATE(zh_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               ah_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1))

      ALLOCATE(zv_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               rrhv_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               irhv_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               kdp_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               adp_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1), &
               zvh_radar_loc(rs_data(ista)%nsmth_above_sfc,1,1))

      CALL calc_dbz_vec_generic (namlist=dbz_meta(ista), &
           ldebug=ldebug_radsim, &
           ydir_lookup_read=TRIM(ydir_mielookup_read), &
           ydir_lookup_write=TRIM(ydir_mielookup_write), &
           rho=rho_loc, &
           t=t_loc, &
           qc=qc_loc, &
           qr=qr_loc, &
           qi=qi_loc, &
           qs=qs_loc, &
           qg=qg_loc, &
           qh=qh_loc, &
           qnc=qnc_loc, &
           qnr=qnr_loc, &
           qni=qni_loc, &
           qns=qns_loc, &
           qng=qng_loc, &
           qnh=qnh_loc, &
           qgl=qgl_loc, &
           qhl=qhl_loc, &
           qnc_s=qnc_s_loc, &
           Tmax_i=Tmax_i_loc, &
           Tmax_s=Tmax_s_loc, &
           Tmax_g=Tmax_g_loc, &
           Tmax_h=Tmax_h_loc, &
           Tmin_g=Tmin_g_loc, &
           Tmin_h=Tmin_h_loc, &
           zh_radar=zh_radar_loc, &
           ah_radar=ah_radar_loc, &
           zv_radar=zv_radar_loc, &
           rrhv_radar=rrhv_radar_loc, &
           irhv_radar=irhv_radar_loc, &
           kdp_radar=kdp_radar_loc, &
           adp_radar=adp_radar_loc, &
           zvh_radar=zvh_radar_loc)

      ! convert dBZ to linear unit mm^6/m^3:
!!$      CALL dbz_to_linear(zh_radar_loc)

      rs_data(ista)%zh_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = zh_radar_loc(:,1,1)
      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        rs_data(ista)%ah_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = ah_radar_loc(:,1,1)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        rs_data(ista)%zv_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = zv_radar_loc(:,1,1)
        rs_data(ista)%rrhv_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = rrhv_radar_loc(:,1,1)
        rs_data(ista)%irhv_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = irhv_radar_loc(:,1,1)
        rs_data(ista)%kdp_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = kdp_radar_loc(:,1,1)
        IF (loutpolall) THEN
          rs_data(ista)%zvh_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = zvh_radar_loc(:,1,1)
        END IF
        IF (lextdbz) THEN
          rs_data(ista)%adp_radar_mod_smth(1:rs_data(ista)%nsmth_above_sfc) = adp_radar_loc(:,1,1)
        END IF
      END IF

      DEALLOCATE(tmp_loc, rho_loc, t_loc, qc_loc, qnc_s_loc, qr_loc)
      IF (ALLOCATED(qi_loc)) DEALLOCATE(qi_loc)
      IF (ALLOCATED(qs_loc)) DEALLOCATE(qs_loc)
      IF (ALLOCATED(qg_loc)) DEALLOCATE(qg_loc)
      IF (ALLOCATED(qh_loc)) DEALLOCATE(qh_loc)
      IF (ALLOCATED(qnc_loc)) DEALLOCATE(qnc_loc)
      IF (ALLOCATED(qnr_loc)) DEALLOCATE(qnr_loc)
      IF (ALLOCATED(qni_loc)) DEALLOCATE(qni_loc)
      IF (ALLOCATED(qns_loc)) DEALLOCATE(qns_loc)
      IF (ALLOCATED(qng_loc)) DEALLOCATE(qng_loc)
      IF (ALLOCATED(qnh_loc)) DEALLOCATE(qnh_loc)
      IF (ALLOCATED(qgl_loc)) DEALLOCATE(qgl_loc)
      IF (ALLOCATED(qhl_loc)) DEALLOCATE(qhl_loc)
      IF (ALLOCATED(Tmax_i_loc)) DEALLOCATE(Tmax_i_loc)
      IF (ALLOCATED(Tmax_s_loc)) DEALLOCATE(Tmax_s_loc)
      IF (ALLOCATED(Tmax_g_loc)) DEALLOCATE(Tmax_g_loc)
      IF (ALLOCATED(Tmax_h_loc)) DEALLOCATE(Tmax_h_loc)
      IF (ALLOCATED(Tmin_g_loc)) DEALLOCATE(Tmin_g_loc)
      IF (ALLOCATED(Tmin_h_loc)) DEALLOCATE(Tmin_h_loc)

      DEALLOCATE(zh_radar_loc, ah_radar_loc)
      DEALLOCATE(zv_radar_loc, rrhv_radar_loc, irhv_radar_loc, &
                 kdp_radar_loc, adp_radar_loc, zvh_radar_loc )

      IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
        WRITE(*,*) 'max rs_data(ista)%zh_radar_mod_smth = ', &
                   MAXVAL(rs_data(ista)%zh_radar_mod_smth), &
                   '  nsmth = ', rs_data(ista)%nsmth
        WRITE(*,*) 'min rs_data(ista)%zh_radar_mod_smth = ', &
                   MINVAL(rs_data(ista)%zh_radar_mod_smth), &
                   '  nsmth = ', rs_data(ista)%nsmth

        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          WRITE(*,*) 'max rs_data(ista)%ah_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%ah_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%ah_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%ah_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
        END IF

        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          WRITE(*,*) 'max rs_data(ista)%zv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%zv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%zv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%zv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%rrhv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%rrhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%rrhv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%rrhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%irhv_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%irhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%irhv_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%irhv_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'max rs_data(ista)%kdp_radar_mod_smth = ', &
                     MAXVAL(rs_data(ista)%kdp_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          WRITE(*,*) 'min rs_data(ista)%kdp_radar_mod_smth = ', &
                     MINVAL(rs_data(ista)%kdp_radar_mod_smth), &
                     '  nsmth = ', rs_data(ista)%nsmth
          IF (loutpolall) THEN
            WRITE(*,*) 'max rs_data(ista)%zvh_radar_mod_smth = ', &
                       MAXVAL(rs_data(ista)%zvh_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
            WRITE(*,*) 'min rs_data(ista)%zvh_radar_mod_smth = ', &
                       MINVAL(rs_data(ista)%zvh_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
          END IF
          IF (lextdbz) THEN
            WRITE(*,*) 'max rs_data(ista)%adp_radar_mod_smth = ', &
                       MAXVAL(rs_data(ista)%adp_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
            WRITE(*,*) 'min rs_data(ista)%adp_radar_mod_smth = ', &
                       MINVAL(rs_data(ista)%adp_radar_mod_smth), &
                       '  nsmth = ', rs_data(ista)%nsmth
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

      CALL finalize_tmax

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_mod_refl_smth_radarbins


  !===================================================================================
  !+ Module procedure in radar_src for the computation of the model radar reflectivity
  !  under consideration of dynamical propagation path
  !-----------------------------------------------------------------------------------
  SUBROUTINE calc_mod_reflectivity_online(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radial reflectivity for each
    !              radar station by
    !              - bilinear interpolation of u,v and w to radar points
    !              - projection onto radar beam
    !
    ! Input and output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method:      Following Jaervinen et al., 2009, Tellus
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                :: ista,igrd,nobsmax,nobs,i,j,k,m,n,o,offset_i,offset_j,iobs,nk,np,no,ko,ke_
    INTEGER                :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)
    REAL(KIND=dp)          :: wl                ,& ! interpolation weight in i-direction
                              wk                   ! interpolation weight in j-direction

    REAL(KIND=dp), ALLOCATABLE, TARGET :: &
         zh_radar_azgrd(:,:,:), ah_radar_azgrd(:,:,:),   &
         zv_radar_azgrd(:,:,:), rrhv_radar_azgrd(:,:,:), irhv_radar_azgrd(:,:,:), &
         kdp_radar_azgrd(:,:,:), adp_radar_azgrd(:,:,:), &
         zvh_radar_azgrd(:,:,:)
    
    REAL(KIND=dp), POINTER, DIMENSION(:,:,:)  :: pzh, pzv, pah, prrhv, pirhv, pkdp, pzvh, padp
    
    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_reflectivity_online
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_reflectivity_online'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model reflectivity & extinction
      ALLOCATE(rs_data(ista)%zh_radar_mod(rs_data(ista)%nobs))
      IF (lextdbz .AND. &
          (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
        ALLOCATE(rs_data(ista)%ah_radar_mod(rs_data(ista)%nobs))
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%irhv_radar_mod(rs_data(ista)%nobs))
        ALLOCATE(rs_data(ista)%kdp_radar_mod(rs_data(ista)%nobs))
        IF (loutpolall) THEN
          ALLOCATE(rs_data(ista)%zvh_radar_mod(rs_data(ista)%nobs))
        END IF
        IF (lextdbz) THEN
          ALLOCATE(rs_data(ista)%adp_radar_mod(rs_data(ista)%nobs))
        END IF
      END IF

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(zh_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                       rs_grid(ista)%nal+1,ke_))
!!$ UB>> Initialize with standard missing value (= point outside model domain):
        zh_radar_azgrd = miss_value
        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          ALLOCATE(ah_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                  rs_grid(ista)%nal+1,ke_))
          ah_radar_azgrd = miss_value
        END IF
        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          ALLOCATE(zv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                  rs_grid(ista)%nal+1,ke_))
          zv_radar_azgrd = miss_value
          ALLOCATE(rrhv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                    rs_grid(ista)%nal+1,ke_))
          rrhv_radar_azgrd = miss_value_rhv
          ALLOCATE(irhv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                    rs_grid(ista)%nal+1,ke_))
          irhv_radar_azgrd = miss_value_rhv
          ALLOCATE(kdp_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                   rs_grid(ista)%nal+1,ke_))
          kdp_radar_azgrd = miss_value
          IF (loutpolall) THEN
            ALLOCATE(zvh_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                     rs_grid(ista)%nal+1,ke_))
            zvh_radar_azgrd = miss_value
          END IF
          IF (lextdbz) THEN
            ALLOCATE(adp_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                     rs_grid(ista)%nal+1,ke_))
            adp_radar_azgrd = miss_value
          END IF
        END IF

        !sort the collected smoothing points data into 3d field according to sorted azimuth, arc length, height.
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, &
                         rs_grid(ista)%nal+1, m, n, k)

          zh_radar_azgrd(m,n,k) = rs_grid(ista)%zh_radar_azgrd(igrd)
          IF (lextdbz .AND. &
              (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
            ah_radar_azgrd(m,n,k) = rs_grid(ista)%ah_radar_azgrd(igrd)
          IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
            zv_radar_azgrd(m,n,k) = rs_grid(ista)%zv_radar_azgrd(igrd)
            rrhv_radar_azgrd(m,n,k) = rs_grid(ista)%rrhv_radar_azgrd(igrd)
            irhv_radar_azgrd(m,n,k) = rs_grid(ista)%irhv_radar_azgrd(igrd)
            kdp_radar_azgrd(m,n,k) = rs_grid(ista)%kdp_radar_azgrd(igrd)
            IF (loutpolall) THEN
              zvh_radar_azgrd(m,n,k) = rs_grid(ista)%zvh_radar_azgrd(igrd)
            END IF
            IF (lextdbz) THEN
              adp_radar_azgrd(m,n,k) = rs_grid(ista)%adp_radar_azgrd(igrd)
            END IF
          END IF

        ENDDO

        rs_data(ista)%zh_radar_mod = miss_value
        IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          rs_data(ista)%ah_radar_mod = miss_value
        END IF
        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          rs_data(ista)%zv_radar_mod = miss_value
          rs_data(ista)%rrhv_radar_mod = miss_value_rhv
          rs_data(ista)%irhv_radar_mod = miss_value_rhv
          rs_data(ista)%kdp_radar_mod = miss_value
          IF (loutpolall) rs_data(ista)%zv_radar_mod = miss_value
          IF (lextdbz) rs_data(ista)%adp_radar_mod = miss_value
        END IF

        ! loop over observation points
        pzh => zh_radar_azgrd  ! for convenience
        pzv => zv_radar_azgrd
        pah => ah_radar_azgrd
        prrhv => rrhv_radar_azgrd
        pirhv => irhv_radar_azgrd
        pkdp => kdp_radar_azgrd
        pzvh => zvh_radar_azgrd
        padp => adp_radar_azgrd
!NEC$ ivdep
        DO iobs = 1, rs_data(ista)%nobs

          ! for each radar point determine the continuous number nk of the grid cube
          ! surrounded by the 8 nearest model grid points
          nk = rs_data(ista)%ind_intp(iobs,1)

          ! for each radar point determine the continuous number np of radar points
          np = rs_data(ista)%ind_intp(iobs,2)

          ! determine indices for interpolation grids in azimuthal, arc length,
          ! vertical direction for current radar point
          CALL ind2sub3D(nk, rs_meta(ista)%naz, rs_grid(ista)%nal+1, m, n, k)
          m = m + nbl_az

          ! determine interpolation weights
          wl = rs_data(ista)%w_intp(iobs,1)
          wk = rs_data(ista)%w_intp(iobs,2)

          no= MIN(n+1,rs_grid(ista)%nal+1)
          ko= one_level_down(k)

          rs_data(ista)%zh_radar_mod(iobs) = &
               boxinterp2D(miss_threshold,pzh(m,n,k),pzh(m,no,k),pzh(m,n,ko),pzh(m,no,ko),wl,wk)
          IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
            rs_data(ista)%ah_radar_mod(iobs) = &
                 boxinterp2D(miss_threshold,pah(m,n,k),pah(m,no,k),pah(m,n,ko),pah(m,no,ko),wl,wk)
          END IF
          IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
            rs_data(ista)%zv_radar_mod(iobs)   = &
                 boxinterp2D(miss_threshold,pzv(m,n,k),pzv(m,no,k),pzv(m,n,ko),pzv(m,no,ko),wl,wk)
            rs_data(ista)%rrhv_radar_mod(iobs) = &
                 boxinterp2D(miss_thresh_rhv,prrhv(m,n,k),prrhv(m,no,k),prrhv(m,n,ko),prrhv(m,no,ko),wl,wk)
            rs_data(ista)%irhv_radar_mod(iobs) = &
                 boxinterp2D(miss_thresh_rhv,pirhv(m,n,k),pirhv(m,no,k),pirhv(m,n,ko),pirhv(m,no,ko),wl,wk)
            rs_data(ista)%kdp_radar_mod(iobs)  = &
                 boxinterp2D(miss_threshold,pkdp(m,n,k),pkdp(m,no,k),pkdp(m,n,ko),pkdp(m,no,ko),wl,wk)
            IF (loutpolall) THEN
              rs_data(ista)%zvh_radar_mod(iobs) = &
                   boxinterp2D(miss_threshold,pzvh(m,n,k),pzvh(m,no,k),pzvh(m,n,ko),pzvh(m,no,ko),wl,wk)
            END IF
            IF (lextdbz) THEN
              rs_data(ista)%adp_radar_mod(iobs) = &
                   boxinterp2D(miss_threshold,padp(m,n,k),padp(m,no,k),padp(m,n,ko),padp(m,no,ko),wl,wk)
            END IF
          END IF

        END DO    ! loop over radar points

        IF (ldebug_radsim .AND. rs_data(ista)%nobs > 0) THEN
          WRITE(*,*) 'max rs_data(ista)%zh_radar_mod = ', &
                     MAXVAL(rs_data(ista)%zh_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          WRITE(*,*) 'min rs_data(ista)%zh_radar_mod = ', &
                     MINVAL(rs_data(ista)%zh_radar_mod), &
                     '  nobs = ', rs_data(ista)%nobs
          IF (lextdbz .AND. &
              (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
            WRITE(*,*) 'max rs_data(ista)%ah_radar_mod = ', &
                       MAXVAL(rs_data(ista)%ah_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%ah_radar_mod = ', &
                       MINVAL(rs_data(ista)%ah_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
          END IF
          IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
            WRITE(*,*) 'max rs_data(ista)%zv_radar_mod = ', &
                       MAXVAL(rs_data(ista)%zv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%zv_radar_mod = ', &
                       MINVAL(rs_data(ista)%zv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'max rs_data(ista)%rrhv_radar_mod = ', &
                       MAXVAL(rs_data(ista)%rrhv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%rrhv_radar_mod = ', &
                       MINVAL(rs_data(ista)%rrhv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'max rs_data(ista)%irhv_radar_mod = ', &
                       MAXVAL(rs_data(ista)%irhv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%irhv_radar_mod = ', &
                       MINVAL(rs_data(ista)%irhv_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'max rs_data(ista)%kdp_radar_mod = ', &
                       MAXVAL(rs_data(ista)%kdp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            WRITE(*,*) 'min rs_data(ista)%kdp_radar_mod = ', &
                       MINVAL(rs_data(ista)%kdp_radar_mod), &
                       '  nobs = ', rs_data(ista)%nobs
            IF (loutpolall) THEN
              WRITE(*,*) 'max rs_data(ista)%zvh_radar_mod = ', &
                         MAXVAL(rs_data(ista)%zvh_radar_mod), &
                         '  nobs = ', rs_data(ista)%nobs
              WRITE(*,*) 'min rs_data(ista)%zvh_radar_mod = ', &
                         MINVAL(rs_data(ista)%zvh_radar_mod), &
                         '  nobs = ', rs_data(ista)%nobs
            END IF
            IF (lextdbz) THEN
              WRITE(*,*) 'max rs_data(ista)%adp_radar_mod = ', &
                         MAXVAL(rs_data(ista)%adp_radar_mod), &
                         '  nobs = ', rs_data(ista)%nobs
              WRITE(*,*) 'min rs_data(ista)%adp_radar_mod = ', &
                         MINVAL(rs_data(ista)%adp_radar_mod), &
                         '  nobs = ', rs_data(ista)%nobs
            END IF
          END IF
        END IF

        DEALLOCATE(zh_radar_azgrd)
        IF (ALLOCATED(ah_radar_azgrd))   DEALLOCATE(ah_radar_azgrd)

        IF (ALLOCATED(zv_radar_azgrd))   DEALLOCATE(zv_radar_azgrd)
        IF (ALLOCATED(rrhv_radar_azgrd)) DEALLOCATE(rrhv_radar_azgrd)
        IF (ALLOCATED(irhv_radar_azgrd)) DEALLOCATE(irhv_radar_azgrd)
        IF (ALLOCATED(kdp_radar_azgrd))  DEALLOCATE(kdp_radar_azgrd)
        IF (ALLOCATED(adp_radar_azgrd))  DEALLOCATE(adp_radar_azgrd)
        IF (ALLOCATED(zvh_radar_azgrd))  DEALLOCATE(zvh_radar_azgrd)

      END IF ! .NOT.lcalc(ista)

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    FUNCTION boxinterp2D(miss_threshold,f11,f21,f12,f22,wl,wk) RESULT(intval)
      
      REAL(KIND=dp), INTENT(in) :: miss_threshold, f11, f12, f21, f22
      REAL(KIND=dp), INTENT(in) :: wl,wk

      REAL(KIND=dp)             :: intval
      
      REAL(KIND=dp)             :: v11, v12, v21, v22, hwk

      ! .. For clarity: these should be the indices n,no,k,ko of the input field elements:
      !    f11 = field(n,k)
      !    f21 = field(no,k)
      !    f12 = field(n,ko)
      !    f22 = field(no,ko)

      v11 = f11
      v21 = f21
      v12 = f12
      v22 = f22

      hwk = wk

      ! .. If the value on a grid point is missing, it will be replaced by value
      !    from the horizontal symmetric grid point;
      !    If both value are missing, neglect these two points while intepolation
      IF (v11 < miss_threshold .AND. v21 >= miss_threshold)     THEN
        v11 = v21
      ELSEIF (v11 >= miss_threshold .AND. v21 < miss_threshold) THEN
        v21 = v11
      ELSEIF (v11 < miss_threshold .AND. v21 < miss_threshold) THEN
        hwk = 1.0
      ELSE
        CONTINUE
      ENDIF

      IF (v12 < miss_threshold .AND. v22 >= miss_threshold)     THEN
        v12 = v22
      ELSEIF (v12 >= miss_threshold .AND. v22 < miss_threshold) THEN
        v22 = v12
      ELSEIF (v12 < miss_threshold .AND. v22 < miss_threshold) THEN
        hwk = 0.0
      ELSE
        CONTINUE
      ENDIF

      intval = (v11*(1.0_dp-wl)+v21*wl)*(1.0_dp - hwk) + (v12*(1.0_dp-wl)+v22*wl)*hwk

    END FUNCTION boxinterp2D

  END SUBROUTINE calc_mod_reflectivity_online


  !=======================================================================================================
  !+ Module procedure in radar_src for the computation of the model radar reflectivity at smoothing points
  !  under consideration of dynamical propagation path
  !-------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mod_reflectivity_onsmth(lcalc)

    !---------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the model radar reflectivity for each
    !              radar station by
    !              - trilinear interpolation of radar reflectivity to smoothing points.
    ! Input and output values are in linear space, not logarithmic (dBZ)!
    !
    !---------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                :: istart(0:num_compute_fwo-1,nradsta),iend(0:num_compute_fwo-1,nradsta)

    INTEGER                :: ista,igrd,nobsmax,ismth,i,j,k,m,n,o,offset_i,offset_j,iobs,nk,mo,no,ko,ke_
    REAL(KIND=dp)          :: wa                ,& ! interpolation weight in i-direction
                              wl                ,& ! interpolation weight in j-direction
                              wk                   ! interpolation weight in j-direction

    REAL(KIND=dp), ALLOCATABLE, TARGET :: &
         zh_radar_azgrd(:,:,:), ah_radar_azgrd(:,:,:), &
         zv_radar_azgrd(:,:,:), rrhv_radar_azgrd(:,:,:), irhv_radar_azgrd(:,:,:), &
         kdp_radar_azgrd(:,:,:), adp_radar_azgrd(:,:,:), &
         zvh_radar_azgrd(:,:,:)

    REAL(KIND=dp), POINTER, DIMENSION(:,:,:)  :: pzh, pzv, pah, prrhv, pirhv, pkdp, pzvh, padp
    
    INTEGER, PARAMETER :: ones(nradsta_max) = 1

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_mod_reflectivity_onsmth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------


    yzroutine(:) = ' '
    yzroutine = 'calc_mod_reflectivity_onsmth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ke_ = nlevels()
    
    CALL para_range_all(ones(1:nradsta),rs_meta(1:nradsta)%naz,nradsta,num_compute_fwo,nbl_az,istart,iend)

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of model radial winds
      ALLOCATE(rs_data(ista)%zh_radar_mod_smth(rs_data(ista)%nsmth))
      IF (lextdbz .AND. &
          (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
        ALLOCATE(rs_data(ista)%ah_radar_mod_smth(rs_data(ista)%nsmth))
      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_data(ista)%zv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%rrhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%irhv_radar_mod_smth(rs_data(ista)%nsmth))
        ALLOCATE(rs_data(ista)%kdp_radar_mod_smth(rs_data(ista)%nsmth))
        IF (loutpolall) &
          ALLOCATE(rs_data(ista)%zvh_radar_mod_smth(rs_data(ista)%nsmth))
        IF (lextdbz) &
          ALLOCATE(rs_data(ista)%adp_radar_mod_smth(rs_data(ista)%nsmth))
      END IF

      IF ( iend(my_cart_id_fwo,ista) > -1 ) THEN

        ALLOCATE(zh_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                rs_grid(ista)%nal+1,ke_))
!!$ UB>> Initialize with standard missing value (= point outside model domain):
        zh_radar_azgrd = miss_value

        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          ALLOCATE(ah_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                  rs_grid(ista)%nal+1,ke_))
          ah_radar_azgrd = miss_value
        END IF
        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          ALLOCATE(zv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                  rs_grid(ista)%nal+1,ke_))
          zv_radar_azgrd = miss_value_rhv
          ALLOCATE(rrhv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                    rs_grid(ista)%nal+1,ke_))
          rrhv_radar_azgrd = miss_value
          ALLOCATE(irhv_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                    rs_grid(ista)%nal+1,ke_))
          irhv_radar_azgrd = miss_value_rhv
          ALLOCATE(kdp_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                   rs_grid(ista)%nal+1,ke_))
          kdp_radar_azgrd = miss_value
          IF (loutpolall) THEN
            ALLOCATE(zvh_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                     rs_grid(ista)%nal+1,ke_))
            zvh_radar_azgrd = miss_value
          END IF
          IF (lextdbz) THEN
            ALLOCATE(adp_radar_azgrd(istart(my_cart_id_fwo,ista):iend(my_cart_id_fwo,ista),&
                                     rs_grid(ista)%nal+1,ke_))
            adp_radar_azgrd = miss_value
          END IF
        END IF

        !sort the collected smoothing points data into 3d filed according to sorted azimuth, arc length, heigh
!NEC$ ivdep
        DO igrd = 1, rs_grid(ista)%nazgrd

          ! determine indices in az, al, hl direction for current grid point
          CALL ind2sub3D(rs_grid(ista)%ind_azgrd(igrd), rs_grid(ista)%naz_nbl, &
                         rs_grid(ista)%nal+1, m, n, k)

          zh_radar_azgrd(m,n,k) = rs_grid(ista)%zh_radar_azgrd(igrd)
          IF (lextdbz .AND. &
              (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
            ah_radar_azgrd(m,n,k) = rs_grid(ista)%ah_radar_azgrd(igrd)
          IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
            zv_radar_azgrd(m,n,k) = rs_grid(ista)%zv_radar_azgrd(igrd)
            rrhv_radar_azgrd(m,n,k) = rs_grid(ista)%rrhv_radar_azgrd(igrd)
            irhv_radar_azgrd(m,n,k) = rs_grid(ista)%irhv_radar_azgrd(igrd)
            kdp_radar_azgrd(m,n,k) = rs_grid(ista)%kdp_radar_azgrd(igrd)
            IF (loutpolall) THEN
              zvh_radar_azgrd(m,n,k) = rs_grid(ista)%zvh_radar_azgrd(igrd)
            END IF
            IF (lextdbz) THEN
              adp_radar_azgrd(m,n,k) = rs_grid(ista)%adp_radar_azgrd(igrd)
            END IF
          END IF

        ENDDO

        rs_data(ista)%zh_radar_mod_smth = miss_value
        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) &
          rs_data(ista)%ah_radar_mod_smth = miss_value
        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          rs_data(ista)%zv_radar_mod_smth = miss_value
          rs_data(ista)%rrhv_radar_mod_smth = miss_value_rhv
          rs_data(ista)%irhv_radar_mod_smth = miss_value_rhv
          rs_data(ista)%kdp_radar_mod_smth = miss_value
          IF (loutpolall) &
            rs_data(ista)%zvh_radar_mod_smth = miss_value
          IF (lextdbz) &
            rs_data(ista)%adp_radar_mod_smth = miss_value
        END IF

        ! loop over observation points
        ! UB>> loop rewritten due to vectorisation:
        pzh => zh_radar_azgrd  ! for convenience
        pzv => zv_radar_azgrd
        pah => ah_radar_azgrd
        prrhv => rrhv_radar_azgrd
        pirhv => irhv_radar_azgrd
        pkdp => kdp_radar_azgrd
        pzvh => zvh_radar_azgrd
        padp => adp_radar_azgrd
!NEC$ ivdep
        DO ismth = 1, rs_data(ista)%nsmth

          ! for each radar point determine the continuous number nk of the grid cube
          ! surrounded by the 8 nearest model grid points
          nk = rs_data(ista)%ind_intp_smth(ismth,1)

          CALL ind2sub3D(nk, rs_grid(ista)%naz_nbl, rs_grid(ista)%nal+1, m, n, k)

          mo= MIN(m+1,iend(my_cart_id_fwo,ista))
          no= MIN(n+1,rs_grid(ista)%nal+1)
          ko= one_level_down(k)

          wa = rs_data(ista)%w_intp_smth(ismth,1)
          wl = rs_data(ista)%w_intp_smth(ismth,2)
          wk = rs_data(ista)%w_intp_smth(ismth,3)

          ! do not interpolate if all of the model grid points are invalid:
          IF (pzh(m, n, k)  >= miss_threshold .OR. &
              pzh(mo,n, k)  >= miss_threshold .OR. &
              pzh(m, no,k)  >= miss_threshold .OR. &
              pzh(mo,no,k)  >= miss_threshold .OR. &
              pzh(m, n, ko) >= miss_threshold .OR. &
              pzh(mo,n, ko) >= miss_threshold .OR. &
              pzh(m, no,ko) >= miss_threshold .OR. &
              pzh(mo,no,ko) >= miss_threshold) THEN

            rs_data(ista)%zh_radar_mod_smth(ismth) = &
                 boxinterp3D(miss_threshold,&
                 &          pzh(m,n,k ),pzh(mo,n,k ),pzh(m,no,k ),pzh(mo,no,k ), &
                 &           pzh(m,n,ko),pzh(mo,n,ko),pzh(m,no,ko),pzh(mo,no,ko), &
                 &           wa,wl,wk)

            IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
              rs_data(ista)%ah_radar_mod_smth(ismth) = &
                   boxinterp3D(miss_threshold,&
                   &           pah(m,n,k ),pah(mo,n,k ),pah(m,no,k ),pah(mo,no,k ), &
                   &           pah(m,n,ko),pah(mo,n,ko),pah(m,no,ko),pah(mo,no,ko), &
                   &           wa,wl,wk)
            END IF
            IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
              rs_data(ista)%zv_radar_mod_smth(ismth) = &
                   boxinterp3D(miss_threshold,&
                   &           pzv(m,n,k ),pzv(mo,n,k ),pzv(m,no,k ),pzv(mo,no,k ), &
                   &           pzv(m,n,ko),pzv(mo,n,ko),pzv(m,no,ko),pzv(mo,no,ko), &
                   &           wa,wl,wk)
              rs_data(ista)%rrhv_radar_mod_smth(ismth) = &
                   boxinterp3D(miss_thresh_rhv,&
                   &           prrhv(m,n,k ),prrhv(mo,n,k ),prrhv(m,no,k ),prrhv(mo,no,k ), &
                   &           prrhv(m,n,ko),prrhv(mo,n,ko),prrhv(m,no,ko),prrhv(mo,no,ko), &
                   &           wa,wl,wk)
              rs_data(ista)%irhv_radar_mod_smth(ismth) = &
                   boxinterp3D(miss_thresh_rhv,&
                   &           pirhv(m,n,k ),pirhv(mo,n,k ),pirhv(m,no,k ),pirhv(mo,no,k ), &
                   &           pirhv(m,n,ko),pirhv(mo,n,ko),pirhv(m,no,ko),pirhv(mo,no,ko), &
                   &           wa,wl,wk)
              rs_data(ista)%kdp_radar_mod_smth(ismth) = &
                   boxinterp3D(miss_threshold,&
                   &           pkdp(m,n,k ),pkdp(mo,n,k ),pkdp(m,no,k ),pkdp(mo,no,k ), &
                   &           pkdp(m,n,ko),pkdp(mo,n,ko),pkdp(m,no,ko),pkdp(mo,no,ko), &
                   &           wa,wl,wk)
              IF (loutpolall) THEN
                rs_data(ista)%zvh_radar_mod_smth(ismth) = &
                     boxinterp3D(miss_threshold,&
                     &           pzvh(m,n,k ),pzvh(mo,n,k ),pzvh(m,no,k ),pzvh(mo,no,k ), &
                     &           pzvh(m,n,ko),pzvh(mo,n,ko),pzvh(m,no,ko),pzvh(mo,no,ko), &
                     &           wa,wl,wk)
              END IF
              IF (lextdbz) THEN
                rs_data(ista)%adp_radar_mod_smth(ismth) = &
                     boxinterp3D(miss_threshold,&
                     &           padp(m,n,k ),padp(mo,n,k ),padp(m,no,k ),padp(mo,no,k ), &
                     &           padp(m,n,ko),padp(mo,n,ko),padp(m,no,ko),padp(mo,no,ko), &
                     &           wa,wl,wk)
              END IF
            END IF
          END IF

        END DO    ! loop over radar points

!!$        IF (ldebug_radsim .AND. rs_data(ista)%nsmth > 0) THEN
!!$          WRITE(*,*) 'max rs_data(ista)%zh_radar_mod = ', &
!!$                     MAXVAL(rs_data(ista)%zh_radar_mod_smth), &
!!$                     ' nsmth = ', rs_data(ista)%nsmth
!!$          WRITE(*,*) 'min rs_data(ista)%zh_radar_mod = ', &
!!$                     MINVAL(rs_data(ista)%zh_radar_mod_smth), &
!!$                     ' nsmth = ', rs_data(ista)%nsmth
!!$          IF (lextdbz .AND. dbz_meta(ista)%itype_refl == 1) THEN
!!$            WRITE(*,*) 'max rs_data(ista)%z_ext_mod = ', &
!!$                       MAXVAL(rs_data(ista)%z_ext_mod_smth), &
!!$                       ' nsmth = ', rs_data(ista)%nsmth
!!$            WRITE(*,*) 'min rs_data(ista)%z_ext_mod = ', &
!!$                       MINVAL(rs_data(ista)%z_ext_mod_smth), &
!!$                       ' nsmth = ', rs_data(ista)%nsmth
!!$          END IF
!!$        END IF

        DEALLOCATE(zh_radar_azgrd)
        IF (ALLOCATED(ah_radar_azgrd))   DEALLOCATE(ah_radar_azgrd)

        IF (ALLOCATED(zv_radar_azgrd))   DEALLOCATE(zv_radar_azgrd)
        IF (ALLOCATED(rrhv_radar_azgrd)) DEALLOCATE(rrhv_radar_azgrd)
        IF (ALLOCATED(irhv_radar_azgrd)) DEALLOCATE(irhv_radar_azgrd)
        IF (ALLOCATED(kdp_radar_azgrd))  DEALLOCATE(kdp_radar_azgrd)
        IF (ALLOCATED(adp_radar_azgrd))  DEALLOCATE(adp_radar_azgrd)
        IF (ALLOCATED(zvh_radar_azgrd))  DEALLOCATE(zvh_radar_azgrd)

      END IF

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  CONTAINS

    FUNCTION boxinterp3D(miss_threshold,f111,f211,f121,f221,f112,f212,f122,f222,wa,wl,wk) RESULT(intval)

      REAL(KIND=dp), INTENT(in) :: miss_threshold
      REAL(KIND=dp), INTENT(in) :: f111, f211, f121, f221, f112, f212, f122, f222
      REAL(KIND=dp), INTENT(in) :: wa,wl,wk
      REAL(KIND=dp)             :: intval

      REAL(KIND=dp)             :: v111, v211, v121, v221, v112, v212, v122, v222, hwl

      
      ! .. For clarity: these should be the indices m,mo,n,no,k,ko of the input field elements:
      !    f111 = field(m,n,k)
      !    f211 = field(mo,n,k)
      !    f121 = field(m,no,k)
      !    f221 = field(mo,no,k)
      !    f112 = field(m,n,ko)
      !    f212 = field(mo,n,ko)
      !    f122 = field(m,no,ko)
      !    f222 = field(mo,no,ko)

      v111 = f111
      v211 = f211
      v121 = f121
      v221 = f221
      v112 = f112
      v212 = f212
      v122 = f122
      v222 = f222
      
      hwl = wl

      ! .. If the value on a grid point is missing miss_value (i.e., is outside the
      !    model domain), it will be replaced by VALUE from the horizontal
      !    symmetric grid point;
      !    If both values are missing, neglect these two points while interpolation
      IF (v111 < miss_threshold .AND. v211 >= miss_threshold) THEN
        v111 = v211
      ELSEIF (v111 >= miss_threshold .AND. v211 < miss_threshold) THEN
        v211 = v111
      ELSEIF (v111 < miss_threshold .AND. v211 < miss_threshold) THEN
        hwl = 1.0
      ELSE
        CONTINUE
      ENDIF

      IF (v121 < miss_threshold .AND. v221 >= miss_threshold) THEN
        v121 = v221
      ELSEIF (v121 >= miss_threshold .AND. v221 < miss_threshold) THEN
        v221 = v121
      ELSEIF (v121 < miss_threshold .AND. v221 < miss_threshold) THEN
        hwl = 0.0
      ELSE
        CONTINUE
      ENDIF

      IF (v112 < miss_threshold .AND. v212 >= miss_threshold) THEN
        v112 = v212
      ELSEIF (v112 >= miss_threshold .AND. v212 < miss_threshold) THEN
        v212 = v112
      ELSEIF (v112 < miss_threshold .AND. v212 < miss_threshold) THEN
        hwl = 1.0
      ELSE
        CONTINUE
      ENDIF

      IF (v122 < miss_threshold .AND. v222 >= miss_threshold) THEN
        v122 = v222
      ELSEIF (v122 >= miss_threshold .AND. v222 < miss_threshold) THEN
        v222 = v122
      ELSEIF (v122 < miss_threshold .AND. v222 < miss_threshold) THEN
        hwl = 0.0
      ELSE
        CONTINUE
      ENDIF

      intval = ((v111*(1.0_dp-wa)+v211*wa)*(1.0_dp - hwl) + &
                (v121*(1.0_dp-wa)+v221*wa)*          hwl ) * (1.0_dp-wk) + &
               ((v112*(1.0_dp-wa)+v212*wa)*(1.0_dp - hwl) + &
                (v122*(1.0_dp-wa)+v222*wa)*          hwl ) * wk

    END FUNCTION boxinterp3D

  END SUBROUTINE calc_mod_reflectivity_onsmth


  !================================================================================================
  !+ Module procedure in radar_src for the computation of the radar reflectivity at auxiliary grids
  !------------------------------------------------------------------------------------------------

  SUBROUTINE calc_grd_reflectivity(lcalc)

    !-------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the radar reflectivity at auxiliary grids for each
    !              radar station. Output values are in linear space, not logarithmic (dBZ)!
    !
    ! Method:
    !
    !-------------------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_grd_reflectivity
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_grd_reflectivity'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
      ! calculate reflectivity on the model grid for this radar station:
      IF (ldebug_radsim) WRITE (*,*)  'calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista, ' ...'
      zh_radar = 0.0_dp
      ah_radar = 0.0_dp

      zv_radar   = 0.0_dp
      rrhv_radar = 0.0_dp
      irhv_radar = 0.0_dp
      kdp_radar  = 0.0_dp
      adp_radar  = 0.0_dp
      zvh_radar  = 0.0_dp

      CALL calc_dbz_vec_modelgrid(itlrad_dyn, itlrad_qx, idom, dbz_meta(ista), &
                                  .TRUE., ldebug_radsim, &
                                  TRIM(ydir_mielookup_read), TRIM(ydir_mielookup_write), &
                                  zh_radar=zh_radar, ah_radar=ah_radar, &
                                  zv_radar=zv_radar, rrhv_radar=rrhv_radar, irhv_radar=irhv_radar, &
                                  kdp_radar=kdp_radar, adp_radar=adp_radar, zvh_radar=zvh_radar)
      ! convert dBZ to linear unit mm^6/m^3:
!!$      CALL dbz_to_linear(zh_radar)

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_compgrid)
#endif
#ifdef __ICON__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif
      IF (ldebug_radsim) WRITE (*,*)  'done with calc_dbz_vec_modelgrid() on proc ', my_radar_id, ' for station ', ista

      ! allocate array of reflectivities at auxiliary grids:
      ALLOCATE(rs_grid(ista)%zh_radar_grd(rs_grid(ista)%ngrd))

      CALL interp_model2azislices_scalar (REAL(zh_radar, kind=wp), rs_grid(ista)%ngrd, &
           rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%zh_radar_grd)

      IF (lextdbz .AND. (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
        ALLOCATE(rs_grid(ista)%ah_radar_grd(rs_grid(ista)%ngrd))
        CALL interp_model2azislices_scalar (REAL(ah_radar, kind=wp), rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%ah_radar_grd)
      END IF

      IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
        ALLOCATE(rs_grid(ista)%zv_radar_grd(rs_grid(ista)%ngrd))
        CALL interp_model2azislices_scalar (REAL(zv_radar, kind=wp), rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%zv_radar_grd)
        ALLOCATE(rs_grid(ista)%rrhv_radar_grd(rs_grid(ista)%ngrd))
        CALL interp_model2azislices_scalar (REAL(rrhv_radar, kind=wp), rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%rrhv_radar_grd)
        ALLOCATE(rs_grid(ista)%irhv_radar_grd(rs_grid(ista)%ngrd))
        CALL interp_model2azislices_scalar (REAL(irhv_radar, kind=wp), rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%irhv_radar_grd)
        ALLOCATE(rs_grid(ista)%kdp_radar_grd(rs_grid(ista)%ngrd))
        CALL interp_model2azislices_scalar (REAL(kdp_radar, kind=wp), rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%kdp_radar_grd)
        IF (loutpolall) THEN
          ALLOCATE(rs_grid(ista)%zvh_radar_grd(rs_grid(ista)%ngrd))
          CALL interp_model2azislices_scalar (REAL(zvh_radar, kind=wp), rs_grid(ista)%ngrd, &
               rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%zvh_radar_grd)
        END IF
        IF (lextdbz) THEN
          ALLOCATE(rs_grid(ista)%adp_radar_grd(rs_grid(ista)%ngrd))
          CALL interp_model2azislices_scalar (REAL(adp_radar, kind=wp), rs_grid(ista)%ngrd, &
               rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, rs_grid(ista)%adp_radar_grd)
        END IF
      END IF

      IF (ldebug_radsim .AND. rs_grid(ista)%ngrd > 0) THEN
        WRITE(*,*) 'max rs_grid(ista)%zh_radar_grd = ', &
                   MAXVAL(rs_grid(ista)%zh_radar_grd), &
                   '  ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%zh_radar_grd = ', &
                   MINVAL(rs_grid(ista)%zh_radar_grd), &
                   '  ngrd = ', rs_grid(ista)%ngrd
        IF (lextdbz .AND. &
            (dbz_meta(ista)%itype_refl == 1 .OR. dbz_meta(ista)%itype_refl > 4)) THEN
          WRITE(*,*) 'max rs_grid(ista)%ah_radar_grd = ', &
                     MAXVAL(rs_grid(ista)%ah_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'min rs_grid(ista)%ah_radar_grd = ', &
                     MINVAL(rs_grid(ista)%ah_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
        END IF
        IF ((loutpolstd .OR. loutpolall) .AND. dbz_meta(ista)%itype_refl > 4) THEN
          WRITE(*,*) 'max rs_grid(ista)%zv_radar_grd = ', &
                     MAXVAL(rs_grid(ista)%zv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'min rs_grid(ista)%zv_radar_grd = ', &
                     MINVAL(rs_grid(ista)%zv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'max rs_grid(ista)%rrhv_radar_grd = ', &
                     MAXVAL(rs_grid(ista)%rrhv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'min rs_grid(ista)%rrhv_radar_grd = ', &
                     MINVAL(rs_grid(ista)%rrhv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'max rs_grid(ista)%irhv_radar_grd = ', &
                     MAXVAL(rs_grid(ista)%irhv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'min rs_grid(ista)%irhv_radar_grd = ', &
                     MINVAL(rs_grid(ista)%irhv_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'max rs_grid(ista)%kdp_radar_grd = ', &
                     MAXVAL(rs_grid(ista)%kdp_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          WRITE(*,*) 'min rs_grid(ista)%kdp_radar_grd = ', &
                     MINVAL(rs_grid(ista)%kdp_radar_grd), &
                     '  ngrd = ', rs_grid(ista)%ngrd
          IF (loutpolall) THEN
            WRITE(*,*) 'max rs_grid(ista)%zvh_radar_grd = ', &
                       MAXVAL(rs_grid(ista)%zvh_radar_grd), &
                       '  ngrd = ', rs_grid(ista)%ngrd
            WRITE(*,*) 'min rs_grid(ista)%zvh_radar_grd = ', &
                       MINVAL(rs_grid(ista)%zvh_radar_grd), &
                       '  ngrd = ', rs_grid(ista)%ngrd
          END IF
          IF (lextdbz) THEN
            WRITE(*,*) 'max rs_grid(ista)%adp_radar_grd = ', &
                       MAXVAL(rs_grid(ista)%adp_radar_grd), &
                       '  ngrd = ', rs_grid(ista)%ngrd
            WRITE(*,*) 'min rs_grid(ista)%adp_radar_grd = ', &
                       MINVAL(rs_grid(ista)%adp_radar_grd), &
                       '  ngrd = ', rs_grid(ista)%ngrd
          END IF
        END IF
      END IF

#ifdef __COSMO__
      CALL get_runtime_timings (i_fwo_comppolar)
#endif

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_grd_reflectivity


  !================================================================================================
  !+ Module procedure in radar_src for the computation of the refraction indices at auxiliary grids
  !------------------------------------------------------------------------------------------------

  SUBROUTINE calc_grd_rfridx(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the refraction indices for each
    !              radar station.
    !
    ! Method:
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):

    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER                 :: ista

    REAL    (KIND=dp), ALLOCATABLE       :: &
         t_gp(:), p_gp(:), e_gp(:), &
         rfridx_grd(:), rfridx_low(:), rfridx_up(:), &
         hfl_gp_low(:), hfl_gp_up(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_grd_rfridx
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_grd_rfridx'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! loop over radar stations
    DO ista = 1,nradsta

      ! decide if current radar station shall be output at this time step
      IF (.NOT.lcalc(ista)) CYCLE

      ! allocate array of refraction indices at auxiliary grids :
      ALLOCATE(rs_grid(ista)%rfridx_grd(rs_grid(ista)%ngrd))

      ! allocate helper vectors:
      ALLOCATE( t_gp(rs_grid(ista)%ngrd) )
      ALLOCATE( p_gp(rs_grid(ista)%ngrd) )
      ALLOCATE( e_gp(rs_grid(ista)%ngrd) )
      ALLOCATE( rfridx_grd(rs_grid(ista)%ngrd) )

      ! temperature at interpolated levels
      CALL interp_model2azislices_scalar ( t, rs_grid(ista)%ngrd, &
           rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, t_gp )

      ! pressure at interpolated levels
      CALL interp_model2azislices_scalar ( p, rs_grid(ista)%ngrd, &
           rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, p_gp )

      ! vapor pressure at interpolated levels
      CALL interp_model2azislices_scalar ( vapor_pres, rs_grid(ista)%ngrd, &
           rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, e_gp )

      rfridx_grd = refr_index_air(t_gp, p_gp, e_gp)

      IF (lsode) THEN

        ALLOCATE( rfridx_low(rs_grid(ista)%ngrd) )
        ALLOCATE( rfridx_up (rs_grid(ista)%ngrd) )
        ALLOCATE( hfl_gp_low(rs_grid(ista)%ngrd) )
        ALLOCATE( hfl_gp_up (rs_grid(ista)%ngrd) )

        ! temperature at one level down:
        CALL interp_model2azislices_scalar ( t, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, t_gp, at_k_lower=.TRUE. )

        ! pressure at one level down:
        CALL interp_model2azislices_scalar ( p, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, p_gp, at_k_lower=.TRUE. )

        ! vapor pressure at one level down:
        CALL interp_model2azislices_scalar ( vapor_pres, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, e_gp, at_k_lower=.TRUE. )

        ! height at one level down:
        CALL interp_model2azislices_scalar ( hfl, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, hfl_gp_low, at_k_lower=.TRUE. )

        rfridx_low = refr_index_air(t_gp, p_gp, e_gp)

        ! temperature at one level up:
        CALL interp_model2azislices_scalar ( t, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, t_gp, at_k_upper=.TRUE. )

        ! pressure at one level up:
        CALL interp_model2azislices_scalar ( p, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, p_gp, at_k_upper=.TRUE. )

        ! vapor pressure at one level up:
        CALL interp_model2azislices_scalar ( vapor_pres, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, e_gp, at_k_upper=.TRUE. )

        ! height at one level up:
        CALL interp_model2azislices_scalar ( hfl, rs_grid(ista)%ngrd, &
             rs_grid(ista)%ind_intp, rs_grid(ista)%w_intp, hfl_gp_up, at_k_upper=.TRUE. )

        rfridx_up = refr_index_air(t_gp, p_gp, e_gp)

        ! for SODE: relative vertical reflectivity gradient
        rs_grid(ista)%rfridx_grd = 1.0/rfridx_grd * ((rfridx_low - rfridx_up)/(hfl_gp_low - hfl_gp_up))

        DEALLOCATE( rfridx_low )
        DEALLOCATE( rfridx_up  )
        DEALLOCATE( hfl_gp_low )
        DEALLOCATE( hfl_gp_up  )

      ELSE

        ! for TORE:
        rs_grid(ista)%rfridx_grd = rfridx_grd

      END IF

      IF (ldebug_radsim .AND. rs_grid(ista)%ngrd > 0) THEN
        WRITE(*,*) 'max rs_grid(ista)%rfridx = ', &
                   MAXVAL(rs_grid(ista)%rfridx_grd), ' ngrd = ', rs_grid(ista)%ngrd
        WRITE(*,*) 'min rs_grid(ista)%rfridx = ', &
                   MINVAL(rs_grid(ista)%rfridx_grd), ' ngrd = ', rs_grid(ista)%ngrd
      END IF

      DEALLOCATE(t_gp, p_gp, e_gp, rfridx_grd)

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_grd_rfridx



  !==============================================================================================
  !+ Module procedure in radar_src for the computation of the refraction indices at radar station
  !----------------------------------------------------------------------------------------------

  SUBROUTINE calc_sta_rfridx(lcalc)

    !------------------------------------------------------------------------------
    !
    ! Description: This subroutine calculates the refraction indices for each
    !              radar station.
    !
    ! Method:
    !
    !------------------------------------------------------------------------------
    !

    ! Scalar arguments with intent(in):
    LOGICAL                              :: lcalc(nradsta_max)

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) :: yzroutine
    CHARACTER (LEN=80) :: yzerrmsg
    CHARACTER (LEN=8)  :: ci

    INTEGER            :: ista, k, ko
    LOGICAL            :: found

    REAL    (KIND=dp)                    :: &
         t_sta, p_sta, e_sta,                     &
         rfridx_sta, rfridx_low, rfridx_up,       &
         hfl_low, hfl_up

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE calc_sta_rfridx
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'calc_sta_rfridx'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    DO ista = 1,nradsta
      rs_grid(ista)%rfridx_sta = 999.99_dp
    ENDDO

    ! loop over radar stations
    DO ista = 1,nradsta

      IF (.NOT.lcalc(ista)) CYCLE

      CALL interp3D_model2geo_scalar (idom, t, rs_meta(ista)%lon, rs_meta(ista)%lat, &
                                      rs_meta(ista)%alt_msl, t_sta, k, found)

      IF (.NOT. found) THEN
        ! Station is not located on this PE, so cycle loop over stations
        CYCLE
      END IF

      CALL interp3D_model2geo_scalar (idom, p, rs_meta(ista)%lon, rs_meta(ista)%lat, &
                                      rs_meta(ista)%alt_msl, p_sta, k, found)

      CALL interp3D_model2geo_scalar (idom, vapor_pres, rs_meta(ista)%lon, &
                                      rs_meta(ista)%lat, rs_meta(ista)%alt_msl, &
                                      e_sta, k, found)

      rfridx_sta = refr_index_air(t_sta, p_sta, e_sta)

      IF (lsode) THEN

        ko = one_level_down(k)
        k  = one_level_up(ko)

        CALL interp2D_model2geo_horiz_scalar(idom, t, k, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, t_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, p, k, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, p_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, vapor_pres, k, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, e_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, hfl, k, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, hfl_up, found)
        rfridx_up = refr_index_air(t_sta, p_sta, e_sta)

        CALL interp2D_model2geo_horiz_scalar(idom, t, ko, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, t_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, p, ko, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, p_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, vapor_pres, ko, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, e_sta, found)
        CALL interp2D_model2geo_horiz_scalar(idom, hfl, ko, rs_meta(ista)%lon, &
                                             rs_meta(ista)%lat, hfl_low, found)
        rfridx_low = refr_index_air(t_sta, p_sta, e_sta)

        ! for lsode (SR ..._onlinenew)
        rs_grid(ista)%rfridx_sta =  1.0/rfridx_sta * ((rfridx_low -  rfridx_up)/(hfl_low - hfl_up))

      ELSE

        ! for online
        rs_grid(ista)%rfridx_sta =  rfridx_sta

      ENDIF

      IF (ldebug_radsim) THEN
        WRITE(*,*) 'Radar station ista=', ista,':  rfridx = ',rs_grid(ista)%rfridx_sta
      END IF
      IF (rs_grid(ista)%rfridx_sta < miss_threshold) THEN
        ci(:) = ' '
        WRITE (ci,'(i8)') ista
        CALL abort_run (my_radar_id, 20739, &
             'ERROR: problem in ' //TRIM(yzroutine)// &
             ': rfridx_sta for radar station ' //TRIM(ADJUSTL(ci))//&
             ', station height probably below model surface!', &
             TRIM(yzroutine)//',calculation of refractive index at radar stations')
      END IF

    END DO   ! loop over radar stations

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_sta_rfridx

  !==============================================================================================
  !+ Module procedure in radar_src for the distribution to all PEs of the refraction indices at radar stations
  !----------------------------------------------------------------------------------------------

  SUBROUTINE distribute_sta_rfridx

    !----------------------------------------------------------------------------------
    ! Description: This subroutine distributes the refraction indices at radar stations
    !              to all processors
    !
    !----------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) :: yzroutine
    CHARACTER (LEN=80) :: yzerrmsg
    CHARACTER (LEN=8)  :: ci

    INTEGER            :: ierr,ista,i

    REAL (KIND=dp)     :: rfridx_sta(nradsta), recvbuf(nradsta*num_compute_fwo)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_sta_rfridx
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distribute_sta_rfridx'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    DO ista = 1, nradsta

      rfridx_sta(ista) = rs_grid(ista)%rfridx_sta

    ENDDO

    recvbuf = 999.99_dp

    IF (num_compute_fwo > 1) THEN

      CALL mpi_gather(rfridx_sta,nradsta,imp_fwo_double,recvbuf,nradsta,imp_fwo_double,0,icomm_cart_fwo,ierr)

      IF (my_cart_id_fwo == 0) THEN

        DO i = 1, nradsta*num_compute_fwo

          IF(recvbuf(i) .NE. 999.99_dp) THEN

            IF (MOD(i,nradsta) == 0) THEN

              IF (rs_grid(nradsta)%rfridx_sta < 900.0_dp .AND. &
                   ABS(rs_grid(nradsta)%rfridx_sta-recvbuf(i)) > 1e-12_dp) THEN
                ! In that case, there must be some error in the
                ! assignment of the radar stations to station indices.
                ci(:) = ' '
                WRITE (ci,'(i8)') nradsta
                CALL abort_run (my_radar_id, 20572, &
                     'ERROR: problem in ' //TRIM(yzroutine)// &
                     ': possible problem with assignment of radar station ' //TRIM(ADJUSTL(ci)), &
                     TRIM(yzroutine)//', distribution of refractive index at radar stations')
              END IF

              rs_grid(nradsta)%rfridx_sta = recvbuf(i)

            ELSE

              IF (rs_grid(MOD(i,nradsta))%rfridx_sta < 900.0_dp .AND. &
                   ABS(rs_grid(MOD(i,nradsta))%rfridx_sta-recvbuf(i)) > 1e-12_dp) THEN
                ! In that case, there must be some error in the
                ! assignment of the radar stations to station indices.
                ci(:) = ' '
                WRITE (ci, '(i8)') nradsta
                CALL abort_run (my_radar_id, 20573, &
                     'ERROR: problem in ' //TRIM(yzroutine)// &
                     ': possible problem with assignment of radar station ' //TRIM(ADJUSTL(ci)), &
                     TRIM(yzroutine)//', distribution of refractive index at radar stations')
              END IF

              rs_grid(MOD(i,nradsta))%rfridx_sta = recvbuf(i)

            ENDIF

          ENDIF

        ENDDO

      ENDIF

    ENDIF

    IF(num_compute_fwo > 1) THEN

      IF(my_cart_id_fwo == 0) THEN

        DO ista = 1, nradsta
          rfridx_sta(ista) = rs_grid(ista)%rfridx_sta
        ENDDO

      ENDIF

      CALL distribute_values_radar(rfridx_sta,nradsta,0,icomm_cart_fwo,ierr)

      IF(my_cart_id_fwo /= 0) THEN

        DO ista = 1, nradsta
          rs_grid(ista)%rfridx_sta = rfridx_sta(ista)
        ENDDO

      ENDIF

    ENDIF


  END SUBROUTINE distribute_sta_rfridx


  SUBROUTINE distribute_onlineinfos_reals(rvector_loc,ipos_loc,size_loc, &
       iaz_start,iaz_end,naz_nbl,rvector_az,icomm,npes)

    !-----------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data for a variable (in type real) from all processors and
    !              sends them back to processors according to the azimuths of data FOR ONE RADAR STATION!
    !
    ! Method:      mpi_alltoall and mpi_alltoallv, based on the azimut distribution from para_range_all().
    !
    ! Remarks:     DEPRECATED! is less efficient than the newer equivalent routines
    !              distribute_onlineinfos_reals() and distr_onlinf_all_reals_all2allv().
    !
    !-----------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-----------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER, INTENT(IN) :: size_loc,                         & ! size of vector
                           icomm,                            & ! MPI-communicator
                           npes,                             & ! number of PEs
                           naz_nbl,                          & ! number of azimuts in total azimut-slice-grid (incl. 2*nbl_az)
                           iaz_start,                        & ! the start index of azimuth
                           iaz_end                             ! the final index of azimuth

    INTEGER,INTENT(IN) :: ipos_loc(size_loc)

    REAL (KIND=dp),      INTENT(IN) :: rvector_loc(size_loc)

    REAL (KIND=dp),      POINTER    :: rvector_az(:)                       ! MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(rvector_az)

    INTEGER            :: ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim

    INTEGER            ::   &
                            iisend(npes),iirecv(npes),         &
                            iscnt(0:npes-1),ircnt(0:npes-1),   &
                            isdsp(0:npes-1),irdsp(0:npes-1),   &
                            istart(0:npes-1),iend(0:npes-1)

    REAL (KIND=dp),      ALLOCATABLE:: isend(:), irecv(:)


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_reals
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distribute_onlineinfos_reals'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor (assuming that each radar has
    ! the same number of azimuts:
    CALL para_range(iaz_start,iaz_end,npes,nbl_az,istart,iend)

    IF(ASSOCIATED(rvector_az)) DEALLOCATE(rvector_az)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0


    IF (num_compute_fwo > 1) THEN

      count1 = 0
      count2 = 0
      dim    = 0

      DO irank = 0, npes-1
        ngrd = 0
!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, size_loc
          ! for each grid point determine the continuous number np of grid points
          np = ipos_loc(igrd)
          ! determine indices in azimuth for current grid point in total azimut-slice-grid
          m  = MOD(np-1,naz_nbl) + 1
          ! check if the current azimuth should be sent to irank-processor
          IF (istart(irank) <= m .AND. m <= iend(irank)) THEN
            ngrd = ngrd + 1
          ENDIF
        ENDDO ! loop over ngrd
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)))
      isend = miss_value

      ngrd = 0
      DO irank = 0, npes-1
!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, size_loc
          np = ipos_loc(igrd)
          m  = MOD(np-1,naz_nbl) + 1
          IF (istart(irank) <= m .AND. m <= iend(irank)) THEN
            ngrd = ngrd + 1
            isend(ngrd) = rvector_loc(igrd)
          ENDIF
        ENDDO ! loop over ngrd
      ENDDO ! loop over npes


      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      !distribute the number of data that each processor may obtain
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      ircnt(0:npes-1) = iirecv(1:npes)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)))
      irecv = miss_value

      !distribute the data to processors
      CALL mpi_alltoallv(isend,iscnt,isdsp,imp_fwo_double,irecv,ircnt,irdsp,imp_fwo_double,icomm_cart_fwo,ierr)


      ALLOCATE(rvector_az(SUM(ircnt)))

      rvector_az(:) = irecv

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      ALLOCATE(rvector_az(size_loc))

      rvector_az(:) = rvector_loc(:)

    END IF   ! num_compute_fwo > 1


    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distribute_onlineinfos_reals


  SUBROUTINE distribute_onlineinfos_integers(ivector_loc,ipos_loc,size_loc,iaz_start, &
       iaz_end,naz_nbl,ivector_az,icomm,npes)

    !-------------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data of a variable (in type integer) from all processors and
    !              sends them back to processors according to the azimuths of data  FOR ONE RADAR STATION!
    !
    ! Method:      mpi_alltoall and mpi_alltoallv, based on the azimut distribution from para_range().
    !
    ! Remarks:     DEPRECATED! is less efficient than the newer equivalent routines
    !              distribute_onlineinfos_integers() and distr_onlinf_all_int_all2allv().
    !
    !-------------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-------------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER, INTENT(IN) ::     size_loc,                         & ! size of vector
                               icomm,                            & ! MPI-communicator
                               npes,                             & ! number of PEs
                               naz_nbl,                          & ! number of azimuths in total azimut-slice-grid (incl. 2*nbl_az)
                               iaz_start,                        & ! the start index of azimuth
                               iaz_end                             ! the final index of azimuth

    INTEGER, INTENT(IN) ::     ipos_loc(size_loc)

    INTEGER, INTENT(IN) ::     ivector_loc(size_loc)

    INTEGER,    POINTER ::     ivector_az(:)                       ! MUST BE INITIALIZED BEFORE INPUT! e.g., nullify(ivector_az)

    INTEGER             ::     ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim

    INTEGER             ::     &
                                                 iisend(npes),iirecv(npes),        &
                                                 iscnt(0:npes-1),ircnt(0:npes-1),  &
                                                 isdsp(0:npes-1),irdsp(0:npes-1),  &
                                                 istart(0:npes-1),iend(0:npes-1)

    INTEGER, ALLOCATABLE::     isend(:), irecv(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_integers
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distribute_onlineinfos_integers'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor, istart saves minimum, iend maximum
    CALL para_range(iaz_start,iaz_end,npes,nbl_az,istart,iend)

    IF(ASSOCIATED(ivector_az)) DEALLOCATE(ivector_az)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0

    IF (num_compute_fwo > 1) THEN

      ierr   = 0
      count1 = 0
      count2 = 0

      DO irank = 0, npes-1
        ngrd = 0
!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, size_loc
          ! for each grid point determine the continuous number np of grid points
          np = ipos_loc(igrd)
          ! determine indices in azimuth for current grid point in total azimut-slice-grid:
          m  = MOD(np-1,naz_nbl) + 1
          ! check if the current azimuth should be sent to irank-processor
          IF (istart(irank) <= m .AND. m <= iend(irank)) THEN
            ngrd = ngrd + 1
          ENDIF
        ENDDO ! loop over ngrd
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)))
      isend = missval_int

      ngrd = 0
      DO irank = 0, npes-1
!CDIR NODEP,VOVERTAKE,VOB
        DO igrd = 1, size_loc
          ! for each grid point determine the continuous number np of grid points
          np = ipos_loc(igrd)
          ! determine indices in azimuth for current grid point in total azimut-slice-grid:
          m  = MOD(np-1,naz_nbl) + 1
          ! check if the current azimuth should be sent to irank-processor
          IF (istart(irank) <= m .AND. m <= iend(irank)) THEN
            ngrd = ngrd + 1
            isend(ngrd) = ivector_loc(igrd)
          ENDIF
        ENDDO ! loop over ngrd
      ENDDO ! loop over npes

      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      !distribute the number of data that each processor may obtain
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      IF (ierr /= 0) THEN
        WRITE(*,*) 'Error in MPI_ALLTOALL'
      ENDIF

      ircnt(0:npes-1) = iirecv(1:npes)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)))
      irecv = missval_int

      !distribute the data to processors
      CALL mpi_alltoallv(isend,iscnt,isdsp,imp_fwo_integers,irecv,ircnt,irdsp,imp_fwo_integers,icomm_cart_fwo,ierr)

      ALLOCATE(ivector_az(SUM(ircnt)))

      ivector_az(:) = irecv

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      ALLOCATE(ivector_az(size_loc))

      ivector_az(:) = ivector_loc(:)

    END IF   ! num_compute_fwo > 1

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distribute_onlineinfos_integers

  SUBROUTINE distr_onlinf_all_reals_all2allv (lcalc,rvector_loc,ipos_loc,nsta, &
       iaz_start,iaz_end,naz_nbl,rvector_az,icomm,npes)

    !-----------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data for a variable (in type real) from all processors and
    !              sends them back to processors according to the combined azimut vector of all radar stations.
    !
    ! Method:      mpi_alltoall and mpi_alltoallv, based on the azimut distribution from para_range_all().
    !
    ! Remarks:    - is more efficient than distribute_onlineinfos_reals()
    !             - might or might not be more efficient than distribute_onlineinfos_all_reals().
    !
    !-----------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-----------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32)  :: yzroutine
    CHARACTER (LEN=80)  :: yzerrmsg

    INTEGER, INTENT(IN) :: nsta,                         & ! number of radar stations
                           icomm,                        & ! MPI-communicator
                           npes,                         & ! number of PEs
                           naz_nbl(nsta),                & ! number of azimuts in total azimut-slice-grid (incl. 2*nbl_az)
                           iaz_start(nsta),              & ! the start index of azimuth
                           iaz_end(nsta)                   ! the final index of azimuth

    LOGICAL,      INTENT(in)  :: lcalc      (nsta)
    TYPE(rpvect), INTENT(in)  :: rvector_loc(nsta)
    TYPE(ipvect), INTENT(in)  :: ipos_loc   (nsta)

    TYPE(rpvect), INTENT(out) :: rvector_az(nsta)   ! Each pointer of the vector must point to a nullified pointer of rs_grid- structure!

    INTEGER            :: ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim,dimp(nsta),&
                          ngrd_sta_send_buf

    INTEGER            :: ngrd_sta_send(nsta,npes),ngrd_sta_recv(nsta,npes), &
                            iisend(npes),iirecv(npes),         &
                            iscnt(0:npes-1),ircnt(0:npes-1),   &
                            isdsp(0:npes-1),irdsp(0:npes-1),   &
                            istart(0:npes-1,nsta),iend(0:npes-1,nsta)  ! ,   &

    REAL (KIND=dp),      ALLOCATABLE:: isend(:), irecv(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_reals
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distr_onlinf_all_reals_all2allv'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor (assuming that each radar has
    ! the same number of azimuts:

    CALL para_range_all(iaz_start,iaz_end,nsta,npes,nbl_az,istart,iend)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0

    IF (num_compute_fwo > 1) THEN

      count1 = 0
      count2 = 0
      dim    = 0

      ngrd_sta_send = 0
      ngrd_sta_recv = 0

      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            ngrd_sta_send_buf = 0
!$omp parallel do private (igrd,np,m) reduction(+:ngrd,ngrd_sta_send_buf)
            DO igrd = 1, rvector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd_sta_send_buf = ngrd_sta_send_buf + 1
                ngrd = ngrd + 1
              ENDIF
            ENDDO ! loop over ngrd
!$omp end parallel do
            ngrd_sta_send(ista,irank+1) = ngrd_sta_send_buf
          END IF
        ENDDO
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)))
      isend = miss_value

!$omp parallel do private (irank,ngrd,ista,igrd,np,m)
      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            DO igrd = 1, rvector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd = ngrd + 1
                isend(isdsp(irank)+ngrd) = rvector_loc(ista)%p(igrd)
              ENDIF
            ENDDO ! loop over ngrd
          END IF
        ENDDO
      ENDDO ! loop over npes
!$omp end parallel do

      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      ! distribute the number of data that each processor may obtain:
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      ircnt(0:npes-1) = iirecv(1:npes)

      ! distribute the number of data per station that each processor may obtain:
      CALL mpi_alltoall(ngrd_sta_send,nsta,imp_fwo_integers,ngrd_sta_recv,nsta,imp_fwo_integers,icomm_cart_fwo,ierr)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)))
      irecv = miss_value

      !distribute the data to processors
      CALL mpi_alltoallv(isend,iscnt,isdsp,imp_fwo_double,irecv,ircnt,irdsp,imp_fwo_double,icomm_cart_fwo,ierr)

      DO ista = 1, nsta
        ngrd = SUM(ngrd_sta_recv(ista,1:npes))
        IF ( ngrd > 0 ) THEN
          ALLOCATE(rvector_az(ista)%p(ngrd))
!$omp parallel do
          DO igrd = 1, ngrd
            rvector_az(ista)%p(igrd) = -987.65_dp
          END DO
!$omp end parallel do
          rvector_az(ista)%n = ngrd
        ELSE
          ALLOCATE(rvector_az(ista)%p(0))
          rvector_az(ista)%n = 0
        END IF
      END DO

      dim = 0
      dimp = 0
      DO irank = 0, npes-1
        DO ista = 1, nsta
          ngrd = ngrd_sta_recv(ista,irank+1)
          IF (ngrd > 0) THEN
!$omp parallel do
            DO igrd = 1, ngrd
              rvector_az(ista)%p(dimp(ista)+igrd) = irecv(dim+igrd)
            END DO
!$omp end parallel do
            dim  = dim  + ngrd
            dimp(ista) = dimp(ista) + ngrd
          END IF
        END DO
      END DO

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      DO ista = 1, nsta
        IF (lcalc(ista)) THEN
          ALLOCATE(rvector_az(ista)%p(rvector_loc(ista)%n))
          rvector_az(ista)%p = rvector_loc(ista)%p
          rvector_az(ista)%n = rvector_loc(ista)%n
        END IF
      END DO

    END IF   ! num_compute_fwo > 1

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distr_onlinf_all_reals_all2allv

  SUBROUTINE distribute_onlineinfos_all_reals(lcalc,rvector_loc,ipos_loc,nsta, &
       iaz_start,iaz_end,naz_nbl,rvector_az,icomm,npes)

    !-----------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data for a variable (in type real) from all processors and
    !              sends them back to processors according to the combined azimut vector of all radar stations.
    !
    ! Method:      mpi_alltoall and non-blocking MPI_ISEND / MPI_IRECV,
    !              based on the azimut distribution from para_range_all().
    !
    ! Remarks:    - is more efficient than distribute_onlineinfos_reals()
    !             - might or might not be more efficient than distr_onlinf_all_reals_all2allv().
    !
    !-----------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-----------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32)  :: yzroutine
    CHARACTER (LEN=80)  :: yzerrmsg

    INTEGER, INTENT(IN) :: nsta,                         & ! number of radar stations
                           icomm,                        & ! MPI-communicator
                           npes,                         & ! number of PEs
                           naz_nbl(nsta),                & ! number of azimuts in total azimut-slice-grid (incl. 2*nbl_az)
                           iaz_start(nsta),              & ! the start index of azimuth
                           iaz_end(nsta)                   ! the final index of azimuth

    LOGICAL,      INTENT(in)  :: lcalc      (nsta)
    TYPE(rpvect), INTENT(in)  :: rvector_loc(nsta)
    TYPE(ipvect), INTENT(in)  :: ipos_loc   (nsta)

    TYPE(rpvect), INTENT(out) :: rvector_az(nsta)   ! Each pointer of the vector must point to a nullified pointer of rs_grid- structure!

    INTEGER            :: ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim,dimp(nsta),&
                          irankr, iranks, ngrd_sta_send_buf

    INTEGER            :: ngrd_sta_send(nsta,npes),ngrd_sta_recv(nsta,npes), &
                          iisend(npes),iirecv(npes),         &
                          iscnt(0:npes-1),ircnt(0:npes-1),   &
                          isdsp(0:npes-1),irdsp(0:npes-1),   &
                          istart(0:npes-1,nsta),iend(0:npes-1,nsta), &
                          itag(0:npes-1,0:npes-1),           &
                          isrequest(0:npes-1), irrequest(0:npes-1), &
                          istatus(MPI_STATUS_SIZE)

    REAL (KIND=dp),      ALLOCATABLE:: isend(:), irecv(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_reals
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distribute_onlineinfos_all_reals'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor (assuming that each radar has
    ! the same number of azimuts:

    CALL para_range_all(iaz_start,iaz_end,nsta,npes,nbl_az,istart,iend)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0

    IF (num_compute_fwo > 1) THEN

      count1 = 0
      count2 = 0
      dim    = 0

      ngrd_sta_send = 0
      ngrd_sta_recv = 0
      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            ngrd_sta_send_buf = 0
!$omp parallel do private (igrd,np,m) reduction(+:ngrd,ngrd_sta_send_buf)
            DO igrd = 1, rvector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd_sta_send_buf = ngrd_sta_send_buf + 1
                ngrd = ngrd + 1
              ENDIF
            ENDDO ! loop over ngrd
!$omp end parallel do
            ngrd_sta_send(ista,irank+1) = ngrd_sta_send_buf
          END IF
        ENDDO
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)))
      isend = miss_value

!$omp parallel do private (irank,ngrd,ista,igrd,np,m)
      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            DO igrd = 1, rvector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd = ngrd + 1
                isend(isdsp(irank)+ngrd) = rvector_loc(ista)%p(igrd)
              ENDIF
            ENDDO ! loop over ngrd
          END IF
        ENDDO
      ENDDO ! loop over npes
!$omp end parallel do

      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      ! distribute the number of data that each processor may obtain:
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      ircnt(0:npes-1) = iirecv(1:npes)

      ! distribute the number of data per station that each processor may obtain:
      CALL mpi_alltoall(ngrd_sta_send,nsta,imp_fwo_integers,ngrd_sta_recv,nsta,imp_fwo_integers,icomm_cart_fwo,ierr)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)))
      irecv = miss_value

      !distribute the data to processors
      !=================================

      ! Create tags for the non-blocking communications:
!$omp parallel do private (irankr,iranks)
      DO irankr=0, npes-1
        DO iranks=0, npes-1
          itag(iranks,irankr) = iranks + irankr*npes
        END DO
      END DO
!$omp end parallel do

      ! Handle the communication from my_cart_id_fwo to itself by a simple copy, not MPI:
      IF (ircnt(my_cart_id_fwo) > 0 ) THEN
        irecv( irdsp(my_cart_id_fwo)+1 : irdsp(my_cart_id_fwo)+ircnt(my_cart_id_fwo) ) = &
             isend( isdsp(my_cart_id_fwo)+1 : isdsp(my_cart_id_fwo)+iscnt(my_cart_id_fwo) )
      END IF

      ! Open non-blocking receive channels: loop through the list of possible send nodes from which to receive data:
      DO iranks=0, npes-1
        IF (ircnt(iranks) > 0 .AND. iranks /= my_cart_id_fwo) THEN
          CALL MPI_IRECV(irecv(irdsp(iranks)+1:irdsp(iranks)+ircnt(iranks)), ircnt(iranks),imp_fwo_double, &
               iranks, itag(iranks,my_cart_id_fwo), icomm_cart_fwo, irrequest(iranks), ierr)
        END IF
      END DO

      ! Open non-blocking send channels: loop through the list of possible receive nodes of data packages from this node:
      DO irankr=0, npes-1
        IF (iscnt(irankr) > 0 .AND. irankr /= my_cart_id_fwo) THEN
          CALL MPI_ISEND(isend(isdsp(irankr)+1:isdsp(irankr)+iscnt(irankr)), iscnt(irankr),imp_fwo_double, &
               irankr, itag(my_cart_id_fwo,irankr), icomm_cart_fwo, isrequest(irankr), ierr)
        END IF
      END DO

      ! Wait for completion of the IREC communication:
      DO irankr=0, npes-1
        IF (iscnt(irankr) > 0 .AND. irankr /= my_cart_id_fwo) THEN
          CALL MPI_WAIT(isrequest(irankr), istatus, ierr)
        END IF
      END DO
      DO iranks=0, npes-1
        IF (ircnt(iranks) > 0 .AND. iranks /= my_cart_id_fwo) THEN
          CALL MPI_WAIT(irrequest(iranks), istatus, ierr)
        END IF
      END DO

      ! Re-distribute the received data to the correct stations:
      !=========================================================

      DO ista = 1, nsta
        ngrd = SUM(ngrd_sta_recv(ista,1:npes))
        IF ( ngrd > 0 ) THEN
          ALLOCATE(rvector_az(ista)%p(ngrd))
!$omp parallel do
          DO igrd = 1, ngrd
            rvector_az(ista)%p(igrd) = -987.65_dp
          END DO
!$omp end parallel do
          rvector_az(ista)%n = ngrd
        ELSE
          ALLOCATE(rvector_az(ista)%p(0))
          rvector_az(ista)%n = 0
        END IF
      END DO

      dim = 0
      dimp = 0
      DO irank = 0, npes-1
        DO ista = 1, nsta
          ngrd = ngrd_sta_recv(ista,irank+1)
          IF (ngrd > 0) THEN
!$omp parallel do
            DO igrd = 1, ngrd
              rvector_az(ista)%p(dimp(ista)+igrd) = irecv(dim+igrd)
            END DO
!$omp end parallel do
            dim  = dim  + ngrd
            dimp(ista) = dimp(ista) + ngrd
          END IF
        END DO
      END DO

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      DO ista = 1, nsta
        IF (lcalc(ista)) THEN
          ALLOCATE(rvector_az(ista)%p(rvector_loc(ista)%n))
          rvector_az(ista)%p = rvector_loc(ista)%p
          rvector_az(ista)%n = rvector_loc(ista)%n
        END IF
      END DO

    END IF   ! num_compute_fwo > 1

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distribute_onlineinfos_all_reals


  SUBROUTINE distr_onlinf_all_int_all2allv(lcalc,ivector_loc,ipos_loc,nsta, &
       iaz_start,iaz_end,naz_nbl,ivector_az,icomm,npes)

    !-----------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data for a variable (in type integer) from all processors and
    !              sends them back to processors according to the combined azimut vector of all radar stations.
    !
    ! Method:      mpi_alltoall and mpi_alltoallv, based on the azimut distribution from para_range_all().
    !
    ! Remarks:    - is more efficient than distribute_onlineinfos_integers()
    !             - might or might not be more efficient than distribute_onlineinfos_all_int().
    !
    !-----------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-----------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32)  :: yzroutine
    CHARACTER (LEN=80)  :: yzerrmsg

    INTEGER, INTENT(IN) :: nsta,                     & ! number of radar stations
                           icomm,                    & ! MPI-communicator
                           npes,                     & ! number of PEs
                           naz_nbl(nsta),            & ! number of azimuts in total azimut-slice-grid (incl. 2*nbl_az)
                           iaz_start(nsta),          & ! the start index of azimuth
                           iaz_end(nsta)               ! the final index of azimuth

    LOGICAL,      INTENT(in)  :: lcalc      (nsta)
    TYPE(ipvect), INTENT(in)  :: ivector_loc(nsta)
    TYPE(ipvect), INTENT(in)  :: ipos_loc   (nsta)

    TYPE(ipvect), INTENT(out) :: ivector_az(nsta)   ! Each pointer of the vector must point to a nullified pointer of rs_grid- structure!

    INTEGER            :: ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim,dimp(nsta),&
                          ngrd_sta_send_buf

    INTEGER            :: ngrd_sta_send(nsta,npes),ngrd_sta_recv(nsta,npes), &
                          iisend(npes),iirecv(npes),         &
                          iscnt(0:npes-1),ircnt(0:npes-1),   &
                          isdsp(0:npes-1),irdsp(0:npes-1),   &
                          istart(0:npes-1,nsta),iend(0:npes-1,nsta)

    INTEGER, ALLOCATABLE :: isend(:), irecv(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_reals
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distr_onlinf_all_int_all2allv'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor (assuming that each radar has
    ! the same number of azimuts:

    CALL para_range_all(iaz_start,iaz_end,nsta,npes,nbl_az,istart,iend)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0

    IF (num_compute_fwo > 1) THEN

      count1 = 0
      count2 = 0
      dim    = 0

      ngrd_sta_send = 0
      ngrd_sta_recv = 0

      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            ngrd_sta_send_buf = 0
!$omp parallel do private (igrd,np,m) reduction(+:ngrd,ngrd_sta_send_buf)
            DO igrd = 1, ivector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd_sta_send_buf = ngrd_sta_send_buf + 1
                ngrd = ngrd + 1
              ENDIF
            ENDDO ! loop over ngrd
!$omp end parallel do
            ngrd_sta_send(ista,irank+1) = ngrd_sta_send_buf
          END IF
        ENDDO
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)),stat=ierr)
      IF (ierr /= 0) THEN
        WRITE (*,*) TRIM(yzroutine)//': Error allocating isend (',SUM(iscnt),')'
        STOP
      END IF
      isend = missval_int

!$omp parallel do private (irank,ngrd,ista,igrd,np,m)
      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            DO igrd = 1, ivector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd = ngrd + 1
                isend(isdsp(irank)+ngrd) = ivector_loc(ista)%p(igrd)
              ENDIF
            ENDDO ! loop over ngrd
          END IF
        ENDDO
      ENDDO ! loop over npes
!$omp end parallel do

      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      ! distribute the number of data that each processor may obtain:
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      ircnt(0:npes-1) = iirecv(1:npes)

      ! distribute the number of data per station that each processor may obtain:
      CALL mpi_alltoall(ngrd_sta_send,nsta,imp_fwo_integers,ngrd_sta_recv,nsta,imp_fwo_integers,icomm_cart_fwo,ierr)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)),stat=ierr)
      IF (ierr /= 0) THEN
        WRITE (*,*) TRIM(yzroutine)//': Error allocating irecv (',SUM(ircnt),')'
        STOP
      END IF
      irecv = missval_int

      !distribute the data to processors
      CALL mpi_alltoallv(isend,iscnt,isdsp,imp_fwo_integers,irecv,ircnt,irdsp,imp_fwo_integers,icomm_cart_fwo,ierr)

      DO ista = 1, nsta
        ngrd = SUM(ngrd_sta_recv(ista,1:npes))
        IF ( ngrd > 0 ) THEN
          ALLOCATE(ivector_az(ista)%p(ngrd),stat=ierr)
          IF (ierr /= 0) THEN
            WRITE (*,*) TRIM(yzroutine)//': Error allocating ivector_az(',ista,')%p(',ngrd,')'
            STOP
          END IF
!$omp parallel do
          DO igrd = 1, ngrd
            ivector_az(ista)%p(igrd) = -987
          END DO
!$omp end parallel do
          ivector_az(ista)%n = ngrd
        ELSE
          ALLOCATE(ivector_az(ista)%p(0),stat=ierr)
          IF (ierr /= 0) THEN
            WRITE (*,*) TRIM(yzroutine)//': Error allocating ivector_az(',ista,')%p(',ngrd,')'
            STOP
          END IF
          ivector_az(ista)%n = 0
        END IF
      END DO

      dim = 0
      dimp = 0
      DO irank = 0, npes-1
        DO ista = 1, nsta
          ngrd = ngrd_sta_recv(ista,irank+1)
          IF (ngrd > 0) THEN
!$omp parallel do
            DO igrd = 1, ngrd
              ivector_az(ista)%p(dimp(ista)+igrd) = irecv(dim+igrd)
            END DO
!$omp end parallel do
            dim  = dim  + ngrd
            dimp(ista) = dimp(ista) + ngrd
          END IF
        END DO
      END DO

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      DO ista = 1, nsta
        IF (lcalc(ista)) THEN
          ALLOCATE(ivector_az(ista)%p(ivector_loc(ista)%n))
          ivector_az(ista)%p = ivector_loc(ista)%p
          ivector_az(ista)%n = ivector_loc(ista)%n
        END IF
      END DO

    END IF   ! num_compute_fwo > 1

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distr_onlinf_all_int_all2allv

  !------------------------------------------------------------------------------

  SUBROUTINE distribute_onlineinfos_all_int(lcalc,ivector_loc,ipos_loc,nsta, &
       iaz_start,iaz_end,naz_nbl,ivector_az,icomm,npes)

    !-----------------------------------------------------------------------------------------------------
    !
    ! Description: This subroutine collects the data for a variable (in type integer) from all processors and
    !              sends them back to processors according to the combined azimut vector of all radar stations.
    !
    ! Method:      mpi_alltoall and non-blocking MPI_ISEND / MPI_IRECV,
    !              based on the azimut distribution from para_range_all().
    !
    ! Remarks:    - is more efficient than distribute_onlineinfos_integers()
    !             - might or might not be more efficient than distr_onlinf_all_int_all2allv().
    !
    !-----------------------------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !-----------------------------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32)  :: yzroutine
    CHARACTER (LEN=80)  :: yzerrmsg

    INTEGER, INTENT(IN) :: nsta,                     & ! number of radar stations
                           icomm,                    & ! MPI-communicator
                           npes,                     & ! number of PEs
                           naz_nbl(nsta),            & ! number of azimuts in total azimut-slice-grid (incl. 2*nbl_az)
                           iaz_start(nsta),          & ! the start index of azimuth
                           iaz_end(nsta)               ! the final index of azimuth

    LOGICAL,      INTENT(in)  :: lcalc      (nsta)
    TYPE(ipvect), INTENT(in)  :: ivector_loc(nsta)
    TYPE(ipvect), INTENT(in)  :: ipos_loc   (nsta)

    TYPE(ipvect), INTENT(out) :: ivector_az(nsta)   ! Each pointer of the vector must point to a nullified pointer of rs_grid- structure!

    INTEGER            :: ierr,ista,irank,igrd,ngrd,np,m,count1,count2,dim,dimp(nsta),&
                            irankr, iranks, ngrd_sta_send_buf

    INTEGER            :: ngrd_sta_send(nsta,npes),ngrd_sta_recv(nsta,npes), &
                            iisend(npes),iirecv(npes),         &
                            iscnt(0:npes-1),ircnt(0:npes-1),   &
                            isdsp(0:npes-1),irdsp(0:npes-1),   &
                            istart(0:npes-1,nsta),iend(0:npes-1,nsta), &
                            itag(0:npes-1,0:npes-1),           &
                            isrequest(0:npes-1), irrequest(0:npes-1), &
                            istatus(MPI_STATUS_SIZE)

    INTEGER, ALLOCATABLE :: isend(:), irecv(:)

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE distribute_onlineinfos_reals
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'distribute_onlineinfos_all_int'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! to estimate the interval of azimuths for each processor:

    CALL para_range_all(iaz_start,iaz_end,nsta,npes,nbl_az,istart,iend)

    iscnt(0:npes-1)= 0
    ircnt(0:npes-1)= 0
    isdsp(0:npes-1)= 0
    irdsp(0:npes-1)= 0

    IF (num_compute_fwo > 1) THEN

      count1 = 0
      count2 = 0
      dim    = 0

      ngrd_sta_send = 0
      ngrd_sta_recv = 0

      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            ngrd_sta_send_buf = 0
!$omp parallel do private (igrd,np,m) reduction(+:ngrd,ngrd_sta_send_buf)
            DO igrd = 1, ivector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd_sta_send_buf = ngrd_sta_send_buf + 1
                ngrd = ngrd + 1
              ENDIF
            ENDDO ! loop over ngrd
!$omp end parallel do
            ngrd_sta_send(ista,irank+1) =  ngrd_sta_send_buf
          END IF
        ENDDO
        iscnt(irank) = ngrd
        isdsp(irank) = count1
        count1 = count1 + iscnt(irank)
      ENDDO ! loop over npes

      ALLOCATE(isend(SUM(iscnt)),stat=ierr)
      IF (ierr /= 0) THEN
        WRITE (*,*) TRIM(yzroutine)//': Error allocating isend (',SUM(iscnt),')'
        STOP
      END IF
      isend = missval_int

!$omp parallel do private (irank,ngrd,ista,igrd,np,m)
      DO irank = 0, npes-1
        ngrd = 0
        DO ista = 1, nsta
          IF (lcalc(ista)) THEN
            DO igrd = 1, ivector_loc(ista)%n
              ! for each grid point determine the continuous number np of grid points
              np = ipos_loc(ista)%p(igrd)
              ! determine indices in azimuth for current grid point in total azimut-slice-grid
              m  = MOD(np-1,naz_nbl(ista)) + 1
              ! check if the current azimuth should be sent to irank-processor
              IF (istart(irank,ista) <= m .AND. m <= iend(irank,ista)) THEN
                ngrd = ngrd + 1
                isend(isdsp(irank)+ngrd) = ivector_loc(ista)%p(igrd)
              ENDIF
            ENDDO ! loop over ngrd
          END IF
        ENDDO
      ENDDO ! loop over npes
!$omp end parallel do

      ! to estimate ircnt
      iisend = iscnt
      iirecv = 0

      ! distribute the number of data that each processor may obtain:
      CALL mpi_alltoall(iisend,1,imp_fwo_integers,iirecv,1,imp_fwo_integers,icomm_cart_fwo,ierr)

      ircnt(0:npes-1) = iirecv(1:npes)

      ! distribute the number of data per station that each processor may obtain:
      CALL mpi_alltoall(ngrd_sta_send,nsta,imp_fwo_integers,ngrd_sta_recv,nsta,imp_fwo_integers,icomm_cart_fwo,ierr)

      ! to estimate irdsp
      DO irank = 0, npes-1

        irdsp(irank) = count2
        count2 = count2 + ircnt(irank)

      ENDDO

      ALLOCATE(irecv(SUM(ircnt)),stat=ierr)
      IF (ierr /= 0) THEN
        WRITE (*,*) TRIM(yzroutine)//': Error allocating irecv (',SUM(ircnt),')'
        STOP
      END IF
      irecv = missval_int


      !distribute the data to processors
      !=================================

      ! Create tags for the non-blocking communications:
!$omp parallel do private (irankr,iranks)
      DO irankr=0, npes-1
        DO iranks=0, npes-1
          itag(iranks,irankr) = iranks + irankr*npes
        END DO
      END DO
!$omp end parallel do

      ! Handle the communication from my_cart_id_fwo to itself by a simple copy, not MPI:
      IF (ircnt(my_cart_id_fwo) > 0 ) THEN
        irecv( irdsp(my_cart_id_fwo)+1 : irdsp(my_cart_id_fwo)+ircnt(my_cart_id_fwo) ) = &
             isend( isdsp(my_cart_id_fwo)+1 : isdsp(my_cart_id_fwo)+iscnt(my_cart_id_fwo) )
      END IF

      ! Open non-blocking receive channels: loop through the list of possible send nodes from which to receive data:
      DO iranks=0, npes-1
        IF (ircnt(iranks) > 0 .AND. iranks /= my_cart_id_fwo) THEN
          CALL MPI_IRECV(irecv(irdsp(iranks)+1:irdsp(iranks)+ircnt(iranks)), ircnt(iranks),imp_fwo_integers, &
               iranks, itag(iranks,my_cart_id_fwo), icomm_cart_fwo, irrequest(iranks), ierr)
        END IF
      END DO

      ! Open non-blocking send channels: loop through the list of possible receive nodes of data packages from this node:
      DO irankr=0, npes-1
        IF (iscnt(irankr) > 0 .AND. my_cart_id_fwo /= irankr) THEN
          CALL MPI_ISEND(isend(isdsp(irankr)+1:isdsp(irankr)+iscnt(irankr)), iscnt(irankr),imp_fwo_integers, &
               irankr, itag(my_cart_id_fwo,irankr), icomm_cart_fwo, isrequest(irankr), ierr)
        END IF
      END DO

      ! Wait for completion of the IREC communication:
      DO irankr=0, npes-1
        IF (iscnt(irankr) > 0 .AND. irankr /= my_cart_id_fwo) THEN
          CALL MPI_WAIT(isrequest(irankr), istatus, ierr)
        END IF
      END DO
      DO iranks=0, npes-1
        IF (ircnt(iranks) > 0 .AND. iranks /= my_cart_id_fwo) THEN
          CALL MPI_WAIT(irrequest(iranks), istatus, ierr)
        END IF
      END DO

      ! Re-distribute the received data to the correct stations:
      !=========================================================

      DO ista = 1, nsta
        ngrd = SUM(ngrd_sta_recv(ista,1:npes))
        IF ( ngrd > 0 ) THEN
          ALLOCATE(ivector_az(ista)%p(ngrd),stat=ierr)
          IF (ierr /= 0) THEN
            WRITE (*,*) TRIM(yzroutine)//': Error allocating ivector_az(',ista,')%p(',ngrd,')'
            STOP
          END IF
!$omp parallel do
          DO igrd = 1, ngrd
            ivector_az(ista)%p(igrd) = -987
          END DO
!$omp end parallel do
          ivector_az(ista)%n = ngrd
        ELSE
          ALLOCATE(ivector_az(ista)%p(0),stat=ierr)
          IF (ierr /= 0) THEN
            WRITE (*,*) TRIM(yzroutine)//': Error allocating ivector_az(',ista,')%p(',ngrd,')'
            STOP
          END IF
          ivector_az(ista)%n = 0
        END IF
      END DO

      dim = 0
      dimp = 0
      DO irank = 0, npes-1
        DO ista = 1, nsta
          ngrd = ngrd_sta_recv(ista,irank+1)
          IF (ngrd > 0) THEN
!$omp parallel do
            DO igrd = 1, ngrd
              ivector_az(ista)%p(dimp(ista)+igrd) = irecv(dim+igrd)
            END DO
!$omp end parallel do
            dim  = dim  + ngrd
            dimp(ista) = dimp(ista) + ngrd
          END IF
        END DO
      END DO

      DEALLOCATE(isend)
      DEALLOCATE(irecv)

    ELSE     ! num_compute_fwo == 1

      DO ista = 1, nsta
        IF (lcalc(ista)) THEN
          ALLOCATE(ivector_az(ista)%p(ivector_loc(ista)%n))
          ivector_az(ista)%p = ivector_loc(ista)%p
          ivector_az(ista)%n = ivector_loc(ista)%n
        END IF
      END DO

    END IF   ! num_compute_fwo > 1

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE distribute_onlineinfos_all_int


  !==============================================================================
  !+ Module procedure in radar_src for the coordinate transformation of radar
  !  radar coordinates to geographical coordinates for constant beam propagation
  !------------------------------------------------------------------------------

  SUBROUTINE rad2geo_const (rsm, lat, lon, h)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation of radar coordinates
    !              (range, azimuth and elevation) to
    !              geographical coordinates (lat, lon and height asl)
    !
    ! Method:      Assumption: 4/3 earth model
    !              Formulas with reference to  Zeng et al. (2014) for 4/3 earth
    !              radius model and Appendix D of Dissertation of Ulrich Blahak
    !              for the geographic coordinates
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    TYPE(radar_meta_type), INTENT(in) :: rsm
    REAL (KIND=dp), ALLOCATABLE, INTENT (inout)     ::        &
         lat(:,:,:),& ! geographical latitude vector       (naz,nra,nel)
         lon(:,:,:),& ! geographical longitude vector      (naz,nra,nel)
         h  (:,:)     ! height above mean sea level vector (nra,nel)  or  (naz,nra) for PRECIP scans with nel=1

    !
    ! Local scalars:
    REAL (KIND=dp)      ::        &
         latsta,  & ! geographical latitude of radar station
         lonsta,  & ! geographical longitude of radar station
         altsta     ! altitude of radar station asl (m)


    REAL (KIND=dp), ALLOCATABLE      ::        &
         ra(:),     & ! array of ranges
         az(:),     & ! array of azimuths
         el(:,:),   & ! array of elevations (nra,nel)  or  (naz,nra) for PRECIP scans with nel=1
         el_precip(:) ! vector of elevations for precip scan

    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine

    INTEGER           :: ira,iaz,iel
    REAL (KIND=dp)    :: h1,s1,re


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE rad2geo_const
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'rad2geo_const'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! allocate and fill up range and azimuth arrays
    CALL get_rangevec (rsm, ra)
    CALL get_azvec    (rsm, az)

    IF (TRIM(rsm%scanname) == 'PRECIP') THEN

      ALLOCATE(el(rsm%naz,1))
      CALL get_elarr_precipscan (rsm, el_precip)
      el(:,1) = el_precip

      IF (.NOT.ALLOCATED(lat)) ALLOCATE( lat(rsm%naz,rsm%nra,rsm%nel))
      IF (.NOT.ALLOCATED(lon)) ALLOCATE( lon(rsm%naz,rsm%nra,rsm%nel))
      IF (.NOT.ALLOCATED(h  )) ALLOCATE( h  (rsm%naz,rsm%nra))

      latsta = rsm%lat
      lonsta = rsm%lon
      altsta = rsm%alt_msl

      re = 4.0_dp/3.0_dp * r_earth_dp

!$omp parallel private(iel,ira,iaz,h1,s1)
      DO iel = 1,rsm%nel
!$omp do
        DO ira = 1,rsm%nra
          DO iaz = 1,rsm%naz

            ! calculate height of radar beam above MSL:
            h1 = SQRT((re+altsta)**2 + ra(ira)**2 + 2*(re+altsta)*ra(ira)*SIN(el(iaz,iel)*degrad)) - re
            ! calculate length of circle at height of MSL:
            s1 = re*ASIN(ra(ira)*COS(el(iaz,iel)*degrad)/(re+h1))

            ! calculate height above mean sea level (only dependent on range and elevation!)
            h(iaz,ira) = h1

            CALL polar2geo_old (lonsta, latsta, r_earth_dp, s1, az(iaz), lon(iaz,ira,iel), lat(iaz,ira,iel))

          END DO
        END DO
!$omp end do
      END DO
!$omp end parallel

    ELSE

      ALLOCATE(el(rsm%naz,rsm%nel))
      DO iel=1, rsm%nel
        el(:,iel) = rsm%el_arr(iel)
      END DO

      IF (.NOT.ALLOCATED(lat)) ALLOCATE( lat(rsm%naz,rsm%nra,rsm%nel))
      IF (.NOT.ALLOCATED(lon)) ALLOCATE( lon(rsm%naz,rsm%nra,rsm%nel))
      IF (.NOT.ALLOCATED(h  )) ALLOCATE( h  (rsm%nra,rsm%nel))

      latsta = rsm%lat
      lonsta = rsm%lon
      altsta = rsm%alt_msl

      re = 4.0_dp/3.0_dp * r_earth_dp

!$omp parallel do private(iel,ira,iaz,h1,s1) collapse(3)
      DO iel = 1,rsm%nel
        DO ira = 1,rsm%nra
          DO iaz = 1,rsm%naz

            ! calculate height of radar beam above MSL:
            h1 = SQRT((re+altsta)**2 + ra(ira)**2 + 2*(re+altsta)*ra(ira)*SIN(el(iaz,iel)*degrad)) - re
            ! calculate length of circle at height of MSL:
            s1 = re*ASIN(ra(ira)*COS(el(iaz,iel)*degrad)/(re+h1))

            ! calculate height above mean sea level (only dependent on range and elevation!)
            h(ira,iel) = h1

            CALL polar2geo_old (lonsta, latsta, r_earth_dp, s1, az(iaz), lon(iaz,ira,iel), lat(iaz,ira,iel))

          END DO
        END DO
      END DO
!$omp end parallel do

    END IF

    ! clean up:
    DEALLOCATE (az, ra, el)
    IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE rad2geo_const

  !==============================================================================
  !+ Module procedure in radar_src for the coordinate transformation of radar
  !  radar coordinates to geographical coordinates for constant beam propagation
  !------------------------------------------------------------------------------

  SUBROUTINE rad2geo_const_vec (rsm,lat,lon,h)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation of radar coordinates
    !              (range, azimuth and elevation) to
    !              geographical coordinates (lat, lon and height asl)
    !
    ! Method:      Assumption: 4/3 earth model
    !              Formulas with reference to  Appendices B and D of Dissertation
    !              of Uli Blahak
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    TYPE(radar_meta_type), INTENT(in) :: rsm
    REAL (KIND=dp), ALLOCATABLE, INTENT (inout)     ::        &
         lat(:),& ! geographical latitude vector (nra*naz*nel)
         lon(:),& ! geographical longitude vector (nra*naz*nel)
         h  (:)   ! height above mean sea level vector (nra*naz*nel)
    !
    ! Local scalars:
    REAL (KIND=dp)      ::        &
         latsta,  & ! geographical latitude of radar station
         lonsta,  & ! geographical longitude of radar station
         altsta     ! altitude of radar station asl (m)

    REAL (KIND=dp), ALLOCATABLE      ::        &
         ra(:),   & ! array of ranges
         az(:),   & ! array of azimuths
         el(:,:), & ! array of elevations
         el_precip(:)  ! vector of elevations for precip scan

    CHARACTER (LEN=32) yzroutine

    INTEGER           :: ira,iaz,iel,irp
    REAL (KIND=dp)    :: h1,s1,re


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE rad2geo_const_vec
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'rad2geo_const_vec'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! allocate and fill up range and azimuth arrays
    CALL get_rangevec (rsm, ra)
    CALL get_azvec    (rsm, az)

    ALLOCATE(el(rsm%naz,rsm%nel))
    IF (TRIM(rsm%scanname) == 'PRECIP') THEN
      CALL get_elarr_precipscan (rsm, el_precip)
      el(:,1) = el_precip
    ELSE
      DO iel=1, rsm%nel
        el(:,iel) = rsm%el_arr(iel)
      END DO
    END IF

    IF (.NOT.ALLOCATED(lat)) ALLOCATE( lat(rsm%naz*rsm%nra*rsm%nel))
    IF (.NOT.ALLOCATED(lon)) ALLOCATE( lon(rsm%naz*rsm%nra*rsm%nel))
    IF (.NOT.ALLOCATED(h  )) ALLOCATE( h  (rsm%naz*rsm%nra*rsm%nel))

    latsta = rsm%lat
    lonsta = rsm%lon
    altsta = rsm%alt_msl

    re = 4.0_dp/3.0_dp * r_earth_dp

!$omp parallel do private(iel,ira,iaz,h1,s1,irp) collapse(3)
    DO iel = 1,rsm%nel
      DO ira = 1,rsm%nra
        DO iaz = 1,rsm%naz

          ! Determine the index of the radar point, based on az,ra,el
          CALL sub2ind3D(iaz,ira,iel, rsm%naz,rsm%nra, irp)

          ! calculate height of radar beam above MSL:
          h1 = SQRT((re+altsta)**2 + ra(ira)**2 + 2*(re+altsta)*ra(ira)*SIN(el(iaz,iel)*degrad)) - re
          h(irp) =  h1
          ! calculate length of circle at height of MSL:
          s1 = re*ASIN(ra(ira)*COS(el(iaz,iel)*degrad)/(re+h1))

          CALL polar2geo_old (lonsta, latsta, r_earth_dp, s1, az(iaz), lon(irp), lat(irp))

        END DO
      END DO
    END DO
!$omp end parallel do

    ! clean up:
    DEALLOCATE (az, ra, el)
    IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE rad2geo_const_vec


  SUBROUTINE rad2geo_const_smth(rsm,lat,lon,h)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation of radar coordinates
    !              (range, azimuth and elevation) to
    !              geographical coordinates (lat, lon and height asl)
    !
    ! Method:      Assumption: 4/3 earth model
    !              Formulas with reference to  Appendices B and D of Dissertation
    !              of Uli Blahak
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    TYPE(radar_meta_type), INTENT(in) :: rsm
    REAL (KIND=dp), ALLOCATABLE, INTENT (inout)     ::        &
         lat(:),& ! geographical latitude vector (nra*naz*nel)
         lon(:),& ! geographical longitude vector (nra*naz*nel)
         h  (:)   ! height above mean sea level vector (nra*naz*nel)

    !
    ! Local scalars:

    REAL (KIND=dp)      ::        &
         latsta,  & ! geographical latitude of radar station
         lonsta,  & ! geographical longitude of radar station
         altsta     ! altitude of radar station asl (m)

    REAL (KIND=dp), ALLOCATABLE      ::        &
         ra(:),   & ! array of ranges
         az(:),   & ! array of azimuths
         el(:,:), & ! array of elevations
         el_precip(:)  ! vector of elevations for precip scan

    CHARACTER (LEN=32) yzroutine

    INTEGER :: ira,iaz,iel,irp,iv,ih,nobsmax
    REAL    (KIND=dp)    :: h1,s1,d1,d2,d3,re
    REAL    (KIND=dp)    ::           &
         el_tmp , & ! elevation of smoothing point
         az_tmp     ! azimuth of smoothing point

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE rad2geo_const_smth
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine = 'rad2geo_const_smth'

    ! allocate and fill up range and azimuth arrays
    CALL get_rangevec (rsm, ra)
    CALL get_azvec    (rsm, az)

    ALLOCATE(el(rsm%naz,rsm%nel))
    IF (TRIM(rsm%scanname) == 'PRECIP') THEN
      CALL get_elarr_precipscan (rsm, el_precip)
      el(:,1) = el_precip
    ELSE
      DO iel=1, rsm%nel
        el(:,iel) = rsm%el_arr(iel)
      END DO
    END IF

    nobsmax = rsm%naz*rsm%nra*rsm%nel*rsm%ngpsm_h*rsm%ngpsm_v
    IF (.NOT.ALLOCATED(lat)) ALLOCATE( lat(nobsmax))
    IF (.NOT.ALLOCATED(lon)) ALLOCATE( lon(nobsmax))
    IF (.NOT.ALLOCATED(h  )) ALLOCATE( h  (nobsmax))

    latsta = rsm%lat
    lonsta = rsm%lon
    altsta = rsm%alt_msl

    re = 4.0_dp/3.0_dp * r_earth_dp

    DO ih = 1,rsm%ngpsm_h
      DO iv = 1,rsm%ngpsm_v
        DO iel = 1,rsm%nel
          DO ira = 1,rsm%nra
            DO iaz = 1,rsm%naz

              ! Determine the index of the radar point, based on az,ra,el
              CALL sub2ind5D(iaz,ira,iel,iv,ih, rsm%naz,rsm%nra,rsm%nel,rsm%ngpsm_v, irp)

              el_tmp = el(iaz,iel) + smth_el_horzscan(rsm%Theta3,rsm%smth_interv_fact,rsm%xabscsm_v(iv))

              ! calculate height of radar beam above MSL:
              h1 = SQRT((re+altsta)**2 + ra(ira)**2 + 2*(re+altsta)*ra(ira)*SIN(el_tmp*degrad)) - re
              h(irp) = h1
              ! calculate length of circle at height of MSL:
              s1 = re*ASIN(ra(ira)*COS(el_tmp*degrad)/(re+h1))

              az_tmp = az(iaz) + smth_az_horzscan( &
                                                   rsm%alpha3_eff_0, &
                                                   rsm%dalpha, &
                                                   rsm%Phi3, &
                                                   rsm%smth_interv_fact, &
                                                   rsm%xabscsm_h(ih), &
                                                   rsm%el_arr(iel) &
                                                   )
              ! Note: for DWD precip scan, rsm%el_arr(1) is a suitable approximation, which is used throughout the code for the eff. beam function

              CALL polar2geo_old (lonsta, latsta, r_earth_dp, s1, az_tmp, lon(irp), lat(irp))

            END DO
          END DO
        END DO
      END DO
    END DO

    ! clean up:
    DEALLOCATE (az, ra, el)
    IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE rad2geo_const_smth


  SUBROUTINE grd2geo(latsta,lonsta,naz,nal,az,al,lat,lon)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation of auxiliary grids coordinates
    !              (azimuth, arc length) to geographical coordinates (lat, lon)
    !
    ! Method:      Compute geographic coordinates from given radar site,
    !              azimut and arc distance at mean sea level
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    REAL (KIND=dp), INTENT (IN)      ::        &
         latsta,  & ! geographical latitude of radar station
         lonsta     ! geographical longitude of radar station

    INTEGER, INTENT (IN) :: &
         naz,nal! number of azimuths, arc distances (array dimensions)

    REAL (KIND=dp),     INTENT (IN) :: &
         az(naz), & ! array of azimuths
         al(nal)    ! array of distance

    REAL (KIND=dp), INTENT (OUT)     ::        &

         lat(naz,nal),& ! geographical latitude
         lon(naz,nal)   ! geographical longitude

    ! Local scalars:
    CHARACTER (LEN=32) yzroutine

    INTEGER           :: iaz,ial
    REAL (KIND=dp)    :: d1,d2,d3


    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE  grd2geo
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine = 'grd2geo'


    DO ial = 1,nal
      DO iaz = 1,naz

        CALL polar2geo_old (lonsta, latsta, r_earth_dp, al(ial), az(iaz), lon(iaz,ial), lat(iaz,ial))

      END DO
    END DO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE grd2geo

  SUBROUTINE online2geo(rsm,al,hl,el_loc,lat,lon)

    !------------------------------------------------------------------------------
    !
    ! Description: Coordinate transformation from azimuth and arc length to
    !              geographical coordinates (lat, lon) from the informations given by
    !              the online beam propagation routines. This routine is written
    !              for input field al (arc distance) as a 3D polar data set.
    !              Blocked ray parts are extrapolated by using the 4/3 earth radius
    !              model, and missing values in the arc distance and height fields are
    !              correspondingly filled.
    !
    ! Method:      Compute geographic coordinates from given radar site,
    !              azimut and arc distance at mean sea level. For blocked parts of
    !              rays, where the arc length is miss_value, the 4/3-earth model
    !              is used to fill in values here, extrapolating from the
    !              first blocked bin minus 1 outwards along the range.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    ! Parameter list:
    TYPE(radar_meta_type), INTENT(in) :: rsm

    REAL (KIND=dp),     INTENT (INOUT) :: &
         hl(:,:,:),   & ! height of ray from online beam propagation (naz,nra,nel))
         al(:,:,:)      ! array of arc distance at sea level         (naz,nra,nel)

    REAL (KIND=dp),     INTENT (IN) :: &
         el_loc(:,:,:)  ! elevation of ray from online beam propagation (naz,nra,nel)

    REAL (KIND=dp), INTENT (OUT)     ::        &
         lat(:,:,:),& ! geographical latitude  (naz,nra,nel)
         lon(:,:,:)   ! geographical longitude (naz,nra,nel)

    ! Local scalars and arrays:
    CHARACTER (LEN=32) :: yzroutine

    INTEGER              :: iaz,ira,iel,ira43(rsm%naz,rsm%nel)
    REAL    (KIND=dp)    :: d1,d2,d3,re,re0,el0,al0,dra

    REAL (KIND=dp), ALLOCATABLE :: &
         az(:),         & ! vector of azimuths
         ra(:),         & ! vector of ranges
         el(:,:),       & ! array of elevations
         el_precip(:)     ! vector of elevations for precip scan

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE  online2geo
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'online2geo'

    re = 4.0_dp/3.0_dp * r_earth_dp

    CALL get_rangevec ( rsm, ra )
    CALL get_azvec ( rsm, az )

    ALLOCATE(el(rsm%naz,rsm%nel))
    IF (TRIM(rsm%scanname) == 'PRECIP') THEN
      CALL get_elarr_precipscan (rsm, el_precip)
      el(:,1) = el_precip
    ELSE
      DO iel=1, rsm%nel
        el(:,iel) = rsm%el_arr(iel)
      END DO
    END IF

    ! to store the first blocked range bin minus 1 for each ray (azi, ele):
    ira43(:,:) = -1

    DO iel = 1,rsm%nel
      DO ira = 1,rsm%nra
        DO iaz = 1,rsm%naz

          IF (al(iaz,ira,iel) < miss_threshold) THEN

            ! Bin is blocked. If no blocking was detected until now, the last range is the valid
            !  range for extrapolation outwards:
            IF (ira43(iaz,iel) == -1) ira43(iaz,iel) = ira - 1

            ! Extrapolate from the first blocked range bin minus 1 using the 4/3 earth radius model:
            IF (ira43(iaz,iel) == 0) THEN
              ! The ray is totally blocked, so extrapolate from the radar station:
              re0 = re + rsm%alt_msl
              el0 = el(iaz,iel)
              al0 = 0.0_dp
              dra = ra(ira)
            ELSE
              ! Extrapolate from the first blocked bin minus 1:
              re0 = re + hl(iaz,ira43(iaz,iel),iel)
              el0 = el_loc (iaz,ira43(iaz,iel),iel)
              al0 = al     (iaz,ira43(iaz,iel),iel)
              dra = ra(ira) - ra(ira43(iaz,iel))
            END IF

            ! Extrapolated height over MSL and arc distance (update INOUT hl and al fields):
            hl(iaz,ira,iel) = SQRT(re0**2 + dra**2 + 2.0*re0*dra*SIN(el0*degrad)) - re
            al(iaz,ira,iel) = al0 + re*ASIN(dra*COS(el0*degrad)/(re+hl(iaz,ira,iel)))

          END IF

          CALL polar2geo_old (rsm%lon, rsm%lat, r_earth_dp, al(iaz,ira,iel), az(iaz), lon(iaz,ira,iel), lat(iaz,ira,iel))

        END DO
      END DO
    END DO

    ! clean up:
    DEALLOCATE (az, ra, el)
    IF (ALLOCATED(el_precip)) DEALLOCATE(el_precip)

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE online2geo

  !=========================================================================================
  !=========================================================================================


  SUBROUTINE calc_vert_weight_online(np,alt,hfl,azarr_idx,naz_loc,istart_az,iend_az,nal,nhl,iaz,ial,weight,kout,wk)

    !------------------------------------------------------------------------------
    !
    ! Description: Calculation of vertical interpolation weight and level index
    !              above the radar point for the interpolation of the auxiliary grids
    !              to the radar points.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    ! Parameter list:

    CHARACTER (LEN=32) yzroutine

    INTEGER, INTENT (IN)     ::        &
         np,istart_az,iend_az,naz_loc,nal,nhl  , &  ! np = no. of radar points, naz_loc = size of azarr
         azarr_idx(naz_loc)                ,&           ! azimut index range of this processor in the aux. grid
         iaz(np)                       ,&
         ial(np)

    INTEGER :: i, kbot, kbotp1, kbotp2, ktop, kincr

    REAL (KIND=dp), INTENT (IN)           ::        &
         alt(np)               ,& !
         weight(np)            ,& ! interpolation weight in i direction
         hfl(istart_az:iend_az,nal,nhl)         ! height a.s.l of radar point

    INTEGER, INTENT (OUT)     ::       &
         kout(np)             ! lowest level index above radar point
                              ! -1 if radar point below topography

    REAL (KIND=dp), INTENT (OUT)     ::             &
         wk(np)               ! vertical interpolation weight
                              ! -1.0 if radar point below topography

    INTEGER :: m,mo,n,no,k,irp,flag(np) ! use flag to check if the radar point has already been done

    REAL    (KIND=dp)    :: hfl_low

    ! Local arrays:
    REAL    (KIND=dp)    :: hfl_up(np) ! full level heights of the 4 grid points
    ! surrounding the radar point
    !- End of header
    !==============================================================================

    yzroutine(:) = ' '
    yzroutine = 'calc_vert_weight_online'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! initialize vertical grid indices:
    kbot   = bottomlevel ()
    kbotp1 = one_level_up(kbot)
    kbotp2 = one_level_up(kbotp1)
    ktop   = toplevel()
    kincr  = levelincr()

    ! initialize vertical interpolation weight and level index
    wk = -1.0_dp
    kout = -1
    flag = 0

   ! Check if the radar point is below topography.
   ! If yes, then flag = 1 (wk,kout = default value); Elseif below the first main level, then flag = 1, wk = 0, kout = bottomlevel().

    DO irp = 1, np

      m = azarr_idx(iaz(irp))
      mo= MIN(m+1,iend_az)
      n = ial(irp)
      no= MIN(n+1,nal)

      IF ( hfl(m,n,kbot)     > miss_threshold .AND. &
           hfl(m,no,kbot)   > miss_threshold .AND. &
           hfl(m,n,kbotp1)   > miss_threshold .AND. &
           hfl(m,no,kbotp1) > miss_threshold       &
           ) THEN

        hfl_low = hfl(m,n,kbot) * (1.0_dp-weight(irp)) + hfl(m,no,kbot) * weight(irp)

        hfl_up(irp) = hfl(m,n,kbotp1) * (1.0_dp-weight(irp)) + hfl(m,no,kbotp1) * weight(irp)


!!$ UB>> ???      IF (alt(irp) < hfl_low .OR. alt(irp) < 0.0 .OR. hfl_low < 0.0 .OR. hfl_up(irp,kbot) < 0.0) THEN
        IF (alt(irp) < hfl_low .OR. alt(irp) < miss_threshold ) THEN

          flag(irp) = 1

        ELSEIF (hfl_up(irp) >= alt(irp)) THEN

          wk(irp) = (alt(irp)-hfl_up(irp))/(hfl_low-hfl_up(irp))

          kout(irp) = kbotp1

          flag(irp) = 1

        ENDIF

      ELSE

        hfl_up(irp) = miss_value

      ENDIF

    ENDDO

    ! Loop over k from second lowest level to the top
    ! Determine the k indice above the radar point
    DO k = kbotp2, ktop, kincr

      DO irp =1, np

        IF (flag(irp) == 0) THEN

          hfl_low = hfl_up(irp)

          m = azarr_idx(iaz(irp))
          mo= MIN(m+1,iend_az)
          n = ial(irp)
          no= MIN(n+1,nal)

          IF ( hfl(m,n,k)     > miss_threshold .AND. &
               hfl(m,no,k)   > miss_threshold       &
               ) THEN

            hfl_up(irp) = hfl(m,n,k) * (1.0_dp-weight(irp)) +  hfl(m,no,k) * weight(irp)

            IF (alt(irp) < miss_threshold) THEN

              flag(irp) = 1

            ELSEIF ( hfl_low >= miss_threshold .AND. &
                     hfl_low <= alt(irp) .AND. alt(irp) < hfl_up(irp) ) THEN

              wk(irp) = (alt(irp)-hfl_up(irp))/(hfl_low-hfl_up(irp))

              kout(irp) = k

              flag(irp) = 1

            END IF

          ELSE

            hfl_up(irp) = miss_value

          ENDIF

        ENDIF ! flag(irp) == 0

      ENDDO

    ENDDO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_vert_weight_online


  SUBROUTINE calc_vert_weight_onsmth(np,alt,hfl,istart_az,iend_az,nal,nhl,iaz,ial,wa,wl,kout,wk)

    !-----------------------------------------------------------------------------------
    !
    ! Description: Calculation of vertical interpolation weight and level index
    !              above the smoothing point for the interpolation of the auxiliary grids
    !              to the smoothing points. All input fields and computation are on the
    !              auxiliary azimutal grid.
    !
    !-----------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    ! Parameter list:

    CHARACTER (LEN=32) yzroutine

    INTEGER, INTENT (IN)     ::        &
         np,istart_az,iend_az,nal,nhl, &
         iaz(np),                      &
         ial(np)

    INTEGER   :: i, kbot, kbotp1, kbotp2, ktop, kincr

    REAL (KIND=dp), INTENT (IN)           ::        &
         alt(np)               ,& !
         wa(np)                ,& ! interpolation weight in i direction
         wl(np)                ,&
         hfl(istart_az:iend_az,nal,nhl)         ! height a.s.l of radar point

    INTEGER, INTENT (OUT)     ::       &
         kout(np)             ! lowest level index above radar point
    ! -1 if radar point below topography

    REAL (KIND=dp), INTENT (OUT)     ::             &
         wk(np)               ! vertical interpolation weight
    ! -1.0 if radar point below topography

    INTEGER :: k,m,mo,n,no,irp,flag(np) ! use flag to check if the radar point has already been done

    REAL    (KIND=dp)    :: w_a,w_l

    REAL    (KIND=dp)    :: hfl_low,hfl_low11,hfl_low21,hfl_low12,hfl_low22

    ! Local arrays:
    REAL    (KIND=dp)    :: hfl_up(np),hfl_up11,hfl_up21,hfl_up12,hfl_up22 ! full level heights of the 4 grid points
    ! surrounding the radar point
    !- End of header
    !==============================================================================

    yzroutine(:) = ' '
    yzroutine = 'calc_vert_weight_onsmth'

    IF (ldebug_radsim) WRITE (*,*) TRIM(yzroutine), ' on proc ', my_radar_id

    ! initialize vertical grid indices:
    kbot   = bottomlevel ()
    kbotp1 = one_level_up(kbot)
    kbotp2 = one_level_up(kbotp1)
    ktop   = toplevel()
    kincr  = levelincr()

    ! initialize vertical interpolation weight and level index
    wk = -1.0_dp
    kout = -1
    flag = 0

    ! Check if the radar point is below topography.
    ! If yes, then flag = 1 (wk,kout = default value); Elseif below the first main level, then flag = 1, wk = 0, kout = kbot.

    DO irp = 1, np

      m = iaz(irp)
      mo= MIN(m+1,iend_az)
      n = ial(irp)
      no= MIN(n+1,nal)

      hfl_low11 = hfl(m,n,kbot)
      hfl_low21 = hfl(mo,n,kbot)
      hfl_low12 = hfl(m,no,kbot)
      hfl_low22 = hfl(mo,no,kbot)

      hfl_up11 = hfl(m,n,kbotp1)
      hfl_up21 = hfl(mo,n,kbotp1)
      hfl_up12 = hfl(m,no,kbotp1)
      hfl_up22 = hfl(mo,no,kbotp1)

      ! Select only azimut-slice-grid intervals which are fully in the model domain:
      IF ( hfl_low11 >= miss_threshold  .AND. &
           hfl_low12 >= miss_threshold  .AND. &
           hfl_low21 >= miss_threshold  .AND. &
           hfl_low22 >= miss_threshold  .AND. &
           hfl_up11  >= miss_threshold  .AND. &
           hfl_up12  >= miss_threshold  .AND. &
           hfl_up21  >= miss_threshold  .AND. &
           hfl_up22  >= miss_threshold        &
           ) THEN

        w_a = wa(irp)
        w_l = wl(irp)

        hfl_low = (hfl_low11*(1.0_dp-w_a)+hfl_low21*w_a)*(1.0_dp-w_l) + &
             (hfl_low12*(1.0_dp-w_a)+hfl_low22*w_a)*w_l

        hfl_up(irp) = (hfl_up11*(1.0_dp-w_a)+hfl_up21*w_a)*(1.0_dp-w_l) + &
             (hfl_up21*(1.0_dp-w_a)+hfl_up22*w_a)*w_l


!!$ UB>>      IF (alt(irp) < hfl_low .or. alt(irp) < 0.0 .or. hfl_low < 0.0 .OR. hfl_up(irp,kbot) < 0.0) THEN
        IF (alt(irp) < hfl_low .OR. alt(irp) < miss_threshold) THEN

          flag(irp) = 1

        ELSEIF(hfl_up(irp) >= alt(irp)) THEN

          wk(irp) = (alt(irp)-hfl_up(irp))/(hfl_low-hfl_up(irp))

          kout(irp) = kbotp1

          flag(irp) = 1

        ENDIF

      ELSE

        ! Do not set flag = 1, because the search should not end here!
        hfl_up(irp) = miss_value

      END IF

    ENDDO

    ! Loop over k from second lowest level to the top
    ! Determine the k indice above the radar point
    DO k = kbotp2, ktop, kincr

      DO irp =1, np

        IF (flag(irp) == 0) THEN

          hfl_low = hfl_up(irp)

          m = iaz(irp)
          mo= MIN(m+1,iend_az)
          n = ial(irp)
          no= MIN(n+1,nal)

          hfl_up11 = hfl(m,n,k)
          hfl_up21 = hfl(mo,n,k)
          hfl_up12 = hfl(m,no,k)
          hfl_up22 = hfl(mo,no,k)

          IF ( hfl_up11 >= miss_threshold .AND. &
               hfl_up21 >= miss_threshold .AND. &
               hfl_up12 >= miss_threshold .AND. &
               hfl_up22 >= miss_threshold &
               ) THEN

            w_a = wa(irp)
            w_l = wl(irp)

            hfl_up(irp) = (hfl_up11*(1.0_dp-w_a)+hfl_up21*w_a)*(1.0_dp-w_l) + &
                 (hfl_up21*(1.0_dp-w_a)+hfl_up22*w_a)*w_l

            IF (alt(irp) < miss_threshold) THEN
              ! This branch should not be necessary (has been
              ! treated alread in the first loop), but for safety ...

              flag(irp) = 1

            ELSEIF ( hfl_low >= miss_threshold .AND. &
                     hfl_low <= alt(irp) .AND. alt(irp) < hfl_up(irp) ) THEN

              wk(irp) = (alt(irp)-hfl_up(irp))/(hfl_low-hfl_up(irp))

              kout(irp) = k

              flag(irp) = 1

            ENDIF ! hfl_up(irp,k) > alt(irp)

          ELSE

            hfl_up(irp) = miss_value

          ENDIF

        ENDIF ! flag(irp) == 0

      ENDDO

    ENDDO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE calc_vert_weight_onsmth



  SUBROUTINE calc_height_rk4(nst,nae,nra,hl1,dhdr1,rfridx,ra_inc,flag,hl2,dhdr2)

    !-----------------------------------------------------------------------------------
    !
    ! Description: Calculation of vertical interpolation weight and level index
    !              above the smoothing point for the interpolation of the auxiliary grids
    !              to the smoothing points
    !-----------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    !------------------------------------------------------------------------------
    ! Parameter list:

    CHARACTER (LEN=32) yzroutine

    INTEGER, INTENT (IN)     :: nst,nae,nra

    INTEGER, INTENT (IN)     :: flag(nae)

    INTEGER   :: ist,iae

    REAL (KIND=dp), INTENT (IN)           :: ra_inc

    REAL (KIND=dp), INTENT (IN)           ::             &
         hl1(nae)               ,& !
         dhdr1(nae)             ,& !
         rfridx(nae)

    REAL (KIND=dp)           :: k1(nae,2),k2(nae,2),k3(nae,2),k4(nae,2),f(nae,2),hltmp(nae),dhdrtmp(nae)

    REAL (KIND=dp)           :: ra_st

    REAL (KIND=dp), INTENT (OUT)          ::             &
         hl2(nae)               ,& !
         dhdr2(nae)

    REAL (kind=dp), PARAMETER ::  inv_3 = 1.0_dp/3.0_dp
    REAL (kind=dp), PARAMETER ::  inv_6 = 1.0_dp/6.0_dp

    !- End of header
    !==============================================================================

    yzroutine(:) = ' '
    yzroutine = 'calc_height_rk4'

    k1 = miss_value
    k2 = miss_value
    k3 = miss_value
    k4 = miss_value

    f  = miss_value

    hltmp   = hl1
    dhdrtmp = dhdr1

    hl2     = miss_value
    dhdr2   = miss_value

    ra_st   = ra_inc/nst

    DO ist = 1, nst

      DO iae = 1, nae

        IF (flag(iae) > 0) THEN

          k1(iae,1) = dhdrtmp(iae)
          k1(iae,2) = (1.0_dp - dhdrtmp(iae)**2)* &
               (rfridx(iae) + 1.0_dp/(r_earth_dp + hltmp(iae)))

          k2(iae,1) = dhdrtmp(iae) + 0.5_dp*k1(iae,2)*ra_st
          k2(iae,2) = (1.0_dp - (dhdrtmp(iae)+0.5_dp*k1(iae,2)*ra_st)**2)* &
               (rfridx(iae) + 1.0_dp/(r_earth_dp + hltmp(iae) + 0.5_dp*k1(iae,1)*ra_st))

          k3(iae,1) = dhdrtmp(iae) + 0.5_dp*k2(iae,2)*ra_st
          k3(iae,2) = (1.0_dp - (dhdrtmp(iae)+0.5_dp*k2(iae,2)*ra_st)**2)* &
               (rfridx(iae) + 1.0_dp/(r_earth_dp + hltmp(iae) + 0.5_dp*k2(iae,1)*ra_st))

          k4(iae,1) = dhdrtmp(iae) + k3(iae,2)*ra_st
          k4(iae,2) = (1.0_dp - (dhdrtmp(iae)+k3(iae,2)*ra_st)**2)* &
               (rfridx(iae) + 1.0_dp/(r_earth_dp + hltmp(iae) + k3(iae,1)*ra_st))

          f(iae,1) = inv_6*k1(iae,1) + inv_3*k2(iae,1) + inv_3*k3(iae,1) + inv_6*k4(iae,1)
          f(iae,2) = inv_6*k1(iae,2) + inv_3*k2(iae,2) + inv_3*k3(iae,2) + inv_6*k4(iae,2)

          hltmp(iae)   = hltmp(iae) + ra_st*f(iae,1)
          dhdrtmp(iae) = dhdrtmp(iae) + ra_st*f(iae,2)

        ENDIF

      ENDDO

    ENDDO

    hl2     = hltmp
    ! dhdr2 = dhdrtmp
    DO iae = 1, nae

      IF (flag(iae) > 0) THEN
        dhdr2(iae) = MIN(MAX(dhdrtmp(iae),-1.0_dp),1.0_dp)
      ENDIF

    ENDDO

  END SUBROUTINE calc_height_rk4


  SUBROUTINE para_range(n1,n2,nprocs,nbl_az,istart,iend)

    !------------------------------------------------------------------------------
    ! Description: This subroutine calculates the range of azimuts of a single station
    !              for each processor
    !              assuming a parallelization model where the azimuts of the
    !              station are evenly divided among the PEs, so that each PE
    !              gets an equal share of azimuts of the station.
    !
    !              If there are no azimuts on a certain PE,
    !              istart and iend will be returned as -1. This happens when
    !              the number of PEs is larger than the number of azimuts.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER, INTENT(IN)       :: n1,n2,nprocs,nbl_az
    INTEGER, INTENT(OUT)      :: istart(0:nprocs-1),iend(0:nprocs-1)
    INTEGER                   :: iwork1,iwork2,irank,ntot
    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE para_range
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'para_range'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! Initialization with missing values:
    istart = -1
    iend   = -1

    ! Computation:
    ntot   = n2 - n1 + 1
    iwork1 = ntot/nprocs
    iwork2 = MOD(ntot,nprocs)

    DO irank = 0, MIN(nprocs,ntot)-1

      istart(irank) = irank*iwork1 + n1 + MIN(irank,iwork2)
      iend(irank)   = istart(irank) + iwork1 - 1
      IF (iwork2 > irank) THEN
        iend(irank) = iend(irank) + 1
      ENDIF

      ! UB>> Append 2*nbl_az overlap at the upper bound of the azimut index range:
      !      (however, physically nbl_az azimut intervals are append both the lower and upper bound)
      iend(irank) = iend(irank) + 2*nbl_az

    ENDDO

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE para_range

  SUBROUTINE para_range_all(n1,n2,nsta,nprocs,nbl_az,istart,iend)

    !------------------------------------------------------------------------------
    ! Description: This subroutine calculates the range of azimuts for each processor
    !              assuming a parallelization strategy where the azimuts of all
    !              stations are combined into one long azimut vector, which is then divided
    !              evenly among the PEs.
    !
    !              If there are no azimuts of a certain station on a certain PE,
    !              istart and iend will be returned as -1. This happens when
    !              the number of PEs is larger than the total number of azimuts
    !              of all stations.
    !
    !------------------------------------------------------------------------------
    !
    ! Subroutine / Function arguments
    ! Scalar arguments with intent(in):

    !------------------------------------------------------------------------------
    !
    ! Local scalars:
    CHARACTER (LEN=32) yzroutine
    CHARACTER (LEN=80) yzerrmsg

    INTEGER, INTENT(IN)       :: nsta,n1(nsta),n2(nsta),nprocs,nbl_az
    INTEGER, INTENT(OUT)      :: istart(0:nprocs-1,nsta),iend(0:nprocs-1,nsta)
    INTEGER                   :: iwork1,iwork2,irank,ista,iop,ntot,irank_start,irank_end,&
                                 nranks_para,istart_tot(0:nprocs-1),iend_tot(0:nprocs-1),&
                                 offs(0:nsta),i,j
    LOGICAL, SAVE :: firstcall = .true.

    !- End of header
    !==============================================================================

    !------------------------------------------------------------------------------
    !- Begin SUBROUTINE para_range
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Section 1: Initializations
    !------------------------------------------------------------------------------

    yzroutine(:) = ' '
    yzroutine = 'para_range_all'

    IF (ldebug_radsim) WRITE (*,*)  TRIM(yzroutine), ' on proc ', my_radar_id

    ! Initialization with missing values:
    ! -----------------------------------

    istart = -1
    iend   = -1

    ! Computation:
    ! ------------

    irank_start = 1
    irank_end = nprocs-1
    nranks_para = irank_end - irank_start + 1

    ! 1) Divide the combined azimut vector of all stations evenly among the PEs:
    !    (The combined vector starts at 1, not n1(1)!)
    ntot   = SUM(n2 - n1 + 1)
    iwork1 = (ntot)/nranks_para
    iwork2 = MOD(ntot,nranks_para)

    IF (irank_start > 0) THEN
      istart_tot(0:irank_start-1) = -1
      iend_tot(0:irank_start-1) = -1
    END IF

    DO irank = 0, MIN(nranks_para,ntot)-1

      istart_tot(irank+irank_start) = irank*iwork1 + 1 + MIN(irank,iwork2)
      iend_tot(irank+irank_start)   = istart_tot(irank+irank_start) + iwork1 - 1
      IF (iwork2 > irank) THEN
        iend_tot(irank+irank_start) = iend_tot(irank+irank_start) + 1
      ENDIF

    ENDDO

    ! 2) Calculate the index range for each station on each PE:

    !   Offset of the index range for each station in the combined azimut vector:
    offs(0) = 0
    DO ista = 1, nsta
      offs(ista) = offs(ista-1) + n2(ista) - n1(ista) + 1
    END DO

    !   now hop through the combined azimut vector from one station boundary (offs)
    !    or PE boundary (iend_tot) to the next station boundary or PE boundary (whichever
    !    comes first), to determine the index range of each station on each PE:

    ista = 1   ! Hop starting value for station (will be increased until nsta)
    iop  = irank_start   ! Hop starting value for PE (will be increased until irank_end)
    DO
      ! UB>> Append 2*nbl_az overlap at the upper bound of the azimut index range:
      !      (however, physically nbl_az azimut intervals are append both the lower and upper bound)
      istart(iop,ista) = MAX( 1, iend(MAX(iop-1,0),ista)-n1(ista)+1-2*nbl_az + 1 ) + n1(ista)-1
      iend(iop,ista)   = MIN( iend_tot(iop)-offs(ista-1), n2(ista)-n1(ista)+1    ) + n1(ista)-1 + 2*nbl_az
      IF (iend_tot(iop) <= offs(ista)) THEN
        iop = iop + 1
        IF (iop > MIN(irank_end+1,ntot)-1) EXIT
        IF (iend(iop-1,ista) == n2(ista)+2*nbl_az) THEN
          ista = ista + 1
          IF (ista > nsta) EXIT
        END IF
      ELSE
        ista = ista + 1
        IF (ista > nsta) EXIT
      END IF
    END DO

!!$IF (my_radar_id == 1 .and. firstcall) THEN
!!$  WRITE (*,'(/,a,i4.4,a,/)') 'ULIproc',my_radar_id,' para_range_all:'
!!$  DO i=0, nprocs-1
!!$    WRITE (*,'(a,i4.4,a,i4,a)') 'ULIproc',my_radar_id,' Proc ',i,'  : '
!!$    DO j=1, nsta
!!$      WRITE (*,'(a,i4.4,a,i3,a,2(1x,i5))') 'ULIproc',my_radar_id,'  ista ',j,'  : ', istart(i,j), iend(i,j)
!!$    END DO
!!$    WRITE (*,*)
!!$  END DO
!!$  firstcall = .FALSE.
!!$END IF

    IF (ldebug_radsim) WRITE (*,*)  'done with '//TRIM(yzroutine)//' on proc ', my_radar_id

  END SUBROUTINE para_range_all

  !=========================================================================
  !
  ! Subroutine for allocating and pre-setting the azimuth vector of a
  !  radar station ista based on the information found in the meta data
  !  of this radar station.
  !
  !=========================================================================

  SUBROUTINE get_azvec ( rsm, azarr )

    TYPE(radar_meta_type), INTENT(in)              :: rsm
    REAL (KIND=dp), INTENT(inout), ALLOCATABLE :: azarr(:)

    IF (.NOT.ALLOCATED(azarr)) ALLOCATE (azarr(rsm%naz))
    CALL get_azvec_values (rsm%naz, rsm%az_inc, rsm%az_start, 0, azarr)

  END SUBROUTINE get_azvec

  SUBROUTINE get_azvec_values ( naz, az_inc, az_start, iaz_offset, azarr )

    REAL (KIND=dp),       INTENT(in)   :: az_inc, az_start
    INTEGER,              INTENT(in)   :: naz, iaz_offset
    REAL (KIND=dp),       INTENT(out)  :: azarr(naz)

    INTEGER :: iaz

    DO iaz = 1, naz
      azarr(iaz) = (iaz-iaz_offset-1) * az_inc + az_start
      ! Folding into the range [0, 360] degrees
      azarr(iaz) = MODULO(azarr(iaz),360.0_dp)
    END DO

  END SUBROUTINE get_azvec_values

  SUBROUTINE get_rangevec ( rsm, raarr )

    TYPE(radar_meta_type), INTENT(in)              :: rsm
    REAL (KIND=dp), INTENT(inout), ALLOCATABLE :: raarr(:)

    IF (.NOT.ALLOCATED(raarr)) ALLOCATE (raarr(rsm%nra))
    CALL get_rangevec_values (rsm%nra, rsm%ra_inc, rsm%ra_inc, raarr)

  END SUBROUTINE get_rangevec

  SUBROUTINE get_rangevec_values ( nra, ra_inc, ra_start, raarr )

    REAL (KIND=dp),       INTENT(in)   :: ra_inc, ra_start
    INTEGER,              INTENT(in)   :: nra
    REAL (KIND=dp),       INTENT(out)  :: raarr(nra)

    INTEGER :: ira

    DO ira = 1, nra
      raarr(ira) = (ira-1) * ra_inc + ra_start
    END DO

  END SUBROUTINE get_rangevec_values




END MODULE radar_model2rays
