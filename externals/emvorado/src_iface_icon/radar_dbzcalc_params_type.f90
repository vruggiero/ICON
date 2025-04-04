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


MODULE radar_dbzcalc_params_type

!-------------------------------------------------------------------------------
!
! Description:
!  This module declares a derived type "dbzcalc_params" which holds the meta data for the
!  reflectivity compuations in calc_dzb_vec()
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

#ifdef __ICON__
#ifdef HAVE_RADARFWO
  USE radar_kind, ONLY : wpfwo  
#else
  USE mo_kind, ONLY : wpfwo => wp  ! only to enable compilation!
#endif
#else
  USE radar_kind, ONLY : wpfwo
#endif
  USE radar_data_mie, ONLY : T0C_fwo ! T[K] at 0 deg C
  
  IMPLICIT NONE

  !==============================================================================

  ! INCLUDE statements

  !==============================================================================

  PUBLIC

  !------------------------------------------------------------------------------

  TYPE t_polMP
    SEQUENCE   
    CHARACTER (LEN=8)      :: ARmodel      ! AR scheme to apply (Reference abbrev., eg 'R11', or 'Poly' for polynomial)
    REAL      (KIND=wpfwo) :: c0           ! AR 0th order polynomial coeff
    REAL      (KIND=wpfwo) :: c1           ! AR 1th order polynomial coeff
    REAL      (KIND=wpfwo) :: ARmin        ! low bound of AR
    REAL      (KIND=wpfwo) :: sig          ! canting angle distribution width
  END TYPE t_polMP

  ! Type to hold the namelist parameters configuring radar reflectivity calculation for
  ! routine calc_dbz_vec() from radar_mie_iface_cosmo_driver.f90:
  
  TYPE t_dbzcalc_params

    SEQUENCE   ! Important: ensures that the following parameters are one block in memory
               !       and that an MPI-distribution by a derived MPI-datatype is possible:
               !       - radar_parallel_utilities.f90, def_mpi_dbzcalc_params_type()
               !       - CALL def_mpi_dbzcalc_params_type (mpi_dbzcalc_params_typ)
               !         MPI_BCAST( <dbzcalc_params>, 1, ..., mpi_dbzcalc_params_typ, ...)
               !       For this, the following data have to be naturally aligned, i.e.,
               !       first the REAL, then the INTEGER, LOGICAL and CHARACTER at the end.

    REAL    (KIND=wpfwo)       :: lambda_radar       ! radar wavelength
    ! melting scheme parameters
    REAL    (KIND=wpfwo)       :: Tmeltbegin_i       ! temperature[K], above which ice is assumed wet
    !REAL    (KIND=wpfwo)       :: Tmin_i             ! Temperature[K] at which smallest ice particles melt
    REAL    (KIND=wpfwo)       :: meltdegTmin_i      ! degree of ice melting at T=Tmin_i
    REAL    (KIND=wpfwo)       :: Tmax_min_i         ! minimum Tmax[K] of ice in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: Tmax_max_i         ! maximum Tmax[K] of ice in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: qthresh_i          ! q threshold to apply dynamic melting scheme for ice
    REAL    (KIND=wpfwo)       :: qnthresh_i         ! qn threshold to apply dynamic melting scheme for ice
    REAL    (KIND=wpfwo)       :: Tmeltbegin_s       ! temperature[K], above which snow is assumed wet
    !REAL    (KIND=wpfwo)       :: Tmin_s             ! Temperature[K] at which smallest snow particles melt
    REAL    (KIND=wpfwo)       :: meltdegTmin_s      ! degree of snow melting at T=Tmin_s
    REAL    (KIND=wpfwo)       :: Tmax_min_s         ! minimum Tmax[K] of snow in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: Tmax_max_s         ! maximum Tmax[K] of snow in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: qthresh_s          ! q threshold to apply dynamic melting scheme for snow
    REAL    (KIND=wpfwo)       :: qnthresh_s         ! qn threshold to apply dynamic melting scheme for snow
    REAL    (KIND=wpfwo)       :: Tmeltbegin_g       ! temperature[K], above which graupel is assumed wet
    !REAL    (KIND=wpfwo)       :: Tmin_g             ! Temperature[K] at which smallest graupel particles melt
    REAL    (KIND=wpfwo)       :: meltdegTmin_g      ! degree of graupel melting at T=Tmin_g
    REAL    (KIND=wpfwo)       :: Tmax_min_g         ! minimum Tmax[K] of graupel in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: Tmax_max_g         ! maximum Tmax[K] of graupel in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: qthresh_g          ! q threshold to apply dynamic melting scheme for graupel
    REAL    (KIND=wpfwo)       :: qnthresh_g         ! qn threshold to apply dynamic melting scheme for graupel
    REAL    (KIND=wpfwo)       :: Tmeltbegin_h       ! temperature[K], above which hail is assumed wet
    !REAL    (KIND=wpfwo)       :: Tmin_h             ! Temperature[K] at which smallest hail particles melt
    REAL    (KIND=wpfwo)       :: meltdegTmin_h      ! degree of hail melting at T=Tmin_h
    REAL    (KIND=wpfwo)       :: Tmax_min_h         ! minimum Tmax[K] of hail in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: Tmax_max_h         ! maximum Tmax[K] of hail in dynamic melting scheme
    REAL    (KIND=wpfwo)       :: qthresh_h          ! q threshold to apply dynamic melting scheme for hail
    REAL    (KIND=wpfwo)       :: qnthresh_h         ! qn threshold to apply dynamic melting scheme for hail

    REAL    (KIND=wpfwo)       :: ext_tune_fac_pure  ! tuning factor for attenuation coefficients of droplets and dry ice species (ah_radar and adp_radar)
    REAL    (KIND=wpfwo)       :: ext_tune_fac_melt  ! tuning factor for attenuation coefficients of melting species (ah_radar and adp_radar)

    INTEGER                    :: station_id         ! unique 6-digit station ID
    INTEGER                    :: itype_refl         ! type of reflectivity calculation (Mie, Tmat or 3 different Rayleigh types)
    INTEGER                    :: isnow_type         ! type of (dry&wet) snow particle model for Mie/Tmat Scattering
    INTEGER                    :: igraupel_type      ! type of melting graupel particle model for Mie/Tmat Scattering
    INTEGER                    :: itype_Dref_fmelt   ! type of defining the Dref reference diameter for the melting degree parameterization
    INTEGER                    :: idummy1            ! dummy for later use, workaraound for effects of memory padding in MPI transfers

    ! Flags for only using certain hydrometeor types for dbz-caldulation:
    !  Order: elem(1) = cloud, (2) = rain, (3) = ice, (4) = snow, (5) = graupel, (6) = hail
    LOGICAL                    :: lhydrom_choice_testing(6) = .TRUE.
    ! Flag for using Mie lookup tables:
    LOGICAL                    :: llookup_mie
    ! if dynamic adaptation of Tmeltbegin_g/h and meltdegTmin_g/h should be applied within grid colums:
    LOGICAL                    :: ldynamic_wetgrowth_gh

    ! strings for defining the effective refractive index of hydrometeors:
    CHARACTER (LEN=12)         :: ctype_dryice_mie     ! dry ice, for Mie/Tmat calculations
    CHARACTER (LEN=12)         :: ctype_wetice_mie     ! wet ice, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_drysnow_mie    ! dry snow, for Mie/Tmat calculations
    CHARACTER (LEN=12)         :: ctype_wetsnow_mie    ! wet snow, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_drygraupel_mie ! dry graupel, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_wetgraupel_mie ! wet graupel, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_dryhail_mie    ! dry hail, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_wethail_mie    ! wet hail, for Mie/Tmat
    CHARACTER (LEN=12)         :: ctype_dryice_ray     ! dry ice, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_wetice_ray     ! wet ice, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_drysnow_ray    ! dry snow, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_wetsnow_ray    ! wet snow, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_drygraupel_ray ! dry graupel, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_wetgraupel_ray ! wet graupel, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_dryhail_ray    ! dry hail, for Rayleigh-theory itype_refl=2
    CHARACTER (LEN=12)         :: ctype_wethail_ray    ! wet hail, for Rayleigh-theory itype_refl=2

    ! polarimetry-relevant microphysics assumptions
    TYPE(t_polMP)              :: polMP_r            ! polarimetric microphysics structure for rain
    TYPE(t_polMP)              :: polMP_i            ! polarimetric microphysics structure for ice
    TYPE(t_polMP)              :: polMP_s            ! polarimetric microphysics structure for snow
    TYPE(t_polMP)              :: polMP_g            ! polarimetric microphysics structure for graupel
    TYPE(t_polMP)              :: polMP_h            ! polarimetric microphysics structure for hail

  END TYPE t_dbzcalc_params

  !.. Global default for namelist parameters for the dbz-calculation
  !   from radar_mie_lm_vec.f90:
  TYPE(t_dbzcalc_params), PARAMETER  :: dbz_namlst_d = t_dbzcalc_params ( &
       0.055_wpfwo    , & ! %lambda_radar         [m]
       273.15_wpfwo   , & ! %Tmeltbegin_i         [K]
       !T0C_fwo        , & ! %Tmin_i               [K]
       0.0_wpfwo      , & ! %meltdegTmin_i        [-]
       275.15_wpfwo   , & ! %Tmax_min_i           [K]
       278.15_wpfwo   , & ! %Tmax_max_i           [K]
       1e-8_wpfwo     , & ! %qthresh_i            [kg/kg]
       1e0_wpfwo      , & ! %qnthresh_i           [#/kg]
       273.15_wpfwo   , & ! %Tmeltbegin_s         [K]
       !T0C_fwo        , & ! %Tmin_s               [K]
       0.0_wpfwo      , & ! %meltdegTmin_s        [-]
       276.15_wpfwo   , & ! %Tmax_min_s           [K]
       283.15_wpfwo   , & ! %Tmax_max_s           [K]
       1e-8_wpfwo     , & ! %qthresh_s            [kg/kg]
       1e0_wpfwo      , & ! %qnthresh_s           [#/kg]
       263.15_wpfwo   , & ! %Tmeltbegin_g         [K]
       !T0C_fwo        , & ! %Tmin_g               [K]
       0.2_wpfwo      , & ! %meltdegTmin_g        [-]
       276.15_wpfwo   , & ! %Tmax_min_g           [K]
       288.15_wpfwo   , & ! %Tmax_max_g           [K]
       1e-8_wpfwo     , & ! %qthresh_g            [kg/kg]
       1e-3_wpfwo     , & ! %qnthresh_g           [#/kg]
       263.15_wpfwo   , & ! %Tmeltbegin_h         [K]
       !T0C_fwo        , & ! %Tmin_h               [K]
       0.2_wpfwo      , & ! %meltdegTmin_h        [-]
       278.15_wpfwo   , & ! %Tmax_min_h           [K]
       303.15_wpfwo   , & ! %Tmax_max_h           [K]
       1e-8_wpfwo     , & ! %qthresh_h            [kg/kg]
       1e-3_wpfwo     , & ! %qnthresh_h           [#/kg]
       1.0_wpfwo      , & ! %ext_tune_fac_pure    [-]
       1.0_wpfwo      , & ! %ext_tune_fac_melt    [-]
       999999         , & ! %station_id           [6-digit number]
       4              , & ! %itype_refl           [1, 2, 3, 4, 5, or 6] ! DO NOT CHANGE! Has to be 4 to enable the hosting model to use calc_dbz_vec() to silenty replace its own dBZ diagnostic
       1              , & ! %isnow_type           [1 or 2]
       1              , & ! %igraupel_type        [1, 2, or 3]
       1              , & ! %itype_Dref_fmelt     [1, or 2]
       0              , & ! %idummy1
       (/ .TRUE.      , & ! %lhydrom_choice_testing(1) = cloud drops
          .TRUE.      , & ! %lhydrom_choice_testing(2) = rain
          .TRUE.      , & ! %lhydrom_choice_testing(3) = cloud ice
          .TRUE.      , & ! %lhydrom_choice_testing(4) = snow
          .TRUE.      , & ! %lhydrom_choice_testing(5) = graupel
          .TRUE.  /)  , & ! %lhydrom_choice_testing(6) = hail
       .TRUE.         , & ! llookup_mie
       .FALSE.        , & ! ldynamic_wetgrowth_gh
       'mis         ' , & ! %ctype_dryice_mie     [12-char string]
       'mawsms      ' , & ! %ctype_wetice_mie     [12-char string]
       'masmas      ' , & ! %ctype_drysnow_mie    [12-char string]
       'mawsasmawsms' , & ! %ctype_wetsnow_mie    [12-char string]
       'mis         ' , & ! %ctype_drygraupel_mie [12-char string]
       'mawsms      ' , & ! %ctype_wetgraupel_mie [12-char string]
       'mis         ' , & ! %ctype_dryhail_mie    [12-char string]
       'mws         ' , & ! %ctype_wethail_mie    [12-char string]
       'mis         ' , & ! %ctype_dryice_ray     [12-char string]
       'mawsms      ' , & ! %ctype_wetice_ray     [12-char string]
       'mas         ' , & ! %ctype_drysnow_ray    [12-char string]
       'mawsms      ' , & ! %ctype_wetsnow_ray    [12-char string]
       'mis         ' , & ! %ctype_drygraupel_ray [12-char string]
       'mawsms      ' , & ! %ctype_wetgraupel_ray [12-char string]
       'mis         ' , & ! %ctype_dryhail_ray    [12-char string]
       'mawsms      ' , & ! %ctype_wethail_ray    [12-char string]
       ! emulate original standard EMVO behaviour (for itype_refl=5, otherwise irrelevant)
       t_polMP('R11     ', 0.0, 0.0, 0.0, -99.0), & ! t_polMP parameters for rain
       t_polMP('R11     ', 0.0, 0.0, 0.0, -99.0), & ! t_polMP parameters for ice
       t_polMP('R11     ', 0.0, 0.0, 0.0, -99.0), & ! t_polMP parameters for snow
       t_polMP('R11     ', 0.0, 0.0, 0.0, -99.0), & ! t_polMP parameters for graupel
       t_polMP('R11     ', 0.0, 0.0, 0.0, -99.0)  & ! t_polMP parameters for hail
       )


  !------------------------------------------------------------------------------


END MODULE radar_dbzcalc_params_type

