!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"

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

#if (defined TWOMOM_SB_NEW || defined TWOMOM_SB_OLD)
#define TWOMOM_SB
#endif

! .. Macro for indexing the Tmax_XX fields (horizontal 2D-fields) in (i,j,k) loops
!     where the outer k-loop goes over the 3rd dimension, j over the seconda and i over the first.
!    In ICON, the horizontal dimensions are the first and third, not the first and second!
#ifdef __ICON__
#define IND_IJ_2D  i,k
#else
#define IND_IJ_2D  i,j
#endif

#ifdef _DACE_
#define __COSMO__
#endif

MODULE radar_mie_iface_cosmo_1mom

!------------------------------------------------------------------------------
!
! Description:
!      Drivers for the interface functions for the COSMO/ICON microphysics
!      to the EMVORADO libraries for polarimetric radar moments computations.
!
!      Calculates polarization parameters of precipitation particles using
!      full Tmatrix theory. An option for simpler Mie theorie for reflectivity and
!      attenuation only is available. All methods include scattering by wet and melting particles.
!
!      Applicable to microwave radiation, wavelength > 1 mm
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:

!!! Now implemented: correct lookup tables for different Mie reflectivity configurations.
!!! Before, the configuration of the first call determined the lookup table, and all
!!! subsequent calls used this table, regardless of the true configuration.
!!!
!!! Rationale: New lookup tables are generated for all hydrometeor types, if the
!!! reflectivity configuration is different from the previous call. Here it can happen,
!!! that a new rain table is generated if a snow reflectivity parameter has changed.
!!! This is for simplicity reasons.
!!! In addition, these tables depend on other parameters unique to each hydrometeor
!!! type (table size, min/max diameter for spectral integration, N0, mu, etc.),
!!! and are recomputed, if something has changed here.
!!! To do this, a hash value is computed based on the overall reflectivity configuration and
!!! these other parameters for each hydrometeor type.
!!! The hash value is also used to generate the file names of binary files to store
!!! the lookup tables on disk. These are re-used at subsequent program calls or
!!! by other PEs in parallel runs for the same hash values to save time.
!!!
!!! NOTE 1: if you change some basic microphysics parameters like ageox, bgeox, n0s, n0g,
!!!       you have to delete old lookup tables on your disk, because these are no
!!!       longer correct. Exception: n0r and mu_rain are part of the hash value for rain
!!!       and can be changed without deleting the lookup table files!
!!!
!!! NOTE 2: if, in generating the lookup tables in parallel, not all PEs call the interface
!!!       routines, the program will hang up due to an exchange in initialize_tmax_1mom_vec_par().
!!!       To avoid this, set the computation neighbourhood "neigh_max" for the Tmax to a value <= 0.0!
!!!       In this case, there will be no exchange because the Tmax is determined
!!!       solely within the single grid columns. These calls should however not lead
!!!       to a valid model output, because Tmax is not determined as intended.


  !------------------------------------------------------------------------------

  USE radar_kind , ONLY :   &
       dp,          &  !           "--"           double precision
       wp              !           "--"           model working precision

  USE radar_dbzcalc_params_type, ONLY : t_polMP

  USE radar_data_mie, ONLY : &
       rho_0,         &  ! 1.225 kg/m^3
       inv_rhow, inv_rhow2, &
       inv_rhoi, inv_rhoi2, &
       T0C_fwo, mw_Tmin, mw_Tmax, mi_Tmin, mi_Tmax, Tmin_f, &
       rain, cloud, snow, ice, graupel, rain_coeffs, &
       q_crit_radar, quasi_zero, &
       ray_const, vt_const, t_mgd_params, n_stuetz, &
       Dmin_r, Dmax_r, Dmin_i, Dmax_i, Dmin_s, Dmax_s, Dmin_g, Dmax_g, &
       nmax_lookup, versionstring_lookup, &
       nq_r, nTa_r,                 nq_i, nTa_id, nTa_iw, nTm_i, &
       nq_s, nTa_sd, nTa_sw, nTm_s, nq_g, nTa_gd, nTa_gw, nTm_g, &
       t_dbzlookuptable, t_tabledef, &
       look_Z_rain, look_Z_ice, look_Z_snow, look_Z_graupel, &
       look_Z_meltice, look_Z_meltsnow, look_Z_meltgraupel, look_Z_meltgraupel_0C

  USE radar_mielib_vec, ONLY : &
       m_complex_water_ray, &
       m_complex_water_ray_vec, m_complex_ice_maetzler_vec

  USE radar_utilities, ONLY : tolower, toupper

  USE radar_mie_utils, ONLY : &
       particle_assign, &
       mgd_1mom, mgd_2mom, &
       nice_mono_1mom, snow_1mom_n0, &
       gamfac_imom_DMGD, &
       D_average_cin_DMGD, D_average_1M_exp_cin, D_of_X_average_1M_exp_cin, &
       hash_radar_lookup, get_hashbase, &
       decode_controlstring_3C, decode_controlstring_2C

  USE radar_mie_meltdegree, ONLY : &
       degree_of_melting_xxx_single, &
       degree_of_melting_xxx_single_Dref

!!$ This division of the radar_mie_specint is necessary because of makedepf90 automatic
!!$  dependency detection for the RADVOP_OFFLINE Makefiles. This tool does not like
!!$  too many continuation lines ...
  USE radar_mie_specint, ONLY : &
       zradar_rain_mie_vec , &
       zradar_ice_mie_vec , &
       zradar_snow_mie_vec , &
       zradar_graupel_mie_vec , &
       zradar_hail_mie_vec , &
       zradar_wetsnow_mie_vec , &
       zradar_wetice_mie_vec , &
       zradar_wetgr_mie_vec , &
       zradar_wetgr_twosph_mie_vec , &
       zradar_wetgr_wsph_mie_vec , &
       zradar_wethail_mie_vec , &
       zradar_rain_1mom_lookupcreate , &
       zradar_ice_1mom_lookupcreate , &
       zradar_snow_1mom_lookupcreate , &
       zradar_graupel_1mom_lookupcreate , &
       zradar_meltice_1mom_lookupcreate , &
       zradar_meltsnow_1mom_lookupcreate , &
       zradar_meltgraupel_1mom_lookupcreate , &
       zradar_triinterp_lookup_add_vec

  USE radar_mie_specint, ONLY : &
       zradar_rayleigh_L_x , &
       zradar_rayleigh_L_Lexp , &
       vtradar_pure_phase , &
       vtradar_mixed_phase , &
       nradar , &
       K_rho_fac_oguchi

  USE radar_mie_iface_cosmo_utils, ONLY : &
       init_radar_rayleigh_consts, init_radar_vt_oguchi, vtradar_normalize 

!===============================================================================
!===============================================================================

  IMPLICIT NONE

!===============================================================================
!==============================================================================

  PUBLIC

!==============================================================================

CONTAINS

!===============================================================================
!===============================================================================

  !==============================================================================
  !==============================================================================
  !
  ! SUBROUTINES FOR GRID POINT REFLECTIVITY/EXTINCTION CALCULATION FOR COSMO/ICON
  !
  !=======================================================================
  !
  ! Interface-Functions for calculation of radar parameters according to Mie or TMatrix
  ! for the 1-moment schemes, which are called in radar_mie_iface_cosmo_driver.f90. These routines use
  ! in turn many subroutines and functions which are defined in other modules.
  ! qx in kg/kg !!!
  ! lambda_radar is rounded to 3 significant digits !!!
  !=======================================================================

  SUBROUTINE radar_mie_1mom_vec(myproc, lambda_radar, itype_gscp_fwo, &
       itype_refl, luse_tmatrix, ldo_nonsphere, &
       isnow_n0temp, igraupel_type, itype_Dref_fmelt, &
       ctype_dryice, ctype_wetice, &
       ctype_drysnow, ctype_wetsnow, &
       ctype_drygraupel, ctype_wetgraupel, ldynamic_wetgrowth_gh, &
       Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
       Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
       Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g, &
       pMPr, pMPi, pMPs, pMPg, &
       rho, t, qc, qr, qi, qs, qg, &
       Tmax_i, Tmax_s, Tmax_g, &
       Tmin_g, &
       ilow, iup, jlow, jup, klow, kup, &
       lalloc_qi, lalloc_qs, lalloc_qg, llookup, &
       impipar_lookupgen, pe_start, pe_end,    &
       linterp_mode_dualpol,                   &
       ydir_lookup_read,ydir_lookup_write, &
       ext_tune_fac_pure, ext_tune_fac_melt, &
       zh_radar, ah_radar, &
       zv_radar, rrhv_radar, irhv_radar, kdp_radar, adp_radar, zvh_radar, &
       lhydrom_choice_testing )

    ! ctype_dryice: 3 chars
    ! ctype_wetice: 6 chars
    ! ctype_drysnow: 3 chars
    ! ctype_wetsnow: 6 chars
    ! ctype_drygraupel: 3 chars
    ! ctype_wetgraupel: 6 chars

    IMPLICIT NONE

    INTEGER,                         INTENT(IN) :: itype_gscp_fwo, igraupel_type, isnow_n0temp, &
                                                   itype_Dref_fmelt, itype_refl

    INTEGER,                         INTENT(IN) :: myproc, ilow, iup, jlow, jup, klow, kup

    REAL(KIND=dp),                   INTENT(IN) :: lambda_radar, &
                                                   Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
                                                   Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
                                                   Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g

    TYPE(t_polMP), INTENT(in)                     :: pMPr, pMPi, pMPs, pMPg

    ! Model fields in model's wp, not dp:
    REAL(KIND=wp), DIMENSION(:,:,:), INTENT(IN) :: rho, t, qc, qr, qi, qs, qg    ! In model wp, not dp
    ! Tmax for degree of melting in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmax_i, Tmax_s, Tmax_g
    ! Tmin for dynamic wet growth of graupel in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmin_g

    CHARACTER(len=*), INTENT(in)                :: ctype_dryice, ctype_wetice, &
                                                   ctype_drysnow, ctype_wetsnow, &
                                                   ctype_drygraupel, ctype_wetgraupel

    LOGICAL, INTENT(in)                         :: llookup, &
                                                   lalloc_qi, lalloc_qs, lalloc_qg, &
                                                   luse_tmatrix, ldo_nonsphere, ldynamic_wetgrowth_gh
    INTEGER, INTENT(in)                         :: impipar_lookupgen, pe_start, pe_end
    LOGICAL, INTENT(in)                         :: linterp_mode_dualpol
    CHARACTER(len=*), INTENT(in)                :: ydir_lookup_read, ydir_lookup_write
    LOGICAL, INTENT(in), OPTIONAL               :: lhydrom_choice_testing(6)

    REAL(KIND=dp),                   INTENT(IN) :: ext_tune_fac_pure, & ! handwaving tuning factors for simulated attenuation coefficients,
                                                   ext_tune_fac_melt    ! to be able to mimick the behaviour of an imperfect, conservative attenuation correction

    REAL(KIND=dp), INTENT(OUT), OPTIONAL        :: zh_radar(:,:,:), ah_radar(:,:,:), &
                                                   zv_radar(:,:,:), &
                                                   rrhv_radar(:,:,:), irhv_radar(:,:,:), &
                                                   kdp_radar(:,:,:), adp_radar(:,:,:), &
                                                   zvh_radar(:,:,:)

    ! Local Variables
    !----------------

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: zh_local, ah_local, &
                                                    zv_local, &
                                                    rrhv_local, irhv_local, &
                                                    kdp_local, adp_local, &
                                                    zvh_local, &
                                                    n0_s

    INTEGER :: i, j, k, ku, ko, ni, nj, nk

    REAL(KIND=dp) :: T_a, T_m,         &
                     fmelt_mean, Dtmp, &
                     q_c,              &
                     q_r,              &
                     q_i, n_i, x_i,    &
                     q_s,              &
                     q_g,              &
                     hlp, &
                     Tmeltbegin_g_loc, meltdegTmin_g_loc  ! for ldynamic_wetgrowth_gh

    REAL(KIND=dp) :: K_w2, K_i2, K2divRho2
    REAL(KIND=dp), PARAMETER :: Dexpo_i = 0.5d0, Dref_i = 100d-6

    REAL(KIND=dp) :: zh, ah, &
                     zv, rrhv, irhv, kdp, adp, zvh

    TYPE(t_tabledef)   :: tableprops
    TYPE(t_mgd_params) :: mgd

    COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: m_i, m_w
    COMPLEX(kind=dp) :: Kw_a, Ki_a

    CHARACTER(len=20) :: &
         mixingrulestring_i_dry, matrixstring_i_dry, inclusionstring_i_dry, &
         mixingrulestring_i, matrixstring_i, inclusionstring_i, &
         hoststring_i, hostmatrixstring_i, hostinclusionstring_i, &
         mixingrulestring_s_shell, matrixstring_s_shell, inclusionstring_s_shell, &
         mixingrulestring_s_core, matrixstring_s_core, inclusionstring_s_core, &
         hoststring_s_shell, hostmatrixstring_s_shell, hostinclusionstring_s_shell, &
         hoststring_s_core, hostmatrixstring_s_core, hostinclusionstring_s_core, &
         mixingrulestring_s_shell_dry, matrixstring_s_shell_dry, inclusionstring_s_shell_dry, &
         mixingrulestring_s_core_dry, matrixstring_s_core_dry, inclusionstring_s_core_dry, &
         mixingrulestring_g, matrixstring_g, inclusionstring_g, &
         mixingrulestring_g_shell, matrixstring_g_shell, inclusionstring_g_shell, &
         hoststring_g, hostmatrixstring_g, hostinclusionstring_g, &
         mixingrulestring_g_dry, matrixstring_g_dry, inclusionstring_g_dry

    LOGICAL, PARAMETER :: mur_is_fixed = .TRUE., n0r_is_fixed = .TRUE.

    ! Pointers to single elements of the lookup table type vectors:
    !  (are necessary as input parameters to the lookup table interpolation
    !   functions to enable correct inlining on the NEC)
    TYPE(t_dbzlookuptable), POINTER :: look_Z_rain_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_ice_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_snow_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_graupel_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_meltice_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_meltsnow_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_meltgraupel_p
    TYPE(t_dbzlookuptable), POINTER :: look_Z_meltgraupel_0C_p

    ! counter for the actual number of wavelengths for which Mie lookup tables exist:
    INTEGER, SAVE              :: nlkup_rain     = 0
    INTEGER, SAVE              :: nlkup_ice      = 0
    INTEGER, SAVE              :: nlkup_meltice  = 0
    INTEGER, SAVE              :: nlkup_snow     = 0
    INTEGER, SAVE              :: nlkup_meltsnow = 0
    INTEGER, SAVE              :: nlkup_grau     = 0
    INTEGER, SAVE              :: nlkup_meltgrau = 0
    INTEGER, SAVE              :: nlkup_meltgrau_0C = 0
    INTEGER                    :: ilkup_rain
    INTEGER                    :: ilkup_ice
    INTEGER                    :: ilkup_meltice
    INTEGER                    :: ilkup_snow
    INTEGER                    :: ilkup_meltsnow
    INTEGER                    :: ilkup_grau
    INTEGER                    :: ilkup_meltgrau
    INTEGER                    :: ilkup_meltgrau_0C
    INTEGER                    :: magicnr_rain, &
                                  magicnr_ice,  magicnr_meltice, &
                                  magicnr_snow, magicnr_meltsnow, &
                                  magicnr_grau, magicnr_meltgrau, magicnr_meltgrau_0C
    INTEGER, PARAMETER         :: unit_lookup = 2430
    REAL(KIND=dp)              :: lambda_radar_3digits
    CHARACTER(len=12)          :: clambda_radar

    CHARACTER(len=500)         :: chydroconfig, ctableconfig, chydroconfig_pMPrain
    CHARACTER(len=3000)        :: hashtext, cbaseconfig, cbaseconfig_g_0C

    LOGICAL :: ldo_qc, ldo_qr, ldo_qi, ldo_qi_rayleigh, ldo_qs, ldo_qg

    REAL(kind=dp), DIMENSION(SIZE(rho, DIM=1)) :: q_x_vec, T_a_vec, T_a_wg_vec, T_m_vec, T0C_vec, a_w_vec, b_w_vec
    LOGICAL, DIMENSION(SIZE(rho, DIM=1)) :: flag_interp, is_melting, T_gt_0C

    CHARACTER (len=*), PARAMETER :: yzroutine = 'emvorado::radar_mie_1mom_vec()'

    ! .. Round lambda_radar to 3 significant digits, to avoid differences to slightly differing numeric
    !     lambda_radar values from namelists and DATA files:
    WRITE(clambda_radar, '(es12.2)') lambda_radar
    READ (clambda_radar, *)          lambda_radar_3digits
    
    ! .. Get field dimensions from input rho (assume that all hydrometeors and t are the same size):
    ni = SIZE(rho, DIM=1)
    nj = SIZE(rho, DIM=2)
    nk = SIZE(rho, DIM=3)

    ! .. Check field dimensions of output fields:
    IF (PRESENT(zh_radar)) THEN
      IF (ni /= SIZE(zh_radar, dim=1) .OR. nj /= SIZE(zh_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and zh_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(ah_radar)) THEN
      IF (ni /= SIZE(ah_radar, dim=1) .OR. nj /= SIZE(ah_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and ah_radar!'
        STOP
      END IF
    END IF

    IF (PRESENT(zv_radar)) THEN
      IF (ni /= SIZE(zv_radar, dim=1) .OR. nj /= SIZE(zv_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and zv_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(rrhv_radar)) THEN
      IF (ni /= SIZE(rrhv_radar, dim=1) .OR. nj /= SIZE(rrhv_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and rrhv_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(irhv_radar)) THEN
      IF (ni /= SIZE(irhv_radar, dim=1) .OR. nj /= SIZE(irhv_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and irhv_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(kdp_radar)) THEN
      IF (ni /= SIZE(kdp_radar, dim=1) .OR. nj /= SIZE(kdp_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and kdp_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(adp_radar)) THEN
      IF (ni /= SIZE(adp_radar, dim=1) .OR. nj /= SIZE(adp_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and adp_radar!'
        STOP
      END IF
    END IF
    IF (PRESENT(zvh_radar)) THEN
      IF (ni /= SIZE(zvh_radar, dim=1) .OR. nj /= SIZE(zvh_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                ': Inconsistent field dimensions of rho and zvh_radar!'
        STOP
      END IF
    END IF

    ! .. Allocate local work arrays:
    ALLOCATE(zh_local(ni,nj,nk), ah_local(ni,nj,nk))
    ALLOCATE(zv_local(ni,nj,nk), &
             rrhv_local(ni,nj,nk), irhv_local(ni,nj,nk), &
             kdp_local(ni,nj,nk), adp_local(ni,nj,nk))
    ALLOCATE(zvh_local(ni,nj,nk))
    ALLOCATE(n0_s(ni,nj,nk))
    ALLOCATE(m_w(ni,nj,nk), m_i(ni,nj,nk))

    ! Generate full control strings for the effective refractive index from input shorthand strings:
    CALL decode_controlstring_2C(ctype_dryice(1:3), &
         mixingrulestring_i_dry, matrixstring_i_dry, inclusionstring_i_dry, &
         'mas')
    CALL decode_controlstring_3C(ctype_wetice(1:6), &
         mixingrulestring_i, matrixstring_i, inclusionstring_i, &
         hoststring_i, hostmatrixstring_i, hostinclusionstring_i, &
         'mawsas')
    CALL decode_controlstring_2C(ctype_drysnow(1:3), &
         mixingrulestring_s_core_dry, matrixstring_s_core_dry, inclusionstring_s_core_dry, &
         'mas')
    CALL decode_controlstring_2C(ctype_drysnow(4:6), &
         mixingrulestring_s_shell_dry, matrixstring_s_shell_dry, inclusionstring_s_shell_dry, &
         'mas')
    CALL decode_controlstring_3C(ctype_wetsnow(1:6), &
         mixingrulestring_s_core, matrixstring_s_core, inclusionstring_s_core, &
         hoststring_s_core, hostmatrixstring_s_core, hostinclusionstring_s_core, &
         'mawsas')
    CALL decode_controlstring_3C(ctype_wetsnow(7:12), &
         mixingrulestring_s_shell, matrixstring_s_shell, inclusionstring_s_shell, &
         hoststring_s_shell, hostmatrixstring_s_shell, hostinclusionstring_s_shell, &
         'mawsas')
    CALL decode_controlstring_2C(ctype_drygraupel(1:3), &
         mixingrulestring_g_dry, matrixstring_g_dry, inclusionstring_g_dry, &
         'mas')
    SELECT CASE (igraupel_type)
    CASE(1)
      CALL decode_controlstring_3C(ctype_wetgraupel(1:6), &
           mixingrulestring_g, matrixstring_g, inclusionstring_g, &
           hoststring_g, hostmatrixstring_g, hostinclusionstring_g, &
           'mawsas')
      mixingrulestring_g_shell(:) = ' '
      matrixstring_g_shell(:)     = ' '
      inclusionstring_g_shell(:)  = ' '
    CASE(2)
      CALL decode_controlstring_2C(ctype_wetgraupel(1:3), &
           mixingrulestring_g, matrixstring_g, inclusionstring_g, &
           'mis')
      CALL decode_controlstring_2C(ctype_wetgraupel(4:6), &
           mixingrulestring_g_shell, matrixstring_g_shell, inclusionstring_g_shell, &
           'mis')
    CASE(3)
      CALL decode_controlstring_2C(ctype_wetgraupel(1:3), &
           mixingrulestring_g, matrixstring_g, inclusionstring_g, &
           'mis')
      mixingrulestring_g_shell(:) = ' '
      matrixstring_g_shell(:)     = ' '
      inclusionstring_g_shell(:)  = ' '
    END SELECT


    zh_local   = 0d0
    ah_local   = 0d0
    zv_local   = 0d0
    rrhv_local = 0d0
    irhv_local = 0d0
    kdp_local  = 0d0
    adp_local  = 0d0
    zvh_local  = 0d0

    IF (PRESENT(zh_radar) .OR. PRESENT(ah_radar)) THEN
      ku = klow
      ko = kup
    ELSE
      ! Simply do nothing
      ku = 2
      ko = 1
    END IF

    IF (PRESENT(lhydrom_choice_testing)) THEN
      ldo_qc = lhydrom_choice_testing(1)
      ldo_qr = lhydrom_choice_testing(2)
      ldo_qi = lhydrom_choice_testing(3) .AND. lalloc_qi
      ldo_qs = lhydrom_choice_testing(4) .AND. lalloc_qs
      ldo_qg = lhydrom_choice_testing(5) .AND. lalloc_qg
    ELSE
      ldo_qc = .TRUE.
      ldo_qr = .TRUE.
      ldo_qi = lalloc_qi
      ldo_qs = lalloc_qs
      ldo_qg = lalloc_qg
    END IF

    ! On its first call or if otherwise necessary, this computes
    ! some constant prefactors needed for cloud drops and cloud ice
    ! (Rayleigh-Oguchi) below
    ! and stores it in global struct "ray_const":
    ! (The values for rain, snow and graupel are not used below)
    CALL init_radar_rayleigh_consts( &
         .TRUE., .FALSE., .TRUE., .FALSE., .FALSE., .FALSE., &
         DBLE(lambda_radar_3digits))

    ! for ice use Rayleigh except for itype_refl=5
    ldo_qi_rayleigh = .NOT.(luse_tmatrix .AND. ldo_nonsphere)

!$OMP PARALLEL DO PRIVATE(j,k)
    DO k = 1, nk
      DO j = 1, nj
        ! Refractive index of water after Ray (1972), only valid for -10 < T < 30 degree C
        m_w(:,j,k) = m_complex_water_ray_vec(DBLE(lambda_radar_3digits),&
           MIN(MAX(t(:,j,k)-T0C_fwo,mw_Tmin),mw_Tmax), ni)
        ! Refractive index of ice after C. Maetzler, only valid for -80 < T < 0 degree C
        m_i(:,j,k) = m_complex_ice_maetzler_vec(DBLE(lambda_radar_3digits),&
           MIN(MAX(t(:,j,k)-T0C_fwo,mi_Tmin),mi_Tmax), ni)
      ENDDO
    ENDDO
!$omp end parallel do

    ! Precalculate NWP-field dependent snow number density field (but only if we
    ! we consider snow at all):
    IF (ldo_qs) &
      CALL snow_1mom_n0(itype_gscp_fwo, isnow_n0temp, &
                        t, rho, qs, snow%a_geo, &
                        n0_s)


    ! .. Reflectivity of cloud droplets and cloud ice with simple Rayleigh-Oguchi-Theory:
!$omp parallel private (i,j,k,T_a,hlp,q_i,n_i,x_i,q_c,&
!$omp&                  Kw_a,K_w2,Ki_a,K_i2,K2divRho2,fmelt_mean,Dtmp,&
!$omp&                  zh,ah)
!$omp do
    DO k = ku, ko
      DO j = jlow, jup
        DO i = ilow, iup

          T_a = t(i,j,k)
          hlp = rho(i,j,k)

          Kw_a = (m_w(i,j,k)**2-1d0)/(m_w(i,j,k)**2+2d0)
          K_w2 = ABS(Kw_a)**2
          Ki_a = (m_i(i,j,k)**2-1d0)/(m_i(i,j,k)**2+2d0)
          K_i2 = ABS(Ki_a)**2

          ! Reflectivity of cloud droplets after Rayleigh-Approximation:
          IF (ldo_qc) THEN

            q_c = qc(i,j,k) * hlp
            IF (q_c < q_crit_radar%rayleigh) q_c = 0d0

            ! JM191030:
            ! set all relevant polarimetric parameters, in the absence of better
            ! assumptions according to relations valid for spheres
            zh = zradar_rayleigh_L_x(q_c, cloud%x_max, K_w2, ray_const%cloud_Zprefac_rho)
            zh_local(i,j,k) = zh_local(i,j,k) + zh
            zv_local(i,j,k) = zv_local(i,j,k) + zh
            ! JM201020: FIXME
            ! Correct to only set rrhv? Basically, for sh==sv: rhv=|rsh^2+ish^2|,
            ! while zh=|rsh+j*ish|^2. Is that indeed identical, hence ok to leave
            ! irhv= 0? Might be for rhv, but probably not for dhv...
            rrhv_local(i,j,k) = rrhv_local(i,j,k) + zh

            ! Rayleigh-Approximation for cloud droplet extinction (only the dominating absorption term):
            ah = ray_const%cloud_extprefac * AIMAG(-Kw_a) * q_c
            ah_local(i,j,k) = ah_local(i,j,k) + ah*ext_tune_fac_pure

          END IF

          ! Reflectivity of cloud ice after Rayleigh:
          IF (ldo_qi .AND. ldo_qi_rayleigh) THEN

            q_i = qi(i,j,k) * hlp
            IF (q_i < q_crit_radar%rayleigh) q_i = 0d0

            n_i = nice_mono_1mom(T_a)
            x_i = q_i / n_i
            x_i = MIN( MAX (x_i, ice%x_min), ice%x_max)

            IF (T_a > Tmeltbegin_i) THEN
              ! m, diameter of particle for calculation of degree of melting
              Dtmp = (x_i/ice%a_geo)**(1d0/ice%b_geo) * 2d0
              CALL degree_of_melting_xxx_single_Dref(&
                      T_a,Dtmp,Dexpo_i,Dref_i,&
                      Tmeltbegin_i,meltdegTmin_i,Tmin_f,&
                      REAL(Tmax_i(IND_IJ_2D),kind=dp),&
                      fmelt_mean)
              K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,&
                      MAX(MIN(fmelt_mean,1d0),0d0))
              zh = zradar_rayleigh_L_x(q_i, x_i, K2divRho2, ray_const%ice1mom_Zprefac)
            ELSE
              zh = zradar_rayleigh_L_x(q_i, x_i, K_i2, ray_const%ice1mom_Zprefac_rho)
            END IF
            ! JM191030:
            ! set all relevant polarimetric parameters, in the absence of better
            ! assumptions according to relations valid for spheres
            zh_local(i,j,k) = zh_local(i,j,k) + zh
            zv_local(i,j,k) = zv_local(i,j,k) + zh
            ! JM201020: FIXME
            ! Correct to only set rrhv? Basically, for sh==sv: rhv=|rsh^2+ish^2|,
            ! while zh=|rsh+j*ish|^2. Is that indeed identical, hence ok to leave
            ! irhv= 0? Might be for rhv, but probably not for dhv...
            rrhv_local(i,j,k) = rrhv_local(i,j,k) + zh

            ! JM191030:
            ! no attenuation at all from cloud ice, not even from melting???

          END IF

        END DO
      END DO
    END DO
!$omp end do
!$omp end parallel

    !==================================================================================
    !
    ! .. Reflectivity of rain, snow and graupel with full Mie/TMatrix-Theory:
    !    (Integration with respect to particle diameter)
    !
    !==================================================================================

    ! .. Integration limits for the integrals of cross sections over the PSDs:
    !    JM210712: Moved to radar_data_mie

    !Dmin_r =  50.0d-6 ! lower limit of integration of reflectivity: 0.05 mm
    !Dmax_r =  10.0d-3 ! upper limit: 10.0 mm

    !Dmin_i =   5.0d-6 ! lower limit of integration of reflectivity: 0.005 mm
    !Dmax_i =  10.0d-3 ! upper limit: 10.0 mm

    !Dmin_s =  50.0d-6 ! lower limit of integration of reflectivity: 0.05 mm
    !Dmax_s =  50.0d-3 ! upper limit: 50.0 mm

    !Dmin_g =  10.0d-6 ! lower limit of integration of reflectivity: 0.01 mm
    !Dmax_g =  30.0d-3 ! upper limit: 30.0 mm

    IF (llookup) THEN

      ! .. Creating lookup tables and calculating reflectivity:

      ! version number for all lookup tables (rain, ice, snow, graupel):
      ! -----------------------------------------------------------
      !
      ! NOTE: adapt "versionstring_lookup" (global module variable, defined in
      !       radar_data_mie) whenever fundamental things of the table are
      !       changed such as:
      !  - Number of nodes and scaling of the table vectors
      !  - Particle models for melting stuff
      !  - Hydrometeor shape and orientation parametrizations
      !  - Adaptions/bugfixes in the basic Mie or TMatrix codes
      !  - Format changes of the lookup table files

      ! magic number to encrypt the actual overall parameter set:
      ! ---------------------------------------------------------
      !
      ! NOTE:
      !  - to keep things simple, we take the entire reflectivity configuration into account,
      !    not only the respective hydrometeor class.
      !    This can lead to the situation, that a change e.g. in the snow configuration might
      !    lead to re-computation of a new table e.g. for rain, although no changes were made
      !    to rain.
      !    This does not do much harm - it takes extra time, but better safe than sorry.

      cbaseconfig = get_hashbase (versionstring_lookup, lambda_radar_3digits, itype_gscp_fwo, &
                                  mur_is_fixed, isnow_n0temp, &
                                  itype_refl, igraupel_type, itype_Dref_fmelt,&
                                  Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
                                  Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
                                  Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g, &
                                  0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                  ctype_dryice, ctype_wetice, &
                                  ctype_drysnow, ctype_wetsnow, &
                                  ctype_drygraupel, ctype_wetgraupel, &
                                  '', '')
      
      IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN
        cbaseconfig_g_0C = get_hashbase (versionstring_lookup, lambda_radar_3digits, itype_gscp_fwo, &
                                         mur_is_fixed, isnow_n0temp, &
                                         itype_refl, igraupel_type, itype_Dref_fmelt,&
                                         Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
                                         Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
                                         T0C_fwo, 0.0_dp,             Tmax_min_g, Tmax_max_g, &
                                         0.0_dp,0.0_dp,0.0_dp,0.0_dp, &
                                         ctype_dryice, ctype_wetice, &
                                         ctype_drysnow, ctype_wetsnow, &
                                         ctype_drygraupel, ctype_wetgraupel, &
                                         '', '')
      END IF
      
      !  - get rain polMP properties once; will later be needed for melt* hashes in addition
      !    to rain itself
      IF (luse_tmatrix) THEN
        IF (toupper(TRIM(ADJUSTL(pMPr%ARmodel))) == 'POLY') THEN
          WRITE (chydroconfig_pMPrain, '(A9,3es13.3,F6.1)') &
                 toupper(TRIM(ADJUSTL(pMPr%ARmodel))), &
                 pMPr%c0, pMPr%c1, pMPr%ARmin, pMPr%sig
        ELSE
          WRITE (chydroconfig_pMPrain, '(A9,F6.1)') &
                 toupper(TRIM(ADJUSTL(pMPr%ARmodel))), pMPr%sig
        END IF
      ELSE
        chydroconfig_pMPrain(:) = ' '
      END IF

      !  - we add parameters which are not part of the reflectivity calculation configuration,
      !    but might affect the table contents, e.g. hydrometeor microphysics and table
      !    dimension info. these are hy<drometeor class specific, henc prepared within each
      !    hydrometeor class loop separately.


      ! JM210712:
      ! definition of table sizes now in radar_data_mie.f90

      IF (ldo_qr) THEN

       ! JM200128:
       ! currently mur_is_fixed and n0r_is_fixed are always true in 1mom scheme
       IF (mur_is_fixed .AND. n0r_is_fixed) THEN

        chydroconfig(:) = ' '
        WRITE (chydroconfig, '(7es16.6,I5,A)') &
               rain%nu, rain%mu, rain%n0_const, rain%a_geo, rain%b_geo, &
               Dmin_r, Dmax_r, n_stuetz, &
               TRIM(chydroconfig_pMPrain)

        ! create lookup table(s) for rain:
        !======================================

        ! table sizes: parameters nTa_r and nq_r now defined in radar_data_mie.f90

        tableprops%nqi = nq_r
        tableprops%nTa = nTa_r
        tableprops%nTm = 0
        tableprops%qilow = q_crit_radar%liquid
        tableprops%qiup  = 0.1_dp
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(2I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_rain, hashtext)

        ! Search lookup table index for rain for the actual parameter set:
        ilkup_rain = -99
        DO i = 1, nlkup_rain
          IF ( look_Z_rain(i)%magicnr == magicnr_rain ) THEN
            ilkup_rain = i
            EXIT
          END IF
        END DO

        IF (ilkup_rain < 0) THEN

          nlkup_rain = nlkup_rain + 1
          IF (nlkup_rain > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for rain! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_rain(nlkup_rain)%magicnr         = magicnr_rain
          look_Z_rain(nlkup_rain)%cversion_lt     = versionstring_lookup
          look_Z_rain(nlkup_rain)%luse_tmatrix    = luse_tmatrix 
          look_Z_rain(nlkup_rain)%ldo_nonsphere   = ldo_nonsphere
          look_Z_rain(nlkup_rain)%itype_refl      = itype_refl   
          look_Z_rain(nlkup_rain)%chydrotype(:)   = ' '
          look_Z_rain(nlkup_rain)%chydrotype      = 'rain-1mom'
          look_Z_rain(nlkup_rain)%chydroconfig(:) = ' '
          look_Z_rain(nlkup_rain)%chydroconfig    = TRIM(chydroconfig)
         
          ilkup_rain = nlkup_rain

          IF ( look_Z_rain(ilkup_rain)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_rain(ilkup_rain) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for rain (only if it has not yet been created for this wavelength)
        IF (.NOT. look_Z_rain(ilkup_rain)%is_initialized) THEN
          CALL zradar_rain_1mom_lookupcreate(&
               look_Z_rain(ilkup_rain),tableprops,&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPr,&
               ray_const%Z_fac,Dmin_r,Dmax_r,&
               impipar_lookupgen, pe_start, pe_end,&
               linterp_mode_dualpol,&
               ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)
        END IF

       ELSE

        WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                ': use of rain Mie lookup tables only '// &
                'implemented for constant mu and constant n0! Stop!'
        STOP

       END IF

      END IF

      IF (ldo_qi .AND. .NOT.ldo_qi_rayleigh) THEN

        chydroconfig(:) = ' '
        IF (luse_tmatrix) THEN
          IF (toupper(TRIM(ADJUSTL(pMPi%ARmodel))) == 'POLY') THEN
            WRITE (chydroconfig, '(7es16.6,I5,A9,3es13.3,F6.1)') &
                   ice%nu, ice%mu, ice%n0_const, ice%a_geo, ice%b_geo, &
                   Dmin_i, Dmax_i, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPi%ARmodel))), &
                   pMPi%c0, pMPi%c1, pMPi%ARmin, pMPi%sig
          ELSE
            WRITE (chydroconfig, '(7es16.6,I5,A9,F6.1)') &
                   ice%nu, ice%mu, ice%n0_const, ice%a_geo, ice%b_geo, &
                   Dmin_i, Dmax_i, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPi%ARmodel))), pMPi%sig
          END IF
        ELSE
          WRITE (chydroconfig, '(7es16.6,I5)') &
                 ice%nu, ice%mu, ice%n0_const, ice%a_geo, ice%b_geo, &
                 Dmin_i, Dmax_i, n_stuetz
        END IF

        ! create lookup table for cloud ice:
        !======================================

        ! table sizes: nTa_id, nq_i in radar_data_mie.f90

        tableprops%nqi = nq_i
        tableprops%nTa = nTa_id
        tableprops%nTm = 0
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.0_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(2I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_ice, hashtext)

        ! Search lookup table index for dry ice for the actual parameter set:
        ilkup_ice = -99
        DO i = 1, nlkup_ice
          IF ( look_Z_ice(i)%magicnr == magicnr_ice ) THEN
            ilkup_ice = i
            EXIT
          END IF
        END DO

        IF (ilkup_ice < 0) THEN

          nlkup_ice = nlkup_ice + 1
          IF (nlkup_ice > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for dry cloud ice! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_ice(nlkup_ice)%magicnr         = magicnr_ice
          look_Z_ice(nlkup_ice)%cversion_lt     = versionstring_lookup
          look_Z_ice(nlkup_ice)%luse_tmatrix    = luse_tmatrix 
          look_Z_ice(nlkup_ice)%ldo_nonsphere   = ldo_nonsphere
          look_Z_ice(nlkup_ice)%itype_refl      = itype_refl   
          look_Z_ice(nlkup_ice)%chydrotype(:)   = ' '
          look_Z_ice(nlkup_ice)%chydrotype      = 'dryice-1mom'
          look_Z_ice(nlkup_ice)%chydroconfig(:) = ' '
          look_Z_ice(nlkup_ice)%chydroconfig    = TRIM(chydroconfig)
          
          ilkup_ice = nlkup_ice

          IF ( look_Z_ice(ilkup_ice)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_ice(ilkup_ice) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for dry ice (only if it has not yet been created)
        IF (.NOT. look_Z_ice(ilkup_ice)%is_initialized ) THEN

          CALL zradar_ice_1mom_lookupcreate(&
               look_Z_ice(ilkup_ice),tableprops,&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPi,&
               ray_const%Z_fac,Dmin_i,Dmax_i,ice,&
               Tmeltbegin_i,mixingrulestring_i_dry,matrixstring_i_dry,inclusionstring_i_dry,&
               impipar_lookupgen,pe_start,pe_end,&
               linterp_mode_dualpol,&
               ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)

        END IF

        ! create lookup table for melting cloud ice:
        !===========================================

        ! melt* tables also depend on rain polMP properties, hence include them
        !   in the melt* hashes in addition to the ice ones
        WRITE (chydroconfig, '(2A)') &
               TRIM(chydroconfig), TRIM(chydroconfig_pMPrain)

        ! table sizes: nTa_iw, nTm_i, nq_i in radar_data_mie.f90
        tableprops%nqi = nq_i
        tableprops%nTa = nTa_iw
        tableprops%nTm = nTm_i
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.5_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(3I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, tableprops%nTm, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_meltice, hashtext)

        ! Search lookup table index for melting ice for the actual parameter set:
        ilkup_meltice = -99
        DO i = 1, nlkup_meltice
          IF ( look_Z_meltice(i)%magicnr == magicnr_meltice ) THEN
            ilkup_meltice = i
            EXIT
          END IF
        END DO

        IF (ilkup_meltice < 0) THEN

          nlkup_meltice = nlkup_meltice + 1
          IF (nlkup_meltice > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for melting cloud ice! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_meltice(nlkup_meltice)%magicnr          = magicnr_meltice
          look_Z_meltice(nlkup_meltice)%cversion_lt      = versionstring_lookup
          look_Z_meltice(nlkup_meltice)%luse_tmatrix     = luse_tmatrix    
          look_Z_meltice(nlkup_meltice)%ldo_nonsphere    = ldo_nonsphere   
          look_Z_meltice(nlkup_meltice)%itype_refl       = itype_refl      
          look_Z_meltice(nlkup_meltice)%itype_Dref_fmelt = itype_Dref_fmelt
          look_Z_meltice(nlkup_meltice)%chydrotype(:)    = ' '
          look_Z_meltice(nlkup_meltice)%chydrotype       = 'meltice-1mom'
          look_Z_meltice(nlkup_meltice)%chydroconfig(:)  = ' '
          look_Z_meltice(nlkup_meltice)%chydroconfig     = TRIM(chydroconfig)
          
          ilkup_meltice = nlkup_meltice

          IF ( look_Z_meltice(ilkup_meltice)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_melticeel(ilkup_meltice) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for melting ice (only if it has not yet been created)
        IF (.NOT. look_Z_meltice(ilkup_meltice)%is_initialized) THEN

          CALL zradar_meltice_1mom_lookupcreate( &
               look_Z_meltice(ilkup_meltice),tableprops,&
               itype_Dref_fmelt,&
               DBLE(Tmeltbegin_i),DBLE(meltdegTmin_i),&
               DBLE(Tmax_min_i),DBLE(Tmax_max_i),&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPi,pMPr,&
               ray_const%Z_fac,Dmin_i,Dmax_i,ice,&
               mixingrulestring_i,matrixstring_i,inclusionstring_i,&
               hoststring_i,hostmatrixstring_i,hostinclusionstring_i,&
               impipar_lookupgen,pe_start,pe_end,&
               linterp_mode_dualpol,&
               ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)

        END IF

      END IF ! ice

      IF (ldo_qs) THEN

        chydroconfig(:) = ' '
        IF (luse_tmatrix) THEN
          IF (toupper(TRIM(ADJUSTL(pMPs%ARmodel))) == 'POLY') THEN
            WRITE (chydroconfig, '(7es16.6,I5,A9,3es13.3,F6.1)') &
                   snow%nu, snow%mu, snow%n0_const, snow%a_geo, snow%b_geo, &
                   Dmin_s, Dmax_s, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPs%ARmodel))), &
                   pMPs%c0, pMPs%c1, pMPs%ARmin, pMPs%sig
          ELSE
            WRITE (chydroconfig, '(7es16.6,I5,A9,F6.1)') &
                   snow%nu, snow%mu, snow%n0_const, snow%a_geo, snow%b_geo, &
                   Dmin_s, Dmax_s, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPs%ARmodel))), pMPs%sig
          END IF
        ELSE
          WRITE (chydroconfig, '(7es16.6,I5)') &
                 snow%nu, snow%mu, snow%n0_const, snow%a_geo, snow%b_geo, &
                 Dmin_s, Dmax_s, n_stuetz
        END IF

        ! create lookup table for dry snow:
        !======================================

        ! table sizes: nTa_sd, nq_s in radar_data_mie.f90

        tableprops%nqi = nq_s
        tableprops%nTa = nTa_sd
        tableprops%nTm = 0
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.5_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(2I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_snow, hashtext)

        ! Search lookup table index for dry snow for the actual parameter set:
        ilkup_snow = -99
        DO i = 1, nlkup_snow
          IF ( look_Z_snow(i)%magicnr == magicnr_snow ) THEN
            ilkup_snow = i
            EXIT
          END IF
        END DO

        IF (ilkup_snow < 0) THEN

          nlkup_snow = nlkup_snow + 1
          IF (nlkup_snow > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for dry snow! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_snow(nlkup_snow)%magicnr         = magicnr_snow
          look_Z_snow(nlkup_snow)%cversion_lt     = versionstring_lookup
          look_Z_snow(nlkup_snow)%luse_tmatrix    = luse_tmatrix 
          look_Z_snow(nlkup_snow)%ldo_nonsphere   = ldo_nonsphere
          look_Z_snow(nlkup_snow)%itype_refl      = itype_refl
          look_Z_snow(nlkup_snow)%chydrotype(:)   = ' '
          look_Z_snow(nlkup_snow)%chydrotype      = 'drysnow-1mom'
          look_Z_snow(nlkup_snow)%chydroconfig(:) = ' '
          look_Z_snow(nlkup_snow)%chydroconfig    = TRIM(chydroconfig)

          ilkup_snow = nlkup_snow

          IF ( look_Z_snow(ilkup_snow)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_snow(ilkup_snow) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for dry snow (only if it has not yet been created for this wavelength)
        IF (.NOT. look_Z_snow(ilkup_snow)%is_initialized) THEN

          CALL zradar_snow_1mom_lookupcreate(&
               look_Z_snow(ilkup_snow),tableprops,&
               isnow_n0temp,&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPs,&
               ray_const%Z_fac,Dmin_s,Dmax_s,snow,DBLE(Tmeltbegin_s),&
               mixingrulestring_s_shell_dry,matrixstring_s_shell_dry,inclusionstring_s_shell_dry,&
               mixingrulestring_s_core_dry,matrixstring_s_core_dry,inclusionstring_s_core_dry,&
               impipar_lookupgen, pe_start, pe_end,    &
               linterp_mode_dualpol,                   &
               ydir_lookup_read,ydir_lookup_write, unit_lookup,hashtext)

        END IF

        ! create lookup table for melting snow:
        !======================================

        ! melt* tables also depend on rain polMP properties, hence include them
        !   in the melt* hashes in addition to the ice ones
        WRITE (chydroconfig, '(2A)') &
               TRIM(chydroconfig), TRIM(chydroconfig_pMPrain)

        ! table sizes: nTa_sw, nTm_s, nq_s in radar_data_mie.f90

        tableprops%nqi = nq_s
        tableprops%nTa = nTa_sw
        tableprops%nTm = nTm_s
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.5_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(3I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, tableprops%nTm, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_meltsnow, hashtext)

        ! Search lookup table index for melting snow for the actual parameter set:
        ilkup_meltsnow = -99
        DO i = 1, nlkup_meltsnow
          IF ( look_Z_meltsnow(i)%magicnr == magicnr_meltsnow ) THEN
            ilkup_meltsnow = i
            EXIT
          END IF
        END DO

        IF (ilkup_meltsnow < 0) THEN

          nlkup_meltsnow = nlkup_meltsnow + 1
          IF (nlkup_meltsnow > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for melting snow! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_meltsnow(nlkup_meltsnow)%magicnr          = magicnr_meltsnow
          look_Z_meltsnow(nlkup_meltsnow)%cversion_lt      = versionstring_lookup
          look_Z_meltsnow(nlkup_meltsnow)%luse_tmatrix     = luse_tmatrix
          look_Z_meltsnow(nlkup_meltsnow)%ldo_nonsphere    = ldo_nonsphere
          look_Z_meltsnow(nlkup_meltsnow)%itype_refl       = itype_refl
          look_Z_meltsnow(nlkup_meltsnow)%itype_Dref_fmelt = itype_Dref_fmelt
          look_Z_meltsnow(nlkup_meltsnow)%chydrotype(:)    = ' '
          look_Z_meltsnow(nlkup_meltsnow)%chydrotype       = 'meltsnow-1mom'
          look_Z_meltsnow(nlkup_meltsnow)%chydroconfig(:)  = ' '
          look_Z_meltsnow(nlkup_meltsnow)%chydroconfig     = TRIM(chydroconfig)

          ilkup_meltsnow = nlkup_meltsnow

          IF ( look_Z_meltsnow(ilkup_meltsnow)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_meltsnow(ilkup_meltsnow) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for melting snow (only if it has not yet been created for this wavelength)
        IF (.NOT. look_Z_meltsnow(ilkup_meltsnow)%is_initialized) THEN

          CALL zradar_meltsnow_1mom_lookupcreate(&
               look_Z_meltsnow(ilkup_meltsnow),tableprops,&
               isnow_n0temp,&
               itype_Dref_fmelt,&
               DBLE(Tmeltbegin_s),DBLE(meltdegTmin_s),&
               DBLE(Tmax_min_s),DBLE(Tmax_max_s),&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPs,pMPr,&
               ray_const%Z_fac,Dmin_s,Dmax_s,snow,&
               mixingrulestring_s_shell,matrixstring_s_shell,inclusionstring_s_shell,&
               mixingrulestring_s_core,matrixstring_s_core,inclusionstring_s_core,&
               hoststring_s_shell,hostmatrixstring_s_shell,hostinclusionstring_s_shell,&
               hoststring_s_core,hostmatrixstring_s_core,hostinclusionstring_s_core,&
               impipar_lookupgen,pe_start,pe_end,&
               linterp_mode_dualpol,&
               ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)

        END IF

      END IF ! snow

      IF (ldo_qg) THEN

        chydroconfig(:) = ' '
        IF (luse_tmatrix) THEN
          IF (toupper(TRIM(ADJUSTL(pMPg%ARmodel))) == 'POLY') THEN
            WRITE (chydroconfig, '(7es16.6,I5,A9,3es13.3,F6.1)') &
                   graupel%nu, graupel%mu, graupel%n0_const, graupel%a_geo, graupel%b_geo, &
                   Dmin_g, Dmax_g, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPg%ARmodel))), &
                   pMPg%c0, pMPg%c1, pMPg%ARmin, pMPg%sig
          ELSE
            WRITE (chydroconfig, '(7es16.6,I5,A9,F6.1)') &
                   graupel%nu, graupel%mu, graupel%n0_const, graupel%a_geo, graupel%b_geo, &
                   Dmin_g, Dmax_g, n_stuetz, &
                   toupper(TRIM(ADJUSTL(pMPg%ARmodel))), pMPg%sig
          END IF
        ELSE
          WRITE (chydroconfig, '(7es16.6,I5)') &
                 graupel%nu, graupel%mu, graupel%n0_const, graupel%a_geo, graupel%b_geo, &
                 Dmin_g, Dmax_g, n_stuetz
        END IF

        ! create lookup table for graupel:
        !======================================

        ! table sizes: nTa_gd, nq_g in radar_data_mie.f90

        tableprops%nqi = nq_g
        tableprops%nTa = nTa_gd
        tableprops%nTm = 0
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.0_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(2I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, &
              tableprops%qilow, tableprops%qiup
        IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN
          CALL hash_radar_lookup ( &
               ctableconfig, chydroconfig, cbaseconfig_g_0C, &
               magicnr_grau, hashtext)
        ELSE
          CALL hash_radar_lookup ( &
               ctableconfig, chydroconfig, cbaseconfig, &
               magicnr_grau, hashtext)
        END IF
        
        ! Search lookup table index for dry graupel for the actual parameter set:
        ilkup_grau = -99
        DO i = 1, nlkup_grau
          IF ( look_Z_graupel(i)%magicnr == magicnr_grau ) THEN
            ilkup_grau = i
            EXIT
          END IF
        END DO

        IF (ilkup_grau < 0) THEN

          nlkup_grau = nlkup_grau + 1
          IF (nlkup_grau > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for dry graupel! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_graupel(nlkup_grau)%magicnr         = magicnr_grau
          look_Z_graupel(nlkup_grau)%cversion_lt     = versionstring_lookup
          look_Z_graupel(nlkup_grau)%luse_tmatrix    = luse_tmatrix 
          look_Z_graupel(nlkup_grau)%ldo_nonsphere   = ldo_nonsphere
          look_Z_graupel(nlkup_grau)%itype_refl      = itype_refl   
          look_Z_graupel(nlkup_grau)%chydrotype(:)   = ' '
          look_Z_graupel(nlkup_grau)%chydrotype      = 'drygraupel-1mom'
          look_Z_graupel(nlkup_grau)%chydroconfig(:) = ' '
          look_Z_graupel(nlkup_grau)%chydroconfig    = TRIM(chydroconfig)
          
          ilkup_grau = nlkup_grau

          IF ( look_Z_graupel(ilkup_grau)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_graupel(ilkup_grau) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for dry graupel (only if it has not yet been created)
        IF (.NOT. look_Z_graupel(ilkup_grau)%is_initialized ) THEN

          IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN
            ! .. Tmeltbegin_g is normally the upper T-limit for the table.
            !    In this case however, the table has to span up to T0C_fwo:
            CALL zradar_graupel_1mom_lookupcreate(&
                 look_Z_graupel(ilkup_grau),tableprops,&
                 DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,&
                 ray_const%Z_fac,Dmin_g,Dmax_g,graupel,&
                 T0C_fwo,mixingrulestring_g_dry,matrixstring_g_dry,inclusionstring_g_dry,&
                 impipar_lookupgen,pe_start,pe_end,&
                 linterp_mode_dualpol,&
                 ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)
          ELSE
            CALL zradar_graupel_1mom_lookupcreate(&
                 look_Z_graupel(ilkup_grau),tableprops,&
                 DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,&
                 ray_const%Z_fac,Dmin_g,Dmax_g,graupel,&
                 Tmeltbegin_g,mixingrulestring_g_dry,matrixstring_g_dry,inclusionstring_g_dry,&
                 impipar_lookupgen,pe_start,pe_end,&
                 linterp_mode_dualpol,&
                 ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)
          END IF

        END IF

        ! create lookup table for melting graupel:
        !=========================================

        ! melt* tables also depend on rain polMP properties, hence include them
        !   in the melt* hashes in addition to the ice ones
        WRITE (chydroconfig, '(2A)') &
               TRIM(chydroconfig), TRIM(chydroconfig_pMPrain)

        ! table sizes: nTa_gw, nTm_g, nq_g in radar_data_mie.f90

        tableprops%nqi = nq_g
        tableprops%nTa = nTa_gw
        tableprops%nTm = nTm_g
        tableprops%qilow = q_crit_radar%frozen
        tableprops%qiup  = 10.0_dp**(-1.5_dp)
        ctableconfig(:) = ' '
        WRITE (ctableconfig, '(3I5,2es16.6)') &
              tableprops%nqi, tableprops%nTa, tableprops%nTm, &
              tableprops%qilow, tableprops%qiup
        CALL hash_radar_lookup ( &
             ctableconfig, chydroconfig, cbaseconfig, &
             magicnr_meltgrau, hashtext)

        ! Search lookup table index for melting graupel for the actual parameter set:
        ilkup_meltgrau = -99
        DO i = 1, nlkup_meltgrau
          IF ( look_Z_meltgraupel(i)%magicnr == magicnr_meltgrau ) THEN
            ilkup_meltgrau = i
            EXIT
          END IF
        END DO

        IF (ilkup_meltgrau < 0) THEN

          nlkup_meltgrau = nlkup_meltgrau + 1
          IF (nlkup_meltgrau > nmax_lookup) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Too many different sets of Mie lookup tables for melting graupel! '// &
                    'Increase parameter nmax_lookup! Stop!'
            STOP
          END IF
          look_Z_meltgraupel(nlkup_meltgrau)%magicnr          = magicnr_meltgrau
          look_Z_meltgraupel(nlkup_meltgrau)%cversion_lt      = versionstring_lookup
          look_Z_meltgraupel(nlkup_meltgrau)%luse_tmatrix     = luse_tmatrix    
          look_Z_meltgraupel(nlkup_meltgrau)%ldo_nonsphere    = ldo_nonsphere   
          look_Z_meltgraupel(nlkup_meltgrau)%itype_refl       = itype_refl      
          look_Z_meltgraupel(nlkup_meltgrau)%itype_Dref_fmelt = itype_Dref_fmelt
          look_Z_meltgraupel(nlkup_meltgrau)%chydrotype(:)    = ' '
          look_Z_meltgraupel(nlkup_meltgrau)%chydrotype       = 'meltgraupel-1mom'
          look_Z_meltgraupel(nlkup_meltgrau)%chydroconfig(:)  = ' '
          look_Z_meltgraupel(nlkup_meltgrau)%chydroconfig     = TRIM(chydroconfig)
          
          ilkup_meltgrau = nlkup_meltgrau

          IF ( look_Z_meltgraupel(ilkup_meltgrau)%is_initialized ) THEN
            WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                    ': Strange error: look_Z_meltgrauel(ilkup_meltgrau) is already '// &
                    'initialized when it should not! Stop!'
            STOP
          END IF

        END IF

        ! create lookup table for melting graupel (only if it has not yet been created)
        IF (.NOT. look_Z_meltgraupel(ilkup_meltgrau)%is_initialized) THEN

          CALL zradar_meltgraupel_1mom_lookupcreate( &
               look_Z_meltgraupel(ilkup_meltgrau),tableprops,&
               igraupel_type,&
               itype_Dref_fmelt,&
               DBLE(Tmeltbegin_g),DBLE(meltdegTmin_g),&
               DBLE(Tmax_min_g),DBLE(Tmax_max_g),&
               DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,pMPr,&
               ray_const%Z_fac,Dmin_g,Dmax_g,graupel,&
               mixingrulestring_g,matrixstring_g,inclusionstring_g,&
               mixingrulestring_g_shell,matrixstring_g_shell,inclusionstring_g_shell,&
               hoststring_g,hostmatrixstring_g,hostinclusionstring_g,&
               impipar_lookupgen,pe_start,pe_end,&
               linterp_mode_dualpol,&
               ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)

        END IF

        IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN

          ! In case of dynamic wet growth: create lookup table with no wet growth to be used together with the "normal" table:
          
          CALL hash_radar_lookup ( &
               ctableconfig, chydroconfig, cbaseconfig_g_0C, &
               magicnr_meltgrau_0C, hashtext)

          ! Search lookup table index for melting graupel for the actual parameter set:
          ilkup_meltgrau_0C = -99
          DO i = 1, nlkup_meltgrau_0C
            IF ( look_Z_meltgraupel_0C(i)%magicnr == magicnr_meltgrau_0C ) THEN
              ilkup_meltgrau_0C = i
              EXIT
            END IF
          END DO

          IF (ilkup_meltgrau_0C < 0) THEN

            nlkup_meltgrau_0C = nlkup_meltgrau_0C + 1
            IF (nlkup_meltgrau_0C > nmax_lookup) THEN
              WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                   ': Too many different sets of Mie lookup tables for melting graupel without wetgrowth! '// &
                   'Increase parameter nmax_lookup! Stop!'
              STOP
            END IF
            ! LUT for melting graupel without wet growth has same metadata as "normal" table
            !  except the magicnr_meltgrau_0C:
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%magicnr          = magicnr_meltgrau_0C
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%cversion_lt      = versionstring_lookup
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%luse_tmatrix     = luse_tmatrix    
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%ldo_nonsphere    = ldo_nonsphere   
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%itype_refl       = itype_refl      
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%itype_Dref_fmelt = itype_Dref_fmelt
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%chydrotype(:)    = ' '
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%chydrotype       = 'meltgraupel-1mom'
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%chydroconfig(:)  = ' '
            look_Z_meltgraupel_0C(nlkup_meltgrau_0C)%chydroconfig     = TRIM(chydroconfig)
          
            ilkup_meltgrau_0C = nlkup_meltgrau_0C

            IF ( look_Z_meltgraupel_0C(ilkup_meltgrau_0C)%is_initialized ) THEN
              WRITE (*,*) 'ERROR in '//TRIM(yzroutine)//&
                   ': Strange error: look_Z_meltgrauel_0C(ilkup_meltgrau_0C) is already '// &
                   'initialized when it should not! Stop!'
              STOP
            END IF

          END IF

          IF (.NOT. look_Z_meltgraupel_0C(ilkup_meltgrau_0C)%is_initialized) THEN

            CALL zradar_meltgraupel_1mom_lookupcreate( &
                 look_Z_meltgraupel_0C(ilkup_meltgrau_0C),tableprops,&
                 igraupel_type,&
                 itype_Dref_fmelt,&
                 T0C_fwo, 0.0_dp,&
                 DBLE(Tmax_min_g),DBLE(Tmax_max_g),&
                 DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,pMPr,&
                 ray_const%Z_fac,Dmin_g,Dmax_g,graupel,&
                 mixingrulestring_g,matrixstring_g,inclusionstring_g,&
                 mixingrulestring_g_shell,matrixstring_g_shell,inclusionstring_g_shell,&
                 hoststring_g,hostmatrixstring_g,hostinclusionstring_g,&
                 impipar_lookupgen,pe_start,pe_end,&
                 linterp_mode_dualpol,&
                 ydir_lookup_read,ydir_lookup_write,unit_lookup,hashtext)

          END IF

        END IF

      END IF ! graupel

      ! Use pointers to single components of the above vectors of lookup table types
      ! These are necessary as input parameters to the lookup table interpolation
      ! functions to enable correct inlining on the NEC:
      IF (ldo_qr) THEN
        look_Z_rain_p        => look_Z_rain(ilkup_rain)
      END IF
      IF (ldo_qi .AND. .NOT.ldo_qi_rayleigh) THEN
        look_Z_ice_p         => look_Z_ice(ilkup_ice)
        look_Z_meltice_p     => look_Z_meltice(ilkup_meltice)
      END IF
      IF (ldo_qs) THEN
        look_Z_snow_p        => look_Z_snow(ilkup_snow)
        look_Z_meltsnow_p    => look_Z_meltsnow(ilkup_meltsnow)
      END IF
      IF (ldo_qg) THEN
        look_Z_graupel_p     => look_Z_graupel(ilkup_grau)
        look_Z_meltgraupel_p => look_Z_meltgraupel(ilkup_meltgrau)
        IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN
          look_Z_meltgraupel_0C_p => look_Z_meltgraupel_0C(ilkup_meltgrau_0C)
        END IF
      END IF


      ! Reflectivity of rain after Mie (Integration with respect to diameter):
      IF (ldo_qr) THEN

!$omp parallel private (i,j,k,T_a_vec,hlp,q_x_vec,T_m_vec,flag_interp)
!$omp do
        DO k = ku, ko
          DO j = jlow, jup

            flag_interp(:) = .FALSE.

            DO i = ilow, iup

              T_a_vec(i) = t(i,j,k)
              hlp        = rho(i,j,k)
              T_m_vec(i) = T_a_vec(i)   ! dummy

              q_x_vec(i) = qr(i,j,k) * hlp
 
              IF (q_x_vec(i) >= q_crit_radar%liquid) THEN
                flag_interp(i) = .TRUE.
              END IF

            END DO

            IF (ANY(flag_interp(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_rain_p,     &
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        & ! for 1-mom rain just a dummy
                   flag_interp    = flag_interp(:),    &
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_pure)
            END IF
            
          END DO
        END DO
!$omp end do
!$omp end parallel

      END IF

      ! Reflectivity of cloud ice after Mie (Integration with respect to diameter):
      IF (ldo_qi .AND. .NOT. ldo_qi_rayleigh) THEN

!$omp parallel private (i,j,k,T_a_vec,hlp,q_x_vec,T_m_vec,flag_interp,is_melting)
!$omp do
        DO k = ku, ko
          DO j = jlow, jup

            flag_interp(:) = .FALSE.
            is_melting(:) = .FALSE.

            DO i = ilow, iup

              T_a_vec(i) = t(i,j,k)
              hlp        = rho(i,j,k)

              q_x_vec(i) = qi(i,j,k) * hlp

              IF (q_x_vec(i) >= q_crit_radar%frozen) THEN
                ! Table lookup:
                IF (T_a_vec(i) > Tmeltbegin_i) THEN
                  T_m_vec(i) = MAX(Tmax_i(IND_IJ_2D), Tmax_min_i)
                  is_melting(i) = .TRUE.
                ELSE
                  T_m_vec(i) = T_a_vec(i)     ! dummy
                  is_melting(i) = .FALSE.
                END IF
                flag_interp(i) = .TRUE.
              END IF

            END DO

            IF (ANY(flag_interp(:) .AND. is_melting(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_meltice_p,&
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        &
                   flag_interp    = (flag_interp(:).AND.is_melting(:)),&
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_melt  )
            END IF

            IF (ANY(flag_interp(:) .AND. .NOT. is_melting(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_ice_p,  &
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        & ! just a dummy for dry cloud ice
                   flag_interp    = (flag_interp(:).AND..NOT.is_melting(:)),&
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_pure  )
            END IF

          END DO
        END DO
!$omp end do
!$omp end parallel

      END IF

      ! Reflectivity of snow after Mie (Integration with respect to diameter):
      IF (ldo_qs) THEN

!$omp parallel private (i,j,k,T_a_vec,hlp,q_x_vec,T_m_vec,flag_interp,is_melting)
!$omp do
        DO k = ku, ko
          DO j = jlow, jup

            flag_interp(:) = .FALSE.
            is_melting(:) = .FALSE.

            DO i = ilow, iup

              T_a_vec(i) = t(i,j,k)
              hlp        = rho(i,j,k)

              q_x_vec(i) = qs(i,j,k) * hlp

              IF (q_x_vec(i) >= q_crit_radar%frozen) THEN
                ! Table lookup:
                IF (T_a_vec(i) > Tmeltbegin_s) THEN
                  T_m_vec(i) = MAX(Tmax_s(IND_IJ_2D), Tmax_min_s)
                  is_melting(i) = .TRUE.
                ELSE
                  T_m_vec(i) = T_a_vec(i)     ! dummy
                  is_melting(i) = .FALSE.
                END IF
                flag_interp(i) = .TRUE.
              END IF

            END DO

            IF (ANY(flag_interp(:) .AND. is_melting(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_meltsnow_p, &
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        &
                   flag_interp    = (flag_interp(:).AND.is_melting(:)),&
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_melt  )
            END IF

            IF (ANY(flag_interp(:) .AND. .NOT. is_melting(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_snow_p,     &
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        & ! just a dummy for dry snow
                   flag_interp    = (flag_interp(:).AND..NOT.is_melting(:)),&
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_pure  )
            END IF

          END DO
        END DO
!$omp end do
!$omp end parallel

      END IF

      ! Reflectivity of graupel after Mie (Integration with respect to diameter):
      IF (ldo_qg) THEN

!$omp parallel private (i,j,k,T_a_vec,T_a_wg_vec,T0C_vec,hlp,q_x_vec,T_m_vec,&
!$omp&                  flag_interp,is_melting,a_w_vec,b_w_vec,T_gt_0C)
!$omp do
        DO k = ku, ko
          DO j = jlow, jup

            flag_interp(:) = .FALSE.
            is_melting(:) = .FALSE.

            DO i = ilow, iup

              T_a_vec(i) = t(i,j,k)
              hlp        = rho(i,j,k)

              q_x_vec(i) = qg(i,j,k) * hlp

              IF (q_x_vec(i) >= q_crit_radar%frozen) THEN
                ! Prepare table lookup:
                IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN
                  T0C_vec(:) = T0C_fwo
                  IF (T_a_vec(i) > Tmin_g(IND_IJ_2D)) THEN
                    T_m_vec(i) = MAX(Tmax_g(IND_IJ_2D), Tmax_min_g)
                    T_a_wg_vec(i) = T_a_vec(i) - Tmin_g(IND_IJ_2D) + Tmeltbegin_g ! shifted T_a
                    a_w_vec(i)    = MAX(Tmin_g(IND_IJ_2D)-Tmeltbegin_g, 0.0_dp) / &
                                    MAX(T0C_fwo-Tmeltbegin_g          , 1e-6_dp)
                    b_w_vec(i)    = MAX(T0C_fwo-T_a_vec(i)       , 0.0_dp) / &
                                    MAX(T0C_fwo-Tmin_g(IND_IJ_2D), 1e-6_dp)
                    is_melting(i) = .TRUE.
                    T_gt_0C(i)    = T_a_vec(i) >= T0C_fwo
                  ELSE
                    T_m_vec(i)    = T_a_vec(i)     ! dummy
                    T_a_wg_vec(i) = T_a_vec(i)     ! dummy
                    a_w_vec(i)    = 1.0_dp         ! dummy
                    b_w_vec(i)    = 1.0_dp         ! dummy
                    is_melting(i) = .FALSE.
                    T_gt_0C(i)    = .FALSE.        ! dummy
                  END IF
                ELSE
                  IF (T_a_vec(i) > Tmeltbegin_g) THEN
                    T_m_vec(i) = MAX(Tmax_g(IND_IJ_2D), Tmax_min_g)
                    is_melting(i) = .TRUE.
                  ELSE
                    T_m_vec(i)    = T_a_vec(i)     ! dummy
                    is_melting(i) = .FALSE.
                  END IF
                END IF
                flag_interp(i) = .TRUE.
              END IF

            END DO

            IF (ANY(flag_interp(:) .AND. is_melting(:))) THEN

              IF (ldynamic_wetgrowth_gh .AND. Tmeltbegin_g < T0C_fwo) THEN

                ! Branch for Tmin_g < T_a < T0C:
                !
                !    z = b*z1 + a*(1-b)*z2 + (1-a)*(1-b)*z3
                !
                ! weights a, b etc. are input via n_i argument:
                                
                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T_a_wg_vec(:),     &  ! shifted T_a
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:).AND..NOT.T_gt_0C(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     n_i            = b_w_vec(:),        &
                     ext_tune_fac   = ext_tune_fac_melt  )

                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_0C_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T0C_vec(:),        & ! T = T0C_fwo
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:).AND..NOT.T_gt_0C(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     n_i            = a_w_vec(:)*(1.0_dp-b_w_vec(:)),  &
                     ext_tune_fac   = ext_tune_fac_melt  )

                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T_a_vec(:),        & ! T = T_a
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:).AND..NOT.T_gt_0C(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     n_i            = (1.0_dp-a_w_vec(:))*(1.0_dp-b_w_vec(:)),  &
                     ext_tune_fac   = ext_tune_fac_melt  )

                ! Branch for T_a >= T0C:
                !
                !    z = a*z1 + (1-a)*z2
                
                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_0C_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T_a_vec(:),        &
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:).AND.T_gt_0C(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     n_i            = a_w_vec(:)      ,  &
                     ext_tune_fac   = ext_tune_fac_melt  )

                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T_a_vec(:),        &
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:).AND.T_gt_0C(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     n_i            = (1.0_dp-a_w_vec(:)),  &
                     ext_tune_fac   = ext_tune_fac_melt  )
                
              ELSE

                CALL zradar_triinterp_lookup_add_vec(    &
                     look_Z_type    = look_Z_meltgraupel_p,&
                     q_i            = q_x_vec(:),        &
                     T_a            = T_a_vec(:),        &
                     T_m            = T_m_vec(:),        &
                     flag_interp    = (flag_interp(:).AND.is_melting(:)),&
                     ilow           = ilow,              &
                     iup            = iup,               &
                     zh_radar_int   = zh_local(:,j,k),   &
                     ah_radar_int   = ah_local(:,j,k),   &
                     zv_radar_int   = zv_local(:,j,k),   &
                     rrhv_radar_int = rrhv_local(:,j,k), &
                     irhv_radar_int = irhv_local(:,j,k), &
                     kdp_radar_int  = kdp_local(:,j,k),  &
                     adp_radar_int  = adp_local(:,j,k),  &
                     zvh_radar_int  = zvh_local(:,j,k),  &
                     ext_tune_fac   = ext_tune_fac_melt  )

              END IF
              
            END IF

            IF (ANY(flag_interp(:) .AND. .NOT. is_melting(:))) THEN
              CALL zradar_triinterp_lookup_add_vec(    &
                   look_Z_type    = look_Z_graupel_p,  &
                   q_i            = q_x_vec(:),        &
                   T_a            = T_a_vec(:),        &
                   T_m            = T_m_vec(:),        & ! just a dummy for dry graupel
                   flag_interp    = (flag_interp(:).AND..NOT.is_melting(:)),&
                   ilow           = ilow,              &
                   iup            = iup,               &
                   zh_radar_int   = zh_local(:,j,k),   &
                   ah_radar_int   = ah_local(:,j,k),   &
                   zv_radar_int   = zv_local(:,j,k),   &
                   rrhv_radar_int = rrhv_local(:,j,k), &
                   irhv_radar_int = irhv_local(:,j,k), &
                   kdp_radar_int  = kdp_local(:,j,k),  &
                   adp_radar_int  = adp_local(:,j,k),  &
                   zvh_radar_int  = zvh_local(:,j,k),  &
                   ext_tune_fac   = ext_tune_fac_pure  )
            END IF

          END DO
        END DO
!$omp end do
!$omp end parallel

      END IF

    ELSE  ! .NOT. llookup

!$omp parallel private (i,j,k,T_a,hlp,q_r,q_s,q_g,&
!$omp&                  zh,ah,zv,rrhv,irhv,kdp,adp,zvh,T_m,mgd,&
!$omp&                  Tmeltbegin_g_loc,meltdegTmin_g_loc)
!$omp do
      DO k = ku, ko
        DO j = jlow, jup
          DO i = ilow, iup

            T_a = t(i,j,k)
            hlp = rho(i,j,k)

            IF (luse_tmatrix .AND. MODULO(i,MIN(10,iup)) == 0 .AND. MODULO(j,MIN(10,jup)) == 0) THEN
              ! computation takes long, so print some progress information to stdout:
              WRITE (*,'("  Tmatrix-computations on proc ",i4," for grid point (i,j,k) ",3i6,' // &
                   '"   from (ni,nj,nk) ",3i6," ...")') &
                   myproc, i,j,k, iup,jup,ko
            END IF

            ! Reflectivity of rain after Mie (Integration with respect to diameter):
            IF (ldo_qr) THEN

              q_r = qr(i,j,k) * hlp

              IF (q_r >= q_crit_radar%liquid) THEN
                ! NOTE: mgd calc functions invalid for q==0.0
                mgd = mgd_1mom(rain,q_r,rain%n0_const)
                CALL zradar_rain_mie_vec(mgd,m_w(i,j,k),&
                     DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPr,&
                     ray_const%Z_fac,Dmin_r,Dmax_r,&
                     zh,ah,zv,rrhv,irhv,kdp,adp,zvh)
                zh_local(i,j,k) = zh_local(i,j,k) + zh
                ah_local(i,j,k) = ah_local(i,j,k) + ah*ext_tune_fac_pure
                zv_local(i,j,k) = zv_local(i,j,k) + zv
                rrhv_local(i,j,k) = rrhv_local(i,j,k) + rrhv
                irhv_local(i,j,k) = irhv_local(i,j,k) + irhv
                kdp_local(i,j,k) = kdp_local(i,j,k) + kdp
                adp_local(i,j,k) = adp_local(i,j,k) + adp*ext_tune_fac_pure
                zvh_local(i,j,k) = zvh_local(i,j,k) + zvh
              END IF  ! q_r(i) >= q_crit_radar

            END IF

            ! Reflectivity of cloud ice after Mie (Integration with respect to diameter):
            IF (ldo_qi .AND. .NOT. ldo_qi_rayleigh) THEN

              q_i = qi(i,j,k) * hlp

              IF (q_i >= q_crit_radar%frozen) THEN

                n_i = nice_mono_1mom(T_a)

                ! Mimic a monodisperse size distribution around a mean mass of q_i/n_i by using
                ! 2-mom mgd method and choosing very large values for ice%mu and ice%nu (in
                ! init_1mom_types(), radar_interface.f90):
                ! NOTE: mgd calc functions invalid for q==0.0
                mgd = mgd_2mom(ice,q_i,n_i)

                IF (T_a > Tmeltbegin_i) THEN
                  CALL zradar_wetice_mie_vec(mgd,T_a,m_i(i,j,k),m_w(i,j,k),&
                       itype_Dref_fmelt,&
                       DBLE(Tmeltbegin_i),DBLE(meltdegTmin_i),&
                       DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPi,pMPr,&
                       ray_const%Z_fac,ice,rain,Dmin_i,Dmax_i,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring_i, matrixstring_i, inclusionstring_i, &
                       hoststring_i, hostmatrixstring_i, hostinclusionstring_i, &
                       REAL(Tmax_i(IND_IJ_2D),kind=dp))
                  ah  = ah  * ext_tune_fac_melt
                  adp = adp * ext_tune_fac_melt
                ELSE
                  CALL zradar_ice_mie_vec(mgd,m_i(i,j,k),&
                       DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPi,&
                       ray_const%Z_fac,ice,Dmin_i,Dmax_i,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring_i_dry, matrixstring_i_dry, inclusionstring_i_dry)
                  ah  = ah  * ext_tune_fac_pure
                  adp = adp * ext_tune_fac_pure
                END IF  ! T_a(i) > Tmeltbegin_i
                zh_local(i,j,k) = zh_local(i,j,k) + zh
                ah_local(i,j,k) = ah_local(i,j,k) + ah
                zv_local(i,j,k) = zv_local(i,j,k) + zv
                rrhv_local(i,j,k) = rrhv_local(i,j,k) + rrhv
                irhv_local(i,j,k) = irhv_local(i,j,k) + irhv
                kdp_local(i,j,k) = kdp_local(i,j,k) + kdp
                adp_local(i,j,k) = adp_local(i,j,k) + adp
                zvh_local(i,j,k) = zvh_local(i,j,k) + zvh
              END IF  ! q_i >= q_crit_radar

            END IF  ! ldo_qi .and. .not.ldo_qi_rayleigh

            ! Reflectivity of snow after Mie (Integration with respect to diameter):
            IF (ldo_qs) THEN

              q_s = qs(i,j,k) * hlp

              IF (q_s >= q_crit_radar%frozen) THEN
                ! NOTE: mgd calc functions invalid for q==0.0
                mgd = mgd_1mom(snow,q_s,n0_s(i,j,k))
                IF (T_a > Tmeltbegin_s) THEN
                  CALL zradar_wetsnow_mie_vec(mgd,T_a,m_i(i,j,k),m_w(i,j,k),&
                       0.5d0,itype_Dref_fmelt,&
                       DBLE(Tmeltbegin_s),DBLE(meltdegTmin_s),&
                       DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPs,pMPr,&
                       ray_const%Z_fac,snow,rain,Dmin_s,Dmax_s,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring_s_shell, matrixstring_s_shell, inclusionstring_s_shell, &
                       mixingrulestring_s_core, matrixstring_s_core, inclusionstring_s_core, &
                       hoststring_s_shell, hostmatrixstring_s_shell, hostinclusionstring_s_shell, &
                       hoststring_s_core, hostmatrixstring_s_core, hostinclusionstring_s_core, &
                       REAL(Tmax_s(IND_IJ_2D),kind=dp))
                  ah  = ah  * ext_tune_fac_melt
                  adp = adp * ext_tune_fac_melt
                ELSE
                  CALL zradar_snow_mie_vec(mgd,m_i(i,j,k),&
                       0.5d0,DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPs,&
                       ray_const%Z_fac,snow,Dmin_s,Dmax_s,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring_s_shell_dry, matrixstring_s_shell_dry, inclusionstring_s_shell_dry, &
                       mixingrulestring_s_core_dry, matrixstring_s_core_dry, inclusionstring_s_core_dry)
                  ah  = ah  * ext_tune_fac_pure
                  adp = adp * ext_tune_fac_pure
                END IF  ! T_a > Tmeltbegin_s
                zh_local(i,j,k) = zh_local(i,j,k) + zh
                ah_local(i,j,k) = ah_local(i,j,k) + ah
                zv_local(i,j,k) = zv_local(i,j,k) + zv
                rrhv_local(i,j,k) = rrhv_local(i,j,k) + rrhv
                irhv_local(i,j,k) = irhv_local(i,j,k) + irhv
                kdp_local(i,j,k) = kdp_local(i,j,k) + kdp
                adp_local(i,j,k) = adp_local(i,j,k) + adp
                zvh_local(i,j,k) = zvh_local(i,j,k) + zvh
              END IF  ! q_s >= q_crit_radar

            END IF   ! ldo_qs

            ! Reflectivity of graupel after Mie (Integration with respect to diameter):
            IF (ldo_qg) THEN

              q_g = qg(i,j,k) * hlp

              IF (q_g >= q_crit_radar%frozen) THEN

                ! NOTE: mgd calc functions invalid for q==0.0
                mgd = mgd_1mom(graupel,q_g,graupel%n0_const)

                Tmeltbegin_g_loc = REAL(Tmin_g(IND_IJ_2D),kind=dp)
                IF (T_a > Tmeltbegin_g_loc) THEN
                  ! In case of dynamic wetgrowth, scale down meltdegTmin_g so that
                  ! it reaches 0 IF Tmin_g approaches Tmin_f = T0C_fwo:
                  meltdegTmin_g_loc = meltdegTmin_g - &
                       meltdegTmin_g / MAX(Tmin_f-Tmeltbegin_g,1e-6_dp) * (Tmeltbegin_g_loc-Tmeltbegin_g)
                  SELECT CASE (igraupel_type)
                  CASE(1)
                    CALL zradar_wetgr_mie_vec(mgd,T_a,m_i(i,j,k),m_w(i,j,k),&
                         itype_Dref_fmelt,&
                         DBLE(Tmeltbegin_g_loc),DBLE(meltdegTmin_g_loc),&
                         DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,pMPr,&
                         ray_const%Z_fac,graupel,rain,Dmin_g,Dmax_g,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                         mixingrulestring_g, matrixstring_g, inclusionstring_g, &
                         hoststring_g, hostmatrixstring_g, hostinclusionstring_g, &
                         REAL(Tmax_g(IND_IJ_2D),kind=dp))
                  CASE(2)
                    CALL zradar_wetgr_twosph_mie_vec(mgd,T_a,m_i(i,j,k),m_w(i,j,k),&
                         itype_Dref_fmelt,&
                         DBLE(Tmeltbegin_g_loc),DBLE(meltdegTmin_g_loc),&
                         DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,pMPr,&
                         ray_const%Z_fac,graupel,rain,Dmin_g,Dmax_g,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                         mixingrulestring_g, matrixstring_g, inclusionstring_g,&
                         mixingrulestring_g_shell, matrixstring_g_shell, inclusionstring_g_shell,&
                         REAL(Tmax_g(IND_IJ_2D),kind=dp))
                  CASE(3)
                    CALL zradar_wetgr_wsph_mie_vec(mgd,T_a,m_i(i,j,k),m_w(i,j,k),&
                         itype_Dref_fmelt,&
                         DBLE(Tmeltbegin_g_loc),DBLE(meltdegTmin_g_loc),&
                         DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,pMPr,&
                         ray_const%Z_fac,graupel,rain,Dmin_g,Dmax_g,&
                         zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                         mixingrulestring_g, matrixstring_g, inclusionstring_g,&
                         REAL(Tmax_g(IND_IJ_2D),kind=dp))
                  END SELECT
                  ah  = ah  * ext_tune_fac_melt
                  adp = adp * ext_tune_fac_melt
                ELSE
                  CALL zradar_graupel_mie_vec(mgd,m_i(i,j,k),&
                       DBLE(lambda_radar_3digits),luse_tmatrix,ldo_nonsphere,pMPg,&
                       ray_const%Z_fac,graupel,Dmin_g,Dmax_g,&
                       zh,ah,zv,rrhv,irhv,kdp,adp,zvh,&
                       mixingrulestring_g_dry, matrixstring_g_dry, inclusionstring_g_dry)
                  ah  = ah  * ext_tune_fac_pure
                  adp = adp * ext_tune_fac_pure
                END IF  ! T_a(i) > Tmeltbegin_g_loc
                zh_local(i,j,k) = zh_local(i,j,k) + zh
                ah_local(i,j,k) = ah_local(i,j,k) + ah
                zv_local(i,j,k) = zv_local(i,j,k) + zv
                rrhv_local(i,j,k) = rrhv_local(i,j,k) + rrhv
                irhv_local(i,j,k) = irhv_local(i,j,k) + irhv
                kdp_local(i,j,k) = kdp_local(i,j,k) + kdp
                adp_local(i,j,k) = adp_local(i,j,k) + adp
                zvh_local(i,j,k) = zvh_local(i,j,k) + zvh
              END IF  ! q_g >= q_crit_radar

            END IF  ! ldo_qg

          END DO  ! i
        END DO  ! j
      END DO  ! k
!$omp end do
!$omp end parallel

    END IF  ! .NOT. llookup

    IF ( PRESENT(zh_radar) ) THEN
      zh_radar = zh_local
    ENDIF
    IF ( PRESENT(ah_radar) ) THEN
      ah_radar = ah_local
    ENDIF

    IF ( PRESENT(zv_radar) ) THEN
      zv_radar = zv_local
    ENDIF
    !rhv_radar = SQRT((rrhv_local*rrhv_local+irhv_local*irhv_local)/(zh_local*zv_local))
    IF ( PRESENT(rrhv_radar) ) THEN
      rrhv_radar = rrhv_local
    ENDIF
    IF ( PRESENT(irhv_radar) ) THEN
      irhv_radar = irhv_local
    ENDIF
    IF ( PRESENT(kdp_radar) ) THEN
      kdp_radar = kdp_local
    ENDIF
    IF ( PRESENT(adp_radar) ) THEN
      adp_radar = adp_local
    ENDIF
    IF ( PRESENT(zvh_radar) ) THEN
      zvh_radar = zvh_local
    ENDIF

    DEALLOCATE (zh_local, ah_local, &
                zv_local, rrhv_local, irhv_local, kdp_local, adp_local, zvh_local, &
                n0_s, m_i, m_w)

    RETURN
  END SUBROUTINE radar_mie_1mom_vec


  ! Rayleigh-approximation, 1-mom-scheme, Oguchi-formulation of refractivity index. Very fast, but not very precisely.
  ! Radar wave lenght only applied in calculation of K_W^2 and K_W_0^2.
  SUBROUTINE radar_rayleigh_oguchi_1mom_vec(myproc, lambda_radar, &
       itype_gscp_fwo, isnow_n0temp, &
       Tmeltbegin_i, meltdegTmin_i, &
       Tmeltbegin_s, meltdegTmin_s, &
       Tmeltbegin_g, meltdegTmin_g, &
       rho, t, qc, qr, qi, qs, qg, &
       Tmax_i, Tmax_s, Tmax_g, &
       Tmin_g, &
       ilow, iup, jlow, jup, klow, kup, &
       lalloc_qi, lalloc_qs, lalloc_qg, &
       zh_radar, lhydrom_choice_testing )

    INTEGER,                 INTENT(IN)         :: itype_gscp_fwo, isnow_n0temp

    INTEGER,                 INTENT(IN)         :: myproc, ilow, iup, jlow, jup, klow, kup

    REAL(KIND=dp),           INTENT(IN)         :: lambda_radar, &
                                                   Tmeltbegin_i, meltdegTmin_i, &
                                                   Tmeltbegin_s, meltdegTmin_s, &
                                                   Tmeltbegin_g, meltdegTmin_g

    ! Model fields in model's wp, not dp:
    REAL(KIND=wp), DIMENSION(:,:,:), INTENT(IN) :: rho, t, qc, qr, qi, qs, qg
    ! Tmax for degree of melting in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmax_i, Tmax_s, Tmax_g
    ! Tmin for dynamic wet growth of graupel in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmin_g

    LOGICAL, INTENT(in)                         ::  lalloc_qi, lalloc_qs, lalloc_qg

    REAL(KIND=dp), INTENT(OUT), OPTIONAL        :: zh_radar(:,:,:)

    LOGICAL, INTENT(in), OPTIONAL               :: lhydrom_choice_testing(6)

    ! Local Variables
    !----------------

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) :: zh_local, n0_s

    INTEGER :: i, j, k, ku, ko, ni, nj, nk
    REAL(KIND=dp) :: T_a, fmelt_mean, Dtmp, &
                     q_r, &
                     q_c, &
                     q_g, &
                     q_s, &
                     q_i, n_i, x_i, &
                     hlp, &
                     Tmeltbegin_g_loc, meltdegTmin_g_loc  ! for ldynamic_wetgrowth_gh
    REAL(KIND=dp) :: K_w2, K_i2, K2divRho2, Dref
    REAL(KIND=dp), PARAMETER :: Dexpo_i = 0.5d0, Dref_i  = 100d-6, &
                                Dexpo_s = 0.5d0, &
                                Dexpo_g = 0.6d0

    COMPLEX(kind=dp) :: Ki_a, Kw_a
    COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: m_i, m_w

    LOGICAL, PARAMETER :: mur_is_fixed = .TRUE.
    !LOGICAL :: mur_is_fixed, n0r_is_fixed
    LOGICAL :: ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg

    TYPE(t_mgd_params) :: mgd_rain

    CHARACTER (len=*), PARAMETER :: yzroutine = 'emvorado::radar_rayleigh_oguchi_1mom_vec()'

    REAL(KIND=dp)              :: lambda_radar_3digits
    CHARACTER(len=12)          :: clambda_radar

    ! .. Round lambda_radar to 3 significant digits, to avoid differences to slightly differing numeric
    !     lambda_radar values from namelists and DATA files:
    WRITE(clambda_radar, '(es12.2)') lambda_radar
    READ (clambda_radar, *)          lambda_radar_3digits
    
    ! .. Get field dimensions from input rho (assume that all hydrometeors and t are the same size):
    ni = SIZE(rho, DIM=1)
    nj = SIZE(rho, DIM=2)
    nk = SIZE(rho, DIM=3)

    ! .. Check field dimensions of output fields:
    IF (PRESENT(zh_radar)) THEN
      IF (ni /= SIZE(zh_radar, dim=1) .OR. nj /= SIZE(zh_radar, dim=2)) THEN
        WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                    '(): Inconsistent field dimensions of qc and zh_radar!'
        STOP
      END IF
    END IF

    ! .. Allocate local work arrays:
    ALLOCATE(zh_local(ni,nj,nk))
    ALLOCATE(n0_s(ni,nj,nk))
    ALLOCATE(m_w(ni,nj,nk), m_i(ni,nj,nk))


    zh_local = 0d0
    IF (PRESENT(zh_radar)) THEN
      ku = klow
      ko = kup
    ELSE
      ! Simply do nothing
      ku = 2
      ko = 1
    END IF

    IF (PRESENT(lhydrom_choice_testing)) THEN
      ldo_qc = lhydrom_choice_testing(1)
      ldo_qr = lhydrom_choice_testing(2)
      ldo_qi = lhydrom_choice_testing(3) .AND. lalloc_qi
      ldo_qs = lhydrom_choice_testing(4) .AND. lalloc_qs
      ldo_qg = lhydrom_choice_testing(5) .AND. lalloc_qg
    ELSE
      ldo_qc = .TRUE.
      ldo_qr = .TRUE.
      ldo_qi = lalloc_qi
      ldo_qs = lalloc_qs
      ldo_qg = lalloc_qg
    END IF

    ! On its first call or if otherwise necessary, this computes
    ! some constant prefactors needed below and stores it on global struct "ray_const":
    CALL init_radar_rayleigh_consts( &
         .TRUE., .TRUE., lalloc_qi, lalloc_qs, lalloc_qg, .FALSE., &
         DBLE(lambda_radar_3digits))

    ! Refrative index of water after Ray (1972), only valid for -10 < T < 30 degree C
    m_w = RESHAPE(m_complex_water_ray_vec(DBLE(lambda_radar_3digits),&
         MIN(MAX(RESHAPE(t(:,:,:),(/ni*nj*nk/))-T0C_fwo,mw_Tmin),mw_Tmax),ni*nj*nk), &
         (/ni,nj,nk/) )
    ! Refrative index of ice after C. Maetzler, only valid for -80 < T < 0 degree C
    m_i = RESHAPE(m_complex_ice_maetzler_vec(DBLE(lambda_radar_3digits),&
         MIN(MAX(RESHAPE(t(:,:,:),(/ni*nj*nk/))-T0C_fwo,mi_Tmin),mi_Tmax), ni*nj*nk), &
         (/ni,nj,nk/) )

    ! Precalculate NWP-field dependent snow number density field (but only if we
    ! we consider snow at all):
    IF (ldo_qs) THEN
      CALL snow_1mom_n0(itype_gscp_fwo, isnow_n0temp, &
                        t, rho, qs, snow%a_geo, &
                        n0_s)
    END IF

    DO k= ku, ko
      DO j = jlow, jup
        DO i = ilow, iup

          T_a = t(i,j,k)
          hlp = rho(i,j,k)

          Kw_a = (m_w(i,j,k)**2-1d0)/(m_w(i,j,k)**2+2d0)
          K_w2 = ABS(Kw_a)**2
          Ki_a = (m_i(i,j,k)**2-1d0)/(m_i(i,j,k)**2+2d0)
          K_i2 = ABS(Ki_a)**2

          ! Reflectivity of cloud droplets after Rayleigh
          IF (ldo_qc) THEN
            q_c = qc(i,j,k) * hlp
            IF (q_c < q_crit_radar%rayleigh) q_c = 0d0

            zh_local(i,j,k) = zh_local(i,j,k) + &
                    zradar_rayleigh_L_x(q_c, cloud%x_max, K_w2, ray_const%cloud_Zprefac_rho)
          END IF

          ! Reflectivity of rain after Rayleigh
          IF (ldo_qr) THEN
            q_r = qr(i,j,k) * hlp
            IF (q_r < q_crit_radar%rayleigh) q_r = 0d0

            zh_local(i,j,k) = zh_local(i,j,k) + &
                    zradar_rayleigh_L_Lexp(q_r,ray_const%rain1mom_expo,K_w2, &
                    ray_const%rain1mom_Zprefac_n0_rho)
          END IF

          IF (ldo_qi) THEN
            q_i = qi(i,j,k) * hlp
            IF (q_i < q_crit_radar%rayleigh) q_i = 0d0

            n_i = nice_mono_1mom(T_a)
            x_i = q_i / n_i
            x_i = MIN( MAX (x_i, ice%x_min), ice%x_max)

            IF (T_a > Tmeltbegin_i) THEN
              ! m, diameter of particle for calculation of degree of melting
              Dtmp = (x_i/ice%a_geo)**(1d0/ice%b_geo) * 2d0
              CALL degree_of_melting_xxx_single_Dref(&
                      T_a,Dtmp,Dexpo_i,Dref_i,&
                      Tmeltbegin_i,meltdegTmin_i,Tmin_f,&
                      REAL(Tmax_i(IND_IJ_2D),kind=dp),&
                      fmelt_mean)
              K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,&
                      MAX(MIN(fmelt_mean,1d0),0d0))
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_x(q_i, x_i, K2divRho2, ray_const%ice1mom_Zprefac)
            ELSE
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_x(q_i, x_i, K_i2, ray_const%ice1mom_Zprefac_rho)
            END IF
          END IF

          IF (ldo_qs) THEN
            q_s = qs(i,j,k) * hlp
            IF (q_s < q_crit_radar%rayleigh) q_s = 0d0

            IF (T_a > Tmeltbegin_s) THEN
              ! m, diameter of particle, accounting for calculation of average degree of melting
              Dtmp = D_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                          ray_const%snow_gam_mub1nu) * 2d0
              Dref = MIN(D_of_X_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                                   ray_const%snow_gam_mub1nu),&
                         6d-3)
              CALL degree_of_melting_xxx_single_Dref(&
                      T_a,Dtmp,Dexpo_s,Dref,&
                      Tmeltbegin_s,meltdegTmin_s,Tmin_f,&
                      REAL(Tmax_s(IND_IJ_2D),kind=dp),&
                      fmelt_mean)
              K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,&
                      MAX(MIN(fmelt_mean,1d0),0d0))
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_Lexp(q_s,ray_const%snow1mom_expo,K2divRho2,&
                      n0_s(i,j,k)**(1d0-ray_const%snow1mom_expo)*&
                      ray_const%snow1mom_Zprefac)
            ELSE
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_Lexp(q_s,ray_const%snow1mom_expo,K_i2,&
                      n0_s(i,j,k)**(1d0-ray_const%snow1mom_expo)*&
                      ray_const%snow1mom_Zprefac_rho)
            END IF
          END IF

          IF (ldo_qg) THEN
            q_g = qg(i,j,k) * hlp
            IF (q_g < q_crit_radar%rayleigh) q_g = 0d0

            Tmeltbegin_g_loc = REAL(Tmin_g(IND_IJ_2D),kind=dp)
            IF (T_a > Tmeltbegin_g_loc) THEN
              ! m, diameter of partical, accounted for calculation of medium degree of melting
              Dtmp = D_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                          ray_const%graupel_gam_mub1nu) * 2d0
              Dref = MIN(D_of_X_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                                   ray_const%graupel_gam_mub1nu),&
                         6d-3)
              ! In case of dynamic wetgrowth, scale down meltdegTmin_g so that
              ! it reaches 0 IF Tmin_g approaches Tmin_f = T0C_fwo:
              meltdegTmin_g_loc = meltdegTmin_g - &
                   meltdegTmin_g / MAX(Tmin_f-Tmeltbegin_g,1e-6_dp) * (Tmeltbegin_g_loc-Tmeltbegin_g)
              CALL degree_of_melting_xxx_single_Dref(&
                      T_a,Dtmp,Dexpo_g,Dref,&
                      Tmeltbegin_g_loc,meltdegTmin_g_loc,Tmin_f,&
                      REAL(Tmax_g(IND_IJ_2D),kind=dp),&
                      fmelt_mean)
              K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,&
                      MAX(MIN(fmelt_mean,1d0),0d0))
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_Lexp(q_g,ray_const%graupel1mom_expo,K2divRho2, &
                      ray_const%graupel1mom_Zprefac_n0)
            ELSE
              zh_local(i,j,k) = zh_local(i,j,k) + &
                      zradar_rayleigh_L_Lexp(q_g,ray_const%graupel1mom_expo,K_i2, &
                      ray_const%graupel1mom_Zprefac_n0_rho)
            END IF
          END IF

        END DO
      END DO
    END DO

    IF ( PRESENT(zh_radar) ) THEN
      zh_radar = zh_local
    ENDIF

    DEALLOCATE (zh_local, n0_s, m_i, m_w)

    RETURN
  END SUBROUTINE radar_rayleigh_oguchi_1mom_vec

  ! Interface routine to compute Vt (reflectivity weighted or number density weighted mean fallspeed of hydrometeors)
  !  on the model grid (version for the COSMO 1-moment schemes):
  SUBROUTINE vtradar_rayleigh_oguchi_1mom_vec(myproc, lambda_radar, &
       itype_gscp_fwo, isnow_n0temp, &
       Tmeltbegin_i, meltdegTmin_i, &
       Tmeltbegin_s, meltdegTmin_s, &
       Tmeltbegin_g, meltdegTmin_g, &
       rho, t, qc, qr, qi, qs, qg, &
       Tmax_i, Tmax_s, Tmax_g, &
       Tmin_g, &
       ilow, iup, jlow, jup, klow, kup, &
       lalloc_qi, lalloc_qs, lalloc_qg, &
       lwdbz, vt_radar, lhydrom_choice_testing)

    INTEGER,                 INTENT(IN)         :: itype_gscp_fwo, isnow_n0temp

    INTEGER,                 INTENT(IN)         :: myproc, ilow, iup, jlow, jup, klow, kup

    REAL(KIND=dp),           INTENT(IN)         :: lambda_radar, &
                                                   Tmeltbegin_i, meltdegTmin_i, &
                                                   Tmeltbegin_s, meltdegTmin_s, &
                                                   Tmeltbegin_g, meltdegTmin_g

    ! Model fields in model's wp, not dp:
    REAL(KIND=wp), DIMENSION(:,:,:), INTENT(IN) :: rho, t, qc, qr, qi, qs, qg
    ! Tmax for degree of melting in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmax_i, Tmax_s, Tmax_g
    ! Tmin for dynamic wet growth of graupel in wp:
    REAL(KIND=wp), DIMENSION(:,:)  , INTENT(IN) :: Tmin_g

    REAL(KIND=dp), INTENT(OUT)                  :: vt_radar(:,:,:)


    LOGICAL,                 INTENT(IN)         :: lwdbz, lalloc_qi, lalloc_qs, lalloc_qg
    LOGICAL, INTENT(in), OPTIONAL               :: lhydrom_choice_testing(6)

    ! Local Variables
    !----------------

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:) ::  vt_local, z_local, n0_s

    INTEGER :: i, j, k, ku, ko, ni, nj, nk
    REAL(KIND=dp) :: T_a, fmelt_mean, fmd, Dtmp, &
                     q_r, &
                     q_c, &
                     q_g, &
                     q_s, &
                     q_i, n_i, x_i, &
                     hlp, &
                     Tmeltbegin_g_loc, meltdegTmin_g_loc  ! for ldynamic_wetgrowth_gh
    REAL(KIND=dp) :: K_w2, K_i2, K2divRho2, Dref
    REAL(KIND=dp), PARAMETER :: Dexpo_i = 0.5d0, Dref_i  = 100d-6, &
                                Dexpo_s = 0.5d0, &
                                Dexpo_g = 0.6d0

    COMPLEX(kind=dp) :: Ki_a, Kw_a
    COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:,:,:) :: m_i, m_w

    LOGICAL :: ldo_qc, ldo_qr, ldo_qi, ldo_qs, ldo_qg

    CHARACTER (len=*), PARAMETER :: yzroutine = 'emvorado::vtradar_rayleigh_oguchi_1mom_vec()'

    REAL(KIND=dp)              :: lambda_radar_3digits
    CHARACTER(len=12)          :: clambda_radar

    ! .. Round lambda_radar to 3 significant digits, to avoid differences to slightly differing numeric
    !     lambda_radar values from namelists and DATA files:
    WRITE(clambda_radar, '(es12.2)') lambda_radar
    READ (clambda_radar, *)          lambda_radar_3digits
    
    ! .. Get field dimensions from input rho (assume that all hydrometeors and t are the same size):
    ni = SIZE(rho, DIM=1)
    nj = SIZE(rho, DIM=2)
    nk = SIZE(rho, DIM=3)

    ! .. Check field dimensions of output fields:
    IF (ni /= SIZE(vt_radar, dim=1) .OR. nj /= SIZE(vt_radar, dim=2)) THEN
      WRITE (*,*) 'ERROR in call to '//TRIM(yzroutine)//&
                  '(): Inconsistent field dimensions of rho and vt_radar!'
      STOP
    END IF

    ! .. Allocate local work arrays:
    ALLOCATE(vt_local(ni,nj,nk), z_local(ni,nj,nk))
    ALLOCATE(n0_s(ni,nj,nk))
    ALLOCATE(m_w(ni,nj,nk), m_i(ni,nj,nk))


    ! .. Initialize the summing variables with 0.0:
    !     vt_local: lambda^4/(pi^5*|K_w_0|^2) * \int sigma_b(D) v(D) N(D) dD  = Vt*Ze
    !      z_local: lambda^4/(pi^5*|K_w_0|^2) * \int sigma_b(D)      N(D) dD  = Ze
    vt_local = 0d0
    z_local = 0d0

    ku = klow
    ko = kup

    ! .. Initializations for PSD-parameters for cloud drops and cloud ice,
    !    since there are no such assumptions in the COSMO microphysics:


    IF (PRESENT(lhydrom_choice_testing)) THEN
      ldo_qc = lhydrom_choice_testing(1)
      ldo_qr = lhydrom_choice_testing(2)
      ldo_qi = lhydrom_choice_testing(3) .AND. lalloc_qi
      ldo_qs = lhydrom_choice_testing(4) .AND. lalloc_qs
      ldo_qg = lhydrom_choice_testing(5) .AND. lalloc_qg
    ELSE
      ldo_qc = .TRUE.
      ldo_qr = .TRUE.
      ldo_qi = lalloc_qi
      ldo_qs = lalloc_qs
      ldo_qg = lalloc_qg
    END IF

    ! On its first call or if otherwise necessary, this computes
    ! some constant prefactors needed below and stores it on global struct "ray_const":
    IF(lwdbz) THEN
      CALL init_radar_rayleigh_consts( &
           .TRUE., .TRUE., lalloc_qi, lalloc_qs, lalloc_qg, .FALSE., &
           DBLE(lambda_radar_3digits))
    END IF
    CALL init_radar_vt_oguchi ( &
           .TRUE., .TRUE., lalloc_qi, lalloc_qs, lalloc_qg, .FALSE., &
           lambda_radar_3digits, lwdbz)

!$OMP PARALLEL DO PRIVATE(j,k)
    DO k = 1, nk
      DO j = 1, nj
        ! Refractive index of water after Ray (1972), only valid for -10 < T < 30 degree C
        m_w(:,j,k) = m_complex_water_ray_vec(DBLE(lambda_radar_3digits),&
           MIN(MAX(t(:,j,k)-T0C_fwo,mw_Tmin),mw_Tmax), ni)
        ! Refractive index of ice after C. Maetzler, only valid for -80 < T < 0 degree C
        m_i(:,j,k) = m_complex_ice_maetzler_vec(DBLE(lambda_radar_3digits),&
           MIN(MAX(t(:,j,k)-T0C_fwo,mi_Tmin),mi_Tmax), ni)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Precalculate NWP-field dependent snow number density field (but only if we
    ! we consider snow at all):
    IF (ldo_qs) &
      CALL snow_1mom_n0(itype_gscp_fwo, isnow_n0temp, &
                        t, rho, qs, snow%a_geo, &
                        n0_s)


    IF (lwdbz) THEN

!$OMP PARALLEL DO PRIVATE(i,j,k,T_a,hlp,Kw_a,K_w2,Ki_a,K_i2,K2divRho2,&
!$OMP&                    q_r,q_c,q_i,n_i,x_i,q_s,q_g,fmelt_mean,fmd,Dtmp,&
!$OMP&                    Tmeltbegin_g_loc,meltdegTmin_g_loc)
      DO k= ku, ko
        DO j = jlow, jup
          DO i = ilow, iup

            T_a = t(i,j,k)
            hlp = rho(i,j,k)

            Kw_a = (m_w(i,j,k)**2-1d0)/(m_w(i,j,k)**2+2d0)
            K_w2 = ABS(Kw_a)**2
            Ki_a = (m_i(i,j,k)**2-1d0)/(m_i(i,j,k)**2+2d0)
            K_i2 = ABS(Ki_a)**2

            q_c = qc (i,j,k) * hlp

            ! Vt*reflectivity and reflectivity of cloud drops after Rayleigh:
            IF (ldo_qc) THEN
              IF (q_c < q_crit_radar%rayleigh) q_c = 0d0

              vt_local(i,j,k) = vt_local(i,j,k) + &
                      vtradar_pure_phase(cloud%x_max, vt_const%cloud_vtexpo,&
                      q_c*K_w2*vt_const%cloud_vtprefac)
              z_local(i,j,k) = z_local(i,j,k) + &
                      zradar_rayleigh_L_x(q_c, cloud%x_max, K_w2, ray_const%cloud_Zprefac_rho)
            END IF  ! ldo_qc

            ! Vt*reflectivity and reflectivity of rain after Rayleigh:
            IF (ldo_qr) THEN
              q_r = qr(i,j,k) * hlp
              IF (q_r < q_crit_radar%rayleigh) q_r = 0d0

              ! weighted with the reflectivity:
              vt_local(i,j,k) = vt_local(i,j,k) + &
                      vtradar_pure_phase(q_r, vt_const%rain_vtexpo, &
                      K_w2*vt_const%rain_vtprefac)
              z_local(i,j,k) = z_local(i,j,k) + &
                      zradar_rayleigh_L_Lexp(q_r,ray_const%rain1mom_expo,K_w2, &
                      ray_const%rain1mom_Zprefac_n0_rho)
            END IF  ! ldo_qr

            ! Vt*reflectivity and reflectivity of cloud ice after Rayleigh:
            IF (ldo_qi) THEN
              q_i = qi(i,j,k) * hlp
              IF (q_i < q_crit_radar%rayleigh) q_i = 0d0

              n_i = nice_mono_1mom(T_a)
              x_i = q_i / n_i
              x_i = MIN( MAX (x_i, ice%x_min), ice%x_max)

              IF (T_a > Tmeltbegin_i) THEN
                ! m, diameter of particle for calculation of degree of melting
                Dtmp = (x_i/ice%a_geo)**(1d0/ice%b_geo) * 2d0
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_i,Dref_i,&
                        Tmeltbegin_i,meltdegTmin_i,Tmin_f,&
                        REAL(Tmax_i(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,fmelt_mean)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase(x_i,&
                        vt_const%ice_vtexpo,vt_const%ice_vtfprefac,&
                        vt_const%ice_vtlexpo,vt_const%ice_vtlprefac,&
                        K2divRho2*q_i,fmd)
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_x(q_i, x_i, K2divRho2, ray_const%ice1mom_Zprefac)
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(x_i, vt_const%ice_vtexpo,&
                        q_i*K_i2*vt_const%ice_vtprefac)
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_x(q_i, x_i, K_i2, ray_const%ice1mom_Zprefac_rho)
              END IF  ! T_a > Tmeltbegin_i
            END IF  ! ldo_qi

            ! Vt*reflectivity and reflectivity of snow after Rayleigh:
            IF (ldo_qs) THEN
              q_s = qs(i,j,k) * hlp
              IF (q_s < q_crit_radar%rayleigh) q_s = 0d0

              IF (T_a > Tmeltbegin_s) THEN
                ! m, diameter of partical, accounted for calculation of medium degree of melting
                Dtmp = D_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                            vt_const%snow_gam_mub1nu) * 2d0
                Dref = MIN(D_of_X_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                                     vt_const%snow_gam_mub1nu),&
                           6d-3)
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_s,Dref,&
                        Tmeltbegin_s,meltdegTmin_s,Tmin_f,&
                        REAL(Tmax_s(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,fmelt_mean)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase( q_s,    &
                        vt_const%snow_vtexpo , vt_const%snow_vtfprefac*n0_s(i,j,k)**(1d0-vt_const%snow_vtexpo ), &
                        vt_const%snow_vtlexpo, vt_const%snow_vtlprefac*n0_s(i,j,k)**(1d0-vt_const%snow_vtlexpo), &
                        K2divRho2, fmd     )
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_Lexp( q_s, ray_const%snow1mom_expo, K2divRho2, &
                        n0_s(i,j,k)**(1d0-ray_const%snow1mom_expo) * ray_const%snow1mom_Zprefac)
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(q_s, vt_const%snow_vtexpo, &
                        K_i2*n0_s(i,j,k)**(1d0-vt_const%snow_vtexpo)*vt_const%snow_vtprefac )
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_Lexp(q_s,ray_const%snow1mom_expo, K_i2, &
                        n0_s(i,j,k)**(1d0-ray_const%snow1mom_expo)*ray_const%snow1mom_Zprefac_rho )
              END IF  ! T_a > Tmeltbegin_s
            END IF  ! ldo_qs

            ! Vt*reflectivity and reflectivity of graupel after Rayleigh:
            IF (ldo_qg) THEN
              q_g = qg(i,j,k) * hlp
              IF (q_g < q_crit_radar%rayleigh) q_g = 0d0
              !IF (q_g >= q_crit_radar%rayleigh) THEN

              Tmeltbegin_g_loc = REAL(Tmin_g(IND_IJ_2D),kind=dp)
              IF (T_a > Tmeltbegin_g_loc) THEN
                ! m, diameter of partical, accounted for calculation of medium degree of melting
                Dtmp = D_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                            vt_const%graupel_gam_mub1nu) * 2d0
                Dref = MIN(D_of_X_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                                     vt_const%graupel_gam_mub1nu),&
                                                     6d-3)
                ! In case of dynamic wetgrowth, scale down meltdegTmin_g so that
                ! it reaches 0 if Tmin_g approaches Tmin_f = T0C_fwo:
                meltdegTmin_g_loc = meltdegTmin_g - &
                     meltdegTmin_g / MAX(Tmin_f-Tmeltbegin_g,1e-6_dp) * (Tmeltbegin_g_loc-Tmeltbegin_g)
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_g,Dref,&
                        Tmeltbegin_g_loc,meltdegTmin_g_loc,Tmin_f,&
                        REAL(Tmax_g(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                K2divRho2 = K_rho_fac_oguchi(Ki_a,Kw_a,inv_rhoi,inv_rhow,fmelt_mean)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase(q_g,&
                        vt_const%graupel_vtexpo, vt_const%graupel_vtfprefac,&
                        vt_const%graupel_vtlexpo, vt_const%graupel_vtlprefac,&
                        K2divRho2,fmd)
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_Lexp(q_g,ray_const%graupel1mom_expo,K2divRho2, &
                        ray_const%graupel1mom_Zprefac_n0)
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(q_g, vt_const%graupel_vtexpo, &
                        K_i2*vt_const%graupel_vtprefac)
                z_local(i,j,k) = z_local(i,j,k) + &
                        zradar_rayleigh_L_Lexp(q_g,ray_const%graupel1mom_expo,K_i2, &
                        ray_const%graupel1mom_Zprefac_n0_rho)
              END IF  ! T_a > Tmeltbegin_g_loc
              !END IF
            END IF  ! ldo_qg

          END DO
        END DO
      END DO

!$OMP END PARALLEL DO
    ELSE    ! lwdbz

!$OMP PARALLEL DO PRIVATE(i,j,k,T_a,hlp,K_w2,K_i2,q_r,q_c,q_i,n_i,x_i,&
!$OMP&                    fmelt_mean,fmd,Dtmp,q_s,q_g,&
!$OMP&                    Tmeltbegin_g_loc,meltdegTmin_g_loc)
      DO k= ku, ko
        DO j = jlow, jup
          DO i = ilow, iup

            T_a = t(i,j,k)
            hlp = rho(i,j,k)

            K_w2 = ABS( (m_w(i,j,k)**2-1d0)/(m_w(i,j,k)**2+2d0) )**2
            K_i2 = ABS( (m_i(i,j,k)**2-1d0)/(m_i(i,j,k)**2+2d0) )**2

            q_c = qc (i,j,k) * hlp

            ! Vt*number dens. and number dens. of cloud drops after Rayleigh:
            IF (ldo_qc) THEN
              IF (q_c < q_crit_radar%rayleigh) q_c = 0d0

              vt_local(i,j,k) = vt_local(i,j,k) + &
                      vtradar_pure_phase(cloud%x_max,&
                      vt_const%cloud_vtexpo, q_c*vt_const%cloud_vtprefac)
              z_local(i,j,k)  = z_local(i,j,k) + &
                      q_c/cloud%x_max
            END IF  ! ldo_qc

            ! Vt*reflectivity and reflectivity of rain after Rayleigh:
            IF (ldo_qr) THEN
              q_r = qr(i,j,k) * hlp
              IF (q_r < q_crit_radar%rayleigh) q_r = 0d0

              ! not reflectivity weighted, simply the spectral mean fall speed:
              vt_local(i,j,k) = vt_local(i,j,k) + &
                      vtradar_pure_phase(q_r, vt_const%rain_vtexpo, &
                      vt_const%rain_vtprefac)
              z_local(i,j,k)  = z_local(i,j,k) + &
                      nradar(q_r, vt_const%rain_nexpo, &
                      vt_const%rain_nprefac)
            END IF  ! ldo_qr

            ! Vt*number dens. and number dens. of cloud ice after Rayleigh:
            IF (ldo_qi) THEN
              q_i = qi(i,j,k) * hlp
              IF (q_i < q_crit_radar%rayleigh) q_i = 0d0

              n_i = nice_mono_1mom(T_a)
              x_i = q_i / n_i
              x_i = MIN( MAX (x_i, ice%x_min), ice%x_max)
              n_i = q_i / x_i

              IF (T_a > Tmeltbegin_i) THEN
                ! m, diameter of particle for calculation of degree of melting
                Dtmp = (x_i/ice%a_geo)**(1d0/ice%b_geo) * 2d0
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_i,Dref_i,&
                        Tmeltbegin_i,meltdegTmin_i,Tmin_f,&
                        REAL(Tmax_i(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase(x_i,&
                        vt_const%ice_vtexpo,vt_const%ice_vtfprefac,&
                        vt_const%ice_vtlexpo,vt_const%ice_vtlprefac,&
                        q_i,fmd)
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(x_i, vt_const%ice_vtexpo, &
                        q_i*vt_const%ice_vtprefac)
              END IF  ! T_a > Tmeltbegin_i
              z_local(i,j,k) = z_local(i,j,k) + n_i
            END IF  ! ldo_qi


            ! Vt*number dens. and number dens. of snow after Rayleigh:
            IF (ldo_qs) THEN
              q_s = qs(i,j,k) * hlp
              IF (q_s < q_crit_radar%rayleigh) q_s = 0d0

              IF (T_a > Tmeltbegin_s) THEN
                ! m, diameter of particle, accounted for calculation of medium degree of melting
                Dtmp = D_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                            vt_const%snow_gam_mub1nu) * 2d0
                Dref = MIN(D_of_X_average_1M_exp_cin(q_s,snow,n0_s(i,j,k),&
                                                     vt_const%snow_gam_mub1nu),&
                           6d-3)
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_s,Dref,&
                        Tmeltbegin_s,meltdegTmin_s,Tmin_f,&
                        REAL(Tmax_s(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase(q_s,&
                        vt_const%snow_vtexpo , vt_const%snow_vtfprefac*n0_s(i,j,k)**(1d0-vt_const%snow_vtexpo ), &
                        vt_const%snow_vtlexpo, vt_const%snow_vtlprefac*n0_s(i,j,k)**(1d0-vt_const%snow_vtlexpo), &
                        1d0, fmd )
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(q_s, vt_const%snow_vtexpo,&
                        n0_s(i,j,k)**(1d0-vt_const%snow_vtexpo)*vt_const%snow_vtprefac)
              END IF  ! T_a > Tmeltbegin_s
              z_local(i,j,k) = z_local(i,j,k) + &
                      nradar(q_s, vt_const%snow_nexpo, &
                      n0_s(i,j,k)**(1d0-vt_const%snow_nexpo)*vt_const%snow_nprefac)
            END IF  ! ldo_qs

            ! Vt*number dens. and number dens. of graupel after Rayleigh:
            IF (ldo_qg) THEN
              q_g = qg(i,j,k) * hlp
              IF (q_g < q_crit_radar%rayleigh) q_g = 0d0

              Tmeltbegin_g_loc = REAL(Tmin_g(IND_IJ_2D),kind=dp)
              IF (T_a > Tmeltbegin_g_loc) THEN
                ! m, diameter of partical, accounted for calculation of medium degree of melting
                Dtmp = D_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                            vt_const%graupel_gam_mub1nu) * 2d0
                Dref = MIN(D_of_X_average_1M_exp_cin(q_g,graupel,graupel%n0_const,&
                                                     vt_const%graupel_gam_mub1nu),&
                           6d-3)
                ! In case of dynamic wetgrowth, scale down meltdegTmin_g so that
                ! it reaches 0 IF Tmin_g approaches Tmin_f = T0C_fwo:
                meltdegTmin_g_loc = meltdegTmin_g - &
                     meltdegTmin_g / MAX(Tmin_f-Tmeltbegin_g,1e-6_dp) * (Tmeltbegin_g_loc-Tmeltbegin_g)
                CALL degree_of_melting_xxx_single_Dref(&
                        T_a,Dtmp,Dexpo_g,Dref,&
                        Tmeltbegin_g_loc,meltdegTmin_g_loc,Tmin_f,&
                        REAL(Tmax_g(IND_IJ_2D),kind=dp),&
                        fmelt_mean)
                fmelt_mean = MAX(MIN(fmelt_mean,1d0),0d0)

                fmd = 0d0
                IF (fmelt_mean > quasi_zero) fmd = EXP(LOG(fmelt_mean)*vt_const%am)

                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_mixed_phase(q_g,&
                        vt_const%graupel_vtexpo, vt_const%graupel_vtfprefac,&
                        vt_const%graupel_vtlexpo, vt_const%graupel_vtlprefac,&
                        1d0,fmd)
              ELSE
                vt_local(i,j,k) = vt_local(i,j,k) + &
                        vtradar_pure_phase(q_g, vt_const%graupel_vtexpo, &
                        vt_const%graupel_vtprefac)
              END IF  ! T_a > Tmeltbegin_g_loc
              z_local(i,j,k) = z_local(i,j,k) + &
                      nradar(q_g, vt_const%graupel_nexpo, &
                      vt_const%graupel_nprefac)
            END IF  ! ldo_qg

          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF  ! lwdbz

    ! .. vt_radar now receives the final values for Vt, including density correction
    !    for large Reynolds numbers
    !    (NOTE: not correct for the contribution of cloud droplets):
    vt_radar = vtradar_normalize(vt_local, z_local, rho, rho_0, &
                                 ni, nj, nk, ilow, iup, jlow, jup, ku, ko)

    DEALLOCATE (vt_local, z_local, n0_s, m_i, m_w)

    RETURN
  END SUBROUTINE vtradar_rayleigh_oguchi_1mom_vec

END MODULE radar_mie_iface_cosmo_1mom
