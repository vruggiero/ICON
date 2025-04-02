!NEC$ options "-finline-max-function-size=200"

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

MODULE radar_mie_utils

!------------------------------------------------------------------------------
!
! Description: Utilitiy functions for the model specific interface(s)
!              (at the moment mainly COSMO) to the EMVORADO libraries for
!              radar reflectivity  calculations based on Mie-Theory.
!
! Method:
!   See subroutines below
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind , ONLY : dp, wp

  USE radar_data_mie, ONLY: &
       particle, particle_rain_coeffs, &
       rain, cloud, snow, ice, graupel, &
       t_mgd_params, rho_w_fwo, T0C_fwo, third, qeps, quasi_zero, &
       t_dbzlookuptable, t_tabparam, i_scal_log, i_scal_fscal, i_scal_lin, i_scal_dbz, &
       scal_crit_lut

  USE radar_gamma_functions_vec, ONLY: gfct

!===============================================================================

  IMPLICIT NONE

!===============================================================================

  PUBLIC

!===============================================================================
!===============================================================================

  PRIVATE :: scale_val_lut_scalar, scale_val_lut_1D, scale_val_lut_2D, scale_val_lut_3D, &
       &     descale_val_lut_scalar, descale_val_lut_1D, descale_val_lut_2D, descale_val_lut_3D
  
  INTERFACE scale_val_lut
    MODULE PROCEDURE           &
         scale_val_lut_scalar, &
         scale_val_lut_1D,     &
         scale_val_lut_2D,     &
         scale_val_lut_3D
  END INTERFACE scale_val_lut
  
  INTERFACE descale_val_lut
    MODULE PROCEDURE             &
         descale_val_lut_scalar, &
         descale_val_lut_1D,     &
         descale_val_lut_2D,     &
         descale_val_lut_3D
  END INTERFACE descale_val_lut
  
  INTERFACE scale_dval_lut
    MODULE PROCEDURE           &
         scale_dval_lut_scalar, &
         scale_dval_lut_1D,     &
         scale_dval_lut_2D,     &
         scale_dval_lut_3D
  END INTERFACE scale_dval_lut
  
!===============================================================================
!===============================================================================

CONTAINS

!===============================================================================
!===============================================================================


!======================================================================
!    Some simple hash functions from the internet ...
!======================================================================

  !========================================================================
  !
  !   Function for hash computation of an input double precision vector:
  !
  !   Commonly referred to as DJB hash (Daniel J. Bernstein) which
  !    computes a hash for character string c by the loop
  !
  !     h(0) = 5381
  !     do i=1, len(c)
  !       h(i) = h(i-1) * 33 + iachar(c(i))
  !     end do
  !
  !   Multiplication with 33 is equivalent to left bit shift by 5.
  !
  !========================================================================

  FUNCTION hash2(text) RESULT(hashed)
    IMPLICIT NONE
    CHARACTER (len=*), INTENT(in)        :: text
    INTEGER                              :: hashed
    INTEGER                              :: hashed_loc, shift_nr     ! 32-bit integers assumed

    INTEGER :: i

    hashed_loc = 5381
    shift_nr = 5
    DO i=1,LEN(text)
      hashed_loc = (ISHFT(hashed_loc,shift_nr) + hashed_loc) + IACHAR(text(i:i))
    END DO
    hashed = hashed_loc

  END FUNCTION hash2

!=========================================================================
!    Functions for computing hashes which represent different
!    reflectivity computation configurations.
!    Needed to create unique names for lookup table files.
!=========================================================================

  SUBROUTINE hash_radar_lookup ( &
       ctableconfig, chydroconfig, cbaseconfig, &
       hashed, hash_text)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in)  :: ctableconfig, chydroconfig, cbaseconfig

    INTEGER, INTENT(OUT)          :: hashed
    CHARACTER(len=*), INTENT(out) :: hash_text

    hash_text(:) = ' '
    WRITE(hash_text, '(3(a,3x))') &
         TRIM(ctableconfig), TRIM(chydroconfig), &
         TRIM(cbaseconfig)
    !WRITE(*,'(3A)') '   hash_text="',hash_text,'"'

    hashed = hash2(TRIM(hash_text))
    !WRITE(*,'(A,I12)') '   initial hashnr=',hashed
    IF (hashed < 0) hashed = hashed + HUGE(hashed)  ! ensure hash to be > 0
    !WRITE(*,'(A,I12)') '   updated hashnr=',hashed

  END SUBROUTINE hash_radar_lookup

  FUNCTION get_hashbase (&
       cversion_lt, lambda_radar, itype_gscp_fwo, &
       mur_is_fixed, isnow_n0temp, &
       itype_refl, igraupel_type, itype_Dref_fmelt,&
       Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
       Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
       Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g, &
       Tmeltbegin_h, meltdegTmin_h, Tmax_min_h, Tmax_max_h, &
       ctype_dryice, ctype_wetice, &
       ctype_drysnow, ctype_wetsnow, &
       ctype_drygraupel, ctype_wetgraupel, &
       ctype_dryhail, ctype_wethail) &
       RESULT (cbaseconfig)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(IN)    :: lambda_radar, &
                                    Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
                                    Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
                                    Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g, &
                                    Tmeltbegin_h, meltdegTmin_h, Tmax_min_h, Tmax_max_h
    INTEGER, INTENT(IN)          :: itype_gscp_fwo, itype_refl, igraupel_type, &
                                    itype_Dref_fmelt, isnow_n0temp
    LOGICAL, INTENT(in)          :: mur_is_fixed
    CHARACTER(len=*), INTENT(in) :: ctype_dryice, ctype_wetice, ctype_drysnow, ctype_wetsnow, &
                                    ctype_drygraupel, ctype_wetgraupel, &
                                    ctype_dryhail, ctype_wethail, cversion_lt

    CHARACTER(len=3000)          :: cbaseconfig

    CHARACTER(len=500)           :: cmiscconfig, cmeltconfig, cemaconfig

    cmiscconfig(:) = ' '
    IF (itype_gscp_fwo < 200) THEN
      WRITE(cmiscconfig, '(a,3x,es12.2,i6,4i3)') &
           TRIM(cversion_lt), lambda_radar, itype_gscp_fwo, &
           isnow_n0temp, itype_refl, igraupel_type, itype_Dref_fmelt
    ELSE
      WRITE(cmiscconfig, '(a,3x,es12.2,i6,L3,3i3)') &
           TRIM(cversion_lt), lambda_radar, itype_gscp_fwo, &
           mur_is_fixed, itype_refl, igraupel_type, itype_Dref_fmelt
    END IF

    cmeltconfig(:) = ' '
    IF (itype_gscp_fwo < 200) THEN
      WRITE(cmeltconfig, '(12es15.5)') &
           Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
           Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
           Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g
    ELSE
      WRITE(cmeltconfig, '(16es15.5)') &
           Tmeltbegin_i, meltdegTmin_i, Tmax_min_i, Tmax_max_i, &
           Tmeltbegin_s, meltdegTmin_s, Tmax_min_s, Tmax_max_s, &
           Tmeltbegin_g, meltdegTmin_g, Tmax_min_g, Tmax_max_g, &
           Tmeltbegin_h, meltdegTmin_h, Tmax_min_h, Tmax_max_h
    END IF

    cemaconfig(:) = ' '
    IF (itype_gscp_fwo < 200) THEN
       WRITE(cemaconfig, '(6(3x,a))') &
           TRIM(ADJUSTL(ctype_dryice)), TRIM(ADJUSTL(ctype_wetice)), &
           TRIM(ADJUSTL(ctype_drysnow)), TRIM(ADJUSTL(ctype_wetsnow)), &
           TRIM(ADJUSTL(ctype_drygraupel)), TRIM(ADJUSTL(ctype_wetgraupel))
    ELSE
      WRITE(cemaconfig, '(8(3x,a))') &
           TRIM(ADJUSTL(ctype_dryice)), TRIM(ADJUSTL(ctype_wetice)), &
           TRIM(ADJUSTL(ctype_drysnow)), TRIM(ADJUSTL(ctype_wetsnow)), &
           TRIM(ADJUSTL(ctype_drygraupel)), TRIM(ADJUSTL(ctype_wetgraupel)), &
           TRIM(ADJUSTL(ctype_dryhail)), TRIM(ADJUSTL(ctype_wethail))
    END IF

    cbaseconfig(:) = ' '
    WRITE(cbaseconfig,'(3(a,3x))') &
          TRIM(cmiscconfig), TRIM(cmeltconfig), TRIM(cemaconfig)

    RETURN
  END FUNCTION get_hashbase

  SUBROUTINE particle_assign(b,a)
    IMPLICIT NONE
    TYPE(particle), INTENT(in)  :: a
    TYPE(particle), INTENT(out) :: b
    b = a
  END SUBROUTINE particle_assign

!================================================================
!
! Size distribution parameters for the 1-moment scheme
!
!================================================================

  SUBROUTINE snow_1mom_n0(itype_gscp_fwo, isnow_n0temp, &
                          t, rho, qs, ageos, &
                          n0_s)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: itype_gscp_fwo, isnow_n0temp
    REAL(kind=wp), INTENT(in) :: t(:,:,:), rho(:,:,:), qs(:,:,:)
    REAL(kind=dp), INTENT(in) :: ageos

    REAL(kind=dp), INTENT(out) :: n0_s(:,:,:)

    INTEGER :: i, j, k, ni, nj, nk

!-----------------------------------------------------------------------------------

    ni = SIZE(rho, DIM=1)
    nj = SIZE(rho, DIM=2)
    nk = SIZE(rho, DIM=3)

    !.. N0-Parameter for snow:
!$OMP PARALLEL DO PRIVATE(i,j,k)
    DO k=1, nk
      IF (itype_gscp_fwo >= 140) THEN
        DO j=1, nj
          DO i=1, ni
            CALL calc_n0_snow(isnow_n0temp, REAL(t(i,j,k), kind=dp), &
                              REAL(qs(i,j,k)*rho(i,j,k), kind=dp), ageos, n0_s(i,j,k))
          END DO
        END DO
      END IF
    END DO
!$OMP END PARALLEL DO

    RETURN
  END SUBROUTINE snow_1mom_n0


  !-----------------------------------------------------------------------------------
  ! Subroutine to calculate n0_s as a function of temperature after Field et al. (2005)
  !
  ! Input:    T_a    Temperature in K
  !           q_s    mass density of snow in kg/m^3
  !           ageos  Parameter of mass_size-relation x=ageos*D^bgeos, x in kg, D in m
  ! Output:   n0_s   Intercept parameter of expon. snow size distrib. in m^-4
  !-----------------------------------------------------------------------------------

  SUBROUTINE calc_n0_snow(isnow_n0temp, T_a, q_s, ageos, n0_s)

    IMPLICIT NONE

    INTEGER, INTENT(in)       :: isnow_n0temp
    REAL(kind=dp), INTENT(in) :: T_a, q_s, ageos

    REAL(kind=dp), INTENT(out) :: n0_s

    REAL (KIND=dp) :: ztc, zn0s, &
                      alf, bet, hlp, m2s, m3s
    INTEGER        :: nn

    REAL (KIND=dp), PARAMETER  :: zn0s1 = 13.5_dp * 5.65E5_dp, & ! parameter in N0S(T)
                                  zn0s2 = -0.107_dp              ! parameter in N0S(T), Field et al

!-----------------------------------------------------------------------------------

    ! Snow: Coeffs for moment relation based on 2nd moment (Field et al. 2005)
    !        (defined as parameters because of better inlining and vectorization properties)
    REAL(KIND=dp), DIMENSION(10), PARAMETER  :: &
         mma = (/5.065339_dp, -0.062659_dp, -3.032362_dp, 0.029469_dp, -0.000285_dp, &
                 0.312550_dp,  0.000204_dp,  0.003199_dp, 0.000000_dp, -0.015952_dp /), &
         mmb = (/0.476221_dp, -0.015896_dp,  0.165977_dp, 0.007468_dp, -0.000141_dp, &
                 0.060366_dp,  0.000079_dp,  0.000594_dp, 0.000000_dp, -0.003577_dp /)

    !.. N0-Parameter for snow:
    n0_s = 1e9_dp  ! initialization

    IF (isnow_n0temp == 1) THEN
      ! Calculate n0s using the temperature-dependent formula of Field et al. (2005)
      ztc = T_a - T0C_fwo
      ztc = MAX(MIN(ztc,0.0_dp),-40.0_dp)
      zn0s = zn0s1*EXP(zn0s2*ztc)
      zn0s = MIN(zn0s,1e9_dp)
      zn0s = MAX(zn0s,1e6_dp)
      n0_s = zn0s
    ELSEIF (isnow_n0temp == 2) THEN
      ! Calculate n0s using the temperature-dependent moment relations of Field et al. (2005)
      ztc = T_a - T0C_fwo
      ztc = MAX(MIN(ztc,0.0_dp),-40.0_dp)

      nn  = 3
      hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
          + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
      alf = 10.0_dp**hlp
      bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
          + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
      IF (q_s >= quasi_zero) THEN
!!$ UB: caution! Here is the exponent bms=2.0 hardwired! not ideal!
        m2s = q_s / ageos
        m3s = alf*EXP(bet*LOG(m2s))
        hlp = zn0s1*EXP(zn0s2*ztc)
!!$ UB: the 13.5 is actually 3^3 / gamma(3) ...
        zn0s = 13.50_dp * m2s**4 / m3s**3
        zn0s = MAX(zn0s,0.5_dp*hlp)
        zn0s = MIN(zn0s,1e2_dp*hlp)
        zn0s = MIN(zn0s,1e9_dp)
        n0_s = MAX(zn0s,8e5_dp)
      ELSE
        n0_s = 1e9_dp
      END IF
    ELSE
      ! Old constant n0s
      n0_s = 8.0e5_dp
    ENDIF

    RETURN
  END SUBROUTINE calc_n0_snow

  !================================================================
  !
  ! Function to compute the gamma DSD shape parameter mu from
  !   the mu-Dm-relation for rain of Axel Seifert (2008).
  ! This is used in the Seifert-Beheng 2-moment scheme outside clouds.
  !
  !   - mu:  N(D) = N0 * D^mu * exp(-lam*D)
  !   - Dm: mean mass diameter
  !
  !================================================================

  FUNCTION mu_d_relation_seifert(x_r,rain,rain_coeffs) RESULT(mue)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)              :: x_r
    TYPE(particle), INTENT(in)             :: rain
    TYPE(particle_rain_coeffs), INTENT(in) :: rain_coeffs

    REAL(KIND=dp) :: mue

    REAL(KIND=dp) :: D_r

    ! use a/b_D instead of a/b_m, adapting to changes in init_twomom_types
    D_r = (x_r/rain%a_geo)**(1d0/rain%b_geo)

    IF (D_r <= rain_coeffs%cmu3) THEN ! see Seifert (2008)
      ! UB: cmu5 is 2.0, so for efficiency we replace it by INT 2 (which can be optimized by the compiler):
      mue = rain_coeffs%cmu0 * &
              TANH((4.0_dp*rain_coeffs%cmu2*(D_r-rain_coeffs%cmu3))**2) &
              + rain_coeffs%cmu4
    ELSE
      mue = rain_coeffs%cmu1 * &
              TANH((1.0_dp*rain_coeffs%cmu2*(D_r-rain_coeffs%cmu3))**2) &
            + rain_coeffs%cmu4
    ENDIF

    RETURN
  END FUNCTION mu_d_relation_seifert

  !================================================================
  !
  ! Function to compute the mean mass diameter Dm based
  !   on particle types in emvorado-notation, i.e.
  !   a/b_geo for x = a_geo * D^b_geo
  !
  !================================================================

  FUNCTION mean_mass_diameter (p, q, n) RESULT (dm)
    
    TYPE(particle), INTENT(in) :: p     ! a/b_geo as x = a_geo * D^b_geo
    REAL(KIND=dp),  INTENT(in) :: q, n  ! either kg/kg or kg/m3

    REAL(kind=dp) :: dm, a, b, x

    REAL(kind=dp), PARAMETER :: eps = 1e-20_dp

    a = (1.0d0/p%a_geo)**(1.0d0/p%b_geo)
    b = (1.0d0/p%b_geo)
    x = MIN( MAX(q/(n+eps), p%x_min), p%x_max)
    ! dm = a*x**b:
    dm = a * EXP(b*LOG(x))
    
  END FUNCTION mean_mass_diameter
  
  !================================================================
  !
  ! Function to compute the parameterized number of ice crystals
  !   as function of Temperature for the ICON 1-moment schemes.
  !
  !================================================================

  FUNCTION nice_mono_1mom (T) RESULT (qni)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: T       ! Temperature in K
    REAL(KIND=dp)              :: qni     ! number density of cloud ice particles

    ! lsuper_coolw is constant in both cloudice and graupel:
    LOGICAL      , PARAMETER   :: lsuper_coolw = .TRUE.
    REAL(KIND=dp), PARAMETER   :: znimax_Thom  = 250.0E3_dp  ! From Greg Thompson, 1/m^3
    REAL(KIND=dp), PARAMETER   :: zThn_1mom    = 236.15_dp ! temperature for hom. freezing of cloud water in K
    REAL(wp)     , PARAMETER   :: tmelt        = 273.15_dp ! melting temperature of ice/snow in K
    REAL(KIND=dp), PARAMETER   :: znimax       = 100.0_dp*EXP(0.2_dp*(Tmelt-zThn_1mom))

    !------------------------------------------------------------------------------

    ! Number of activate ice crystals:
    IF (lsuper_coolw) THEN
      qni    = MIN( 5.0_dp  *EXP(0.304_dp*MAX((Tmelt-T),0.0_dp) ), znimax_Thom)  ! 1/m^3
    ELSE
      qni    = MIN( 100.0_dp*EXP(0.2_dp  *(Tmelt-T) )            , znimax)       ! 1/m^3
    END IF

  END FUNCTION nice_mono_1mom

  !=======================================================================
  !
  ! some help functions
  !
  !=======================================================================

  ! NOT USED, BUT KEPT FOR REFERENCE:
  !
  ! Function to compute the imom'th-moment M_imom from the 0'th and 1'st moments N and L
  !  of a generalized gamma distribution, depending on the shape parameters:
  !
  !  M_imom =  gamfac_imom * L * (L/N)^(imom-1)
  !
  FUNCTION gamfac_imom (parti,imom) RESULT(zg)
    IMPLICIT NONE
    TYPE(particle) :: parti
    REAL(KIND=dp) :: imom
    REAL(KIND=dp) :: zg
    zg =     gfct((imom+parti%nu+1.0_dp)/parti%mu)        &
         & / gfct((     parti%nu+1.0_dp)/parti%mu)        &
         & * ( gfct((parti%nu+1.0_dp)/parti%mu)      &
         &   / gfct((parti%nu+2.0_dp)/parti%mu) )**imom
  END FUNCTION gamfac_imom

  ! THE FOLLOWING IS ACTUALLY USED:
  !
  ! Function to compute the (m(or x)-based) imom'th-moment M_imom from the
  ! 0'th and 1'st moments N and L of a generalized gamma distribution given in m(or x)-space,
  ! depending on the shape parameters given in D(!)-space
  !
  !  M_imom =  gamfac_imom * L * (L/N)^(imom-1)
  !
  FUNCTION gamfac_imom_DMGD (parti,imom) RESULT(zg)

    IMPLICIT NONE

    TYPE(particle), INTENT(in) :: parti
    REAL(KIND=dp), INTENT(in)  :: imom

    REAL(KIND=dp) :: zg, gfct_tmp

    gfct_tmp = gfct( (parti%mu+1.0_dp) / parti%nu )
    zg       = gfct((parti%mu+imom*parti%b_geo+1.0_dp)/parti%nu) / gfct_tmp * &
             ( gfct_tmp / gfct((parti%mu+parti%b_geo+1.0_dp)/parti%nu) )**imom
  END FUNCTION gamfac_imom_DMGD

  ! Medial diameter of hydro meteors with gen. Gamma distribution
  FUNCTION D_average_DMGD (q, qn, parti) RESULT (D)

    ! valid for D-based parti parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, qn
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: bm, fac, x

    bm = 1d0/parti%b_geo
    fac = gamfac_imom_DMGD(parti,bm)

    !x = MIN(MAX((q/(qn+1d-20)),1d-12),1.0_dp)
    x = MIN(MAX(q/(qn+quasi_zero),parti%x_min),parti%x_max)

    D = fac * (x/parti%a_geo)**bm

  END FUNCTION D_average_DMGD

  FUNCTION D_average_cin (q, qn, parti, fac) RESULT (D)

    ! valid for m(or x)-based parti parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, qn, fac
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: x

    !x = MIN(MAX((q/(qn+1d-20)),1d-12),1.0_dp)
    x = MIN(MAX(q/(qn+quasi_zero),parti%x_min),parti%x_max)

    D = fac * parti%a_geo * x ** parti%b_geo

  END FUNCTION D_average_cin

  FUNCTION D_average_cin_DMGD (x, parti, fac) RESULT (D)

    ! valid for D-based parti parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: x, fac
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D

    D = fac * (x/parti%a_geo)**(1d0/parti%b_geo)

  END FUNCTION D_average_cin_DMGD

  ! vectorized version:
  FUNCTION D_average_vec_DMGD (q, qn, parti, anz) RESULT (D)

    ! valid for D-based parti parameters

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: anz
    REAL(KIND=dp), INTENT(in)  :: q(anz), qn(anz)
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D(anz)
    REAL(KIND=dp) :: bm, fac, x(anz)

    bm = 1d0/parti%b_geo
    fac = gamfac_imom_DMGD(parti,bm)

    !x = MIN(MAX((q/(qn+1d-20)),1d-12),1.0_dp)
    x = MIN(MAX(q/(qn+quasi_zero),parti%x_min),parti%x_max)

    D = fac * (x/parti%a_geo)**bm

  END FUNCTION D_average_vec_DMGD

  FUNCTION D_average_1M_exp (q, parti, n0) RESULT (D)

    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, n0
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: fac

    fac = 1.0_dp / (parti%a_geo * n0 * gfct(parti%b_geo+1.0_dp))

    D = (MAX(q,qeps) * fac) ** (1.0_dp/(1.0_dp+parti%b_geo))

  END FUNCTION D_average_1M_exp

  FUNCTION D_average_1M_exp_cin (q, parti, n0, gfc_bgeop1) RESULT (D)

    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, n0, gfc_bgeop1
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: fac

!!$   gfc_bgeop1 =  gfct(b_geo+1d0))
    fac = 1.0_dp / (parti%a_geo * n0 * gfc_bgeop1)
    D = EXP(LOG(MAX(q,qeps) * fac ) * (1.0_dp/(1.0_dp+parti%b_geo)))

  END FUNCTION D_average_1M_exp_cin

  ! vectorized version:
  FUNCTION D_average_1M_exp_vec (q, parti, n0, anz) RESULT (D)

    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: anz
    REAL(KIND=dp), INTENT(in)  :: q(anz), n0
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D(anz)
    REAL(KIND=dp) :: fac

    fac = 1.0_dp / (parti%a_geo * n0 * gfct(parti%b_geo+1.0_dp))

    D = (MAX(q,qeps) * fac) ** (1.0_dp/(1.0_dp+parti%b_geo))

  END FUNCTION D_average_1M_exp_vec

  FUNCTION D_average_1M_exp_vec_n0s (q, parti, n0, anz) RESULT (D)

    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: anz
    REAL(KIND=dp), INTENT(in)  :: q(anz), n0(anz)
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D(anz)
    REAL(KIND=dp) :: fac(anz)

    fac = 1.0_dp / (parti%a_geo * n0 * gfct(parti%b_geo+1.0_dp))

    D = (MAX(q,qeps) * fac) ** (1.0_dp/(1.0_dp+parti%b_geo))

  END FUNCTION D_average_1M_exp_vec_n0s

  FUNCTION D_of_X_average_1M_exp (q, parti, n0) RESULT (D)

    ! Diameter of mean mass particle
    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, n0
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: fac, expo

    expo = 1.0_dp/(parti%b_geo*(parti%b_geo+1.0_dp))
    fac = gfct(parti%b_geo+1.0_dp) ** expo / (parti%a_geo*n0)**(parti%b_geo*expo)
    D = MAX(q,qeps)**(parti%b_geo*expo) * fac

  END FUNCTION D_of_X_average_1M_exp

  FUNCTION D_of_X_average_1M_exp_cin (q, parti, n0, gfc_bgeop1) RESULT (D)

    ! Diameter of mean mass particle
    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: q, n0, gfc_bgeop1
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D
    REAL(KIND=dp) :: fac, expo

!!$ gfc_bgeop1 = gfct(b_geo+1.0_dp)

    expo = 1.0_dp/(parti%b_geo*(parti%b_geo+1.0_dp))
    fac = gfc_bgeop1 ** expo / (parti%a_geo*n0)**(parti%b_geo*expo)
    D = MAX(q,qeps)**(parti%b_geo*expo) * fac

  END FUNCTION D_of_X_average_1M_exp_cin

  ! vectorized version:
  FUNCTION D_of_X_average_1M_exp_vec (q, parti, n0, anz) RESULT (D)

    ! Diameter of mean mass particle
    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: anz
    REAL(KIND=dp), INTENT(in)  :: q(anz), n0
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D(anz)
    REAL(KIND=dp) :: fac, expo

    expo = 1.0_dp/(parti%b_geo*(parti%b_geo+1.0_dp))
    fac = gfct(parti%b_geo+1.0_dp) ** expo / (parti%a_geo*n0)**(parti%b_geo*expo)
    D = MAX(q,qeps)**(parti%b_geo*expo) * fac

  END FUNCTION D_of_X_average_1M_exp_vec

  FUNCTION D_of_X_average_1M_exp_vec_n0s (q, parti, n0, anz) RESULT (D)

    ! Diameter of mean mass particle
    ! NOTE: ONLY valid for standard exponential distribution (mu=0, nu=1)
    ! valid for D-based microphysical parameters

    IMPLICIT NONE

    INTEGER, INTENT(in)        :: anz
    REAL(KIND=dp), INTENT(in)  :: q(anz), n0(anz)
    TYPE(PARTICLE), INTENT(in) :: parti

    REAL(KIND=dp) :: D(anz)
    REAL(KIND=dp) :: fac(anz), expo

    expo = 1.0_dp/(parti%b_geo*(parti%b_geo+1.0_dp))
    fac = gfct(parti%b_geo+1.0_dp) ** expo / (parti%a_geo*n0)**(parti%b_geo*expo)
    D = MAX(q,qeps)**(parti%b_geo*expo) * fac

  END FUNCTION D_of_X_average_1M_exp_vec_n0s

  ! Coefficient vector for Simpson rule
  FUNCTION simpson_weights(n_stuetz) RESULT(simpson)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n_stuetz

    REAL(KIND=dp) :: simpson(n_stuetz+1)

    REAL(KIND=dp), PARAMETER, DIMENSION(3) :: basis = (/third, 4d0*third, third/)
    INTEGER :: i

    simpson = 0.0_dp
    DO i=1,n_stuetz-1,2
      simpson(i:i+2) = simpson(i:i+2) + basis
    END DO

  END FUNCTION simpson_weights

  ! PDF value vector of MGD
  FUNCTION fd_mgd(mgd,D,nD) RESULT (fd)

    IMPLICIT NONE

    INTEGER,            INTENT(in) :: nD
    REAL(KIND=dp),      INTENT(in) :: D(nD)
    TYPE(t_mgd_params), INTENT(in) :: mgd

    REAL(KIND=dp) :: fd(nD)

    !f_d = n_0 * D_w**mu * EXP(-lam * D_w**nu)
    fd = mgd%n0 * D**mgd%mu * EXP(-mgd%lam * D**mgd%nu)

    RETURN
  END FUNCTION fd_mgd

  ! Integration weights vector of MGD
  FUNCTION integweights_mgd(mgd,D,dD,simpson,nD) RESULT (itw)

    IMPLICIT NONE

    INTEGER,            INTENT(in) :: nD
    REAL(KIND=dp),      INTENT(in) :: dD, D(nD), simpson(nD)
    TYPE(t_mgd_params), INTENT(in) :: mgd

    REAL(KIND=dp) :: itw(nD)

    !eta = SUM(f_d * CBACK * simpson * dD)
    itw = fd_mgd(mgd,D,nD) * simpson * dD

    RETURN
  END FUNCTION integweights_mgd

  ! Calculate MGD parameters for a 1-mom particle and additional n0
  FUNCTION  mgd_1mom(parti,L,n_0) RESULT (mgd)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)  :: L, n_0
    TYPE(particle), INTENT(in) :: parti

    TYPE(t_mgd_params) :: mgd

    REAL(kind=dp) :: tmp1, tmp2

    mgd%mu  = parti%mu
    mgd%nu  = parti%nu

    mgd%n0  = n_0

    tmp1 = (mgd%mu+1.0_dp) / mgd%nu
    tmp2 = (mgd%mu+parti%b_geo+1.0_dp) / mgd%nu

    mgd%lam = ( (parti%a_geo*mgd%n0*gfct(tmp2)) / (L*mgd%nu) )**(1.0_dp/tmp2)

    mgd%q = L
    mgd%qn = mgd%n0 * gfct(tmp1) / (mgd%nu * mgd%lam**tmp1)

    RETURN
  END FUNCTION mgd_1mom

  ! Calculate MGD parameters for a 2-mom particle
  FUNCTION mgd_2mom(parti,L,N,muD,nuD) RESULT (mgd)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)           :: L, N
    REAL(KIND=dp), OPTIONAL, INTENT(in) :: muD, nuD
    TYPE(particle), INTENT(in)          :: parti

    TYPE(t_mgd_params) :: mgd

    REAL(kind=dp) :: tmp1, tmp2, gfct_tmp1, gfct_tmp2

    IF (PRESENT(nuD) .AND. PRESENT(muD)) THEN
      mgd%mu = muD
      mgd%nu = nuD
    ELSE
      mgd%mu = parti%mu
      mgd%nu = parti%nu
    END IF

    tmp1 = (mgd%mu+1.0_dp) / mgd%nu
    tmp2 = (mgd%mu+parti%b_geo+1.0_dp) / mgd%nu
    gfct_tmp1 = gfct(tmp1)
    gfct_tmp2 = gfct(tmp2)

    mgd%lam = ( (parti%a_geo*N*gfct_tmp2) / (L*gfct_tmp1) )**(mgd%nu/parti%b_geo)
    mgd%n0  = N * mgd%nu * mgd%lam**tmp1 / gfct_tmp1

    mgd%q  = L
    mgd%qn = N

    RETURN
  END FUNCTION mgd_2mom

  SUBROUTINE set_dummy_mgd(mgd)
    IMPLICIT NONE

    TYPE(t_mgd_params), INTENT(inout) :: mgd

    mgd%n0  = 0.0_dp
    mgd%lam = 1.0_dp ! a non-zero number

    mgd%q   = 0.0_dp
    mgd%qn  = 0.0_dp

    RETURN
  END SUBROUTINE set_dummy_mgd

  !==================================================================================
  !==================================================================================
  !
  ! Generating full control strings for the degree of melting computations
  !  from the shorthand notation from the namelist input
  !
  !==================================================================================
  !==================================================================================

  SUBROUTINE decode_controlstring_3C(ctype_relrefind, &
       mixingrulestring, matrixstring, inclusionstring, &
       hoststring, hostmatrixstring, hostinclusionstring, &
       ctype_default)

    IMPLICIT NONE

    CHARACTER(len=6), INTENT(in) :: ctype_relrefind, ctype_default
    CHARACTER(len=*), INTENT(out) :: &
         mixingrulestring, matrixstring, inclusionstring, &
         hoststring, hostmatrixstring, hostinclusionstring

    INTEGER                 :: i
    CHARACTER(len=6) :: buf(2)

    buf(1:2) = REPEAT(' ',6)
    buf(1) = ctype_default(1:6)
    buf(2) = ctype_relrefind(1:6)

    ! Ausgabestrings initialisieren:
    mixingrulestring    = REPEAT(' ',LEN(mixingrulestring))
    matrixstring        = REPEAT(' ',LEN(matrixstring))
    inclusionstring     = REPEAT(' ',LEN(inclusionstring))
    hoststring          = REPEAT(' ',LEN(hoststring))
    hostmatrixstring    = REPEAT(' ',LEN(hostmatrixstring))
    hostinclusionstring = REPEAT(' ',LEN(hostinclusionstring))

    ! Im ersten Durchlauf defaults setzen, im 2. Durchlauf die eigentlichen Strings:
    DO i=1,2

      SELECT CASE (buf(i)(1:1))
      CASE('m','M')
        mixingrulestring = 'maxwellgarnett'
      CASE('o','O','w','W')
        mixingrulestring = 'oguchi'
      CASE('b','B')
        mixingrulestring = 'bruggemann'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 1. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(2:2))
      CASE('n','N')
        hoststring = 'none'
      CASE('a','A')
        hoststring = 'air'
      CASE('i','I')
        hoststring = 'ice'
      CASE('w','W')
        hoststring = 'water'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 2. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(3:3))
      CASE('a','A')
        matrixstring = 'air'
      CASE('i','I')
        matrixstring = 'ice'
      CASE('w','W')
        matrixstring = 'water'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 3. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(4:4))
      CASE('s','S')
        inclusionstring = 'spheroidal'
      CASE('k','K')
        inclusionstring = 'spherical'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 4. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(5:5))
      CASE('a','A')
        hostmatrixstring = 'air'
      CASE('i','I')
        hostmatrixstring = 'ice'
      CASE('w','W')
        hostmatrixstring = 'water'
      CASE('m','M')
        hostmatrixstring = 'icewater'
      CASE('s','S')
        hostmatrixstring = 'airice'
      CASE('r','R')
        hostmatrixstring = 'airwater'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 5. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(6:6))
      CASE('s','S')
        hostinclusionstring = 'spheroidal'
      CASE('k','K')
        hostinclusionstring = 'spherical'
      CASE default
        WRITE (*,*) '  decode_controlstrings_3C: invalid ctype_relrefind, 6. character: ', TRIM(buf(i))
      END SELECT

    END DO

    RETURN
  END SUBROUTINE decode_controlstring_3C

  SUBROUTINE decode_controlstring_2C(ctype_relrefind, &
       mixingrulestring, matrixstring, inclusionstring, &
       ctype_default)

    IMPLICIT NONE

    CHARACTER(len=3), INTENT(in) :: ctype_relrefind, ctype_default
    CHARACTER(len=*), INTENT(out) :: &
         mixingrulestring, matrixstring, inclusionstring

    INTEGER                 :: i
    CHARACTER(len=3) :: buf(2)

    buf(1:2) = REPEAT(' ',3)
    buf(1) = ctype_default(1:3)
    buf(2) = ctype_relrefind(1:3)

    ! initilize output strings
    mixingrulestring    = REPEAT(' ',LEN(mixingrulestring))
    matrixstring        = REPEAT(' ',LEN(matrixstring))
    inclusionstring     = REPEAT(' ',LEN(inclusionstring))

    ! first set default values, but correct values within second loop
    DO i=1,2

      SELECT CASE (buf(i)(1:1))
      CASE('m','M')
        mixingrulestring = 'maxwellgarnett'
      CASE('o','O','w','W')
        mixingrulestring = 'oguchi'
      CASE('b','B')
        mixingrulestring = 'bruggemann'
      CASE default
        WRITE (*,*) '  decode_controlstrings_2C: invalid ctype_relrefind, 1. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(2:2))
      CASE('a','A')
        matrixstring = 'air'
      CASE('i','I')
        matrixstring = 'ice'
      CASE('w','W')
        matrixstring = 'water'
      CASE default
        WRITE (*,*) '  decode_controlstrings_2C: invalid ctype_relrefind, 2. character: ', TRIM(buf(i))
      END SELECT

      SELECT CASE (buf(i)(3:3))
      CASE('s','S')
        inclusionstring = 'spheroidal'
      CASE('k','K')
        inclusionstring = 'spherical'
      CASE default
        WRITE (*,*) '  decode_controlstrings_2C: invalid ctype_relrefind, 3. character: ', TRIM(buf(i))
      END SELECT

    END DO

    RETURN
  END SUBROUTINE decode_controlstring_2C

  
  !================================================================================================
  !================================================================================================
  ! Procedures for polarimetric lookup tables
  !================================================================================================
  !================================================================================================

  
  !================================================================================================
  !
  ! Constructor for type dbzlookuptable:
  !
  ! Allocate mem space, associate parameter pointers, define table vectors from table boundaries,
  ! and define scaling for the equidistant qi table nodes.
  !
  ! However, does not compute table values and does not set any scaling flags for actual interpolations!
  ! This has to be done by subsequent calls to init_scaling_tabparams() for each parameter.
  !
  !================================================================================================

  SUBROUTINE init_dbzlookuptable (lut, nqi, nTa, nTm, qilow, qiup, Talow, Taup, Tmlow, Tmup, &
       &                          flag_qi_scal_eq, f_eq, ierr)

    TYPE(t_dbzlookuptable), INTENT(inout) :: lut  ! The table to initialize
    INTEGER,  INTENT(in) :: nqi, nTa, nTm       ! The dimensions of the table
    REAL(dp), INTENT(in) :: qilow, qiup, Talow, Taup, Tmlow, Tmup  ! min/max values of table vectors. low- and up-value should be different, except for Tm, where equality is tolerated
    INTEGER,  INTENT(in) :: flag_qi_scal_eq     ! flag for scaling type of equidistant qi_axis (i_scal_log, i_scal_fscal, or i_scal_lin, i_scal_dbz)
    REAL(dp), INTENT(in) :: f_eq                ! scaling exponent for the equidistant qi-axis, > 0.0!
    INTEGER,  INTENT(out):: ierr                ! error flag on return, a value /= 0 inicates an error

    CHARACTER(len=*), PARAMETER :: yzroutine = 'emvorado::init_dbzlookuptable'
    INTEGER :: i

    ierr = 0

    ! catch some obvious errors or problems:
    IF (lut%is_initialized) THEN
      WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': lookup table for '//TRIM(lut%chydrotype)//&
           ' already initialized when it should not be! Programming error!'
      ierr = 1
      RETURN
    END IF

    ! wrong input for Ta-boundaries of the table:
    IF (Taup-Talow < -quasi_zero) THEN
      WRITE (*,'(a,f0.2,a,f0.2,a)') 'ERROR '//TRIM(yzroutine)//' initializing lookup table for '//&
           TRIM(lut%chydrotype)//': Taup=',Taup,' (upper T-value of table) is smaller than Talow=',Talow,' (lower T-value)!'
      ierr = 3
      RETURN
    END IF
    
    ! wrong input for Tmax-boundaries for melting scheme:
    IF (Tmup-Tmlow < -quasi_zero) THEN
      WRITE (*,'(a,f0.2,a,f0.2,a)') 'ERROR '//TRIM(yzroutine)//' initializing lookup table for '//&
           TRIM(lut%chydrotype)//': Tmup=',Tmup,' (''Tmax_max'') is smaller than Tmlow=',Tmlow,' (''Tmax_min'')!'
      ierr = 2
      RETURN
    END IF

    ! number of table nodes (one more to include initial values qi0, Ta0 and Tm0)
    lut%nqi = nqi + 1
    IF (Taup-Talow < quasi_zero .AND. nTa > 1) THEN
      ! reset table node number to 1, because there is effectively only one Ta-value:
      lut%nTa = 1
      WRITE (*,'(a,f0.2,a,f0.2,a)') 'WARNING '//TRIM(yzroutine)//' initializing lookup table for '//&
           TRIM(lut%chydrotype)//': upper and lower T-value of table are equal (',Talow,'), so lut%nTa is set to 1!'
    ELSE
      lut%nTa = nTa + 1
    END IF
    IF (Tmup-Tmlow < quasi_zero .AND. nTm > 1) THEN
      ! reset table node number to 1, because there is effectively only one Tm-value:
      lut%nTm = 1
      WRITE (*,'(a,f0.2,a,f0.2,a)') 'WARNING '//TRIM(yzroutine)//' initializing lookup table for '//&
           TRIM(lut%chydrotype)//': Tmup is equal to Tmlow (',Tmlow,'), so lut%nTm is set to 1!'
    ELSE
      lut%nTm = nTm + 1
    END IF
    
    ! flag and exponent for scaling of the equidistant q_i-axis:
    lut%flag_qi_scal_eq = flag_qi_scal_eq
    lut%f_eq  = f_eq
    lut%if_eq = 1d0 / f_eq

    IF (.NOT.ALLOCATED(lut%q_i)) THEN
      ALLOCATE(lut%q_i(lut%nqi))
      ALLOCATE(lut%q_i_lin(lut%nqi))
      ALLOCATE(lut%T_a(lut%nTa))
      ALLOCATE(lut%T_m(lut%nTm))
    ELSE
      ierr = 2
      lut%is_initialized = .FALSE.
      WRITE (*,*) 'ERROR: q_i in lookup table for '//TRIM(lut%chydrotype)//&
           ' already allocated when it should not be! Programming error!'
      RETURN
    END IF

    ! memory blocks for table values and their derivatives w.r.t. q_i
    IF (.NOT.ASSOCIATED(lut%mem)) THEN
      ALLOCATE(lut%mem (lut%nqi,lut%nTa,lut%nTm,lut%nparams))
    ELSE
      ierr = 3
      lut%is_initialized = .FALSE.
      WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': mem in lookup table for '//TRIM(lut%chydrotype)//&
           ' already allocated when it should not be! Programming error!'
      RETURN
    END IF
    IF (.NOT.ASSOCIATED(lut%dmem)) THEN
      ALLOCATE(lut%dmem(lut%nqi,lut%nTa,lut%nTm,lut%nparams))
    ELSE
      ierr = 4
      lut%is_initialized = .FALSE.
      WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': dmem in lookup table for '//TRIM(lut%chydrotype)//&
           ' already allocated when it should not be! Programming error!'
      RETURN
    END IF
    IF (.NOT.ASSOCIATED(lut%qmem)) THEN
      ALLOCATE(lut%qmem(lut%nqi,lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%flag_qi_scal)) THEN
      ALLOCATE(lut%flag_qi_scal(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%f_scal_qi)) THEN
      ALLOCATE(lut%f_scal_qi(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%if_scal_qi)) THEN
      ALLOCATE(lut%if_scal_qi(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%flag_mem_scal)) THEN
      ALLOCATE(lut%flag_mem_scal(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%f_scal_mem)) THEN
      ALLOCATE(lut%f_scal_mem(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%if_scal_mem)) THEN
      ALLOCATE(lut%if_scal_mem(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF
    IF (.NOT.ASSOCIATED(lut%flag_mem_interp_qi)) THEN
      ALLOCATE(lut%flag_mem_interp_qi(lut%nparams))
      ! don't bother about already associated here. Has been already checked for mem and dmem.
    END IF

    ! Associate pointers in the t_tabparam types to the memory blocks for each parameter:
    lut%zh%q_i     => lut%qmem(:,lut%i_zh  )
    lut%ah%q_i     => lut%qmem(:,lut%i_ah  )
    lut%zv%q_i     => lut%qmem(:,lut%i_zv  )
    lut%rrhv%q_i   => lut%qmem(:,lut%i_rrhv)
    lut%irhv%q_i   => lut%qmem(:,lut%i_irhv)
    lut%kdp%q_i    => lut%qmem(:,lut%i_kdp )
    lut%adp%q_i    => lut%qmem(:,lut%i_adp )
    lut%zvh%q_i    => lut%qmem(:,lut%i_zvh )

    lut%zh%val    => lut%mem(:,:,:,lut%i_zh  )
    lut%ah%val    => lut%mem(:,:,:,lut%i_ah  )
    lut%zv%val    => lut%mem(:,:,:,lut%i_zv  )
    lut%rrhv%val  => lut%mem(:,:,:,lut%i_rrhv)
    lut%irhv%val  => lut%mem(:,:,:,lut%i_irhv)
    lut%kdp%val   => lut%mem(:,:,:,lut%i_kdp )
    lut%adp%val   => lut%mem(:,:,:,lut%i_adp )
    lut%zvh%val   => lut%mem(:,:,:,lut%i_zvh )

    lut%zh%dval   => lut%dmem(:,:,:,lut%i_zh  )
    lut%ah%dval   => lut%dmem(:,:,:,lut%i_ah  )
    lut%zv%dval   => lut%dmem(:,:,:,lut%i_zv  )
    lut%rrhv%dval => lut%dmem(:,:,:,lut%i_rrhv)
    lut%irhv%dval => lut%dmem(:,:,:,lut%i_irhv)
    lut%kdp%dval  => lut%dmem(:,:,:,lut%i_kdp )
    lut%adp%dval  => lut%dmem(:,:,:,lut%i_adp )
    lut%zvh%dval  => lut%dmem(:,:,:,lut%i_zvh )

    lut%zh%qi_scal   => lut%flag_qi_scal(lut%i_zh  )
    lut%ah%qi_scal   => lut%flag_qi_scal(lut%i_ah  )
    lut%zv%qi_scal   => lut%flag_qi_scal(lut%i_zv  )
    lut%rrhv%qi_scal => lut%flag_qi_scal(lut%i_rrhv)
    lut%irhv%qi_scal => lut%flag_qi_scal(lut%i_irhv)
    lut%kdp%qi_scal  => lut%flag_qi_scal(lut%i_kdp )
    lut%adp%qi_scal  => lut%flag_qi_scal(lut%i_adp )
    lut%zvh%qi_scal  => lut%flag_qi_scal(lut%i_zvh )
    
    lut%zh%f_qi   => lut%f_scal_qi(lut%i_zh  )
    lut%ah%f_qi   => lut%f_scal_qi(lut%i_ah  )
    lut%zv%f_qi   => lut%f_scal_qi(lut%i_zv  )
    lut%rrhv%f_qi => lut%f_scal_qi(lut%i_rrhv)
    lut%irhv%f_qi => lut%f_scal_qi(lut%i_irhv)
    lut%kdp%f_qi  => lut%f_scal_qi(lut%i_kdp )
    lut%adp%f_qi  => lut%f_scal_qi(lut%i_adp )
    lut%zvh%f_qi  => lut%f_scal_qi(lut%i_zvh )
    
    lut%zh%if_qi   => lut%if_scal_qi(lut%i_zh  )
    lut%ah%if_qi   => lut%if_scal_qi(lut%i_ah  )
    lut%zv%if_qi   => lut%if_scal_qi(lut%i_zv  )
    lut%rrhv%if_qi => lut%if_scal_qi(lut%i_rrhv)
    lut%irhv%if_qi => lut%if_scal_qi(lut%i_irhv)
    lut%kdp%if_qi  => lut%if_scal_qi(lut%i_kdp )
    lut%adp%if_qi  => lut%if_scal_qi(lut%i_adp )
    lut%zvh%if_qi  => lut%if_scal_qi(lut%i_zvh )

    lut%zh%f_val   => lut%f_scal_mem(lut%i_zh  )
    lut%ah%f_val   => lut%f_scal_mem(lut%i_ah  )
    lut%zv%f_val   => lut%f_scal_mem(lut%i_zv  )
    lut%rrhv%f_val => lut%f_scal_mem(lut%i_rrhv)
    lut%irhv%f_val => lut%f_scal_mem(lut%i_irhv)
    lut%kdp%f_val  => lut%f_scal_mem(lut%i_kdp )
    lut%adp%f_val  => lut%f_scal_mem(lut%i_adp )
    lut%zvh%f_val  => lut%f_scal_mem(lut%i_zvh )
    
    lut%zh%if_val   => lut%if_scal_mem(lut%i_zh  )
    lut%ah%if_val   => lut%if_scal_mem(lut%i_ah  )
    lut%zv%if_val   => lut%if_scal_mem(lut%i_zv  )
    lut%rrhv%if_val => lut%if_scal_mem(lut%i_rrhv)
    lut%irhv%if_val => lut%if_scal_mem(lut%i_irhv)
    lut%kdp%if_val  => lut%if_scal_mem(lut%i_kdp )
    lut%adp%if_val  => lut%if_scal_mem(lut%i_adp )
    lut%zvh%if_val  => lut%if_scal_mem(lut%i_zvh )
    
    lut%zh  %val_scal => lut%flag_mem_scal(lut%i_zh  )
    lut%ah  %val_scal => lut%flag_mem_scal(lut%i_ah  )
    lut%zv  %val_scal => lut%flag_mem_scal(lut%i_zv  )
    lut%rrhv%val_scal => lut%flag_mem_scal(lut%i_rrhv)
    lut%irhv%val_scal => lut%flag_mem_scal(lut%i_irhv)
    lut%kdp %val_scal => lut%flag_mem_scal(lut%i_kdp )
    lut%adp %val_scal => lut%flag_mem_scal(lut%i_adp )
    lut%zvh %val_scal => lut%flag_mem_scal(lut%i_zvh )
    
    lut%zh  %val_interp => lut%flag_mem_interp_qi(lut%i_zh  )
    lut%ah  %val_interp => lut%flag_mem_interp_qi(lut%i_ah  )
    lut%zv  %val_interp => lut%flag_mem_interp_qi(lut%i_zv  )
    lut%rrhv%val_interp => lut%flag_mem_interp_qi(lut%i_rrhv)
    lut%irhv%val_interp => lut%flag_mem_interp_qi(lut%i_irhv)
    lut%kdp %val_interp => lut%flag_mem_interp_qi(lut%i_kdp )
    lut%adp %val_interp => lut%flag_mem_interp_qi(lut%i_adp )
    lut%zvh %val_interp => lut%flag_mem_interp_qi(lut%i_zvh )
    
    ! initial and final value for lookup table vector of q_i in the desired equidistant scaling:
    SELECT CASE (lut%flag_qi_scal_eq)
    CASE (i_scal_fscal,i_scal_log,i_scal_lin,i_scal_dbz)
      lut%qi0          = scale_val_lut (qilow, lut%flag_qi_scal_eq, lut%f_eq)
      lut%q_i(lut%nqi) = scale_val_lut (qiup , lut%flag_qi_scal_eq, lut%f_eq)
    CASE default
      WRITE (*,*) 'ERROR '//TRIM(yzroutine)//': wrong flag_qi_scal_eq for '//TRIM(lut%chydrotype)//&
           ', programming error! STOP!'
      STOP
    END SELECT
    ! initial and final value for lookup table vectors of T_a and T_m in linear scaling:
    lut%Ta0          = Talow
    lut%T_a(lut%nTa) = Taup
    lut%Tm0          = Tmlow
    lut%T_m(lut%nTm) = Tmup

    ! calculate increments in the equidistant scaling:
    lut%dqi = (lut%q_i(lut%nqi) - lut%qi0) / DBLE(lut%nqi-1)
    lut%dTa = (lut%T_a(lut%nTa) - lut%Ta0) / DBLE(MAX(lut%nTa-1,1))
    lut%dTm = (lut%T_m(lut%nTm) - lut%Tm0) / DBLE(MAX(lut%nTm-1,1))
    ! invert these increments:
    lut%idqi = 1d0 / lut%dqi
    IF (ABS(lut%dTa) >= quasi_zero) THEN
      lut%idTa = 1d0 / lut%dTa
    ELSE
      lut%idTa = 0d0
    END IF
    IF (ABS(lut%dTm) >= quasi_zero) THEN
      lut%idTm = 1d0 / lut%dTm
    ELSE
      lut%idTm = 0d0
    END IF

    ! create lookup table vectors in the equidistant scaling:
    DO i=1, lut%nqi
      ! q_i axis for computing the scaled equidistant interpolation interval:
      lut%q_i(i) = lut%qi0 + DBLE(i-1)*lut%dqi
    END DO
    ! store q_i also in linear space:
    lut%q_i_lin(:) = descale_val_lut (lut%q_i, lut%flag_qi_scal_eq, lut%if_eq)

    ! T_a and T_m are always in linear scaling:
    DO i=1, lut%nTa
      lut%T_a(i) = lut%Ta0 + DBLE(i-1)*lut%dTa
    END DO
    DO i=1, lut%nTm
      lut%T_m(i) = lut%Tm0 + DBLE(i-1)*lut%dTm
    END DO

  END SUBROUTINE init_dbzlookuptable

  !================================================================================================
  !
  ! Procedure to initialize scaling and interpolation flags for one lookup table parameter in type t_tabparam.
  !
  ! Must to be called only after init_dbzlookuptable()!
  !
  !================================================================================================

  SUBROUTINE init_scaling_tabparams (tpar, lut, flag_scal_qi, f_scal_qi, flag_scal_val, f_scal_val, flag_interp_val_qi)

    TYPE (t_tabparam),       INTENT(inout) :: tpar
    TYPE (t_dbzlookuptable), INTENT(in)    :: lut
    INTEGER,                 INTENT(in)    :: flag_scal_qi, flag_scal_val, flag_interp_val_qi
    REAL(dp),                INTENT(in)    :: f_scal_qi, f_scal_val

    ! Define interpolation scaling for each parameter w.r.t. qi (w.r.t. Ta will be always linear):
    ! --------------------------------------------------------------------------------------------

    !  1) scaling of the q-axis and store the scaled qi-axis values:
    !      (%f_qi and %if_qi ONLY relevant if %qi_scal = i_scal_fscal)

    tpar % qi_scal = flag_scal_qi
    tpar % f_qi    = f_scal_qi
    tpar % if_qi   = 1.0_dp / f_scal_qi
    tpar % q_i     = scale_val_lut (lut%q_i_lin, flag_scal_qi, f_scal_qi)
    
    !  2) scaling of the parameter itself:
    !      (%f_val and %if_val ONLY relevant if %val_scal = i_scal_fscal)
    
    tpar % val_scal = flag_scal_val
    tpar % f_val    = f_scal_val
    tpar % if_val   = 1.0_dp / f_scal_val

    ! 3) flag for type of interpolation of the scaled parameter with respect to the scaled qi:
    tpar % val_interp = flag_interp_val_qi
    
  END SUBROUTINE init_scaling_tabparams
  
  !================================================================================================
  !
  ! Module procedure functions to scale and descale table values based on scaling flags
  !  i_scal_log, i_scal_fscal, or i_scal_lin.
  !
  !================================================================================================

  FUNCTION scale_val_lut_scalar (val, flag_scal, f_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: f_scal
    REAL(dp)             :: outval

    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      outval = LOG10(MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = val**f_scal
    CASE (i_scal_lin)
      outval = val
    CASE (i_scal_dbz)
      outval = 10.0_dp*LOG10(MAX(val,quasi_zero))
    END SELECT

  END FUNCTION scale_val_lut_scalar

  FUNCTION scale_val_lut_1D (val, flag_scal, f_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: f_scal
    REAL(dp)             :: outval(SIZE(val,dim=1))
    
    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      outval = LOG10(MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = val**f_scal
    CASE (i_scal_lin)
      outval = val
    CASE (i_scal_dbz)
      outval = 10.0_dp*LOG10(MAX(val,quasi_zero))
    END SELECT

  END FUNCTION scale_val_lut_1D

  FUNCTION scale_val_lut_2D (val, flag_scal, f_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: f_scal
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2))
    
    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      outval = LOG10(MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = val**f_scal
    CASE (i_scal_lin)
      outval = val
    CASE (i_scal_dbz)
      outval = 10.0_dp*LOG10(MAX(val,quasi_zero))
    END SELECT

  END FUNCTION scale_val_lut_2D

  FUNCTION scale_val_lut_3D (val, flag_scal, f_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:,:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: f_scal
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2),SIZE(val,dim=3))
    
    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      outval = LOG10(MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = val**f_scal
    CASE (i_scal_lin)
      outval = val
    CASE (i_scal_dbz)
      outval = 10.0_dp*LOG10(MAX(val,quasi_zero))
    END SELECT

  END FUNCTION scale_val_lut_3D

  FUNCTION descale_val_lut_scalar (val, flag_scal, if_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: if_scal
    REAL(dp)             :: outval

    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)

    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      IF (val >= scal_crit_lut(flag_scal)) THEN
        outval = EXP(val*ln10)
      ELSE
        outval = 0.0_dp
      END IF
    CASE (i_scal_fscal)
      IF (val >= scal_crit_lut(flag_scal)) THEN
        outval = EXP(if_scal*LOG(val)) ! val**if_scal
      ELSE
        outval = 0.0_dp
      END IF
    CASE (i_scal_lin)
      outval = val
    CASE (i_scal_dbz)
      IF (val >= scal_crit_lut(flag_scal)) THEN
        outval = EXP(val*0.1_dp*ln10)
      ELSE
        outval = 0.0_dp
      END IF
    END SELECT

  END FUNCTION descale_val_lut_scalar

  FUNCTION descale_val_lut_1D (val, flag_scal, if_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: if_scal
    REAL(dp)             :: outval(SIZE(val,dim=1))

    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    INTEGER              :: i

    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      DO i=1, SIZE(val,dim=1)
        IF (val(i) >= scal_crit_lut(flag_scal)) THEN
          outval(i) = EXP(val(i)*ln10)
        ELSE
          outval(i) = 0.0_dp
        END IF
      END DO
    CASE (i_scal_fscal)
      DO i=1, SIZE(val,dim=1)
        IF (val(i) >= scal_crit_lut(flag_scal)) THEN
          outval(i) = EXP(if_scal*LOG(val(i))) ! val**if_scal
        ELSE
          outval(i) = 0.0_dp
        END IF
      END DO
    CASE (i_scal_lin)
      outval(:) = val(:)
    CASE (i_scal_dbz)
      DO i=1, SIZE(val,dim=1)
        IF (val(i) >= scal_crit_lut(flag_scal)) THEN
          outval(i) = EXP(val(i)*0.1_dp*ln10)
        ELSE
          outval(i) = 0.0_dp
        END IF
      END DO
    END SELECT

  END FUNCTION descale_val_lut_1D

  FUNCTION descale_val_lut_2D (val, flag_scal, if_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: if_scal
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2))

    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    INTEGER              :: i, j

    SELECT CASE (flag_scal)
    CASE (i_scal_log)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          IF (val(i,j) >= scal_crit_lut(flag_scal)) THEN
            outval(i,j) = EXP(val(i,j)*ln10)
          ELSE
            outval(i,j) = 0.0_dp
          END IF
        END DO
      END DO
    CASE (i_scal_fscal)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          IF (val(i,j) >= scal_crit_lut(flag_scal)) THEN
            outval(i,j) = EXP(if_scal*LOG(val(i,j))) ! val**if_scal
          ELSE
            outval(i,j) = 0.0_dp
          END IF
        END DO
      END DO
    CASE (i_scal_lin)
      outval(:,:) = val(:,:)
    CASE (i_scal_dbz)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          IF (val(i,j) >= scal_crit_lut(flag_scal)) THEN
            outval(i,j) = EXP(val(i,j)*0.1_dp*ln10)
          ELSE
            outval(i,j) = 0.0_dp
          END IF
        END DO
      END DO
    END SELECT

  END FUNCTION descale_val_lut_2D

  FUNCTION descale_val_lut_3D (val, flag_scal, if_scal) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:,:)
    INTEGER,  INTENT(in) :: flag_scal
    REAL(dp), INTENT(in) :: if_scal
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2),SIZE(val,dim=3))

    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    INTEGER              :: i, j, k

    SELECT CASE (flag_scal)
    CASE (i_scal_log)
!$omp parallel do private(i,j,k)
     DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            IF (val(i,j,k) >= scal_crit_lut(flag_scal)) THEN
              outval(i,j,k) = EXP(val(i,j,k)*ln10)
            ELSE
              outval(i,j,k) = 0.0_dp
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do
    CASE (i_scal_fscal)
!$omp parallel do private(i,j,k)
      DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            IF (val(i,j,k) >= scal_crit_lut(flag_scal)) THEN
              outval(i,j,k) = EXP(if_scal*LOG(val(i,j,k))) ! val**if_scal
            ELSE
              outval(i,j,k) = 0.0_dp
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do
    CASE (i_scal_lin)
      outval(:,:,:) = val(:,:,:)
    CASE (i_scal_dbz)
!$omp parallel do private(i,j,k)
      DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            IF (val(i,j,k) >= scal_crit_lut(flag_scal)) THEN
              outval(i,j,k) = EXP(val(i,j,k)*0.1_dp*ln10)
            ELSE
              outval(i,j,k) = 0.0_dp
            END IF
          END DO
        END DO
      END DO
!$omp end parallel do
    END SELECT

  END FUNCTION descale_val_lut_3D

  !================================================================================================
  !
  ! Module procedure functions to scale and descale derivatives of table values based on scaling flags
  !  i_scal_log, i_scal_fscal, or i_scal_lin.
  !
  ! Scaled p : f(p)   = val
  ! Scaled qi: m(qi)  = qi_scaled
  ! dval     : unscaled dp/dqi
  ! if_scal_qi : inverse of f_scal_qi
  !
  ! Compute df(p(qi(m)))/dm = df/dp * dp/dqi * dqi/dm = df/dp * dval * dqi/dm
  ! where df/dp may depend on val itself.
  !
  ! For the 1D, 2D and 3D versions: the first dimension must be qi!
  !
  !================================================================================================

  FUNCTION scale_dval_lut_scalar (val, dval, qi_scaled, flag_scal_qi, if_scal_qi, flag_scal_val, f_scal_val) RESULT (outval)

    REAL(dp), INTENT(in) :: val, dval, qi_scaled
    INTEGER,  INTENT(in) :: flag_scal_qi, flag_scal_val
    REAL(dp), INTENT(in) :: if_scal_qi, f_scal_val
    REAL(dp)             :: outval

    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    REAL(dp), PARAMETER  :: ln10_10 = LOG(10.0_dp) * 0.1_dp

    ! .. compute df/dp * dval: 
    SELECT CASE (flag_scal_val)
    CASE (i_scal_log)
      outval = dval / (ln10*MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = dval*f_scal_val*val**(f_scal_val-1.0_dp)
    CASE (i_scal_lin)
      outval = dval
    CASE (i_scal_dbz)
      outval = dval / (ln10_10*MAX(val,quasi_zero))
    END SELECT

    ! .. compute dqi/dm:
    SELECT CASE (flag_scal_qi)
    CASE (i_scal_log)
      outval = outval * ln10*EXP(ln10*qi_scaled)
    CASE (i_scal_fscal)
      outval = outval * if_scal_qi * qi_scaled**(if_scal_qi-1.0_dp)
    CASE (i_scal_lin)
      outval = outval
    CASE (i_scal_dbz)
      outval = outval * ln10_10*EXP(ln10_10*qi_scaled)
    END SELECT

  END FUNCTION scale_dval_lut_scalar

  FUNCTION scale_dval_lut_1D (val, dval, qi_scaled, flag_scal_qi, if_scal_qi, flag_scal_val, f_scal_val) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:), dval(:)  ! vector dimension must be qi
    REAL(dp), INTENT(in) :: qi_scaled(:) 
    INTEGER,  INTENT(in) :: flag_scal_qi, flag_scal_val
    REAL(dp), INTENT(in) :: if_scal_qi, f_scal_val
    REAL(dp)             :: outval(SIZE(val,dim=1))
    
    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    REAL(dp), PARAMETER  :: ln10_10 = LOG(10.0_dp) * 0.1_dp

    SELECT CASE (flag_scal_val)
    CASE (i_scal_log)
      outval = dval / (ln10*MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = dval*f_scal_val*val**(f_scal_val-1.0_dp)
    CASE (i_scal_lin)
      outval = dval
    CASE (i_scal_dbz)
      outval = dval / (ln10_10*MAX(val,quasi_zero))
    END SELECT

    ! .. compute dqi/dm: (assuming that qi is the vector dimension of val and dval:
    SELECT CASE (flag_scal_qi)
    CASE (i_scal_log)
      outval = outval * ln10*EXP(ln10*qi_scaled)
    CASE (i_scal_fscal)
      outval = outval * if_scal_qi * qi_scaled**(if_scal_qi-1.0_dp)
    CASE (i_scal_lin)
      outval = outval
    CASE (i_scal_dbz)
      outval = outval * ln10_10*EXP(ln10_10*qi_scaled)
    END SELECT

  END FUNCTION scale_dval_lut_1D

  FUNCTION scale_dval_lut_2D (val, dval, qi_scaled, flag_scal_qi, if_scal_qi, flag_scal_val, f_scal_val) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:), dval(:,:)  ! first dimension must be qi
    REAL(dp), INTENT(in) :: qi_scaled(:)
    INTEGER,  INTENT(in) :: flag_scal_qi, flag_scal_val
    REAL(dp), INTENT(in) :: if_scal_qi, f_scal_val
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2))
    
    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    REAL(dp), PARAMETER  :: ln10_10 = LOG(10.0_dp) * 0.1_dp

    INTEGER              :: i, j

    SELECT CASE (flag_scal_val)
    CASE (i_scal_log)
      outval = dval / (ln10*MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = dval*f_scal_val*val**(f_scal_val-1.0_dp)
    CASE (i_scal_lin)
      outval = dval
    CASE (i_scal_dbz)
      outval = dval / (ln10_10*MAX(val,quasi_zero))
    END SELECT

    ! .. compute dqi/dm: (assuming that qi is the vector dimension of val and dval:
    SELECT CASE (flag_scal_qi)
    CASE (i_scal_log)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          outval(i,j) = outval(i,j) * ln10*EXP(ln10*qi_scaled(i))
        END DO
      END DO
    CASE (i_scal_fscal)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          outval(i,j) = outval(i,j) * if_scal_qi * qi_scaled(i)**(if_scal_qi-1.0_dp)
        END DO
      END DO
    CASE (i_scal_lin)
      CONTINUE
    CASE (i_scal_dbz)
      DO j=1, SIZE(val,dim=2)
        DO i=1, SIZE(val,dim=1)
          outval(i,j) = outval(i,j) * ln10_10*EXP(ln10_10*qi_scaled(i))
        END DO
      END DO
    END SELECT

  END FUNCTION scale_dval_lut_2D

  FUNCTION scale_dval_lut_3D (val, dval, qi_scaled, flag_scal_qi, if_scal_qi, flag_scal_val, f_scal_val) RESULT (outval)

    REAL(dp), INTENT(in) :: val(:,:,:), dval(:,:,:)  ! first dimension must be qi
    REAL(dp), INTENT(in) :: qi_scaled(:)
    INTEGER,  INTENT(in) :: flag_scal_qi, flag_scal_val
    REAL(dp), INTENT(in) :: if_scal_qi, f_scal_val
    REAL(dp)             :: outval(SIZE(val,dim=1),SIZE(val,dim=2),SIZE(val,dim=3))
    
    REAL(dp), PARAMETER  :: ln10 = LOG(10.0_dp)
    REAL(dp), PARAMETER  :: ln10_10 = LOG(10.0_dp) * 0.1_dp

    INTEGER              :: i, j, k

    SELECT CASE (flag_scal_val)
    CASE (i_scal_log)
      outval = dval / (ln10*MAX(val,quasi_zero))
    CASE (i_scal_fscal)
      outval = dval*f_scal_val*val**(f_scal_val-1.0_dp)
    CASE (i_scal_lin)
      outval = dval
    CASE (i_scal_dbz)
      outval = dval / (ln10_10*MAX(val,quasi_zero))
    END SELECT

    ! .. compute dqi/dm: (assuming that qi is the vector dimension of val and dval:
    SELECT CASE (flag_scal_qi)
    CASE (i_scal_log)
!$omp parallel do private(k,j,i)
      DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            outval(i,j,k) = outval(i,j,k) * ln10*EXP(ln10*qi_scaled(i))
          END DO
        END DO
      END DO
!$omp end parallel do
    CASE (i_scal_fscal)
!$omp parallel do private(k,j,i)
      DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            outval(i,j,k) = outval(i,j,k) * if_scal_qi * qi_scaled(i)**(if_scal_qi-1.0_dp)
          END DO
        END DO
      END DO
!$omp end parallel do
    CASE (i_scal_lin)
      CONTINUE
    CASE (i_scal_dbz)
!$omp parallel do private(k,j,i)
      DO k=1, SIZE(val,dim=3)
        DO j=1, SIZE(val,dim=2)
          DO i=1, SIZE(val,dim=1)
            outval(i,j,k) = outval(i,j,k) * ln10_10*EXP(ln10_10*qi_scaled(i))
          END DO
        END DO
      END DO
!$omp end parallel do
    END SELECT

  END FUNCTION scale_dval_lut_3D

END MODULE radar_mie_utils


