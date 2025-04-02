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

! JM changes (adapting to JSnyder offline changes):
!
! - added S1, S2 as output parameters to subroutines
!     ...
! - in soubroutines
!     SPHERE_SCATTER_BH_VEC,
!     SPHERE_SCATTER_BH,
!     COATEDSPHERE_SCATTER_BH_VEC,
!     COATEDSPHERE_SCATTER_BH,
!     COATEDSPHERE_SCATTER_VEC,
!     COATEDSPHERE_SCATTER
!     RAYLEIGH_SOAK_WETGR_VEC
!     RAYLEIGH_SOAK_WETGR
!     RAYLEIGH_DRYGR_VEC
!     RAYLEIGH_DRYGR
!       - replaced C* (xsecs) with f* (scattering amplitudes) as output parameters
!       - downgrade C* to local variables
!       - added S1,S2 as local variables
! - in subroutines
!     MIE_WATERSPH_WETHAIL_VEC
!     MIE_WATERSPH_WETHAIL
!     MIE_WATERSPH_WETHAIL_BH_VEC
!     MIE_WATERSPH_WETHAIL_BH
!     MIE_SPONGY_WETHAIL_VEC
!     MIE_SPONGY_WETHAIL
!     MIE_SPONGY_WETHAIL_BH_VEC
!     MIE_SPONGY_WETHAIL_BH
!     MIE_DRYHAIL_VEC
!     MIE_DRYHAIL
!     MIE_WATERSPH_WETGR_VEC
!     MIE_WATERSPH_WETGR
!     MIE_WATERSPH_WETGR_BH_VEC
!     MIE_WATERSPH_WETGR_BH
!     MIE_SOAK_TWOSPH_WETGR_VEC
!     MIE_SOAK_TWOSPH_WETGR
!     MIE_SOAK_TWOSPH_WETGR_BH_VEC
!     MIE_SOAK_TWOSPH_WETGR_BH
!     MIE_MEAN_WETGR_VEC
!     MIE_MEAN_WETGR
!     MIE_SOAK_WETGR_VEC
!     MIE_SOAK_WETGR
!     MIE_DRYSNOW_VEC
!     MIE_DRYSNOW
!     MIE_DRYSNOW_TWOSPH_VEC
!     MIE_DRYSNOW_TWOSPH
!     MIE_DRYGR_VEC
!     MIE_DRYGR
!     MIE_WETSNOW_TWOSPH_VEC
!     MIE_WETSNOW_TWOSPH
!       - replaced C* (xsecs) with f* (scattering amplitudes) as output parameters
!       - downgrade C* to local variables
!       - added T-matrix branch
!       - added aspectratio as input parameter
!
! - added subroutines
!     SPHEROID_SCATTER_TMAT_VEC,
!     SPHEROID_SCATTER_TMAT,
!     COATEDSPHEROID_SCATTER_TMAT_VEC,
!     COATEDSPHEROID_SCATTER_TMAT

MODULE radar_mielib_vec

!------------------------------------------------------------------------------
!
! Description:  Library for Radar Reflectivity Calculations based on Mie-Theory
!               (calculates radar reflectivity of precipitation particles using
!               full mie theorie which includes scattering by wet and melting particles)
!
!               Applicable to microwave radiation, wavelength > 1 mm
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind, ONLY : dp

  USE radar_data, ONLY : miss_value, miss_value_rhv

  USE radar_sphbessel_vec, ONLY: CSPHJYVEC, CSPHJY

  USE radar_data_mie, ONLY :   &
       rho_w_fwo,   &  ! density of liquid water
       rho_ice_fwo, &  ! density of ice          (kg/m^3)
       T0C_fwo,     &  ! T[K] at 0°C
       pi_dp, pih => pih_dp, pi6 => pi6_dp, inv_pi6 => inv_pi6_dp, &
       nci => nci_dp, degrad_dp, &
       m_air, &
       third, quasi_zero, Deps, xeps

  USE radar_dualpol_t_matrix_mod, ONLY : tmatrix

  USE radar_dualpol_t_matrix2_mod, ONLY : tmatrix2

  IMPLICIT NONE

  PUBLIC
  INTEGER, PARAMETER, PRIVATE :: NANG = 2

CONTAINS

  !----------------------------------------------------------------------------------------!
  ! 1) Routinen zur Berechnung des Komplexen Brechungsindex fuer Wasser und Eis
  !    als Funktion von Wellenlaenge und Temperatur nach verschiedenen
  !    Parametrisierungen aus der Literatur.
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! 1a) Komplexer Brechungsindex für Flüssigwasser nach Ray (1972)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_water_ray_vec(lambda,T,anz)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Wasser als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.001,1.0] m; T aus [-10.0,30.0] Grad C
    ! nach Ray (1972)

    INTEGER, INTENT(in)       :: anz
    REAL(KIND=dp), INTENT(in) :: T(anz), lambda

    COMPLEX(kind=dp) :: m_complex_water_ray_vec(anz)

    INTEGER                     :: k
    REAL(KIND=dp)               :: epsinf,epss,epsr,epsi,alpha,lambdas,nenner,Toff

    DO k=1,anz
      Toff = T(k) - 25d0
      epsinf  = 5.27137d0 + 0.02164740d0 * T(k)  &
           &              - 0.00131198d0 * T(k)*T(k)
      epss    = 78.54d+0 * (1d0 - 4.579d-3 * Toff  &
           &                      + 1.190d-5 * Toff*Toff  &
           &                      - 2.800d-8 * Toff*Toff*Toff)
      alpha   = -16.8129d0/(T(k)+T0C_fwo) + 0.0609265d0
      lambdas = 0.00033836d0 * EXP(2513.98d0/(T(k)+T0C_fwo)) * 1d-2

      nenner  = 1d0 + 2d0*(lambdas/lambda)**(1d0-alpha)*SIN(alpha*PIh) &
           &         +        (lambdas/lambda)**(2d0-2d0*alpha)

      epsr   = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)    &
           &                           * SIN(alpha*PIh)+1d0)) / nenner
      epsi   =          ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)    &
           &                           * COS(alpha*PIh)+0d0)) / nenner  &
           & + lambda*1.25664d0/1.88496d0
      m_complex_water_ray_vec(k) = SQRT(CMPLX(epsr,-epsi,kind=dp))
    END DO

    RETURN
  END FUNCTION m_complex_water_ray_vec

  FUNCTION m_complex_water_ray(lambda,T)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Wasser als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.001,1.0] m; T aus [-10.0,30.0] Grad C
    ! nach Ray (1972)

    REAL(KIND=dp), INTENT(in) :: T, lambda

    COMPLEX(kind=dp) :: m_complex_water_ray

    REAL(KIND=dp)               :: epsinf,epss,alpha,lambdas,nenner,epsr,epsi,Toff

    Toff = T - 25d0
    epsinf  = 5.27137d0 + 0.02164740d0 * T  &
         &              - 0.00131198d0 * T*T
    epss    = 78.54d+0 * (1d0 - 4.579d-3 * Toff  &
         &                      + 1.190d-5 * Toff*Toff  &
         &                      - 2.800d-8 * Toff*Toff*Toff)
    alpha   = -16.8129d0/(T+T0C_fwo) + 0.0609265d0
    lambdas = 0.00033836d0 * exp(2513.98d0/(T+T0C_fwo)) * 1d-2

    nenner  = 1d0 + 2d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PIh) &
         &         +      (lambdas/lambda)**(2d0-2d0*alpha)

    epsr   = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)    &
         &                           * sin(alpha*PIh)+1d0)) / nenner
    epsi   =          ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)    &
         &                           * cos(alpha*PIh)+0d0)) / nenner  &
         & + lambda*1.25664d0/1.88496d0

    m_complex_water_ray = SQRT(CMPLX(epsr,-epsi,kind=dp))

    RETURN
  END FUNCTION m_complex_water_ray

  !----------------------------------------------------------------------------------------!
  ! 1b) Komplexer Brechungsindex für Eis nach Ray (1972)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_ice_ray_vec(lambda,T,anz)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.001,10**7] m; T aus [-20.0,0.0] Grad C
    ! nach Ray (1972)

    !  ;;; Vorsicht: Erzeugt falsche Werte im Vergleich zur Abb. in Ray
    !  ;;; 1972  -- muss an den Formeln im Paper liegen!
    !  ;;; mit 1d-3 bei lambdas scheints zu stimmen ...

    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_ice_ray_vec(anz)
    REAL(KIND=dp)     :: T(anz),lambda
    REAL(KIND=dp)     :: epsinf,epss,alpha,lambdas,sigma,nenner,epsr,epsi
    INTEGER :: k

    DO k=1,anz
      epsinf  = 3.168d0
      epss    = 203.168d0 + 2.5d0*T(k) + 0.15d0*T(k)*T(k)
      alpha   = 0.288d0 + 0.00520d0 * T(k) &
           &            + 0.00023d0 * T(k)*T(k)
      lambdas = 9.990288d-4 * EXP(13200d0/(T(k)+T0C_fwo)/1.9869d0) * 1d-3
      sigma   = 1.26d0 * EXP(-12500d0/(T(k)+T0C_fwo)/1.9869d0)

      nenner = 1d0 + 2d0*(lambdas/lambda)**(1d0-alpha)*SIN(alpha*PIh)  &
           & + (lambdas/lambda)**(2d0-2d0*alpha)

      epsr   = epsinf + ((epss-epsinf) * (1d0+(lambdas/lambda)**(1d0-alpha) &
           &                           * SIN(alpha*PIh))) / nenner
      epsi   =          ((epss-epsinf) * (lambdas/lambda)**(1d0-alpha) &
           &                           * COS(alpha*PIh)) / nenner &
           & +   lambda*sigma/1.88496d9
      m_complex_ice_ray_vec(k) = SQRT(CMPLX(epsr,-epsi,kind=dp))
    END DO


  END FUNCTION m_complex_ice_ray_vec

  FUNCTION m_complex_ice_ray(lambda,T)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.001,10**7] m; T aus [-20.0,0.0] Grad C
    ! nach Ray (1972)

    !  ;;; Vorsicht: Erzeugt falsche Werte im Vergleich zur Abb. in Ray
    !  ;;; 1972  -- muss an den Formeln im Paper liegen!
    !  ;;; mit 1d-3 bei lambdas scheints zu stimmen ...

    COMPLEX(kind=dp) :: m_complex_ice_ray
    REAL(KIND=dp)     :: T,lambda
    REAL(KIND=dp)     :: epsinf,epss,alpha,lambdas,sigma,nenner,epsr,epsi

    epsinf  = 3.168d0
    epss    = 203.168d0 + 2.5d0*T + 0.15d0*T*T
    alpha   = 0.288d0 + 0.00520d0 * T &
         &            + 0.00023d0 * T*T
    lambdas = 9.990288d-4 * exp(13200d0/(T+T0C_fwo)/1.9869d0) * 1d-3
    sigma   = 1.26d0 * exp(-12500d0/(T+T0C_fwo)/1.9869d0)

    nenner = 1d0 + 2d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PIh)  &
         & + (lambdas/lambda)**(2d0-2d0*alpha)

    epsr   = epsinf + ((epss-epsinf) * (1d0+(lambdas/lambda)**(1d0-alpha) &
         &                           * sin(alpha*PIh))) / nenner
    epsi   =          ((epss-epsinf) * (lambdas/lambda)**(1d0-alpha) &
         &                           * cos(alpha*PIh)) / nenner &
         & +   lambda*sigma/1.88496d9

    m_complex_ice_ray = SQRT(CMPLX(epsr,-epsi,kind=dp))

  END FUNCTION m_complex_ice_ray

  !----------------------------------------------------------------------------------------!
  ! 1c) Komplexer Brechungsindex für Flüssigwasser nach Liebe (1991)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_water_liebe_vec(lambda,T,anz)

    ! Komplexer Brechungsindex von Wasser als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.0003,0.3] m; T aus [-3.0,30.0] Grad C
    ! nach Liebe (1991) (aus den Matlab-Programmen von Prof. Maetzler, Juni 2003)

    IMPLICIT NONE

    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_water_liebe_vec(anz)
    REAL(KIND=dp)     :: T(anz),lambda
    REAL(KIND=dp)     :: c,TK,fGHz,TETA,e0,e1,f1,e2,f2
    INTEGER :: k

    c = 2.99d8 ! m/s, Lichtgeschw. im Vakuum
    fGHz = c / lambda * 1d-9

    DO k=1,anz
      TK = T(k) + T0C_fwo ! Kelvin
      TETA = 1d0 - 300d0/TK
      e0 = 77.66d0 - 103.3d0*TETA
      e1 = 0.0671d0*e0
      f1 = 20.2d0 + 146.4d0*TETA + 316d0*TETA*TETA
      e2 = 3.52d0 + 7.52d0*TETA
      ! version of Liebe MPM 1993 uses: e2=3.52
      f2 = 39.8d0*f1
      ! Liebe verwendet die Konvention mit positivem Imaginaerteil; auf
      ! negativen umstellen.
      ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
      ! positivem Realteil (Hauptzweig), was in diesem Falle richtig ist.
      m_complex_water_liebe_vec(k) = SQRT(CONJG(e2 + (e1-e2)/CMPLX(1d0,-fGHz/f2,kind=dp) &
           + (e0-e1)/CMPLX(1d0,-fGHz/f1,kind=dp)))
    END DO

  END FUNCTION m_complex_water_liebe_vec

  FUNCTION m_complex_water_liebe(lambda,T)

    ! Komplexer Brechungsindex von Wasser als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.0003,0.3] m; T aus [-3.0,30.0] Grad C
    ! nach Liebe (1991) (aus den Matlab-Programmen von Prof. Maetzler, Juni 2003)

    IMPLICIT NONE

    COMPLEX(kind=dp) :: m_complex_water_liebe
    REAL(KIND=dp)     :: T,lambda
    REAL(KIND=dp)     :: c,TK,fGHz,TETA,e0,e1,f1,e2,f2


    c = 2.99d8 ! m/s, Lichtgeschw. im Vakuum
    TK = T + T0C_fwo ! Kelvin
    fGHz = c / lambda * 1d-9

    TETA = 1d0 - 300d0/TK
    e0 = 77.66d0 - 103.3d0*TETA
    e1 = 0.0671d0*e0
    f1 = 20.2d0 + 146.4d0*TETA + 316d0*TETA*TETA
    e2 = 3.52d0 + 7.52d0*TETA
    ! version of Liebe MPM 1993 uses: e2=3.52
    f2 = 39.8d0*f1
    ! Liebe verwendet die Konvention mit positivem Imaginaerteil; auf
    ! negativen umstellen.
    ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
    ! positivem Realteil (Hauptzweig), was in diesem Falle richtig ist.
    m_complex_water_liebe = SQRT(CONJG(e2 + (e1-e2)/CMPLX(1d0,-fGHz/f2,kind=dp) + &
         (e0-e1)/CMPLX(1d0,-fGHz/f1,kind=dp)))

  END FUNCTION m_complex_water_liebe

  !----------------------------------------------------------------------------------------!
  ! 1d) Komplexer Brechungsindex für Eis nach Maetzler (1998)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_ice_maetzler_vec(lambda,T,anz)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.0001,30] m; T aus [-250.0,0.0] Grad C
    ! nach Maetzler (1998)

    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_ice_maetzler_vec(anz)
    REAL(KIND=dp)     :: T(anz),lambda
    REAL(KIND=dp)     :: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa
    INTEGER :: k

!!! Originalkommentar aus der Matlab-Routine von Prof. Maetzler:
!!! ============================================================
    ! Function for calculating the relative permittivity of pure ice in the
    ! microwave region, according to C. Maetzler, "Microwave properties of ice
    ! and snow", in B. Schmitt et al. (eds.) Solar System Ices, Astrophys.
    ! and Space Sci. Library, Vol. 227, Kluwer Academic Publishers,
    ! Dordrecht, pp. 241-257 (1998). Input:
    ! TK = temperature (K), range 20 to 273.15
    ! f = frequency in GHz, range 0.01 to 3000

    c = 2.99d8 ! m/s, Lichtgeschw. im Vakuum
    f = c / lambda * 1d-9

    DO k= 1,anz
      TK = T(k) + T0C_fwo ! Kelvin
      B1 = 0.0207d0
      B2 = 1.16d-11
      b = 335d0
      deltabeta = EXP(-10.02d0 + 0.0364d0*(TK-T0C_fwo))
      betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1d0)**2d0) ) + B2*f*f
      beta = betam + deltabeta
      theta = 300d0 / TK - 1d0
      alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
      m_complex_ice_maetzler_vec(k) = 3.1884d0 + 9.1d-4*(TK-T0C_fwo)
      m_complex_ice_maetzler_vec(k) = m_complex_ice_maetzler_vec(k) + &
           CMPLX(0d0, (alfa/f + beta*f), kind=dp)
      ! Maetzler verwendet die Konvention mit positivem Imaginaerteil; auf
      ! negativen umstellen.
      ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
      ! positivem Realteil (Hauptzweig), was in diesem Falle richtig ist.
      m_complex_ice_maetzler_vec(k) = SQRT(CONJG(m_complex_ice_maetzler_vec(k)))
    END DO

  END FUNCTION m_complex_ice_maetzler_vec

  FUNCTION m_complex_ice_maetzler(lambda,T)
    IMPLICIT NONE

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [0.0001,30] m; T aus [-250.0,0.0] Grad C
    ! nach Maetzler (1998)

    COMPLEX(kind=dp) :: m_complex_ice_maetzler
    REAL(KIND=dp)     :: T,lambda
    REAL(KIND=dp)     :: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa


!!! Originalkommentar aus der Matlab-Routine von Prof. Maetzler:
!!! ============================================================
    ! Function for calculating the relative permittivity of pure ice in the
    ! microwave region, according to C. Maetzler, "Microwave properties of ice
    ! and snow", in B. Schmitt et al. (eds.) Solar System Ices, Astrophys.
    ! and Space Sci. Library, Vol. 227, Kluwer Academic Publishers,
    ! Dordrecht, pp. 241-257 (1998). Input:
    ! TK = temperature (K), range 20 to 273.15
    ! f = frequency in GHz, range 0.01 to 3000

    c = 2.99d8 ! m/s, Lichtgeschw. im Vakuum
    TK = T + T0C_fwo ! Kelvin
    f = c / lambda * 1d-9

    B1 = 0.0207d0
    B2 = 1.16d-11
    b = 335d0
    deltabeta = EXP(-10.02d0 + 0.0364d0*(TK-T0C_fwo))
    betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1d0)**2d0) ) + B2*f*f
    beta = betam + deltabeta
    theta = 300d0 / TK - 1d0
    alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
    m_complex_ice_maetzler = 3.1884d0 + 9.1d-4*(TK-T0C_fwo)
    m_complex_ice_maetzler = m_complex_ice_maetzler + CMPLX(0d0, (alfa/f + beta*f), kind=dp)
    ! Maetzler verwendet die Konvention mit positivem Imaginaerteil; auf
    ! negativen umstellen.
    ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
    ! positivem Realteil (Hauptzweig), was in diesem Falle richtig ist.
    m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))

  END FUNCTION m_complex_ice_maetzler

  !----------------------------------------------------------------------------------------!
  ! 1e) Komplexer Brechungsindex für Eis nach Warren (1983)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_ice_warren_vec(lambda,T,anz)

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [45*1d-9,8.6] m; T aus [-60.0,0.0] Grad C
    ! nach Warren (1983) - seine Fortran-routine wird benutzt
    ! Brechungsindex m_i wird in der Konvention mit negativem
    ! Imaginaerteil zurueckgegeben.

    IMPLICIT NONE

    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_ice_warren_vec(anz)
    REAL(KIND=dp)     :: T(anz),lambda
    REAL(KIND=dp)     :: m_r(anz), m_c(anz), absind(anz), abscof(anz)

    CALL REFICEVEC(1,lambda*1000d0,T+T0C_fwo,m_r,m_c,absind,abscof,anz)

    m_complex_ice_warren_vec = CMPLX(m_r,-m_c,kind=dp)

  CONTAINS

    SUBROUTINE REFICEVEC(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF,anz)

      IMPLICIT NONE
      ! Arguments:
      INTEGER :: anz
      INTEGER :: IUNIT
      REAL(KIND=dp) :: ABSCOF(anz),ABSIND(anz),CN(anz),XLAM,RN(anz),T(anz)
      ! Parameters:
      INTEGER :: I,LT1(anz),LT2(anz),NWL,NWLT
      PARAMETER(NWL=468,NWLT=62)
      ! Local variables:
      REAL(KIND=dp) :: &
           ALAM,CUTICE,T1,T2,TK(anz),WLMAX,WLMIN, &
           X,X1,X2,Y,Y1,Y2,YLO,YHI

      REAL(KIND=dp) :: &
           TABIM(NWL),TABIMT(NWLT,4),TABRE(NWL),TABRET(NWLT,4),TEMREF(4),&
           WL(NWL),WLT(NWLT)

      INTEGER :: k, I1(1)
      INTEGER :: Ihlp(NWL), IhlpT(NWLT)

      !     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
      !     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
      !     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.

      !     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
      !                               RN  VS.        T
      !                           LOG(CN) VS. LOG(XLAM)
      !                           LOG(CN) VS.        T

      !     STEPHEN G. WARREN - 1983
      !     DEPT. OF ATMOSPHERIC SCIENCES
      !     UNIVERSITY OF WASHINGTON
      !     SEATTLE, WA  98195

      !     BASED ON

      !        WARREN,S.G.,1984.
      !        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
      !        APPLIED OPTICS,23,1206-1225

      !     INPUT PARAMETERS

      !     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
      !           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
      !           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
      !           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
      !     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
      !     T = TEMPERATURE ( DEGREES KELVIN )

      !     OUTPUT PARAMETERS

      !     RN = REAL PORTION ( SCATTERING )
      !     CN = COMPLEX PORTION ( ABSORPTION )
      !     ABSIND = ABSORPTIVE INDEX ( CN/RN )
      !     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )

      !      DIMENSION WL(NWL),WLT(NWLT)
      !      DIMENSION TABRE(NWL),TABRET(NWLT,4),TABIM(NWL),TABIMT(NWLT,4)
      !      DIMENSION TEMREF(4)

      !     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA

      TEMREF = (/272.16d0,268.16d0,253.16d0,213.16d0/)

      WLMIN = 0.045d0
      WLMAX = 8.6d6
      CUTICE = 167d0

      WL(1:114) = (/ &
           0.4430d-1,0.4510d-1,0.4590d-1,0.4680d-1,0.4770d-1,0.4860d-1,&
           0.4960d-1,0.5060d-1,0.5170d-1,0.5280d-1,0.5390d-1,0.5510d-1,&
           0.5640d-1,0.5770d-1,0.5900d-1,0.6050d-1,0.6200d-1,0.6360d-1,&
           0.6530d-1,0.6700d-1,0.6890d-1,0.7080d-1,0.7290d-1,0.7380d-1,&
           0.7510d-1,0.7750d-1,0.8000d-1,0.8270d-1,0.8550d-1,0.8860d-1,&
           0.9180d-1,0.9300d-1,0.9540d-1,0.9920d-1,0.1033d0,0.1078d0,&
           0.1100d0,0.1127d0,0.1140d0,0.1181d0,0.1210d0,0.1240d0,&
           0.1272d0,0.1295d0,0.1305d0,0.1319d0,0.1333d0,0.1348d0,&
           0.1362d0,0.1370d0,0.1378d0,0.1387d0,0.1393d0,0.1409d0,&
           0.1425d0,0.1435d0,0.1442d0,0.1450d0,0.1459d0,0.1468d0,&
           0.1476d0,0.1480d0,0.1485d0,0.1494d0,0.1512d0,0.1531d0,&
           0.1540d0,0.1550d0,0.1569d0,0.1580d0,0.1589d0,0.1610d0,&
           0.1625d0,0.1648d0,0.1669d0,0.1692d0,0.1713d0,0.1737d0,&
           0.1757d0,0.1779d0,0.1802d0,0.1809d0,0.1821d0,0.1833d0,&
           0.1843d0,0.1850d0,0.1860d0,0.1870d0,0.1880d0,0.1890d0,&
           0.1900d0,0.1910d0,0.1930d0,0.1950d0,0.2100d0,0.2500d0,&
           0.3000d0,0.3500d0,0.4000d0,0.4100d0,0.4200d0,0.4300d0,&
           0.4400d0,0.4500d0,0.4600d0,0.4700d0,0.4800d0,0.4900d0,&
           0.5000d0,0.5100d0,0.5200d0,0.5300d0,0.5400d0,0.5500d0/)
      WL(115:228) = (/ &
           0.5600d0,0.5700d0,0.5800d0,0.5900d0,0.6000d0,0.6100d0,&
           0.6200d0,0.6300d0,0.6400d0,0.6500d0,0.6600d0,0.6700d0,&
           0.6800d0,0.6900d0,0.7000d0,0.7100d0,0.7200d0,0.7300d0,&
           0.7400d0,0.7500d0,0.7600d0,0.7700d0,0.7800d0,0.7900d0,&
           0.8000d0,0.8100d0,0.8200d0,0.8300d0,0.8400d0,0.8500d0,&
           0.8600d0,0.8700d0,0.8800d0,0.8900d0,0.9000d0,0.9100d0,&
           0.9200d0,0.9300d0,0.9400d0,0.9500d0,0.9600d0,0.9700d0,&
           0.9800d0,0.9900d0,0.1000d1,0.1010d1,0.1020d1,0.1030d1,&
           0.1040d1,0.1050d1,0.1060d1,0.1070d1,0.1080d1,0.1090d1,&
           0.1100d1,0.1110d1,0.1120d1,0.1130d1,0.1140d1,0.1150d1,&
           0.1160d1,0.1170d1,0.1180d1,0.1190d1,0.1200d1,0.1210d1,&
           0.1220d1,0.1230d1,0.1240d1,0.1250d1,0.1260d1,0.1270d1,&
           0.1280d1,0.1290d1,0.1300d1,0.1310d1,0.1320d1,0.1330d1,&
           0.1340d1,0.1350d1,0.1360d1,0.1370d1,0.1380d1,0.1390d1,&
           0.1400d1,0.1410d1,0.1420d1,0.1430d1,0.1440d1,0.1449d1,&
           0.1460d1,0.1471d1,0.1481d1,0.1493d1,0.1504d1,0.1515d1,&
           0.1527d1,0.1538d1,0.1563d1,0.1587d1,0.1613d1,0.1650d1,&
           0.1680d1,0.1700d1,0.1730d1,0.1760d1,0.1800d1,0.1830d1,&
           0.1840d1,0.1850d1,0.1855d1,0.1860d1,0.1870d1,0.1890d1/)
      WL(229:342) = (/ &
           0.1905d1,0.1923d1,0.1942d1,0.1961d1,0.1980d1,0.2000d1,&
           0.2020d1,0.2041d1,0.2062d1,0.2083d1,0.2105d1,0.2130d1,&
           0.2150d1,0.2170d1,0.2190d1,0.2220d1,0.2240d1,0.2245d1,&
           0.2250d1,0.2260d1,0.2270d1,0.2290d1,0.2310d1,0.2330d1,&
           0.2350d1,0.2370d1,0.2390d1,0.2410d1,0.2430d1,0.2460d1,&
           0.2500d1,0.2520d1,0.2550d1,0.2565d1,0.2580d1,0.2590d1,&
           0.2600d1,0.2620d1,0.2675d1,0.2725d1,0.2778d1,0.2817d1,&
           0.2833d1,0.2849d1,0.2865d1,0.2882d1,0.2899d1,0.2915d1,&
           0.2933d1,0.2950d1,0.2967d1,0.2985d1,0.3003d1,0.3021d1,&
           0.3040d1,0.3058d1,0.3077d1,0.3096d1,0.3115d1,0.3135d1,&
           0.3155d1,0.3175d1,0.3195d1,0.3215d1,0.3236d1,0.3257d1,&
           0.3279d1,0.3300d1,0.3322d1,0.3345d1,0.3367d1,0.3390d1,&
           0.3413d1,0.3436d1,0.3460d1,0.3484d1,0.3509d1,0.3534d1,&
           0.3559d1,0.3624d1,0.3732d1,0.3775d1,0.3847d1,0.3969d1,&
           0.4099d1,0.4239d1,0.4348d1,0.4387d1,0.4444d1,0.4505d1,&
           0.4547d1,0.4560d1,0.4580d1,0.4719d1,0.4904d1,0.5000d1,&
           0.5100d1,0.5200d1,0.5263d1,0.5400d1,0.5556d1,0.5714d1,&
           0.5747d1,0.5780d1,0.5814d1,0.5848d1,0.5882d1,0.6061d1,&
           0.6135d1,0.6250d1,0.6289d1,0.6329d1,0.6369d1,0.6410d1/)
      WL(343:456) = (/ &
           0.6452d1,0.6494d1,0.6579d1,0.6667d1,0.6757d1,0.6897d1,&
           0.7042d1,0.7143d1,0.7246d1,0.7353d1,0.7463d1,0.7576d1,&
           0.7692d1,0.7812d1,0.7937d1,0.8065d1,0.8197d1,0.8333d1,&
           0.8475d1,0.8696d1,0.8929d1,0.9091d1,0.9259d1,0.9524d1,&
           0.9804d1,0.1000d2,0.1020d2,0.1031d2,0.1042d2,0.1053d2,&
           0.1064d2,0.1075d2,0.1087d2,0.1100d2,0.1111d2,0.1136d2,&
           0.1163d2,0.1190d2,0.1220d2,0.1250d2,0.1282d2,0.1299d2,&
           0.1316d2,0.1333d2,0.1351d2,0.1370d2,0.1389d2,0.1408d2,&
           0.1429d2,0.1471d2,0.1515d2,0.1538d2,0.1563d2,0.1613d2,&
           0.1639d2,0.1667d2,0.1695d2,0.1724d2,0.1818d2,0.1887d2,&
           0.1923d2,0.1961d2,0.2000d2,0.2041d2,0.2083d2,0.2222d2,&
           0.2260d2,0.2305d2,0.2360d2,0.2460d2,0.2500d2,0.2600d2,&
           0.2857d2,0.3100d2,0.3333d2,0.3448d2,0.3564d2,0.3700d2,&
           0.3824d2,0.3960d2,0.4114d2,0.4276d2,0.4358d2,0.4458d2,&
           0.4550d2,0.4615d2,0.4671d2,0.4736d2,0.4800d2,0.4878d2,&
           0.5003d2,0.5128d2,0.5275d2,0.5350d2,0.5424d2,0.5500d2,&
           0.5574d2,0.5640d2,0.5700d2,0.5746d2,0.5840d2,0.5929d2,&
           0.6000d2,0.6100d2,0.6125d2,0.6250d2,0.6378d2,0.6467d2,&
           0.6558d2,0.6655d2,0.6760d2,0.6900d2,0.7053d2,0.7300d2/)
      WL(457:468) = (/ &
           0.7500d2,0.7629d2,0.8000d2,0.8297d2,0.8500d2,0.8680d2,&
           0.9080d2,0.9517d2,0.1000d3,0.1200d3,0.1500d3,0.1670d3/)
      WLT = (/ &
           0.1670d3,0.1778d3,0.1884d3,&
           0.1995d3,0.2113d3,0.2239d3,0.2371d3,0.2512d3,0.2661d3,&
           0.2818d3,0.2985d3,0.3162d3,0.3548d3,0.3981d3,0.4467d3,&
           0.5012d3,0.5623d3,0.6310d3,0.7943d3,0.1000d4,0.1259d4,&
           0.2500d4,0.5000d4,0.1000d5,0.2000d5,0.3200d5,0.3500d5,&
           0.4000d5,0.4500d5,0.5000d5,0.6000d5,0.7000d5,0.9000d5,&
           0.1110d6,0.1200d6,0.1300d6,0.1400d6,0.1500d6,0.1600d6,&
           0.1700d6,0.1800d6,0.2000d6,0.2500d6,0.2900d6,0.3200d6,&
           0.3500d6,0.3800d6,0.4000d6,0.4500d6,0.5000d6,0.6000d6,&
           0.6400d6,0.6800d6,0.7200d6,0.7600d6,0.8000d6,0.8400d6,&
           0.9000d6,0.1000d7,0.2000d7,0.5000d7,0.8600d7/)
      TABRE(1:114) = (/ &
           0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,&
           0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,&
           0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,&
           0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,&
           1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,&
           1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,&
           1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,&
           1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,&
           1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,&
           1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,&
           1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,&
           1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,&
           1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,&
           1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,&
           1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,&
           1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,&
           1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,&
           1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,&
           1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/)
      TABRE(115:228) = (/ &
           1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,&
           1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,&
           1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,&
           1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,&
           1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,&
           1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,&
           1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,&
           1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,&
           1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,&
           1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,&
           1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,&
           1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,&
           1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,&
           1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,&
           1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,&
           1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,&
           1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,&
           1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,&
           1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/)
      TABRE(229:342) = (/ &
           1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,&
           1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,&
           1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,&
           1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,&
           1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,&
           1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,&
           1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,&
           1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,&
           0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,&
           1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,&
           1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,&
           1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,&
           1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,&
           1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,&
           1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,&
           1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,&
           1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,&
           1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,&
           1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/)
      TABRE(343:456) = (/ &
           1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,&
           1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,&
           1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,&
           1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,&
           1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,&
           1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,&
           1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,&
           1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,&
           1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,&
           1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,&
           1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,&
           1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,&
           1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,&
           1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,&
           1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,&
           1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,&
           1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,&
           1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,&
           1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/)
      TABRE(457:468) = (/ &
           1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,&
           1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/)
      TABRET(1:NWLT,1) = (/ &
           1.82961,   1.83258,   1.83149,&
           1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,&
           1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,&
           1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,&
           1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,&
           1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78800/)
      TABRET(1:NWLT,2) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,&
           1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,&
           1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,&
           1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,&
           1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78720/)
      TABRET(1:NWLT,3) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,   1.82224,&
           1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,&
           1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,&
           1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,&
           1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,&
           1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78400,   1.78450/)
      TABRET(1:NWLT,4) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,&
           1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,&
           1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,&
           1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77800/)
      TABIM(1:114)  = (/ &
           0.1640d0,0.1730d0,0.1830d0,0.1950d0,0.2080d0,0.2230d0,&
           0.2400d0,0.2500d0,0.2590d0,0.2680d0,0.2790d0,0.2970d0,&
           0.3190d0,0.3400d0,0.3660d0,0.3920d0,0.4160d0,0.4400d0,&
           0.4640d0,0.4920d0,0.5170d0,0.5280d0,0.5330d0,0.5340d0,&
           0.5310d0,0.5240d0,0.5100d0,0.5000d0,0.4990d0,0.4680d0,&
           0.3800d0,0.3600d0,0.3390d0,0.3180d0,0.2910d0,0.2510d0,&
           0.2440d0,0.2390d0,0.2390d0,0.2440d0,0.2470d0,0.2240d0,&
           0.1950d0,0.1740d0,0.1720d0,0.1800d0,0.1940d0,0.2130d0,&
           0.2430d0,0.2710d0,0.2890d0,0.3340d0,0.3440d0,0.3820d0,&
           0.4010d0,0.4065d0,0.4050d0,0.3890d0,0.3770d0,0.3450d0,&
           0.3320d0,0.3150d0,0.2980d0,0.2740d0,0.2280d0,0.1980d0,&
           0.1720d0,0.1560d0,0.1100d0,0.8300d-1,0.5800d-1,0.2200d-1,&
           0.1000d-1,0.3000d-2,0.1000d-2,0.3000d-3,0.1000d-3,0.3000d-4,&
           0.1000d-4,0.3000d-5,0.1000d-5,0.7000d-6,0.4000d-6,0.2000d-6,&
           0.1000d-6,0.6377d-7,0.3750d-7,0.2800d-7,0.2400d-7,0.2200d-7,&
           0.1900d-7,0.1750d-7,0.1640d-7,0.1590d-7,0.1325d-7,0.8623d-8,&
           0.5504d-8,0.3765d-8,0.2710d-8,0.2510d-8,0.2260d-8,0.2080d-8,&
           0.1910d-8,0.1540d-8,0.1530d-8,0.1550d-8,0.1640d-8,0.1780d-8,&
           0.1910d-8,0.2140d-8,0.2260d-8,0.2540d-8,0.2930d-8,0.3110d-8/)
      TABIM(115:228)  = (/ &
           0.3290d-8,0.3520d-8,0.4040d-8,0.4880d-8,0.5730d-8,0.6890d-8,&
           0.8580d-8,0.1040d-7,0.1220d-7,0.1430d-7,0.1660d-7,0.1890d-7,&
           0.2090d-7,0.2400d-7,0.2900d-7,0.3440d-7,0.4030d-7,0.4300d-7,&
           0.4920d-7,0.5870d-7,0.7080d-7,0.8580d-7,0.1020d-6,0.1180d-6,&
           0.1340d-6,0.1400d-6,0.1430d-6,0.1450d-6,0.1510d-6,0.1830d-6,&
           0.2150d-6,0.2650d-6,0.3350d-6,0.3920d-6,0.4200d-6,0.4440d-6,&
           0.4740d-6,0.5110d-6,0.5530d-6,0.6020d-6,0.7550d-6,0.9260d-6,&
           0.1120d-5,0.1330d-5,0.1620d-5,0.2000d-5,0.2250d-5,0.2330d-5,&
           0.2330d-5,0.2170d-5,0.1960d-5,0.1810d-5,0.1740d-5,0.1730d-5,&
           0.1700d-5,0.1760d-5,0.1820d-5,0.2040d-5,0.2250d-5,0.2290d-5,&
           0.3040d-5,0.3840d-5,0.4770d-5,0.5760d-5,0.6710d-5,0.8660d-5,&
           0.1020d-4,0.1130d-4,0.1220d-4,0.1290d-4,0.1320d-4,0.1350d-4,&
           0.1330d-4,0.1320d-4,0.1320d-4,0.1310d-4,0.1320d-4,0.1320d-4,&
           0.1340d-4,0.1390d-4,0.1420d-4,0.1480d-4,0.1580d-4,0.1740d-4,&
           0.1980d-4,0.2500d-4,0.5400d-4,0.1040d-3,0.2030d-3,0.2708d-3,&
           0.3511d-3,0.4299d-3,0.5181d-3,0.5855d-3,0.5899d-3,0.5635d-3,&
           0.5480d-3,0.5266d-3,0.4394d-3,0.3701d-3,0.3372d-3,0.2410d-3,&
           0.1890d-3,0.1660d-3,0.1450d-3,0.1280d-3,0.1030d-3,0.8600d-4,&
           0.8220d-4,0.8030d-4,0.8500d-4,0.9900d-4,0.1500d-3,0.2950d-3/)
      TABIM(229:342)  = (/ &
           0.4687d-3,0.7615d-3,0.1010d-2,0.1313d-2,0.1539d-2,0.1588d-2,&
           0.1540d-2,0.1412d-2,0.1244d-2,0.1068d-2,0.8414d-3,0.5650d-3,&
           0.4320d-3,0.3500d-3,0.2870d-3,0.2210d-3,0.2030d-3,0.2010d-3,&
           0.2030d-3,0.2140d-3,0.2320d-3,0.2890d-3,0.3810d-3,0.4620d-3,&
           0.5480d-3,0.6180d-3,0.6800d-3,0.7300d-3,0.7820d-3,0.8480d-3,&
           0.9250d-3,0.9200d-3,0.8920d-3,0.8700d-3,0.8900d-3,0.9300d-3,&
           0.1010d-2,0.1350d-2,0.3420d-2,0.7920d-2,0.2000d-1,0.3800d-1,&
           0.5200d-1,0.6800d-1,0.9230d-1,0.1270d0,0.1690d0,0.2210d0,&
           0.2760d0,0.3120d0,0.3470d0,0.3880d0,0.4380d0,0.4930d0,&
           0.5540d0,0.6120d0,0.6250d0,0.5930d0,0.5390d0,0.4910d0,&
           0.4380d0,0.3720d0,0.3000d0,0.2380d0,0.1930d0,0.1580d0,&
           0.1210d0,0.1030d0,0.8360d-1,0.6680d-1,0.5400d-1,0.4220d-1,&
           0.3420d-1,0.2740d-1,0.2200d-1,0.1860d-1,0.1520d-1,0.1260d-1,&
           0.1060d-1,0.8020d-2,0.6850d-2,0.6600d-2,0.6960d-2,0.9160d-2,&
           0.1110d-1,0.1450d-1,0.2000d-1,0.2300d-1,0.2600d-1,0.2900d-1,&
           0.2930d-1,0.3000d-1,0.2850d-1,0.1730d-1,0.1290d-1,0.1200d-1,&
           0.1250d-1,0.1340d-1,0.1400d-1,0.1750d-1,0.2400d-1,0.3500d-1,&
           0.3800d-1,0.4200d-1,0.4600d-1,0.5200d-1,0.5700d-1,0.6900d-1,&
           0.7000d-1,0.6700d-1,0.6500d-1,0.6400d-1,0.6200d-1,0.5900d-1/)
      TABIM(343:456)  = (/ &
           0.5700d-1,0.5600d-1,0.5500d-1,0.5700d-1,0.5800d-1,0.5700d-1,&
           0.5500d-1,0.5500d-1,0.5400d-1,0.5200d-1,0.5200d-1,0.5200d-1,&
           0.5200d-1,0.5000d-1,0.4700d-1,0.4300d-1,0.3900d-1,0.3700d-1,&
           0.3900d-1,0.4000d-1,0.4200d-1,0.4400d-1,0.4500d-1,0.4600d-1,&
           0.4700d-1,0.5100d-1,0.6500d-1,0.7500d-1,0.8800d-1,0.1080d0,&
           0.1340d0,0.1680d0,0.2040d0,0.2480d0,0.2800d0,0.3410d0,&
           0.3790d0,0.4090d0,0.4220d0,0.4220d0,0.4030d0,0.3890d0,&
           0.3740d0,0.3540d0,0.3350d0,0.3150d0,0.2940d0,0.2710d0,&
           0.2460d0,0.1980d0,0.1640d0,0.1520d0,0.1420d0,0.1280d0,&
           0.1250d0,0.1230d0,0.1160d0,0.1070d0,0.7900d-1,0.7200d-1,&
           0.7600d-1,0.7500d-1,0.6700d-1,0.5500d-1,0.4500d-1,0.2900d-1,&
           0.2750d-1,0.2700d-1,0.2730d-1,0.2890d-1,0.3000d-1,0.3400d-1,&
           0.5300d-1,0.7550d-1,0.1060d0,0.1350d0,0.1761d0,0.2229d0,&
           0.2746d0,0.3280d0,0.3906d0,0.4642d0,0.5247d0,0.5731d0,&
           0.6362d0,0.6839d0,0.7091d0,0.6790d0,0.6250d0,0.5654d0,&
           0.5433d0,0.5292d0,0.5070d0,0.4883d0,0.4707d0,0.4203d0,&
           0.3771d0,0.3376d0,0.3056d0,0.2835d0,0.3170d0,0.3517d0,&
           0.3902d0,0.4509d0,0.4671d0,0.4779d0,0.4890d0,0.4899d0,&
           0.4873d0,0.4766d0,0.4508d0,0.4193d0,0.3880d0,0.3433d0/)
      TABIM(457:468)  = (/ &
           0.3118d0,0.2935d0,0.2350d0,0.1981d0,0.1865d0,0.1771d0,&
           0.1620d0,0.1490d0,0.1390d0,0.1200d0,0.9620d-1,0.8300d-1/)


      TABIMT(1:NWLT,1)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,&
           0.4560d-1,0.3790d-1,0.3140d-1,0.2620d-1,0.2240d-1,0.1960d-1,&
           0.1760d-1,0.1665d-1,0.1620d-1,0.1550d-1,0.1470d-1,0.1390d-1,&
           0.1320d-1,0.1250d-1,0.1180d-1,0.1060d-1,0.9540d-2,0.8560d-2,&
           0.6210d-2,0.4490d-2,0.3240d-2,0.2340d-2,0.1880d-2,0.1740d-2,&
           0.1500d-2,0.1320d-2,0.1160d-2,0.8800d-3,0.6950d-3,0.4640d-3,&
           0.3400d-3,0.3110d-3,0.2940d-3,0.2790d-3,0.2700d-3,0.2640d-3,&
           0.2580d-3,0.2520d-3,0.2490d-3,0.2540d-3,0.2640d-3,0.2740d-3,&
           0.2890d-3,0.3050d-3,0.3150d-3,0.3460d-3,0.3820d-3,0.4620d-3,&
           0.5000d-3,0.5500d-3,0.5950d-3,0.6470d-3,0.6920d-3,0.7420d-3,&
           0.8200d-3,0.9700d-3,0.1950d-2,0.5780d-2,0.9700d-2/)
      TABIMT(1:NWLT,2)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4560d-1,&
           0.3790d-1,0.3140d-1,0.2620d-1,0.2240d-1,0.1960d-1,0.1760d-1,&
           0.1665d-1,0.1600d-1,0.1500d-1,0.1400d-1,0.1310d-1,0.1230d-1,&
           0.1150d-1,0.1080d-1,0.9460d-2,0.8290d-2,0.7270d-2,0.4910d-2,&
           0.3300d-2,0.2220d-2,0.1490d-2,0.1140d-2,0.1060d-2,0.9480d-3,&
           0.8500d-3,0.7660d-3,0.6300d-3,0.5200d-3,0.3840d-3,0.2960d-3,&
           0.2700d-3,0.2520d-3,0.2440d-3,0.2360d-3,0.2300d-3,0.2280d-3,&
           0.2250d-3,0.2200d-3,0.2160d-3,0.2170d-3,0.2200d-3,0.2250d-3,&
           0.2320d-3,0.2390d-3,0.2600d-3,0.2860d-3,0.3560d-3,0.3830d-3,&
           0.4150d-3,0.4450d-3,0.4760d-3,0.5080d-3,0.5400d-3,0.5860d-3,&
           0.6780d-3,0.1280d-2,0.3550d-2,0.5600d-2/)
      TABIMT(1:NWLT,3)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4560d-1,0.3790d-1,&
           0.3140d-1,0.2620d-1,0.2190d-1,0.1880d-1,0.1660d-1,0.1540d-1,&
           0.1470d-1,0.1350d-1,0.1250d-1,0.1150d-1,0.1060d-1,0.9770d-2,&
           0.9010d-2,0.7660d-2,0.6520d-2,0.5540d-2,0.3420d-2,0.2100d-2,&
           0.1290d-2,0.7930d-3,0.5700d-3,0.5350d-3,0.4820d-3,0.4380d-3,&
           0.4080d-3,0.3500d-3,0.3200d-3,0.2550d-3,0.2120d-3,0.2000d-3,&
           0.1860d-3,0.1750d-3,0.1660d-3,0.1560d-3,0.1490d-3,0.1440d-3,&
           0.1350d-3,0.1210d-3,0.1160d-3,0.1160d-3,0.1170d-3,0.1200d-3,&
           0.1230d-3,0.1320d-3,0.1440d-3,0.1680d-3,0.1800d-3,0.1900d-3,&
           0.2090d-3,0.2160d-3,0.2290d-3,0.2400d-3,0.2600d-3,0.2920d-3,&
           0.6100d-3,0.1020d-2,0.1810d-2/)
      TABIMT(1:NWLT,4)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4450d-1,0.3550d-1,0.2910d-1,&
           0.2440d-1,0.1970d-1,0.1670d-1,0.1400d-1,0.1235d-1,0.1080d-1,&
           0.8900d-2,0.7340d-2,0.6400d-2,0.5600d-2,0.5000d-2,0.4520d-2,&
           0.3680d-2,0.2990d-2,0.2490d-2,0.1550d-2,0.9610d-3,0.5950d-3,&
           0.3690d-3,0.2670d-3,0.2510d-3,0.2290d-3,0.2110d-3,0.1960d-3,&
           0.1730d-3,0.1550d-3,0.1310d-3,0.1130d-3,0.1060d-3,0.9900d-4,&
           0.9300d-4,0.8730d-4,0.8300d-4,0.7870d-4,0.7500d-4,0.6830d-4,&
           0.5600d-4,0.4960d-4,0.4550d-4,0.4210d-4,0.3910d-4,0.3760d-4,&
           0.3400d-4,0.3100d-4,0.2640d-4,0.2510d-4,0.2430d-4,0.2390d-4,&
           0.2370d-4,0.2380d-4,0.2400d-4,0.2460d-4,0.2660d-4,0.4450d-4,&
           0.8700d-4,0.1320d-3/)

      !     ZERO PARAMETERS

      RN=0d0
      CN=0d0
      ABSIND=0d0
      ABSCOF=0d0

      !     CONVERT WAVELENGTH TO MICRONS

      ALAM=XLAM
      IF (IUNIT.EQ.1) ALAM=1000d0*ALAM
      IF (IUNIT.EQ.2) ALAM=10000d0*ALAM
      IF (IUNIT.EQ.3) ALAM=10000d0*(1d0/ALAM)

      IF (ALAM.GE.WLMIN.AND.ALAM.LE.WLMAX) THEN
        IF (ALAM.LE.CUTICE) THEN
          !     REGION FROM 0.045 MICRONS TO 167.0 MICRONS - NO TEMPERATURE DEPEND
          !        DO I=2,NWL
          !          IF (ALAM.LT.WL(I)) EXIT
          !        END DO
          ! alternatively: vectorized search for interpolation interval:
          Ihlp = 0
          DO k=2,NWL
            IF (ALAM.LE.WL(k).AND.ALAM.GE.WL(k-1)) THEN
              Ihlp(k) = 1
            END IF
          END DO
          I1 = MAXLOC(Ihlp)
          I = I1(1)

          X1=LOG(WL(I-1))
          X2=LOG(WL(I))
          Y1=TABRE(I-1)
          Y2=TABRE(I)
          X=LOG(ALAM)
          Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
          RN=Y
          Y1=LOG(ABS(TABIM(I-1)))
          Y2=LOG(ABS(TABIM(I)))
          Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
          CN=EXP(Y)
        ELSE
          !     REGION FROM 167.0 MICRONS TO 8.6 METERS - TEMPERATURE DEPENDENCE
          TK=T
          WHERE (TK.GT.TEMREF(1))
            TK=TEMREF(1)
          END WHERE
          WHERE (TK.LT.TEMREF(4))
            TK=TEMREF(4)
          END WHERE

          LT1 = 2
          DO I=2,4
            DO k=1,anz
              IF(TK(k).GE.TEMREF(I).AND.TK(k).LE.TEMREF(I-1)) THEN
                LT1(k) = I
              END IF
            END DO
          END DO
          LT2=LT1-1

          !        DO I=2,NWLT
          !          IF (ALAM.LT.WLT(I)) EXIT
          !        END DO
          ! alternatively: vectorized search for interpolation interval:
          IhlpT = 0
          DO k=2,NWLT
            IF (ALAM.LE.WLT(k).AND.ALAM.GE.WLT(k-1)) THEN
              IhlpT(k) = 1
            END IF
          END DO
          I1 = MAXLOC(IhlpT)
          I = I1(1)

          X1=LOG(WLT(I-1))
          X2=LOG(WLT(I))

          X=LOG(ALAM)
          DO k=1,anz
            T1=TEMREF(LT1(k))
            T2=TEMREF(LT2(k))
            Y1=TABRET(I-1,LT1(k))
            Y2=TABRET(I,LT1(k))
            YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
            Y1=TABRET(I-1,LT2(k))
            Y2=TABRET(I,LT2(k))
            YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
            Y=((TK(k)-T1)*(YHI-YLO)/(T2-T1))+YLO
            RN(k)=Y
            Y1=LOG(ABS(TABIMT(I-1,LT1(k))))
            Y2=LOG(ABS(TABIMT(I,LT1(k))))
            YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
            Y1=LOG(ABS(TABIMT(I-1,LT2(k))))
            Y2=LOG(ABS(TABIMT(I,LT2(k))))
            YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
            Y=((TK(k)-T1)*(YHI-YLO)/(T2-T1))+YLO
            CN(k)=EXP(Y)
          END DO
        END IF

        !     ABSORPTIVE QUANITIES
        ABSIND=CN/RN
        ABSCOF=4.0*PI_dp*CN/ALAM

      END IF

      RETURN

    END SUBROUTINE REFICEVEC

  END FUNCTION m_complex_ice_warren_vec

  FUNCTION m_complex_ice_warren(lambda,T)

    ! Komplexer Brechungsindex von Eis als Funktion der Temperatur
    ! [Grad C] und der Wellenlaenge [m]; Gueltig fuer
    ! lambda aus [45*1d-9,8.6] m; T aus [-60.0,0.0] Grad C
    ! nach Warren (1983) - seine Fortran-routine wird benutzt
    ! Brechungsindex m_i wird in der Konvention mit negativem
    ! Imaginaerteil zurueckgegeben.

    IMPLICIT NONE

    COMPLEX(kind=dp) :: m_complex_ice_warren
    REAL(KIND=dp)     :: T,lambda,TT
    REAL(KIND=dp)     :: m_r, m_c, absind, abscof

    TT = T + T0C_fwo
    CALL REFICE(1,lambda*1000d0,TT,m_r,m_c,absind,abscof)

    m_complex_ice_warren = CMPLX(m_r,-m_c,kind=dp)

  CONTAINS

    SUBROUTINE REFICE(IUNIT,XLAM,T,RN,CN,ABSIND,ABSCOF)

      IMPLICIT NONE
      ! Arguments:
      INTEGER :: IUNIT
      REAL(KIND=dp) :: ABSCOF,ABSIND,CN,XLAM,RN,T
      ! Parameters:
      INTEGER :: I,LT1,LT2,NWL,NWLT
      PARAMETER(NWL=468,NWLT=62)
      ! Local variables:
      REAL(KIND=dp) :: &
           ALAM,CUTICE,T1,T2,TK,WLMAX,WLMIN, &
           X,X1,X2,Y,Y1,Y2,YLO,YHI

      REAL(KIND=dp) :: &
           TABIM(NWL),TABIMT(NWLT,4),TABRE(NWL),TABRET(NWLT,4),TEMREF(4),&
           WL(NWL),WLT(NWLT)

      !     DEFINES WAVELENGTH DEPENDENT COMPLEX INDEX OF REFRACTION FOR ICE.
      !     ALLOWABLE WAVELENGTH RANGE EXTENDS FROM 0.045 MICRONS TO 8.6 METER
      !     TEMPERATURE DEPENDENCE ONLY CONSIDERED BEYOND 167 MICRONS.

      !     INTERPOLATION IS DONE     RN  VS. LOG(XLAM)
      !                               RN  VS.        T
      !                           LOG(CN) VS. LOG(XLAM)
      !                           LOG(CN) VS.        T

      !     STEPHEN G. WARREN - 1983
      !     DEPT. OF ATMOSPHERIC SCIENCES
      !     UNIVERSITY OF WASHINGTON
      !     SEATTLE, WA  98195

      !     BASED ON

      !        WARREN,S.G.,1984.
      !        OPTICAL CONSTANTS OF ICE FROM THE ULTRAVIOLET TO THE MICROWAVE.
      !        APPLIED OPTICS,23,1206-1225

      !     INPUT PARAMETERS

      !     IUNIT = 0 FOR WAVELENGTH SPECIFIED IN MICRONS
      !           = 1 FOR WAVELENGTH SPECIFIED IN MILLIMETERS
      !           = 2 FOR WAVELENGTH SPECIFIED IN CENTIMETERS
      !           = 3 FOR WAVELENGTH SPECIFIED IN INVERSE CENTIMETERS ( WAVE N
      !     XLAM = WAVELENGTH ( MICRONS OR MM OR CM OR CM**-1 )
      !     T = TEMPERATURE ( DEGREES KELVIN )

      !     OUTPUT PARAMETERS

      !     RN = REAL PORTION ( SCATTERING )
      !     CN = COMPLEX PORTION ( ABSORPTION )
      !     ABSIND = ABSORPTIVE INDEX ( CN/RN )
      !     ABSCOF = ABORPTION COEFFICIENT ( 4*PI*CN/XLAM )

      !      DIMENSION WL(NWL),WLT(NWLT)
      !      DIMENSION TABRE(NWL),TABRET(NWLT,4),TABIM(NWL),TABIMT(NWLT,4)
      !      DIMENSION TEMREF(4)

      !     REFERENCE TEMPERATURES ARE -1.0,-5.0,-20.0, AND -60.0 DEG CENTIGRA

      TEMREF = (/272.16d0,268.16d0,253.16d0,213.16d0/)

      WLMIN = 0.045d0
      WLMAX = 8.6D6
      CUTICE = 167d0

      WL(1:114) = (/ &
           0.4430d-1,0.4510d-1,0.4590d-1,0.4680d-1,0.4770d-1,0.4860d-1,&
           0.4960d-1,0.5060d-1,0.5170d-1,0.5280d-1,0.5390d-1,0.5510d-1,&
           0.5640d-1,0.5770d-1,0.5900d-1,0.6050d-1,0.6200d-1,0.6360d-1,&
           0.6530d-1,0.6700d-1,0.6890d-1,0.7080d-1,0.7290d-1,0.7380d-1,&
           0.7510d-1,0.7750d-1,0.8000d-1,0.8270d-1,0.8550d-1,0.8860d-1,&
           0.9180d-1,0.9300d-1,0.9540d-1,0.9920d-1,0.1033d0,0.1078d0,&
           0.1100d0,0.1127d0,0.1140d0,0.1181d0,0.1210d0,0.1240d0,&
           0.1272d0,0.1295d0,0.1305d0,0.1319d0,0.1333d0,0.1348d0,&
           0.1362d0,0.1370d0,0.1378d0,0.1387d0,0.1393d0,0.1409d0,&
           0.1425d0,0.1435d0,0.1442d0,0.1450d0,0.1459d0,0.1468d0,&
           0.1476d0,0.1480d0,0.1485d0,0.1494d0,0.1512d0,0.1531d0,&
           0.1540d0,0.1550d0,0.1569d0,0.1580d0,0.1589d0,0.1610d0,&
           0.1625d0,0.1648d0,0.1669d0,0.1692d0,0.1713d0,0.1737d0,&
           0.1757d0,0.1779d0,0.1802d0,0.1809d0,0.1821d0,0.1833d0,&
           0.1843d0,0.1850d0,0.1860d0,0.1870d0,0.1880d0,0.1890d0,&
           0.1900d0,0.1910d0,0.1930d0,0.1950d0,0.2100d0,0.2500d0,&
           0.3000d0,0.3500d0,0.4000d0,0.4100d0,0.4200d0,0.4300d0,&
           0.4400d0,0.4500d0,0.4600d0,0.4700d0,0.4800d0,0.4900d0,&
           0.5000d0,0.5100d0,0.5200d0,0.5300d0,0.5400d0,0.5500d0/)
      WL(115:228) = (/ &
           0.5600d0,0.5700d0,0.5800d0,0.5900d0,0.6000d0,0.6100d0,&
           0.6200d0,0.6300d0,0.6400d0,0.6500d0,0.6600d0,0.6700d0,&
           0.6800d0,0.6900d0,0.7000d0,0.7100d0,0.7200d0,0.7300d0,&
           0.7400d0,0.7500d0,0.7600d0,0.7700d0,0.7800d0,0.7900d0,&
           0.8000d0,0.8100d0,0.8200d0,0.8300d0,0.8400d0,0.8500d0,&
           0.8600d0,0.8700d0,0.8800d0,0.8900d0,0.9000d0,0.9100d0,&
           0.9200d0,0.9300d0,0.9400d0,0.9500d0,0.9600d0,0.9700d0,&
           0.9800d0,0.9900d0,0.1000d1,0.1010d1,0.1020d1,0.1030d1,&
           0.1040d1,0.1050d1,0.1060d1,0.1070d1,0.1080d1,0.1090d1,&
           0.1100d1,0.1110d1,0.1120d1,0.1130d1,0.1140d1,0.1150d1,&
           0.1160d1,0.1170d1,0.1180d1,0.1190d1,0.1200d1,0.1210d1,&
           0.1220d1,0.1230d1,0.1240d1,0.1250d1,0.1260d1,0.1270d1,&
           0.1280d1,0.1290d1,0.1300d1,0.1310d1,0.1320d1,0.1330d1,&
           0.1340d1,0.1350d1,0.1360d1,0.1370d1,0.1380d1,0.1390d1,&
           0.1400d1,0.1410d1,0.1420d1,0.1430d1,0.1440d1,0.1449d1,&
           0.1460d1,0.1471d1,0.1481d1,0.1493d1,0.1504d1,0.1515d1,&
           0.1527d1,0.1538d1,0.1563d1,0.1587d1,0.1613d1,0.1650d1,&
           0.1680d1,0.1700d1,0.1730d1,0.1760d1,0.1800d1,0.1830d1,&
           0.1840d1,0.1850d1,0.1855d1,0.1860d1,0.1870d1,0.1890d1/)
      WL(229:342) = (/ &
           0.1905d1,0.1923d1,0.1942d1,0.1961d1,0.1980d1,0.2000d1,&
           0.2020d1,0.2041d1,0.2062d1,0.2083d1,0.2105d1,0.2130d1,&
           0.2150d1,0.2170d1,0.2190d1,0.2220d1,0.2240d1,0.2245d1,&
           0.2250d1,0.2260d1,0.2270d1,0.2290d1,0.2310d1,0.2330d1,&
           0.2350d1,0.2370d1,0.2390d1,0.2410d1,0.2430d1,0.2460d1,&
           0.2500d1,0.2520d1,0.2550d1,0.2565d1,0.2580d1,0.2590d1,&
           0.2600d1,0.2620d1,0.2675d1,0.2725d1,0.2778d1,0.2817d1,&
           0.2833d1,0.2849d1,0.2865d1,0.2882d1,0.2899d1,0.2915d1,&
           0.2933d1,0.2950d1,0.2967d1,0.2985d1,0.3003d1,0.3021d1,&
           0.3040d1,0.3058d1,0.3077d1,0.3096d1,0.3115d1,0.3135d1,&
           0.3155d1,0.3175d1,0.3195d1,0.3215d1,0.3236d1,0.3257d1,&
           0.3279d1,0.3300d1,0.3322d1,0.3345d1,0.3367d1,0.3390d1,&
           0.3413d1,0.3436d1,0.3460d1,0.3484d1,0.3509d1,0.3534d1,&
           0.3559d1,0.3624d1,0.3732d1,0.3775d1,0.3847d1,0.3969d1,&
           0.4099d1,0.4239d1,0.4348d1,0.4387d1,0.4444d1,0.4505d1,&
           0.4547d1,0.4560d1,0.4580d1,0.4719d1,0.4904d1,0.5000d1,&
           0.5100d1,0.5200d1,0.5263d1,0.5400d1,0.5556d1,0.5714d1,&
           0.5747d1,0.5780d1,0.5814d1,0.5848d1,0.5882d1,0.6061d1,&
           0.6135d1,0.6250d1,0.6289d1,0.6329d1,0.6369d1,0.6410d1/)
      WL(343:456) = (/ &
           0.6452d1,0.6494d1,0.6579d1,0.6667d1,0.6757d1,0.6897d1,&
           0.7042d1,0.7143d1,0.7246d1,0.7353d1,0.7463d1,0.7576d1,&
           0.7692d1,0.7812d1,0.7937d1,0.8065d1,0.8197d1,0.8333d1,&
           0.8475d1,0.8696d1,0.8929d1,0.9091d1,0.9259d1,0.9524d1,&
           0.9804d1,0.1000d2,0.1020d2,0.1031d2,0.1042d2,0.1053d2,&
           0.1064d2,0.1075d2,0.1087d2,0.1100d2,0.1111d2,0.1136d2,&
           0.1163d2,0.1190d2,0.1220d2,0.1250d2,0.1282d2,0.1299d2,&
           0.1316d2,0.1333d2,0.1351d2,0.1370d2,0.1389d2,0.1408d2,&
           0.1429d2,0.1471d2,0.1515d2,0.1538d2,0.1563d2,0.1613d2,&
           0.1639d2,0.1667d2,0.1695d2,0.1724d2,0.1818d2,0.1887d2,&
           0.1923d2,0.1961d2,0.2000d2,0.2041d2,0.2083d2,0.2222d2,&
           0.2260d2,0.2305d2,0.2360d2,0.2460d2,0.2500d2,0.2600d2,&
           0.2857d2,0.3100d2,0.3333d2,0.3448d2,0.3564d2,0.3700d2,&
           0.3824d2,0.3960d2,0.4114d2,0.4276d2,0.4358d2,0.4458d2,&
           0.4550d2,0.4615d2,0.4671d2,0.4736d2,0.4800d2,0.4878d2,&
           0.5003d2,0.5128d2,0.5275d2,0.5350d2,0.5424d2,0.5500d2,&
           0.5574d2,0.5640d2,0.5700d2,0.5746d2,0.5840d2,0.5929d2,&
           0.6000d2,0.6100d2,0.6125d2,0.6250d2,0.6378d2,0.6467d2,&
           0.6558d2,0.6655d2,0.6760d2,0.6900d2,0.7053d2,0.7300d2/)
      WL(457:468) = (/ &
           0.7500d2,0.7629d2,0.8000d2,0.8297d2,0.8500d2,0.8680d2,&
           0.9080d2,0.9517d2,0.1000d3,0.1200d3,0.1500d3,0.1670d3/)
      WLT = (/ &
           0.1670d3,0.1778d3,0.1884d3,&
           0.1995d3,0.2113d3,0.2239d3,0.2371d3,0.2512d3,0.2661d3,&
           0.2818d3,0.2985d3,0.3162d3,0.3548d3,0.3981d3,0.4467d3,&
           0.5012d3,0.5623d3,0.6310d3,0.7943d3,0.1000d4,0.1259d4,&
           0.2500d4,0.5000d4,0.1000d5,0.2000d5,0.3200d5,0.3500d5,&
           0.4000d5,0.4500d5,0.5000d5,0.6000d5,0.7000d5,0.9000d5,&
           0.1110d6,0.1200d6,0.1300d6,0.1400d6,0.1500d6,0.1600d6,&
           0.1700d6,0.1800d6,0.2000d6,0.2500d6,0.2900d6,0.3200d6,&
           0.3500d6,0.3800d6,0.4000d6,0.4500d6,0.5000d6,0.6000d6,&
           0.6400d6,0.6800d6,0.7200d6,0.7600d6,0.8000d6,0.8400d6,&
           0.9000d6,0.1000d7,0.2000d7,0.5000d7,0.8600d7/)
      TABRE(1:114) = (/ &
           0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,&
           0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,&
           0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,&
           0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,&
           1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,&
           1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,&
           1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,&
           1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,&
           1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,&
           1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,&
           1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,&
           1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,&
           1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,&
           1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,&
           1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,&
           1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,&
           1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,&
           1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,&
           1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/)
      TABRE(115:228) = (/ &
           1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,&
           1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,&
           1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,&
           1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,&
           1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,&
           1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,&
           1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,&
           1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,&
           1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,&
           1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,&
           1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,&
           1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,&
           1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,&
           1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,&
           1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,&
           1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,&
           1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,&
           1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,&
           1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/)
      TABRE(229:342) = (/ &
           1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,&
           1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,&
           1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,&
           1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,&
           1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,&
           1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,&
           1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,&
           1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,&
           0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,&
           1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,&
           1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,&
           1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,&
           1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,&
           1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,&
           1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,&
           1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,&
           1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,&
           1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,&
           1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/)
      TABRE(343:456) = (/ &
           1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,&
           1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,&
           1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,&
           1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,&
           1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,&
           1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,&
           1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,&
           1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,&
           1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,&
           1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,&
           1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,&
           1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,&
           1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,&
           1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,&
           1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,&
           1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,&
           1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,&
           1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,&
           1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/)
      TABRE(457:468) = (/ &
           1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,&
           1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/)
      TABRET(1:NWLT,1) = (/ &
           1.82961,   1.83258,   1.83149,&
           1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,&
           1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,&
           1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,&
           1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,&
           1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,&
           1.78720,   1.78720,   1.78720,   1.78720,   1.78800/)
      TABRET(1:NWLT,2) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,&
           1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,&
           1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,&
           1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,&
           1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,&
           1.78650,   1.78650,   1.78650,   1.78720/)
      TABRET(1:NWLT,3) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,   1.82224,&
           1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,&
           1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,&
           1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,&
           1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,&
           1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,&
           1.78370,   1.78400,   1.78450/)
      TABRET(1:NWLT,4) = (/ &
           1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,&
           1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,&
           1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,&
           1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,&
           1.77720,   1.77800/)
      TABIM(1:114)  = (/ &
           0.1640d0,0.1730d0,0.1830d0,0.1950d0,0.2080d0,0.2230d0,&
           0.2400d0,0.2500d0,0.2590d0,0.2680d0,0.2790d0,0.2970d0,&
           0.3190d0,0.3400d0,0.3660d0,0.3920d0,0.4160d0,0.4400d0,&
           0.4640d0,0.4920d0,0.5170d0,0.5280d0,0.5330d0,0.5340d0,&
           0.5310d0,0.5240d0,0.5100d0,0.5000d0,0.4990d0,0.4680d0,&
           0.3800d0,0.3600d0,0.3390d0,0.3180d0,0.2910d0,0.2510d0,&
           0.2440d0,0.2390d0,0.2390d0,0.2440d0,0.2470d0,0.2240d0,&
           0.1950d0,0.1740d0,0.1720d0,0.1800d0,0.1940d0,0.2130d0,&
           0.2430d0,0.2710d0,0.2890d0,0.3340d0,0.3440d0,0.3820d0,&
           0.4010d0,0.4065d0,0.4050d0,0.3890d0,0.3770d0,0.3450d0,&
           0.3320d0,0.3150d0,0.2980d0,0.2740d0,0.2280d0,0.1980d0,&
           0.1720d0,0.1560d0,0.1100d0,0.8300d-1,0.5800d-1,0.2200d-1,&
           0.1000d-1,0.3000d-2,0.1000d-2,0.3000d-3,0.1000d-3,0.3000d-4,&
           0.1000d-4,0.3000d-5,0.1000d-5,0.7000d-6,0.4000d-6,0.2000d-6,&
           0.1000d-6,0.6377d-7,0.3750d-7,0.2800d-7,0.2400d-7,0.2200d-7,&
           0.1900d-7,0.1750d-7,0.1640d-7,0.1590d-7,0.1325d-7,0.8623d-8,&
           0.5504d-8,0.3765d-8,0.2710d-8,0.2510d-8,0.2260d-8,0.2080d-8,&
           0.1910d-8,0.1540d-8,0.1530d-8,0.1550d-8,0.1640d-8,0.1780d-8,&
           0.1910d-8,0.2140d-8,0.2260d-8,0.2540d-8,0.2930d-8,0.3110d-8/)
      TABIM(115:228)  = (/ &
           0.3290d-8,0.3520d-8,0.4040d-8,0.4880d-8,0.5730d-8,0.6890d-8,&
           0.8580d-8,0.1040d-7,0.1220d-7,0.1430d-7,0.1660d-7,0.1890d-7,&
           0.2090d-7,0.2400d-7,0.2900d-7,0.3440d-7,0.4030d-7,0.4300d-7,&
           0.4920d-7,0.5870d-7,0.7080d-7,0.8580d-7,0.1020d-6,0.1180d-6,&
           0.1340d-6,0.1400d-6,0.1430d-6,0.1450d-6,0.1510d-6,0.1830d-6,&
           0.2150d-6,0.2650d-6,0.3350d-6,0.3920d-6,0.4200d-6,0.4440d-6,&
           0.4740d-6,0.5110d-6,0.5530d-6,0.6020d-6,0.7550d-6,0.9260d-6,&
           0.1120d-5,0.1330d-5,0.1620d-5,0.2000d-5,0.2250d-5,0.2330d-5,&
           0.2330d-5,0.2170d-5,0.1960d-5,0.1810d-5,0.1740d-5,0.1730d-5,&
           0.1700d-5,0.1760d-5,0.1820d-5,0.2040d-5,0.2250d-5,0.2290d-5,&
           0.3040d-5,0.3840d-5,0.4770d-5,0.5760d-5,0.6710d-5,0.8660d-5,&
           0.1020d-4,0.1130d-4,0.1220d-4,0.1290d-4,0.1320d-4,0.1350d-4,&
           0.1330d-4,0.1320d-4,0.1320d-4,0.1310d-4,0.1320d-4,0.1320d-4,&
           0.1340d-4,0.1390d-4,0.1420d-4,0.1480d-4,0.1580d-4,0.1740d-4,&
           0.1980d-4,0.2500d-4,0.5400d-4,0.1040d-3,0.2030d-3,0.2708d-3,&
           0.3511d-3,0.4299d-3,0.5181d-3,0.5855d-3,0.5899d-3,0.5635d-3,&
           0.5480d-3,0.5266d-3,0.4394d-3,0.3701d-3,0.3372d-3,0.2410d-3,&
           0.1890d-3,0.1660d-3,0.1450d-3,0.1280d-3,0.1030d-3,0.8600d-4,&
           0.8220d-4,0.8030d-4,0.8500d-4,0.9900d-4,0.1500d-3,0.2950d-3/)
      TABIM(229:342)  = (/ &
           0.4687d-3,0.7615d-3,0.1010d-2,0.1313d-2,0.1539d-2,0.1588d-2,&
           0.1540d-2,0.1412d-2,0.1244d-2,0.1068d-2,0.8414d-3,0.5650d-3,&
           0.4320d-3,0.3500d-3,0.2870d-3,0.2210d-3,0.2030d-3,0.2010d-3,&
           0.2030d-3,0.2140d-3,0.2320d-3,0.2890d-3,0.3810d-3,0.4620d-3,&
           0.5480d-3,0.6180d-3,0.6800d-3,0.7300d-3,0.7820d-3,0.8480d-3,&
           0.9250d-3,0.9200d-3,0.8920d-3,0.8700d-3,0.8900d-3,0.9300d-3,&
           0.1010d-2,0.1350d-2,0.3420d-2,0.7920d-2,0.2000d-1,0.3800d-1,&
           0.5200d-1,0.6800d-1,0.9230d-1,0.1270d0,0.1690d0,0.2210d0,&
           0.2760d0,0.3120d0,0.3470d0,0.3880d0,0.4380d0,0.4930d0,&
           0.5540d0,0.6120d0,0.6250d0,0.5930d0,0.5390d0,0.4910d0,&
           0.4380d0,0.3720d0,0.3000d0,0.2380d0,0.1930d0,0.1580d0,&
           0.1210d0,0.1030d0,0.8360d-1,0.6680d-1,0.5400d-1,0.4220d-1,&
           0.3420d-1,0.2740d-1,0.2200d-1,0.1860d-1,0.1520d-1,0.1260d-1,&
           0.1060d-1,0.8020d-2,0.6850d-2,0.6600d-2,0.6960d-2,0.9160d-2,&
           0.1110d-1,0.1450d-1,0.2000d-1,0.2300d-1,0.2600d-1,0.2900d-1,&
           0.2930d-1,0.3000d-1,0.2850d-1,0.1730d-1,0.1290d-1,0.1200d-1,&
           0.1250d-1,0.1340d-1,0.1400d-1,0.1750d-1,0.2400d-1,0.3500d-1,&
           0.3800d-1,0.4200d-1,0.4600d-1,0.5200d-1,0.5700d-1,0.6900d-1,&
           0.7000d-1,0.6700d-1,0.6500d-1,0.6400d-1,0.6200d-1,0.5900d-1/)
      TABIM(343:456)  = (/ &
           0.5700d-1,0.5600d-1,0.5500d-1,0.5700d-1,0.5800d-1,0.5700d-1,&
           0.5500d-1,0.5500d-1,0.5400d-1,0.5200d-1,0.5200d-1,0.5200d-1,&
           0.5200d-1,0.5000d-1,0.4700d-1,0.4300d-1,0.3900d-1,0.3700d-1,&
           0.3900d-1,0.4000d-1,0.4200d-1,0.4400d-1,0.4500d-1,0.4600d-1,&
           0.4700d-1,0.5100d-1,0.6500d-1,0.7500d-1,0.8800d-1,0.1080d0,&
           0.1340d0,0.1680d0,0.2040d0,0.2480d0,0.2800d0,0.3410d0,&
           0.3790d0,0.4090d0,0.4220d0,0.4220d0,0.4030d0,0.3890d0,&
           0.3740d0,0.3540d0,0.3350d0,0.3150d0,0.2940d0,0.2710d0,&
           0.2460d0,0.1980d0,0.1640d0,0.1520d0,0.1420d0,0.1280d0,&
           0.1250d0,0.1230d0,0.1160d0,0.1070d0,0.7900d-1,0.7200d-1,&
           0.7600d-1,0.7500d-1,0.6700d-1,0.5500d-1,0.4500d-1,0.2900d-1,&
           0.2750d-1,0.2700d-1,0.2730d-1,0.2890d-1,0.3000d-1,0.3400d-1,&
           0.5300d-1,0.7550d-1,0.1060d0,0.1350d0,0.1761d0,0.2229d0,&
           0.2746d0,0.3280d0,0.3906d0,0.4642d0,0.5247d0,0.5731d0,&
           0.6362d0,0.6839d0,0.7091d0,0.6790d0,0.6250d0,0.5654d0,&
           0.5433d0,0.5292d0,0.5070d0,0.4883d0,0.4707d0,0.4203d0,&
           0.3771d0,0.3376d0,0.3056d0,0.2835d0,0.3170d0,0.3517d0,&
           0.3902d0,0.4509d0,0.4671d0,0.4779d0,0.4890d0,0.4899d0,&
           0.4873d0,0.4766d0,0.4508d0,0.4193d0,0.3880d0,0.3433d0/)
      TABIM(457:468)  = (/ &
           0.3118d0,0.2935d0,0.2350d0,0.1981d0,0.1865d0,0.1771d0,&
           0.1620d0,0.1490d0,0.1390d0,0.1200d0,0.9620d-1,0.8300d-1/)


      TABIMT(1:NWLT,1)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,&
           0.4560d-1,0.3790d-1,0.3140d-1,0.2620d-1,0.2240d-1,0.1960d-1,&
           0.1760d-1,0.1665d-1,0.1620d-1,0.1550d-1,0.1470d-1,0.1390d-1,&
           0.1320d-1,0.1250d-1,0.1180d-1,0.1060d-1,0.9540d-2,0.8560d-2,&
           0.6210d-2,0.4490d-2,0.3240d-2,0.2340d-2,0.1880d-2,0.1740d-2,&
           0.1500d-2,0.1320d-2,0.1160d-2,0.8800d-3,0.6950d-3,0.4640d-3,&
           0.3400d-3,0.3110d-3,0.2940d-3,0.2790d-3,0.2700d-3,0.2640d-3,&
           0.2580d-3,0.2520d-3,0.2490d-3,0.2540d-3,0.2640d-3,0.2740d-3,&
           0.2890d-3,0.3050d-3,0.3150d-3,0.3460d-3,0.3820d-3,0.4620d-3,&
           0.5000d-3,0.5500d-3,0.5950d-3,0.6470d-3,0.6920d-3,0.7420d-3,&
           0.8200d-3,0.9700d-3,0.1950d-2,0.5780d-2,0.9700d-2/)
      TABIMT(1:NWLT,2)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4560d-1,&
           0.3790d-1,0.3140d-1,0.2620d-1,0.2240d-1,0.1960d-1,0.1760d-1,&
           0.1665d-1,0.1600d-1,0.1500d-1,0.1400d-1,0.1310d-1,0.1230d-1,&
           0.1150d-1,0.1080d-1,0.9460d-2,0.8290d-2,0.7270d-2,0.4910d-2,&
           0.3300d-2,0.2220d-2,0.1490d-2,0.1140d-2,0.1060d-2,0.9480d-3,&
           0.8500d-3,0.7660d-3,0.6300d-3,0.5200d-3,0.3840d-3,0.2960d-3,&
           0.2700d-3,0.2520d-3,0.2440d-3,0.2360d-3,0.2300d-3,0.2280d-3,&
           0.2250d-3,0.2200d-3,0.2160d-3,0.2170d-3,0.2200d-3,0.2250d-3,&
           0.2320d-3,0.2390d-3,0.2600d-3,0.2860d-3,0.3560d-3,0.3830d-3,&
           0.4150d-3,0.4450d-3,0.4760d-3,0.5080d-3,0.5400d-3,0.5860d-3,&
           0.6780d-3,0.1280d-2,0.3550d-2,0.5600d-2/)
      TABIMT(1:NWLT,3)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4560d-1,0.3790d-1,&
           0.3140d-1,0.2620d-1,0.2190d-1,0.1880d-1,0.1660d-1,0.1540d-1,&
           0.1470d-1,0.1350d-1,0.1250d-1,0.1150d-1,0.1060d-1,0.9770d-2,&
           0.9010d-2,0.7660d-2,0.6520d-2,0.5540d-2,0.3420d-2,0.2100d-2,&
           0.1290d-2,0.7930d-3,0.5700d-3,0.5350d-3,0.4820d-3,0.4380d-3,&
           0.4080d-3,0.3500d-3,0.3200d-3,0.2550d-3,0.2120d-3,0.2000d-3,&
           0.1860d-3,0.1750d-3,0.1660d-3,0.1560d-3,0.1490d-3,0.1440d-3,&
           0.1350d-3,0.1210d-3,0.1160d-3,0.1160d-3,0.1170d-3,0.1200d-3,&
           0.1230d-3,0.1320d-3,0.1440d-3,0.1680d-3,0.1800d-3,0.1900d-3,&
           0.2090d-3,0.2160d-3,0.2290d-3,0.2400d-3,0.2600d-3,0.2920d-3,&
           0.6100d-3,0.1020d-2,0.1810d-2/)
      TABIMT(1:NWLT,4)  = (/ &
           0.8300d-1,0.6900d-1,0.5700d-1,0.4450d-1,0.3550d-1,0.2910d-1,&
           0.2440d-1,0.1970d-1,0.1670d-1,0.1400d-1,0.1235d-1,0.1080d-1,&
           0.8900d-2,0.7340d-2,0.6400d-2,0.5600d-2,0.5000d-2,0.4520d-2,&
           0.3680d-2,0.2990d-2,0.2490d-2,0.1550d-2,0.9610d-3,0.5950d-3,&
           0.3690d-3,0.2670d-3,0.2510d-3,0.2290d-3,0.2110d-3,0.1960d-3,&
           0.1730d-3,0.1550d-3,0.1310d-3,0.1130d-3,0.1060d-3,0.9900d-4,&
           0.9300d-4,0.8730d-4,0.8300d-4,0.7870d-4,0.7500d-4,0.6830d-4,&
           0.5600d-4,0.4960d-4,0.4550d-4,0.4210d-4,0.3910d-4,0.3760d-4,&
           0.3400d-4,0.3100d-4,0.2640d-4,0.2510d-4,0.2430d-4,0.2390d-4,&
           0.2370d-4,0.2380d-4,0.2400d-4,0.2460d-4,0.2660d-4,0.4450d-4,&
           0.8700d-4,0.1320d-3/)

      !     ZERO PARAMETERS

      RN=0d0
      CN=0d0
      ABSIND=0d0
      ABSCOF=0d0

      !     CONVERT WAVELENGTH TO MICRONS

      ALAM=XLAM
      IF (IUNIT.EQ.1) ALAM=1000d0*ALAM
      IF (IUNIT.EQ.2) ALAM=10000d0*ALAM
      IF (IUNIT.EQ.3) ALAM=10000d0*(1d0/ALAM)
      IF (ALAM.LT.WLMIN.OR.ALAM.GT.WLMAX) RETURN
      IF (ALAM.LE.CUTICE) THEN
        !     REGION FROM 0.045 MICRONS TO 167.0 MICRONS - NO TEMPERATURE DEPEND
        DO I=2,NWL
          IF (ALAM.LT.WL(I)) EXIT
        END DO
        X1=LOG(WL(I-1))
        X2=LOG(WL(I))
        Y1=TABRE(I-1)
        Y2=TABRE(I)
        X=LOG(ALAM)
        Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        RN=Y
        Y1=LOG(ABS(TABIM(I-1)))
        Y2=LOG(ABS(TABIM(I)))
        Y=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        CN=EXP(Y)
      ELSE
        !     REGION FROM 167.0 MICRONS TO 8.6 METERS - TEMPERATURE DEPENDENCE

        TK=T
        IF (TK.GT.TEMREF(1)) TK=TEMREF(1)
        IF (TK.LT.TEMREF(4)) TK=TEMREF(4)
        DO I=2,4
          IF(TK.GE.TEMREF(I)) EXIT
        END DO
        LT1=I
        LT2=I-1
        DO I=2,NWLT
          IF (ALAM.LE.WLT(I)) EXIT
        END DO
        X1=LOG(WLT(I-1))
        X2=LOG(WLT(I))
        Y1=TABRET(I-1,LT1)
        Y2=TABRET(I,LT1)
        X=LOG(ALAM)
        YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        Y1=TABRET(I-1,LT2)
        Y2=TABRET(I,LT2)
        YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        T1=TEMREF(LT1)
        T2=TEMREF(LT2)
        Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
        RN=Y
        Y1=LOG(ABS(TABIMT(I-1,LT1)))
        Y2=LOG(ABS(TABIMT(I,LT1)))
        YLO=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        Y1=LOG(ABS(TABIMT(I-1,LT2)))
        Y2=LOG(ABS(TABIMT(I,LT2)))
        YHI=((X-X1)*(Y2-Y1)/(X2-X1))+Y1
        Y=((TK-T1)*(YHI-YLO)/(T2-T1))+YLO
        CN=EXP(Y)
      END IF

      !     ABSORPTIVE QUANITIES
      ABSIND=CN/RN
      ABSCOF=4.0*PI_dp*CN/ALAM
      RETURN

    END SUBROUTINE REFICE

  END FUNCTION m_complex_ice_warren


!*********************************************************************************************

  !----------------------------------------------------------------------------------------!
  ! 2) Routinen zum effektiven Brechungsindex von Mischmaterialien nach verschiedenen
  !    Ansaetzen (es werden jeweils immer 3 Mischungspartner beruecksichtigt).
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! 2a) Homogene Mischung 3er Substanzen nach Oguchi (1983),
  !     erwaehnt in Wiener (1912), zurueckgehend auf Mosotti, Clausius, Poisson
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_oguchi_vec(vol1, vol2, vol3, m1, m2, m3, anz, fehler)

    ! vol1, vol2, vol3 sind die Volumenanteile der drei Substanzen mit dem
    ! komplexen Brechungsindex m1, m2, m3 an der homogenen Mischung (am Gesamtvolumen).
    ! Hier wird vol1 bzw. m1 als Matrix und die anderen beiden als
    ! kugelfoermige mikroskopische Einschluesse (formfakt = 2.0) behandelt.

    ! Formel ist symetrisch in der Wahl, welcher Stoff die Matrix und
    ! welche beiden die Einschluesse sind.

    IMPLICIT NONE
    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_oguchi_vec(anz)
    COMPLEX(kind=dp) :: m1(anz), m2(anz), m3(anz)
    REAL(KIND=dp) :: vol1(anz), vol2(anz), vol3(anz)
    COMPLEX(kind=dp) :: m1t(anz), m2t(anz), m3t(anz)
    REAL(KIND=dp) :: formfakt
    INTEGER, INTENT(out) :: fehler

    fehler = 0

    IF (ANY(ABS(vol1+vol2+vol3-1d0) > 1d-6)) THEN
      WRITE (*,*) 'M_COMPLEX_OGUCHI_VEC: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_oguchi_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    formfakt = 2d0  ! Einschluesse kugelfoermig

    m1t = m1*m1
    m2t = m2*m2
    m3t = m3*m3

    m_complex_oguchi_vec = SQRT((m1t*vol1+(m1t+formfakt)/(m2t+formfakt)*m2t*vol2+(m1t+formfakt)/&
         (m3t+formfakt)*m3t*vol3) / &
         (vol1+(m1t+formfakt)/(m2t+formfakt)*vol2+(m1t+formfakt)/(m3t+formfakt)*vol3))
    ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
    ! positivem Realteil, was in diesem Falle richtig ist. Das
    ! Vorzeichen des Imagin"arteils bleibt dann in der Konvention der
    ! Brechungsindices m1, m2 und m3.
    RETURN

  END FUNCTION m_complex_oguchi_vec

  FUNCTION m_complex_oguchi(vol1, vol2, vol3, m1, m2, m3, fehler)

    ! vol1, vol2, vol3 sind die Volumenanteile der drei Substanzen mit dem
    ! komplexen Brechungsindex m1, m2, m3 an der homogenen Mischung (am Gesamtvolumen).
    ! Hier wird vol1 bzw. m1 als Matrix und die anderen beiden als
    ! kugelfoermige mikroskopische Einschluesse (formfakt = 2.0) behandelt.

    ! Formel ist symetrisch in der Wahl, welcher Stoff die Matrix und
    ! welche beiden die Einschluesse sind.

    IMPLICIT NONE
    COMPLEX(kind=dp) :: m_complex_oguchi
    COMPLEX(kind=dp) :: m1, m2, m3
    REAL(KIND=dp) :: vol1, vol2, vol3
    COMPLEX(kind=dp) :: m1t, m2t, m3t
    REAL(KIND=dp) :: formfakt
    INTEGER, INTENT(out) :: fehler

    fehler = 0

    IF (ABS(vol1+vol2+vol3-1d0) > 1d-6) THEN
      WRITE (*,*) 'M_COMPLEX_OGUCHI: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_oguchi = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    formfakt = 2d0  ! Einschluesse kugelfoermig

    m1t = m1*m1
    m2t = m2*m2
    m3t = m3*m3

    m_complex_oguchi = SQRT((m1t*vol1+(m1t+formfakt)/(m2t+formfakt)*m2t*vol2+(m1t+formfakt)/(m3t+formfakt)*m3t*vol3) / &
         (vol1+(m1t+formfakt)/(m2t+formfakt)*vol2+(m1t+formfakt)/(m3t+formfakt)*vol3))
    ! Die komplexe Wurzelfunktion liefert diejenige Wurzel mit
    ! positivem Realteil, was in diesem Falle richtig ist. Das
    ! Vorzeichen des Imagin"arteils bleibt dann in der Konvention der
    ! Brechungsindices m1, m2 und m3.
    RETURN

  END FUNCTION m_complex_oguchi

  !----------------------------------------------------------------------------------------!
  ! 2b) Homogene Mischung 3er Substanzen nach Maxwell-Garnett (1904),
  !     entnommen Bohren und Huffman (1983)
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_maxwellgarnett_vec(vol1, vol2, vol3, m1, m2, m3, inclusionstring, anz, fehler)

    ! vol1, vol2, vol3 sind die Volumenanteile der drei Substanzen mit dem
    ! komplexen Brechungsindex m1, m2, m3 an der homogenen Mischung (am Gesamtvolumen).
    ! Hier wird vol1 bzw. m1 als Matrix und die anderen beiden als
    ! kugelfoermige bzw. ellipsoide mikroskopische Einschluesse behandelt.

    ! Der Stoff vol1 ist die Matrix, vol2 und vol3 sind die
    ! Einschluesse.

    ! Der inclusionstring legt fest, wie die inclusion zu behandeln
    ! ist: 'shperical' oder 'spheroidal'

    ! Diese Formel ist nicht symetrisch bezueglich der Wahl der
    ! Matrix. Ist jedoch die Matrix festgelegt, koennen die anderen
    ! beiden Stoffe (kugelf"ormige Einschluesse) beliebig vertauscht werden.

    ! Bohren, Huffman 1983; Bohren 1982

    ! Ergebnisse eines Profilings mit gprof: Die Auswertung des komplexen Logarithmus ist bei
    ! weitem die aufwendigste Funktion im untenstehenden Code. Fast 2/3 der gesamten fuer
    ! die Funktion benoetigten Rechenzeit geht hierfuer drauf. Danach kommt die
    ! komplexe Wurzelfunktion. Die select-Anweisung spielt nur eine sehr untergeordnete
    ! Rolle.

    IMPLICIT NONE
    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_maxwellgarnett_vec(anz)
    COMPLEX(kind=dp) :: m1(anz), m2(anz), m3(anz)
    REAL(KIND=dp) :: vol1(anz), vol2(anz), vol3(anz)
    CHARACTER(len=*) :: inclusionstring

    COMPLEX(kind=dp) :: beta2(anz), beta3(anz), m1t(anz), m2t(anz), m3t(anz)
    INTEGER, INTENT(out) :: fehler
    CHARACTER(len=20) :: csi

    fehler = 0

    IF (ANY(ABS(vol1+vol2+vol3-1d0) > 1d-6)) THEN
      WRITE (*,*) 'M_COMPLEX_MAXWELLGARNETT_VEC: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_maxwellgarnett_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    m1t = m1*m1
    m2t = m2*m2
    m3t = m3*m3

    csi(:) = ' '
    csi = inclusionstring(1:LEN_TRIM(inclusionstring))
    csi = ADJUSTL(csi)

!!$ here, "TRIM()" obstructs inlining, therefore changed to something which does inline:
!!$    SELECT CASE (TRIM(ADJUSTL(inclusionstring)))
    SELECT CASE (csi(1:10))
!!$    CASE ('spherical','SPHERICAL','Spherical')
    CASE ('spherical ','SPHERICAL ','Spherical ')
      beta2 = 3d0*m1t/(m2t+2d0*m1t)
      beta3 = 3d0*m1t/(m3t+2d0*m1t)
    CASE ('spheroidal','SPHEROIDAL','Spheroidal')
      ! Die LOG-Funktion liefert den Hauptzweig
      ! der Komplexen Funktion, wie gefordert:
      beta2 = 2d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1d0)
      beta3 = 2d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1d0)
    CASE DEFAULT
      WRITE (*,*) 'M_COMPLEX_MAXWELLGARNETT_VEC: Unbekannter inclusionstring: ', inclusionstring
      m_complex_maxwellgarnett_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END SELECT

    m_complex_maxwellgarnett_vec = SQRT(((1d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
         (1d0-vol2-vol3+vol2*beta2+vol3*beta3))

    ! Die komplexe Wurzelfunktion liefert diejenige der beiden Wurzeln mit
    ! positivem Realteil, was in diesem Falle richtig ist. Das
    ! Vorzeichen des Imagin"arteils bleibt dann in der Konvention der
    ! Brechungsindices m1, m2 und m3.
    RETURN
  END FUNCTION m_complex_maxwellgarnett_vec

  FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3, m1, m2, m3, inclusionstring, fehler)

    ! vol1, vol2, vol3 sind die Volumenanteile der drei Substanzen mit dem
    ! komplexen Brechungsindex m1, m2, m3 an der homogenen Mischung (am Gesamtvolumen).
    ! Hier wird vol1 bzw. m1 als Matrix und die anderen beiden als
    ! kugelfoermige bzw. ellipsoide mikroskopische Einschluesse behandelt.

    ! Der Stoff vol1 ist die Matrix, vol2 und vol3 sind die
    ! Einschluesse.

    ! Der inclusionstring legt fest, wie die inclusion zu behandeln
    ! ist: 'shperical' oder 'spheroidal'

    ! Diese Formel ist nicht symetrisch bezueglich der Wahl der
    ! Matrix. Ist jedoch die Matrix festgelegt, koennen die anderen
    ! beiden Stoffe (kugelf"ormige Einschluesse) beliebig vertauscht werden.

    ! Bohren, Huffman 1983; Bohren 1982

    ! Ergebnisse eines Profilings mit gprof: Die Auswertung des komplexen Logarithmus ist bei
    ! weitem die aufwendigste Funktion im untenstehenden Code. Fast 2/3 der gesamten fuer
    ! die Funktion benoetigten Rechenzeit geht hierfuer drauf. Danach kommt die
    ! komplexe Wurzelfunktion. Die select-Anweisung spielt nur eine sehr untergeordnete
    ! Rolle.

    IMPLICIT NONE

    COMPLEX(kind=dp) :: m_complex_maxwellgarnett
    COMPLEX(kind=dp) :: m1, m2, m3
    REAL(KIND=dp) :: vol1, vol2, vol3
    CHARACTER(len=*) :: inclusionstring

    COMPLEX(kind=dp) :: beta2, beta3, m1t, m2t, m3t
    INTEGER, INTENT(out) :: fehler

    fehler = 0

    IF (ABS(vol1+vol2+vol3-1d0) > 1d-6) THEN
      WRITE (*,*) 'M_COMPLEX_MAXWELLGARNETT: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_maxwellgarnett = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    m1t = m1*m1
    m2t = m2*m2
    m3t = m3*m3

    SELECT CASE (TRIM(ADJUSTL(inclusionstring)))
    CASE ('spherical','SPHERICAL','Spherical')
      beta2 = 3d0*m1t/(m2t+2d0*m1t)
      beta3 = 3d0*m1t/(m3t+2d0*m1t)
    CASE ('spheroidal','SPHEROIDAL','Spheroidal')
      ! Die LOG-Funktion liefert den Hauptzweig
      ! der Komplexen Funktion, wie gefordert:
      beta2 = 2d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1d0)
      beta3 = 2d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1d0)
    CASE DEFAULT
      WRITE (*,*) 'M_COMPLEX_MAXWELLGARNETT: Unbekannter inclusionstring: ', inclusionstring
      m_complex_maxwellgarnett = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END SELECT

    m_complex_maxwellgarnett = SQRT(((1d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
         (1d0-vol2-vol3+vol2*beta2+vol3*beta3))

    ! Die komplexe Wurzelfunktion liefert diejenige der beiden Wurzeln mit
    ! positivem Realteil, was in diesem Falle richtig ist. Das
    ! Vorzeichen des Imagin"arteils bleibt dann in der Konvention der
    ! Brechungsindices m1, m2 und m3.
    RETURN
  END FUNCTION m_complex_maxwellgarnett

  !----------------------------------------------------------------------------------------!
  ! 2c) Homogene Mischung 3er Substanzen nach Bruggemann
  !----------------------------------------------------------------------------------------!

  FUNCTION m_complex_bruggemann_vec (volair, volice, volwater, mair, mice, mwater, anz, fehler)

    ! volwater, volice, volair sind die Volumenanteile von Wasser, Eis
    ! und Luft an der homogenen Mischung (am Gesamtvolumen).

    IMPLICIT NONE
    INTEGER :: anz
    COMPLEX(kind=dp) :: m_complex_bruggemann_vec(anz)
    REAL(KIND=dp), INTENT(in) :: volair(anz), volice(anz), volwater(anz)
    COMPLEX(kind=dp), INTENT(in) :: mair(anz), mice(anz), mwater(anz)
    INTEGER, INTENT(out) :: fehler
    COMPLEX(kind=dp) :: m1t(anz), m2t(anz), m3t(anz), koeff(anz,4), roots(anz,3), sroots(anz,3), &
         p(anz), q(anz), u(anz), v(anz), diskri(anz), y(anz,3), epsilon1, epsilon2

    INTEGER :: i, j, count(anz)

    fehler = 0

    IF (ANY(ABS(volair+volice+volwater-1d0) > 1d-6)) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN_VEC: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_bruggemann_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    m1t = mair*mair
    m2t = mice*mice
    m3t = mwater*mwater

    ! die m's sind jetzt die epsilons!

    ! Die Bruggemann-Formel fuer N=3 auf ein Polynom 3. Grades
    ! umschreiben, von dem dann die 0-Stellen gesucht werden:
    ! hier kommen die Koeffizienten:

    koeff(:,1) = m1t*m2t*m3t
    koeff(:,2) = m1t*m2t*(2d0*volair+2d0*volice-volwater) + &
         m1t*m3t*(2d0*volair-volice+2d0*volwater) + &
         m2t*m3t*(-volair+2d0*volice+2d0*volwater)
    koeff(:,3) = m1t*(4d0*volair-2d0*volice-2d0*volwater) + &
         m2t*(-2d0*volair+4d0*volice-2d0*volwater) + &
         m3t*(-2d0*volair-2d0*volice+4d0*volwater)
    koeff(:,4) = (-4d0,0d0)

    ! Laguerre-Verfahren zur iterativen Nullstellenbestimmung
    ! (in dieser vektorisierten Fassung leider noch nicht implementiert...):
!    call zrootsvec(koeff,3,anz,roots,.true.)

    !=============================================================
    ! alternativ: Cardanische Formel (exakte Loesungen):

    ! Normalform des Polynoms:
    koeff(:,1) = koeff(:,1) / koeff(:,4)
    koeff(:,2) = koeff(:,2) / koeff(:,4)
    koeff(:,3) = koeff(:,3) / koeff(:,4)
    koeff(:,4) = koeff(:,4) / koeff(:,4)

    ! reduzierte Form:
    p = (3d0*koeff(:,2)-koeff(:,3)*koeff(:,3)) * third
    q = 2d0*koeff(:,3)*koeff(:,3)*koeff(:,3)/27d0 - koeff(:,2)*koeff(:,3)*third + koeff(:,1)

    ! Cardanische Formel fuer Polynom 3. Grades:
    diskri = p*p*p*third*third*third + q*q*0.25d0
    u = (-q/2d0+SQRT(diskri))**third
    v = -p/(3d0*u)

    epsilon1 = CMPLX(-0.5d0,SQRT(3d0)*0.5d0,kind=dp)
    epsilon2 = CMPLX(-0.5d0,-SQRT(3d0)*0.5d0,kind=dp)
    y(:,1) = u + v
    y(:,2) = epsilon1*u + epsilon2*v
    y(:,3) = epsilon2*u + epsilon1*v

    roots(:,1) = y(:,1) - koeff(:,3)*third
    roots(:,2) = y(:,2) - koeff(:,3)*third
    roots(:,3) = y(:,3) - koeff(:,3)*third
    !===============================================================

    sroots = SQRT(roots)

    ! herumspielen mit einem IDL-testdriver-Programm (darstellen der
    ! berechneten 3 sroots in der komplexen Zahlenebene) hat
    ! ergeben, dass die physikalisch sinnvolle Loesung einen Realteil
    ! hat, der zwischen 1.0 (Luft) und dem Wert fuer reines Wasser
    ! liegt. Der Imaginaerteil ist dann kleiner als 0.0 und kleiner als der Wert fuer
    ! Wasser. Mit diesen Bedingungen kann die richtige Wurzel
    ! absepariert werden, weil die anderen 2 nie in diesem Bereich liegen.

    ! Aus numerischen Gruenden nicht genau 1.0, sondern 0.98!
    count = 0
    DO i=1,3
      DO j=1,anz
        IF (REAL(sroots(j,i)) >= 0.98d0 .AND. AIMAG(sroots(j,i)) < 0d0) THEN
          count(j) = count(j) + 1
          m_complex_bruggemann_vec(j) = sroots(j,i)
        END IF
      END DO
    END DO
    IF (ANY(COUNT == 0)) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN_VEC: Keine passende Loesung gefunden! Abbruch!'
      m_complex_bruggemann_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF
    IF (ANY(COUNT > 1)) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN_VEC: Mehrere passende Loesungen gefunden! Da is was faul! Abbruch!'
      m_complex_bruggemann_vec = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    RETURN

  END FUNCTION m_complex_bruggemann_vec

  FUNCTION m_complex_bruggemann (volair, volice, volwater, mair, mice, mwater, fehler)

    ! volwater, volice, volair sind die Volumenanteile von Wasser, Eis
    ! und Luft an der homogenen Mischung (am Gesamtvolumen).

    IMPLICIT NONE
    COMPLEX(kind=dp) :: m_complex_bruggemann
    REAL(KIND=dp), INTENT(in) :: volair, volice, volwater
    COMPLEX(kind=dp), INTENT(in) :: mair, mice, mwater
    INTEGER, INTENT(out) :: fehler
    COMPLEX(kind=dp) :: m1t, m2t, m3t, koeff(4), roots(3), sroots(3), mm, check, &
         p, q, u, v, diskri, y(3), epsilon1, epsilon2

    INTEGER :: i, count

    fehler = 0

    IF (ABS(volair+volice+volwater-1d0) > 1d-6) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN: Summe der rel. Partialvolumen ergibt nicht 1! Abbruch!'
      m_complex_bruggemann = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    m1t = mair*mair
    m2t = mice*mice
    m3t = mwater*mwater

    ! die m's sind jetzt die epsilons!

    ! Die Bruggemann-Formel fuer N=3 auf ein Polynom 3. Grades
    ! umschreiben, von dem dann die 0-Stellen gesucht werden:
    ! hier kommen die Koeffizienten:

    koeff(1) = m1t*m2t*m3t
    koeff(2) = m1t*m2t*(2d0*volair+2d0*volice-volwater) + &
         m1t*m3t*(2d0*volair-volice+2d0*volwater) + &
         m2t*m3t*(-volair+2d0*volice+2d0*volwater)
    koeff(3) = m1t*(4d0*volair-2d0*volice-2d0*volwater) + &
         m2t*(-2d0*volair+4d0*volice-2d0*volwater) + &
         m3t*(-2d0*volair-2d0*volice+4d0*volwater)
    koeff(4) = (-4d0,0d0)

    ! Laguerre-Verfahren zur iterativen Nullstellenbestimmung:
!    call zroots(koeff,3,roots,.true.)

    !=============================================================
    ! alternativ: Cardanische Formel (exakte Loesungen):

    ! Normalform des Polynoms:
    koeff = koeff / koeff(4)

    ! reduzierte Form:
    p = (3d0*koeff(2)-koeff(3)*koeff(3)) * third
    q = 2d0*koeff(3)*koeff(3)*koeff(3)/27d0 - koeff(2)*koeff(3)*third + koeff(1)

    ! Cardanische Formel fuer Polynom 3. Grades:
    diskri = p*p*p*third*third*third + q*q*0.25d0
    u = (-q/2d0+SQRT(diskri))**third
    v = -p/(3d0*u)

    epsilon1 = CMPLX(-0.5d0,SQRT(3d0)*0.5d0,kind=dp)
    epsilon2 = CMPLX(-0.5d0,-SQRT(3d0)*0.5d0,kind=dp)
    y(1) = u + v
    y(2) = epsilon1*u + epsilon2*v
    y(3) = epsilon2*u + epsilon1*v

    roots = y - koeff(3)*third
    !===============================================================

    sroots = SQRT(roots)

    ! herumspielen mit einem IDL-testdriver-Programm (darstellen der
    ! berechneten 3 sroots in der komplexen Zahlenebene) hat
    ! ergeben, dass die physikalisch sinnvolle Loesung einen Realteil
    ! hat, der zwischen 1.0 (Luft) und dem Wert fuer reines Wasser
    ! liegt. Der Imaginaerteil ist dann kleiner als 0.0 und kleiner als der Wert fuer
    ! Wasser. Mit diesen Bedingungen kann die richtige Wurzel
    ! absepariert werden, weil die anderen 2 nie in diesem Bereich liegen.

    ! Aus numerischen Gruenden nicht genau 1.0, sondern 0.98!
    count = 0
    DO i=1,3
      IF (REAL(sroots(i)) >= 0.98d0 .AND. AIMAG(sroots(i)) < 0d0) THEN
        count = count + 1
        m_complex_bruggemann = sroots(i)
        mm = roots(i)
      END IF
    END DO
    IF (count == 0) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN: Keine passende Loesung gefunden! Abbruch!'
      m_complex_bruggemann = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF
    IF (count > 1) THEN
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN: Mehrere passende Loesungen gefunden! Da is was faul! Abbruch!'
      m_complex_bruggemann = CMPLX(miss_value,miss_value,kind=dp)
      fehler = 1
      RETURN
    END IF

    ! Einsetzen der Loesung in die Bruggemann-Formel muss 0 ergeben: (nur bei Anwendung des Laguerre-Verfahrens pruefen!)
!    check = volair*(m1t-mm)/(m1t+2d0*mm)+ &
!         volice*(m2t-mm)/(m2t+2d0*mm)+ &
!         volwater*(m3t-mm)/(m3t+2d0*mm)
!    IF (ABS(check) > 1d-2) THEN
!      WRITE (*,'(a,1x,d12.5)') 'M_COMPLEX_BRUGGEMANN: Die Loesung fuer m_eff hat ein Residuum von ', ABS(check)
!    END IF

    RETURN

  CONTAINS

!!! 2 Routinen fuer das Laguerre-Verfahren zur Bestimmung der Wurzeln eines komplexen Polynoms,
!!! siehe Numerical Recipes in Fortran 77, p. 365 ff.

    SUBROUTINE laguer(a,m,x,its)
      IMPLICIT NONE
      INTEGER :: m,its,MAXIT,MR,MT
      REAL(KIND=dp) :: EPSS
      COMPLEX(kind=dp) :: a(m+1),x
      PARAMETER (EPSS=2.d-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER :: iter,j
      REAL(KIND=dp) :: abx,abp,abm,err,frac(MR)
      COMPLEX(kind=dp) :: dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      frac = (/.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1d0/)
      DO iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=CMPLX(0.,0.,kind=dp)
        f=CMPLX(0.,0.,kind=dp)
        abx=abs(x)
        DO j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
        END DO
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2d0*f/b
          sq=sqrt(DBLE(m-1)*(DBLE(m)*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0d0) then
            dx=DBLE(m)/gp
          else
            dx=EXP(CMPLX(LOG(1d0+abx),dble(iter),kind=dp))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
      END DO
      WRITE (*,*) 'M_COMPLEX_BRUGGEMANN: too many iterations in laguer (Laguerre-Method)'
      STOP
      RETURN
    END SUBROUTINE laguer

    SUBROUTINE zroots(a,m,roots,polish)
      INTEGER :: m,MAXM
      REAL(KIND=dp) :: EPS
      COMPLEX(kind=dp) :: a(m+1),roots(m)
      LOGICAL :: polish
      PARAMETER (EPS=1.d-6,MAXM=101)
      !    USES laguer
      INTEGER :: i,j,jj,its
      COMPLEX(kind=dp) :: ad(MAXM),x,b,c
      DO j=1,m+1
        ad(j)=a(j)
      END DO
      DO j=m,1,-1
        x=CMPLX(0d0,0d0,kind=dp)
        CALL laguer(ad,j,x,its)
        IF(ABS(AIMAG(x)).LE.2d0*EPS*EPS*ABS(REAL(x))) x=CMPLX(DBLE(x),0d0,kind=dp)
        roots(j)=x
        b=ad(j+1)
        DO jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
        END DO
      END DO
      if (polish) then
        DO j=1,m
          call laguer(a,m,roots(j),its)
        END DO
      endif
      DO j=2,m
        x=roots(j)
        DO i=j-1,1,-1
          if(real(roots(i)).le.real(x)) exit
          roots(i+1)=roots(i)
        END DO
        !        i=0
        roots(i+1)=x
      END DO
      RETURN
    END SUBROUTINE zroots

  END FUNCTION m_complex_bruggemann

  !----------------------------------------------------------------------------------------!
  ! 2d) Interfaces zu den Brechungsindexroutinen, gesteuert durch Schluesselworstrings:
  !----------------------------------------------------------------------------------------!

  ! Fuer Eis-Wasser-Luft-Mischung in einem Durchgang:
  FUNCTION get_m_mix_vec(m_a, m_i, m_w, volair, volice, volwater, &
       mixingrulestring, matrixstring, inclusionstring, anz, fehler) RESULT(get_m_mix)

    IMPLICIT NONE

    INTEGER, INTENT(in)  :: anz
    REAL(KIND=dp), INTENT(in) :: volice(anz), volair(anz), volwater(anz)
    COMPLEX(kind=dp), INTENT(in) :: m_a(anz), m_i(anz), m_w(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring
    INTEGER, INTENT(out) :: fehler
    COMPLEX(kind=dp) :: get_m_mix(anz)

    fehler = 0
    get_m_mix = CMPLX(1d0,0d0,kind=dp)

    SELECT CASE (TRIM(ADJUSTL(mixingrulestring)))

    CASE ('maxwellgarnett','MAXWELLGARNETT','Maxwellgarnett')
      ! Als erstes Argument kommt der Volumenanteil der Matrix, dann die anderen beiden Konstituenten
      ! (Luft kommt in der Mischung nicht vor und wird mit volair = 0.0 belegt)
      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE ('ice','ICE','Ice')
!CDIR NEXPAND
        get_m_mix = m_complex_maxwellgarnett_vec(volice, volair, volwater, m_i, m_a, m_w, inclusionstring, anz, fehler)
      CASE ('water','WATER','Water')
!CDIR NEXPAND
        get_m_mix = m_complex_maxwellgarnett_vec(volwater, volair, volice, m_w, m_a, m_i, inclusionstring, anz, fehler)
      CASE ('air','AIR','Air')
!CDIR NEXPAND
        get_m_mix = m_complex_maxwellgarnett_vec(volair, volwater, volice, m_a, m_w, m_i, inclusionstring, anz, fehler)
      CASE DEFAULT
        WRITE (*,*) 'GET_M_MIX_VEC: Unbekannter matrixstring: ', matrixstring
        fehler = 1
      END SELECT
    CASE ('bruggemann','BRUGGEMANN','Bruggemann')
!CDIR NEXPAND
      get_m_mix = m_complex_bruggemann_vec(volair, volice, volwater, m_a, m_i, m_w, anz, fehler)
    CASE ('oguchi','OGUCHI','Oguchi')
!CDIR NEXPAND
      get_m_mix = m_complex_oguchi_vec(volair, volice, volwater, m_a, m_i, m_w, anz, fehler)
    CASE default
      WRITE (*,*) 'GET_M_MIX_VEC: Unbekannter mixingrulestring: ', mixingrulestring
      fehler = 2
    END SELECT
    IF (fehler /= 0) THEN
      WRITE (*,*) 'GET_M_MIX_VEC: Fehler beim eff. Brechungsindex aufgetreten!'
      STOP
      RETURN
    END IF

    RETURN
  END FUNCTION get_m_mix_vec

  COMPLEX(kind=dp) FUNCTION get_m_mix(m_a, m_i, m_w, volair, volice, volwater, &
       mixingrulestring, matrixstring, inclusionstring, fehler)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: volice, volair, volwater
    COMPLEX(kind=dp), INTENT(in) :: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring
    INTEGER, INTENT(out) :: fehler

    fehler = 0
    get_m_mix = CMPLX(1d0,0d0,kind=dp)

    SELECT CASE (TRIM(ADJUSTL(mixingrulestring)))

    CASE ('maxwellgarnett','MAXWELLGARNETT','Maxwellgarnett')
      ! Als erstes Argument kommt der Volumenanteil der Matrix, dann die anderen beiden Konstituenten
      ! (Luft kommt in der Mischung nicht vor und wird mit volair = 0.0 belegt)
      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE ('ice','ICE','Ice')
        get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater, m_i, m_a, m_w, inclusionstring, fehler)
      CASE ('water','WATER','Water')
        get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice, m_w, m_a, m_i, inclusionstring, fehler)
      CASE ('air','AIR','Air')
        get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice, m_a, m_w, m_i, inclusionstring, fehler)
      CASE DEFAULT
        WRITE (*,*) 'GET_M_MIX: Unbekannter matrixstring: ', matrixstring
        fehler = 1
      END SELECT
    CASE ('bruggemann','BRUGGEMANN','Bruggemann')
      get_m_mix = m_complex_bruggemann(volair, volice, volwater, m_a, m_i, m_w, fehler)
    CASE ('oguchi','OGUCHI','Oguchi')
      get_m_mix = m_complex_oguchi(volair, volice, volwater, m_a, m_i, m_w, fehler)
    CASE default
      WRITE (*,*) 'GET_M_MIX: Unbekannter mixingrulestring: ', mixingrulestring
      fehler = 2
    END SELECT
    IF (fehler /= 0) THEN
      WRITE (*,*) 'GET_M_MIX: Fehler beim eff. Brechungsindex aufgetreten!'
      STOP
      RETURN
    END IF

    RETURN
  END FUNCTION get_m_mix

  ! Ineinandergeschachtelt: ( (m1 + m2) + m3), wobei m1,m2,m3 jeweils Luft, Eis oder Wasser sein koennen!
  FUNCTION get_m_mix_nested_vec(m_a, m_i, m_w, volair, volice, volwater, &
       mixingrulestring, hoststring, matrixstring, inclusionstring, hostmatrixstring, &
       hostinclusionstring, anz, kumfehler) &
       RESULT (get_m_mix_nested)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: anz
    REAL(KIND=dp), INTENT(in) :: volice(anz), volair(anz), volwater(anz)
    COMPLEX(kind=dp), INTENT(in) :: m_a(anz), m_i(anz), m_w(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, hoststring, matrixstring, &
         inclusionstring, hostmatrixstring, hostinclusionstring
    INTEGER, INTENT(out) :: kumfehler
    COMPLEX(kind=dp) :: get_m_mix_nested(anz)

    REAL(KIND=dp) :: vol1(anz), vol2(anz), nullv(anz)
    COMPLEX(kind=dp) :: mtmp(anz)
    INTEGER :: fehler

    kumfehler = 0
    get_m_mix_nested = CMPLX(1d0,0d0,kind=dp)
    nullv = 0d0

    SELECT CASE (TRIM(ADJUSTL(hoststring)))

    CASE('air','AIR','Air')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('air','AIR','Air')
        WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix_vec(m_a, m_i, m_w, nullv, vol1, vol2, &
             mixingrulestring, matrixstring, inclusionstring, anz, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('air','AIR','Air')
          get_m_mix_nested = get_m_mix_vec(m_a, mtmp, 2d0*m_a, volair, (1d0-volair), nullv, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE('icewater','ICEWATER','Icewater')
          get_m_mix_nested = get_m_mix_vec(m_a, mtmp, 2d0*m_a, volair, (1d0-volair), nullv, &
               mixingrulestring, 'ice', hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE('ice','ICE','Ice')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('ice','ICE','Ice')
        WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix_vec(m_a, m_i, m_w, vol1, nullv, vol2, &
             mixingrulestring, matrixstring, inclusionstring, anz, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('ice','ICE','Ice')
          get_m_mix_nested = get_m_mix_vec(mtmp, m_i, 2d0*m_a, (1d0-volice), volice, nullv, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE('airwater','AIRWATER','Airwater')
          get_m_mix_nested = get_m_mix_vec(mtmp, m_i, 2d0*m_a, (1d0-volice), volice, nullv, &
               mixingrulestring, 'air', hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE('water','WATER','Water')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('water','WATER','Water')
        WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix_vec(m_a, m_i, m_w, vol1, vol2, nullv, &
             mixingrulestring, matrixstring, inclusionstring, anz, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('water','WATER','Water')
          get_m_mix_nested = get_m_mix_vec(2d0*m_a, mtmp, m_w, nullv, (1d0-volwater), volwater, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE('airice','AIRICE','Airice')
          get_m_mix_nested = get_m_mix_vec(2d0*m_a, mtmp, m_w, nullv, (1d0-volwater), volwater, &
               mixingrulestring, 'ice', hostinclusionstring, anz, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE ('none','NONE','None')

      get_m_mix_nested = get_m_mix_vec(m_a, m_i, m_w, volair, volice, volwater, &
           mixingrulestring, matrixstring, inclusionstring, anz, fehler)
      kumfehler = kumfehler + fehler

    CASE default
      WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Unbekannter hoststring: ', mixingrulestring
      kumfehler = kumfehler + 1
    END SELECT


    IF (kumfehler /= 0) THEN
      WRITE (*,*) 'GET_M_MIX_NESTED_VEC: Fehler beim eff. Brechungsindex aufgetreten!'
      get_m_mix_nested = CMPLX(1d0,0d0,kind=dp)
      STOP
      RETURN
    END IF

    RETURN
  END FUNCTION get_m_mix_nested_vec

  COMPLEX(kind=dp) FUNCTION get_m_mix_nested(m_a, m_i, m_w, volair, volice, volwater, &
       mixingrulestring, hoststring, matrixstring, inclusionstring, hostmatrixstring, hostinclusionstring, kumfehler)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in) :: volice, volair, volwater
    COMPLEX(kind=dp), INTENT(in) :: m_a, m_i, m_w
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, hoststring, matrixstring, &
         inclusionstring, hostmatrixstring, hostinclusionstring
    INTEGER, INTENT(out) :: kumfehler

    REAL(KIND=dp) :: vol1, vol2
    COMPLEX(kind=dp) :: mtmp
    INTEGER :: fehler

    kumfehler = 0
    get_m_mix_nested = CMPLX(1d0,0d0,kind=dp)

    SELECT CASE (TRIM(ADJUSTL(hoststring)))

    CASE('air','AIR','Air')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('air','AIR','Air')
        WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix(m_a, m_i, m_w, 0d0, vol1, vol2, &
             mixingrulestring, matrixstring, inclusionstring, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('air','AIR','Air')
          get_m_mix_nested = get_m_mix(m_a, mtmp, 2d0*m_a, volair, (1d0-volair), 0d0, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE('icewater','ICEWATER','Icewater')
          get_m_mix_nested = get_m_mix(m_a, mtmp, 2d0*m_a, volair, (1d0-volair), 0d0, &
               mixingrulestring, 'ice', hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE('ice','ICE','Ice')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('ice','ICE','Ice')
        WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix(m_a, m_i, m_w, vol1, 0d0, vol2, &
             mixingrulestring, matrixstring, inclusionstring, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('ice','ICE','Ice')
          get_m_mix_nested = get_m_mix(mtmp, m_i, 2d0*m_a, (1d0-volice), volice, 0d0, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE('airwater','AIRWATER','Airwater')
          get_m_mix_nested = get_m_mix(mtmp, m_i, 2d0*m_a, (1d0-volice), volice, 0d0, &
               mixingrulestring, 'air', hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE('water','WATER','Water')

      SELECT CASE (TRIM(ADJUSTL(matrixstring)))
      CASE('water','WATER','Water')
        WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger matrixstring: ', matrixstring
        kumfehler = kumfehler + 1
      CASE default
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1d0 - vol1
        mtmp = get_m_mix(m_a, m_i, m_w, vol1, vol2, 0d0, &
             mixingrulestring, matrixstring, inclusionstring, fehler)
        kumfehler = kumfehler + fehler

        SELECT CASE (TRIM(ADJUSTL(hostmatrixstring)))
        CASE('water','WATER','Water')
          get_m_mix_nested = get_m_mix(2d0*m_a, mtmp, m_w, 0d0, (1d0-volwater), volwater, &
               mixingrulestring, hostmatrixstring, hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE('airice','AIRICE','Airice')
          get_m_mix_nested = get_m_mix(2d0*m_a, mtmp, m_w, 0d0, (1d0-volwater), volwater, &
               mixingrulestring, 'ice', hostinclusionstring, fehler)
          kumfehler = kumfehler + fehler
        CASE default
          WRITE (*,*) 'GET_M_MIX_NESTED: Ungueltiger hostmatrixstring: ', hostmatrixstring
          kumfehler = kumfehler + 1
        END SELECT
      END SELECT

    CASE ('none','NONE','None')

      get_m_mix_nested = get_m_mix(m_a, m_i, m_w, volair, volice, volwater, &
           mixingrulestring, matrixstring, inclusionstring, fehler)
      kumfehler = kumfehler + fehler

    CASE default
      WRITE (*,*) 'GET_M_MIX_NESTED: Unbekannter hoststring: ', mixingrulestring
      kumfehler = kumfehler + 1
    END SELECT


    IF (kumfehler /= 0) THEN
      WRITE (*,*) 'GET_M_MIX_NESTED: Fehler beim eff. Brechungsindex aufgetreten!'
      get_m_mix_nested = CMPLX(1d0,0d0,kind=dp)
      STOP
      RETURN
    END IF

    RETURN
  END FUNCTION get_m_mix_nested

!*********************************************************************************************

  !----------------------------------------------------------------------------------------!
  ! 3) Routinen zur Berechnung der Streufunktionen nach Mie
  !----------------------------------------------------------------------------------------!

  !----------------------------------------------------------------------------------------!
  ! 3a) Rekursionsformeln fuer Qext, Qback und Qseit nach
  !     Bohren und Huffman (1983) fuer einschalige Kugeln
  !----------------------------------------------------------------------------------------!

  SUBROUTINE BHMIESEITVEC(X,REFREL,QBACK,QEXT,QSCA,S1,S2,anz)

    ! Es werden nur Streuwinkel von 180 Grad und 0 Grad (Extinktion!) berechnet.
    ! Erweiterungen auf bel. Streurichtung aber sehr einfach moeglich!!!

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0d0), NANG = 2

    INTEGER :: anz
    REAL(KIND=dp) :: GSCA(anz),QBACK(anz),QEXT(anz),QSCA(anz),X(anz)
    COMPLEX(kind=dp) :: REFREL(anz)
    COMPLEX(kind=dp) :: S1(anz,NANG),S2(anz,NANG)

    INTEGER :: J,N,NSTOP(anz),NMX(anz),NN,K
    REAL(KIND=dp) :: APSI(anz),APSI1(anz),CHI(anz),CHI0(anz),CHI1(anz),FN,XSTOP(anz),YMOD(anz)
    REAL(KIND=dp) ::  AMU(NANG),PI(anz,NANG),PI0(anz,NANG),PI1(anz,NANG),TAU(anz,NANG),THETA(NANG)
    REAL(KIND=dp) :: PSI0(anz),PSI1(anz),PSI(anz),DN
    COMPLEX(kind=dp) :: AN(anz),AN1(anz),BN(anz),BN1(anz),XI(anz),XI1(anz),Y(anz)
    COMPLEX(kind=dp), ALLOCATABLE :: D(:,:)
!===========================================================================
!
! Singularity for X = 0.0! The calling program has to take care of this!
!
!===========================================================================
    !***********************************************************************
    ! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
    !    to calculate scattering and absorption by a homogenous isotropic
    !    sphere.
    ! Given:
    !    X = 2*pi*a/lambda
    !    REFREL = (complex refr. index of sphere)/(real index of medium)
    !    NANG = number of angles between 0 and 90 degrees
    !           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
    !           if called with NANG<2, will set NANG=2 and will compute
    !           scattering for theta=0,90,180.
    ! Returns:
    !    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
    !                                scatt. E perp. to scatt. plane)
    !    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
    !                               scatt. E parr. to scatt. plane)
    !   QEXT = C_ext/pi*a**2 = efficiency factor for extinction
    !   QSCA = C_sca/pi*a**2 = efficiency factor for scattering
    !   QBACK = (dC_sca/domega)/pi*a**2
    !         = backscattering efficiency
    !   GSCA = <cos(theta)> for scattering
    !
    !Original program taken from Bohren and Huffman (1983), Appendix A
    !Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
    !in order to compute <cos(theta)>
    !91/05/07 (BTD): Modified to allow NANG=1
    !91/08/15 (BTD): Corrected error (failure to initialize P)
    !91/08/15 (BTD): Modified to enhance vectorizability.
    !91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
    !91/08/15 (BTD): Changed definition of QBACK.
    !92/01/08 (BTD): Note that this version has been superceded by
    !                fully double precision version = bhmie.f which,
    !                unfortunately, is not standard f77.
    !                However, retain this in case standard f77 version
    !                is required for porting to some other system.
    !2005/11/07 U. Blahak: transferred to F90 and full double precision,
    !                major modifications
    !***********************************************************************

    Y=X*REFREL
    YMOD=ABS(Y)
    !
    !*** Series expansion terminated after NSTOP terms
    !   Logarithmic derivatives calculated from NMX on down
    XSTOP=FLOOR(X+4d0*X**third+2d0)+1d0
    !*** Original code:
    !     NMX=AMAX1(XSTOP,YMOD)+15
    !     NSTOP=XSTOP
    !*** Experimental code:
    NMX=MAX(XSTOP,YMOD)+15
    NSTOP=XSTOP
    ALLOCATE(D(anz,MAXVAL(NMX)))

    !*** Require NANG.GE.1 in order to calculate scattering intensities
!!! Hier werden die Streuwinkel theta (zenitwinkel zur Vorwaertsstreurichtung)
!!! gesetzt:
    THETA(1) = 0d0 * degrad_dp
!!! THETA(...) = ... ! THETA prinzipiell erweiterbar auf beliebige Winkel. NANG entsprechend erhoehen!!!
    THETA(2) = 180d0 * degrad_dp

    DO j=1,NANG
      AMU(J)=COS(THETA(J))
    END DO
    DO J=1,NANG
      DO k=1,anz
        PI0(k,J)=0d0
        PI1(k,J)=1d0
      END DO
    END DO
    NN=NANG
    DO J=1,NANG
      DO k=1,anz
        S1(k,J)=CMPLX(0d0,0d0,kind=dp)
        S2(k,J)=CMPLX(0d0,0d0,kind=dp)
      END DO
    END DO

! Maybe simply set NMX to the scalar value maxval(NMX)?

    !*** Logarithmic derivative D(J) calculated by downward recurrence
    !   beginning with initial value (0.,0.) at J=NMX
!    DO k=1,anz
!      D(k,NMX(k))=CMPLX(0d0,0d0,kind=dp)
!    END DO
    D = CMPLX(0d0, 0d0, kind=dp)
    NN=MAXVAL(NMX)-1
    DO N=1,NN
      DO k=1,anz
        IF (NMX(k)-N >= 1) THEN
          DN=NMX(k)-N+1
          D(k,NMX(k)-N)=(DN/Y(k))-(1d0/(D(k,NMX(k)-N+1)+DN/Y(k)))
        END IF
      END DO
    END DO

    !*** Riccati-Bessel functions with real argument X
    !   calculated by upward recurrence

    PSI0=COS(X)
    PSI1=SIN(X)
    CHI0=-SIN(X)
    CHI1=COS(X)
    !APSI0 never used, so this line removed from program:
    !     APSI0=PSI0
    APSI1=PSI1
    !XI0 never used, so this line removed from program:
    !     In the next two lines, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
    !     XI0=CMPLX(APSI0,-CHI0)
    !      XI1=CMPLX(APSI1,-CHI1)
    XI1=CMPLX(APSI1,CHI1,kind=dp)
    QSCA=0d0
    GSCA=0d0
    AN = (0d0, 0d0)
    BN = (0d0, 0d0)
    DO N=1,MAXVAL(NSTOP)
      DN=N
      FN=(2d0*DN+1d0)/(DN*(DN+1d0))
      DO k=1,anz
        IF (N <= NSTOP(k)) THEN
          PSI(k)=(2d0*DN-1d0)*PSI1(k)/X(k)-PSI0(k)
          APSI(k)=PSI(k)
          CHI(k)=(2d0*DN-1d0)*CHI1(k)/X(k)-CHI0(k)

      !     In the next line, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
      !          XI(k)=CMPLX(APSI(k),-CHI(k))
          XI(k)=CMPLX(APSI(k),CHI(k),kind=dp)

      !*** Store previous values of AN and BN for use
      !   in computation of g=<cos(theta)>
          AN1(k)=AN(k)
          BN1(k)=BN(k)


      !*** Compute AN and BN:
          AN(k)=(D(k,N)/REFREL(k)+DN/X(k))*APSI(k)-APSI1(k)
          AN(k)=AN(k)/((D(k,N)/REFREL(k)+DN/X(k))*XI(k)-XI1(k))
          BN(k)=(REFREL(k)*D(k,N)+DN/X(k))*APSI(k)-APSI1(k)
          BN(k)=BN(k)/((REFREL(k)*D(k,N)+DN/X(k))*XI(k)-XI1(k))

      !*** Augment sums for Qsca and g=<cos(theta)>
          QSCA(k)=QSCA(k)+(2d0*DN+1d0)*(ABS(AN(k))*ABS(AN(k))+ABS(BN(k))*ABS(BN(k)))
          GSCA(k)=GSCA(k)+((2d0*DN+1d0)/(DN*(DN+1d0))) * &
               (DBLE(AN(k))*DBLE(BN(k))+AIMAG(AN(k))*AIMAG(BN(k)))
        END IF
      END DO

      IF(N.GT.1)THEN
        DO k=1,anz
          IF (N <= NSTOP(k)) THEN
            GSCA(k)=GSCA(k)+((DN-1d0)*(DN+1d0)/DN)* &
                 (DBLE(AN1(k))*DBLE(AN(k))+AIMAG(AN1(k))*AIMAG(AN(k))+ &
                 DBLE(BN1(k))*DBLE(BN(k))+AIMAG(BN1(k))*AIMAG(BN(k)))
          END IF
        END DO
      ENDIF

      !*** Now calculate scattering intensity pattern
      DO J=1,NANG
        DO k=1,anz
          IF (N <= NSTOP(k)) THEN
            PI(k,J)=PI1(k,J)
            TAU(k,J)=DN*AMU(J)*PI(k,J)-(DN+1d0)*PI0(k,J)
            S1(k,J)=S1(k,J)+FN*(AN(k)*PI(k,J)+BN(k)*TAU(k,J))
            S2(k,J)=S2(k,J)+FN*(AN(k)*TAU(k,J)+BN(k)*PI(k,J))

      !*** Compute pi_n for next value of n
      !   For each angle J, compute pi_n+1
      !   from PI = pi_n , PI0 = pi_n-1
            PI1(k,J)=((2.*DN+1.)*AMU(J)*PI(k,J)-(DN+1d0)*PI0(k,J))/DN
            PI0(k,J)=PI(k,J)

          END IF
        END DO
      END DO


      DO k=1,anz
        IF (N <= NSTOP(k)) THEN
          PSI0(k)=PSI1(k)
          PSI1(k)=PSI(k)
          APSI1(k)=PSI1(k)
          CHI0(k)=CHI1(k)
          CHI1(k)=CHI(k)

          !     In the next line, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
          !          XI1=CMPLX(APSI1,-CHI1)
          XI1(k)=CMPLX(APSI1(k),CHI1(k),kind=dp)
        END IF
      END DO

    END DO

    !*** Have summed sufficient terms.
    !   Now compute QSCA,QEXT,QBACK,and GSCA (efficiency factors)
    GSCA=2d0*GSCA/QSCA
    QSCA=(2d0/(X*X))*QSCA
    QEXT=(4d0/(X*X))*DBLE(S1(:,1))
    ! QSEIT ohne Ber. der Polarisation, es wird vorausgesetzt, dass das Messgeraet alle Polarisationsrichtungen empfangen kann.
    !  QSEIT=(4d0/(X*X))*(ABS(S1(2:NANG-1))**2*SIN(PHI(2:NANG-1)*degrad_dp)**2 + ABS(S1(NANG(2:NANG-1)))**2*COS(PHI(2:NANG-1)*degrad_dp)**2)
    QBACK=(4d0/(X*X))*ABS(S1(:,NANG))*ABS(S1(:,NANG))

    DEALLOCATE(D)

    RETURN
  END SUBROUTINE BHMIESEITVEC

  SUBROUTINE BHMIESEIT(X,REFREL,QBACK,QEXT,QSCA,S1,S2)

    ! Es werden nur Streuwinkel von 180 Grad und 0 Grad (Extinktion!) berechnet.
    ! Erweiterungen auf bel. Streurichtung aber sehr einfach moeglich!!!
    ! Siehe unten, Variablen theta() und NANG erweitern.

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0d0), NANG = 2

    REAL(KIND=dp) :: GSCA,QBACK,QEXT,QSCA,X
    COMPLEX(kind=dp) :: REFREL
    COMPLEX(kind=dp) :: S1(NANG),S2(NANG)

    INTEGER :: J,N,NSTOP,NMX,NN
    REAL(KIND=dp) :: APSI,APSI1,CHI,CHI0,CHI1,FN,XSTOP,YMOD
    REAL(KIND=dp) ::  AMU(NANG),PI(NANG),PI0(NANG),PI1(NANG),TAU(NANG),THETA(NANG)
    REAL(KIND=dp) :: PSI0,PSI1,PSI,DN
    COMPLEX(kind=dp) :: AN,AN1,BN,BN1,XI,XI1,Y
    COMPLEX(kind=dp), ALLOCATABLE :: D(:)
!===========================================================================
!
! Singularity for X = 0.0! The calling program has to take care of this!
!
!===========================================================================
    !***********************************************************************
    ! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
    !    to calculate scattering and absorption by a homogenous isotropic
    !    sphere.
    ! Given:
    !    X = 2*pi*a/lambda
    !    REFREL = (complex refr. index of sphere)/(real index of medium)
    !    NANG = number of angles between 0 and 90 degrees
    !           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
    !           if called with NANG<2, will set NANG=2 and will compute
    !           scattering for theta=0,90,180.
    ! Returns:
    !    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
    !                                scatt. E perp. to scatt. plane)
    !    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
    !                               scatt. E parr. to scatt. plane)
    !   QEXT = C_ext/pi*a**2 = efficiency factor for extinction
    !   QSCA = C_sca/pi*a**2 = efficiency factor for scattering
    !   QBACK = (dC_sca/domega)/pi*a**2
    !         = backscattering efficiency
    !   GSCA = <cos(theta)> for scattering
    !
    !Original program taken from Bohren and Huffman (1983), Appendix A
    !Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
    !in order to compute <cos(theta)>
    !91/05/07 (BTD): Modified to allow NANG=1
    !91/08/15 (BTD): Corrected error (failure to initialize P)
    !91/08/15 (BTD): Modified to enhance vectorizability.
    !91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
    !91/08/15 (BTD): Changed definition of QBACK.
    !92/01/08 (BTD): Note that this version has been superceded by
    !                fully double precision version = bhmie.f which,
    !                unfortunately, is not standard f77.
    !                However, retain this in case standard f77 version
    !                is required for porting to some other system.
    !2005/11/07 U. Blahak: transferred to F90 and full double precision,
    !                major modifications
    !***********************************************************************

    Y=X*REFREL
    YMOD=ABS(Y)
    !
    !*** Series expansion terminated after NSTOP terms
    !   Logarithmic derivatives calculated from NMX on down
    XSTOP=FLOOR(X+4d0*X**third+2d0)+1d0
    !*** Original code:
    !     NMX=AMAX1(XSTOP,YMOD)+15
    !     NSTOP=XSTOP
    !*** Experimental code:
    NMX=MAX(XSTOP,YMOD)+15
    NSTOP=XSTOP
    ALLOCATE(D(NMX))

    !*** Require NANG.GE.1 in order to calculate scattering intensities
!!! Hier werden die Streuwinkel theta (zenitwinkel zur Vorwaertsstreurichtung)
!!! gesetzt:
    THETA(1) = 0d0 * degrad_dp
!!! THETA(...) = ... ! THETA prinzipiell erweiterbar auf beliebige Winkel. NANG entsprechend erhoehen!!!
    THETA(2) = 180d0 * degrad_dp
    DO j=1,NANG
      AMU(J)=COS(THETA(J))
    END DO
    DO J=1,NANG
      PI0(J)=0d0
      PI1(J)=1d0
    END DO
    NN=NANG
    DO J=1,NANG
      S1(J)=CMPLX(0d0,0d0,kind=dp)
      S2(J)=CMPLX(0d0,0d0,kind=dp)
    END DO

    !*** Logarithmic derivative D(J) calculated by downward recurrence
    !   beginning with initial value (0.,0.) at J=NMX

    D(NMX)=CMPLX(0d0,0d0,kind=dp)
    NN=NMX-1
    DO N=1,NN
      DN=NMX-N+1
      D(NMX-N)=(DN/Y)-(1d0/(D(NMX-N+1)+DN/Y))
    END DO

    !*** Riccati-Bessel functions with real argument X
    !   calculated by upward recurrence

    PSI0=COS(X)
    PSI1=SIN(X)
    CHI0=-SIN(X)
    CHI1=COS(X)
    !APSI0 never used, so this line removed from program:
    !     APSI0=PSI0
    APSI1=PSI1
    !XI0 never used, so this line removed from program:
    !     In the next two lines, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
    !     XI0=CMPLX(APSI0,-CHI0)
    !      XI1=CMPLX(APSI1,-CHI1)
    XI1=CMPLX(APSI1,CHI1,kind=dp)
    QSCA=0d0
    GSCA=0d0
    DO N=1,NSTOP
      DN=N
      FN=(2d0*DN+1d0)/(DN*(DN+1d0))
      PSI=(2d0*DN-1d0)*PSI1/X-PSI0
      APSI=PSI
      CHI=(2d0*DN-1d0)*CHI1/X-CHI0

      !     In the next line, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
      !          XI=CMPLX(APSI,-CHI)
      XI=CMPLX(APSI,CHI,kind=dp)

      !*** Store previous values of AN and BN for use
      !   in computation of g=<cos(theta)>
      IF(N.GT.1)THEN
        AN1=AN
        BN1=BN
      ENDIF

      !*** Compute AN and BN:
      AN=(D(N)/REFREL+DN/X)*APSI-APSI1
      AN=AN/((D(N)/REFREL+DN/X)*XI-XI1)
      BN=(REFREL*D(N)+DN/X)*APSI-APSI1
      BN=BN/((REFREL*D(N)+DN/X)*XI-XI1)

      !*** Augment sums for Qsca and g=<cos(theta)>
      QSCA=QSCA+(2d0*DN+1d0)*(ABS(AN)*ABS(AN)+ABS(BN)*ABS(BN))
      GSCA=GSCA+((2d0*DN+1d0)/(DN*(DN+1d0))) * &
           (DBLE(AN)*DBLE(BN)+AIMAG(AN)*AIMAG(BN))
      IF(N.GT.1)THEN
        GSCA=GSCA+((DN-1d0)*(DN+1d0)/DN)* &
             (DBLE(AN1)*DBLE(AN)+AIMAG(AN1)*AIMAG(AN)+ &
             DBLE(BN1)*DBLE(BN)+AIMAG(BN1)*AIMAG(BN))
      ENDIF

      !*** Now calculate scattering intensity pattern
      DO J=1,NANG
        PI(J)=PI1(J)
        TAU(J)=DN*AMU(J)*PI(J)-(DN+1d0)*PI0(J)
        S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
        S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
      END DO

      PSI0=PSI1
      PSI1=PSI
      APSI1=PSI1
      CHI0=CHI1
      CHI1=CHI

      !     In the next line, the minus-sign is changed to switch to convention of REFREL with negative imaginary part
      !          XI1=CMPLX(APSI1,-CHI1)
      XI1=CMPLX(APSI1,CHI1,kind=dp)

      !*** Compute pi_n for next value of n
      !   For each angle J, compute pi_n+1
      !   from PI = pi_n , PI0 = pi_n-1
      DO J=1,NANG
        PI1(J)=((2.*DN+1.)*AMU(J)*PI(J)-(DN+1.)*PI0(J))/DN
        PI0(J)=PI(J)
      END DO

    END DO

    !*** Have summed sufficient terms.
    !   Now compute QSCA,QEXT,QBACK,and GSCA (efficiency factors)
    GSCA=2d0*GSCA/QSCA
    QSCA=(2d0/(X*X))*QSCA
    QEXT=(4d0/(X*X))*DBLE(S1(1))
    ! QSEIT ohne Ber. der Polarisation, es wird vorausgesetzt, dass das Messgeraet alle Polarisationsrichtungen empfangen kann.
    !  QSEIT=(4d0/(X*X))*(ABS(S1(2:NANG-1))**2*SIN(PHI(2:NANG-1)*degrad_dp)**2 + ABS(S1(NANG(2:NANG-1)))**2*COS(PHI(2:NANG-1)*degrad_dp)**2)
    QBACK=(4d0/(X*X))*ABS(S1(NANG))*ABS(S1(NANG))

    DEALLOCATE(D)

    RETURN
  END SUBROUTINE BHMIESEIT



  !----------------------------------------------------------------------------------------!
  ! 3b) Rekursionsformeln nach Bohren und Huffman (1983) fuer zweischalige Kugeln
  !----------------------------------------------------------------------------------------!

  SUBROUTINE BHCOATBACKVEC(X,Y,RRFRL1,RRFRL2,QBACK,QEXT,QSCA,anz)

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0d0)

    INTEGER :: anz
    COMPLEX(kind=dp) :: RRFRL1(anz),RRFRL2(anz)

    INTEGER :: IFLAG(anz),N,NSTOP(anz),K
    REAL(KIND=dp) :: CHI0Y(anz),CHI1Y(anz),CHIY(anz),DEL,PSI0Y(anz),&
         PSI1Y(anz),PSIY(anz),QEXT(anz),QBACK(anz),QSCA(anz),X(anz),Y(anz),YSTOP(anz),RN
    COMPLEX(kind=dp) :: AMESS1,AMESS2,AMESS3,AMESS4,AN(anz),ANCAP, &
         BN(anz),BNCAP,BRACK(anz), &
         CHI0X2(anz),CHI0Y2(anz),CHI1X2(anz),CHI1Y2(anz),&
         CHIX2(anz),CHIPX2(anz),CHIPY2(anz),CHIY2(anz),CRACK(anz), &
         D0X1(anz),D0X2(anz),D0Y2(anz),D1X1(anz),D1X2(anz),D1Y2(anz),DNBAR,GNBAR,II, &
         REFREL(anz),RFREL1(anz),RFREL2(anz), &
         XBACK(anz),XI0Y(anz),XI1Y(anz),XIY(anz), &
         X1(anz),X2(anz),Y2(anz)
!===========================================================================
!
! - Singularity for X = 0.0, Y = 0.0. The calling program has to deal with this.
!
! - Use of optimization may slightly change the results, i.e., ifort -O2
!   or all optimization options for sxf90.
!
!===========================================================================
    !***********************************************************************
    !
    ! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back for coated sphere.
    ! All bessel functions computed by upward recurrence.
    ! Input:
    !        X = 2*PI*RCORE*REFMED/WAVEL
    !        Y = 2*PI*RMANT*REFMED/WAVEL
    !        RFREL1 = REFCOR/REFMED
    !        RFREL2 = REFMAN/REFMED
    ! where  REFCOR = complex refr.index of core)
    !        REFMAN = complex refr.index of mantle)
    !        REFMED = real refr.index of medium)
    !        RCORE = radius of core
    !        RMANT = radius of mantle
    !        WAVEL = wavelength of light in ambient medium
    !
    ! Routine BHCOAT is taken from Bohren & Huffman (1983)
    ! Obtained from C.L.Joseph
    !
    ! history:
    ! 92/11/24 (BTD) Explicit declaration of all variables
    ! 2005/11/07 U.Blahak completely rewritten in F90
    !***********************************************************************

    DEL = 1d-8
    II = CMPLX(0.D0,1.D0,kind=dp)
    RFREL1=RRFRL1
    RFREL2=RRFRL2
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    ystop = y + 4d0*y**third + 2d0
    refrel = rfrel2/rfrel1
    nstop = ystop
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    d0x1 = COS(x1)/SIN(x1)
    d0x2 = COS(x2)/SIN(x2)
    d0y2 = COS(y2)/SIN(y2)
    psi0y = COS(y)
    psi1y = SIN(y)
    chi0y = -SIN(y)
    chi1y = COS(y)
    !     The minus-sign in the next 2 lines was changed to + to switch to the refractive index convention with negative imaginary part
    !     xi0y = psi0y-II*chi0y
    !     xi1y = psi1y-II*chi1y
    xi0y = psi0y+II*chi0y
    xi1y = psi1y+II*chi1y
    chi0y2 = -SIN(y2)
    chi1y2 = COS(y2)
    chi0x2 = -SIN(x2)
    chi1x2 = COS(x2)
    qsca = 0d0
    qext = 0d0
    xback = CMPLX(0d0,0d0,kind=dp)

    iflag = 0

    DO n=1,MAXVAL(nstop)

      rn = n

      DO k=1,anz
        IF (n <= nstop(k)) THEN
          psiy(k) = (2d0*rn-1d0)*psi1y(k)/y(k) - psi0y(k)
          chiy(k) = (2d0*rn-1d0)*chi1y(k)/y(k) - chi0y(k)
          xiy(k) = psiy(k)+II*chiy(k)
          d1y2(k) = 1d0/(rn/y2(k)-d0y2(k)) - rn/y2(k)
        END IF
      END DO

      DO k=1,anz
        IF (n <= nstop(k)) THEN
          IF (iflag(k) .EQ. 0) THEN
            d1x1(k) = 1d0/(rn/x1(k)-d0x1(k)) - rn/x1(k)
            d1x2(k) = 1d0/(rn/x2(k)-d0x2(k)) - rn/x2(k)
            chix2(k) = (2d0*rn - 1d0)*chi1x2(k)/x2(k) - chi0x2(k)
            chiy2(k) = (2d0*rn - 1d0)*chi1y2(k)/y2(k) - chi0y2(k)
            chipx2(k) = chi1x2(k) - rn*chix2(k)/x2(k)
            chipy2(k) = chi1y2(k) - rn*chiy2(k)/y2(k)
            ancap = refrel(k)*d1x1(k) - d1x2(k)
            ancap = ancap/(refrel(k)*d1x1(k)*chix2(k) - chipx2(k))
            ancap = ancap/(chix2(k)*d1x2(k) - chipx2(k))
            brack(k) = ancap*(chiy2(k)*d1y2(k) - chipy2(k))
            bncap = refrel(k)*d1x2(k) - d1x1(k)
            bncap = bncap/(refrel(k)*chipx2(k) - d1x1(k)*chix2(k))
            bncap = bncap/(chix2(k)*d1x2(k) - chipx2(k))
            crack(k) = bncap*(chiy2(k)*d1y2(k) - chipy2(k))
            amess1 = brack(k)*chipy2(k)
            amess2 = brack(k)*chiy2(k)
            amess3 = crack(k)*chipy2(k)
            amess4 = crack(k)*chiy2(k)
            IF (ABS(amess1) .LE. del*ABS(d1y2(k)) .AND. &
                 ABS(amess2) .LE. del .AND. &
                 ABS(amess3) .LE. del*ABS(d1y2(k)) .AND. &
                 ABS(amess4) .LE. del) THEN
              brack(k) = CMPLX(0d0,0d0,kind=dp)
              crack(k) = CMPLX(0d0,0d0,kind=dp)
              iflag(k) = 1
            END IF
          END IF
          dnbar = d1y2(k) - brack(k)*chipy2(k)
          dnbar = dnbar/(1d0-brack(k)*chiy2(k))
          gnbar = d1y2(k) - crack(k)*chipy2(k)
          gnbar = gnbar/(1d0-crack(k)*chiy2(k))
          an(k) = (dnbar/rfrel2(k) + rn/y(k))*psiy(k) - psi1y(k)
          an(k) = an(k)/((dnbar/rfrel2(k)+rn/y(k))*xiy(k)-xi1y(k))
          bn(k) = (rfrel2(k)*gnbar + rn/y(k))*psiy(k) - psi1y(k)
          bn(k) = bn(k)/((rfrel2(k)*gnbar+rn/y(k))*xiy(k)-xi1y(k))
          qsca(k) = qsca(k) + 0.5d0*(2d0*rn+1d0)*(ABS(an(k))*ABS(an(k))+ABS(bn(k))*ABS(bn(k)))
          xback(k) = xback(k) + 0.5d0*(2d0*rn+1d0)*(-1)**n*(an(k)-bn(k))
          qext(k) = qext(k) + 0.5d0*(2d0*rn + 1d0)*DBLE(an(k)+bn(k))
          psi0y(k) = psi1y(k)
          psi1y(k) = psiy(k)
          chi0y(k) = chi1y(k)
          chi1y(k) = chiy(k)
          xi1y(k) = psi1y(k)+II*chi1y(k)
          chi0x2(k) = chi1x2(k)
          chi1x2(k) = chix2(k)
          chi0y2(k) = chi1y2(k)
          chi1y2(k) = chiy2(k)
          d0x1(k) = d1x1(k)
          d0x2(k) = d1x2(k)
          d0y2(k) = d1y2(k)
        END IF
      END DO
    END DO

    QSCA = (4d0/(y*y))*QSCA
    QEXT = (4d0/(y*y))*QEXT
    qback = (4d0/(y*y))*ABS(XBACK)*ABS(XBACK)
    RETURN
  END SUBROUTINE BHCOATBACKVEC

  SUBROUTINE BHCOATBACK(X,Y,RRFRL1,RRFRL2,QBACK,QEXT,QSCA)

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0d0)

    COMPLEX(kind=dp) :: RRFRL1,RRFRL2

    INTEGER :: IFLAG,N,NSTOP
    REAL(KIND=dp) :: CHI0Y,CHI1Y,CHIY,DEL,PSI0Y,PSI1Y,PSIY,QEXT,QBACK,RN,QSCA,X,Y,YSTOP
    COMPLEX(kind=dp) :: AMESS1,AMESS2,AMESS3,AMESS4,AN,ANCAP, &
         BN,BNCAP,BRACK, &
         CHI0X2,CHI0Y2,CHI1X2,CHI1Y2,CHIX2,CHIPX2,CHIPY2,CHIY2,CRACK, &
         D0X1,D0X2,D0Y2,D1X1,D1X2,D1Y2,DNBAR,GNBAR,II, &
         REFREL,RFREL1,RFREL2, &
         XBACK,XI0Y,XI1Y,XIY, &
         X1,X2,Y2
    !***********************************************************************
    !
    ! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back for coated sphere.
    ! All bessel functions computed by upward recurrence.
    ! Input:
    !        X = 2*PI*RCORE*REFMED/WAVEL
    !        Y = 2*PI*RMANT*REFMED/WAVEL
    !        RFREL1 = REFCOR/REFMED
    !        RFREL2 = REFMAN/REFMED
    ! where  REFCOR = complex refr.index of core)
    !        REFMAN = complex refr.index of mantle)
    !        REFMED = real refr.index of medium)
    !        RCORE = radius of core
    !        RMANT = radius of mantle
    !        WAVEL = wavelength of light in ambient medium
    !
    ! Routine BHCOAT is taken from Bohren & Huffman (1983)
    ! Obtained from C.L.Joseph
    !
    ! history:
    ! 92/11/24 (BTD) Explicit declaration of all variables
    ! 2005/11/07 U.Blahak completely rewritten in F90
    !***********************************************************************

    DEL = 1d-8
    II = CMPLX(0.D0,1.D0,kind=dp)
    RFREL1=RRFRL1
    RFREL2=RRFRL2
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    ystop = y + 4d0*y**third + 2d0
    refrel = rfrel2/rfrel1
    nstop = ystop
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    d0x1 = COS(x1)/SIN(x1)
    d0x2 = COS(x2)/SIN(x2)
    d0y2 = COS(y2)/SIN(y2)
    psi0y = COS(y)
    psi1y = SIN(y)
    chi0y = -SIN(y)
    chi1y = COS(y)
    !     The minus-sign in the next 2 lines was changed to + to switch to the refractive index convention with negative imaginary part
    !     xi0y = psi0y-II*chi0y
    !     xi1y = psi1y-II*chi1y
    xi0y = psi0y+II*chi0y
    xi1y = psi1y+II*chi1y
    chi0y2 = -SIN(y2)
    chi1y2 = COS(y2)
    chi0x2 = -SIN(x2)
    chi1x2 = COS(x2)
    qsca = 0d0
    qext = 0d0
    xback = CMPLX(0d0,0d0,kind=dp)

    iflag = 0
    DO n=1,nstop

      rn = n
      psiy = (2d0*rn-1d0)*psi1y/y - psi0y
      chiy = (2d0*rn-1d0)*chi1y/y - chi0y
      xiy = psiy+II*chiy
      d1y2 = 1d0/(rn/y2-d0y2) - rn/y2
      IF (iflag .EQ. 0) THEN
        d1x1 = 1d0/(rn/x1-d0x1) - rn/x1
        d1x2 = 1d0/(rn/x2-d0x2) - rn/x2
        chix2 = (2d0*rn - 1d0)*chi1x2/x2 - chi0x2
        chiy2 = (2d0*rn - 1d0)*chi1y2/y2 - chi0y2
        chipx2 = chi1x2 - rn*chix2/x2
        chipy2 = chi1y2 - rn*chiy2/y2
        ancap = refrel*d1x1 - d1x2
        ancap = ancap/(refrel*d1x1*chix2 - chipx2)
        ancap = ancap/(chix2*d1x2 - chipx2)
        brack = ancap*(chiy2*d1y2 - chipy2)
        bncap = refrel*d1x2 - d1x1
        bncap = bncap/(refrel*chipx2 - d1x1*chix2)
        bncap = bncap/(chix2*d1x2 - chipx2)
        crack = bncap*(chiy2*d1y2 - chipy2)
        amess1 = brack*chipy2
        amess2 = brack*chiy2
        amess3 = crack*chipy2
        amess4 = crack*chiy2
        IF (ABS(amess1) .LE. del*ABS(d1y2) .AND. &
             ABS(amess2) .LE. del .AND. &
             ABS(amess3) .LE. del*ABS(d1y2) .AND. &
             ABS(amess4) .LE. del) THEN
          brack = CMPLX(0d0,0d0,kind=dp)
          crack = CMPLX(0d0,0d0,kind=dp)
          iflag = 1
        END IF
      END IF
      dnbar = d1y2 - brack*chipy2
      dnbar = dnbar/(1d0-brack*chiy2)
      gnbar = d1y2 - crack*chipy2
      gnbar = gnbar/(1d0-crack*chiy2)
      an = (dnbar/rfrel2 + rn/y)*psiy - psi1y
      an = an/((dnbar/rfrel2+rn/y)*xiy-xi1y)
      bn = (rfrel2*gnbar + rn/y)*psiy - psi1y
      bn = bn/((rfrel2*gnbar+rn/y)*xiy-xi1y)
      qsca = qsca + 0.5d0*(2d0*rn+1d0)*(ABS(an)*ABS(an)+ABS(bn)*ABS(bn))
      xback = xback + 0.5d0*(2d0*rn+1d0)*(-1)**n*(an-bn)
      qext = qext + 0.5d0*(2d0*rn + 1d0)*DBLE(an+bn)
      psi0y = psi1y
      psi1y = psiy
      chi0y = chi1y
      chi1y = chiy
      xi1y = psi1y+II*chi1y
      chi0x2 = chi1x2
      chi1x2 = chix2
      chi0y2 = chi1y2
      chi1y2 = chiy2
      d0x1 = d1x1
      d0x2 = d1x2
      d0y2 = d1y2
    END DO

    QSCA = (4d0/(y*y))*QSCA
    QEXT = (4d0/(y*y))*QEXT
    qback = (4d0/(y*y))*ABS(XBACK)*ABS(XBACK)
    RETURN
  END SUBROUTINE BHCOATBACK

  SUBROUTINE BHCOATBACKVECS(X,Y,RRFRL1,RRFRL2,QBACK,QEXT,QSCA,S1,S2,anz)

    !***********************************************************************
    !
    ! JM190705:
    ! This is SUBROUTINE BHCOATBACKVEC extended to also return scattering
    ! amplitudes (needed for polarimetric parameter calculations).
    !
    ! Extention is based on theory as shown in Guzzi98: calculate S1 and S2
    ! as done for homogenous spheres, just applying the coated sphere
    ! modified an and bn (ie the an and bn as calculated and applied for
    ! QBACK,QEXT,QSCA calculation in the original BHCOATBACKVEC routine).
    ! The code required in addition to what was already in BHCOATBACKVEC (ie
    ! the angular dependency and S1/2 sum up parts) have been copied from
    ! the homogeneous sphere code in subroutine BHMIESEITVEC.
    !
    !***********************************************************************

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0.0d0)

    INTEGER :: anz
    COMPLEX(kind=dp) :: RRFRL1(anz),RRFRL2(anz)

    INTEGER :: IFLAG(anz),N,NSTOP(anz),K
    REAL(KIND=dp) :: CHI0Y(anz),CHI1Y(anz),CHIY(anz),DEL,PSI0Y(anz),&
         PSI1Y(anz),PSIY(anz),QEXT(anz),QBACK(anz),QSCA(anz),X(anz),Y(anz),YSTOP(anz),RN
    COMPLEX(kind=dp) :: AMESS1,AMESS2,AMESS3,AMESS4,AN(anz),ANCAP, &
         BN(anz),BNCAP,BRACK(anz), &
         CHI0X2(anz),CHI0Y2(anz),CHI1X2(anz),CHI1Y2(anz),&
         CHIX2(anz),CHIPX2(anz),CHIPY2(anz),CHIY2(anz),CRACK(anz), &
         D0X1(anz),D0X2(anz),D0Y2(anz),D1X1(anz),D1X2(anz),D1Y2(anz),DNBAR,GNBAR,II, &
         REFREL(anz),RFREL1(anz),RFREL2(anz), &
         XBACK(anz),XI0Y(anz),XI1Y(anz),XIY(anz), &
         X1(anz),X2(anz),Y2(anz)

    INTEGER, PARAMETER :: NANG = 2
    INTEGER :: J
    COMPLEX(kind=dp) :: S1(anz,NANG),S2(anz,NANG)
    REAL(KIND=dp) :: AMU(NANG),PI(anz,NANG),PI0(anz,NANG),PI1(anz,NANG),TAU(anz,NANG),THETA(NANG)
    REAL(KIND=dp) :: DN,FN

!===========================================================================
!
! - Singularity for X = 0.0, Y = 0.0. The calling program has to deal with this.
!
! - Use of optimization may slightly change the results, i.e., ifort -O2
!   or all optimization options for sxf90.
!
!===========================================================================
    !***********************************************************************
    !
    ! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back for coated sphere.
    ! All bessel functions computed by upward recurrence.
    ! Input:
    !        X = 2*PI*RCORE*REFMED/WAVEL
    !        Y = 2*PI*RMANT*REFMED/WAVEL
    !        RFREL1 = REFCOR/REFMED
    !        RFREL2 = REFMAN/REFMED
    ! where  REFCOR = complex refr.index of core)
    !        REFMAN = complex refr.index of mantle)
    !        REFMED = real refr.index of medium)
    !        RCORE = radius of core
    !        RMANT = radius of mantle
    !        WAVEL = wavelength of light in ambient medium
    !
    ! Routine BHCOAT is taken from Bohren & Huffman (1983)
    ! Obtained from C.L.Joseph
    !
    ! History:
    ! 92/11/24 (BTD) Explicit declaration of all variables
    ! 2005/11/07 U.Blahak completely rewritten in F90
    !***********************************************************************

    DEL = 1d-8
    II = CMPLX(0.D0,1.D0,kind=dp)
    RFREL1=RRFRL1
    RFREL2=RRFRL2
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    ystop = y + 4.*y**0.3333 + 2.0
    refrel = rfrel2/rfrel1
    nstop = ystop
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    d0x1 = COS(x1)/SIN(x1)
    d0x2 = COS(x2)/SIN(x2)
    d0y2 = COS(y2)/SIN(y2)
    psi0y = COS(y)
    psi1y = SIN(y)
    chi0y = -SIN(y)
    chi1y = COS(y)
    !     The minus-sign in the next 2 lines was changed to + to switch to the refractive index convention with negative imaginary part
    !     xi0y = psi0y-II*chi0y
    !     xi1y = psi1y-II*chi1y
    xi0y = psi0y+II*chi0y
    xi1y = psi1y+II*chi1y
    chi0y2 = -SIN(y2)
    chi1y2 = COS(y2)
    chi0x2 = -SIN(x2)
    chi1x2 = COS(x2)
    qsca = 0.0
    qext = 0.0
    xback = CMPLX(0.0d0,0.0d0,kind=dp)

    THETA(1) = 0.0d0 * degrad_dp
    THETA(2) = 180.0d0 * degrad_dp
    DO J=1,NANG
      AMU(J)=COS(THETA(J))
    END DO
    DO J=1,NANG
      DO k=1,anz
        PI0(k,J)=0.d0
        PI1(k,J)=1.d0
      END DO
    END DO
    DO J=1,NANG
      DO k=1,anz
        S1(k,J)=CMPLX(0.d0,0.d0,kind=dp)
        S2(k,J)=CMPLX(0.d0,0.d0,kind=dp)
      END DO
    END DO

    iflag = 0

    DO n=1,MAXVAL(nstop)

      rn = n

!$NEC novector
      DO k=1,anz
        IF (n <= nstop(k)) THEN
          psiy(k) = (2.0*rn-1.)*psi1y(k)/y(k) - psi0y(k)
          chiy(k) = (2.0*rn-1.)*chi1y(k)/y(k) - chi0y(k)
          xiy(k) = psiy(k)+II*chiy(k)
          d1y2(k) = 1.0/(rn/y2(k)-d0y2(k)) - rn/y2(k)
        END IF
      END DO

      DO k=1,anz
        IF (n <= nstop(k)) THEN
          IF (iflag(k) .EQ. 0) THEN
            d1x1(k) = 1.0/(rn/x1(k)-d0x1(k)) - rn/x1(k)
            d1x2(k) = 1.0/(rn/x2(k)-d0x2(k)) - rn/x2(k)
            chix2(k) = (2.0*rn - 1.0)*chi1x2(k)/x2(k) - chi0x2(k)
            chiy2(k) = (2.0*rn - 1.0)*chi1y2(k)/y2(k) - chi0y2(k)
            chipx2(k) = chi1x2(k) - rn*chix2(k)/x2(k)
            chipy2(k) = chi1y2(k) - rn*chiy2(k)/y2(k)
            ancap = refrel(k)*d1x1(k) - d1x2(k)
            ancap = ancap/(refrel(k)*d1x1(k)*chix2(k) - chipx2(k))
            ancap = ancap/(chix2(k)*d1x2(k) - chipx2(k))
            brack(k) = ancap*(chiy2(k)*d1y2(k) - chipy2(k))
            bncap = refrel(k)*d1x2(k) - d1x1(k)
            bncap = bncap/(refrel(k)*chipx2(k) - d1x1(k)*chix2(k))
            bncap = bncap/(chix2(k)*d1x2(k) - chipx2(k))
            crack(k) = bncap*(chiy2(k)*d1y2(k) - chipy2(k))
            amess1 = brack(k)*chipy2(k)
            amess2 = brack(k)*chiy2(k)
            amess3 = crack(k)*chipy2(k)
            amess4 = crack(k)*chiy2(k)
            IF (ABS(amess1) .LE. del*ABS(d1y2(k)) .AND. &
                 ABS(amess2) .LE. del .AND. &
                 ABS(amess3) .LE. del*ABS(d1y2(k)) .AND. &
                 ABS(amess4) .LE. del) THEN
              brack(k) = CMPLX(0.d0,0.d0,kind=dp)
              crack(k) = CMPLX(0.d0,0.d0,kind=dp)
              iflag(k) = 1
            END IF
          END IF
          dnbar = d1y2(k) - brack(k)*chipy2(k)
          dnbar = dnbar/(1.0-brack(k)*chiy2(k))
          gnbar = d1y2(k) - crack(k)*chipy2(k)
          gnbar = gnbar/(1.0-crack(k)*chiy2(k))
          an(k) = (dnbar/rfrel2(k) + rn/y(k))*psiy(k) - psi1y(k)
          an(k) = an(k)/((dnbar/rfrel2(k)+rn/y(k))*xiy(k)-xi1y(k))
          bn(k) = (rfrel2(k)*gnbar + rn/y(k))*psiy(k) - psi1y(k)
          bn(k) = bn(k)/((rfrel2(k)*gnbar+rn/y(k))*xiy(k)-xi1y(k))
          qsca(k) = qsca(k) + 0.5*(2.0*rn+1.0)*(ABS(an(k))*ABS(an(k))+ABS(bn(k))*ABS(bn(k)))
          xback(k) = xback(k) + 0.5*(2.0*rn+1.0)*(-1)**n*(an(k)-bn(k))
          qext(k) = qext(k) + 0.5*(2.0*rn + 1.0)*DBLE(an(k)+bn(k))
          psi0y(k) = psi1y(k)
          psi1y(k) = psiy(k)
          chi0y(k) = chi1y(k)
          chi1y(k) = chiy(k)
          xi1y(k) = psi1y(k)+II*chi1y(k)
          chi0x2(k) = chi1x2(k)
          chi1x2(k) = chix2(k)
          chi0y2(k) = chi1y2(k)
          chi1y2(k) = chiy2(k)
          d0x1(k) = d1x1(k)
          d0x2(k) = d1x2(k)
          d0y2(k) = d1y2(k)
        END IF
      END DO

      DN = n
      FN=(2.d0*DN+1.d0)/(DN*(DN+1.d0))
      DO J=1,NANG
        DO k=1,anz
          IF (N <= NSTOP(k)) THEN
            PI(k,J)=PI1(k,J)
            TAU(k,J)=DN*AMU(J)*PI(k,J)-(DN+1.d0)*PI0(k,J)
            S1(k,J)=S1(k,J)+FN*(AN(k)*PI(k,J)+BN(k)*TAU(k,J))
            S2(k,J)=S2(k,J)+FN*(AN(k)*TAU(k,J)+BN(k)*PI(k,J))

            PI1(k,J)=((2.*DN+1.)*AMU(J)*PI(k,J)-(DN+1.)*PI0(k,J))/DN
            PI0(k,J)=PI(k,J)
          END IF
        END DO
      END DO

    END DO

    QSCA = (4.0/(y*y))*QSCA
    QEXT = (4.0/(y*y))*QEXT
    qback = (4.0/(y*y))*(ABS(XBACK))**2

    RETURN
  END SUBROUTINE BHCOATBACKVECS

  SUBROUTINE BHCOATBACKS(X,Y,RRFRL1,RRFRL2,QBACK,QEXT,QSCA,S1,S2)

    !***********************************************************************
    !
    ! JM190705:
    ! This is SUBROUTINE BHCOATBACK extended to also return scattering
    ! amplitudes (needed for polarimetric parameter calculations).
    ! For debugging it also (temporarily) provides Qext and Qback calculated
    ! from S1.
    !
    ! Extention is based on theory as shown in Guzzi98: calculate S1 and S2
    ! as done for homogenous spheres, just applying the coated sphere
    ! modified an and bn (ie the an and bn as calculated and applied for
    ! QBACK,QEXT,QSCA calculation in the original BHCOATBACK routine).
    ! The code required in addition to what was already in BHCOATBACK (ie
    ! the angular dependency and S1/2 sum up parts) has been copied from
    ! the homogeneous sphere code in subroutine BHMIESEIT.
    !
    !***********************************************************************

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0.0d0)

    COMPLEX(kind=dp) :: RRFRL1,RRFRL2

    INTEGER :: IFLAG,N,NSTOP
    REAL(KIND=dp) :: CHI0Y,CHI1Y,CHIY,DEL,PSI0Y,PSI1Y,PSIY,QEXT,QBACK,RN,QSCA,X,Y,YSTOP
    COMPLEX(kind=dp) :: AMESS1,AMESS2,AMESS3,AMESS4,AN,ANCAP, &
         BN,BNCAP,BRACK, &
         CHI0X2,CHI0Y2,CHI1X2,CHI1Y2,CHIX2,CHIPX2,CHIPY2,CHIY2,CRACK, &
         D0X1,D0X2,D0Y2,D1X1,D1X2,D1Y2,DNBAR,GNBAR,II, &
         REFREL,RFREL1,RFREL2, &
         XBACK,XI0Y,XI1Y,XIY, &
         X1,X2,Y2

    INTEGER, PARAMETER :: NANG = 2
    INTEGER :: J
    COMPLEX(kind=dp) :: S1(NANG),S2(NANG)
    REAL(KIND=dp) :: AMU(NANG),PI(NANG),PI0(NANG),PI1(NANG),TAU(NANG),THETA(NANG)
    REAL(KIND=dp) :: DN,FN

    !***********************************************************************
    !
    ! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back for coated sphere.
    ! All bessel functions computed by upward recurrence.
    ! Input:
    !        X = 2*PI*RCORE*REFMED/WAVEL
    !        Y = 2*PI*RMANT*REFMED/WAVEL
    !        RFREL1 = REFCOR/REFMED
    !        RFREL2 = REFMAN/REFMED
    ! where  REFCOR = complex refr.index of core)
    !        REFMAN = complex refr.index of mantle)
    !        REFMED = real refr.index of medium)
    !        RCORE = radius of core
    !        RMANT = radius of mantle
    !        WAVEL = wavelength of light in ambient medium
    !
    ! Routine BHCOAT is taken from Bohren & Huffman (1983)
    ! Obtained from C.L.Joseph
    !
    ! History:
    ! 92/11/24 (BTD) Explicit declaration of all variables
    ! 2005/11/07 U.Blahak completely rewritten in F90
    !***********************************************************************

    DEL = 1d-8
    II = CMPLX(0.D0,1.D0,kind=dp)
    RFREL1=RRFRL1
    RFREL2=RRFRL2
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    ystop = y + 4.*y**0.3333 + 2.0
    refrel = rfrel2/rfrel1
    nstop = ystop
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    d0x1 = COS(x1)/SIN(x1)
    d0x2 = COS(x2)/SIN(x2)
    d0y2 = COS(y2)/SIN(y2)
    psi0y = COS(y)
    psi1y = SIN(y)
    chi0y = -SIN(y)
    chi1y = COS(y)
    !     The minus-sign in the next 2 lines was changed to + to switch to the refractive index convention with negative imaginary part
    !     xi0y = psi0y-II*chi0y
    !     xi1y = psi1y-II*chi1y
    xi0y = psi0y+II*chi0y
    xi1y = psi1y+II*chi1y
    chi0y2 = -SIN(y2)
    chi1y2 = COS(y2)
    chi0x2 = -SIN(x2)
    chi1x2 = COS(x2)
    qsca = 0.0
    qext = 0.0
    xback = CMPLX(0.0d0,0.0d0,kind=dp)

    THETA(1) = 0.0d0 * degrad_dp
    THETA(2) = 180.0d0 * degrad_dp
    DO J=1,NANG
      AMU(J)=COS(THETA(J))
    END DO
    DO J=1,NANG
      PI0(J)=0.d0
      PI1(J)=1.d0
    END DO
    DO J=1,NANG
      S1(J)=CMPLX(0.d0,0.d0,kind=dp)
      S2(J)=CMPLX(0.d0,0.d0,kind=dp)
    END DO

    iflag = 0
    DO n=1,nstop

      rn = n
      psiy = (2.0*rn-1.)*psi1y/y - psi0y
      chiy = (2.0*rn-1.)*chi1y/y - chi0y
      xiy = psiy+II*chiy
      d1y2 = 1.0/(rn/y2-d0y2) - rn/y2
      IF (iflag .EQ. 0) THEN
        d1x1 = 1.0/(rn/x1-d0x1) - rn/x1
        d1x2 = 1.0/(rn/x2-d0x2) - rn/x2
        chix2 = (2.0*rn - 1.0)*chi1x2/x2 - chi0x2
        chiy2 = (2.0*rn - 1.0)*chi1y2/y2 - chi0y2
        chipx2 = chi1x2 - rn*chix2/x2
        chipy2 = chi1y2 - rn*chiy2/y2
        ancap = refrel*d1x1 - d1x2
        ancap = ancap/(refrel*d1x1*chix2 - chipx2)
        ancap = ancap/(chix2*d1x2 - chipx2)
        brack = ancap*(chiy2*d1y2 - chipy2)
        bncap = refrel*d1x2 - d1x1
        bncap = bncap/(refrel*chipx2 - d1x1*chix2)
        bncap = bncap/(chix2*d1x2 - chipx2)
        crack = bncap*(chiy2*d1y2 - chipy2)
        amess1 = brack*chipy2
        amess2 = brack*chiy2
        amess3 = crack*chipy2
        amess4 = crack*chiy2
        IF (ABS(amess1) .LE. del*ABS(d1y2) .AND. &
             ABS(amess2) .LE. del .AND. &
             ABS(amess3) .LE. del*ABS(d1y2) .AND. &
             ABS(amess4) .LE. del) THEN
          brack = CMPLX(0.d0,0.d0,kind=dp)
          crack = CMPLX(0.d0,0.d0,kind=dp)
          iflag = 1
        END IF
      END IF
      dnbar = d1y2 - brack*chipy2
      dnbar = dnbar/(1.0-brack*chiy2)
      gnbar = d1y2 - crack*chipy2
      gnbar = gnbar/(1.0-crack*chiy2)
      an = (dnbar/rfrel2 + rn/y)*psiy - psi1y
      an = an/((dnbar/rfrel2+rn/y)*xiy-xi1y)
      bn = (rfrel2*gnbar + rn/y)*psiy - psi1y
      bn = bn/((rfrel2*gnbar+rn/y)*xiy-xi1y)
      qsca = qsca + 0.5*(2.0*rn+1.0)*(ABS(an)*ABS(an)+ABS(bn)*ABS(bn))
      xback = xback + 0.5*(2.0*rn+1.0)*(-1)**n*(an-bn)
      qext = qext + 0.5*(2.0*rn + 1.0)*DBLE(an+bn)

      DN = n
      FN=(2.d0*DN+1.d0)/(DN*(DN+1.d0))
      DO J=1,NANG
        PI(J)=PI1(J)
        TAU(J)=DN*AMU(J)*PI(J)-(DN+1.d0)*PI0(J)
        S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
        S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
      END DO
      DO J=1,NANG
        PI1(J)=((2.*DN+1.)*AMU(J)*PI(J)-(DN+1.)*PI0(J))/DN
        PI0(J)=PI(J)
      END DO

      psi0y = psi1y
      psi1y = psiy
      chi0y = chi1y
      chi1y = chiy
      xi1y = psi1y+II*chi1y
      chi0x2 = chi1x2
      chi1x2 = chix2
      chi0y2 = chi1y2
      chi1y2 = chiy2
      d0x1 = d1x1
      d0x2 = d1x2
      d0y2 = d1y2
    END DO

    QSCA = (4.0/(y*y))*QSCA
    QEXT = (4.0/(y*y))*QEXT
    qback = (4.0/(y*y))*(ABS(XBACK))**2

    RETURN
  END SUBROUTINE BHCOATBACKS

  !----------------------------------------------------------------------------------------!
  ! 3c) Ansteuern von BHMIESEIT()
  !----------------------------------------------------------------------------------------!

  SUBROUTINE SPHERE_SCATTER_BH_VEC(D,m,m_medium,lambda,&
      fa,fb,fa0,fb0,anz)
    IMPLICIT NONE

    !..Berechnung des Rueckstreuquerschnitts nach Bohren and Huffman (1983)

    !..Eingabevariablen
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D(anz),lambda,m_medium   ! in [m] angeben, dann C_back in [m^2]
    COMPLEX(kind=dp), INTENT(in) :: m(anz)

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa(anz), fb(anz), fa0(anz), fb0(anz)

    !..Lokale Variablen
    REAL(KIND=dp)    :: alf(anz)
    REAL(KIND=dp)    :: Q_back(anz), Q_ext(anz), Q_sca(anz)
    COMPLEX(kind=dp) :: S1(anz,NANG), S2(anz,NANG)

    !..Dimensionslose Radius
    alf = pi_dp / lambda * D * m_medium

    !..Berechung der Streueffizienz Q_back:
    CALL BHMIESEITVEC(MAX(alf,xeps),m/m_medium,Q_back,Q_ext,Q_sca,S1,S2,anz)

    ! Rueckstreuquerschnitt in der Einheit [Laenge^2]
    WHERE (alf >= xeps)
      fa  = (lambda/(2d0*pi_dp*m_medium))*S1(:,NANG)*nci
      fa0 = (lambda/(2d0*pi_dp*m_medium))*S1(:,1)*nci
   ELSEWHERE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fa0 = fa
    END WHERE

    fb  = fa
    fb0 = fa0

  END SUBROUTINE SPHERE_SCATTER_BH_VEC

  SUBROUTINE SPHERE_SCATTER_BH(D,m,m_medium,lambda,fa,fb,fa0,fb0)
    IMPLICIT NONE

    !..Berechnung des Rueckstreuquerschnitts nach Bohren and Huffman (1983)

    !..Eingabevariablen
    REAL(KIND=dp), INTENT(in)    :: D,lambda,m_medium   ! in [m] angeben, dann C_back in [m^2]
    COMPLEX(kind=dp), INTENT(in) :: m

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa, fb, fa0, fb0

    !..Lokale Variablen
    REAL(KIND=dp)    :: alf
    REAL(KIND=dp)    :: Q_back, Q_ext, Q_sca
    COMPLEX(kind=dp) :: S1(NANG), S2(NANG)

    !..Dimensionslose Radius
    alf = pi_dp / lambda * D * m_medium

    IF (alf >= xeps) THEN
      !..Berechung der Streueffizienz Q_back:
      CALL BHMIESEIT(alf,m/m_medium,Q_back,Q_ext,Q_sca,S1,S2)

      ! Rueckstreuquerschnitt in der Einheit [Laenge^2]
      !C_back = pi_dp * 0.25d0*D**2 *  Q_back !!! blahak, 4.3.2003
      !C_ext = pi_dp * 0.25d0*D**2 *  Q_ext !!! blahak, 4.3.2003
      !C_sca = pi_dp * 0.25d0*D**2 *  Q_sca !!! blahak, 4.3.2003
      fa  = (lambda/(2*pi_dp*m_medium))*S1(NANG)*nci
      fa0 = (lambda/(2*pi_dp*m_medium))*S1(1)*nci

    ELSE

      !C_back = 0d0
      !C_ext = 0d0
      !C_sca = 0d0
      fa  = CMPLX(0.d0,0.d0,kind=dp)
      fa0 = fa

    END IF

    !fa  = sqrt(C_back/(4*pi_dp))
    !fa  = (lambda/(2*pi_dp*m_medium))*S1(NANG)*nci
    !fa0 = (lambda/(2*pi_dp*m_medium))*S1(1)*nci
    fb  = fa
    fb0 = fa0

  END SUBROUTINE SPHERE_SCATTER_BH

  !----------------------------------------------------------------------------------------!
  ! 3c.1) JCS - TMATRIX SCATTERING
  !----------------------------------------------------------------------------------------!

  SUBROUTINE SPHEROID_SCATTER_TMAT_VEC(&
      Dveq,m,m_medium,aspectratio,lambda,fa,fb,fa0,fb0,anz,do_tmat)
    IMPLICIT NONE

    !..Eingabevariablen
    INTEGER, INTENT(in)           :: anz
    REAL(KIND=dp), INTENT(in)     :: Dveq(anz),aspectratio(anz),lambda,m_medium
    COMPLEX(kind=dp), INTENT(in)  :: m(anz)
    LOGICAL, OPTIONAL, INTENT(in) :: do_tmat(anz)

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out)    :: fa(anz),fb(anz),fa0(anz),fb0(anz)

    INTEGER :: i, tmerr(anz)
    LOGICAL :: do_tmat_tmp(anz)

    IF (PRESENT(do_tmat)) THEN
      do_tmat_tmp = do_tmat
    ELSE
      do_tmat_tmp = .TRUE.
    END IF

    tmerr = -1
    DO i=1,anz
      IF (do_tmat_tmp(i)) THEN
        CALL tmatrix(Dveq(i),lambda,m(i),aspectratio(i),fa(i),fb(i),fa0(i),fb0(i),tmerr(i))
      END IF
    END DO

    ! JM190725:
    ! Convert from tmatrix sign conventions to current EMVORADO one.
    ! FIXME: we need to check however, whether this is really
    !        consistent with calc_radar_vars orientation conversion
    !        and for other polarimetric moments.
    !        Also check, whether this is original tmatrix convention or
    !        modified by Jeff (and forgot fb0?)
    ! JM191021:
    ! Removed any sign changings in tmatrix(), instead do ALL of them here.
    WHERE (tmerr .EQ. 0)
      fa  = CONJG(fa)*(-1.d0)  ! fa*(-1.d0)
      fb  = CONJG(fb)          ! fb*(-1.d0)
      fa0 = CONJG(fa0)*(-1.d0) ! fa0*(-1.d0)
      fb0 = CONJG(fb0)*(-1.d0) ! !fb0 = fb0*(-1.d0)
    ELSEWHERE
      fa  = miss_value
      fb  = miss_value
      fa0 = miss_value
      fb0 = miss_value
    END WHERE

  END SUBROUTINE SPHEROID_SCATTER_TMAT_VEC
!!!
  SUBROUTINE SPHEROID_SCATTER_TMAT(Dveq,m,m_medium,aspectratio,lambda,fa,fb,fa0,fb0)
    IMPLICIT NONE

    !..Eingabevariablen
    REAL(KIND=dp), INTENT(in)    :: Dveq,aspectratio,lambda,m_medium
    COMPLEX(kind=dp), INTENT(in) :: m

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa,fb,fa0,fb0

    INTEGER :: tmerr

    CALL tmatrix(Dveq,lambda,m,aspectratio,fa,fb,fa0,fb0,tmerr)

    ! JM190725:
    ! Convert from tmatrix sign conventions to current EMVORADO one.
    ! FIXME: we need to check however, whether this is really
    !        consistent with calc_radar_vars orientation conversion
    !        and for other polarimetric moments.
    ! JM191021:
    ! Removed any sign changings in tmatrix(), instead do ALL of them here.
    IF (tmerr .EQ. 0) THEN
      fa  = CONJG(fa)*(-1.d0)  ! fa*(-1.d0)
      fb  = CONJG(fb)          ! fb*(-1.d0)
      fa0 = CONJG(fa0)*(-1.d0) ! fa0*(-1.d0)
      fb0 = CONJG(fb0)*(-1.d0) ! !fb0 = fb0*(-1.d0)
    ELSE
      fa  = miss_value
      fb  = miss_value
      fa0 = miss_value
      fb0 = miss_value
    END IF

  END SUBROUTINE SPHEROID_SCATTER_TMAT
  !----------------------------------------------------------------------------------------!
  ! 3d) Ansteuern von BHCOATBACKVEC()
  !----------------------------------------------------------------------------------------!

  SUBROUTINE COATEDSPHERE_SCATTER_BH_VEC(D,D_core,m_core,m_shell,m_medium,lambda,&
      fa,fb,fa0,fb0,anz)
    IMPLICIT NONE

    !..Berechnung des Rueckstreuquerschnitts, des Wirkungsquerschnitts
    !  und des effektiven Dielektrizitaetsfaktors nach Bohren and Huffman (1983)

    !..Eingabevariablen
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D(anz),D_core(anz),&
                                    lambda,m_medium
    COMPLEX(kind=dp), INTENT(in) :: m_core(anz),m_shell(anz)

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa(anz),fb(anz),fa0(anz),fb0(anz)

    !..Lokale Variablen
    REAL(KIND=dp)    :: alf(anz),nue(anz),alfhlp(anz),nuehlp(anz)
    REAL(KIND=dp)    :: Q_sca(anz),Q_ext(anz),Q_back(anz),&
                        Q_scahlp(anz),Q_exthlp(anz),Q_backhlp(anz)
    COMPLEX(kind=dp) :: S1(anz,NANG),S2(anz,NANG),S1hlp(anz,NANG),S2hlp(anz,NANG)

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = MAX(pi_dp / lambda * D_core * m_medium , 1d-14)
    nue = pi_dp / lambda * D * m_medium

    alfhlp = MAX(alf,1d-7)
    nuehlp = MAX(MAX(nue,2d-7),alf)

    !CALL BHCOATBACKVEC(alfhlp,nuehlp,m_core/m_medium,m_shell/m_medium,&
    !                   Q_back,Q_ext,Q_sca,anz)
    CALL BHCOATBACKVECS(alfhlp,nuehlp,m_core/m_medium,m_shell/m_medium,&
                        Q_back,Q_ext,Q_sca,S1,S2,anz)

    IF (ANY(alf < 1d-7 .AND. nue >= 2d-7)) THEN
      CALL BHMIESEITVEC(nuehlp,m_shell/m_medium,&
                        Q_backhlp,Q_exthlp,Q_scahlp,S1hlp,S2hlp,anz)
      !WHERE(alf < 1d-7 .AND. nue >= 2d-7)
      !  Q_back = Q_backhlp
      !  Q_ext  = Q_exthlp
      !  Q_sca  = Q_scahlp
      !END WHERE
    END IF

    WHERE (nue >= 2d-7 .AND. nue >= alf)
      WHERE(alf < 1d-7 .AND. nue >= 2d-7)
        fa  = (lambda/(2*pi_dp*m_medium))*S1hlp(:,NANG)*nci
        fa0 = (lambda/(2*pi_dp*m_medium))*S1hlp(:,1)*nci
      ELSEWHERE
        fa  = (lambda/(2*pi_dp*m_medium))*S1(:,NANG)*nci
        fa0 = (lambda/(2*pi_dp*m_medium))*S1(:,1)*nci
      END WHERE
    ELSEWHERE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fa0 = fa
   END WHERE

    fb  = fa
    fb0 = fa0

    RETURN
  END SUBROUTINE COATEDSPHERE_SCATTER_BH_VEC

  SUBROUTINE COATEDSPHERE_SCATTER_BH(D,D_core,m_core,m_shell,m_medium,lambda,&
      fa,fb,fa0,fb0)
    IMPLICIT NONE

    !..Berechnung des Rueckstreuquerschnitts, des Wirkungsquerschnitts
    !  und des effektiven Dielektrizitaetsfaktors nach Bohren and Huffman (1983)

    !..Eingabevariablen
    REAL(KIND=dp), INTENT(in)    :: D,D_core,lambda
    COMPLEX(kind=dp), INTENT(in) :: m_core,m_shell
    REAL(KIND=dp), INTENT(in)    :: m_medium

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa,fb,fa0,fb0

    !..Lokale Variablen
    REAL(KIND=dp)    :: alf,nue
    REAL(KIND=dp)    :: Q_sca,Q_ext,Q_back
    COMPLEX(kind=dp) :: S1(NANG),S2(NANG)

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = MAX(pi_dp / lambda * D_core * m_medium , 1d-14)
    nue = pi_dp / lambda * D * m_medium

    IF (nue >= 2d-7 .AND. nue >= alf) THEN

      !..Berechung der div. Querschnitte:
      IF (alf >= 1d-7) THEN
        !CALL BHCOATBACK(alf,nue,m_core/m_medium,m_shell/m_medium,&
        !                Q_back,Q_ext,Q_sca)
        CALL BHCOATBACKS(alf,nue,m_core/m_medium,m_shell/m_medium,&
                        Q_back,Q_ext,Q_sca,S1,S2)
      ELSE
        CALL BHMIESEIT(nue,m_shell/m_medium,&
                       Q_back,Q_ext,Q_sca,S1,S2)
      END IF

      !C_back =  pi_dp * 0.25d0*D**2 * Q_back
      !C_ext =  pi_dp * 0.25d0*D**2 * Q_ext
      !C_sca =  pi_dp * 0.25d0*D**2 * Q_sca
      fa  = (lambda/(2*pi_dp*m_medium))*S1(NANG)*nci
      fa0 = (lambda/(2*pi_dp*m_medium))*S1(1)*nci

    ELSE

      !C_back = 0d0
      !C_ext = 0d0
      !C_sca = 0d0
      fa  = CMPLX(0.d0,0.d0,kind=dp)
      fa0 = fa

    END IF

    !fa  = (lambda/(2*pi_dp*m_medium))*S1(NANG)*nci
    !fa0 = (lambda/(2*pi_dp*m_medium))*S1(1)*nci
    fb  = fa
    fb0 = fa0

    RETURN
  END SUBROUTINE COATEDSPHERE_SCATTER_BH

  !!!!!
  ! JCS T-matrix coated sphere
  !!!!!
  SUBROUTINE COATEDSPHEROID_SCATTER_TMAT_VEC(&
      D,D_core,m_core,m_shell,m_medium,aspectratio,lambda,&
      fa,fb,fa0,fb0,anz,do_tmat)
    IMPLICIT NONE

    INTEGER, INTENT(in)           :: anz
    REAL(KIND=dp), INTENT(in)     :: D(anz),D_core(anz),aspectratio(anz),lambda,m_medium
    COMPLEX(kind=dp), INTENT(in)  :: m_core(anz),m_shell(anz)
    LOGICAL, OPTIONAL, INTENT(in) :: do_tmat(anz)

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa(anz),fb(anz),fa0(anz),fb0(anz)

    !..Lokale Variablen
    REAL(KIND=dp)    :: alf(anz),nue(anz)
    COMPLEX(kind=dp) :: fahlp(anz),fbhlp(anz),fa0hlp(anz),fb0hlp(anz)
    INTEGER          :: i, tmerr(anz)
    LOGICAL          :: do_tmat_tmp(anz), do_1lay(anz)

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = MAX(pi_dp / lambda * D_core * m_medium , 1d-14)
    nue = pi_dp / lambda * D * m_medium

    do_1lay = (alf < 1d-7 .AND. nue >= 2d-7)

    do_tmat_tmp = .NOT. do_1lay
    IF (PRESENT(do_tmat)) &
      do_tmat_tmp = (do_tmat_tmp .AND. do_tmat)

    ! tmatrix2 doesn't like if size of entire particle or just the core is
    ! exactly 0.0. 
    ! JM210213: D[_core] > Deps ensured by calling routines
    !
    ! JM190724:
    ! in case of small r_core, wouldn't it be better to call tmatrix
    ! instead right from the start?
    ! particularly since we do the calc per individual size anyways and
    ! call tmatrix below anyways for small alf (that is, where the core
    ! is small?!?), then overwriting what was done here.
    tmerr = -1
    DO i=1,anz
      IF (do_tmat_tmp(i)) THEN
        ! FIX'ED(?): should unify tmatrix2 and tmatrix2_qp regarding:
        ! (a) unit of input params lambda and/or D,D_core
        ! (b) unit of output params f*
!!$        IF (low_prec) THEN
        CALL tmatrix2(lambda,D(i),D_core(i),&
                      aspectratio(i),aspectratio(i),m_shell(i),m_core(i),&
                      fa(i),fb(i),fa0(i),fb0(i),tmerr(i))
!!$        ELSE
!!$          CALL tmatrix2_qp(lambda,D,D_core,&
!!$                           aspectratio,aspectratio,m_shell,m_core,fa,fb,fa0,fb0,tmerr)
!!$        END IF
      END IF
    END DO
    WHERE (tmerr .NE. 0)
      fa  = miss_value
      fb  = miss_value
      fa0 = miss_value
      fb0 = miss_value
    END WHERE

    IF (ANY(do_1lay)) THEN
      do_tmat_tmp = do_1lay
      IF (PRESENT(do_tmat)) &
        do_tmat_tmp = (do_tmat_tmp .AND. do_tmat)

      tmerr = -1
      DO i=1,anz
        IF (do_tmat_tmp(i)) THEN
          ! JM210213
          ! changed tmatrix to make use of (vol equiv) diameter, not radius
          CALL tmatrix(MAX(D(i),Deps),lambda,m_shell(i),aspectratio(i),&
                       fahlp(i),fbhlp(i),fa0hlp(i),fb0hlp(i),tmerr(i))
        END IF
      END DO
      ! JM190725:
      ! Convert from tmatrix sign conventions to current EMVORADO one.
      ! FIXME: we need to check however, whether this is really
      !        consistent with calc_radar_vars orientation conversion
      !        and for other polarimetric moments.
      ! JM191021:
      ! Removed any sign changings in tmatrix(), instead do ALL of them here.
      WHERE (tmerr .EQ. 0)
        fahlp  = CONJG(fahlp)*(-1.d0)   ! fahlp*(-1.d0)
        fbhlp  = CONJG(fbhlp)           ! fbhlp*(-1.d0)
        fa0hlp = CONJG(fa0hlp)*(-1.d0)  ! fa0hlp*(-1.d0)
        fb0hlp = CONJG(fb0hlp)*(-1.d0)  ! !fb0hlp = fb0hlp*(-1.d0)
      ELSEWHERE
        fahlp  = miss_value
        fbhlp  = miss_value
        fa0hlp = miss_value
        fb0hlp = miss_value
      END WHERE

      WHERE(do_tmat_tmp)
        fa  = fahlp
        fb  = fbhlp
        fa0 = fa0hlp
        fb0 = fb0hlp
      END WHERE
    END IF

    RETURN
  END SUBROUTINE COATEDSPHEROID_SCATTER_TMAT_VEC

  SUBROUTINE COATEDSPHEROID_SCATTER_TMAT(&
      D,D_core,m_core,m_shell,m_medium,aspectratio,lambda,&
      fa,fb,fa0,fb0)
    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)    :: D,D_core,aspectratio,lambda,m_medium
    COMPLEX(kind=dp), INTENT(in) :: m_core,m_shell

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa,fb,fa0,fb0

    !..Lokale Variablen
    INTEGER       :: tmerr
    REAL(KIND=dp) :: alf,nue

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = MAX(pi_dp / lambda * D_core * m_medium , 1d-14)
    nue = pi_dp / lambda * D * m_medium

    IF (nue >= 2d-7 .AND. nue >= alf) THEN

      !..Berechung der div. Querschnitte:
      IF (alf >= 1d-7) THEN
        CALL tmatrix2(lambda,D,D_core,&
                      aspectratio,aspectratio,m_shell,m_core,fa,fb,fa0,fb0,tmerr)
        IF (tmerr .NE. 0) THEN
          fa  = miss_value
          fb  = miss_value
          fa0 = miss_value
          fb0 = miss_value
        END IF
      ELSE
        CALL tmatrix(MAX(D,Deps),lambda,m_shell,aspectratio,fa,fb,fa0,fb0,tmerr)
        IF (tmerr .EQ. 0) THEN
          fa  = CONJG(fa)*(-1.d0)  ! fa*(-1.d0)
          fb  = CONJG(fb)          ! fb*(-1.d0)
          fa0 = CONJG(fa0)*(-1.d0) ! fa0*(-1.d0)
          fb0 = CONJG(fb0)*(-1.d0) ! !fb0 = fb0*(-1.d0)
        ELSE
          fa  = miss_value
          fb  = miss_value
          fa0 = miss_value
          fb0 = miss_value
        END IF
      END IF
    ELSE
      fa  = 0.0d0
      fb  = 0.0d0
      fa0 = 0.0d0
      fb0 = 0.0d0
    END IF
    RETURN
  END SUBROUTINE COATEDSPHEROID_SCATTER_TMAT


!----------------------------------------------------------------------------------------!
!
! Alternativ fuer BHCOATBACK:
! Mie-Theorie fuer zweischalige Kugel nach Kerker (1969)
! (etwas langsamer, aber nachvollziehbar)
!
!----------------------------------------------------------------------------------------!

!==================================================================================
!
! Mie-Koeffizienten an und bn der Ordnung n einer Kugel mit
! Brechungsindex m1 (=m_kugel/m_medium), die mit einer Kugelschale
! des Materials mit Brechungsindex m2 (=m_schale/m_medium)
! ueberzogen ist, als Funktion des Mie-Parameters der Kugel x1
! (=2*pi*r1/lambda_medium = 2*pi*r1*m_medium/lambda) und des Mie-Parameters
! von Kugel+Schale x2 (=2*pi*(r1+r2)/lambda_medium) -- beide reell,
! d.h. das aeussere Medium absorbiert nicht: nach Kerker (1969): The scattering
! of light, S. 189 ff.
! Das ganze wird in double-Arithmetik gerechnet!
!
!     ---> hat fuer x1 = 0.0 eine Singularitaet!
!
!==================================================================================

  ! Vektorisierte Form: Es werden die Streukoeffizienten an und bn von der
  ! Ordnung 0 bis zur Ordnung NMAX zurueckgegeben, an(1:NMAX), bn(1:NMAX)
  !
  ! Davon die nochmals ueber verschiedene x-elemente vektorisierte Fassung:
  SUBROUTINE COATEDSPHERE_SCATTERCOEFFS_VEC(anz,NMAX,m1,m2,x1,x2,an,bn,maske)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: anz, NMAX
    COMPLEX(kind=dp), INTENT(in) :: m1(anz), m2(anz)
    REAL(KIND=dp), INTENT(in) :: x1(anz), x2(anz)
    LOGICAL, OPTIONAL, INTENT(in) :: maske(anz)

    COMPLEX(kind=dp), INTENT(out) :: an(anz,1:NMAX), bn(anz,1:NMAX)
    COMPLEX(kind=dp), PARAMETER   :: iimag = (0d0,1d0), czero = (0d0,0d0)

    COMPLEX(kind=dp) :: y1(anz), y2(anz), z2(anz)
    COMPLEX(kind=dp), DIMENSION(1:anz,0:NMAX) :: sphj_y1, sphdj_y1, sphy_y1, sphdy_y1, &
         sphj_y2, sphdj_y2, sphy_y2, sphdy_y2, &
         sphj_z2, sphdj_z2, sphy_z2, sphdy_z2, &
         sphj_x2, sphdj_x2, sphy_x2, sphdy_x2
    COMPLEX(kind=dp) :: psiy1, psiy2, psiz2, &
         dpsiy1, dpsiy2, dpsiz2, &
         xsiy2, xsiz2, dxsiy2, dxsiz2, zetax2, dzetax2
    COMPLEX(kind=dp) :: det1, det2, det3, det4, det5, det6, nenner
    REAL(KIND=dp)     :: psix2, dpsix2
    INTEGER :: NM(anz), i, k
    LOGICAL :: maskeloc(anz)

    IF (ANY(x1 > x2)) THEN
      WRITE (*,*) 'COATEDSPHERE_SCATTERCOEFFS_VEC: Fehler: Innerer Kugelradius > Aeusserer Radius!'
    ENDIF

    IF (.NOT.PRESENT(maske)) THEN
      maskeloc = .TRUE.
    ELSE
      maskeloc = maske
    END IF

    y1 = m1 * x1
    y2 = m2 * x1
    z2 = m2 * x2

    ! Zuerst die weiter unten benoetigten Bessel-Funktionen und deren Ableitungen
    ! bis zur Ordnung NMAX berechnen (backward recurrence):
    CALL CSPHJYVEC(NMAX,y1,NM,sphj_y1,sphdj_y1,sphy_y1,sphdy_y1,anz)
    CALL CSPHJYVEC(NMAX,y2,NM,sphj_y2,sphdj_y2,sphy_y2,sphdy_y2,anz)
    CALL CSPHJYVEC(NMAX,z2,NM,sphj_z2,sphdj_z2,sphy_z2,sphdy_z2,anz)
    CALL CSPHJYVEC(NMAX,CMPLX(x2,0d0,kind=dp),NM,sphj_x2,sphdj_x2,sphy_x2,sphdy_x2,anz)

    an = czero
    bn = czero
    DO i=1,NMAX
      DO k=1,anz
        IF (maskeloc(k)) THEN
          IF (x1(k) > 0d0) THEN
            ! Die Riccati-Bessel-Funktionen der benoetigten Argumente:

            psiy1 = y1(k) * sphj_y1(k,i)
            dpsiy1 = sphj_y1(k,i) + y1(k) * sphdj_y1(k,i)

            psiy2 = y2(k) * sphj_y2(k,i)
            dpsiy2 = sphj_y2(k,i) + y2(k) * sphdj_y2(k,i)

            psiz2 = z2(k) * sphj_z2(k,i)
            dpsiz2 = sphj_z2(k,i) + z2(k) * sphdj_z2(k,i)

            psix2 = x2(k) * DBLE(sphj_x2(k,i))
            dpsix2 = DBLE(sphj_x2(k,i)) + x2(k) * DBLE(sphdj_x2(k,i))

            xsiy2 = -y2(k) * sphy_y2(k,i)
            dxsiy2 = -sphy_y2(k,i) - y2(k) * sphdy_y2(k,i)

            xsiz2 = -z2(k) * sphy_z2(k,i)
            dxsiz2 = -sphy_z2(k,i) - z2(k) * sphdy_z2(k,i)

            zetax2 = x2(k) * (DBLE(sphj_x2(k,i)) - iimag * DBLE(sphy_x2(k,i)))
            dzetax2 = x2(k) * (DBLE(sphdj_x2(k,i)) - iimag * DBLE(sphdy_x2(k,i))) + &
                 (DBLE(sphj_x2(k,i)) - iimag * DBLE(sphy_x2(k,i)))

            ! Einige vorbereitende Terme auf die Streukoeffizienten:

            det1 = m1(k)*dpsiy2*psiy1 - m2(k)*psiy2*dpsiy1
            det2 = dxsiz2*psix2 - m2(k)*xsiz2*dpsix2
            det3 = m1(k)*dxsiy2*psiy1 - m2(k)*xsiy2*dpsiy1
            det4 = dpsiz2*psix2 - m2(k)*psiz2*dpsix2
            det5 = dxsiz2*zetax2 - m2(k)*xsiz2*dzetax2
            det6 = dpsiz2*zetax2 - m2(k)*psiz2*dzetax2

            nenner = det1*det5-det3*det6

            IF (nenner /= czero) THEN
              an(k,i) = (det1*det2-det3*det4)/nenner
            ELSE
              an(k,i) = czero
            END IF

            det1 = m2(k)*dpsiy2*psiy1 - m1(k)*psiy2*dpsiy1
            det2 = m2(k)*dxsiz2*psix2 - xsiz2*dpsix2
            det3 = m2(k)*dxsiy2*psiy1 - m1(k)*xsiy2*dpsiy1
            det4 = m2(k)*dpsiz2*psix2 - psiz2*dpsix2
            det5 = m2(k)*dxsiz2*zetax2 - xsiz2*dzetax2
            det6 = m2(k)*dpsiz2*zetax2 - psiz2*dzetax2

            nenner = det1*det5-det3*det6

            IF (nenner /= czero) THEN
              bn(k,i) = (det1*det2-det3*det4)/nenner
            ELSE
              bn(k,i) = czero
            END IF

          ELSE

            ! Die Riccati-Bessel-Funktionen der benoetigten Argumente:

            psiz2 = z2(k) * sphj_z2(k,i)
            dpsiz2 = sphj_z2(k,i) + z2(k) * sphdj_z2(k,i)

            psix2 = x2(k) * DBLE(sphj_x2(k,i))
            dpsix2 = DBLE(sphj_x2(k,i)) + x2(k) * DBLE(sphdj_x2(k,i))

            zetax2 = x2(k) * (DBLE(sphj_x2(k,i)) - iimag * DBLE(sphy_x2(k,i)))
            dzetax2 = x2(k) * (DBLE(sphdj_x2(k,i)) - iimag * DBLE(sphdy_x2(k,i))) + &
                 (DBLE(sphj_x2(k,i)) - iimag * DBLE(sphy_x2(k,i)))

            det4 = dpsiz2*psix2 - m2(k)*psiz2*dpsix2
            det6 = dpsiz2*zetax2 - m2(k)*psiz2*dzetax2

            IF (det6 /= czero) THEN
              an(k,i) = det4 / det6
            ELSE
              an(k,i) = czero
            END IF

            det4 = m2(k)*dpsiz2*psix2 - psiz2*dpsix2
            det6 = m2(k)*dpsiz2*zetax2 - psiz2*dzetax2

            IF (det6 /= czero) THEN
              bn(k,i) = det4 / det6
            ELSE
              bn(k,i) = czero
            END IF
          END IF
        END IF
      END DO
    END DO

    RETURN

  END SUBROUTINE COATEDSPHERE_SCATTERCOEFFS_VEC

  SUBROUTINE COATEDSPHERE_SCATTER_VEC(D,D_core,m_core,m_shell,m_medium,lambda,&
      fa,fb,fa0,fb0,anz,maske)
    IMPLICIT NONE

    !..Berechnung der Rueckstreueffizienz
    !  r = Radius der Kugel (Gesamtradius, r = r_shell + r_core)
    !  r_shell = Dicke der auesseren Schicht
    !  lambda  = Wellenlaenge


    !..Eingabevariablen
    INTEGER :: anz
    REAL(KIND=dp), INTENT(in)     :: D(anz),D_core(anz),lambda,m_medium
    COMPLEX(kind=dp),INTENT(in)   :: m_core(anz),m_shell(anz)
    LOGICAL, OPTIONAL, INTENT(in) :: maske(anz)

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa(anz),fb(anz),fa0(anz),fb0(anz)

    !..Lokale Variablen
    INTEGER, PARAMETER :: NMAX = 50
    INTEGER          :: n,k,m1expn
    REAL(KIND=dp)    :: alf(anz),nue(anz),nuehlp(anz),&
                        chlp_sca(anz),csum_sca(anz)
    REAL(KIND=dp)    :: C_back(anz) !,C_ext(anz),C_sca(anz)
    COMPLEX(kind=dp) :: an(anz,NMAX),bn(anz,NMAX),&
                        chlp_back(anz),chlp_ext(anz),&
                        csum_back(anz),csum_ext(anz)
    REAL(KIND=dp), PARAMETER :: eps = 1D-6
    LOGICAL :: maskeloc(anz)

    IF (.NOT.PRESENT(maske)) THEN
      maskeloc = .true.
    ELSE
      maskeloc = maske
    END IF

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = pi_dp / lambda * D_core * m_medium
    nue = pi_dp / lambda * D * m_medium

    nuehlp = MAX(nue, xeps)

    !.. Streukoeffizienten bis zur NMAX-ten Ordnung berechnen:
    CALL COATEDSPHERE_SCATTERCOEFFS_VEC(anz,NMAX,m_core/m_medium,m_shell/m_medium,alf,nuehlp,an,bn,maske=maskeloc)

    !..Berechung des Rueckstreuquerschnitts
    !  nach van de Hulst, S. 284
    csum_back = (0d0,0d0)
    csum_ext  = (0d0,0d0)
    csum_sca  = 0d0
    chlp_back = (1d0,0d0)
    chlp_ext  = (1d0,0d0)
    chlp_sca  = 1d0
    m1expn    = 1

    DO n=1,NMAX
      m1expn = m1expn*(-1)
      DO k=1,anz
        IF ( maskeloc(k) .AND. (ABS(chlp_back(k)) >= eps * ABS(csum_back(k)) .OR. &
             ABS(chlp_ext(k)) >= eps * ABS(csum_ext(k)) .OR. &
             ABS(chlp_sca(k)) >= eps * ABS(csum_sca(k))) ) THEN

          chlp_back(k) = (n + 0.5d0) * m1expn * ( an(k,n) - bn(k,n) )
          csum_back(k) = csum_back(k) + chlp_back(k)
          chlp_ext(k) = (n + 0.5d0) * ( an(k,n) + bn(k,n) )
          csum_ext(k) = csum_ext(k) + chlp_ext(k)
          chlp_sca(k) = (n + 0.5d0) * (ABS(an(k,n))*ABS(an(k,n)) + ABS(bn(k,n))*ABS(bn(k,n)) )
          csum_sca(k) = csum_sca(k) + chlp_sca(k)

        END IF
      END DO
    END DO


    WHERE(nue >= xeps .AND. maskeloc)

      C_back = lambda**2 / pi_dp * ABS(csum_back)*ABS(csum_back)
      !C_ext = lambda**2 / pi_dp * DBLE(csum_ext)
      !C_sca = lambda**2 / pi_dp * csum_sca

    ELSEWHERE

      C_back = 0d0
      !C_ext = 0d0
      !C_sca = 0d0

    END WHERE

    !FIXME
    !JM190710: Not really convinced this is correct >:-/
    !  this means fa = fa0 for non-zero fa. odd because:
    !  (a) implies backscatt == extinction or
    !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
    !  both alternatives (a) and (b) seem fishy.
    fa  = sqrt(C_back/(4*pi_dp))
    fa0 = fa
    fb  = fa
    fb0 = fa

  END SUBROUTINE COATEDSPHERE_SCATTER_VEC

  SUBROUTINE COATEDSPHERE_SCATTERCOEFFS(NMAX,m1,m2,x1,x2,an,bn)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: NMAX
    COMPLEX(kind=dp), INTENT(in) :: m1, m2
    REAL(KIND=dp), INTENT(in) :: x1, x2

    COMPLEX(kind=dp), INTENT(out) :: an(1:NMAX), bn(1:NMAX)
    COMPLEX(kind=dp), PARAMETER :: iimag = (0d0,1d0), czero = (0d0,0d0)

    COMPLEX(kind=dp) :: y1, y2, z2
    COMPLEX(kind=dp), DIMENSION(0:NMAX) :: sphj_y1, sphdj_y1, sphy_y1, sphdy_y1, &
         sphj_y2, sphdj_y2, sphy_y2, sphdy_y2, &
         sphj_z2, sphdj_z2, sphy_z2, sphdy_z2, &
         sphj_x2, sphdj_x2, sphy_x2, sphdy_x2
    COMPLEX(kind=dp) :: psiy1, psiy2, psiz2, &
                            dpsiy1, dpsiy2, dpsiz2, &
                            xsiy2, xsiz2, dxsiy2, dxsiz2, zetax2, dzetax2, &
                            det1, det2, det3, det4, det5, det6, nenner
    REAL(KIND=dp) :: psix2, dpsix2
    INTEGER :: NM, i

    IF (x1 > x2) THEN
      WRITE (*,*) 'COATEDSPHERE_SCATTERCOEFFS: Fehler: Innerer Kugelradius > Aeusserer Radius!'
    ENDIF

    y1 = m1 * x1
    y2 = m2 * x1
    z2 = m2 * x2


    IF (x1 > 0d0) THEN

      ! Zuerst die weiter unten benoetigten Bessel-Funktionen und deren Ableitungen
      ! bis zur Ordnung NMAX berechnen (backward recurrence):
      CALL CSPHJY(NMAX,y1,NM,sphj_y1,sphdj_y1,sphy_y1,sphdy_y1)
      CALL CSPHJY(NMAX,y2,NM,sphj_y2,sphdj_y2,sphy_y2,sphdy_y2)
      CALL CSPHJY(NMAX,z2,NM,sphj_z2,sphdj_z2,sphy_z2,sphdy_z2)
      CALL CSPHJY(NMAX,CMPLX(x2,0d0,kind=dp),NM,sphj_x2,sphdj_x2,sphy_x2,sphdy_x2)

      DO i=1,NMAX

        ! Die Riccati-Bessel-Funktionen der benoetigten Argumente:

        psiy1 = y1 * sphj_y1(i)
        dpsiy1 = sphj_y1(i) + y1 * sphdj_y1(i)

        psiy2 = y2 * sphj_y2(i)
        dpsiy2 = sphj_y2(i) + y2 * sphdj_y2(i)

        psiz2 = z2 * sphj_z2(i)
        dpsiz2 = sphj_z2(i) + z2 * sphdj_z2(i)

        psix2 = x2 * DBLE(sphj_x2(i))
        dpsix2 = DBLE(sphj_x2(i)) + x2 * DBLE(sphdj_x2(i))

        xsiy2 = -y2 * sphy_y2(i)
        dxsiy2 = -sphy_y2(i) - y2 * sphdy_y2(i)

        xsiz2 = -z2 * sphy_z2(i)
        dxsiz2 = -sphy_z2(i) - z2 * sphdy_z2(i)

        zetax2 = x2 * (DBLE(sphj_x2(i)) - iimag * DBLE(sphy_x2(i)))
        dzetax2 = x2 * (DBLE(sphdj_x2(i)) - iimag * DBLE(sphdy_x2(i))) + &
                  (DBLE(sphj_x2(i)) - iimag * DBLE(sphy_x2(i)))

        ! Einige vorbereitende Terme auf die Streukoeffizienten:

        det1 = m1*dpsiy2*psiy1 - m2*psiy2*dpsiy1
        det2 = dxsiz2*psix2 - m2*xsiz2*dpsix2
        det3 = m1*dxsiy2*psiy1 - m2*xsiy2*dpsiy1
        det4 = dpsiz2*psix2 - m2*psiz2*dpsix2
        det5 = dxsiz2*zetax2 - m2*xsiz2*dzetax2
        det6 = dpsiz2*zetax2 - m2*psiz2*dzetax2

        nenner = det1*det5-det3*det6

        IF (nenner /= czero) THEN
          an(i) = (det1*det2-det3*det4)/nenner
        ELSE
          an(i) = czero
        END IF

        det1 = m2*dpsiy2*psiy1 - m1*psiy2*dpsiy1
        det2 = m2*dxsiz2*psix2 - xsiz2*dpsix2
        det3 = m2*dxsiy2*psiy1 - m1*xsiy2*dpsiy1
        det4 = m2*dpsiz2*psix2 - psiz2*dpsix2
        det5 = m2*dxsiz2*zetax2 - xsiz2*dzetax2
        det6 = m2*dpsiz2*zetax2 - psiz2*dzetax2

        nenner = det1*det5-det3*det6

        IF (nenner /= czero) THEN
          bn(i) = (det1*det2-det3*det4)/nenner
        ELSE
          bn(i) = czero
        END IF


      END DO

    ELSE

      ! Zuerst die weiter unten benoetigten Bessel-Funktionen und deren Ableitungen
      ! bis zur Ordnung NMAX berechnen (backward recurrence):
      CALL CSPHJY(NMAX,z2,NM,sphj_z2,sphdj_z2,sphy_z2,sphdy_z2)
      CALL CSPHJY(NMAX,CMPLX(x2,0d0,kind=dp),NM,sphj_x2,sphdj_x2,sphy_x2,sphdy_x2)

      DO i=1,NMAX

        ! Die Riccati-Bessel-Funktionen der benoetigten Argumente:

        psiz2 = z2 * sphj_z2(i)
        dpsiz2 = sphj_z2(i) + z2 * sphdj_z2(i)

        psix2 = x2 * DBLE(sphj_x2(i))
        dpsix2 = DBLE(sphj_x2(i)) + x2 * DBLE(sphdj_x2(i))

        zetax2 = x2 * (DBLE(sphj_x2(i)) - iimag * DBLE(sphy_x2(i)))
        dzetax2 = x2 * (DBLE(sphdj_x2(i)) - iimag * DBLE(sphdy_x2(i))) + &
                 (DBLE(sphj_x2(i)) - iimag * DBLE(sphy_x2(i)))

        det4 = dpsiz2*psix2 - m2*psiz2*dpsix2
        det6 = dpsiz2*zetax2 - m2*psiz2*dzetax2

        IF (det6 /= czero) THEN
          an(i) = det4 / det6
        ELSE
          an(i) = czero
        END IF

        det4 = m2*dpsiz2*psix2 - psiz2*dpsix2
        det6 = m2*dpsiz2*zetax2 - psiz2*dzetax2

        IF (det6 /= czero) THEN
          bn(i) = det4 / det6
        ELSE
          bn(i) = czero
        END IF

      END DO

    END IF

    RETURN

  END SUBROUTINE COATEDSPHERE_SCATTERCOEFFS

  SUBROUTINE COATEDSPHERE_SCATTER(D,D_core,m_core,m_shell,m_medium,lambda,&
      fa,fb,fa0,fb0)
    IMPLICIT NONE

    !..Berechnung der Rueckstreueffizienz
    !  r = Radius der Kugel (Gesamtradius, r = r_shell + r_core)
    !  r_shell = Dicke der auesseren Schicht
    !  lambda  = Wellenlaenge


    !..Eingabevariablen
    REAL(KIND=dp), INTENT(in)   :: D,D_core,lambda,m_medium
    COMPLEX(kind=dp),INTENT(in) :: m_core,m_shell

    !..Ausgabevariablen
    COMPLEX(kind=dp), INTENT(out) :: fa,fb,fa0,fb0

    !..Lokale Variablen
    INTEGER, PARAMETER :: NMAX = 50
    INTEGER          :: n,m1expn
    REAL(KIND=dp)    :: alf,nue,chlp_sca,csum_sca,dbln
    REAL(KIND=dp)    :: C_back !,C_ext,C_sca
    COMPLEX(kind=dp) :: an(NMAX),bn(NMAX),&
                        chlp_back,chlp_ext,csum_back,csum_ext
    REAL(KIND=dp), PARAMETER :: eps = 1D-6

    !..Dimensionslose Radien (vgl. Kerker, S. 192, Gl. (5.1.25)
    alf = pi_dp / lambda * D_core * m_medium
    nue = pi_dp / lambda * D * m_medium

    IF (nue >= xeps) THEN

      !.. Streukoeffizienten bis zur NMAX-ten Ordnung berechnen:
      CALL COATEDSPHERE_SCATTERCOEFFS(NMAX,m_core/m_medium,m_shell/m_medium,alf,nue,an,bn)

      !..Berechung des Rueckstreuquerschnitts
      !  nach van de Hulst, S. 284
      csum_back = (0d0,0d0)
      csum_ext = (0d0,0d0)
      csum_sca = 0d0
      chlp_back = (1d0,0d0)
      chlp_ext = (1d0,0d0)
      chlp_sca = 1d0
      m1expn   = 1

      DO n=1,NMAX

        IF ( ABS(chlp_back) >= eps * ABS(csum_back) .OR. &
             ABS(chlp_ext) >= eps * ABS(csum_ext) .OR. &
             ABS(chlp_sca) >= eps * ABS(csum_sca) ) THEN

          dbln = DBLE(n)
          m1expn = m1expn*(-1)
          chlp_back = (dbln + 0.5d0) * m1expn * ( an(n) - bn(n) )
          csum_back = csum_back + chlp_back
          chlp_ext = (dbln + 0.5d0) * ( an(n) + bn(n) )
          csum_ext = csum_ext + chlp_ext
          chlp_sca = (dbln + 0.5d0) * (ABS(an(n))*ABS(an(n)) + ABS(bn(n))*ABS(bn(n)) )
          csum_sca = csum_sca + chlp_sca

        ELSE

          EXIT

        END IF

      ENDDO

      C_back = lambda**2 / pi_dp * ABS(csum_back)*ABS(csum_back)
      !C_ext = lambda**2 / pi_dp * DBLE(csum_ext)
      !C_sca = lambda**2 / pi_dp * csum_sca

    ELSE

      C_back = 0d0
      !C_ext = 0d0
      !C_sca = 0d0

   END IF

    !FIXME
    !JM190710: Not really convinced this is correct >:-/
    !  this means fa = fa0 for non-zero fa. odd because:
    !  (a) implies backscatt == extinction or
    !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
    !  both alternatives (a) and (b) seem fishy.
    fa  = sqrt(C_back/(4*pi_dp))
    fa0 = fa
    fb  = fa
    fb0 = fa

  END SUBROUTINE COATEDSPHERE_SCATTER



!*********************************************************************************************

  !----------------------------------------------------------------------------------------!
  ! 4) Routinen zur Berechnung von Q_back, Q_ext und Q_sca von verschiedenen
  !    einzelnen Hydrometeoren
  !----------------------------------------------------------------------------------------!


  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für nassen Hagel: Eiskugel mit Wasserhaut
  !
  ! Schwachpunkt: Shedding wird nicht beachtet --- die Wasserhaut ist in ihrer Dicke
  !               somit nicht begrenzt. (noch einbauen!!)
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_WATERSPH_WETHAIL_VEC(D_h,a_geo,b_geo,fmelt,&
      m_i,m_w,aspectratio,lambda,&
      fa,fb,fa0,fb0,anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_h(anz), fmelt(anz), &
                                    a_geo, b_geo, lambda, &
                                    aspectratio(anz)
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)

    COMPLEX(kind=dp), INTENT(out) :: fa(anz), fb(anz), fa0(anz), fb0(anz)

    REAL(KIND=dp), DIMENSION(anz) :: x_h, x_w, x_core, x_shell, &
                                     D_large, D_core, &
                                     rho_core, rho_shell, &
                                     fm
    COMPLEX(kind=dp)              :: m_core(anz), m_shell(anz)

    fm = MAX(MIN(fmelt, 1d0), 0d0)
    x_h = a_geo * D_h**b_geo
    x_w = x_h * fm

    m_core = m_i
    m_shell = m_w

    rho_core = rho_ice_fwo
    rho_shell = rho_w_fwo

    x_shell = x_w
    x_core  = x_h-x_w

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-10)
      D_large = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    IF (luse_tmatrix) THEN
      !..T-Matrix
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-scattering
      CALL COATEDSPHERE_SCATTER_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                    m_core,m_shell,DBLE(m_air),lambda,&
                                    fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_large < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_WATERSPH_WETHAIL_VEC

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_WATERSPH_WETHAIL_BH_VEC(D_h,a_geo,b_geo,fmelt,&
      m_i,m_w,aspectratio,lambda,&
      fa,fb,fa0,fb0,anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_h(anz), fmelt(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)

    COMPLEX(kind=dp), INTENT(out) :: fa(anz), fb(anz), fa0(anz), fb0(anz)

    REAL(KIND=dp), DIMENSION(anz) :: x_h, x_w, x_core, x_shell, &
                                     D_large, D_core, &
                                     rho_core, rho_shell, &
                                     fm
    COMPLEX(kind=dp)              :: m_core(anz), m_shell(anz)

    fm = MAX(MIN(fmelt, 1d0), 0d0)
    x_h = a_geo * D_h**b_geo
    x_w = x_h * fm

    m_core = m_i
    m_shell = m_w

    rho_core = rho_ice_fwo
    rho_shell = rho_w_fwo

    x_shell = x_w
    x_core  = x_h-x_w

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-10)
      D_large = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    IF (luse_tmatrix) THEN
      !..T-Matrix
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                       m_core,m_shell,DBLE(m_air),lambda,&
                                       fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_large < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_WATERSPH_WETHAIL_BH_VEC

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für spongy Hagel (Eiskern + Wasser/Eis-Schale mit Schmelzgrad fmelt;
  ! Das Schmelzwasser verteilt sich auf eine Schale, die ihrerseits einen Wasser-Volumenanteil
  ! von f_water hat. f_water darf nicht 0.0 oder kleiner sein!!!
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_SPONGY_WETHAIL_VEC(D_h,a_geo,b_geo,fmelt,f_water,&
      m_i,m_w,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_h(anz), fmelt(anz), f_water(anz), &
                                    aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    LOGICAL :: maske(anz)
    REAL(KIND=dp), DIMENSION(anz) :: volair, volice, volwater, volh, &
                                     rho_core, rho_shell, &
                                     x_h, x_w, x_core, x_shell, &
                                     D_large, D_core, &
                                     f_mass, f_ice, fwater, &
                                     fm
    COMPLEX(kind=dp)              :: m_core(anz), m_shell(anz), m_airv(anz)
    COMPLEX(kind=dp), DIMENSION(anz) :: fal, fbl, fa0l, fb0l

    fehler = 0

    m_airv = m_air

    !..Total mass of melted water:
    fm = MAX(MIN(fmelt, 1d0), 0d0)
    x_h = a_geo * D_h**b_geo
    x_w = x_h * fm

    !..Volume fraction of water in spongy shell:
    fwater = MAX(MIN(f_water, 1d0), 0d0)
    !..Volume fraction of ice in spongy shell:
    f_ice = 1d0 - fwater

    m_core = m_i

    rho_shell = rho_w_fwo * fwater + rho_ice_fwo * f_ice
    rho_core  = rho_ice_fwo

    f_mass = rho_w_fwo / rho_shell * fwater

    ! Mass of shell: distinction of cases needed, since in case of f_water < 1
    ! already for melting degree < 1 it happens that the complete particle
    ! consists of an ice-water-mixture such that the core completely disappears
    ! and the complete particle is shell. With further melting the relative
    ! volume fraction of water, f_water, changes.

    x_shell = x_w / f_mass

    !.. 1) Case x_shell < x_h
    maske = (x_shell < x_h)
    WHERE(maske)
      x_core  = x_h - x_shell
    ELSEWHERE
      x_core  = 0d0
    END WHERE

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-10)
      D_large = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    !..Complex refractive index of water-ice-mixture of shell:
    ! Volume fraction of water in ice-water-mixture
    volwater = fwater
    volice = 1d0 - volwater
    volair = 0d0

    m_shell = get_m_mix_vec(m_airv, m_i, m_w, volair, volice, volwater, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) RETURN

    maske = (x_shell < x_h .AND. D_large >= Deps)
    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                    m_core,m_shell,DBLE(m_air),lambda,&
                                    fa,fb,fa0,fb0,anz,maske=maske)
    ENDIF

    WHERE (.NOT.maske)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

    !.. 2) Case x_shell = x_h

    volh = (x_h - x_w) / rho_ice_fwo + x_w / rho_w_fwo

    D_large  = (inv_pi6 * volh) ** third

    volice = (x_h - x_w) / (volh * rho_ice_fwo)
    volwater = 1d0 - volice
    volair = 0d0

    m_core = get_m_mix_vec(m_airv, m_i, m_w, volair, volice, volwater, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) RETURN

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),&
                                     m_core,DBLE(m_air),aspectratio,lambda,&
                                     fal,fbl,fa0l,fb0l,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),&
                                 m_core,DBLE(m_air),lambda,&
                                 fal,fbl,fa0l,fb0l,anz)
    ENDIF

    WHERE (D_large >= Deps .AND. x_shell >= x_h)
      fa  = fal
      fb  = fbl
      fa0 = fa0l
      fb0 = fb0l
    END WHERE

  END SUBROUTINE MIE_SPONGY_WETHAIL_VEC

  ! Version with calculation of scattering function according to recursion formulas
  ! by Bohren und Huffman (1983):

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_SPONGY_WETHAIL_BH_VEC(D_h,a_geo,b_geo,fmelt,f_water,&
      m_i,m_w,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_h(anz), fmelt(anz), f_water(anz), &
                                    aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    LOGICAL :: do_tmat(anz)
    REAL(KIND=dp), DIMENSION(anz) :: volair, volice, volwater, volh, &
                                     rho_core, rho_shell, &
                                     x_h, x_w, x_core, x_shell, &
                                     D_large_wcore, D_large_ncore, D_core, &
                                     f_mass, f_ice, fwater, &
                                     fm
    COMPLEX(kind=dp)              :: m_core(anz), m_shell(anz), m_airv(anz)
    COMPLEX(kind=dp), DIMENSION(anz) :: fa1, fb1, fa01, fb01

    fehler = 0

    m_airv = m_air

    !..Total mass of melted water:
    fm = MAX(MIN(fmelt, 1d0), 0d0)
    x_h = a_geo * D_h**b_geo
    x_w = x_h * fm

    !..Volume fraction of water in spongy shell:
    fwater = MAX(MIN(f_water, 1d0), 0d0)
    !..Volume fraction of ice in spongy shell:
    f_ice = 1d0 - fwater

    m_core = m_i

    rho_shell = rho_w_fwo * fwater + rho_ice_fwo * f_ice
    rho_core  = rho_ice_fwo

    f_mass = rho_w_fwo / rho_shell * fwater

    ! Mass of shell: distinction of cases needed, since in case of f_water < 1
    ! already for melting degree < 1 it happens that the complete particle
    ! consists of an ice-water-mixture such that the core completely disappears
    ! and the complete particle is shell. With further melting the relative
    ! volume fraction of water, f_water, changes.

    x_shell = x_w / f_mass

    !.. 1.0) Prep for case x_shell < x_h

    WHERE (x_shell < x_h)
      x_core = x_h - x_shell
    ELSEWHERE
      x_core = 0d0
    END WHERE

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-10)
      D_large_wcore = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large_wcore = D_core
    END WHERE

    !.. 2.0) Prep for case x_shell = x_h

    volh = (x_h - x_w) / rho_ice_fwo + x_w / rho_w_fwo
    D_large_ncore  = (inv_pi6 * volh) ** third


    !.. 1) Case x_shell < x_h

    do_tmat = .TRUE.
    WHERE (x_shell >= x_h .OR. D_large_wcore < Deps)
      do_tmat = .FALSE.
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    ELSEWHERE (D_large_ncore >= Deps .AND. x_shell >= x_h)
      ! JM200719: these will be handled in case 2)
      do_tmat = .FALSE.
    END WHERE

    !..Complex refractive index of water-ice-mixture of the shell:
    ! Volume fraction of water in ice-water-mixture
    volwater = fwater
    volice   = 1d0 - volwater
    volair   = 0d0

    m_shell = get_m_mix_vec(m_airv, m_i, m_w, volair, volice, volwater, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      WHERE (do_tmat)
        fa  = CMPLX(0d0,0d0,kind=dp)
        fb  = CMPLX(0d0,0d0,kind=dp)
        fa0 = CMPLX(0d0,0d0,kind=dp)
        fb0 = CMPLX(0d0,0d0,kind=dp)
      END WHERE
    ELSE

      IF (luse_tmatrix) THEN
        CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large_wcore,Deps),MAX(D_core,Deps),&
                                             m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                             fa,fb,fa0,fb0,anz,do_tmat)
      ELSE
        !..Mie-Streuung
        CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_large_wcore,Deps),MAX(D_core,Deps),&
                                         m_core,m_shell,DBLE(m_air),lambda,&
                                         fa,fb,fa0,fb0,anz)
      END IF
    END IF


    !.. 2) Case x_shell = x_h

    do_tmat = .FALSE.
    WHERE (D_large_ncore >= Deps .AND. x_shell >= x_h)
      do_tmat = .TRUE.
    ENDWHERE

    volice   = (x_h - x_w) / (volh * rho_ice_fwo)
    volwater = 1d0 - volice
    volair   = 0d0

    m_core = get_m_mix_vec(m_airv, m_i, m_w, volair, volice, volwater, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0)  THEN
      WHERE (do_tmat)
        fa  = CMPLX(0d0,0d0,kind=dp)
        fb  = CMPLX(0d0,0d0,kind=dp)
        fa0 = CMPLX(0d0,0d0,kind=dp)
        fb0 = CMPLX(0d0,0d0,kind=dp)
      END WHERE
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_large_ncore,Deps),&
                                     m_core,DBLE(m_air),aspectratio,lambda,&
                                     fa1,fb1,fa01,fb01,anz,do_tmat)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_large_ncore,Deps),&
                                 m_core,DBLE(m_air),lambda,&
                                 fa1,fb1,fa01,fb01,anz)
    ENDIF

    WHERE (do_tmat)
        fa  = fa1
        fb  = fb1
        fa0 = fa01
        fb0 = fb01
    END WHERE

  END SUBROUTINE MIE_SPONGY_WETHAIL_BH_VEC

  !----------------------------------------------------------------------------------------!
  ! Backscattering cross section for dry hail: ice sphere
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_DRYHAIL_VEC(D_h,m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_h(anz),aspectratio(anz),lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz)

    COMPLEX(kind=dp), INTENT(out) :: fa(anz),fb(anz),fa0(anz),fb0(anz)

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_h,Deps),&
                                     m_i,DBLE(m_air),aspectratio,lambda,&
                                     fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_h,Deps),&
                                 m_i,DBLE(m_air),lambda,&
                                 fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_h < 1d-8)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_DRYHAIL_VEC

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für nassen Graupel/Hagel: Zweischalenmodell, wobei sich
  ! das am Aussenrand des Eis/Luft-Kernes schmelzende Wasser am Aussenrand des Partikels
  ! als Wasserhaut ansammelt. Der Brechungsindex des Kerns wird nach Mischungsregeln
  ! berechnet, die durch die Schluesselworte mixingrulestring, matrixstring, inclusionstring
  ! definiert werden.
  !
  ! Schwachpunkt: Shedding wird nicht beachtet --- die Wasserhaut ist in ihrer Dicke
  !               somit nicht begrenzt. (noch einbauen!!)
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_WATERSPH_WETGR_VEC(D_g,a_geo,b_geo,fmelt,&
      m_w,m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring, matrixstring, inclusionstring,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz) :: x_g, x_w, x_core, x_shell, D_large, D_core, &
                                     volice, volair, &
                                     rho_core, rho_shell, &
                                     f_core, fzero, fm
    COMPLEX(kind=dp), DIMENSION(anz) :: m_core, m_airv

    m_airv = m_air
    fzero = 0d0

    fm = MAX(MIN(fmelt, 1d0), 0d0)

    x_g = a_geo * D_g**b_geo
    x_w = x_g * fm
    rho_shell = rho_w_fwo

    rho_core  = MAX(MIN(x_g * inv_pi6 / (D_g*D_g*D_g), rho_ice_fwo), 10d0)

    f_core = rho_core/rho_ice_fwo

    x_shell = x_w
    x_core  = x_g - x_w

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-6)
      D_large = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    !..Complex refractive index of air-ice-mixture of core:
    m_core = get_m_mix_vec(m_airv, m_i, m_w, (1d0-f_core), f_core, fzero, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_w,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                    m_core,m_w,DBLE(m_air),lambda,&
                                    fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_g < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_WATERSPH_WETGR_VEC

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_WATERSPH_WETGR_BH_VEC(D_g,a_geo,b_geo,fmelt,&
      m_w,m_i,aspectratio,ardry,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring, matrixstring, inclusionstring,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), ardry(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz) :: x_g, x_w, x_core, x_shell, D_large, D_core, &
                                     volice, volair, &
                                     v_core, rho_core, rho_shell, &
                                     f_core, fzero, fm
    COMPLEX(kind=dp), DIMENSION(anz) :: m_core, m_airv

    m_airv = m_air
    fzero = 0d0

    fm = MAX(MIN(fmelt, 1d0), 0d0)

    x_g = a_geo * D_g**b_geo
    x_w = x_g * fm
    rho_shell = rho_w_fwo

    v_core    = pi6 * D_g**3 * ardry
    rho_core  = MAX(MIN(x_g / v_core, rho_ice_fwo), 10d0)

    f_core = rho_core/rho_ice_fwo

    x_shell = x_w
    x_core  = x_g - x_w

    D_core  = (inv_pi6 * x_core/rho_core)**third
    WHERE (fm >= 1d-6)
      D_large = (inv_pi6 * x_shell/rho_shell + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    !..Complex refractive index of air-ice-mixture of core:
    m_core = get_m_mix_vec(m_airv, m_i, m_w, (1d0-f_core), f_core, fzero, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_w,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                       m_core,m_w,DBLE(m_air),lambda,&
                                       fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_g < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_WATERSPH_WETGR_BH_VEC

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für nassen Graupel/Hagel: Zweischalenmodell, wobei das
  ! aussen schmelzende Wasser sukzessive die Poren in der aeusseren Kugelschale vollstaendig
  ! fuellt. Die aeussere Schale besteht somit aus Eis und Wasser, waehrend die
  ! innere Schale aus Eis und Luft besteht. Solange nicht alle Hohlraeume mit Schmelzwasser
  ! gefuellt sind, schmilzt das Partikel nur von aussen, d.h. das verbleibende Eisskellet
  ! im inneren schmilzt zunaechst nicht mit. Erst wenn alle Hohlraeume voll sind, schmilzt
  ! auch das Eisskellet.
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_SOAK_TWOSPH_WETGR_VEC(D_g,a_geo,b_geo,fmelt,&
      m_w,m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,&
      mixingrulestring_icewater,matrixstring_icewater,inclusionstring_icewater,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring_iceair, matrixstring_iceair, &
                                    inclusionstring_iceair, &
                                    mixingrulestring_icewater, matrixstring_icewater, &
                                    inclusionstring_icewater

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    LOGICAL :: maske(anz)
    REAL(KIND=dp), DIMENSION(anz) :: x_g, x_w, D_large, D_core, D_c3, &
                                     rho_core, rho_shell, &
                                     f_ice, f_water, fzero, fm
    COMPLEX(kind=dp), DIMENSION(anz) :: m_core, m_shell, m_airv, &
                                        fal, fbl, fa0l, fb0l

    fzero = 0d0
    m_airv = m_air

    fm = MAX(MIN(fmelt, 1d0), 0d0)

    x_g = a_geo * D_g**b_geo
    x_w = x_g * fm

    rho_core = MAX(MIN(x_g * inv_pi6 / (D_g*D_g*D_g), rho_ice_fwo), 10d0)
    f_ice    = rho_core / rho_ice_fwo

    f_water   = 1d0 - f_ice
    rho_shell = rho_ice_fwo * f_ice + rho_w_fwo * f_water

    D_c3 = inv_pi6 * x_g/rho_core * (1d0 - (rho_shell*x_w) / (f_water*rho_w_fwo*x_g))

    !.. 1) Case D_c3 >= 0d0:
    !      (Not all cavities filled yet with water. Partial density of ice remains
    !       constant, remaining ice skeleton is not yet melting.)

    ! UB: D_c3 could be negative, so we need to care about it:
    WHERE (D_c3 >= 0d0)
      D_core  = D_c3**third
    ELSEWHERE
      D_core = 0d0
    END WHERE

    WHERE (fm >= 1d-6)
      ! UB: we have to use D_core**3 instead of D_c3, because the latter
      !     can be negative and crash the "**third"
      D_large = (inv_pi6 * x_w/rho_w_fwo/f_water + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    !..Complex refractive index of air-ice-mixture of core:
    m_core = get_m_mix_vec(m_airv, m_i, m_w, (1d0-f_ice), f_ice, fzero, &
         mixingrulestring_iceair, matrixstring_iceair, inclusionstring_iceair, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    !..Complex refractive index of water-ice-mixture of shell:
    m_shell = get_m_mix_vec(m_airv, m_i, m_w, fzero, f_ice, (1d0-f_ice), &
         mixingrulestring_icewater, matrixstring_icewater, inclusionstring_icewater, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    maske = (D_c3 >= 0d0 .AND. D_g >= Deps)
    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                    m_core,m_shell,DBLE(m_air),lambda,&
                                    fa,fb,fa0,fb0,anz,maske=maske)
    ENDIF

    WHERE (.NOT.maske)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

    !.. 2) Case D_c3 < 0d0:
    !      All cavities water filled, ice skeleton starts to melt:

    f_ice   = (1d0 - fm) / (rho_ice_fwo * (fm/rho_w_fwo + (1d0-fm)/rho_ice_fwo))
    f_water = 1d0 - f_ice

    rho_shell = rho_ice_fwo * f_ice + rho_w_fwo * f_water

    D_large = (inv_pi6 * x_g/rho_shell)**third

    !..Complex refractive index of water-ice-mixture of shell:
    m_shell = get_m_mix_vec(m_airv, m_i, m_w, fzero, f_ice, f_water, &
         mixingrulestring_icewater, matrixstring_icewater, inclusionstring_icewater, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),&
                                     m_shell,DBLE(m_air),aspectratio,lambda,&
                                     fal,fbl,fa0l,fb0l,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),&
                                 m_shell,DBLE(m_air),lambda,&
                                 fal,fbl,fa0l,fb0l,anz)
    ENDIF

    WHERE (D_large >= Deps .AND. D_c3 < 0d0)
      fa  = fal
      fb  = fbl
      fa0 = fa0l
      fb0 = fb0l
    END WHERE

  END SUBROUTINE MIE_SOAK_TWOSPH_WETGR_VEC

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_SOAK_TWOSPH_WETGR_BH_VEC(D_g,a_geo,b_geo,fmelt,&
      m_w,m_i,aspectratio,ardry,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,&
      mixingrulestring_icewater,matrixstring_icewater,inclusionstring_icewater,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), ardry(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring_iceair, matrixstring_iceair, &
                                    inclusionstring_iceair, &
                                    mixingrulestring_icewater, matrixstring_icewater, &
                                    inclusionstring_icewater

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz) :: x_g, x_w, D_large, D_core, D_c3, &
                                     v_core, rho_core, rho_shell, &
                                     f_ice, f_water, fzero, fm
    COMPLEX(kind=dp), DIMENSION(anz) :: m_core, m_shell, m_airv, &
                                        fal, fbl, fa0l, fb0l

    fzero = 0d0
    m_airv = m_air

    fm = MAX(MIN(fmelt, 1d0), 0d0)

    x_g = a_geo * D_g**b_geo
    x_w = x_g * fm

    v_core   = pi6 * D_g**3 * ardry
    rho_core = MAX(MIN(x_g / v_core, rho_ice_fwo), 10d0)
    f_ice    = rho_core / rho_ice_fwo

    f_water   = 1d0 - f_ice
    rho_shell = rho_ice_fwo * f_ice + rho_w_fwo * f_water

    D_c3 = inv_pi6 * x_g/rho_core * (1d0 - (rho_shell*x_w) / (f_water*rho_w_fwo*x_g))

    !.. 1) Case D_c3 >= 0d0:
    !      (Not all cavities filled yet with water. Partial density of ice remains
    !       constant, remaining ice skeleton is not yet melting.)

    WHERE (D_c3 >= 0d0)
      D_core  = D_c3**third
    ELSEWHERE
      D_core = 0d0
    END WHERE

    WHERE (fm >= 1d-6)
      ! UB: we have to use D_core**3 instead of D_c3, because D_c3 can be negative and crash the "**third"
      D_large = (inv_pi6 * x_w/rho_w_fwo/f_water + D_core**3)**third
    ELSEWHERE
      D_large = D_core
    END WHERE

    !..Complex refractive index of air-ice-mixture of core:
    m_core = get_m_mix_vec(m_airv,m_i,m_w,f_water,f_ice,fzero,&
                           mixingrulestring_iceair,matrixstring_iceair,&
                           inclusionstring_iceair,&
                           anz,fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    !..Complex refractive index of water-ice-mixture of shell:
    m_shell = get_m_mix_vec(m_airv,m_i,m_w,fzero,f_ice,f_water,&
                            mixingrulestring_icewater,matrixstring_icewater,&
                            inclusionstring_icewater,&
                            anz,fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                           m_core,m_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),MAX(D_core,Deps),&
                                       m_core,m_shell,DBLE(m_air),lambda,&
                                       fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_c3 < 0.0d0 .OR. D_large < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

    !.. 2) Case D_c3 < 0d0:
    !      All cavities water filled, ice skeleton starts to melt:

    f_ice   = (1d0 - fm) / (rho_ice_fwo * (fm/rho_w_fwo + (1d0-fm)/rho_ice_fwo))
    f_water = 1d0 - f_ice

    rho_shell = rho_ice_fwo * f_ice + rho_w_fwo * f_water

    D_large = (inv_pi6 * x_g/rho_shell)**third

    !..Complex refractive index of water-ice-mixture of shell:
    m_shell = get_m_mix_vec(m_airv,m_i,m_w,fzero,f_ice,f_water,&
                            mixingrulestring_icewater,matrixstring_icewater,&
                            inclusionstring_icewater,&
                            anz,fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_large,Deps),&
                                     m_shell,DBLE(m_air),aspectratio,lambda,&
                                     fal,fbl,fa0l,fb0l,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_large,Deps),&
                                 m_shell,DBLE(m_air),lambda,&
                                 fal,fbl,fa0l,fb0l,anz)
    ENDIF

    WHERE (D_large >= Deps .AND. D_c3 < 0d0 )
      fa  = fal
      fb  = fbl
      fa0 = fa0l
      fb0 = fb0l
    END WHERE

  END SUBROUTINE MIE_SOAK_TWOSPH_WETGR_BH_VEC

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für nassen Graupel/Hagel
  ! als massengewichtetes Mittel aus MIE_SOAK_TWOSPHERE_WETGRAUPEL und MIE_WATERSPHERE_WETGRAUPEL
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_MEAN_WETGR_VEC(D_g,a_geo,b_geo,fmelt,&
      m_w,m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,&
      mixingrulestring_icewater,matrixstring_icewater,inclusionstring_icewater,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring_iceair, matrixstring_iceair, &
                                    inclusionstring_iceair, &
                                    mixingrulestring_icewater, matrixstring_icewater, &
                                    inclusionstring_icewater

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    REAL(KIND=dp), DIMENSION(anz) :: C_soak, C_watersphere, &
                                     C_soak_ext, C_watersphere_ext, &
                                     fm
    COMPLEX(kind=dp), DIMENSION(anz) :: fa_soak, fb_soak, &
                                        fa0_soak, fb0_soak, &
                                        fa_watersphere, fb_watersphere, &
                                        fa0_watersphere, fb0_watersphere

    fm = MAX(MIN(fmelt, 1d0), 0d0)

    CALL MIE_SOAK_TWOSPH_WETGR_VEC(D_g,a_geo,b_geo,fmelt,&
         m_w,m_i,aspectratio,lambda,&
         fa_soak,fb_soak,fa0_soak,fb0_soak,&
         mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,&
         mixingrulestring_icewater,matrixstring_icewater,inclusionstring_icewater,&
         anz,luse_tmatrix)
    CALL MIE_WATERSPH_WETGR_VEC(D_g,a_geo,b_geo,fmelt,&
         m_w,m_i,aspectratio,lambda,&
         fa_watersphere,fb_watersphere,fa0_watersphere,fb0_watersphere,&
         mixingrulestring_iceair,matrixstring_iceair,inclusionstring_iceair,&
         anz,luse_tmatrix)

  ! JCS -- NEED TO COME BACK AND ADDRESS THIS!!!!!
    fa  = fm * fa_watersphere  + (1d0 - fm) * fa_soak
    fb  = fm * fb_watersphere  + (1d0 - fm) * fb_soak
    fa0 = fm * fa0_watersphere + (1d0 - fm) * fa0_soak
    fb0 = fm * fb0_watersphere + (1d0 - fm) * fb0_soak

  END SUBROUTINE MIE_MEAN_WETGR_VEC


  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für nassen Graupel/Hagel/Schnee: Einschalenmodell, wobei das
  ! aussen schmelzende Wasser sukzessive homogen die Poren in der ganzen Kugel
  ! fuellt. Das Eisskelett schmilzt dabei nur von aussen, so dass die Partialdichte
  ! des Eises zunaechst konstant bleibt. Dies bildet das durch Kapillarkraefte
  ! verursachte Einsaugen von Schmelzwasser in die zunaechst intakt bleibende
  ! Eisstruktur nach, wie es auch in Schmelzexperimenten bei Graupel geringerer bis
  ! mittlerer Dichte und Schneeflocken beobachtet wird.
  ! Erst wenn alle Hohlraeume mit Schmelzwasser
  ! gefuellt sind, dann schmilzt auch das nun disintegrierende Eisskelett.
  ! In diesem Stadium wird angenommen, dass die Reste des Eisskelettes homogen
  ! im resultierenden Tropfen verteilt umherschwimmen.
  !
  ! Als Alternative kann man mittels des Schalters schmelzstring auf ein gaenzlich
  ! anderes Schmelzverhalten umschalten: Bei schmelzstring = 'volumentreu' wird
  ! im Gegensatz zu oben (schmelzstring = 'volumenverkleinernd') angenommen,
  ! dass das Eisskelett im Inneren schmilzt und
  ! somit das Partikelvolumen waehrend des kompletten Schmelzvorganges konstant bleibt.
  ! Erst im allerletzten Moment kollabiert das Teilchen zu Wassertropfen.
  ! Vor allem im letzten Schmelzstadium ist dies jedoch hoechst unrealistisch.
  !
  ! NEU: Bisher nur Schmelzen entweder vollstaendig an der Aussenseite oder vollstaendig
  !      im Inneren implementiert. Nun kann mittels des Verhaeltnisses
  !      0 <= meltratio_outside <= 1 derjenige Anteil der geschmolzenen Masse
  !      vorgegeben werden, der von der Aussenseite geschmolzen und in das
  !      Partikel aufgesogen worden ist. (1 - meltratio_outside)
  !      ist dann der Anteil der im Inneren geschmolzenen Masse an der gesamten
  !      Schmelzmasse.
  !
  ! Das Schluesselwort schmelzstring daraufhin wieder entfernt!!!
  ! "Volumentreues" Schmelzen kann nun mittels meltratio_outside = 0.0 erreicht werden,
  ! und "volumenverkleinerndes" Schmelzen mittels meltratio_outside = 1.0
  !
  ! Wenn man unterstellt, dass die Partikel das Schmelzwasser durch Kapillarkraefte
  ! vollstaendig aufsaugen und homogen im Inneren verteilen, dann duerfte die "Wahrheit"
  ! in etwa zwischen den obigen beiden Annahmen zu suchen sein. Man koennte also
  ! als Mittelweg festlegen, dass ein bestimmter Anteil des Schmelzwassers an der Aussenseite und
  ! der andere Teil im Inneren schmilzt, wobei dieser Anteil entweder fest vorgegeben
  ! oder abhaengig vom Schmelzgrad definiert werden koennte. Dies wuerde aber wiederum
  ! noch mehr Freiheitsgrade einfuehren und ist noch nicht implementiert.
  !
  ! Doch, jetzt implementiert: meltratio_outside, Werte zwischen 0 und 1.
  !
  ! Wenn meltratio_outside < 1, dann kann es vorkommen (je nach Dichte des Teilchens), dass
  ! das Eisskelett geschmolzen ist, bevor noch alle Hohlraeume mit Schmelzwasser gefuellt
  ! worden sind. Dies wuerde ein instabiles und sehr schnell platzendes Teilchen bedeuten.
  ! Um dies zu vermeiden, wird meltratio_outside als Grenzwert bei fmelt = 0 angenommen
  ! und waechst linear mit fmelt an, bis es fuer fmelt = 1 ebenfalls den Wert 1 annimmt.
  ! Auf diese Weise geht das Schmelzpartikel asymptotisch gegen einen Wassertropfen.
  !
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_SOAK_WETGR_VEC(&
      D_g,a_geo,b_geo,fmelt,meltratio_outside,&
      m_w,m_i,aspectratio,ardry,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring,&
      anz,luse_tmatrix)
    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), fmelt(anz), aspectratio(anz), ardry(anz), &
                                    meltratio_outside, a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring, &
                                    hoststring, hostmatrixstring, hostinclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fa0, fb, fb0

    INTEGER :: fehler, i
    REAL(KIND=dp), DIMENSION(anz) :: volair, volice, volwater, volg, vg, &
                                     rho_g, Dm, x_g, x_w, meltrat_outs_max, &
                                     fm, fm_max, mro
    COMPLEX(kind=dp)              :: m_core(anz), m_airv(anz)

    m_airv = m_air

    fm = MAX(MIN(fmelt, 1d0), 0d0)
    mro = MAX(MIN(meltratio_outside, 1d0), 0d0)
    ! Outside melting fraction shall increase from given value (between 0 and 1)
    ! to 1, when melting degree fm approaches 1, such that melting particle
    ! asymptotically appproaches a water drop.
    ! Most simple approach linear:
    mro = mro + (1d0-mro)*fm

    Dm =  MAX(D_g,quasi_zero)
    x_g = a_geo * Dm**b_geo
    x_w = x_g * fm

    ! Volume of dry (ie still completely frozen) particle
    vg = pi6 * Dm**3 * ardry
    ! Bulk density of dry ice-air mixture (adjusted to 10kg/m3 <= rho_g <= rho_solid_ice)
    rho_g = MAX(MIN(x_g / vg, rho_ice_fwo), 10d0)
    ! rho-range adjusted volume of dry particle
    vg = x_g / rho_g

    meltrat_outs_max = 1d0 - rho_g / rho_w_fwo

    DO i=1,anz
      ! Core cavities won't get completely filled with meltwater
      IF (mro(i) <= meltrat_outs_max(i)) THEN
        volg(i) = vg(i) * (1d0 - mro(i) * fm(i))

      ! Core cavities might get completely filled with meltwater if melting
      ! degree is sufficiently high
      ELSE
        fm_max(i) = (rho_ice_fwo-rho_g(i)) / &
                    (mro(i)*rho_ice_fwo - rho_g(i) + rho_ice_fwo*rho_g(i)/rho_w_fwo)

        ! Not all cavities filled yet
        IF (fm(i) <= fm_max(i)) THEN
          volg(i) = (1d0 - mro(i) * fm(i)) * vg(i)
        ! All cavities water filled, ice skeleton starts to melt
        ELSE
          volg(i) = (x_g(i) - x_w(i)) / rho_ice_fwo + x_w(i) / rho_w_fwo
        END IF
      END IF
    END DO

    ! Vol-equiv particle diameter (hence not AR-scaling!):
    Dm  = (volg*inv_pi6) ** third

    volice = (x_g - x_w) / (volg * rho_ice_fwo)
    volwater = x_w / (rho_w_fwo * volg)
    volair = 1d0 - volice - volwater

    !..Complex refractive indices of air-ice-mixture of core
    m_core = get_m_mix_nested_vec(m_airv, m_i, m_w, &
                                  volair, volice, volwater, &
                                  mixingrulestring, hoststring, matrixstring, inclusionstring, &
                                  hostmatrixstring, hostinclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(Dm,Deps),&
                                     m_core,DBLE(m_air),aspectratio,lambda,&
                                     fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(Dm,Deps),&
                                 m_core,DBLE(m_air),lambda,&
                                 fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (Dm < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_SOAK_WETGR_VEC

  ! Backscattering cross section of melting graupel/hail/snow/ice accoring to Rayleigh
  ! approximation:   (for clarity, we keep lambda in the interface instead of Cfac = p^5/lambda^4)
  ! Vectorized version
  SUBROUTINE RAYLEIGH_SOAK_WETGR_VEC(x_g,a_geo,b_geo,fmelt,meltratio_outside,&
      m_w,m_i,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring,&
      anz)

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: x_g(anz), fmelt(anz), &
                                    meltratio_outside, a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w(anz), m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring
    CHARACTER(len=*), INTENT(in) :: hoststring, hostmatrixstring, hostinclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    COMPLEX(kind=dp) :: m_core(anz), m_airv(anz)
    REAL(KIND=dp), DIMENSION(anz) :: D_large, D_g, rho_g, x_w, xw_a, fm, fmgrenz,volg, vg, &
                                     volair, volice, volwater, meltratio_outside_grenz, mra
    INTEGER :: fehler, i

    m_airv = m_air

    fm = MAX(MIN(fmelt, 1d0), 0d0)
    mra = MAX(MIN(meltratio_outside, 1d0), 0d0)
    ! Der aussen schmelzende Anteil soll von dem eingegebenen Wert (zwischen 0 und 1)
    ! auf den Wert 1 ansteigen, wenn der Schmelzgrad fm gegen 1 geht,
    ! damit das Schmelzpartikel asymptotisch gegen einen Wassertropfen
    ! geht. Einfachster Ansatz linear:
    mra = mra + (1d0-mra)*fm

    x_w = x_g * fm

    D_g = MAX(a_geo * x_g**b_geo, quasi_zero)


    vg = pi6 * (D_g*D_g*D_g)
    rho_g = MAX(MIN(x_g / vg, rho_ice_fwo), 10d0)
    vg = x_g / rho_g

    meltratio_outside_grenz = 1d0 - rho_g / rho_w_fwo

    DO i=1,anz
      IF (mra(i) <= meltratio_outside_grenz(i)) THEN
        ! Bei diesem Partikel kann es nicht vorkommen, dass die Poren vollstaendig
        ! mit Schmelzwasser gefuellt werden.

        volg(i) = vg(i) * (1d0 - mra(i) * fm(i))

      ELSE

        ! Bei diesem Partikel kann es vorkommen, dass die Poren vollstaendig
        ! mit Schmelzwasser gefuellt werden, wenn der Schmelzgrad fm entsprechend hoch ist.

        fmgrenz(i) = (rho_ice_fwo-rho_g(i)) / (mra(i)*rho_ice_fwo - rho_g(i) + rho_ice_fwo*rho_g(i)/rho_w_fwo)

        IF (fm(i) <= fmgrenz(i)) THEN
          ! noch sind nicht alle Hohlraeume mit Wasser gefuellt.
          volg(i) = (1d0 - mra(i) * fm(i)) * vg(i)
        ELSE
          ! Alle Hohlraeume sind gefuellt, nun schmilzt auch das Eisskellet mit:
          volg(i) = (x_g(i) - x_w(i)) / rho_ice_fwo + x_w(i) / rho_w_fwo
        END IF

      END IF
    END DO

    D_large  = (inv_pi6 * volg) ** third

    volice = (x_g - x_w) / (volg * rho_ice_fwo)
    volwater = x_w / (rho_w_fwo * volg)
    volair = 1d0 - volice - volwater

    !..Komplexer Brechungindex der Luft-Eis-Mischung des Kerns:
    m_core = get_m_mix_nested_vec(m_airv, m_i, m_w, volair, volice, volwater, &
         mixingrulestring, hoststring, matrixstring, inclusionstring, &
         hostmatrixstring, hostinclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    WHERE (D_g >= 0.5d0*Deps)
      !..Rayleigh-Streuung
      ! UB: keep original formulation for clarity:
      !C_back = (ABS((m_core**2-1d0)/(m_core**2+2d0)))**2 * pi_dp**5 * D_large**6 / lambda**4

      !FIXME
      !JM190710: Not really convinced this is correct >:-/
      !  this means fa = fa0 for non-zero fa. odd because:
      !  (a) implies backscatt == extinction or
      !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
      !  both alternatives (a) and (b) seem fishy.
        fa  = (pi_dp**2 * D_large**3) / (2d0 * lambda**2) * sqrt(m_core)
        fb  = fa
        fa0 = fa
        fb0 = fb
    ELSEWHERE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE RAYLEIGH_SOAK_WETGR_VEC

  ! Originalfassung:   (for clarity, we keep lambda in the interface instead of Cfac = p^5/lambda^4)
  SUBROUTINE RAYLEIGH_SOAK_WETGR(x_g,a_geo,b_geo,fmelt,meltratio_outside,&
      m_w,m_i,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      hoststring,hostmatrixstring,hostinclusionstring)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)    :: x_g, a_geo, b_geo, fmelt, meltratio_outside, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_w, m_i
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring
    CHARACTER(len=*), INTENT(in) :: hoststring, hostmatrixstring, hostinclusionstring

    COMPLEX(kind=dp), INTENT(out) :: fa, fb, fa0, fb0

    COMPLEX(kind=dp) :: m_core
    REAL(KIND=dp) :: D_large, D_g, rho_g, x_w, xw_a, fm, fmgrenz,volg, vg, &
                     volair, volice, volwater, meltratio_outside_grenz, mra
    INTEGER :: fehler

    fm = MAX(MIN(fmelt, 1d0), 0d0)
    mra = MAX(MIN(meltratio_outside, 1d0), 0d0)
    ! Der aussen schmelzende Anteil soll von dem eingegebenen Wert (zwischen 0 und 1)
    ! auf den Wert 1 ansteigen, wenn der Schmelzgrad fm gegen 1 geht,
    ! damit das Schmelzpartikel asymptotisch gegen einen Wassertropfen
    ! geht. Einfachster Ansatz linear:
    mra = mra + (1d0-mra)*fm

    x_w = x_g * fm

    D_g = a_geo * x_g**b_geo

    IF (D_g >= 0.5d0*Deps) THEN

      vg = pi6 * (D_g*D_g*D_g)
      rho_g = MAX(MIN(x_g / vg, rho_ice_fwo), 10d0)
      vg = x_g / rho_g

      meltratio_outside_grenz = 1d0 - rho_g / rho_w_fwo

      IF (mra <= meltratio_outside_grenz) THEN
        ! Bei diesem Partikel kann es nicht vorkommen, dass die Poren vollstaendig
        ! mit Schmelzwasser gefuellt werden.

        volg = vg * (1d0 - mra * fm)

      ELSE

        ! Bei diesem Partikel kann es vorkommen, dass die Poren vollstaendig
        ! mit Schmelzwasser gefuellt werden, wenn der Schmelzgrad fm entsprechend hoch ist.

        fmgrenz = (rho_ice_fwo-rho_g) / (mra*rho_ice_fwo - rho_g + rho_ice_fwo*rho_g/rho_w_fwo)

        IF (fm <= fmgrenz) THEN
          ! noch sind nicht alle Hohlraeume mit Wasser gefuellt.
          volg = (1d0 - mra * fm) * vg
        ELSE
          ! Alle Hohlraeume sind gefuellt, nun schmilzt auch das Eisskellet mit:
          volg = (x_g - x_w) / rho_ice_fwo + x_w / rho_w_fwo
        END IF

      END IF

      D_large  = (inv_pi6 * volg) ** third

      volice = (x_g - x_w) / (volg * rho_ice_fwo)
      volwater = x_w / (rho_w_fwo * volg)
      volair = 1d0 - volice - volwater

      !..Komplexer Brechungindex der Luft-Eis-Mischung des Kerns:
      m_core = get_m_mix_nested(m_air, m_i, m_w, volair, volice, volwater, &
           mixingrulestring, hoststring, matrixstring, inclusionstring, &
           hostmatrixstring, hostinclusionstring, fehler)
      IF (fehler /= 0) THEN
        fa  = CMPLX(0d0,0d0,kind=dp)
        fb  = CMPLX(0d0,0d0,kind=dp)
        fa0 = CMPLX(0d0,0d0,kind=dp)
        fb0 = CMPLX(0d0,0d0,kind=dp)
        RETURN
      END IF

      !..Rayleigh-Streuung
      ! UB: keep original formulation for clarity:
      !C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2 * pi_dp**5 * D_large**6 / lambda**4

      !FIXME
      !JM190710: Not really convinced this is correct >:-/
      !  this means fa = fa0 for non-zero fa. odd because:
      !  (a) implies backscatt == extinction or
      !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
      !  both alternatives (a) and (b) seem fishy.
      fa  = (pi_dp**2 * D_large**3) / (2d0 * lambda**2) * sqrt(m_core)
      fb  = fa
      fa0 = fa
      fb0 = fb
    ELSE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)

    END IF

  END SUBROUTINE RAYLEIGH_SOAK_WETGR


  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für trockenen Schnee: Einschalenmodell
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_DRYSNOW_VEC(Ds,a_geo,b_geo,&
      m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: Ds(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz)    :: volice, volair, volwater, &
                                        rho_s, x_s, D_s
    COMPLEX(kind=dp), DIMENSION(anz) :: m_s, m_wdummy, m_airv

    m_airv = m_air
    m_wdummy = m_air * 5d0

    D_s    = MAX(Ds,quasi_zero)
    x_s   = a_geo * D_s**b_geo

    rho_s = MAX(MIN(x_s * inv_pi6 / (D_s*D_s*D_s), rho_ice_fwo), 5d0)
    D_s   = (x_s*inv_pi6/rho_s) ** third

    !..Complex refractive indices of air-ice-mixture:
    volice = rho_s / rho_ice_fwo
    volair = 1d0 - volice
    volwater = 0d0

    m_s = get_m_mix_vec(m_airv, m_i, m_wdummy, volair, volice, volwater, &
           mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_s,Deps),&
                                     m_s,DBLE(m_air),aspectratio,lambda,&
                                     fa,fb,fa0,fb0,anz)
    ELSE
    !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(D_s,Deps),&
                                 m_s,DBLE(m_air),lambda,&
                                 fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_s < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_DRYSNOW_VEC

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für trockenen Schnee: Zweischalenmodell
  ! Das Radienverhaeltnis der inneren zur aeusseren Schale wird hier als
  ! freier Parameter angesehen und muss vorgegeben werden. Erlaubt sind
  ! Werte 0 < radienverh < 1
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_DRYSNOW_TWOSPH_VEC(D_s,a_geo,b_geo,&
      m_i,Rratio,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_s(anz), aspectratio(anz), &
                                    Rratio, a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                    mixingrulestring_core, matrixstring_core, inclusionstring_core

    COMPLEX(kind=dp),INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz)    :: volice, volair, volwater, &
                                        volair_core, volice_core, volwater_core, &
                                        volair_shell, volice_shell, volwater_shell, &
                                        v_s_tot, v_s_core, &
                                        rho_s, rho_s_core, rho_s_shell, &
                                        D_s_tot, x_s_tot, D_s_core, x_s_core
    COMPLEX(kind=dp), DIMENSION(anz) :: m_s_shell, m_s_core, m_wdummy, m_airv

    m_airv = m_air
    m_wdummy = m_air * 5d0

    x_s_tot = a_geo * D_s**b_geo
    ! Volume of dry particle
    v_s_tot = pi6 * D_s**3 * aspectratio
    ! Bulk density of dry ice-air mixture (adjusted to 5kg/m3 <= rho_g <= rho_solid_ice)
    rho_s   = MAX(MIN(x_s_tot / v_s_tot, rho_ice_fwo), 5d0)
    ! rho-range adjusted volume and Dmax of dry ice-air particle
    v_s_tot = x_s_tot / rho_s
    D_s_tot = (v_s_tot*inv_pi6 / aspectratio) ** third

    D_s_core = D_s_tot * MAX(MIN(Rratio,1d0),0d0)
    x_s_core = a_geo * D_s_core**b_geo
    D_s_core = MAX(D_s_core,Deps*0.05d0)
    v_s_core = pi6 * D_s_core**3 * aspectratio
    rho_s_core = MAX(MIN(x_s_core / v_s_core, rho_ice_fwo), 5d0)
    ! The rho-max-min adjustment might have changed the radius ratio (Rratio) ---> adapt:
    v_s_core = x_s_core / rho_s_core
    D_s_core = (v_s_core*inv_pi6 / aspectratio) ** third

    rho_s_shell = MIN( (x_s_tot - x_s_core) / &
                         MAX((v_s_tot - v_s_core), Deps**3), &
                       rho_ice_fwo)

    !..Volume fractions of ice, water, air of core and shell:
    ! volice = v_ice/v_tot (with v_tot==v_airice)
    ! v_ice = m / rho_ice = a*D**b / rho_ice
    ! v_tot = v_airice = m / rho_airice => m = v_airice * rho_airice
    ! volice = (m / rho_ice) / (m / rho_airice) = rho_airice / rho_ice
    volice_core = rho_s_core / rho_ice_fwo
    ! dry particle
    volwater_core = 0d0
    volair_core = 1d0 - volice_core - volwater_core

    volice_shell = rho_s_shell / rho_ice_fwo
    volwater_shell = 0d0
    volair_shell = 1d0 - volice_shell - volwater_shell

    !..Complex refractive indices of air-ice-mixture of core and shell:
    m_s_core = get_m_mix_vec(m_airv, m_i, m_wdummy, &
                             volair_core, volice_core, volwater_core, &
                             mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                             anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    m_s_shell = get_m_mix_vec(m_airv, m_i, m_wdummy, &
                              volair_shell, volice_shell, volwater_shell, &
                              mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                              anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      !CALL SPHEROID_SCATTER_TMAT_VEC(MAX(D_s_tot,Deps),&
      !                               m_s_core,DBLE(m_air),aspectratio,lambda,&
      !                               fa,fb,fa0,fb0,anz)
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_s_tot,Deps),MAX(D_s_core,Deps),&
                                           m_s_core,m_s_shell,DBLE(m_air),aspectratio,lambda,&
                                           fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_s_tot,Deps),MAX(D_s_core,Deps),&
                                       m_s_core,m_s_shell,DBLE(m_air),lambda,&
                                       fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (D_s < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

    RETURN
  END SUBROUTINE MIE_DRYSNOW_TWOSPH_VEC


  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für trockenen Graupel/Hagel/Schnee: Einschalenmodell
  ! Effektiver Brechungsindex der Eis-Luft-Mischung wird durch die Schluesselwortparameter
  ! mixingrulestring, matrixstring und inclusionstring.
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_DRYGR_VEC(D_g,a_geo,b_geo,&
      m_i,aspectratio,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      anz,luse_tmatrix)
    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_g(anz), aspectratio(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out) :: fa(anz), fb(anz), fa0(anz), fb0(anz)

    INTEGER :: fehler
    REAL(KIND=dp), DIMENSION(anz)    :: volice, volair, volwater, rho, Dm, x_g, v_g
    COMPLEX(kind=dp), DIMENSION(anz) :: m_g, m_wdummy, m_airv

    m_airv = m_air
    m_wdummy = m_air * 5d0

    ! Adjust D, such that particle has a valid and sane density
    Dm = MAX(D_g,quasi_zero)
    x_g = a_geo * Dm**b_geo
    v_g = pi6 * Dm**3 * aspectratio
    rho = MAX(MIN( x_g / v_g, rho_ice_fwo), 10d0)
    ! volume of a particle with mass x_g and density rho
    v_g = x_g / rho
    ! Dveq of a particle with v_g
    Dm  = (v_g*inv_pi6) ** third

    !..Complex refractive index of air-ice-mixture according to Bohren and Battan
    volice = rho / rho_ice_fwo
    volair = 1d0 - volice
    volwater = 0d0

    m_g = get_m_mix_vec(m_airv, m_i, m_wdummy, &
                        volair, volice, volwater, &
                        mixingrulestring, matrixstring, inclusionstring, &
                        anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF(luse_tmatrix) THEN
      CALL SPHEROID_SCATTER_TMAT_VEC(MAX(Dm,Deps),&
                                     m_g,DBLE(m_air),aspectratio,lambda,&
                                     fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL SPHERE_SCATTER_BH_VEC(MAX(Dm,Deps),&
                                 m_g,DBLE(m_air),lambda,&
                                 fa,fb,fa0,fb0,anz)
    ENDIF

    WHERE (Dm < 1d-6)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE MIE_DRYGR_VEC

  ! Rayleigh-Naeherung des Rueckstreuquerschnitts von trockenem Graupel/Hagel/Schnee/Eis
  ! mit anderer Brechungsindexbeschreibung als Oguchi-Formel:
  ! Vektorisierte Ffassung:  (for clarity, we keep lambda in the interface instead of Cfac = p^5/lambda^4)
  SUBROUTINE RAYLEIGH_DRYGR_VEC(x_g,a_geo,b_geo,&
      m_i,lambda,&
     fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring,&
      anz)

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: x_g(anz), &
                                    a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    COMPLEX(kind=dp), DIMENSION(anz) :: m_g, m_airv
    REAL(KIND=dp), DIMENSION(anz)     :: volice, volair, volwater, rho_g, D_g
    INTEGER :: fehler

    m_airv = m_air

    D_g   = MAX(a_geo * x_g**b_geo, quasi_zero)

    rho_g = MAX(MIN(x_g * inv_pi6 / (D_g*D_g*D_g), rho_ice_fwo), 10d0)
    D_g   = (x_g*inv_pi6/rho_g) ** third

    !..Komplexer Brechungindex der Luft-Eis-Mischung nach Bohren und Battan
    volice = rho_g / rho_ice_fwo
    volair = 1d0 - volice
    volwater = 0d0

    m_g = get_m_mix_vec(m_airv, m_i, 5d0*m_airv, volair, volice, volwater, &
         mixingrulestring, matrixstring, inclusionstring, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    !..Mie-Streuung
    WHERE (D_g >= 1d-12)
      ! UB: keep original formulation for clarity:
      !C_back = pi_dp**5 * (ABS((m_g**2-1.0d0)/(m_g**2+2.0d0)))**2 * D_g**6 / lambda**4

      !FIXME
      !JM190710: Not really convinced this is correct >:-/
      !  this means fa = fa0 for non-zero fa. odd because:
      !  (a) implies backscatt == extinction or
      !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
      !  both alternatives (a) and (b) seem fishy.
      fa  = (pi_dp**2 * D_g**3) / (2d0 * lambda**2) * sqrt(m_g)
      fb  = fa
      fa0 = fa
      fb0 = fb
    ELSEWHERE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

  END SUBROUTINE RAYLEIGH_DRYGR_VEC

  ! Originalfassung:  (for clarity, we keep lambda in the interface instead of Cfac = p^5/lambda^4)
  SUBROUTINE RAYLEIGH_DRYGR(x_g,a_geo,b_geo,&
      m_i,lambda,&
      fa,fb,fa0,fb0,&
      mixingrulestring,matrixstring,inclusionstring)

    IMPLICIT NONE

    REAL(KIND=dp), INTENT(in)    :: x_g, a_geo, b_geo, lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i
    CHARACTER(len=*), INTENT(in) :: mixingrulestring, matrixstring, inclusionstring

    COMPLEX(kind=dp), INTENT(out) :: fa, fb, fa0, fb0

    COMPLEX(kind=dp) :: m_g
    REAL(KIND=dp)     :: volice, volair, volwater, rho_g, D_g
    INTEGER :: fehler

    D_g   = a_geo * x_g**b_geo

    IF (D_g > 0.5d0*Deps) THEN

      rho_g = MAX(MIN(x_g * inv_pi6 / (D_g*D_g*D_g), rho_ice_fwo), 10d0)
      D_g   = (x_g*inv_pi6/rho_g) ** third

      !..Komplexer Brechungindex der Luft-Eis-Mischung nach Bohren und Battan
      volice = rho_g / rho_ice_fwo
      volair = 1d0 - volice
      volwater = 0d0

      m_g = get_m_mix(m_air, m_i, 5d0*m_air, volair, volice, volwater, &
           mixingrulestring, matrixstring, inclusionstring, fehler)
      IF (fehler /= 0) THEN
        fa  = CMPLX(0d0,0d0,kind=dp)
        fb  = CMPLX(0d0,0d0,kind=dp)
        fa0 = CMPLX(0d0,0d0,kind=dp)
        fb0 = CMPLX(0d0,0d0,kind=dp)
        RETURN
      END IF

      !..Mie-Streuung
      ! UB: keep original formulation for clarity:
      !C_back = pi_dp**5 * (ABS((m_g**2-1.0d0)/(m_g**2+2.0d0)))**2 * D_g**6 / lambda**4

      !FIXME
      !JM190710: Not really convinced this is correct >:-/
      !  this means fa = fa0 for non-zero fa. odd because:
      !  (a) implies backscatt == extinction or
      !  (b) considering that extinction will be from IMAG(fa0), implies zero extinction.
      !  both alternatives (a) and (b) seem fishy.
      fa  = (pi_dp**2 * D_g**3) / (2d0 * lambda**2) * sqrt(m_g)
      fb  = fa
      fa0 = fa
      fb0 = fb
    ELSE
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END IF

  END SUBROUTINE RAYLEIGH_DRYGR

  !----------------------------------------------------------------------------------------!
  ! Rückstreuquerschnitt für schmelzenden Schnee: Zweischalenmodell
  ! Das Radienverhaeltnis der inneren zur aeusseren Schale wird hier als
  ! freier Parameter angesehen und muss vorgegeben werden. Erlaubt sind
  ! Werte 0 < radienverh < 1
  !----------------------------------------------------------------------------------------!

  ! Vectorized version with input in D-space
  SUBROUTINE MIE_WETSNOW_TWOSPH_VEC(&
      D_s,a_geo,b_geo,fmelt,meltingratio_outside,&
      m_w,m_i,aspectratio,ardry,lambda,Rratio,&
      fa,fb,fa0,fb0,&
      mixingrulestring_shell,matrixstring_shell,inclusionstring_shell,&
      mixingrulestring_core,matrixstring_core,inclusionstring_core,&
      hoststring_shell,hostmatrixstring_shell,hostinclusionstring_shell,&
      hoststring_core,hostmatrixstring_core,hostinclusionstring_core,&
      anz,luse_tmatrix)

    IMPLICIT NONE

    LOGICAL, INTENT(in)          :: luse_tmatrix
    INTEGER, INTENT(in)          :: anz
    REAL(KIND=dp), INTENT(in)    :: D_s(anz), fmelt(anz), aspectratio(anz), ardry(anz), &
                                    Rratio, meltingratio_outside, a_geo, b_geo,lambda
    COMPLEX(kind=dp), INTENT(in) :: m_i(anz),m_w(anz)
    CHARACTER(len=*), INTENT(in) :: mixingrulestring_shell, matrixstring_shell, inclusionstring_shell, &
                                    mixingrulestring_core, matrixstring_core, inclusionstring_core, &
                                    hoststring_shell, hostmatrixstring_shell, hostinclusionstring_shell, &
                                    hoststring_core, hostmatrixstring_core, hostinclusionstring_core

    COMPLEX(kind=dp), INTENT(out), DIMENSION(anz) :: fa, fb, fa0, fb0

    INTEGER :: fehler, i
    REAL(KIND=dp), DIMENSION(anz)    :: volice, volair, volwater, &
                                        volair_core, volice_core, volwater_core, &
                                        volair_shell, volice_shell, volwater_shell, &
                                        vol_tot, vol_core, vol_shell, &
                                        vol_wet_core, vol_wet_shell, &
                                        v_s_tot, v_s_core, v_s_shell, &
                                        rho_s, rho_s_core, rho_s_shell, &
                                        D_s_tot, x_s_tot, D_s_core, x_s_core, &
                                        D_s_core_melting, D_s_melting, &
                                        meltrat_outs_max_core, meltrat_outs_max_shell, &
                                        fm, fm_max_core, fm_max_shell, mro
    COMPLEX(kind=dp), DIMENSION(anz) :: m_s_shell, m_s_core, m_airv

    REAL(KIND=dp), PARAMETER :: Deps3 = Deps*Deps*Deps

    m_airv = m_air

    fm = MAX(MIN(fmelt,1d0),0d0)
    mro = MAX(MIN(meltingratio_outside,1d0),0d0)
    ! Outside melting fraction shall increase from given value (between 0 and 1)
    ! to 1, when melting degree fm approaches 1, such that melting particle
    ! asymptotically appproaches a water drop.
    ! Most simple approach linear:
    mro = mro + (1d0-mro)*fm

    x_s_tot = a_geo * D_s**b_geo
    ! Volume of dry particle
    vol_tot = pi6 * D_s**3 * ardry
    ! Bulk density of dry ice-air mixture (adjusted to 5kg/m3 <= rho_g <= rho_solid_ice)
    rho_s   = MAX(MIN(x_s_tot / vol_tot, rho_ice_fwo), 5d0)
    ! rho-range adjusted volume and Dmax of dry ice-air particle
    vol_tot = x_s_tot / rho_s
    D_s_tot = (vol_tot*inv_pi6 / ardry) ** third

    D_s_core = D_s_tot * MAX(MIN(Rratio,1d0),0d0)
    x_s_core = a_geo * D_s_core**b_geo
    D_s_core = MAX(D_s_core,Deps*0.05d0)
    vol_core = pi6 * D_s_core**3 * ardry
    rho_s_core = MAX(MIN(x_s_core / vol_core, rho_ice_fwo), 5d0)
    ! The rho-max-min adjustment might have changed the radius ratio (Rratio) ---> adapt:
    vol_core = x_s_core / rho_s_core
    !D_s_core = (vol_core*inv_pi6 / aspectratio) ** third  ! commented out because not needed below, but keep for clarity

    vol_shell = vol_tot - vol_core

    rho_s_shell = MIN( (x_s_tot - x_s_core) / &
                         MAX(vol_shell, Deps3), &
                       rho_ice_fwo)

    meltrat_outs_max_core = 1d0 - rho_s_core / rho_w_fwo
    meltrat_outs_max_shell = 1d0 - rho_s_shell / rho_w_fwo

    DO i=1,anz

      ! First calculate volume of core:

      ! Core cavities won't get completely filled with meltwater
      IF (mro(i) < meltrat_outs_max_core(i)) THEN
        vol_wet_core(i) = vol_core(i) * (1d0 - mro(i) * fm(i))

      ! Core cavities might get completely filled with meltwater if melting
      ! degree is sufficiently high
      ELSE
        fm_max_core(i) = (rho_ice_fwo-rho_s_core(i)) / &
                         (mro(i)*rho_ice_fwo - rho_s_core(i) + &
                          rho_ice_fwo*rho_s_core(i)/rho_w_fwo)

        ! Not all cavities filled yet
        IF (fm(i) <= fm_max_core(i)) THEN
          vol_wet_core(i) = (1d0 - mro(i) * fm(i)) * vol_core(i)
        ! All cavities water filled, ice skeleton starts to melt
        ELSE
          vol_wet_core(i) = x_s_core(i) * (1d0-fm(i)) / rho_ice_fwo + &
                            x_s_core(i) * fm(i) / rho_w_fwo
        END IF
      END IF

      ! Now calculate volume of shell:

      ! Shell cavities won't get completely filled with meltwater
      IF (mro(i) < meltrat_outs_max_shell(i)) THEN
        vol_wet_shell(i) = vol_shell(i) * (1d0 - mro(i) * fm(i))

      ! Shell cavities might get completely filled with meltwater if melting
      ! degree is sufficiently high
      ELSE
        fm_max_shell(i) = (rho_ice_fwo-rho_s_shell(i)) / (mro(i)*rho_ice_fwo - rho_s_shell(i) + &
                          rho_ice_fwo*rho_s_shell(i)/rho_w_fwo)

        ! Not all cavities filled yet
        IF (fm(i) <= fm_max_shell(i)) THEN
          vol_wet_shell(i) = (1d0 - mro(i) * fm(i)) * vol_shell(i)
        ! All cavities water filled, ice skeleton starts to melt
        ELSE
          vol_wet_shell(i) = (x_s_tot(i)-x_s_core(i)) * (1d0-fm(i)) / rho_ice_fwo + &
                             (x_s_tot(i)-x_s_core(i)) * fm(i) / rho_w_fwo
        END IF
      END IF

    END DO

    ! Vol-equiv particle diameters (hence not AR-scaling!) of core and shell:
    D_s_core_melting = (inv_pi6*vol_wet_core) ** third
    D_s_melting = (inv_pi6*(vol_wet_core+vol_wet_shell)) ** third

    !..Complex refractive indices of air-ice-mixture of core and shell
    ! 1) Core:
    volice_core = x_s_core * (1d0-fm)/ (rho_ice_fwo * MAX(vol_wet_core, Deps3))
    volwater_core = x_s_core * fm / (rho_w_fwo * MAX(vol_wet_core, Deps3))
    volair_core = 1d0 - volice_core - volwater_core

    m_s_core = get_m_mix_nested_vec(m_airv, m_i, m_w, volair_core, volice_core, volwater_core, &
         mixingrulestring_core, hoststring_core, matrixstring_core, inclusionstring_core, &
         hostmatrixstring_core, hostinclusionstring_core, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    ! 2) Shell:
    volice_shell = (x_s_tot-x_s_core) * (1d0 - fm) / (rho_ice_fwo * MAX(vol_wet_shell,Deps3))
    volwater_shell = (x_s_tot-x_s_core) * fm / (rho_w_fwo * MAX(vol_wet_shell,Deps3))
    volair_shell = 1d0 - volice_shell - volwater_shell

    m_s_shell = get_m_mix_nested_vec(m_airv, m_i, m_w, volair_shell, volice_shell, volwater_shell, &
         mixingrulestring_shell, hoststring_shell, matrixstring_shell, inclusionstring_shell, &
         hostmatrixstring_shell, hostinclusionstring_shell, anz, fehler)
    IF (fehler /= 0) THEN
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
      RETURN
    END IF

    IF (luse_tmatrix) THEN
      CALL COATEDSPHEROID_SCATTER_TMAT_VEC(MAX(D_s_melting,Deps),&
                                           MAX(D_s_core_melting,Deps),&
                                           m_s_core,m_s_shell,DBLE(m_air),aspectratio,&
                                           lambda,fa,fb,fa0,fb0,anz)
    ELSE
      !..Mie-Streuung
      CALL COATEDSPHERE_SCATTER_BH_VEC(MAX(D_s_melting,Deps),&
                                       MAX(D_s_core_melting,Deps),&
                                       m_s_core,m_s_shell,DBLE(m_air),lambda,&
                                       fa,fb,fa0,fb0,anz)
    END IF

    WHERE (D_s < Deps)
      fa  = CMPLX(0d0,0d0,kind=dp)
      fb  = CMPLX(0d0,0d0,kind=dp)
      fa0 = CMPLX(0d0,0d0,kind=dp)
      fb0 = CMPLX(0d0,0d0,kind=dp)
    END WHERE

    RETURN
  END SUBROUTINE MIE_WETSNOW_TWOSPH_VEC


END MODULE radar_mielib_vec




