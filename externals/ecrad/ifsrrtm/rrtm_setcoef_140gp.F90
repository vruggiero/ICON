! This file has been modified for the use in ICON

SUBROUTINE RRTM_SETCOEF_140GP (KIDIA,KFDIA,KLEV,P_COLDRY,P_WBROAD,P_WKL,&
 & P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,&
 & P_COLH2O,P_COLCO2,P_COLO3,P_COLN2O,P_COLCH4, P_COLO2,P_CO2MULT, P_COLBRD, &
 & K_LAYTROP,K_LAYSWTCH,K_LAYLOW,PAVEL,P_TAVEL,P_SELFFAC,P_SELFFRAC,K_INDSELF,&
 & K_INDMINOR,P_SCALEMINOR,P_SCALEMINORN2,P_MINORFRAC,&
 & PRAT_H2OCO2, PRAT_H2OCO2_1, PRAT_H2OO3, PRAT_H2OO3_1, &
 & PRAT_H2ON2O, PRAT_H2ON2O_1, PRAT_H2OCH4, PRAT_H2OCH4_1, &
 & PRAT_N2OCO2, PRAT_N2OCO2_1, PRAT_O3CO2, PRAT_O3CO2_1)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714
!        NEC           25-Oct-2007 Optimisations
!     201305 ABozzo updated to rrtmg_lw_v4.85
!     201507 RHogan Bug fix: swapped P_COLO2 & P_CO2MULT in argument list


!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.
!     Also calculate the values of the integrated Planck functions 
!     for each band at the level and layer temperatures.

USE PARKIND1 , ONLY : JPIM, JPRB
USE ecradhook  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM  , ONLY : JPINPX
USE YOERRTRF , ONLY : PREFLOG   ,TREF, CHI_MLS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLDRY(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_WBROAD(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLBRD(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_WKL(KIDIA:KFDIA,JPINPX,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FAC11(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_FORFRAC(KIDIA:KFDIA,KLEV)
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_JT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLCO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLO3(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLN2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLCH4(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_COLO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_CO2MULT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYTROP(KIDIA:KFDIA) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYSWTCH(KIDIA:KFDIA) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_LAYLOW(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAVEL(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_INDSELF(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_INDFOR(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: K_INDMINOR(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SCALEMINOR(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_SCALEMINORN2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_MINORFRAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: &                 !
                     & PRAT_H2OCO2(KIDIA:KFDIA,KLEV),PRAT_H2OCO2_1(KIDIA:KFDIA,KLEV), &
                     & PRAT_H2OO3(KIDIA:KFDIA,KLEV) ,PRAT_H2OO3_1(KIDIA:KFDIA,KLEV), & 
                     & PRAT_H2ON2O(KIDIA:KFDIA,KLEV),PRAT_H2ON2O_1(KIDIA:KFDIA,KLEV), &
                     & PRAT_H2OCH4(KIDIA:KFDIA,KLEV),PRAT_H2OCH4_1(KIDIA:KFDIA,KLEV), &
                     & PRAT_N2OCO2(KIDIA:KFDIA,KLEV),PRAT_N2OCO2_1(KIDIA:KFDIA,KLEV), &
                     & PRAT_O3CO2(KIDIA:KFDIA,KLEV) ,PRAT_O3CO2_1(KIDIA:KFDIA,KLEV)
!- from INTFAC      
!- from INTIND
!- from PROFDATA             
!- from PROFILE             
!- from SELF             
INTEGER(KIND=JPIM) :: JP1, JLAY
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) :: Z_CO2REG, Z_COMPFP, Z_FACTOR, Z_FP, Z_FT, Z_FT1, Z_PLOG, Z_SCALEFAC, Z_STPFAC, Z_WATER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('RRTM_SETCOEF_140GP',0,ZHOOK_HANDLE)

!$ACC PARALLEL DEFAULT(NONE) PRESENT(P_COLDRY, P_WBROAD, P_COLBRD, P_WKL, P_FAC00, P_FAC01, P_FAC10, P_FAC11, &
!$ACC   P_FORFAC, P_FORFRAC, K_JP, K_JT, K_JT1, P_COLH2O, P_COLCO2, P_COLO3, P_COLN2O, P_COLCH4, P_COLO2, P_CO2MULT, &
!$ACC   K_LAYTROP, K_LAYSWTCH, K_LAYLOW, PAVEL, P_TAVEL, P_SELFFAC, P_SELFFRAC, K_INDSELF, K_INDFOR, K_INDMINOR, &
!$ACC   P_SCALEMINOR, P_SCALEMINORN2, P_MINORFRAC, PRAT_H2OCO2, PRAT_H2OCO2_1, PRAT_H2OO3, PRAT_H2OO3_1, PRAT_H2ON2O, &
!$ACC   PRAT_H2ON2O_1, PRAT_H2OCH4, PRAT_H2OCH4_1, PRAT_N2OCO2, PRAT_N2OCO2_1, PRAT_O3CO2, PRAT_O3CO2_1) ASYNC(1)
!$ACC LOOP GANG VECTOR
DO JLON = KIDIA, KFDIA
  Z_STPFAC = 296._JPRB/1013._JPRB

  K_LAYTROP(JLON)  = 0
  K_LAYSWTCH(JLON) = 0
  K_LAYLOW(JLON)   = 0
  !$ACC LOOP SEQ
  DO JLAY = 1, KLEV
!        Find the two reference pressures on either side of the
!        layer pressure.  Store them in JP and JP1.  Store in FP the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.
    Z_PLOG = LOG(PAVEL(JLON,JLAY))
    K_JP(JLON,JLAY) = INT(36._JPRB - 5*(Z_PLOG+0.04_JPRB))
    IF (K_JP(JLON,JLAY)  <  1) THEN
      K_JP(JLON,JLAY) = 1
    ELSEIF (K_JP(JLON,JLAY)  >  58) THEN
      K_JP(JLON,JLAY) = 58
    ENDIF
    JP1 = K_JP(JLON,JLAY) + 1
    Z_FP = 5._JPRB * (PREFLOG(K_JP(JLON,JLAY)) - Z_PLOG)
! bound Z_FP in case Z_PLOG is outside range of ref. pressure PREFLOG
! (in LVERTFE, pressure at last full level is known, but not in finite diff (NH)
    Z_FP = MAX(-1.0_JPRB, MIN(1.0_JPRB, Z_FP))
!        Determine, for each reference pressure (JP and JP1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  Store these indices in JT and JT1, resp.
!        Store in FT (resp. FT1) the fraction of the way between JT
!        (JT1) and the next highest reference temperature that the 
!        layer temperature falls.

    K_JT(JLON,JLAY) = INT(3._JPRB + (P_TAVEL(JLON,JLAY)-TREF(K_JP(JLON,JLAY)))/15._JPRB)
    IF (K_JT(JLON,JLAY)  <  1) THEN
      K_JT(JLON,JLAY) = 1
    ELSEIF (K_JT(JLON,JLAY)  >  4) THEN
      K_JT(JLON,JLAY) = 4
    ENDIF
    Z_FT = ((P_TAVEL(JLON,JLAY)-TREF(K_JP(JLON,JLAY)))/15._JPRB) - REAL(K_JT(JLON,JLAY)-3)
    K_JT1(JLON,JLAY) = INT(3._JPRB + (P_TAVEL(JLON,JLAY)-TREF(JP1))/15._JPRB)
    IF (K_JT1(JLON,JLAY)  <  1) THEN
      K_JT1(JLON,JLAY) = 1
    ELSEIF (K_JT1(JLON,JLAY)  >  4) THEN
      K_JT1(JLON,JLAY) = 4
    ENDIF
    Z_FT1 = ((P_TAVEL(JLON,JLAY)-TREF(JP1))/15._JPRB) - REAL(K_JT1(JLON,JLAY)-3)

    Z_WATER = P_WKL(JLON,1,JLAY)/P_COLDRY(JLON,JLAY)
    Z_SCALEFAC = PAVEL(JLON,JLAY) * Z_STPFAC / P_TAVEL(JLON,JLAY)

!        If the pressure is less than ~100mb, perform a different
!        set of species interpolations.
!         IF (PLOG .LE. 4.56) GO TO 5300
!--------------------------------------         
    IF (Z_PLOG  >  4.56_JPRB) THEN
      K_LAYTROP(JLON) =  K_LAYTROP(JLON) + 1
!        For one band, the "switch" occurs at ~300 mb. 
!      IF (Z_PLOG  >=  5.76_JPRB) K_LAYSWTCH(JLON) = K_LAYSWTCH(JLON) + 1
!      IF (Z_PLOG  >=  6.62_JPRB) K_LAYLOW(JLON) = K_LAYLOW(JLON) + 1


!        water vapor foreign continuum
      P_FORFAC(JLON,JLAY) = Z_SCALEFAC / (1.0_JPRB+Z_WATER)
      Z_FACTOR = (332.0_JPRB-P_TAVEL(JLON,JLAY))/36.0_JPRB
      K_INDFOR(JLON,JLAY) = MIN(2, MAX(1, INT(Z_FACTOR)))
      P_FORFRAC(JLON,JLAY) = Z_FACTOR - REAL(K_INDFOR(JLON,JLAY))

!        Set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.
!C           SELFFAC(LAY) = WATER * SCALEFAC / (1.+WATER)
      P_SELFFAC(JLON,JLAY) = Z_WATER * P_FORFAC(JLON,JLAY)
      Z_FACTOR = (P_TAVEL(JLON,JLAY)-188.0_JPRB)/7.2_JPRB
      K_INDSELF(JLON,JLAY) = MIN(9, MAX(1, INT(Z_FACTOR)-7))
      P_SELFFRAC(JLON,JLAY) = Z_FACTOR - REAL(K_INDSELF(JLON,JLAY) + 7)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         P_SCALEMINOR(JLON,JLAY) = PAVEL(JLON,JLAY)/P_TAVEL(JLON,JLAY)
         P_SCALEMINORN2(JLON,JLAY) = (PAVEL(JLON,JLAY)/P_TAVEL(JLON,JLAY)) &
           &  *(P_WBROAD(JLON,JLAY)/(P_COLDRY(JLON,JLAY)+P_WKL(JLON,1,JLAY)))
         Z_FACTOR = (P_TAVEL(JLON,JLAY)-180.8_JPRB)/7.2_JPRB
         K_INDMINOR(JLON,JLAY) = MIN(18, MAX(1, INT(Z_FACTOR)))
         P_MINORFRAC(JLON,JLAY) = Z_FACTOR - REAL(K_INDMINOR(JLON,JLAY))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in lower atmosphere.
         PRAT_H2OCO2(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY))/CHI_MLS(2,K_JP(JLON,JLAY))
         PRAT_H2OCO2_1(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY)+1)/CHI_MLS(2,K_JP(JLON,JLAY)+1)

         PRAT_H2OO3(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY))/CHI_MLS(3,K_JP(JLON,JLAY))
         PRAT_H2OO3_1(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY)+1)/CHI_MLS(3,K_JP(JLON,JLAY)+1)

         PRAT_H2ON2O(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY))/CHI_MLS(4,K_JP(JLON,JLAY))
         PRAT_H2ON2O_1(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY)+1)/CHI_MLS(4,K_JP(JLON,JLAY)+1)

         PRAT_H2OCH4(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY))/CHI_MLS(6,K_JP(JLON,JLAY))
         PRAT_H2OCH4_1(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY)+1)/CHI_MLS(6,K_JP(JLON,JLAY)+1)

         PRAT_N2OCO2(JLON,JLAY)=CHI_MLS(4,K_JP(JLON,JLAY))/CHI_MLS(2,K_JP(JLON,JLAY))
         PRAT_N2OCO2_1(JLON,JLAY)=CHI_MLS(4,K_JP(JLON,JLAY)+1)/CHI_MLS(2,K_JP(JLON,JLAY)+1)



!        Calculate needed column amounts.
      P_COLH2O(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,1,JLAY)
      P_COLCO2(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,2,JLAY)
      P_COLO3(JLON,JLAY)  = 1.E-20_JPRB * P_WKL(JLON,3,JLAY)
      P_COLN2O(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,4,JLAY)
      P_COLCH4(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,6,JLAY)
      P_COLO2(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,7,JLAY)
      P_COLBRD(JLON,JLAY) = 1.E-20_JPRB * P_WBROAD(JLON,JLAY)
      IF (P_COLCO2(JLON,JLAY)  ==  0.0_JPRB) P_COLCO2(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
      IF (P_COLN2O(JLON,JLAY)  ==  0.0_JPRB) P_COLN2O(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
      IF (P_COLCH4(JLON,JLAY)  ==  0.0_JPRB) P_COLCH4(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
!        Using E = 1334.2 cm-1.
      Z_CO2REG = 3.55E-24_JPRB * P_COLDRY(JLON,JLAY)
      P_CO2MULT(JLON,JLAY)= (P_COLCO2(JLON,JLAY) - Z_CO2REG) *&
       & 272.63_JPRB*EXP(-1919.4_JPRB/P_TAVEL(JLON,JLAY))/(8.7604E-4_JPRB*P_TAVEL(JLON,JLAY))  
!         GO TO 5400
!------------------
    ELSE
!        Above LAYTROP.
! 5300    CONTINUE

!        Calculate needed column amounts.
      P_FORFAC(JLON,JLAY) = Z_SCALEFAC / (1.0_JPRB+Z_WATER)
      Z_FACTOR = (P_TAVEL(JLON,JLAY)-188.0_JPRB)/36.0_JPRB
      K_INDFOR(JLON,JLAY) = 3
      P_FORFRAC(JLON,JLAY) = Z_FACTOR - 1.0_JPRB

!  Set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
      P_SELFFAC(JLON,JLAY) = Z_WATER * P_FORFAC(JLON,JLAY)

!  Set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
      P_SCALEMINOR(JLON,JLAY) = PAVEL(JLON,JLAY)/P_TAVEL(JLON,JLAY)         
      P_SCALEMINORN2(JLON,JLAY) = (PAVEL(JLON,JLAY)/P_TAVEL(JLON,JLAY)) &
        &    * (P_WBROAD(JLON,JLAY)/(P_COLDRY(JLON,JLAY)+P_WKL(JLON,1,JLAY)))
      Z_FACTOR = (P_TAVEL(JLON,JLAY)-180.8_JPRB)/7.2_JPRB
      K_INDMINOR(JLON,JLAY) = MIN(18, MAX(1, INT(Z_FACTOR)))
      P_MINORFRAC(JLON,JLAY) = Z_FACTOR - REAL(K_INDMINOR(JLON,JLAY))

!  Setup reference ratio to be used in calculation of binary
!  species parameter in upper atmosphere.
      PRAT_H2OCO2(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY))/CHI_MLS(2,K_JP(JLON,JLAY))
      PRAT_H2OCO2_1(JLON,JLAY)=CHI_MLS(1,K_JP(JLON,JLAY)+1)/CHI_MLS(2,K_JP(JLON,JLAY)+1)         

      PRAT_O3CO2(JLON,JLAY)=CHI_MLS(3,K_JP(JLON,JLAY))/CHI_MLS(2,K_JP(JLON,JLAY))
      PRAT_O3CO2_1(JLON,JLAY)=CHI_MLS(3,K_JP(JLON,JLAY)+1)/CHI_MLS(2,K_JP(JLON,JLAY)+1)         


!  Calculate needed column amounts.
      P_COLH2O(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,1,JLAY)
      P_COLCO2(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,2,JLAY)
      P_COLO3(JLON,JLAY)  = 1.E-20_JPRB * P_WKL(JLON,3,JLAY)
      P_COLN2O(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,4,JLAY)
      P_COLCH4(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,6,JLAY)
      P_COLO2(JLON,JLAY) = 1.E-20_JPRB * P_WKL(JLON,7,JLAY)
      P_COLBRD(JLON,JLAY) = 1.E-20_JPRB * P_WBROAD(JLON,JLAY)
      IF (P_COLCO2(JLON,JLAY)  ==  0.0_JPRB) P_COLCO2(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
      IF (P_COLN2O(JLON,JLAY)  ==  0.0_JPRB) P_COLN2O(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
      IF (P_COLCH4(JLON,JLAY)  ==  0.0_JPRB) P_COLCH4(JLON,JLAY) = 1.E-32_JPRB * P_COLDRY(JLON,JLAY)
      Z_CO2REG = 3.55E-24_JPRB * P_COLDRY(JLON,JLAY)
      P_CO2MULT(JLON,JLAY)= (P_COLCO2(JLON,JLAY) - Z_CO2REG) *&
       & 272.63_JPRB*EXP(-1919.4_JPRB/P_TAVEL(JLON,JLAY))/(8.7604E-4_JPRB*P_TAVEL(JLON,JLAY))  
!----------------     
    ENDIF
! 5400    CONTINUE

!        We have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  We multiply the pressure 
!        fraction FP with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines TAUGBn for band n).

    Z_COMPFP = 1.0_JPRB - Z_FP
    P_FAC10(JLON,JLAY) = Z_COMPFP * Z_FT
    P_FAC00(JLON,JLAY) = Z_COMPFP * (1.0_JPRB - Z_FT)
    P_FAC11(JLON,JLAY) = Z_FP * Z_FT1
    P_FAC01(JLON,JLAY) = Z_FP * (1.0_JPRB - Z_FT1)

!  Rescale selffac and forfac for use in taumol
    P_SELFFAC(JLON,JLAY) = P_COLH2O(JLON,JLAY)*P_SELFFAC(JLON,JLAY)
    P_FORFAC(JLON,JLAY) = P_COLH2O(JLON,JLAY)*P_FORFAC(JLON,JLAY)


  ENDDO

! MT 981104 
!-- Set LAYLOW for profiles with surface pressure less than 750 hPa. 
  IF (K_LAYLOW(JLON) == 0) K_LAYLOW(JLON)=1
ENDDO
!$ACC END PARALLEL

IF (LHOOK) CALL DR_HOOK('RRTM_SETCOEF_140GP',1,ZHOOK_HANDLE)

END SUBROUTINE RRTM_SETCOEF_140GP
