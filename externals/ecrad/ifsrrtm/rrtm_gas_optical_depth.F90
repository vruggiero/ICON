! This file has been modified for the use in ICON

!option! -pvctl no_on_adb
!option! -pvctl nocollapse
SUBROUTINE RRTM_GAS_OPTICAL_DEPTH(KIDIA,KFDIA,KLEV,POD,PAVEL, PCOLDRY,PCOLBRD,PWX,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,PCOLO3,PCOLN2O,PCOLCH4,PCOLO2,P_CO2MULT,&
 & KLAYTROP,KLAYSWTCH,KLAYLOW,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & KINDMINOR,PSCALEMINOR,PSCALEMINORN2,PMINORFRAC,&
 &  PRAT_H2OCO2, PRAT_H2OCO2_1, PRAT_H2OO3, PRAT_H2OO3_1, &
 &  PRAT_H2ON2O, PRAT_H2ON2O_1, PRAT_H2OCH4, PRAT_H2OCH4_1, &
 &  PRAT_N2OCO2, PRAT_N2OCO2_1, PRAT_O3CO2, PRAT_O3CO2_1)  

! R. J. Hogan   20160831  Adapted from rrtm_gasabs1a_140gp.F90

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE ecradhook   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND   ,JPXSEC
USE YOERRTM  , ONLY : JPGPT

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POD(JPGPT,KLEV,KIDIA:KFDIA) ! Optical depth: note different orientation
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLDRY(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWX(KIDIA:KFDIA,JPXSEC,KLEV) ! Amount of trace gases
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUAERL(KIDIA:KFDIA,KLEV,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC00(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC01(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC10(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAC11(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KJP(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KJT1(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PONEMINUS

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLH2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO3(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLN2O(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLCH4(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLO2(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_CO2MULT(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAYTROP(KIDIA:KFDIA) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAYSWTCH(KIDIA:KFDIA) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAYLOW(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSELFFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDSELF(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRAC(KIDIA:KFDIA,JPGPT,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORFRAC(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDFOR(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMINORFRAC(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCALEMINOR(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCALEMINORN2(KIDIA:KFDIA,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDMINOR(KIDIA:KFDIA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLBRD(KIDIA:KFDIA,KLEV)  
REAL(KIND=JPRB)  , INTENT(IN) :: &                  !
                    &   PRAT_H2OCO2(KIDIA:KFDIA,KLEV),PRAT_H2OCO2_1(KIDIA:KFDIA,KLEV), &
                    &   PRAT_H2OO3(KIDIA:KFDIA,KLEV),PRAT_H2OO3_1(KIDIA:KFDIA,KLEV), & !    DIMENSIONS: (NLAYERS)
                    &   PRAT_H2ON2O(KIDIA:KFDIA,KLEV),PRAT_H2ON2O_1(KIDIA:KFDIA,KLEV), &
                    &   PRAT_H2OCH4(KIDIA:KFDIA,KLEV),PRAT_H2OCH4_1(KIDIA:KFDIA,KLEV), &
                    &   PRAT_N2OCO2(KIDIA:KFDIA,KLEV),PRAT_N2OCO2_1(KIDIA:KFDIA,KLEV), &
                    &   PRAT_O3CO2(KIDIA:KFDIA,KLEV),PRAT_O3CO2_1(KIDIA:KFDIA,KLEV)
           
REAL(KIND=JPRB) :: ZTAU   (KIDIA:KFDIA,JPGPT,KLEV)

INTEGER(KIND=JPIM) :: JI, JLEV
INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "rrtm_taumol1.intfb.h"
#include "rrtm_taumol10.intfb.h"
#include "rrtm_taumol11.intfb.h"
#include "rrtm_taumol12.intfb.h"
#include "rrtm_taumol13.intfb.h"
#include "rrtm_taumol14.intfb.h"
#include "rrtm_taumol15.intfb.h"
#include "rrtm_taumol16.intfb.h"
#include "rrtm_taumol2.intfb.h"
#include "rrtm_taumol3.intfb.h"
#include "rrtm_taumol4.intfb.h"
#include "rrtm_taumol5.intfb.h"
#include "rrtm_taumol6.intfb.h"
#include "rrtm_taumol7.intfb.h"
#include "rrtm_taumol8.intfb.h"
#include "rrtm_taumol9.intfb.h"

IF (LHOOK) CALL DR_HOOK('RRTM_GAS_OPTICAL_DEPTH',0,ZHOOK_HANDLE)

!$ACC KERNELS DEFAULT(NONE) PRESENT(PFRAC) ASYNC(1)
PFRAC(:,:,:) = 0.0_jprb
!$ACC END KERNELS

ASSOCIATE(NFLEVG=>KLEV)

!$ACC DATA CREATE(ZTAU) PRESENT(POD)

CALL RRTM_TAUMOL1  (KIDIA,KFDIA,KLEV,ZTAU,PAVEL,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, PMINORFRAC, &
 & KINDMINOR,PSCALEMINORN2,PCOLBRD)  
CALL RRTM_TAUMOL2  (KIDIA,KFDIA,KLEV,ZTAU,PAVEL,PCOLDRY,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC)  
CALL RRTM_TAUMOL3  (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,PCOLN2O,PCOLDRY,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCO2, PRAT_H2OCO2_1,PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL4  (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,PCOLO3,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCO2, PRAT_H2OCO2_1, PRAT_O3CO2, PRAT_O3CO2_1)  
CALL RRTM_TAUMOL5  (KIDIA,KFDIA,KLEV,ZTAU,PWX,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,PCOLO3,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCO2, PRAT_H2OCO2_1, PRAT_O3CO2, PRAT_O3CO2_1,PMINORFRAC,KINDMINOR)   
CALL RRTM_TAUMOL6  (KIDIA,KFDIA,KLEV,ZTAU,PWX,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,PCOLCO2,PCOLDRY,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC,PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL7  (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLO3,PCOLCO2,PCOLDRY,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OO3, PRAT_H2OO3_1,PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL8  (KIDIA,KFDIA,KLEV,ZTAU,PWX,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,PCOLO3,PCOLN2O,PCOLCO2,PCOLDRY,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL9  (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLN2O,PCOLCH4,PCOLDRY,KLAYTROP,KLAYSWTCH,KLAYLOW,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCH4,PRAT_H2OCH4_1,PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL10 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC)  
CALL RRTM_TAUMOL11 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLH2O,PCOLO2,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC,PMINORFRAC,KINDMINOR,PSCALEMINOR)  
CALL RRTM_TAUMOL12 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCO2, PRAT_H2OCO2_1)  
CALL RRTM_TAUMOL13 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLN2O,PCOLCO2,PCOLO3,PCOLDRY,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2ON2O, PRAT_H2ON2O_1,PMINORFRAC,KINDMINOR)  
CALL RRTM_TAUMOL14 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,&
 & PCOLCO2,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC)  
CALL RRTM_TAUMOL15 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCO2,PCOLN2O,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_N2OCO2, PRAT_N2OCO2_1,PMINORFRAC,KINDMINOR,PSCALEMINOR,PCOLBRD)  
CALL RRTM_TAUMOL16 (KIDIA,KFDIA,KLEV,ZTAU,&
 & PTAUAERL,PFAC00,PFAC01,PFAC10,PFAC11,PFORFAC,PFORFRAC,KINDFOR,KJP,KJT,KJT1,PONEMINUS,&
 & PCOLH2O,PCOLCH4,KLAYTROP,PSELFFAC,PSELFFRAC,KINDSELF,PFRAC, &
 & PRAT_H2OCH4,PRAT_H2OCH4_1)   

!TO CHECK TOTAL OD FOR EACH BAND
    ! print*,'ZTAU2= ',sum(ZTAU(:,11:22,:),2)
    ! print*,'ZTAU3= ',sum(ZTAU(:,23:38,:),2)
    ! print*,'ZTAU4= ',sum(ZTAU(:,39:52,:),2)
    ! print*,'ZTAU5= ',sum(ZTAU(:,53:68,:),2)
    ! print*,'ZTAU6= ',sum(ZTAU(:,69:76,:),2)
    ! print*,'ZTAU7= ',sum(ZTAU(:,77:88,:),2)
    ! print*,'ZTAU8= ',sum(ZTAU(:,89:96,:),2)
    ! print*,'ZTAU9= ',sum(ZTAU(:,97:108,:),2)
    ! print*,'ZTAU10= ',sum(ZTAU(:,109:114,:),2)
    ! print*,'ZTAU11= ',sum(ZTAU(:,115:122,:),2)
    ! print*,'ZTAU12= ',sum(ZTAU(:,123:130,:),2)
    ! print*,'ZTAU13= ',sum(ZTAU(:,131:134,:),2)
    ! print*,'ZTAU14= ',sum(ZTAU(:,135:136,:),2)
    ! print*,'ZTAU15= ',sum(ZTAU(:,137:138,:),2)
    ! print*,'ZTAU16= ',sum(ZTAU(:,139:140,:),2)


!- Loop over g-channels.
!$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
!$ACC LOOP GANG VECTOR TILE(2, 64, 1)
DO JLEV = 1, KLEV
!cdir unroll=4
  DO JI = 1, JPGPT
    DO JLON = KIDIA, KFDIA
      POD(JI,JLEV,JLON) = ZTAU(JLON,JI,JLEV)
    ENDDO
  ENDDO
ENDDO
!$ACC END PARALLEL
!     -----------------------------------------------------------------

!$ACC WAIT
!$ACC END DATA

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('RRTM_GAS_OPTICAL_DEPTH',1,ZHOOK_HANDLE)

END SUBROUTINE RRTM_GAS_OPTICAL_DEPTH
