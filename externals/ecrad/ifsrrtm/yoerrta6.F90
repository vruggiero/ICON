! This file has been modified for the use in ICON

MODULE YOERRTA6

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
!     ABozzo 201306 updaten to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG6  = 8

REAL(KIND=JPRB) , DIMENSION(NG6) :: FRACREFA

REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC11ADJ
REAL(KIND=JPRB) , DIMENSION(NG6) :: CFC12


REAL(KIND=JPRB) :: KA(5,13,NG6),ABSA(65,NG6)
REAL(KIND=JPRB) :: SELFREF(10,NG6)
REAL(KIND=JPRB) :: KA_MCO2(19,NG6)
REAL(KIND=JPRB) :: FORREF(4,NG6)

EQUIVALENCE (KA(1,1,1),ABSA(1,1))

!$ACC DECLARE CREATE(FRACREFA, CFC11ADJ, CFC12, KA, ABSA, SELFREF, KA_MCO2, &
!$ACC                FORREF)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE *

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSCO2  : REAL     absorption coefficient for CO2
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! CFC11ADJ: REAL     absorption coefficient for CFC-11 (adjusted)
! CFC12   : REAL     absorption coefficient for CFC-12
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! SELFREF : REAL     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
END MODULE YOERRTA6
