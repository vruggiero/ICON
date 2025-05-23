! This file has been modified for the use in ICON

MODULE YOERRTA16

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA16* - RRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG16 = 2

REAL(KIND=JPRB) :: FRACREFA(NG16,9)
REAL(KIND=JPRB) , DIMENSION(NG16) :: FRACREFB

REAL(KIND=JPRB) :: KA(9,5,13,NG16) ,ABSA(585,NG16)
REAL(KIND=JPRB) :: KB(5,13:59,NG16), ABSB(235,NG16)
REAL(KIND=JPRB) :: SELFREF(10,NG16)
REAL(KIND=JPRB) :: FORREF(4,NG16)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))

!$ACC DECLARE CREATE(FRACREFA, FRACREFB, KA, ABSA, KB, ABSB, SELFREF, FORREF)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRAT  : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA16
