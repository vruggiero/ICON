! This file has been modified for the use in ICON

MODULE YOESRTA21

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA21* - SRTM COEFFICIENTS FOR INTERVAL 21
!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG21 = 16

REAL(KIND=JPRB) :: KA(9,5,13,JPG)   
REAL(KIND=JPRB) :: KB(5,5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(9,5,13,JPG)   
REAL(KIND=JPRD) :: KB_D(5,5,13:59,JPG)
REAL(KIND=JPRB) :: SELFREF(10,JPG)  ,FORREF(4,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG,9)
REAL(KIND=JPRB) :: RAYL              ,STRRAT
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(9,5,13,NG21)   ,ABSA(585,NG21)
REAL(KIND=JPRB) :: KBC(5,5,13:59,NG21),ABSB(1175,NG21)
REAL(KIND=JPRB) :: SELFREFC(10,NG21)  ,FORREFC(4,NG21)
REAL(KIND=JPRB) :: SFLUXREFC(NG21,9)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

!$ACC DECLARE CREATE(KAC, ABSA, KBC, ABSB, SELFREFC, FORREFC, SFLUXREFC)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SELFREF : REAL     self brodening coefficient for water vapour
! FORREF  : REAL     foreign broadening coefficient for water vapour
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! STRRAT  : REAL     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
END MODULE YOESRTA21

