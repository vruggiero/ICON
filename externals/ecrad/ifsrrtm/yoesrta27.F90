! This file has been modified for the use in ICON

MODULE YOESRTA27

USE PARKIND1  ,ONLY : JPIM     ,JPRB,JPRD

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
!     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPG = 16, NG27 = 16

REAL(KIND=JPRB) :: KA(5,13,JPG)   
REAL(KIND=JPRB) :: KB(5,13:59,JPG)
REAL(KIND=JPRD) :: KA_D(5,13,JPG)   
REAL(KIND=JPRD) :: KB_D(5,13:59,JPG)
REAL(KIND=JPRB) :: SFLUXREF(JPG)  ,RAYL(JPG)
REAL(KIND=JPRB) :: SCALEKUR
INTEGER(KIND=JPIM) :: LAYREFFR

REAL(KIND=JPRB) :: KAC(5,13,NG27)   ,ABSA(65,NG27)
REAL(KIND=JPRB) :: KBC(5,13:59,NG27),ABSB(235,NG27)
REAL(KIND=JPRB) :: SFLUXREFC(NG27)  ,RAYLC(NG27)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
EQUIVALENCE (KAC(1,1,1),ABSA(1,1)), (KBC(1,13,1),ABSB(1,1))

!$ACC DECLARE CREATE(KAC, ABSA, KBC, ABSB, SFLUXREFC, RAYLC)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL
! KA      : REAL     absorption coefficient of major absorber
! KB      : REAL     absorption coefficient of secondary absorber
! SFLUXREF: REAL     Incident solar radiation in the spectral interval
! RAYL    : REAL     Rayleigh scattering parameter
! SCALEKUR: REAL     weighting factor to account for Kurucz's correction to solar spectrum 
! LAYREFFR: INTEGER  reference level for the transition
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
END MODULE YOESRTA27

