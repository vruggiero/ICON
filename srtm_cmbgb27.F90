! This file has been modified for the use in ICON

SUBROUTINE SRTM_CMBGB27

!     BAND 27:  29000-38000 cm-1 (low - O3; high - O3)
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM , JPRB
USE ecradhook   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOESRTM  , ONLY : NGN
USE YOESRTWN , ONLY : NGC, NGS, RWGT
!USE YOESRTWN , ONLY : NGC, NGS, NGN, RWGT
USE YOESRTA27, ONLY : KA, KB, SFLUXREF, RAYL, &
                    & KAC, KBC, SFLUXREFC, RAYLC

IMPLICIT NONE

! Local variables
INTEGER(KIND=JPIM) :: JT, JP, IGC, IPR, IPRSM
REAL(KIND=JPRB)    :: ZSUMK, ZSUMF1, ZSUMF2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB27',0,ZHOOK_HANDLE)

DO JT = 1,5
  DO JP = 1,13
    IPRSM = 0
    DO IGC = 1,NGC(12)
      ZSUMK = 0.
      DO IPR = 1, NGN(NGS(11)+IGC)
        IPRSM = IPRSM + 1
        ZSUMK = ZSUMK + KA(JT,JP,IPRSM)*RWGT(IPRSM+176)
      ENDDO
      KAC(JT,JP,IGC) = ZSUMK
    ENDDO
  ENDDO

  DO JP = 13,59
    IPRSM = 0
    DO IGC = 1,NGC(12)
      ZSUMK = 0.
      DO IPR = 1, NGN(NGS(11)+IGC)
        IPRSM = IPRSM + 1
        ZSUMK = ZSUMK + KB(JT,JP,IPRSM)*RWGT(IPRSM+176)
      ENDDO
      KBC(JT,JP,IGC) = ZSUMK
    ENDDO
  ENDDO
ENDDO

IPRSM = 0
DO IGC = 1,NGC(12)
  ZSUMF1 = 0.
  ZSUMF2 = 0.
  DO IPR = 1, NGN(NGS(11)+IGC)
    IPRSM = IPRSM + 1
    ZSUMF1 = ZSUMF1 + SFLUXREF(IPRSM)
    ZSUMF2 = ZSUMF2 + RAYL(IPRSM)*RWGT(IPRSM+176)
  ENDDO
  SFLUXREFC(IGC) = ZSUMF1
  RAYLC(IGC) = ZSUMF2
ENDDO

!$ACC UPDATE DEVICE(KAC, KBC, SFLUXREFC, RAYLC)

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_CMBGB27',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_CMBGB27

