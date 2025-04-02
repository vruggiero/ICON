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

MODULE radar_sphbessel_vec

!------------------------------------------------------------------------------
!
! Description:  This module provides basic subroutines and functions for the
!               computation of spherical Bessel-functions, which
!               are used in the Karlsruhe Libraries for Radar Reflectivity
!               Calculations based on Mie-Theory.
!
!------------------------------------------------------------------------------
!
! Declarations:
!
! Modules used:
!

  USE radar_kind, ONLY : dp

  !=====================================================================
  !
  ! The following module radar_sphbessel_vec provides
  ! spherical Bessel-functionen 1. and 2. kind of complexe
  ! arguments and of integer order and there derivatives
  !
  !=====================================================================

  IMPLICIT NONE

  PUBLIC

  PRIVATE MSTA1VEC, MSTA2VEC, ENVJ, ENVJVEC

CONTAINS

!=====================================================================
!
! some help functions for calculation of spherical
! Bessel-functionen in CSPHJY() (see below)
!
!=====================================================================

  FUNCTION ENVJVEC(N,X,anz) RESULT (Y)
    INTEGER :: anz
    INTEGER :: N(anz)
    REAL(KIND=dp) :: X(anz)
    REAL(KIND=dp) :: Y(anz)
    Y = 0.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
    RETURN
  END FUNCTION ENVJVEC

  FUNCTION ENVJ(N,X) RESULT (Y)
    INTEGER :: N
    REAL(KIND=dp) :: X
    REAL(KIND=dp) :: Y
    Y = 0.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
    RETURN
  END FUNCTION ENVJ

  FUNCTION MSTA1VEC(X,MP,anz) RESULT (Y)

    !       ===================================================
    !       Purpose: Determine the starting point for backward
    !                recurrence such that the magnitude of
    !                Jn(x) at that point is about 10^(-MP)
    !       Input :  x     --- Argument of Jn(x)
    !                MP    --- Value of magnitude
    !       Output:  MSTA1 --- Starting point
    !       ===================================================

    IMPLICIT NONE

    INTEGER :: anz
    INTEGER :: MP
    REAL(KIND=dp) :: X(anz)
    REAL(KIND=dp) :: Y(anz)
    REAL(KIND=dp) :: A0(anz), F(anz), F0(anz), F1(anz)
    INTEGER :: N0(anz), N1(anz), IT, NN(anz), K

    A0 = ABS(X)
    N0 = INT(1.1*A0) + 1
    F0 = ENVJVEC(N0,A0,anz) - MP
    N1 = N0 + 5
    F1 = ENVJVEC(N1,A0,anz) - MP
    NN = N1-(N1-N0)/(1.0D0-F0/F1)
    F  = ENVJVEC(NN,A0,anz)-MP
    DO IT=1,20
      DO k=1,anz
        IF(ABS(NN(k)-N1(k)) >= 1) THEN
          N0(k)=N1(k)
          F0(k)=F1(k)
          N1(k)=NN(k)
          F1(k)=F(k)
          NN(k)=N1(k)-(N1(k)-N0(k))/(1.0D0-F0(k)/F1(k))
          F(k) =ENVJ(NN(k),A0(k))-MP
        END IF
      END DO
    END DO
    Y=NN

    RETURN
  END FUNCTION MSTA1VEC

  FUNCTION MSTA2VEC(X,N,MP,anz) RESULT (Y)

    !       ===================================================
    !       Purpose: Determine the starting point for backward
    !                recurrence such that all Jn(x) has MP
    !                significant digits
    !       Input :  x  --- Argument of Jn(x)
    !                n  --- Order of Jn(x)
    !                MP --- Significant digit
    !       Output:  MSTA2 --- Starting point
    !       ===================================================

    IMPLICIT NONE

    INTEGER :: anz
    INTEGER :: N, MP
    REAL(KIND=dp) :: X(anz)
    REAL(KIND=dp) :: Y(anz)

    INTEGER :: N0(anz), N1(anz), IT, NN(anz), K
    REAL(KIND=dp) :: EJN, OBJ(anz), A0(anz), HMP, F(anz), F0(anz), F1(anz)

    A0 = ABS(X)
    HMP = 0.5D0*MP
    DO K=1,anz
      EJN = ENVJ(N,A0(K))
      IF (EJN < HMP) THEN
        OBJ(K) = MP
        N0(K) = INT(1.1*A0(K))
      ELSE
        OBJ(K) = HMP + EJN
        N0(K) = N
      END IF
    END DO
    F0 = ENVJVEC(N0,A0,anz)-OBJ
    N1 = N0 + 5
    F1 = ENVJVEC(N1,A0,anz) - OBJ
    NN = N1-(N1-N0)/(1.0D0-F0/F1)
    F = ENVJVEC(NN,A0,anz) - OBJ
    DO IT=1,20
      DO k=1,anz
        IF (ABS(NN(k)-N1(k)) >= 1) THEN
          N0(k)=N1(k)
          F0(k)=F1(k)
          N1(k)=NN(k)
          F1(k)=F(k)
          NN(k) = N1(k)-(N1(k)-N0(k))/(1.0D0-F0(k)/F1(k))
          F(k) = ENVJ(NN(k),A0(k)) - OBJ(k)
        END IF
      END DO
    END DO
    Y = NN + 10

    RETURN
  END FUNCTION MSTA2VEC

  SUBROUTINE CSPHJYVEC(N,Z,NM,CSJ,CDJ,CSY,CDY,anz)
!*********************************************************************
!
! Vectorized version!
!
!*********************************************************************
    !       ==========================================================
    !       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
    !                and their derivatives for a complex argument up
    !                to order N.
    !       Input :  z --- Complex argument
    !                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
    !
    !       Output:  Vectors CSJ(0:N), CDJ(0:N),CSY(0:N),CDY(0:N):
    !                CSJ(n) --- jn(z)
    !                CDJ(n) --- jn'(z)
    !                CSY(n) --- yn(z)
    !                CDY(n) --- yn'(z)
    !                NM --- Highest order computed
    !                --> Higher orders are virtually 0.0 and need not
    !                be computed!!!
    !       Routines called:
    !                MSTA1 and MSTA2 for computing the starting
    !                point for backward recurrence
    !
    !       From: "Computation of Special Functions",  Shanjie Zhang, Jianming Jin, 1996, Wiley,
    !              ISBN: 0-471-11963-6
    !
    !       ==========================================================

    IMPLICIT NONE

    INTEGER :: anz
    INTEGER :: N, NM(anz)
    COMPLEX(kind=dp) :: CSJ(anz,0:N),CDJ(anz,0:N),CSY(anz,0:N),CDY(anz,0:N), Z(anz)

    COMPLEX(kind=dp) :: CSA(anz), CSB(anz), CF0(anz), CF1(anz), CF(anz), CS(anz)

    INTEGER :: K, M(anz), J, Mhelp(anz)

    REAL(KIND=dp) :: A0(anz)

    A0=ABS(Z)
    NM=N

    ! Preset values for abs(A0) < 1d-60, which will
    ! be overwritten later in case of abs(A0) >= 1d-60.
    CSJ(:,0:N)=0.0D0
    CDJ(:,0:N)=0.0D0
    CSY(:,0:N)=-1.0D+300
    CDY(:,0:N)=1.0D+300
    CSJ(:,0)=(1.0D0,0.0D0)
    CDJ(:,1)=(0.333333333333333D0,0.0D0)

    DO j=1,anz
      IF (A0(j) >= 1.0D-60) THEN
        CSJ(j,0)=SIN(Z(j))/Z(j)
        CSJ(j,1)=(CSJ(j,0)-COS(Z(j)))/Z(j)
      END IF
    END DO

    IF (N >= 2) THEN
      CSA=CSJ(:,0)
      CSB=CSJ(:,1)
      M=MSTA1VEC(A0,200,anz)
      WHERE (M < N)
        NM=M
      END WHERE
      Mhelp=MSTA2VEC(A0,N,15,anz)
      WHERE (M >= N)
        M = Mhelp
      END WHERE

      CF0=0.0D0
      CF1=1.0D-100
      DO K=MAXVAL(M),0,-1
        DO j=1,anz
          IF (A0(j) >= 1.0D-60 .AND. K <= M(j)) THEN
            CF(j)=(2.0D0*K+3.0D0)*CF1(j)/Z(j)-CF0(j)
            CF0(j)=CF1(j)
            CF1(j)=CF(j)
          END IF
          IF (A0(j) >= 1.0D-60 .AND. K <= M(j) .AND. K <= NM(j)) THEN
            CSJ(j,K)=CF(j)
          END IF
        END DO
      END DO
      WHERE (ABS(CSA) > ABS(CSB))
        CS=CSA/CF
      ELSEWHERE
        CS=CSB/CF0
      END WHERE
      DO K=0,MAXVAL(NM)
        DO j=1,anz
          IF (A0(j) >= 1.0D-60 .AND. K <= NM(j)) THEN
            CSJ(j,K)=CS(j)*CSJ(j,K)
          END IF
        END DO
      END DO
    ENDIF


    DO j=1,anz
      IF (A0(j) >= 1.0D-60) THEN
        CDJ(j,0)=(COS(Z(j))-SIN(Z(j))/Z(j))/Z(j)
      END IF
    END DO
    DO K=1,MAXVAL(NM)
      DO j=1,anz
        IF (A0(j) >= 1.0D-60 .AND. K <= NM(j)) THEN
          CDJ(j,K)=CSJ(j,K-1)-(K+1.0D0)*CSJ(j,K)/Z(j)
        END IF
      END DO
    END DO
    DO j=1,anz
      IF (A0(j) >= 1.0D-60) THEN
        CSY(j,0)=-COS(Z(j))/Z(j)
        CSY(j,1)=(CSY(j,0)-SIN(Z(j)))/Z(j)
        CDY(j,0)=(SIN(Z(j))+COS(Z(j))/Z(j))/Z(j)
        CDY(j,1)=(2.0D0*CDY(j,0)-COS(Z(j)))/Z(j)
      END IF
    END DO

    DO K=2,MAXVAL(NM)
      DO j=1,anz
        IF (A0(j) >= 1.0D-60 .AND. K <= NM(j) .AND. ABS(CSJ(j,K-1)) > ABS(CSJ(j,K-2))) THEN
          CSY(j,K)=(CSJ(j,K)*CSY(j,K-1)-1.0D0/(Z(j)*Z(j)))/CSJ(j,K-1)
        ELSE IF (A0(j) >= 1.0D-60 .AND. K <= NM(j)) THEN
          CSY(j,K)=(CSJ(j,K)*CSY(j,K-2)-(2.0D0*K-1.0D0)/Z(j)**3)/CSJ(j,K-2)
        ENDIF
      END DO
    END DO
    DO K=2,MAXVAL(NM)
      DO j=1,anz
        IF (A0(j) >= 1.0D-60 .AND. K <= NM(j)) THEN
          CDY(j,K)=CSY(j,K-1)-(K+1.0D0)*CSY(j,K)/Z(j)
        END IF
      END DO
    END DO

    RETURN
  END SUBROUTINE CSPHJYVEC

  INTEGER FUNCTION MSTA1(X,MP)

    !       ===================================================
    !       Purpose: Determine the starting point for backward
    !                recurrence such that the magnitude of
    !                Jn(x) at that point is about 10^(-MP)
    !       Input :  x     --- Argument of Jn(x)
    !                MP    --- Value of magnitude
    !       Output:  MSTA1 --- Starting point
    !       ===================================================

    IMPLICIT NONE

    REAL(KIND=dp) :: X
    INTEGER :: MP
    REAL(KIND=dp) :: A0, F, F0, F1
    INTEGER :: N0, N1, IT, NN

    A0 = ABS(X)
    N0 = INT(1.1*A0) + 1
    F0 = ENVJ(N0,A0) - MP
    N1 = N0 + 5
    F1 = ENVJ(N1,A0) - MP
    DO IT=1,20
      NN=N1-(N1-N0)/(1.0D0-F0/F1)
      F=ENVJ(NN,A0)-MP
      IF(ABS(NN-N1) < 1) exit
      N0=N1
      F0=F1
      N1=NN
      F1=F
    END DO
    MSTA1=NN

    RETURN
  END FUNCTION MSTA1

  INTEGER FUNCTION MSTA2(X,N,MP)

    !       ===================================================
    !       Purpose: Determine the starting point for backward
    !                recurrence such that all Jn(x) has MP
    !                significant digits
    !       Input :  x  --- Argument of Jn(x)
    !                n  --- Order of Jn(x)
    !                MP --- Significant digit
    !       Output:  MSTA2 --- Starting point
    !       ===================================================

    IMPLICIT NONE

    INTEGER :: N, MP
    REAL(KIND=dp) :: X

    INTEGER :: N0, N1, IT, NN
    REAL(KIND=dp) :: EJN, OBJ, A0, HMP, F, F0, F1

    A0 = ABS(X)
    HMP = 0.5D0*MP
    EJN = ENVJ(N,A0)
    IF (EJN < HMP) THEN
      OBJ = MP
      N0 = INT(1.1*A0)
    ELSE
      OBJ = HMP + EJN
      N0 = N
    ENDIF
    F0 = ENVJ(N0,A0)-OBJ
    N1 = N0 + 5
    F1 = ENVJ(N1,A0) - OBJ
    DO IT=1,20
      NN = N1-(N1-N0)/(1.0D0-F0/F1)
      F = ENVJ(NN,A0) - OBJ
      IF (ABS(NN-N1) < 1) EXIT
      N0=N1
      F0=F1
      N1=NN
      F1=F
    END DO
    MSTA2 = NN + 10

    RETURN
  END FUNCTION MSTA2

  SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)

    !       ==========================================================
    !       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
    !                and their derivatives for a complex argument up
    !                to order N.
    !       Input :  z --- Complex argument
    !                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
    !
    !       Output:  Vectors CSJ(0:N), CDJ(0:N),CSY(0:N),CDY(0:N):
    !                CSJ(n) --- jn(z)
    !                CDJ(n) --- jn'(z)
    !                CSY(n) --- yn(z)
    !                CDY(n) --- yn'(z)
    !                NM --- Highest order computed
    !                --> Higher orders are virtually 0.0 and need not
    !                be computed!!!
    !       Routines called:
    !                MSTA1 and MSTA2 for computing the starting
    !                point for backward recurrence
    !
    !       From: "Computation of Special Functions",  Shanjie Zhang, Jianming Jin, 1996, Wiley,
    !              ISBN: 0-471-11963-6
    !
    !       ==========================================================

    IMPLICIT NONE

    INTEGER :: N, NM
    COMPLEX(kind=dp) :: CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N), Z

    COMPLEX(kind=dp) :: CSA, CSB, CF0, CF1, CF, CS

    INTEGER :: K, M

    REAL(KIND=dp) :: A0

    A0=ABS(Z)
    NM=N
    IF (A0 < 1.0D-60) THEN
      DO  K=0,N
        CSJ(K)=0.0D0
        CDJ(K)=0.0D0
        CSY(K)=-1.0D+300
        CDY(K)=1.0D+300
      END DO
      CSJ(0)=(1.0D0,0.0D0)
      CDJ(1)=(0.333333333333333D0,0.0D0)
      RETURN
    ENDIF

    CSJ(0)=SIN(Z)/Z
    CSJ(1)=(CSJ(0)-COS(Z))/Z
    IF (N >= 2) THEN
      CSA=CSJ(0)
      CSB=CSJ(1)
      M=MSTA1(A0,200)
      IF (M < N) THEN
        NM=M
      ELSE
        M=MSTA2(A0,N,15)
      ENDIF
      CF0=0.0D0
      CF1=1.0D-100
      DO K=M,0,-1
        CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
        IF (K <= NM) CSJ(K)=CF
        CF0=CF1
        CF1=CF
      END DO
      IF (ABS(CSA) > ABS(CSB)) THEN
        CS=CSA/CF
      ELSE
        CS=CSB/CF0
      END IF
      DO K=0,NM
        CSJ(K)=CS*CSJ(K)
      END DO
    ENDIF
    CDJ(0)=(COS(Z)-SIN(Z)/Z)/Z
    DO K=1,NM
      CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
    END DO
    CSY(0)=-COS(Z)/Z
    CSY(1)=(CSY(0)-SIN(Z))/Z
    CDY(0)=(SIN(Z)+COS(Z)/Z)/Z
    CDY(1)=(2.0D0*CDY(0)-COS(Z))/Z
    DO K=2,NM
      IF (ABS(CSJ(K-1)) > ABS(CSJ(K-2))) THEN
        CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
      ELSE
        CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
      ENDIF
    END DO
    DO K=2,NM
      CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
    END DO

    RETURN
  END SUBROUTINE CSPHJY

END MODULE radar_sphbessel_vec
