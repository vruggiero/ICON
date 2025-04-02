!
!+ Some public domain SLATEC routines adapted to Fortran 90
!
MODULE slatec_module
!
! Description:
!   Some public domain SLATEC routines adapted to Fortran 90
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  public: xrs  ! eigenvalues + vectors
! V1_5         2009/05/25 Harald Anlauf
!  subroutine xRS_opt (optimized for NEC SX).
!  sort: replace dpsort, ipsort by a modified quicksort
! V1_8         2009/12/09 Harald Anlauf
!  sortix, sortrx: add workaround for xlf bug
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  remove unused variables
! V1_19        2012-04-16 Harald Anlauf
!  revert workaround (rev 5366) for xlf v12.1, it produces bad results
! V1_22        2013-02-13 Harald Anlauf
!  new include file (Template for Hybrid QuickSort), replace cpsort by cpsortx
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2000  original routines gathered
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_kind   ,only: wp, sp, dp, i8
  use mo_system ,only: abort
  implicit none

  !----------------
  ! Public entities
  !----------------
  private
  public :: rs      ! eigenvalues + vectors
  public :: xrs     ! eigenvalues + vectors
  public :: xrs_opt ! eigenvalues + vectors (opt. version for SX9)
  public :: sort    ! sorting
!------------------------------------------------------------------------------
  interface sort
    module procedure spsort             ! SLATEC quicksort for default real
!   module procedure spsortx            ! Hybrid quicksort for default real
#if defined(__ICON__)
    module procedure cpsort             ! SLATEC quicksort for character(*)
!   module procedure cpsortx            ! INTERNAL COMPILER ERROR ???
#else
!   module procedure cpsort             ! SLATEC quicksort for character(*)
    module procedure cpsortx            ! Hybrid quicksort for character(*)
#endif
#if defined(__ibm__)
!   module procedure dpsort,  ipsort    ! SLATEC quicksort
    module procedure dpsortx, ipsortx   ! Modified quicksort (xlf V12.1 bug)
#else
!   module procedure dpsort,  ipsort    ! SLATEC quicksort
    module procedure dpsortx, ipsortx   ! Modified quicksort
#endif
    module procedure i8sortx            ! Hybrid quicksort for integer(8)
  end interface
!------------------------------------------------------------------------------
  interface rmach
    module procedure spmach, dpmach
  end interface
!==============================================================================
contains
!==============================================================================
  subroutine rs (a,w,z,ierr)
  real(wp) ,intent(in)            :: a(:,:) ! real symmetric matrix
  real(wp) ,intent(out)           :: w(:)   ! eigenvalues
  real(wp) ,intent(out) ,optional :: z(:,:) ! eigenvectors
  integer  ,intent(out) ,optional :: ierr   !  0=ok
                                            ! =j: jth eigenvalue not determined
  !-----------------------------------------------------------------------
  ! Fortran90 Wrapper for SLATEC routine rs
  !
  ! This subroutine calls the recommended sequence of
  ! subroutines from the eigensystem subroutine package (EISPACK)
  ! to find the eigenvalues and eigenvectors (if desired)
  ! of a REAL SYMMETRIC matrix.
  !
  ! On Input
  !
  !    A contains the real symmetric matrix.  A is a two-dimensional
  !      REAL array, dimensioned A(N,N).
  !
  ! On Output
  !
  !    A is unaltered.
  !
  !    W contains the eigenvalues in ascending order.  W is a one-
  !      dimensional REAL array, dimensioned W(N).
  !
  !    Z contains the eigenvectors if present.  The
  !      eigenvectors are orthonormal.  Z is a two-dimensional
  !      REAL array, dimensioned Z(NM,N).
  !
  !    IERR is an INTEGER flag set to
  !      Zero       for normal return,
  !      J          if the J-th eigenvalue has not been
  !                 determined after 30 iterations.
  !                 The eigenvalues, and eigenvectors if requested,
  !                 should be correct for indices 1, 2, ..., IERR-1.
  !
  !      The routine aborts if an error occures and IERR is not present
  !-----------------------------------------------------------------------
    integer               :: n, ie
    real(wp) ,allocatable :: zz(:,:), fv1(:), fv2(:)

    n = size(a,1)
    if (    size(a,2) /= n) call abort ('rs: a')
    if (    size(w)   /= n) call abort ('rs: w')
    allocate (fv1(n))
    allocate (fv2(n))
    if(present(z)) then
      if (any(shape(z)/= n)) call abort ('rs: z')
      call xrs (n,n,a,w,1,z,fv1,fv2,ie)
    else
      allocate (zz(n,n))
      call xrs (n,n,a,w,0,zz,fv1,fv2,ie)
      deallocate (zz)
    endif
    if(present(ierr)) then
      ierr = ie
    else
      if(ie/=0) call abort ('rs: ierr')
    endif
    deallocate (fv1)
    deallocate (fv2)
  end subroutine rs
!==============================================================================
      SUBROUTINE xRS (NM, N, A, W, MATZ, Z, FV1, FV2, IERR)
!***BEGIN PROLOGUE  RS
!***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
!            of a real symmetric matrix.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A1
!***TYPE      SINGLE PRECISION (RS-S, CH-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a REAL SYMMETRIC matrix.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the real symmetric matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        A is unaltered.
!
!        W contains the eigenvalues in ascending order.  W is a one-
!          dimensional REAL array, dimensioned W(N).
!
!        Z contains the eigenvectors if MATZ is not zero.  The
!          eigenvectors are orthonormal.  Z is a two-dimensional
!          REAL array, dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues, and eigenvectors if requested,
!                     should be correct for indices 1, 2, ..., IERR-1.
!
!        FV1 and FV2 are one-dimensional REAL arrays used for temporary
!          storage, dimensioned FV1(N) and FV2(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  TQL2, TQLRAT, TRED1, TRED2
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RS
!
      INTEGER N,NM,IERR,MATZ
      REAL(wp) A(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)
!
!***FIRST EXECUTABLE STATEMENT  RS
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10 IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END subroutine xrs

!------------------------------------------------------------------------------
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
!***BEGIN PROLOGUE  TQL2
!***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
!            tridiagonal matrix.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (TQL2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQL2,
!     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!     Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!     The eigenvectors of a FULL SYMMETRIC matrix can also
!     be found if  TRED2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        Z contains the transformation matrix produced in the
!          reduction by  TRED2, if performed.  If the eigenvectors
!          of the tridiagonal matrix are desired, Z must contain
!          the identity matrix.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., IERR-1.
!
!        E has been destroyed.
!
!        Z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  If an error exit is made,
!          Z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TQL2
!
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      REAL(WP) D(*),E(*),Z(NM,*)
      REAL(WP) B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
!
!***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO I = 2, N
        E(I-1) = E(I)
      end do
!
      F = 0.0_wp
      B = 0.0_wp
      E(N) = 0.0_wp
!
      DO 240 L = 1, N
        J = 0
        H = ABS(D(L)) + ABS(E(L))
        IF (B .LT. H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
        DO 110 M = L, N
          IF (B + ABS(E(M)) .EQ. B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110   CONTINUE
!
  120   IF (M .EQ. L) GO TO 220
  130   IF (J .EQ. 30) GO TO 1000
        J = J + 1
!     .......... FORM SHIFT ..........
        L1 = L + 1
        L2 = L1 + 1
        G = D(L)
        P = (D(L1) - G) / (2.0_wp * E(L))
        R = PYTHAG(P,1.0_wp)
        D(L) = E(L) / (P + SIGN(R,P))
        D(L1) = E(L) * (P + SIGN(R,P))
        DL1 = D(L1)
        H = G - D(L)
        IF (L2 .GT. N) GO TO 145
!
        DO I = L2, N
          D(I) = D(I) - H
        end do
!
  145   F = F + H
!     .......... QL TRANSFORMATION ..........
        P = D(M)
        C = 1.0_wp
        C2 = C
        EL1 = E(L1)
        S = 0.0_wp
        MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
        DO 200 II = 1, MML
          C3 = C2
          C2 = C
          S2 = S
          I = M - II
          G = C * E(I)
          H = C * P
          IF (ABS(P) .LT. ABS(E(I))) GO TO 150
          C = E(I) / P
          R = SQRT(C*C+1.0_wp)
          E(I+1) = S * P * R
          S = C / R
          C = 1.0_wp / R
          GO TO 160
  150     C = P / E(I)
          R = SQRT(C*C+1.0_wp)
          E(I+1) = S * E(I) * R
          S = 1.0_wp / R
          C = C * S
  160     P = C * D(I) - S * G
          D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
          DO 180 K = 1, N
            H = Z(K,I+1)
            Z(K,I+1) = S * Z(K,I) + C * H
            Z(K,I) = C * Z(K,I) - S * H
  180     CONTINUE
!
  200   CONTINUE
!
        P = -S * S2 * C3 * EL1 * E(L) / DL1
        E(L) = S * P
        D(L) = C * P
        IF (B + ABS(E(L)) .GT. B) GO TO 130
  220   D(L) = D(L) + F
  240 END DO
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
        I = II - 1
        K = I
        P = D(I)
!
        DO 260 J = II, N
          IF (D(J) .GE. P) GO TO 260
          K = J
          P = D(J)
  260   CONTINUE
!
        IF (K .EQ. I) GO TO 300
        D(K) = D(I)
        D(I) = P
!
        DO 280 J = 1, N
          P = Z(J,I)
          Z(J,I) = Z(J,K)
          Z(J,K) = P
  280   CONTINUE
!
  300 END DO
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END subroutine tql2
!------------------------------------------------------------------------------
      SUBROUTINE TRED1 (NM, N, A, D, E, E2)
!***BEGIN PROLOGUE  TRED1
!***PURPOSE  Reduce a real symmetric matrix to symmetric tridiagonal
!            matrix using orthogonal similarity transformations.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (TRED1-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED1,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a REAL SYMMETRIC matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the real symmetric input matrix.  Only the lower
!          triangle of the matrix need be supplied.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!     On Output
!
!        A contains information about the orthogonal transformations
!          used in the reduction in its strict lower triangle.  The
!          full upper triangle of A is unaltered.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is set
!          to zero.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2 may coincide with E if the squares are not needed.
!          E2 is a one-dimensional REAL array, dimensioned E2(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TRED1
!
      INTEGER I,J,K,L,N,II,NM,JP1
      real(wp) A(NM,*),D(*),E(*),E2(*)
      real(wp) F,G,H,SCALE
!
!***FIRST EXECUTABLE STATEMENT  TRED1
      DO I = 1, N
        D(I) = A(I,I)
      end do
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
        I = N + 1 - II
        L = I - 1
        H = 0.0_wp
        SCALE = 0.0_wp
        IF (L .LT. 1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
        DO K = 1, L
          SCALE = SCALE + ABS(A(I,K))
        end do
!
        IF (SCALE .NE. 0.0_wp) GO TO 140
  130   E(I) = 0.0_wp
        E2(I) = 0.0_wp
        GO TO 290
!
  140   DO 150 K = 1, L
          A(I,K) = A(I,K) / SCALE
          H = H + A(I,K) * A(I,K)
  150   CONTINUE
!
        E2(I) = SCALE * SCALE * H
        F = A(I,L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE * G
        H = H - F * G
        A(I,L) = F - G
        IF (L .EQ. 1) GO TO 270
        F = 0.0_wp
!
        DO 240 J = 1, L
          G = 0.0_wp
!     .......... FORM ELEMENT OF A*U ..........
          DO K = 1, J
            G = G + A(J,K) * A(I,K)
          end do
!
          JP1 = J + 1
          IF (L .LT. JP1) GO TO 220
!
          DO K = JP1, L
            G = G + A(K,J) * A(I,K)
          end do
!     .......... FORM ELEMENT OF P ..........
  220     E(J) = G / H
          F = F + E(J) * A(I,J)
  240   CONTINUE
!
        H = F / (H + H)
!     .......... FORM REDUCED A ..........
        DO J = 1, L
          F = A(I,J)
          G = E(J) - H * F
          E(J) = G
!
          DO K = 1, J
            A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
          end do
        end do
!
  270   DO K = 1, L
          A(I,K) = SCALE * A(I,K)
        end do
!
  290   H = D(I)
        D(I) = A(I,I)
        A(I,I) = H
  300 END DO
!
      RETURN
      END subroutine tred1
!------------------------------------------------------------------------------
      SUBROUTINE TRED2 (NM, N, A, D, E, Z)
!***BEGIN PROLOGUE  TRED2
!***PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal
!            matrix using and accumulating orthogonal transformations.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B1
!***TYPE      SINGLE PRECISION (TRED2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRED2,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a REAL SYMMETRIC matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the real symmetric input matrix.  Only the lower
!          triangle of the matrix need be supplied.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!     On Output
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is set
!          to zero.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        Z contains the orthogonal transformation matrix produced in
!          the reduction.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!        A and Z may coincide.  If distinct, A is unaltered.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TRED2
!
      INTEGER I,J,K,L,N,II,NM,JP1
      real(wp) A(NM,*),D(*),E(*),Z(NM,*)
      real(wp) F,G,H,HH,SCALE
!
!***FIRST EXECUTABLE STATEMENT  TRED2
      DO I = 1, N
!
        DO J = 1, I
          Z(I,J) = A(I,J)
        end do
      end do
!
      IF (N .EQ. 1) GO TO 320
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
        I = N + 2 - II
        L = I - 1
        H = 0.0_wp
        SCALE = 0.0_wp
        IF (L .LT. 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
        DO K = 1, L
          SCALE = SCALE + ABS(Z(I,K))
        end do
!
        IF (SCALE .NE. 0.0_wp) GO TO 140
  130   E(I) = Z(I,L)
        GO TO 290
!
  140   DO 150 K = 1, L
          Z(I,K) = Z(I,K) / SCALE
          H = H + Z(I,K) * Z(I,K)
  150   CONTINUE
!
        F = Z(I,L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE * G
        H = H - F * G
        Z(I,L) = F - G
        F = 0.0_wp
!
        DO 240 J = 1, L
          Z(J,I) = Z(I,J) / H
          G = 0.0_wp
!     .......... FORM ELEMENT OF A*U ..........
          DO K = 1, J
            G = G + Z(J,K) * Z(I,K)
          end do
!
          JP1 = J + 1
          IF (L .LT. JP1) GO TO 220
!
          DO K = JP1, L
            G = G + Z(K,J) * Z(I,K)
          end do
!     .......... FORM ELEMENT OF P ..........
  220     E(J) = G / H
          F = F + E(J) * Z(I,J)
  240   CONTINUE
!
        HH = F / (H + H)
!     .......... FORM REDUCED A ..........
        DO J = 1, L
          F = Z(I,J)
          G = E(J) - HH * F
          E(J) = G
!
          DO K = 1, J
            Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
          end do
        end do
!
  290   D(I) = H
  300 END DO
!
  320 D(1) = 0.0_wp
      E(1) = 0.0_wp
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
        L = I - 1
        IF (D(I) .EQ. 0.0_wp) GO TO 380
!
        DO J = 1, L
          G = 0.0_wp
!
          DO K = 1, L
            G = G + Z(I,K) * Z(K,J)
          end do
!
          DO K = 1, L
            Z(K,J) = Z(K,J) - G * Z(K,I)
          end do
        end do
!
  380   D(I) = Z(I,I)
        Z(I,I) = 1.0_wp
        IF (L .LT. 1) GO TO 500
!
        DO 400 J = 1, L
          Z(I,J) = 0.0_wp
          Z(J,I) = 0.0_wp
  400   CONTINUE
!
  500 END DO
!
      RETURN
      END subroutine tred2
!------------------------------------------------------------------------------
      SUBROUTINE TQLRAT (N, D, E2, IERR)
!***BEGIN PROLOGUE  TQLRAT
!***PURPOSE  Compute the eigenvalues of symmetric tridiagonal matrix
!            using a rational variant of the QL method.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (TQLRAT-S)
!***KEYWORDS  EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX, EISPACK,
!             QL METHOD
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQLRAT.
!
!     This subroutine finds the eigenvalues of a SYMMETRIC
!     TRIDIAGONAL matrix by the rational QL method.
!
!     On Input
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E2 contains the squares of the subdiagonal elements of the
!          symmetric tridiagonal matrix in its last N-1 positions.
!          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
!          dimensioned E2(N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1, 2, ..., IERR-1, but may not be
!          the smallest eigenvalues.
!
!        E2 has been destroyed.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!               C. H. Reinsch, Eigenvalues of a real, symmetric, tri-
!                 diagonal matrix, Algorithm 464, Communications of the
!                 ACM 16, 11 (November 1973), pp. 689.
!***ROUTINES CALLED  PYTHAG, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  TQLRAT
!
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      real(wp) D(*),E2(*)
      real(wp) B,C,F,G,H,P,R,S,MACHEP
      LOGICAL FIRST
!
      SAVE FIRST, MACHEP
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  TQLRAT
      IF (FIRST) THEN
        MACHEP = RMACH(4, 1._wp)
      ENDIF
      FIRST = .FALSE.
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO I = 2, N
        E2(I-1) = E2(I)
      end do
!
      F = 0.0_wp
      B = 0.0_wp
      E2(N) = 0.0_wp
!
      DO 290 L = 1, N
        J = 0
        H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
        IF (B .GT. H) GO TO 105
        B = H
        C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105   DO 110 M = L, N
          IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110   CONTINUE
!
  120   IF (M .EQ. L) GO TO 210
  130   IF (J .EQ. 30) GO TO 1000
        J = J + 1
!     .......... FORM SHIFT ..........
        L1 = L + 1
        S = SQRT(E2(L))
        G = D(L)
        P = (D(L1) - G) / (2.0_wp * S)
        R = PYTHAG(P,1.0_wp)
        D(L) = S / (P + SIGN(R,P))
        H = G - D(L)
!
        DO I = L1, N
          D(I) = D(I) - H
        end do
!
        F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
        G = D(M)
        IF (G .EQ. 0.0_wp) G = B
        H = G
        S = 0.0_wp
        MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
        DO 200 II = 1, MML
          I = M - II
          P = G * H
          R = P + E2(I)
          E2(I+1) = S * R
          S = E2(I) / R
          D(I+1) = H + S * (H + D(I))
          G = D(I) - E2(I) / G
          IF (G .EQ. 0.0_wp) G = B
          H = G * P / R
  200   CONTINUE
!
        E2(L) = S * G
        D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
        IF (H .EQ. 0.0_wp) GO TO 210
        IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
        E2(L) = H * E2(L)
        IF (E2(L) .NE. 0.0_wp) GO TO 130
  210   P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
        IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
        DO 230 II = 2, L
          I = L + 2 - II
          IF (P .GE. D(I-1)) GO TO 270
          D(I) = D(I-1)
  230   CONTINUE
!
  250   I = 1
  270   D(I) = P
  290 END DO
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END subroutine tqlrat
!------------------------------------------------------------------------------
  function spmach(i,r)
  integer  :: i
  real(sp) :: r, spmach, z
  select case (i)
  case (1)
    spmach = tiny(r)
  case (2)
    spmach = huge(r)
  case (3)
    spmach = radix(r)**(-digits(r))
  case (4)
    spmach = epsilon(r)
  case (5)
    z = radix(r)
    spmach = log10(z)
  case default
    call abort ('spmach: i')
  end select
  end function spmach
!------------------------------------------------------------------------------
  function dpmach(i,r)
  integer  :: i
  real(dp) :: r, dpmach, z
  select case (i)
  case (1)
    dpmach = tiny(r)
  case (2)
    dpmach = huge(r)
  case (3)
#if defined (__PGI) || defined (__FLANG)
    dpmach = epsilon(r)/real(radix(r),dp)
#else
    dpmach = radix(r)**(-digits(r))
#endif
  case (4)
    dpmach = epsilon(r)
  case (5)
    z = radix(r)
    dpmach = log10(z)
  case default
    call abort ('dpmach: i')
  end select
  end function dpmach
!------------------------------------------------------------------------------
      REAL(wp) FUNCTION PYTHAG (A, B)
!***BEGIN PROLOGUE  PYTHAG
!***SUBSIDIARY
!***PURPOSE  Compute the complex square root of a complex number without
!            destructive overflow or underflow.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (PYTHAG-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Finds sqrt(A**2+B**2) without overflow or destructive underflow
!
!***SEE ALSO  EISDOC
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PYTHAG
      real(wp) A,B
!
      real(wp) P,Q,R,S,T
!***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0_wp) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0_wp + R
         IF (T .EQ. 4.0_wp) GO TO 20
         S = R/T
         P = P + 2.0_wp*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END function pythag
!------------------------------------------------------------------------------
      SUBROUTINE IPSORT (IX, IPERM, KFLAG, IER)
!***BEGIN PROLOGUE  IPSORT
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1A, N6A2A
!***TYPE      INTEGER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   IPSORT returns the permutation vector IPERM generated by sorting
!   the array IX and, optionally, rearranges the values in IX.  IX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that IX(IPERM(I)) is the Ith value in the
!   rearrangement of IX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between IPSORT and its active sorting equivalent
!   ISORT is that the data are referenced indirectly rather than
!   directly.  Therefore, IPSORT should require approximately twice as
!   long to execute as ISORT.  However, IPSORT is more general.
!
!   Description of Parameters
!      IX - input/output -- integer array of values to be sorted.
!           If ABS(KFLAG) = 2, then the values in IX will be
!           rearranged on output; otherwise, they are unchanged.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              IX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting IX in increasing order and sort IX also.
!            =  1  means return the permutation vector resulting from
!                  sorting IX in increasing order and do not sort IX.
!            = -1  means return the permutation vector resulting from
!                  sorting IX in decreasing order and do not sort IX.
!            = -2  means return the permutation vector resulting from
!                  sorting IX in decreasing order and sort IX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   810801  Further modified by David K. Kahaner.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  IPSORT
!     .. Scalar Arguments ..
      INTEGER IER, KFLAG
!     .. Array Arguments ..
      INTEGER IPERM(:), IX(:)
!     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, INDX, INDX0, ISTRT, ITEMP, J, K, KK, L, LM, LMT, M,&
     &        NN, N
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!!    EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  IPSORT
      n = size (ix)
      if (n==0) return
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
        IER = 1
!!         CALL XERMSG ('SLATEC', 'IPSORT',                               &
!!     &    'The number of values to be sorted, N, is not positive.',     &
!!     &    IER, 1)
        RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        IER = 2
!!         CALL XERMSG ('SLATEC', 'IPSORT',                               &
!!     &    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!!     &    IER, 1)
        RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=1,NN
        IPERM(I) = I
   10 END DO
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. 1) RETURN
!
!     Alter array IX to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
        DO 20 I=1,NN
          IX(I) = -IX(I)
   20   CONTINUE
      ENDIF
!
!     Sort IX only
!
      M = 1
      I = 1
      J = NN
      R = .375E0
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437E0) THEN
        R = R+3.90625E-2
      ELSE
        R = R-0.21875E0
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (IX(IPERM(I)) .GT. IX(LM)) THEN
        IPERM(IJ) = IPERM(I)
        IPERM(I) = LM
        LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (IX(IPERM(J)) .LT. IX(LM)) THEN
        IPERM(IJ) = IPERM(J)
        IPERM(J) = LM
        LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
        IF (IX(IPERM(I)) .GT. IX(LM)) THEN
          IPERM(IJ) = IPERM(I)
          IPERM(I) = LM
          LM = IPERM(IJ)
        ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
      IF (IX(IPERM(L)) .GT. IX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
      IF (IX(IPERM(K)) .LT. IX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M+1
      ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M+1
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
!
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (IX(IPERM(I)) .LE. IX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+1) = IPERM(K)
      K = K-1
!
      IF (IX(LM) .LT. IX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -1) THEN
        DO 130 I=1,NN
          IX(I) = -IX(I)
  130   CONTINUE
      ENDIF
!
!     Rearrange the values of IX if desired
!
      IF (KK .EQ. 2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
        DO 150 ISTRT=1,NN
          IF (IPERM(ISTRT) .GE. 0) THEN
            INDX = ISTRT
            INDX0 = INDX
            ITEMP = IX(ISTRT)
  140       IF (IPERM(INDX) .GT. 0) THEN
              IX(INDX) = IX(IPERM(INDX))
              INDX0 = INDX
              IPERM(INDX) = -IPERM(INDX)
              INDX = ABS(IPERM(INDX))
              GO TO 140
            ENDIF
            IX(INDX0) = ITEMP
          ENDIF
  150   CONTINUE
!
!        Revert the signs of the IPERM values
!
        DO 160 I=1,NN
          IPERM(I) = -IPERM(I)
  160   CONTINUE
!
      ENDIF
!
      RETURN
      END subroutine ipsort
!------------------------------------------------------------------------------
      SUBROUTINE CPSORT (CX, IPERM, KFLAG, IER)
!***BEGIN PROLOGUE  CPSORT
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1A, N6A2A
!***TYPE      INTEGER (SPSORT-S, DPSORT-D, CPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   CPSORT returns the permutation vector IPERM generated by sorting
!   the array CX and, optionally, rearranges the values in CX.  CX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that CX(IPERM(I)) is the Ith value in the
!   rearrangement of CX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between CPSORT and its active sorting equivalent
!   ISORT is that the data are referenced indirectly rather than
!   directly.  Therefore, CPSORT should require approximately twice as
!   long to execute as ISORT.  However, CPSORT is more general.
!
!   Description of Parameters
!      CX - input/output -- integer array of values to be sorted.
!           If ABS(KFLAG) = 2, then the values in CX will be
!           rearranged on output; otherwise, they are unchanged.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              CX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting CX in increasing order and sort CX also.
!            =  1  means return the permutation vector resulting from
!                  sorting CX in increasing order and do not sort CX.
!            = -1  means return the permutation vector resulting from
!                  sorting CX in decreasing order and do not sort CX.
!            = -2  means return the permutation vector resulting from
!                  sorting CX in decreasing order and sort CX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   810801  Further modified by David K. Kahaner.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  CPSORT
!     .. Scalar Arguments ..
      INTEGER IER, KFLAG
!     .. Array Arguments ..
      character(len=*) ::  CX(:)
      integer          ::  IPERM(:)
!     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M,        &
     &        NN, N
      character :: ctemp
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!!    EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  CPSORT
      n = size (cx)
      if (n==0) return
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
        IER = 1
!!         CALL XERMSG ('SLATEC', 'CPSORT',                               &
!!     &    'The number of values to be sorted, N, is not positive.',     &
!!     &    IER, 1)
        RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        IER = 2
!!         CALL XERMSG ('SLATEC', 'CPSORT',                               &
!!     &    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!!     &    IER, 1)
        RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=1,NN
        IPERM(I) = I
   10 END DO
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. 1) RETURN
!
!     Alter array CX to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
        CX = cx( (/(i,i=nn,1,-1)/) )
!        DO 20 I=1,NN
!          CX(I) = -CX(I)
!   20   CONTINUE
      ENDIF
!
!     Sort CX only
!
      M = 1
      I = 1
      J = NN
      R = .375E0
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437E0) THEN
        R = R+3.90625E-2
      ELSE
        R = R-0.21875E0
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (CX(IPERM(I)) .GT. CX(LM)) THEN
        IPERM(IJ) = IPERM(I)
        IPERM(I) = LM
        LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (CX(IPERM(J)) .LT. CX(LM)) THEN
        IPERM(IJ) = IPERM(J)
        IPERM(J) = LM
        LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
        IF (CX(IPERM(I)) .GT. CX(LM)) THEN
          IPERM(IJ) = IPERM(I)
          IPERM(I) = LM
          LM = IPERM(IJ)
        ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
      IF (CX(IPERM(L)) .GT. CX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
      IF (CX(IPERM(K)) .LT. CX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M+1
      ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M+1
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
!
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (CX(IPERM(I)) .LE. CX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+1) = IPERM(K)
      K = K-1
!
      IF (CX(LM) .LT. CX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -1) THEN
        CX = cx( (/(i,i=nn,1,-1)/) )
!        DO 130 I=1,NN
!          CX(I) = -CX(I)
!  130   CONTINUE
      ENDIF
!
!     Rearrange the values of CX if desired
!
      IF (KK .EQ. 2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
        DO 150 ISTRT=1,NN
          IF (IPERM(ISTRT) .GE. 0) THEN
            INDX = ISTRT
            INDX0 = INDX
            CTEMP = CX(ISTRT)
  140       IF (IPERM(INDX) .GT. 0) THEN
              CX(INDX) = CX(IPERM(INDX))
              INDX0 = INDX
              IPERM(INDX) = -IPERM(INDX)
              INDX = ABS(IPERM(INDX))
              GO TO 140
            ENDIF
            CX(INDX0) = CTEMP
          ENDIF
  150   CONTINUE
!
!        Revert the signs of the IPERM values
!
        DO 160 I=1,NN
          IPERM(I) = -IPERM(I)
  160   CONTINUE
!
      ENDIF
!
      RETURN
      END subroutine cpsort
!------------------------------------------------------------------------------
      SUBROUTINE DPSORT (DX, IPERM, KFLAG, IER)
!***BEGIN PROLOGUE  DPSORT
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1B, N6A2B
!***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   DPSORT returns the permutation vector IPERM generated by sorting
!   the array DX and, optionally, rearranges the values in DX.  DX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that DX(IPERM(I)) is the Ith value in the
!   rearrangement of DX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between DPSORT and its active sorting equivalent
!   DSORT is that the data are referenced indirectly rather than
!   directly.  Therefore, DPSORT should require approximately twice as
!   long to execute as DSORT.  However, DPSORT is more general.
!
!   Description of Parameters
!      DX - input/output -- double precision array of values to be
!           sorted.  If ABS(KFLAG) = 2, then the values in DX will be
!           rearranged on output; otherwise, they are unchanged.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              DX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting DX in increasing order and sort DX also.
!            =  1  means return the permutation vector resulting from
!                  sorting DX in increasing order and do not sort DX.
!            = -1  means return the permutation vector resulting from
!                  sorting DX in decreasing order and do not sort DX.
!            = -2  means return the permutation vector resulting from
!                  sorting DX in decreasing order and sort DX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890619  Double precision version of SPSORT created by D. W. Lozier.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  DPSORT
!     .. Scalar Arguments ..
      INTEGER IER, KFLAG
!     .. Array Arguments ..
      real(dp) DX(:)
      INTEGER IPERM(:)
!     .. Local Scalars ..
      real(dp) R, TEMP
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN, N
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!!    EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  DPSORT
      n = size (dx)
      if (n==0) return
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
        IER = 1
!!         CALL XERMSG ('SLATEC', 'DPSORT',                               &
!!     &    'The number of values to be sorted, N, is not positive.',     &
!!     &    IER, 1)
        RETURN
      ENDIF
!
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        IER = 2
!!         CALL XERMSG ('SLATEC', 'DPSORT',                               &
!!     &    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!!     &    IER, 1)
        RETURN
      ENDIF
!
!     Initialize permutation vector
!
      DO 10 I=1,NN
        IPERM(I) = I
   10 END DO
!
!     Return if only one value is to be sorted
!
      IF (NN .EQ. 1) RETURN
!
!     Alter array DX to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
        DO 20 I=1,NN
          DX(I) = -DX(I)
   20   CONTINUE
      ENDIF
!
!     Sort DX only
!
      M = 1
      I = 1
      J = NN
      R = .375_dp
!
   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437_dp) THEN
        R = R+3.90625E-2_dp
      ELSE
        R = R-0.21875_dp
      ENDIF
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
      IF (DX(IPERM(I)) .GT. DX(LM)) THEN
        IPERM(IJ) = IPERM(I)
        IPERM(I) = LM
        LM = IPERM(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than LM, interchange with LM
!
      IF (DX(IPERM(J)) .LT. DX(LM)) THEN
        IPERM(IJ) = IPERM(J)
        IPERM(J) = LM
        LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
        IF (DX(IPERM(I)) .GT. DX(LM)) THEN
          IPERM(IJ) = IPERM(I)
          IPERM(I) = LM
          LM = IPERM(IJ)
        ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
      IF (DX(IPERM(L)) .GT. DX(LM)) GO TO 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
      IF (DX(IPERM(K)) .LT. DX(LM)) GO TO 70
!
!     Interchange these elements
!
      IF (K .LE. L) GO TO 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M+1
      ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M+1
      ENDIF
      GO TO 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)
!
   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1
!
  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (DX(IPERM(I)) .LE. DX(LM)) GO TO 100
      K = I
!
  110 IPERM(K+1) = IPERM(K)
      K = K-1
      IF (DX(LM) .LT. DX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
!
!     Clean up
!
  120 IF (KFLAG .LE. -1) THEN
        DO 130 I=1,NN
          DX(I) = -DX(I)
  130   CONTINUE
      ENDIF
!
!     Rearrange the values of DX if desired
!
      IF (KK .EQ. 2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
        DO 150 ISTRT=1,NN
          IF (IPERM(ISTRT) .GE. 0) THEN
            INDX = ISTRT
            INDX0 = INDX
            TEMP = DX(ISTRT)
  140       IF (IPERM(INDX) .GT. 0) THEN
              DX(INDX) = DX(IPERM(INDX))
              INDX0 = INDX
              IPERM(INDX) = -IPERM(INDX)
              INDX = ABS(IPERM(INDX))
              GO TO 140
            ENDIF
            DX(INDX0) = TEMP
          ENDIF
  150   CONTINUE
!
!        Revert the signs of the IPERM values
!
        DO 160 I=1,NN
          IPERM(I) = -IPERM(I)
  160   CONTINUE
!
      ENDIF
!
      RETURN
      END subroutine dpsort
!==============================================================================
      SUBROUTINE SPSORT (DX, IPERM, KFLAG, IER)
!------------------------------------------------------------------------------
!***BEGIN PROLOGUE  SPSORT
!***PURPOSE  Return the permutation vector generated by sorting a given
!            array and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A1B, N6A2B
!***TYPE      DOUBLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SPSORT returns the permutation vector IPERM generated by sorting
!   the array DX and, optionally, rearranges the values in DX.  DX may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that DX(IPERM(I)) is the Ith value in the
!   rearrangement of DX.  IPERM may be applied to another array by
!   calling IPPERM, SPPERM, DPPERM or HPPERM.
!
!   The main difference between SPSORT and its active sorting equivalent
!   DSORT is that the data are referenced indirectly rather than
!   directly.  Therefore, SPSORT should require approximately twice as
!   long to execute as DSORT.  However, SPSORT is more general.
!
!   Description of Parameters
!      DX - input/output -- double precision array of values to be
!           sorted.  If ABS(KFLAG) = 2, then the values in DX will be
!           rearranged on output; otherwise, they are unchanged.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              DX array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting DX in increasing order and sort DX also.
!            =  1  means return the permutation vector resulting from
!                  sorting DX in increasing order and do not sort DX.
!            = -1  means return the permutation vector resulting from
!                  sorting DX in decreasing order and do not sort DX.
!            = -2  means return the permutation vector resulting from
!                  sorting DX in decreasing order and sort DX also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890619  Double precision version of SPSORT created by D. W. Lozier.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG.LT.0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  SPSORT
!------------------------------------------------------------------------------
      !-----------------
      ! Scalar Arguments
      !-----------------
      INTEGER IER, KFLAG
      !----------------
      ! Array Arguments
      !----------------
      real(sp) DX(:)
      INTEGER IPERM(:)
      !--------------
      ! Local Scalars
      !--------------
      real(sp) R, TEMP
      INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN, N
      !-------------
      ! Local Arrays
      !-------------
      INTEGER IL(21), IU(21)
      !----------------------------------
      ! External Subroutine Calls removed
      ! EXTERNAL XERMSG
      !----------------------------------
      !--------------------
      ! Intrinsic Functions
      !--------------------
      INTRINSIC ABS, INT
!-------------------------------------
!***FIRST EXECUTABLE STATEMENT  SPSORT
!-------------------------------------
      n = size (dx)
      if (n==0) return
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
        IER = 1
!!         CALL XERMSG ('SLATEC', 'SPSORT',                               &
!!     &    'The number of values to be sorted, N, is not positive.',     &
!!     &    IER, 1)
        RETURN
      ENDIF

      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        IER = 2
!!         CALL XERMSG ('SLATEC', 'SPSORT',                               &
!!     &    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
!!     &    IER, 1)
        RETURN
      ENDIF
      !------------------------------
      ! Initialize permutation vector
      !------------------------------
      DO I=1,NN
        IPERM(I) = I
      END DO
      !-----------------------------------------
      ! Return if only one value is to be sorted
      !-----------------------------------------
      IF (NN .EQ. 1) RETURN
      !-------------------------------------------------
      ! Alter array DX to get decreasing order if needed
      !-------------------------------------------------
      IF (KFLAG .LE. -1) THEN
        DO I=1,NN
          DX(I) = -DX(I)
        end do
      ENDIF
      !-------------
      ! Sort DX only
      !-------------
      M = 1
      I = 1
      J = NN
      R = .375_sp

   30 IF (I .EQ. J) GO TO 80
      IF (R .LE. 0.5898437_sp) THEN
        R = R+3.90625E-2_sp
      ELSE
        R = R-0.21875_sp
      ENDIF

   40 K = I
      !----------------------------------------------------------------
      ! Select a central element of the array and save it in location L
      !----------------------------------------------------------------
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
      !------------------------------------------------------------------
      ! If first element of array is greater than LM, interchange with LM
      !------------------------------------------------------------------
      IF (DX(IPERM(I)) .GT. DX(LM)) THEN
        IPERM(IJ) = IPERM(I)
        IPERM(I) = LM
        LM = IPERM(IJ)
      ENDIF
      L = J
      !--------------------------------------------------------------
      ! If last element of array is less than LM, interchange with LM
      !--------------------------------------------------------------
      IF (DX(IPERM(J)) .LT. DX(LM)) THEN
        IPERM(IJ) = IPERM(J)
        IPERM(J) = LM
        LM = IPERM(IJ)
        !----------------------------------------------------------
        ! If first element of array is greater than LM, interchange
        ! with LM
        !----------------------------------------------------------
        IF (DX(IPERM(I)) .GT. DX(LM)) THEN
          IPERM(IJ) = IPERM(I)
          IPERM(I) = LM
          LM = IPERM(IJ)
        ENDIF
      ENDIF
      GO TO 60
   50 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
      !-----------------------------------------------------------------
      ! Find an element in the second half of the array which is smaller
      ! than LM
      !-----------------------------------------------------------------
   60 L = L-1
      IF (DX(IPERM(L)) .GT. DX(LM)) GO TO 60
      !----------------------------------------------------------------
      ! Find an element in the first half of the array which is greater
      ! than LM
      !----------------------------------------------------------------
   70 K = K+1
      IF (DX(IPERM(K)) .LT. DX(LM)) GO TO 70
      !---------------------------
      ! Interchange these elements
      !---------------------------
      IF (K .LE. L) GO TO 50
      !--------------------------------------------------------------
      ! Save upper and lower subscripts of the array yet to be sorted
      !--------------------------------------------------------------
      IF (L-I .GT. J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M+1
      ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M+1
      ENDIF
      GO TO 90
      !-----------------------------------------------------
      ! Begin again on another portion of the unsorted array
      !-----------------------------------------------------
   80 M = M-1
      IF (M .EQ. 0) GO TO 120
      I = IL(M)
      J = IU(M)

   90 IF (J-I .GE. 1) GO TO 40
      IF (I .EQ. 1) GO TO 30
      I = I-1

  100 I = I+1
      IF (I .EQ. J) GO TO 80
      LM = IPERM(I+1)
      IF (DX(IPERM(I)) .LE. DX(LM)) GO TO 100
      K = I

  110 IPERM(K+1) = IPERM(K)
      K = K-1
      IF (DX(LM) .LT. DX(IPERM(K))) GO TO 110
      IPERM(K+1) = LM
      GO TO 100
      !---------
      ! Clean up
      !---------
  120 IF (KFLAG .LE. -1) THEN
        DO I=1,NN
          DX(I) = -DX(I)
        end do
      ENDIF
      !--------------------------------------
      ! Rearrange the values of DX if desired
      !--------------------------------------
      IF (KK .EQ. 2) THEN
        !------------------------------------------------------------
        ! Use the IPERM vector as a flag.
        ! If IPERM(I) < 0, then the I-th value is in correct location
        !------------------------------------------------------------
        DO 150 ISTRT=1,NN
          IF (IPERM(ISTRT) .GE. 0) THEN
            INDX = ISTRT
            INDX0 = INDX
            TEMP = DX(ISTRT)
  140       IF (IPERM(INDX) .GT. 0) THEN
              DX(INDX) = DX(IPERM(INDX))
              INDX0 = INDX
              IPERM(INDX) = -IPERM(INDX)
              INDX = ABS(IPERM(INDX))
              GO TO 140
            ENDIF
            DX(INDX0) = TEMP
          ENDIF
  150   CONTINUE
        !-------------------------------------
        ! Revert the signs of the IPERM values
        !-------------------------------------
        DO I=1,NN
          IPERM(I) = -IPERM(I)
        end do

      ENDIF

      END subroutine spsort

!==============================================================================
  !----------------------------
  ! Emulate DPSORT using SORTRX
  !----------------------------
  subroutine dpsortx (DX, IPERM, KFLAG, IER)
    !-----------------
    ! Scalar Arguments
    !-----------------
    integer  ,intent(in)    :: KFLAG
    integer  ,intent(out)   :: IER
    !----------------
    ! Array Arguments
    !----------------
    real(dp)                :: DX(:)
    integer  ,intent(out)   :: IPERM(:)
    !--------------
    ! Local Scalars
    !--------------
    integer  :: n

    ier = 0
    n = size (dx)
    select case (kflag)
    case ( 1, 2)
       call SORTRX (n,  DX, IPERM)
       if (kflag ==  2) DX = DX(IPERM(1:n))
    case (-1,-2)
       call SORTRX (n, -DX, IPERM)
       if (kflag == -2) DX = DX(IPERM(1:n))
    case default
       ier = 2
    end select
  end subroutine dpsortx
!------------------------------------------------------------------------------
  !-------------------------------
  ! Emulate SPSORT using SORT_real
  !-------------------------------
  subroutine spsortx (DX, IPERM, KFLAG, IER)
    !-----------------
    ! Scalar Arguments
    !-----------------
    integer  ,intent(in)    :: KFLAG
    integer  ,intent(out)   :: IER
    !----------------
    ! Array Arguments
    !----------------
    real(sp)                :: DX(:)
    integer  ,intent(out)   :: IPERM(:)
    !--------------
    ! Local Scalars
    !--------------
    integer  :: n

    ier = 0
    n = size (dx)
    select case (kflag)
    case ( 1, 2)
       call SORT_real (n,  DX, IPERM)
       if (kflag ==  2) DX = DX(IPERM(1:n))
    case (-1,-2)
       call SORT_real (n, -DX, IPERM)
       if (kflag == -2) DX = DX(IPERM(1:n))
    case default
       ier = 2
    end select
  end subroutine spsortx
!------------------------------------------------------------------------------
  !----------------------------
  ! Emulate IPSORT using SORTIX
  !----------------------------
  subroutine ipsortx (DX, IPERM, KFLAG, IER)
    !-----------------
    ! Scalar Arguments
    !-----------------
    integer  ,intent(in)    :: KFLAG
    integer  ,intent(out)   :: IER
    !----------------
    ! Array Arguments
    !----------------
    integer                 :: DX(:)
    integer  ,intent(out)   :: IPERM(:)
    !--------------
    ! Local Scalars
    !--------------
    integer  :: n

    ier = 0
    n = size (dx)
    select case (kflag)
    case ( 1, 2)
       call SORTIX (n,  DX, IPERM)
       if (kflag ==  2) DX = DX(IPERM(1:n))
    case (-1,-2)
       call SORTIX (n, -DX, IPERM)
       if (kflag == -2) DX = DX(IPERM(1:n))
    case default
       ier = 2
    end select
  end subroutine ipsortx
!------------------------------------------------------------------------------
  subroutine i8sortx (DX, IPERM, KFLAG, IER)
    !-----------------
    ! Scalar Arguments
    !-----------------
    integer  ,intent(in)    :: KFLAG
    integer  ,intent(out)   :: IER
    !----------------
    ! Array Arguments
    !----------------
    integer(i8)             :: DX(:)
    integer  ,intent(out)   :: IPERM(:)
    !--------------
    ! Local Scalars
    !--------------
    integer  :: n

    ier = 0
    n = size (dx)
    select case (kflag)
    case ( 1, 2)
       call SORT_int8 (n,  DX, IPERM)
       if (kflag ==  2) DX = DX(IPERM(1:n))
    case (-1,-2)
       call SORT_int8 (n, -DX, IPERM)
       if (kflag == -2) DX = DX(IPERM(1:n))
    case default
       ier = 2
    end select
  end subroutine i8sortx
!------------------------------------------------------------------------------
  !-------------------------------
  ! Emulate CPSORT using SORT_char
  !-------------------------------
  subroutine cpsortx (DX, IPERM, KFLAG, IER)
    !-----------------
    ! Scalar Arguments
    !-----------------
    integer  ,intent(in)    :: KFLAG
    integer  ,intent(out)   :: IER
    !----------------
    ! Array Arguments
    !----------------
    character(len=*)        :: DX(:)
    integer  ,intent(out)   :: IPERM(:)
    !--------------
    ! Local Scalars
    !--------------
    integer  :: n

    ier = 0
    n = size (dx)
    call SORT_char (n, DX, IPERM)
    select case (kflag)
    case ( 1, 2)
       if (kflag ==  2) DX = DX(IPERM(1:n))
    case (-1,-2)
       ! Beware: reverse sorting of strings may provide unstable results!
       IPERM = IPERM(n:1:-1)
       if (kflag == -2) DX = DX(IPERM(1:n))
    case default
       ier = 2
    end select
  end subroutine cpsortx
!==============================================================================
! From Leonard J. Moss of SLAC:

! Here's a hybrid QuickSort I wrote a number of years ago.  It's
! based on suggestions in Knuth, Volume 3, and performs much better
! than a pure QuickSort on short or partially ordered input arrays.

! SUBROUTINE SORTRX (N, DATA, INDEX)
!===================================================================
!
!     SORTRX -- SORT, Real input, indeX output
!
!
!     Input:  N     INTEGER
!             DATA  REAL
!
!     Output: INDEX INTEGER (DIMENSION N)
!
! This routine performs an in-memory sort of the first N elements of
! array DATA, returning into array INDEX the indices of elements of
! DATA arranged in ascending order.  Thus,
!
!    DATA(INDEX(1)) will be the smallest number in array DATA;
!    DATA(INDEX(N)) will be the largest number in DATA.
!
! The original data is not physically rearranged.  The original order
! of equal input values is not necessarily preserved.
!
!===================================================================
!
! SORTRX uses a hybrid QuickSort algorithm, based on several
! suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
! "pivot key" [my term] for dividing each subsequence is chosen to be
! the median of the first, last, and middle values of the subsequence;
! and the QuickSort is cut off when a subsequence has 9 or fewer
! elements, and a straight insertion sort of the entire array is done
! at the end.  The result is comparable to a pure insertion sort for
! very short arrays, and very fast for very large arrays (of order 12
! micro-sec/element on the 3081K for arrays of 10K elements).  It is
! also not subject to the poor performance of the pure QuickSort on
! partially ordered data.
!
! Created:  15 Jul 1986  Len Moss
!
!-------------------------------------------------------------------
!
! The present version is a minor rewrite in Fortran 95 and takes into
! account the initial order through a modified final insertion sort.
! It performs sufficiently well on NEC SX-9.
!
! Harald Anlauf  (2009).
!
!===================================================================
! The actual implementation was moved to the template ljmqsort.incf.
!==============================================================================
! Define qsort for double precision reals
#define SORT_MYTYPE  SORTRX
#define MYTYPE       REAL(dp)
#undef  TYPE_IS_CHAR_ARRAY
#include "ljmqsort.incf"
#undef  SORT_MYTYPE
#undef  MYTYPE
!==============================================================================
! Define qsort for default integer
#define SORT_MYTYPE  SORTIX
#define MYTYPE       INTEGER
#undef  TYPE_IS_CHAR_ARRAY
#include "ljmqsort.incf"
#undef  SORT_MYTYPE
#undef  MYTYPE
!==============================================================================
! Define qsort for integer(8)
#define SORT_MYTYPE  sort_int8
#define MYTYPE       INTEGER(i8)
#undef  TYPE_IS_CHAR_ARRAY
#include "ljmqsort.incf"
#undef  SORT_MYTYPE
#undef  MYTYPE
!==============================================================================
! Define qsort for default reals
#define SORT_MYTYPE  sort_real
#define MYTYPE       REAL(sp)
#undef  TYPE_IS_CHAR_ARRAY
#include "ljmqsort.incf"
#undef  SORT_MYTYPE
#undef  MYTYPE
!==============================================================================
! Define qsort for character(*) arrays
#define SORT_MYTYPE  sort_char
#define MYTYPE       CHARACTER(len=*)
#define TYPE_IS_CHAR_ARRAY
#include "ljmqsort.incf"
#undef  SORT_MYTYPE
#undef  MYTYPE
#undef  TYPE_IS_CHAR_ARRAY
!==============================================================================
  SUBROUTINE xRS_opt (NM, Z, W, IERR)
! Compute the eigenvalues and the eigenvectors of multiple
! real symmetric matrices in one step.
! "Vectorized" version of the SLATEC (EISPACK) routine xRS by loop-pushing.
!
!     On Input
!
!        NM must be set to the total number of matrices to be diagonlized
!
!        Z contains the collection of real symmetric matrices.
!          Z is a three-dimensional REAL array, dimensioned Z(NM,N,N).
!
!     On Output
!
!        Z contains the collection of eigenvectors.  The
!          eigenvectors are orthonormal.  Z is a three-dimensional
!          REAL array, dimensioned Z(NM,N,N).
!
!        W contains the collection of eigenvalues in ascending order.
!          W is a two-dimensional REAL array, dimensioned W(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues, and eigenvectors if requested,
!                     should be correct for indices 1, 2, ..., IERR-1.
!
    INTEGER,  intent(in)     :: NM
    INTEGER,  intent(out)    :: IERR
    REAL(wp), intent(inout)  :: Z(:,:,:)
    REAL(wp), intent(out)    :: W(:,:)

    ! Local variables
    INTEGER :: N, n1, n2
    REAL(wp), allocatable :: E(:,:)

    n1 = size (Z, dim=1)
    n2 = size (Z, dim=2)
    n  = size (Z, dim=3)     ! Dimension of eigenproblems

 !    Check valid value for number of eigenproblems (NM)
    if (NM > n1) then
       IERR = 10 * NM
       return
    end IF

    allocate (E(n1,n))
    call tred2_opt (NM, Z(:,:,:), W(:,:), E(:,:))
    call tql2_opt(NM,W(:,:),E(:,:),Z(:,:,:),IERR)
    if (ierr /= 0) then
       write (0,*) "*** FATAL ***", ierr
       stop
    end if
  END subroutine xRS_opt

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine tred2_opt (NM, Z, D, E)
!
!     This subroutine reduces a collection of REAL SYMMETRIC matrices to
!     symmetric tridiagonal form using and accumulating
!     orthogonal similarity transformations.  This routine is intended to
!     be called only by xRS_opt.
!
!    On Input
!
!        NM must be set to the total number of matrices to be reduced.
!        NM is an INTEGER variable.
!
!        Z contains the real symmetric input matrices.  Only the lower
!          triangle of the matrix need be supplied.  Z is a three-
!          dimensional REAL array, dimensioned Z(NM,N,N).
!
!     On Output
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrices.  D is a two-dimensional REAL array, dimensioned D(NM,N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrices in their last N-1 positions.  E(1) is set
!          to zero.  E is a two-dimensional REAL array, dimensioned
!          E(NM,N).
!
!        Z contains the orthogonal transformation matrices produced in
!          the reduction.  Z is a three-dimensional REAL array,
!          dimensioned Z(NM,N,N).
!

    INTEGER,  intent(in)     :: NM
    real(wp), intent(inout)  :: Z(:,:,:)
    real(wp), intent(out)    :: D(:,:), E(:,:)

    INTEGER :: I,J,K,L,N, nn, ii, kk
    real(wp), dimension(NM) :: F, G, H, HH, SCALE

    logical,  dimension(NM) :: mask
    integer,  allocatable   :: ix(:)
!
    n = size (Z, dim=3)

!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
    DO 300 I = N, 2, -1
       L = I - 1
       H(:NM) = 0.0_wp
       SCALE(:NM) = 0.0_wp
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
       DO K = 1, L
          SCALE(:NM) = SCALE(:NM) + ABS(Z(:NM,I,K))
       end do
!
       mask(1:NM) = (SCALE(:NM) /= 0.0_wp)
       nn = count (mask)
       allocate (ix(nn))
       ix(1:nn) = pack ((/ (k, k=1,NM) /), mask)

       where (.not. mask)
          E(:NM,I) = Z(:NM,I,L)
       end where
!
       DO K = 1, L
!NEC$ ivdep
          do kk = 1, nn
             ii = ix(kk)
             Z(ii,I,K) = Z(ii,I,K) / SCALE(ii)
             H(ii)     = H(ii) + Z(ii,I,K) ** 2
          end do
       END DO
!
!NEC$ ivdep
       do kk = 1, nn
          ii = ix(kk)
          F(ii) = Z(ii,I,L)
          G(ii) = -SIGN(SQRT(H(ii)),F(ii))
          E(ii,I) = SCALE(ii) * G(ii)
          H(ii) = H(ii) - F(ii) * G(ii)
          Z(ii,I,L) = F(ii) - G(ii)
          F(ii) = 0.0_wp
       end do
!
       DO J = 1, L
!NEC$ ivdep
          do kk = 1, nn
             ii = ix(kk)
             Z(ii,J,I) = Z(ii,I,J) / H(ii)
             G(ii) = 0.0_wp
          end do
!     .......... FORM ELEMENT OF A*U ..........
          DO K = 1, J
!NEC$ ivdep
             do kk = 1, nn
                ii = ix(kk)
                G(ii) = G(ii) + Z(ii,J,K) * Z(ii,I,K)
             end do
          end do
!
          DO K = J+1, L
!NEC$ ivdep
             do kk = 1, nn
                ii = ix(kk)
                G(ii) = G(ii) + Z(ii,K,J) * Z(ii,I,K)
             end do
          end do
!     .......... FORM ELEMENT OF P ..........
!NEC$ ivdep
          do kk = 1, nn
             ii = ix(kk)
             E(ii,J) = G(ii) / H(ii)
             F(ii) = F(ii) + E(ii,J) * Z(ii,I,J)
          end do
       END DO
!
       HH(ix(:)) = F(ix(:)) / (2 * H(ix(:)))
!     .......... FORM REDUCED A ..........
       DO J = 1, L
!NEC$ ivdep
          do kk = 1, nn
             ii = ix(kk)
             F(ii) = Z(ii,I,J)
             G(ii) = E(ii,J) - HH(ii) * F(ii)
             E(ii,J) = G(ii)
          end do
!
          DO K = 1, J
!NEC$ ivdep
             do kk = 1, nn
                ii = ix(kk)
                Z(ii,J,K) = Z(ii,J,K) - F(ii) * E(ii,K) - G(ii) * Z(ii,I,K)
             end do
          end do
       end do
!
       D(:NM,I) = H(:NM)

       deallocate (ix)
300 END DO
!
    D(NM,1) = 0.0_wp
    E(NM,1) = 0.0_wp
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
    DO I = 1, N
       L = I - 1
!!! CHECK !!!
!!!        IF (D(I) .EQ. 0.0_wp) GO TO 380
!
       DO J = 1, L
          G(:NM) = 0.0_wp
!
          DO K = 1, L
             G(:NM) = G(:NM) + Z(:NM,I,K) * Z(:NM,K,J)
          end do
!
          DO K = 1, L
             Z(:NM,K,J) = Z(:NM,K,J) - G(:NM) * Z(:NM,K,I)
          end do
       end do
!
380    D(:NM,I) = Z(:NM,I,I)
       Z(:NM,I,I) = 1.0_wp
!
       DO J = 1, L
          Z(:NM,I,J) = 0.0_wp
          Z(:NM,J,I) = 0.0_wp
       END DO
!
    END DO
!
  END subroutine tred2_opt
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE tql2_opt (NM, D, E, Z, IERR)
!     This subroutine is a translation of the ALGOL procedure TQL2,
!     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!     Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a collection of SYMMETRIC TRIDIAGONAL matrices by the
!     QL method.
!     The eigenvectors of FULL SYMMETRIC matrices can also
!     be found if  TRED2_OPT  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     On Input
!
!        NM must be set to the total number of tridiagonal
!        matrices to be processed. NM is an INTEGER variable.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrices.  D is a two-dimensional REAL array, dimensioned D(NM,N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrices in their last N-1 positions.  E(1) is
!          arbitrary.  E is a two-dimensional REAL array, dimensioned
!          E(NM,N).
!
!        Z contains the transformation matrices produced in the
!          reduction by  TRED2, if performed.  If the eigenvectors
!          of the tridiagonal matrices are desired, Z must contain
!          the identity matrix.  Z is a three-dimensional REAL array,
!          dimensioned Z(NM,N,N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., IERR-1.
!
!        E has been destroyed.
!
!        Z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrices. If an error exit is made,
!          Z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!
!


      integer,  intent(in)     :: NM
      real(wp), intent(inout)  :: D(:,:), E(:,:)
      real(wp), intent(inout)  :: Z(:,:,:)
      integer  :: I,J,K,L,M,N,II,L1,L2,IERR,N1,N2,KK,MMAX
      integer,  dimension(NM) :: MM
      real(wp), dimension(NM) :: B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      logical,  dimension(NM) :: mask1, mask2
      integer, allocatable    :: idx1(:)
      integer                 :: idx2(NM)

      n = size(Z,dim=3)
!
!***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) return
!
      DO I = 2, N
        E(:NM,I-1) = E(:NM,I)
      end do
!
      F(:NM) = 0.0_wp
      B(:NM) = 0.0_wp
      E(:NM,N) = 0.0_wp
!
      DO L = 1, N
        H(:NM) = ABS(D(:NM,L)) + ABS(E(:NM,L))
        where (B(:NM) .LT. H(:NM))
           B(:NM) = H(:NM)
        end where
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
        DO K=1, NM
           DO M = L, N
              IF (B(K) + ABS(E(K,M)) .EQ. B(K)) EXIT
              !     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
              !                THROUGH THE BOTTOM OF THE LOOP ..........
           END DO
           MM(K) = M
        END DO
        J = 0
        ! start j-loop
        DO
           ! update of mask1,idx1
           mask1=(MM /= L)
           n1=count(mask1)
           if(n1==0) exit
           IF (J .EQ. 30) EXIT
           J = J + 1
!#ifdef _FTRACE
!           CALL ftrace_region_begin('idx1 and form shift')
!#endif
           allocate(idx1(1:n1))
           idx1(1:n1) = pack((/ (k,k=1,NM) /),mask1)
           !     .......... FORM SHIFT ..........
           L1 = L + 1
           L2 = L1 + 1
!NEC$ ivdep
           do kk = 1, n1
              ii = idx1(kk)
              G(ii) =  D(ii,L)
              P(ii) = (D(ii,L1) - G(ii)) / (2.0_wp * E(ii,L))
              R(ii) = sqrt(P(ii)**2 + 1.0_wp)
              D(ii,L)  = E(ii,L) / (P(ii) + SIGN(R(ii),P(ii)))
              D(ii,L1) = E(ii,L) * (P(ii) + SIGN(R(ii),P(ii)))
              DL1(ii) = D(ii,L1)
              H(ii) = G(ii) - D(ii,L)
              !
           end do ! kk
!NEC$ ivdep
           do kk = 1, n1
              ii = idx1(kk)
              DO I = L2, N
                 D(ii,I) = D(ii,I) - H(ii)
              end do
           end do ! kk
           !

!#ifdef _FTRACE
!           CALL ftrace_region_end('idx1 and form shift')
!#endif
!#ifdef _FTRACE
!           CALL ftrace_region_begin('QL Trafo')
!#endif
!NEC$ ivdep
           do kk=1,n1
              k=idx1(kk)
              F(k) = F(k) + H(k)
           !     .......... QL TRANSFORMATION ..........
              P(k) = D(k,MM(k))
              C(k) = 1.0_wp
              C2(k) = C(k)
              EL1(k) = E(k,L1)
              S(k) = 0.0_wp
           end do
           !     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
           MMAX = maxval(MM)
           DO I = MMAX-1,L,-1
              mask2 = (I .le. MM-1)
              n2 = count(mask2)
!             allocate(idx2(1:n2))
              idx2(:n2) = pack((/ (k,k=1,NM) /),mask2)
!NEC$ ivdep
              DO KK=1,n2
                 k = idx2(KK)
                 C3(k) = C2(k)
                 C2(k) = C(k)
                 S2(k) = S(k)
                 G(k) = C(k) * E(k,I)
                 H(k) = C(k) * P(k)
                 IF (ABS(P(k)) .LT. ABS(E(k,I))) THEN
                    C(k) = P(k) / E(k,I)
                    R(k) = SQRT(C(k)*C(k)+1.0_wp)
                    E(k,I+1) = S(k) * E(k,I) * R(k)
                    S(k) = 1.0_wp / R(k)
                    C(k) = C(k) * S(k)
                 ELSE
                    C(k) = E(k,I) / P(k)
                    R(k) = SQRT(C(k)*C(k)+1.0_wp)
                    E(k,I+1) = S(k) * P(k) * R(k)
                    S(k) = C(k) / R(k)
                    C(k) = 1.0_wp / R(k)
                 END IF
              P(k) = C(k) * D(k,I) - S(k) * G(k)
              D(k,I+1) = H(k) + S(k) * (C(k) * G(k) + S(k) * D(k,I))
              END DO
              !     .......... FORM VECTOR ..........
              DO K = 1, N
!NEC$ ivdep
                 do kk = 1, n2
                    ii = idx2(kk)
                    H(ii) = Z(ii,K,I+1)
                    Z(ii,K,I+1) = S(ii) * Z(ii,K,I) + C(ii) * H(ii)
                    Z(ii,K,I)   = C(ii) * Z(ii,K,I) - S(ii) * H(ii)
                 end do
              END DO
              !
!              deallocate(idx2)
           END DO
           !
!NEC$ ivdep
           do kk = 1, n1
              ii = idx1(kk)
              P(ii)   = -S(ii) * S2(ii) * C3(ii) * EL1(ii) * E(ii,L) / DL1(ii)
              E(ii,L) =  S(ii) * P(ii)
              D(ii,L) =  C(ii) * P(ii)
           end do
!#ifdef _FTRACE
!           CALL ftrace_region_end('QL Trafo')
!#endif
           where (B(:NM) + ABS(E(:NM,L)) .EQ. B(:NM))
              MM(:NM) = L
           end where
           deallocate(idx1)
           !end of j-loop
        END DO
        IF(J .EQ. 30) THEN
           IERR = L
           RETURN
        END IF
        D(:NM,L) = D(:NM,L) + F(:NM)
     END DO

!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
!#ifdef _FTRACE
!     CALL ftrace_region_begin('order ew and ev')
!#endif
     DO kk=1,NM
      DO 300 II = 2, N
        I = II - 1
        K = I
        P(kk) = D(kk,I)
!
        DO 260 J = II, N
          IF (D(kk,J) .GE. P(kk)) GO TO 260
          K = J
          P(kk) = D(kk,J)
260     CONTINUE
!
        IF (K .EQ. I) GO TO 300
        D(kk,K) = D(kk,I)
        D(kk,I) = P(kk)
!
        DO J = 1, N
          P(kk) = Z(kk,J,I)
          Z(kk,J,I) = Z(kk,J,K)
          Z(kk,J,K) = P(kk)
        END DO
!
300  END DO
     END DO
!#ifdef _FTRACE
!     CALL ftrace_region_end('order ew and ev')
!#endif
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END subroutine tql2_opt

!------------------------------------------------------------------------------
end module slatec_module
