!
!+ GNSS Radio occultation operator: Matrix inversion,quasi-inversion,regression
!
MODULE Matrix
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Matrix inversion, quasi-inversion, regression.
!
! Current Code Owner: DWD, Andreas Rhodin
!    phone: +49 69 8062 2722
!    fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Reference:
!   Michael E. Gorbunov and Luis Kornblueh
!   Principles of variational assimilation of GNSS radio occultation data.
!   Max-Planck-Institut fuer Meteorologie, Hamburg, Report No. 350 (2003)
!
! Author:
! Michael E. Gorbunov  2004  original code
! Changes:
! Andreas Rhodin             adapted to DWD 3D-VAR
!==============================================================================
!
! Module Matrix
!
! Matrix inversion, quasi-inversion, regression.
!----------------------------------------------------------
! (C) Copyright 1998-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 29 Sep 1998 | Original code.
!   2.0   | 05 Nov 1998 | Diag.
!   3.0   | 25 Feb 1999 | Fast_Invert_Matrix.
!   4.0   | 24 Apr 2000 | Tensor_Product.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
Private
Public :: Regression, Basic_Polynomials, Quasi_Invert, Invert_Matrix, &
          Basic_Splines, Diag, Tensor_Product
!----------------------------------------------------------
!
Contains


!==========================================================
Function Diag &
  (A)       ! <-- Array of diagonal elements
            ! --> Diagonal matrix
!
! Generation of diagonal matrix.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 05 Nov 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   A(1:)                 ! Array of diagonal elements
!
! Function result:
!
Real(Double) :: &
   Diag(Size(A),Size(A)) ! Diagonal matrix
!----------------------------------------------------------
! Local Scalars:
!
Integer :: i   ! Diagonal index
!----------------------------------------------------------


Diag(:,:) = 0

Do i=1,Size(A)
   Diag(i,i) = A(i)
End Do


End Function Diag




!==========================================================
Function Tensor_Product &
  (X,     & ! <-- First vector
   Y)     & ! <-- Second vector
Result(T)   ! --> Tensor product T(i,j) = X(i)*Y(j)
!
! Tensor product of two vectors.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 24 Apr 2000 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   X(:)     ! First vector
Real(Double), Intent(In) :: &
   Y(:)     ! Second vector
!
! Function result:
!
Real(Double) :: &
   T(Size(X),Size(Y)) ! Tensor product T(i,j) = X(i)*Y(j)
!----------------------------------------------------------
! Local Scalars:
!
Integer :: i   ! Dimension index
!----------------------------------------------------------


Do i=1,Size(X)
   T(i,:) = X(i)*Y(:)
End Do


End Function Tensor_Product




!==========================================================
Subroutine Invert_Matrix &
  (A,       & ! <-> Matrix to invert
   B,       & ! --> Inverted matrix
   ErrCode)   ! --> Error code
!
! Matrix inversion.
!----------------------------------------------------------
! Method:
!   Gauss elimination.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Sep 1998 | Original code.
!   1.1   | 28 Sep 1998 | Fewer error codes.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Inout arguments:
!
Real(Double), Intent(Inout) :: &
   A(1:,1:) ! Matrix to invert.
            ! A is converted to diagonal form:
            ! Det(A) = A(1,1)*...*A(N,N)
!
! Output arguments:
!
Real(Double), Intent(Out)   :: &
   B(1:,1:) ! Inverted matrix
Integer,      Intent(Out)   :: &
   ErrCode  ! Error code:
            !   0 - no error
            !   1 - non-square matrix A
            !   2 - Shape(A) /= Shape(B)
            !   3 - matrix may be degenerated
            !   4 - degenerated matrix
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N          ! Matrix dimension
Integer      :: i, k       ! Matrix indeces
Real(Double) :: Alpha      ! Row combination coefficien
Real(Double) :: Det        ! Matrix determinant
Real(Double) :: Amin, Amax ! Max and min diagonal elements
!
! Local Arrays:
!
Integer      :: m(1)       ! Search index
Real(Double) :: &
   F(Size(A,1))            ! Exchange work space
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ARRAY SHAPE CHECK
!----------------------------------------------------------

If (Size(A,1) /= Size(A,2)) then
   ErrCode = 1
   Return
End If
If (Any(Shape(B) /= Shape(B))) then
   ErrCode = 2
   Return
End If


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

N      = Size(A,1)
B(:,:) = 0

Do i=1,N
   B(i,i) = 1
End Do


!----------------------------------------------------------
! 2. ELIMINATION OF SUB-DIAGONAL ELEMENTS
!----------------------------------------------------------

Lower: Do k=1,N-1

   If (A(k,k) == 0) then
      m(:) = MaxLoc(A(k+1:N,k), A(k+1:N,k) /= 0)
      If (m(1) == 0) then
         Exit Lower
      End If
      i      = k + m(1)
      F(:)   = A(k,:)
      A(k,:) = A(i,:)
      A(i,:) = F(:)
      F(:)   = B(k,:)
      B(k,:) = B(i,:)
      B(i,:) = F(:)
   End If

   Do i=k+1,N
      Alpha  = A(i,k)/A(k,k)
      A(i,:) = A(i,:) - Alpha*A(k,:)
      B(i,:) = B(i,:) - Alpha*B(k,:)
   End Do

End Do Lower


!----------------------------------------------------------
! 3. CHECKING FOR DEGENERATED MATRIX
!----------------------------------------------------------

Det  = 1
Amin = Abs (A(1,1))
Amax = Abs (A(1,1))
Diagonal: Do i=1,N
   Amin = Min(Abs (A(i,i)), Amin)
   Amax = Max(Abs (A(i,i)), Amax)
   Det  = Det*A(i,i)
End Do Diagonal

If (Det == 0) then
   ErrCode = 4
   Return
Else If (Amin/Amax <= 1e-12) then
   ErrCode = 3
Else
   ErrCode = 0
End If


!----------------------------------------------------------
! 4. ELIMINATION OF SUPER-DIAGONAL ELEMENTS
!----------------------------------------------------------

Upper: Do k=N,2,-1
   Do i=k-1,1,-1
      Alpha  = A(i,k)/A(k,k)
      A(i,:) = A(i,:) - Alpha*A(k,:)
      B(i,:) = B(i,:) - Alpha*B(k,:)
   End Do
End Do Upper

Do i=1,N
   B(i,:) = B(i,:)/A(i,i)
End Do


End Subroutine Invert_Matrix



!==========================================================
Subroutine Fast_Invert_Matrix &
  (A,       & ! <-> Matrix to invert
   B)         ! --> Inverted matrix
!
! Matrix inversion without pivoting and error checking.
!----------------------------------------------------------
! Method:
!   Gauss elimination.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Feb 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Inout arguments:
!
Real(Double), Intent(Inout) :: &
   A(1:,1:) ! Matrix to invert.
            ! A is converted to diagonal form:
            ! Det(A) = A(1,1)*...*A(N,N)
!
! Output arguments:
!
Real(Double), Intent(Out)   :: &
   B(1:,1:) ! Inverted matrix
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N          ! Matrix dimension
Integer      :: i, k       ! Matrix indeces
Real(Double) :: Alpha      ! Row combination coefficien
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

N      = Size(A,1)
B(:,:) = 0

Do i=1,N
   B(i,i) = 1
End Do


!----------------------------------------------------------
! 2. ELIMINATION OF SUB-DIAGONAL ELEMENTS
!----------------------------------------------------------

Do k=1,N-1
   Do i=k+1,N
      Alpha  = A(i,k)/A(k,k)
      A(i,:) = A(i,:) - Alpha*A(k,:)
      B(i,:) = B(i,:) - Alpha*B(k,:)
   End Do
End Do


!----------------------------------------------------------
! 3. ELIMINATION OF SUPER-DIAGONAL ELEMENTS
!----------------------------------------------------------

Do k=N,2,-1
   Do i=k-1,1,-1
      Alpha  = A(i,k)/A(k,k)
      A(i,:) = A(i,:) - Alpha*A(k,:)
      B(i,:) = B(i,:) - Alpha*B(k,:)
   End Do
End Do

Do i=1,N
   B(i,:) = B(i,:)/A(i,i)
End Do


End Subroutine Fast_Invert_Matrix



!==========================================================
Subroutine Quasi_Invert &
  (K,      & ! <-- Matrix to quasi-invert
   Q,      & ! --> Quasi-inverse matrix
   ErrCode)  ! --> Error code
!
! Quasi-inverse of a matrix K(dimY, dimX) of system:
! Kx = y
! 1. dimY >= dimX:
!    x = Qy - vector minimizing ||Kx - y||
!    QK = E; Q is left inverse operator.
! 1. dimY <= dimX:
!    x = Qy - solution minimizing ||x||
!    KQ = E; Q is right inverse operator.
!----------------------------------------------------------
! Method:
!                        T  -1 T
! 1. dimY >= dimX: Q = (K K)  K
!
!                       T   T -1
! 2. dimY <= dimX: Q = K (KK )
!
! (See Yu. A. Pyt'ev "Mathematical methods of
!  interpretation of experimental data".)
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 28 Sep 1998 | Original code
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   K(1:,1:)  ! Matrix to quasi-invert
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Q(1:,1:)  ! Quasi-inverse
             ! If Shape(K) = (DimY, DimX),
             ! then Shape(Q) must be (DimX, DimY)
!
Integer, Intent(Out)      :: &
   ErrCode   ! Error code:
             ! 0 - no error
             ! 2 - Q shape misfit
             ! 3 - Rank(K) may be < Min(DimX, DimY)
             ! 4 - Rank(K) < Min(DimX, DimY)
!----------------------------------------------------------
! Local Scalars:
!
Integer :: DimX, DimY  ! Dimension of x and y spaces
Integer :: N           ! Work matrix dimension
!
! Local Arrays:
!
Real(Double), Allocatable :: &
   W(:,:), WI(:,:)     ! Work arrays for inversion
!----------------------------------------------------------


!----------------------------------------------------------
! 0. DIMENSION CHECK, INITIALIZATION
!----------------------------------------------------------

DimX = Size(K,2)
DimY = Size(K,1)

If (Any(Shape(Q) /= (/ DimX, DimY /))) then
   ErrCode = 2
   Return
End If

N = Min(DimX, DimY)
Allocate (W(N,N), WI(N,N))


!----------------------------------------------------------
! 1. QUASI-INVERSION
!----------------------------------------------------------

If (DimY >= DimX) then

   !----------------------------------------------------------
   ! 1.1. LEFT INVERSION
   !----------------------------------------------------------

   !      T          T  -1
   ! W = K K, WI = (K K)

   W = MatMul(Transpose(K), K)
   Call Invert_Matrix(W, WI, ErrCode)

   !       T  -1 T
   ! Q = (K K)  K

   Q = MatMul(WI, Transpose(K))


Else

   !----------------------------------------------------------
   ! 1.2. RIGHT INVERSION
   !----------------------------------------------------------

   !       T          T -1
   ! W = KK , WI = (KK )

   W = MatMul(K, Transpose(K))
   Call Invert_Matrix(W, WI, ErrCode)

   !       T   T -1
   ! Q =  K (KK )

   Q = MatMul(Transpose(K), WI)


End If


!----------------------------------------------------------
! 2. FREEING DYNAMIC MEMORY AND EXIT
!----------------------------------------------------------

Deallocate (W, WI)


End Subroutine Quasi_Invert



!==========================================================
Subroutine Regression &
  (K,      & ! <-- Matrix of basic functions f_j(x_i)
   y,      & ! <-- The data
   a,      & ! --> Regression coefficients
   ErrCode)  ! --> Error code
!
! Linear regression.
!----------------------------------------------------------
! Method:
!   Quasi-inversion of matrix of basic functions:
!   ||sum a  f (x ) - y || = min
!      j   j  j  i     i
!   The solution is a = Q y,
!   where Q is left inverse of K  = f (x )
!                               ij   j  i
!   KQ is the projector to the space spanned on f
!                                                j
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 28 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   K(1:,1:)   ! Matrix of basic functions
              ! K  = f (x )
              !  ij   j  i
!
Real(Double), Intent(In)  :: &
   y(1:)      ! Regression data
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   a(1:)      ! Regression coeffients
!
Integer, Intent(Out)      :: &
   ErrCode    ! Error code:
              !   0 - no error
              !   1 - a-size misfit
              !   2 - y-size misfit
              !   3 - basic functions not independent
!----------------------------------------------------------
! Local Arrays:
!
Real(Double) :: &
   Q(Size(K,2),Size(K,1)) ! Quasi-inverse matrix
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ARRAY SHAPE CHECK
!----------------------------------------------------------

If (Size(K,2) /= Size(a)) then
   ErrCode = 1
   Return
End If
If (Size(K,1) /= Size(y)) then
   ErrCode = 2
   Return
End If


!----------------------------------------------------------
! 1. REGRESSION
!----------------------------------------------------------

Call Quasi_Invert(K, Q, ErrCode)
a = MatMul(Q, y)


End Subroutine Regression



!==========================================================
Subroutine Basic_Polynomials &
  (x,      & ! <-- Argument x grid
   K,      & ! --> Polynomial matrix
   ErrCode)  ! --> Error code
!
! Generation of matrix of basic polynomials for
! polynomial regression.
!----------------------------------------------------------
! Method:
!                  j
!   K  = f (x ) = x  , j=0..UBound(K,2)
!    ij   j  i     i
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 29 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   x(1:)      ! Grid of argument x
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   K(1:,0:)   ! Matrix of basic polynomials
              ! Size(K,1) must be equal to Size(x)
!
Integer, Intent(Out)      :: &
   ErrCode    ! Error code:
              !   0 - no error
              !   1 - x-size misfit
!----------------------------------------------------------
! Local Scalars:
!
Integer     :: j  ! Polynomial power
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ARRAY SIZE CHECK
!----------------------------------------------------------

If (Size(K,1) /= Size(X)) then
   ErrCode = 1
   Return
End If

!----------------------------------------------------------
! 2. GENERATION OF BASIC FUNCTION MATRIX
!----------------------------------------------------------

K(:,0) = 1
Do j=1,UBound(K,2)
   K(:,j) = x(:)**j
End Do


End Subroutine Basic_Polynomials



!==========================================================
Subroutine Basic_Splines &
  (x,      & ! <-- Argument x grid
   xs,     & ! <-- x-grid of delta-splines
   K,      & ! --> Delta-spline matrix
   ErrCode)  ! --> Error code
!
! Generation of matrix of basic polynomials for
! polynomial regression.
!----------------------------------------------------------
! Method:
!   K  = f (x ) = S (x ), j=0..UBound(K,2)
!    ij   j  i     j  i
!   S  - Delta-splines: S (xs ) = Delta
!    j                   j   k         jk
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 29 Sep 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Init_Spline, Spline
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   x(1:)      ! Grid of argument x
!
Real(Double), Intent(In)  :: &
   xs(1:)     ! x-grid of delta-splines
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   K(1:,1:)   ! Matrix of basic polynomials
              ! Size(K,1) must be equal to Size(x)
              ! Size(K,2) must be equal to Size(xs)
!
Integer, Intent(Out)      :: &
   ErrCode    ! Error code:
              !   0 - no error
              !   1 - xs is not monotonously increasing
              !   2 - x-size misfit
              !   3 - xs-size misfit
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: j   ! Delta-spline number
Integer      :: i   ! x index
!
! Local Arrays:
!
Real(Double) :: &
   S(Size(xs)),   & ! Delta-spline
   D2S(Size(xs))    ! Delta-spline second derivative
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ARRAY SIZE CHECK
!----------------------------------------------------------

If (Size(K,1) /= Size(x)) then
   ErrCode = 1
   Return
End If
If (Size(K,2) /= Size(xs)) then
   ErrCode = 2
   Return
End If


!----------------------------------------------------------
! 1. GENERATION OF BASIC FUNCTION MATRIX
!----------------------------------------------------------

S(:) = 0

Do j = 1, Size(K,2)
   S(j) = 1
   Call Init_Spline(xs, S, D2S, ErrCode)
   Do i = 1, Size(x)
      Call Spline (xs, S, D2S, x(i),  K(i,j))
   End Do
   S(j) = 0
End Do


End Subroutine Basic_Splines



End Module Matrix


