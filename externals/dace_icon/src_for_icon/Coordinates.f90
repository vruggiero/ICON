!
!+ GNSS Radio occultation operator: cartesian and geocentric coordinates
!
MODULE Coordinates
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Definition and transformations of cartesian
!   and geocentric coordinates, vetor analysis.
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
! Language: Fortran 95.
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
! Module Coordinates
!
! Definition and transformations of cartesian
! and geocentric coordinates, vetor analysis.
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | Definition of type Cartesian.
!   1.2   | 17 Sep 1998 | Use of new versions of Defauls
!         |             |   and Earth modules.
!   2.0   | 18 Sep 1998 | Earth-related subprograms
!         |             |   in module Earth.
!   2.1   | 18 Sep 1998 | Operator overloadings.
!   2.2   | 24 Sep 1998 | Unary vector minus.
!   2.3   | 22 Oct 1998 | Vector_Normed, Vector_Norm
!   2.4   | 18 Feb 1999 ! Integer_x_Vector
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
Private
Public :: Cartesian, Spherical, Spher_from_Cart, Rotate, Vector_Angle, &
          Vector_Norm, Vector_Normed, &
          operator(.x.), operator (*), operator (-), operator (/),     &
          operator (+), assignment(=)
!----------------------------------------------------------
! Public Type Definitions:
!
Type Cartesian            ! Cartesian coordinates
  Real(Double) :: X(3)      ! XYZ-coordinates
End Type

Type Spherical            ! Spherical coordinates:
  Real(Double) :: R         ! Radius (km)
  Real(Double) :: Phi       ! Latitude from North Pole (radian)
  Real(Double) :: Lambda    ! Longitude (radian)
End Type
!
! Operator definitions:
!
Interface Assignment(=)
   Module Procedure Vector_to_Array
   Module Procedure Array_to_Vector
   Module Procedure Real_to_Vector
   Module Procedure Integer_to_Vector
End Interface
!
Interface Operator(.x.)
   Module Procedure Vector_Product
End Interface
!
Interface Operator(*)
   Module Procedure Scalar_Product
End Interface
!
Interface Operator(+)
   Module Procedure Vector_Sum
End Interface
!
Interface Operator(-)
   Module Procedure Vector_Dif
End Interface
!
Interface Operator(-)
   Module Procedure Vector_Minus
End Interface
!
Interface Operator(*)
   Module Procedure Real_x_Vector
   Module Procedure Vector_x_Real
   Module Procedure Integer_x_Vector
End Interface
!
Interface Operator(/)
   Module Procedure Vector_over_Scalar
End Interface
!----------------------------------------------------------
Contains



!==========================================================
Subroutine Real_to_Vector &
  (V,   & ! --> Vector
   A)     ! <-- Scalar
!
! Conversion of scalar to vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Outut arguments:
!
Type(Cartesian), Intent(Out) :: &
   V     ! Resulting vector
!
! Input arguments:
!
Real(Double), Intent(In)     :: &
   A     ! Scalar
!----------------------------------------------------------


V%X(:) = A


End Subroutine Real_to_Vector



!==========================================================
Subroutine Integer_to_Vector &
  (V,   & ! --> Vector
   A)     ! <-- Scalar
!
! Conversion of scalar to vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Mar 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Outut arguments:
!
Type(Cartesian), Intent(Out) :: &
   V     ! Resulting vector
!
! Input arguments:
!
Integer, Intent(In)     :: &
   A     ! Scalar
!----------------------------------------------------------


V%X(:) = A


End Subroutine Integer_to_Vector



!==========================================================
Subroutine Vector_to_Array &
  (A,   & ! --> Array
   V)     ! <-- Vector
!
! Conversion of vector to array(3).
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Outut arguments:
!
Real(Double), Intent(Out)   :: &
   A(3)  ! Resulting array
!
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   V     ! Vector
!----------------------------------------------------------


A(:) = V%X(:)


End Subroutine Vector_to_Array



!==========================================================
Subroutine Array_to_Vector &
  (V,   & ! --> Vector
   A)     ! <-- Array
!
! Conversion of array(3) to vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Outut arguments:
!
Type(Cartesian), Intent(Out) :: &
   V     ! Resulting vector
!
! Input arguments:
!
Real(Double), Intent(In)     :: &
   A(3)  ! Array
!----------------------------------------------------------


V%X(:) = A(:)


End Subroutine Array_to_Vector



!==========================================================
Function Vector_Sum &
  (X,   & ! <-- First vector
   Y)     ! <-- Second vector
          ! --> Sum of vectors
!
! Calculation of sum of two vectors.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X, Y    ! Vectors to summate
!
! Function value:
!
Type(Cartesian)            :: &
   Vector_Sum  ! Sum of vectors
!----------------------------------------------------------


Vector_Sum%X(:) = X%X(:) + Y%X(:)


End Function Vector_Sum



!==========================================================
Function Vector_Dif &
  (X,   & ! <-- First vector
   Y)     ! <-- Second vector
          ! --> Sum of vectors
!
! Calculation of difference of two vectors.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X, Y    ! Vectors to find difference
!
! Function value:
!
Type(Cartesian)            :: &
   Vector_Dif  ! Difference of vectors
!----------------------------------------------------------


Vector_Dif%X(:) = X%X(:) - Y%X(:)


End Function Vector_Dif



!==========================================================
Function Vector_Minus &
  (X)     ! <-- Vector
          ! --> Negated vector
!
! Calculation of negated vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 24 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X    ! Vector to negate
!
! Function value:
!
Type(Cartesian)   :: &
   Vector_Minus  ! Difference of vectors
!----------------------------------------------------------


Vector_Minus%X(:) = -X%X(:)


End Function Vector_Minus



!==========================================================
Function Real_x_Vector &
  (A,  &  ! <-- Scalar
   X)     ! <-- Vector
          ! --> Product of scalar and vector
!
! Calculation of product of scalar and vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)    :: &
   A    ! Scalar
Type(Cartesian), Intent(In) :: &
   X    ! Vector
!
! Function value:
!
Type(Cartesian)            :: &
   Real_x_Vector  ! Product of scalar and vector
!----------------------------------------------------------


Real_x_Vector%X(:) = A*X%X(:)


End Function Real_x_Vector



!==========================================================
Function Integer_x_Vector &
  (A,  &  ! <-- Scalar
   X)     ! <-- Vector
          ! --> Product of scalar and vector
!
! Calculation of product of scalar and vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)         :: &
   A    ! Scalar
Type(Cartesian), Intent(In) :: &
   X    ! Vector
!
! Function value:
!
Type(Cartesian)            :: &
   Integer_x_Vector  ! Product of scalar and vector
!----------------------------------------------------------


Integer_x_Vector%X(:) = A*X%X(:)


End Function Integer_x_Vector



!==========================================================
Function Vector_x_Real &
  (X,   & ! <-- Vector
   A)     ! <-- Scalar
          ! --> Product of vector and scalar
!
! Calculation of product of vector and scalar.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 19 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X    ! Vector
Real(Double), Intent(In)    :: &
   A    ! Scalar
!
! Function value:
!
Type(Cartesian)            :: &
   Vector_x_Real  ! Product of vector and scalar
!----------------------------------------------------------


Vector_x_Real%X(:) = X%X(:)*A


End Function Vector_x_Real



!==========================================================
Function Vector_over_Scalar &
  (X,  & ! <-- Vector
   A)    ! <-- Scalar
         ! --> Vector divided by scalar
!
! Division of vector by scalar.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X    ! Vector
Real(Double), Intent(In)    :: &
   A    ! Scalar
!
! Function value:
!
Type(Cartesian)            :: &
   Vector_over_Scalar  ! Vector divided by scalar
!----------------------------------------------------------


Vector_over_Scalar%X(:) = X%X(:)/A


End Function Vector_over_Scalar



!==========================================================
Function Rotate &
  (X,   &     ! <-- The vector to rotate
   A,   &     ! <-- The rotation axis
   Phi)       ! <-- The rotation angle
              ! --> The rotated vector
!
! Rotation of a vector in cartesian coordinates
! around a given axis by a given angle.
!----------------------------------------------------------
! Method:
!   N*(N,X) + [N,X]*Sin(Phi) + (X-N*(N,X))*Cos(Phi),
!   where N=A/|A|
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Oct 1997 | Original code.
!   2.0   | 17 Sep 1998 | Function.
!   2.1   | 19 Sep 1998 | Simplified formula.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In) :: &
   X    ! Vector to rotate
Type (Cartesian), Intent(In) :: &
   A    ! Rotation axis
Real(Double), Intent(In)     :: &
   Phi  ! Rotation angle (rad)
!
! Function result:
!
Type (Cartesian)            :: &
   Rotate  ! Rotated vector
!----------------------------------------------------------
! Local Scalars:
!
Type(Cartesian) :: N  ! Normed rotation axis
!----------------------------------------------------------


N = A/Sqrt(A*A)

Rotate = N*(N*X) + (N.x.X)*Sin(Phi) + (X-N*(N*X))*Cos(Phi)


End Function Rotate



!==========================================================
Function Scalar_Product &
  (X,   & ! <-- First vector
   Y)     ! <-- Second vector
          ! --> Scalar product
!
! Calculation of scalar product of two vectors.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 17 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X, Y    ! Vectors to find scalar product of
!
! Function value:
!
Real(Double)               :: &
   Scalar_Product  ! Scalar product
!----------------------------------------------------------


Scalar_Product = Dot_Product(X%X(:), Y%X(:))


End Function Scalar_Product



!==========================================================
Function Vector_Product &
  (X, &    ! <-- First vector
   Y)      ! <-- Second vector
           ! --> Vector product
!
! Calculation of the vector product of two vectors.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorubnov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Oct 1997 | Original code.
!   2.0   | 17 Sep 1998 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In) :: &
   X,  Y  !     Vectors to multiply
!
! Function result:
!
Type (Cartesian)            :: &
   Vector_Product  ! Vector product of X and Y
!----------------------------------------------------------


Vector_Product%X(:) =                       &
         (/ X%X(2)*Y%X(3) - X%X(3)*Y%X(2),  &
            X%X(3)*Y%X(1) - X%X(1)*Y%X(3),  &
            X%X(1)*Y%X(2) - X%X(2)*Y%X(1) /)


End Function Vector_Product



!==========================================================
Function Vector_Normed &
  (X)      ! <-- Vector
           ! --> Normed vector
!
! Calculation of normed vector.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorubnov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 Oct 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In) :: &
   X  !     Vectors to norm
!
! Function result:
!
Type (Cartesian)            :: &
   Vector_Normed  ! Normed vector X
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: N  ! Vector norm
!----------------------------------------------------------

N = Sqrt(Sum(X%X(:)**2))

Vector_Normed%X(:) = X%X(:)/N


End Function Vector_Normed



!==========================================================
Function Vector_Norm &
  (X)      ! <-- Vector
           ! --> Vector norm
!
! Calculation of vector norm.
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorubnov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 22 Oct 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In) :: &
   X  !     Vectors to find norm of
!
! Function result:
!
Real(Double)               :: &
   Vector_Norm  ! Vector norm
!----------------------------------------------------------


Vector_Norm = Sqrt(Sum(X%X(:)**2))


End Function Vector_Norm



!==========================================================
Function Vector_Angle &
  (X, Y,   &  ! <-- Vectors to find the angle between
   A)         ! <~~ Orientation axis
              ! --> The angle between the vectors
!
! Calculation of the angle between two vectors.
!----------------------------------------------------------
! Method:
!   N = A/|A|
!   Alpha = [N,X], Beta = X-N*(N,X), Gamma=Y-N*(N,Y)
!   Sin(Phi)=(Alpha,Gamma)/(Alpha,Alpha)
!   Cos(Phi)=(Beta,Gamma)/(Beta,Beta)
!   Vector_Angle=ATan2(Sin(Phi),Cos(Phi))
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Oct 1997 | Original code.
!   2.0   | 18 Sep 1998 | Orientation axis.
!   2.1   | 18 Sep 1998 | Optimal solution for Sin&Cos
!   2.2   | 26 Feb 1999 | Correction for zero angle
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In)            :: &
   X, Y  ! Vectors to find the angle between
Type (Cartesian), Intent(In), Optional  :: &
   A     ! Orientation axis
         ! Default: A = [X,Y]/|[X,Y]|
!
! Function result
!
Real(Double)                           :: &
   Vector_Angle  ! Angle between X and Y.
!                ! The sign of the angle is the sign of
!                ! rotation direction from X to Y
!                ! around axis A.
!----------------------------------------------------------
! Local Scalars:
!
Type(Cartesian) :: N, Alpha, Beta, Gamma
Real(Double)    :: NN
!----------------------------------------------------------


If (Present(A)) then
   N = A
Else
   N = X .x. Y
End If

NN = N*N
If (NN == 0) then
   Vector_Angle = 0
   Return
End If
N  = N/Sqrt(NN)

Alpha = N .x. X
Beta  = X - (N*X)*N
Gamma = Y - (N*Y)*N

Vector_Angle = ATan2(Alpha*Gamma, Beta*Gamma)


End Function Vector_Angle



!==========================================================
Function Spher_from_Cart &
  (X)    & ! <-- Cartesian coordinates
Result(S)  ! --> Spherical coordinates
!
! Conversion from Cartesian to spherical coordinates.
!----------------------------------------------------------
! Method:
!    Standard formulas
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 08 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 17 Sep 1998 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In) :: &
   X   ! Cartesian coordinates of a point
!
! Function result:
!
Type (Spherical)            :: &
   S   ! Spherical coordinates of the point
!----------------------------------------------------------


S%R      = Sqrt(Sum(X%X(:)**2))
S%Phi    = ACos(X%X(3)/S%R)
S%Lambda = ATan2(X%X(2),X%X(1))


End Function Spher_from_Cart



!==========================================================
Function Cart_from_Spher &
  (S)    & ! <-- Spherical coordinates
Result(C)  ! --> Cartesian coordinates
!
! Conversion from spherical to Cartesian coordinates
!----------------------------------------------------------
! Method:
!    Standard formulas
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 08 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Spherical), Intent(In) :: &
   S   ! Spherical coordinates of a point
!
! Function result:
!
Type (Cartesian)            :: &
   C   ! Cartesian coordinates of the point
!----------------------------------------------------------


C%X(3) = S%R*Cos(S%Phi)
C%X(2) = S%R*Sin(S%Phi)*Sin(S%Lambda)
C%X(1) = S%R*Sin(S%Phi)*Cos(S%Lambda)


End Function Cart_from_Spher



!==========================================================
Function Cart_from_Radial &
  (RGrad,   & ! <-- Radial gradient df/d|x|
   X)       & ! <-- Cartesian coordinates of point x
Result(Grad)  ! --> Cartesian gradient df/dx
!
! Conversion of the radial gradient of a radial scalar
! field to a gradient vector in Cartesian coordinates.
!----------------------------------------------------------
! Method:
!   Standard formulas
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)     :: &
   RGrad  ! Radial gradient gradient of a scalar field
Type (Cartesian), Intent(In) :: &
   X      ! Cartesian cooradinates of the point
!
! Function result:
!
Type (Cartesian)            :: &
   Grad  ! Cartesian gradient vector
!        ! of the field at the point
!----------------------------------------------------------


Grad%X(:) = RGrad*X%X(:)/Sqrt(Sum(X%X(:)**2))


End Function Cart_from_Radial



End Module Coordinates

