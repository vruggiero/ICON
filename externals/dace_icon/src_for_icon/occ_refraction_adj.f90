!
!+ GNSS Radio occultation operator: Adjoint version of Occ_Refraction module
!
MODULE Occ_refraction_adj
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Adjoint version of Occ_Refraction module.
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
! Module Occ_refraction_adj
!
! Adjoint version of Occ_Refraction module.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 03 May 2000 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
Private
Public :: Ray_Derivatives, Ray_Refraction_adj
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Ray_Refraction_adj &
  (Z0,        & ! <-- Initial ray point (X0,U0)
   ZN,        & ! <-- Final ray point (XN,UN)
   V0,        & ! <-- Transmitter velocity [km/s]
   VN,        & ! <-- Receiver velocity [km/s]
   Eps,       & ! --> Refraction angle [rad]
   P,         & ! --> Impact parameter [km]
   E_Z0,      & ! --> d(Eps)/d(Z0)
   E_ZN,      & ! --> d(Eps)/d(ZN)
   P_Z0,      & ! --> d(P)/d(Z0)
   P_ZN)        ! --> d(P)/d(ZN)
!
! Calculation of refraction angle and impact parameter
! as a functional of a ray: Adjoint version.
!----------------------------------------------------------
! Method:
!   Calculation of Doppler frequency shift from ray end
!   points and direction, calculation of refraction angle
!   and impact parameter from Doppler shift under the
!   spherical symmetry approximation.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 28 Apr 2000 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,      &
! Imported Operators:
    Operator(-),    &
    Operator(/),    &
    Operator(+),    &
    Operator(*),    &
    Operator(.x.),  &
    Assignment(=)
!
!Use Occ_Refraction, only: &
! Imported Routines:
!    Doppler_to_Refraction
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)    :: &
   Z0(6)        ! Initial ray point (X0,U0)
Real(Double), Intent(In)    :: &
   ZN(6)        ! Final ray point (XN,UN)
Type(Cartesian), Intent(In) :: &
   V0           ! Transmitter velocity [km/s]
Type(Cartesian), Intent(In) :: &
   VN           ! Receiver velocity [km/s]
!
! Output arguments:
!
Real(Double), Intent(Out)   :: &
   Eps          ! Refraction angle [rad]
Real(Double), Intent(Out)   :: &
   P            ! Impact parameter [km]
Real(Double), Intent(Out)   :: &
   E_Z0(6)      ! d(Eps)/d(Z0)
Real(Double), Intent(Out)   :: &
   E_ZN(6)      ! d(Eps)/d(ZN)
Real(Double), Intent(Out)   :: &
   P_Z0(6)      ! d(P)/d(Z0)
Real(Double), Intent(Out)   :: &
   P_ZN(6)      ! d(P)/d(ZN)
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   c = 300000.0_Double     ! Light velocity (km/s)
! Local Scalars:
!
Type(Cartesian) :: X0     ! Transmitter position
Type(Cartesian) :: U0     ! Ray direction at transmitter
Type(Cartesian) :: XN     ! Receiver position
Type(Cartesian) :: UN     ! Ray direction at receiver
Real(Double)    :: D      ! Doppler freaquency shift
Real(Double)    :: P_d    ! d(P)/d(d)
Real(Double)    :: E_d    ! d(Eps)/d(d)
!
! Local Arrays:
!
Real(Double)  :: &
   P_XRT(6),   & ! d(P)/d(XN,X0)
   E_XRT(6),   & ! d(Eps)/d(XN,X0)
   D_U0(3),    & ! d(D)/d(U0)
   D_UN(3)       ! d(D)/d(UN)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF DOPPLER SHIFT
!----------------------------------------------------------

X0 = Cartesian(Z0(1:3))
U0 = Cartesian(Z0(4:6))
XN = Cartesian(ZN(1:3))
UN = Cartesian(ZN(4:6))

D = ((V0*U0) - (VN*UN))/(c - (V0*U0))

D_UN(:) = -VN/(c - (VN*UN))
D_U0(:) =  V0*(c - (VN*UN))/(c - (V0*U0))**2


!----------------------------------------------------------
! 2. CALCULATION OF REFRACTION ANGLE AND IMPACT PARAMETER
!----------------------------------------------------------

Call Doppler_to_Refraction_adj &
  (X0,     & ! <-- Transmitter position
   V0,     & ! <-- Transmitter velocity
   XN,     & ! <-- Receiver position
   VN,     & ! <-- Receiver velocity
   D,      & ! <-- Relative Doppler frequency shift
   P,      & ! --> Impact parameter
   Eps,    & ! --> Refraction angle
   P_D,    & ! --> d(P)/d(D)
   P_XRT,  & ! --> d(P)/d(XR,XT)
   E_D,    & ! --> d(Eps)/d(D)
   E_XRT)    ! --> d(Eps)/d(XR,XT)


E_Z0(1:3) = E_XRT(4:6)
E_Z0(4:6) = E_D*D_U0(:)

E_ZN(1:3) = E_XRT(1:3)
E_ZN(4:6) = E_D*D_UN(:)

P_Z0(1:3) = P_XRT(4:6)
P_Z0(4:6) = P_D*D_U0(:)

P_ZN(1:3) = P_XRT(1:3)
P_ZN(4:6) = P_D*D_UN(:)


End Subroutine Ray_Refraction_adj



!==========================================================
Subroutine Doppler_to_Refraction_adj &
  (XT,     & ! <-- Transmitter position
   VT,     & ! <-- Transmitter velocity
   XR,     & ! <-- Receiver position
   VR,     & ! <-- Receiver velocity
   D,      & ! <-- Relative Doppler frequency shift
   P,      & ! --> Impact parameter
   Eps,    & ! --> Refraction angle
   P_d,    & ! --> d(P)/d(d)
   P_XRT,  & ! --> d(P)/d(XR,XT)
   E_d,    & ! --> d(Eps)/d(d)
   E_XRT)    ! --> d(Eps)/d(XR,XT)
!
! Calculation of refraction angle and impact parameter
! from relative Doppler frequency shift: Adjoint version.
!----------------------------------------------------------
! Method:
!   Iterative solution of system:
!   (c - (VR, UR))/(c - (VT, UT)) - 1 = D
!   [XR, UR] - [XT, UT] = 0
!   (UR, UR) = 1
!   (UT, UT) = 1
!   where UR and UT are the ray directions
!   at the receiver and transmitter respectively.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Sep 1998 | Basic non-adjoint version.
!   1.0   | 28 Apr 2000 | Adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,    &
! Imported Routines:
    Vector_Normed,  &
! Imported Operators:
    Operator(-),    &
    Operator(/),    &
    Operator(+),    &
    Operator(*),    &
    Operator(.x.),  &
    Assignment(=)
!
Use Matrix, only: &
! Imported Routines:
    Invert_Matrix
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   XT      ! Position of transmitter (km)
Type(Cartesian), Intent(In) :: &
   VT      ! Velocity of transmitter (km/s)
Type(Cartesian), Intent(In) :: &
   XR      ! Position of receiver (km)
Type(Cartesian), Intent(In) :: &
   VR      ! Velocity of receiver (km/s)
Real(Double), Intent(In)    :: &
   d       ! Relative Doppler frequency shift
!
! Output arguments:
!
Real(Double), Intent(Out)    :: &
   P         ! Impact parameter (km)
Real(Double), Intent(Out)    :: &
   Eps       ! Signed refraction angle (rad)
Real(Double), Intent(Out)    :: &
   P_d       ! d(P)/d(d)
Real(Double), Intent(Out)    :: &
   P_XRT(6)  ! d(P)/d(XR,XT)
Real(Double), Intent(Out)    :: &
   E_d       ! d(Eps)/d(d)
Real(Double), Intent(Out)    :: &
   E_XRT(6)  ! d(Eps)/d(XR,XT)
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter      :: &
   N_it = 10               ! Number of iterations
Real(Double), Parameter :: &
   c = 300000.0_Double     ! Light velocity (km/s)
Real(Double), Parameter :: &
   E(3,3,3) =   &          ! Antisymmetrical tensor
     Reshape  &
      ((/0,  0,  0,  0,  0,  1,  0, -1,  0,   &
         0,  0, -1,  0,  0,  0,  1,  0,  0,   &
         0,  1,  0, -1,  0,  0,  0,  0,  0/), &
         Shape = (/3,3,3/), Order = (/3,2,1/))
!
! Local Scalars:
!
Type(Cartesian) :: UR        ! Ray direction at receiver
Type(Cartesian) :: UT        ! Ray direction at transmitter
Integer         :: ErrCode   ! Inversion error code
Integer         :: i,j       ! Matrix indices
Integer         :: K         ! Iteration number
!
! Local Arrays:
!
Real(Double)    :: DZ(6)        ! Perturbation of (UR, UT)
Real(Double)    :: DF(6)        ! Discreapancy vector
Real(Double)    :: A(6,6)       ! Matrix of linearized system
                                ! A*DZ = DF
Real(Double)    :: AI(6,6)      ! Inverted system matrix
Real(Double)    :: E_URT(6)     ! d(Eps)/d(UR,UT)
Real(Double)    :: P_UR(3)      ! d(P)/d(UR)
Real(Double)    :: P_XR(3)      ! d(P)/d(XR)
Real(Double)    :: URT_D(6)     ! d(UR)/d(D)
Real(Double)    :: BXR(3,3)     ! Influence of XR
Real(Double)    :: BXT(3,3)     ! Influence of XT
Real(Double)    :: B(6,6)       ! Influence of XR and XT
Real(Double)    :: URT_XRT(6,6) ! d(UR,UT)/d(XR,XT)
!----------------------------------------------------------



!----------------------------------------------------------
! 1. INITIAL APPROXIMATION
!----------------------------------------------------------

! Straight-line XT --> XR direction

UR = Vector_Normed(XR-XT)
UT = UR


!----------------------------------------------------------
! 2. ITERATIVE SOLUTION
!----------------------------------------------------------

Iterations: Do K=1,N_it

   !--- 2.1. Calculation of discreapancy

   DF(1)   = D - ((VT*UT) - (VR*UR))/(c - (VT*UT))
   DF(2:4) = (XT.x.UT) - (XR.x.UR)
   DF(5)   = 1 - (UR*UR)
   DF(6)   = 1 - (UT*UT)

   !--- 2.2 Calculation of the matrix of the linearized system

   A(:,:)   = 0
   A(1,1:3) = -VR/(c - (VR*UR))
   A(1,4:6) =  VT*(c - (VR*UR))/(c - (VT*UT))**2
   Do i=1,3
      Do j=1,3
         A(i+1,j)   =  Sum(E(i,:,j)*XR%X(:))
         A(i+1,j+3) = -Sum(E(i,:,j)*XT%X(:))
      End Do
   End Do
   A(5,1:3) = 2.0_Double*UR
   A(6,4:6) = 2.0_Double*UT

   !--- 2.3. Solution of the linearized system

   Call Invert_Matrix(A, AI, ErrCode)
   DZ = MatMul(AI, DF)

   !--- 2.4. Calculation of next approximation

   UR = UR + Cartesian(DZ(1:3))
   UT = UT + Cartesian(DZ(4:6))

End Do Iterations


!----------------------------------------------------------
! 3. CALCULATION OF REFRACTION ANGLE AND IMPACT PARAMETER
!----------------------------------------------------------


Call Vector_Angle_adj &
  (UT, UR,      & ! <-- Vectors to find the angle between
   XT.x.XR,     & ! <-- Orientation axis
   Eps,         & ! --> The angle between the vectors
   E_URT(4:6),  & ! --> d(AXY)/d(X)
   E_URT(1:3))    ! --> d(AXY)/d(Y)

P = Sqrt((XR.x.UR)*(XR.x.UR))


!----------------------------------------------------------
! 4. ADJOINT CALCULATIONS
!----------------------------------------------------------


!--- 4.1. Derivatives of UR and UT

URT_D(:) = AI(:,1)

Do i=1,3
   Do j=1,3
      BXR(i,j) = -Sum(E(i,j,:)*UR%X(:))
      BXT(i,j) =  Sum(E(i,j,:)*UT%X(:))
   End Do
End Do

B(:,:)       = 0
B(2:4,1:3)   = BXR(:,:)
B(2:4,4:6)   = BXT(:,:)
URT_XRT(:,:) = MatMul(AI(:,:),B(:,:))


!--- 4.2. Derivatives of Eps

E_D      = Sum(E_URT(:)*URT_D(:))
E_XRT(:) = MatMul(E_URT(:),URT_XRT(:,:))


!--- 4.3. Derivatives of P

P_UR(:)    = -(XR.x.(XR.x.UR))/P
P_XR(:)    = -(UR.x.(UR.x.XR))/P

P_D        = Sum(P_UR(:)*URT_D(1:3))

P_XRT(:)   = MatMul(P_UR(:),URT_XRT(1:3,:))
P_XRT(1:3) = P_XRT(1:3) + P_XR(:)


End Subroutine Doppler_to_Refraction_adj



!==========================================================
Subroutine Vector_Angle_adj &
  (X, Y,    & ! <-- Vectors to find the angle between
   A,       & ! <-- Orientation axis
   AXY,     & ! --> The angle between the vectors
   AXY_X,   & ! --> d(AXY)/d(X)
   AXY_Y)     ! --> d(AXY)/d(Y)
!
! Calculation of the angle between two vectors:
! Ajodint version.
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
!   2.2   | 26 Feb 1999 | Basic non-adjoint version.
!   1.0   | 28 Apr 2000 | Adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In)  :: &
   X, Y       ! Vectors to find the angle between
Type (Cartesian), Intent(In)  :: &
   A          ! Orientation axis
!
! Output arguments:
!
Real(Double), Intent(Out)     :: &
   AXY        ! Angle between X and Y.
!             ! The sign of the angle is the sign of
!             ! rotation direction from X to Y
!             ! around axis A.
Real(Double), Intent(Out)     :: &
   AXY_X(3)   ! d(AXY)/d(X)
Real(Double), Intent(Out)     :: &
   AXY_Y(3)   ! d(AXY)/d(Y)
!----------------------------------------------------------
! Local Scalars:
!
Type(Cartesian) :: N, Alpha, Beta, Gamma
Real(Double)    :: NN
Real(Double)    :: AG, BG
!----------------------------------------------------------


!----------------------------------------------------------
! 1. NON-ADJOINT CALCULATIONS
!----------------------------------------------------------

NN = A*A
If (NN == 0) then
   AXY = 0
   Return
End If
N  = A/Sqrt(NN)

Alpha = N .x. X
Beta  = X - (N*X)*N
Gamma = Y - (N*Y)*N

AG = Alpha*Gamma
BG = Beta*Gamma

AXY = ATan2(AG, BG)


!----------------------------------------------------------
! 2. ADJOINT CALCULATIONS
!----------------------------------------------------------

AXY_X = -(AG*Gamma - BG*(Y.x.N))/(AG**2 + BG**2)
AXY_Y = -(AG*Beta  - BG*Alpha)/(AG**2 + BG**2)


End Subroutine Vector_Angle_adj



!==========================================================
Subroutine Ray_Derivatives &
  (Z0,        & ! <-- Initial ray point (X0,U0)
   ZN,        & ! <-- Final ray point (XN,UN)
   ZN_Z0,     & ! <-- d(ZN)/d(Z0)
   P_Z0,      & ! <-- d(P)/d(Z0)
   P_ZN,      & ! <-- d(P)/d(ZN)
   E_Z0,      & ! <-- d(E)/d(Z0)
   E_ZN,      & ! <-- d(D)/d(ZN)
   FZ0_P,     & ! --> D(Z0)/D(P)
   FE_ZN)       ! --> D(E)/D(ZN)
!
! Calculation of full derivatives D(Z0)/D(P) and
! D(E)/D(ZN).
!----------------------------------------------------------
! Method:
!   Solution for vertical variations of initial ray
!   direction.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Apr 2000 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Vector_Normed, &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Matrix, only: &
! Imported Routines:
    Invert_Matrix
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   Z0(6)        ! Initial ray point (X0,U0)
Real(Double), Intent(In) :: &
   ZN(6)        ! Final ray point (XN,UN)
Real(Double), Intent(In) :: &
   ZN_Z0(6,6)   ! d(ZN)/d(Z0)
Real(Double), Intent(In) :: &
   P_Z0(6)      ! d(P)/d(Z0)
Real(Double), Intent(In) :: &
   P_ZN(6)      ! d(P)/d(ZN)
Real(Double), Intent(In) :: &
   E_Z0(6)      ! d(E)/d(Z0)
Real(Double), Intent(In) :: &
   E_ZN(6)      ! d(E)/d(ZN)
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   FZ0_P(6)     ! D(Z0)/D(P)
Real(Double), Intent(Out) :: &
   FE_ZN(6)     ! D(E)/D(ZN)
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   E3(3,3) = Reshape((/1, 0, 0,  &
                       0, 1, 0,  &
                       0, 0, 1 /),  &
                     (/3,3/))     ! Unit 3x3 matrix
!
! Local Scalars:
!
Type(Cartesian) :: X0     ! Initial ray point
Type(Cartesian) :: U0     ! Initial ray direction
Type(Cartesian) :: XN     ! Final ray point
Type(Cartesian) :: UN     ! Final ray direction
Integer         :: Stat   ! Matrix inversion status
Real(Double)    :: FE_P   ! D(E)/D(P)
!
! Local Arrays:
!
Real(Double)    :: FE_Z0(6)  ! D(E)/D(Z0)
Real(Double)    :: A(6,6)    ! System matrix
Real(Double)    :: AI(6,6)   ! Inverted system matrix
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

X0 = Z0(1:3)
U0 = Z0(4:6)
XN = ZN(1:3)
UN = ZN(4:6)


!----------------------------------------------------------
! 2. CALCULATION OF SYSTEM MATRIX
!----------------------------------------------------------

A(1,:)     = P_Z0(:) + MatMul(P_ZN(:),ZN_Z0(:,:))
A(2:4,1:3) = E3(:,:)
A(2:4,4:6) = 0
A(5,1:3)   = 0
A(5,4:6)   = Vector_Normed(X0.x.XN)
A(6,1:3)   = 0
A(6,4:6)   = U0


!----------------------------------------------------------
! 3. SOLUTION OF SYSTEM
!----------------------------------------------------------

Call Invert_Matrix(A, AI, Stat)

FZ0_P(:) = AI(:,1)


!----------------------------------------------------------
! 4. CALCULATION OF D(E)/D(ZN)
!----------------------------------------------------------

FE_Z0(:) = E_Z0(:) + MatMul(E_ZN(:), ZN_Z0(:,:))
FE_P     = Sum(FE_Z0(:)*FZ0_P(:))

FE_ZN(:) = E_ZN(:) - FE_P*P_ZN(:)


End Subroutine Ray_Derivatives



End Module Occ_refraction_adj


