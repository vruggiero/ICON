!
!+ GNSS Radio occultation observation operator: Coordinate transforms
!
MODULE Occ_Coordinates
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Coordinate transforms used in processing of
!   a radio occultation and calculation.
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
! Module Occ_Coordinates
!
! Coordinate transforms used in processing of
! a radio occultation and calculation.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 19 Mar 1999 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, Pi
!----------------------------------------------------------
Implicit None
Private
Public :: Plane_Coordinates, Satellite_Velocities, Perigee

!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Occ_Geometry &
  (Year,      & ! <-- Occultation year
   Month,     & ! <-- Occultation month
   Day,       & ! <-- Occultation day
   Hour,      & ! <-- Occultation hour
   Minute,    & ! <-- Occultation begin minute
   Second,    & ! <-- Occultation begin second
   TR,        & ! <-- Array of relative time of samples [sec]
   RLEO,      & ! <-- LEO coordinates (Cartesian)
   RGPS,      & ! <-- GPS coordinates (Cartesian)
   ERLEO,     & ! --> LEO coordinates in Earth frame
   ERGPS,     & ! --> GPS coordinates in Earth frame
   GP,        & ! --> Geodetic coordinates of occultation point
   ERLC,      & ! --> Curvature center in Earth frame
   RLC,       & ! --> Curvature center (Cartesian)
   RE,        & ! --> Local curvature radius
   Stat)        ! ~~> Error status
!
! Determination of occultation geometry.
!----------------------------------------------------------
! Method:
!   Coordinate frame rotation, finding occultation point,
!   finding local curvature center.
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 18 Oct 1998 | Original code.
!   1.1.  | 01 Nov 1998 | Determination of occultation
!         |             | in a separate subroutine.
!   2.0   | 17 Mar 1999 | Error status check.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, &
! Imported Routines:
    Rotate,       &
    Vector_Angle, &
! Imported Operators:
    Operator(.x.),  &
    Operator(*)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic, &
! Imported Routines:
    GAST, Curvature,  &
    Cart_from_Geod
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)          :: &
   Year      ! Occultation year
!
Integer, Intent(In)          :: &
   Month     ! Occultation month
!
Integer, Intent(In)          :: &
   Day       ! Occultation day
!
Integer, Intent(In)          :: &
   Hour      ! Occultation hour
!
Integer, Intent(In)          :: &
   Minute    ! Occultation begin minute
!
Real(Double), Intent(In)     :: &
   Second    ! Occultation begin second
!
Real(Double), Intent(In)     :: &
   TR(1:)    ! Array of relative time of samples [sec]
!
Type(Cartesian), Intent(In)  :: &
   RLEO(1:)  ! LEO coordinates (Cartesian)
!
Type(Cartesian), Intent(In)  :: &
   RGPS(1:)  ! GPS coordinates (Cartesian)
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   ERLEO(1:) ! LEO coordinates in Earth frame
!
Type(Cartesian), Intent(Out) :: &
   ERGPS(1:) ! GPS coordinates in Earth frame
!
Type(Geodetic), Intent(Out)  :: &
   GP        ! Geodetic coordinates of occultation point
!
Type(Cartesian), Intent(Out) :: &
   ERLC      ! Curvature center in Earth frame
!
Type(Cartesian), Intent(Out) :: &
   RLC       ! Curvature center (Cartesian)
!
Real(Double), Intent(Out)    :: &
   RE        ! Local curvature radius
!
Integer, Optional, Intent(Out)    :: &
   Stat      ! Error status:
             !    0 - no error
             !   -1 - array size mismatch
!----------------------------------------------------------
! Local Parameters:
!
Type(Cartesian), Parameter :: &
   PA =  &  ! Polar axis
      Cartesian((/0,0,1/))
!
! Local Scalars:
!
Integer      :: N      ! Number of data
Integer      :: i      ! Data point index
Real(Double) :: Theta  ! Cross section azimuth
Integer      :: Iocc   ! Occultation point index
Type(Cartesian) :: CP  ! Perigee position (cartesian)
!
! Local Arrays:
!
Real(Double) :: &
   Phi(Size(ERLEO))    ! GPS frame rotation
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ERROR CHECK
!----------------------------------------------------------

N = Size(TR)

If (.not. All((/         &
      Size(RLEO)  == N,  &
      Size(RGPS)  == N,  &
      Size(ERLEO) == N,  &
      Size(ERGPS) == N   &
      /))) then
   If (Present(Stat)) then
      Stat = -1
   End If
   Return
Else
   If (Present(Stat)) then
      Stat = 0
   End If
End If


!----------------------------------------------------------
! 1. FRAME ROTATION
!----------------------------------------------------------

Do i=1,N
   Phi(i) =  GAST &
     (Year,   & ! <-- Year of occultation
      Month,  & ! <-- Month
      Day,    & ! <-- Day
      Hour,   & ! <-- Hour
      Minute, & ! <-- Minute
      Second, & ! <-- Second of beginning of occultation
      TR(i))    ! <-- Second of occultation duration
   ERLEO(i) = Rotate(RLEO(i), PA, -Phi(i))
   ERGPS(i) = Rotate(RGPS(i), PA, -Phi(i))
End Do


!----------------------------------------------------------
! 2. DETERMINATION OF LOCAL CURVATURE CENTER
!----------------------------------------------------------


!--- 2.1. Determination of the occultation point

Call Occ_Point &
  (ERLEO,     & ! <-- LEO coordinates in Earth frame
   ERGPS,     & ! <-- GPS coordinates in Earth frame
   GP,        & ! --> Occultation point (Geodetic)
   Iocc)        ! --> Occultation point index


!--- 2.2. Determination of curvature center and radius

CP    = Cart_from_Geod(GP)
Theta = Vector_Angle(ERGPS(Iocc).x.ERLEO(Iocc), PA.x.CP)
Call Curvature(GP, Theta, ERLC, RE)
RLC   = Rotate(ERLC, PA, Phi(Iocc))


End Subroutine Occ_Geometry



!==========================================================
Function Perigee  &
  (X0,     & ! <-- Ray beginning point
   X1,     & ! <-- Ray ending point
   Eps)    & ! <-- Refraction angle
Result(XP)   ! --> Ray Perigee point
!
! Calculation of perigee point in approximation
! of spherical symmetry.
!----------------------------------------------------------
! Method:
!   Trigonometrical formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 24 Sep 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Rotate,        &
! Imported Operators:
    Operator(*),  Operator(.x.)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian) :: &
   X0  ! Ray beginning point
Type(Cartesian) :: &
   X1  ! Ray ending point
Real(Double)    :: &
   Eps ! Refraction angle
!
! Function result:
!
Type(Cartesian) :: &
   XP  ! Pergiee point of ray connecting X0 and X1 calculated
!      ! in approximation of spherical symmetry.
!----------------------------------------------------------
! Local Scalars:
!
Real(Double)    :: P     ! Ray impact parameter
Real(Double)    :: R0    ! Length of X0
Real(Double)    :: Alpha ! Angle X0^XP
!----------------------------------------------------------


P     = Impact_Parameter(X0, X1, Eps)
R0    = Sqrt(X0*X0)
Alpha = Eps/2 + ACos(P/R0)
XP    = Rotate(X0, X0.x.X1, Alpha)*(P/R0)


End Function Perigee



!==========================================================
Function Impact_Parameter &
  (X0,    & ! <-- Ray beginning
   X1,    & ! <-- Ray end
   Eps)   & ! <-- Refraction angle
Result(P)   ! --> Impact parameter
!
! Calculation of the impact parameter connecting two
! given points, in the assumption of the spherical
! symmetry.
!----------------------------------------------------------
! Method:
!   Trigonometric formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 24 Sep 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Routines:
    Vector_Angle,  &
! Imported Operators:
    Operator(*)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian) :: &
   X0    ! Vector of ray beginnig point
Type(Cartesian) :: &
   X1    ! Vector of ray ending point
Real(Double)    :: &
   Eps   ! Refraction angle
!
! Function result:
!
Real(Double) :: &
   P     ! Ray impact parameter calculated
!        ! in the assumption of spherical symmetry
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: R0, R1  ! Lengths of X0 and X1
Real(Double) :: Omega   ! Complementary to X0^X1 - Eps
Real(Double) :: TanAlfa ! Tan(X0^(X0-X1))
!----------------------------------------------------------


R0      = Sqrt(X0*X0)
R1      = Sqrt(X1*X1)
Omega   = Pi - Vector_Angle(X0, X1) + Eps
TanAlfa = R1*Sin(Omega)/(R0 + R1*Cos(Omega))

P       = R0*TanAlfa/Sqrt(1 + TanAlfa**2)


End Function Impact_Parameter



!==========================================================
Subroutine Occ_Point &
  (ERLEO,     & ! <-- LEO coordinates in Earth frame
   ERGPS,     & ! <-- GPS coordinates in Earth frame
   GP,        & ! --> Occultation point (Geodetic)
   Iocc)        ! --> Occultation point index
!
! Determination of the occultation point
!----------------------------------------------------------
! Method:
!   Lowest occultation perigee point projected to the
!   Earth's surface
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Nov 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic, &
! Imported Routines:
    Geod_from_Cart
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   ERLEO(1:) ! LEO coordinates in Earth frame
!
Type(Cartesian), Intent(In) :: &
   ERGPS(1:) ! GPS coordinates in Earth frame
!
! Output arguments:
!
Type(Geodetic), Intent(Out)  :: &
   GP        ! Geodetic coordinates of occultation point
!
Integer, Intent(Out)         :: &
   Iocc      ! Occultation point index
!----------------------------------------------------------
! Local Scalars:
!
Integer :: N  ! Number of data
Integer :: i  ! Data sample index
!
! Local Arrays:
!
Type(Geodetic) :: &
   GPX(Size(ERLEO))       ! Perigee positions (geodetic)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. DETERMINATION OF RAY PERIGEES
!----------------------------------------------------------

N = Size(ERLEO)

Do i=1,N
   GPX(i) = Geod_from_Cart(Perigee(ERLEO(i), ERGPS(i), 0.0_Double))
End Do


!----------------------------------------------------------
! 2. DETERMINATION OF OCCULTATION POINT
!----------------------------------------------------------


!--- 2.1. Finding the lowest perigee

Iocc = Sum(MinLoc(GPX(:)%H))
GP   = GPX(Iocc)


!--- 2.2. Projecting to the surface

GP%H = 0


End Subroutine Occ_Point



!==========================================================
Subroutine Satellite_Velocities &
  (TR,       & ! <-- Array of relative time of samples [sec]
   RLEO,     & ! <-- LEO coordinates (cartesian)
   RGPS,     & ! <-- GPS coordinates (cartesian)
   XLEO,     & ! --> LEO coordinates from regression (cartesian)
   VLEO,     & ! --> LEO velocities from regression (cartesian)
   XGPS,     & ! --> GPS coordinates from regression (cartesian)
   VGPS,     & ! --> GPS velocities from regression (cartesian)
   Stat)       ! ~~> Error status
!
! Calculation of satellite velocities
!----------------------------------------------------------
! Method:
!   Polynomial regression
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Nov 1998 | Original code.
!   2.0   | 17 Mar 1999 | Error status check.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian,     &
! Imported Operators:
    Operator(*),   &
    Operator(-),   &
    Operator(/)
!
Use Matrix, only: &
! Imported Routines:
    Regression,         &
    Basic_Polynomials
!
Use Interpolation, only: &
! Imported Routines:
    Polynomial
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)     :: &
   TR(1:)     ! Array of relative time of samples [sec]
!
Type(Cartesian), Intent(In)  :: &
   RLEO(1:)   ! LEO coordinates
!
Type(Cartesian), Intent(In)  :: &
   RGPS(1:)   ! GPS coordinates
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   XLEO(1:)   ! LEO coordinates from regression
!
Type(Cartesian), Intent(Out) :: &
   VLEO(1:)   ! LEO velocities from regression
!
Type(Cartesian), Intent(Out) :: &
   XGPS(1:)   ! GPS coordinates from regression
!
Type(Cartesian), Intent(Out) :: &
   VGPS(1:)   ! GPS velocities from regression
!
Integer, Optional, Intent(Out) :: &
   Stat       ! Error status:
              !    0 - no error
              !   -1 - array size mismatch
!----------------------------------------------------------
! Local Parameters:
!
Integer, Parameter      :: &
   NV = 5       ! Polynomial degree for calculation of velocity
!
! Local Scalars:
!
Integer         :: N       ! Number of data
Integer         :: i       ! Sample number
Integer         :: m       ! Dimension index
Integer         :: ErrCode ! Error code
!
! Local Arrays:
!
Real(Double)    :: &
   TRV(Size(TR))        ! Normed TR
Real(Double)    :: &
   KV(Size(TR),0:NV)    ! Regression matrix for velocities
Real(Double)    :: &
   BG(0:NV,3),   &      ! Regression coefficients for VGPS
   BL(0:NV,3)           ! Regression coefficients for VLEO
!----------------------------------------------------------


!----------------------------------------------------------
! 0. ERROR CHECK
!----------------------------------------------------------

N = Size(TR)

If (.not. All((/         &
      Size(RLEO) == N,   &
      Size(RGPS) == N,   &
      Size(XLEO) == N,   &
      Size(XGPS) == N,   &
      Size(VLEO) == N,   &
      Size(VGPS) == N    &
      /))) then
   If (Present(Stat)) then
      Stat = -1
   End If
   Return
Else
   If (Present(Stat)) then
      Stat = 0
   End If
End If


!----------------------------------------------------------
! 1. RENORMING TIME
!----------------------------------------------------------

TRV(:) = (TR(:) - TR(1))/(TR(N) - TR(1))


!----------------------------------------------------------
! 2. CALCULATION OF REGRESSION COEFFICIENTS
!----------------------------------------------------------

Call Basic_Polynomials(TRV, KV, ErrCode)
Do m=1,3
   Call Regression(KV, RGPS(:)%X(m), BG(:,m), ErrCode)
   Call Regression(KV, RLEO(:)%X(m), BL(:,m), ErrCode)
End Do


!----------------------------------------------------------
! 3. REGRESSION
!----------------------------------------------------------

Do i=1,N
   Do m=1,3
      Call Polynomial(BG(:,m), TRV(i), XGPS(i)%X(m), VGPS(i)%X(m))
      Call Polynomial(BL(:,m), TRV(i), XLEO(i)%X(m), VLEO(i)%X(m))
   End Do
   VGPS(i) = VGPS(i)/(TR(N) - TR(1))
   VLEO(i) = VLEO(i)/(TR(N) - TR(1))
End Do


End Subroutine Satellite_Velocities



!==========================================================
Subroutine Plane_Coordinates &
  (RLEO,      & ! <-- LEO coordinates (Cartesian)
   RGPS,      & ! <-- GPS coordinates (Cartesian)
   RLC,       & ! <-- Curvature center (Cartesian)
   RE,        & ! <-- Local curvature radius [km]
   XLEO,      & ! --> X coordinates of LEO
   YLEO,      & ! --> Y coordinates of LEO
   XGPS,      & ! --> X coordinates of GPS
   YGPS,      & ! --> Y coordinates of GPS
   AX,        & ! --> Occultation plane X basis vector
   AY)          ! --> Occultation plane Y basis vector
!
! Calculation of occultation plane coordinates of GPS and LEO.
!----------------------------------------------------------
! Method:
!   Calculation of coordinates in GPS-O-LEO plane,
!   transform  to unmoving GPS. (See: Report No. 211,
!   MPIM, Hamburg, 1996).
!----------------------------------------------------------
! (C) Copyright 1998-99, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 31 Oct 1998 | Original code.
!   1.1   | 31 Oct 1998 | Both rising and settings
!         |             | acceptable.
!   1.2.  | 01 Nov 1998 | Use of Occ_Point.
!   2.0   | 09 Nov 1998 | Use of local curvature center.
!   2.1   | 25 Feb 1999 | Excessive variables excluded.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, &
! Imported Routines:
    Vector_Normed,   &
    Vector_Norm,     &
! Imported Operators:
    Operator(-),   &
    Operator(+),   &
    Operator(*),   &
    Operator(/),   &
    Assignment(=)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   RLEO(1:)   ! LEO coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RGPS(1:)   ! GPS coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RLC        ! Cartesian coordinates of curvature center
!
Real(Double), Intent(In)   :: &
   RE         ! Local curvature radius [km]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   XLEO(1:)   ! X coordinates of LEO
!
Real(Double), Intent(Out) :: &
   YLEO(1:)   ! Y coordinates of LEO
!
Real(Double), Intent(Out) :: &
   XGPS(1:)   ! X coordinates of GPS
!
Real(Double), Intent(Out) :: &
   YGPS(1:)   ! Y coordinates of GPS
!
Type(Cartesian), Intent(Out) :: &
   AX         ! Occultation plane X basis vector
!
Type(Cartesian), Intent(Out) :: &
   AY         ! Occultation plane Y basis vector
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: N       ! Number of data
Type(Geodetic)  :: GP      ! Occultation point
Integer         :: i       ! Array index
Integer         :: Iocc    ! Occultation point index
Type(Cartesian) :: NG, NL  ! Orthonormalized (RGPS, RLEO) basis
Real(Double)    :: CT, ST  ! Cos(RLEO^AY) and Sin(RLEO^AY)
Real(Double)    :: F       ! Coordinate invariant
!
! Local Array:
!
Type(Cartesian) :: &
   AXS(Size(RLEO)), & ! Current X basis vector
   AYS(Size(RLEO))    ! Current Y basis vector
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF OCCULTATION POINT
!----------------------------------------------------------

N = Size(RLEO)

Call Occ_Point &
  (RLEO,     & ! <-- LEO coordinates
   RGPS,     & ! <-- GPS coordinates
   GP,       & ! --> Occultation point (Geodetic)
   Iocc)       ! --> Occultation point index


!----------------------------------------------------------
! 2. TRANSFORM TO COORDINATES IN OCCULTATION PLANE
!----------------------------------------------------------

Do i=1,N
   NG      = Vector_Normed(RGPS(i)-RLC)
   NL      = Vector_Normed(RLEO(i)-RLC - ((RLEO(i)-RLC)*NG)*NG)
   CT      = RE/Vector_Norm(RGPS(i)-RLC)
   ST      = Sqrt(1 - CT**2)
   AXS(i)  = CT*NL - ST*NG
   AYS(i)  = ST*NL + CT*NG
   XLEO(i) = AXS(i)*(RLEO(i)-RLC)
   YLEO(i) = AYS(i)*(RLEO(i)-RLC)
   XGPS(i) = AXS(i)*(RGPS(i)-RLC)
   YGPS(i) = AYS(i)*(RGPS(i)-RLC)
End Do

AX = AXS(Iocc)
AY = AYS(Iocc)


!----------------------------------------------------------
! 3. TRANSFORM TO UNMOVING GPS
!----------------------------------------------------------

Do i=1,N
   F       = XLEO(i)*XGPS(i)/(XLEO(i)-XGPS(i))
   XGPS(i) = XGPS(Iocc)
   XLEO(i) = F*XGPS(i)/(F-XGPS(i))
End Do


End Subroutine Plane_Coordinates



!==========================================================
Subroutine Plane_Basis &
  (RLEO,      & ! <-- LEO coordinates (Cartesian)
   RGPS,      & ! <-- GPS coordinates (Cartesian)
   RLC,       & ! <-- Curvature center (Cartesian)
   RE,        & ! <-- Local curvature radius [km]
   AX,        & ! --> Occultation plane X basis vector
   AY)          ! --> Occultation plane Y basis vector
!
! Calculation of basis vectors of occultation plane.
!----------------------------------------------------------
! Method:
!   See: Report No. 211, MPIM, Hamburg, 1996.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Feb 1999 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, &
! Imported Routines:
    Vector_Normed,   &
    Vector_Norm,     &
! Imported Operators:
    Operator(-),   &
    Operator(+),   &
    Operator(*),   &
    Operator(/),   &
    Assignment(=)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   RLEO(1:)   ! LEO coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RGPS(1:)   ! GPS coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RLC        ! Cartesian coordinates of curvature center
!
Real(Double), Intent(In)   :: &
   RE         ! Local curvature radius [km]
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   AX         ! Occultation plane X basis vector
!
Type(Cartesian), Intent(Out) :: &
   AY         ! Occultation plane Y basis vector
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: Iocc    ! Occultation point index
Type(Geodetic)  :: GP      ! Occultation point
Type(Cartesian) :: NG, NL  ! Orthonormalized (RGPS, RLEO) basis
Real(Double)    :: CT, ST  ! Cos(RLEO^AY) and Sin(RLEO^AY)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF OCCULTATION POINT
!----------------------------------------------------------

Call Occ_Point &
  (RLEO,     & ! <-- LEO coordinates
   RGPS,     & ! <-- GPS coordinates
   GP,       & ! --> Occultation point (Geodetic)
   Iocc)       ! --> Occultation point index


!----------------------------------------------------------
! 2. CACLULATION OF BASIS VECTORS
!----------------------------------------------------------

NG = Vector_Normed(RGPS(Iocc)-RLC)
NL = Vector_Normed(RLEO(Iocc)-RLC - ((RLEO(Iocc)-RLC)*NG)*NG)
CT = RE/Vector_Norm(RGPS(Iocc)-RLC)
ST = Sqrt(1 - CT**2)
AX = CT*NL - ST*NG
AY = ST*NL + CT*NG


End Subroutine Plane_Basis



!==========================================================
Function BPP_Position &
  (RLEO,      & ! <-- LEO coordinates (Cartesian)
   RGPS,      & ! <-- GPS coordinates (Cartesian)
   BRLEO,     & ! <-- BP LEO coordinates (Cartesian)
   BRGPS,     & ! <-- BP GPS coordinates (Cartesian)
   RLC,       & ! <-- Curvature center (Cartesian)
   RE)        & ! <-- Local curvature radius [km]
Result(XB)      ! --> BP plane X-coordinate
!
! Calculation of position of BP plane.
!----------------------------------------------------------
! Method:
!   Calculation of occultation plane coordinates.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Nov 1998 | Original code.
!   2.0   | 09 Nov 1998 | Use of local curvature center.
!   2.1   | 10 Dec 1999 | Allocatable arrays, N and NB.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, &
! Imported Operators:
    Operator(*),  &
    Operator(-)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   RLEO(1:)   ! LEO coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RGPS(1:)   ! GPS coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   BRLEO(1:)  ! BP LEO coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   BRGPS(1:)  ! BP GPS coordinates (Cartesian)
!
Type(Cartesian), Intent(In) :: &
   RLC        ! Cartesian coordinates of curvature center
!
Real(Double), Intent(In)   :: &
   RE         ! Local curvature radius [km]
!
! Function result:
!
Real(Double) :: &
  XB          ! BP plane X coordinate
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: N    ! Number of GPS data
Integer         :: NB   ! Number of BP data
Integer         :: i    ! Array index
Type(Cartesian) :: &
   AX,       & ! Occultation plane X basis vector
   AY          ! Occultation plane Y basis vector
!
! Local Arrays:
!
Real(Double), Allocatable    :: &
   XLEO(:),    & ! X coordinates of LEO
   YLEO(:),    & ! Y coordinates of LEO
   XGPS(:),    & ! X coordinates of GPS
   YGPS(:)       ! Y coordinates of GPS
Real(Double), Allocatable    :: &
   BXLEO(:),   & ! BP X coordinates of LEO
   BYLEO(:),   & ! BP Y coordinates of LEO
   BXGPS(:),   & ! BP X coordinates of GPS
   BYGPS(:)      ! BP Y coordinates of GPS
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF OCCULTATION PLANE COORDINATES
!----------------------------------------------------------


!--- 1.1. Array allocation

N = Size(RLEO)

Allocate(XLEO(1:N), YLEO(1:N), XGPS(1:N), YGPS(1:N))


!--- 1.2. Coordinate calculation

Call Plane_Coordinates &
  (RLEO,     & ! <-- LEO coordinates (Cartesian)
   RGPS,     & ! <-- GPS coordinates (Cartesian)
   RLC,      & ! <-- Curvature center (Cartesian)
   RE,       & ! <-- Local curvature radius [km]
   XLEO,     & ! --> X coordinates of LEO
   YLEO,     & ! --> Y coordinates of LEO
   XGPS,     & ! --> X coordinates of GPS
   YGPS,     & ! --> Y coordinates of GPS
   AX,       & ! --> Occultation plane X basis vector
   AY)         ! --> Occultation plane Y basis vector


!--- 1.3. Array deallocation

Deallocate(XLEO, YLEO, XGPS, YGPS)


!----------------------------------------------------------
! 2. CALCULATION OF XY COORDINATES OF BP GPS/LEO
!----------------------------------------------------------


!--- 2.1. Array allocation

NB = Size(BRLEO)

Allocate(BXLEO(1:NB), BYLEO(1:NB), BXGPS(1:NB), BYGPS(1:NB))


!--- 2.2. Coordinate calculation

Do i=1,NB
   BXLEO(i) = AX*(BRLEO(i)-RLC)
   BYLEO(i) = AY*(BRLEO(i)-RLC)
   BXGPS(i) = AX*(BRGPS(i)-RLC)
   BYGPS(i) = AY*(BRGPS(i)-RLC)
End Do


!----------------------------------------------------------
! 3. CALCULATION OF BP PLANE POSITION
!----------------------------------------------------------


!--- 3.1. Calculation of BPP position

XB = Sum(BXLEO(:))/NB


!--- 3.2. Array deallocation

Deallocate(BXLEO, BYLEO, BXGPS, BYGPS)


End Function BPP_Position



!==========================================================
Subroutine Find_Positions &
  (PN,        & ! <-- Impact parameter
   EN,        & ! <-- Refraction angle
   TR,        & ! <-- Array of relative time of samples [sec]
   ERLEO,     & ! <-- LEO coordinates in Earth frame
   ERGPS,     & ! <-- GPS coordinates in Earth frame
   RLC,       & ! <-- Cartesian coordinates of curvature center
   DRLEO,     & ! --> LEO coordinates for EN(PN)
   DVLEO,     & ! --> LEO velocities for EN(PN)
   DRGPS,     & ! --> GPS coordinates for EN(PN)
   DVGPS)       ! --> GPS velocities for EN(PN)
!
! Finding satellite positions and velocities corresponding
! to refraction angle profile.
!----------------------------------------------------------
! Method:
!   Solution of equation:
!   XGPS^XLEO - acos(P/|XGPS|) - acos(P/|XLEO|) = Eps
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Nov 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, &
! Imported Routines:
    Vector_Angle,    &
    Vector_Norm,     &
! Imported Operators:
    Operator(-),   &
    Operator(+),   &
    Operator(*),   &
    Operator(/),   &
    Assignment(=)
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   PN(1:)      ! Impact parameter
!
Real(Double), Intent(In) :: &
   EN(1:)      ! Refraction angle
!
Real(Double), Intent(In)     :: &
   TR(1:)      ! Array of relative time of samples [sec]
!
Type(Cartesian), Intent(In) :: &
   ERLEO(1:)   ! LEO coordinates in Earth frame
!
Type(Cartesian), Intent(In) :: &
   ERGPS(1:)   ! GPS coordinates in Earth frame
!
Type(Cartesian), Intent(In) :: &
   RLC         ! Cartesian coordinates of curvature center
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   DRLEO(1:)   ! LEO coordinates for EN(PN)
!
Type(Cartesian), Intent(Out) :: &
   DVLEO(1:)   ! LEO velocities for EN(PN)
!
Type(Cartesian), Intent(Out) :: &
   DRGPS(1:)   ! GPS coordinates for EN(PN)
!
Type(Cartesian), Intent(Out) :: &
   DVGPS(1:)   ! GPS velocities for EN(PN)
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: Imin  ! Index of maximum p
Integer      :: Imax  ! Index of minimum p
Integer      :: N     ! Number of GPS/MET data
Integer      :: NB    ! Number of refraction angles
Integer      :: i     ! Sample index
Integer      :: Jmin  ! Lowper index limit for dichotomy
Integer      :: Jmax  ! Upper index limit for dichotomy
Integer      :: J     ! Dichotomically found index
Real(Double) :: EGL   ! Estimation of refraction angle
!
! Local Arrays:
!
Type(Cartesian)  :: &
   XGPS(Size(ERLEO)),  & ! GPS positions from regression
   XLEO(Size(ERLEO)),  & ! LEO positions from regression
   VGPS(Size(ERLEO)),  & ! GPS velocities from regression
   VLEO(Size(ERLEO))     ! LEO velocities from regression
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF VELOCITIES
!----------------------------------------------------------

Call Satellite_Velocities &
  (TR,       & ! <-- Array of relative time of samples [sec]
   ERLEO,    & ! <-- LEO coordinates (cartesian)
   ERGPS,    & ! <-- GPS coordinates (cartesian)
   XLEO,     & ! --> LEO coordinates from regression (cartesian)
   VLEO,     & ! --> LEO velocities from regression (cartesian)
   XGPS,     & ! --> GPS coordinates from regression (cartesian)
   VGPS)       ! --> GPS velocities from regression (cartesian)


!----------------------------------------------------------
! 2. FINDING SATELLITE POSITIONS FOR EPS(P)
!----------------------------------------------------------

NB = Size(EN)
N  = Size(XLEO)

If (PN(NB) > PN(1)) then
   Imin = 1
   Imax = N
Else
   Imin = N
   Imax = 1
End If

Do i=1,NB
   Jmin = Imax
   Jmax = Imin
   Find: Do
      J   = (JMin + JMax)/2
      EGL = Vector_Angle (XLEO(J)-RLC, XGPS(J)-RLC) - &
                ACos (PN(i)/Vector_Norm(XLEO(J)-RLC))  - &
                ACos (PN(i)/Vector_Norm(XGPS(J)-RLC))
      If (EGL > EN(i)) then
         Jmax = J
      Else
         Jmin = J
      End if
      If (Abs(Jmax-Jmin) == 1) then
         Exit Find
      End If
   End Do Find
   DRLEO(i) = XLEO(J)
   DRGPS(i) = XGPS(J)
   DVLEO(i) = VLEO(J)
   DVGPS(i) = VGPS(J)
End Do


End Subroutine Find_Positions



End Module Occ_Coordinates

