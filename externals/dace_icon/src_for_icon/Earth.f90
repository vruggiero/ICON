!
!+ GNSS Radio occultation operator: Earth's parameters, transformations
!
MODULE Earth
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Earth's parameters, transformations of geodetic coordinates,
!   gravity field, and GPS/MET frame rotation.
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
! V1_6         2009/06/10 Harald Anlauf
!  Fix precision of FP constants
! V1_9         2010/04/20 Harald Anlauf
!  Change value of C2 (0.37) to the one given in the literature (0.373)
!  Gravity, Geop_from_Alt,Alt_from_Geop,Alt_from_Geop_adj: save on sine calls
!  Fix digit transposition (Zahlendreher) in Earth's flatness parameter
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  Revise coefficients for the 3-term refractivity formula (C1,C2,C3)
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_28        2014/02/26 Harald Anlauf
!  Dry air gas constant: GME /= ICON, use #ifdef __ICON__
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
! Module Earth
!
! Earth's parameters, transformations of geodetic coordinates,
! gravity field, and GPS/MET frame rotation.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Sep 1998 | Original code.
!   2.0   | 18 Sep 1998 | Geodetic coordinates,
!         |             | gravity field and rotation.
!   3.0   | 10 Jun 1999 | Consecutive use of Earth constants.
!   4.0   | 06 Jul 1999 | Rd, Rnu, C1, and C2.
!   4.1   | 20 Apr 2000 | Jacobian_GC_adj, Hessian_GC.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian, Spherical
!----------------------------------------------------------
Implicit None
Private
Public :: Geodetic, Rd, C1, Gclat_from_Gdlat, Gravity, Geop_from_Alt, &
          Gast, Curvature, Cart_from_Geod, Geod_from_Cart, R_Earth,   &
          Rnu, C2, C3, H_Atm, G_Ave, Alt_from_Geop, Hessian_GC,       &
          Jacobian_GC, Alt_from_Geop_adj
!----------------------------------------------------------
! Public Type Definitions:
!
Type Geodetic                ! Geodetic coordinates:
   Real(Double) :: H         ! Height above reference ellipsoid [km]
   Real(Double) :: Phi       ! Latitude from equator [degree]
   Real(Double) :: Lambda    ! Longitude [degree]
End Type
!
! Public Parameters:
!
Real(Double), Parameter :: &
   R_Earth = 6370.0_Double,  & ! Average Earth radius [km]
   H_atm   = 130.0_Double,   & ! Conventional height of the atmosphere [km]
   g_ave   = 9.80665_Double    ! Average gravity acceleration [m/s^2]
!
Real(Double), Parameter :: &
#ifdef __ICON__
   Rd      = 287.04_Double,  & ! Dry air gas constant [J/kg*K] (ICON)
#else
   Rd      = 287.05_Double,  & ! Dry air gas constant [J/kg*K] (GME)
#endif
   Rnu     = 461.51_Double     ! Water vapour gas constant [J/kg*K]
!
#ifdef _SMITH_WEINTRAUB_TWOTERM
! Coefficients for the classic Smith & Weintraub 2-term refractivity formula:
Real(Double), Parameter :: &
   C1      = 77.6e-6_Double,  & !  N = C1*P/T + C2*Pw/T^2
   C2      = 0.373_Double,    &
   C3      = 0
#else
! Revised coefficients for the 3-term formula from:
! Huw Lewis, Refractivity calculations in ROPP, GRAS SAF Report 05 (2008).
Real(Double), Parameter :: &
   C1      = 77.6890e-6_Double, & !  N = C1*P/T + C2*Pw/T^2 + C3*Pw/T
   C2      =  0.375463_Double,  &
   C3      = -6.3938e-6_Double
#endif
!
!
! Private parameters:
!
! --- Earth's shape
!
Real(Double), Private, Parameter :: &
   f      = 1/298.257_Double,       & ! Flatness
   Re     = 6378.136_Double,        & ! Equatorial semiaxis [km]
   Rp     = Re*(1-f),               & ! Polar semiaxis [km]
   flatfn = (2 - f)*f,              &
   funsq  = (1 - f)**2,             &
   e      = (Re**2 - Rp**2)/(Re**2)   ! Squared eccentricity
!
! --- Earth's gravity
!
Real(Double), Private, Parameter :: &
   GM = 3.9860044e14_Double,        & ! Earth's gravity constant [m^3/s^2]
   m  = 0.00345_Double,             & ! Centrifugal force/gravity at equator
   ge = GM/(1e6_Double*Re**2*       &
       (1-f+3*m/2-15*m*f/14)),      & ! Gravity at equator [m/s^2]
   f2 = -f + 5*m/2 -      &
        17*f*m/14 +15*m**2/4,       &
   f4 = -f**2/2 + 5*f*m/2
!----------------------------------------------------------
Contains


!==========================================================
Function Cart_from_Geod &
  (G)     &  ! <-- Geodetic coordinates
Result(X)    ! --> Cartesion coordinates
!
! Conversion from geodetic to Cartesian coordinates
!----------------------------------------------------------
! Method:
!   Escobal, Methods of orbit determination
!   1965, Wiley & sons, Inc. Pp. 27 - 29.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    dtr
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In)  :: &
   G  ! The geodetic coordinates of a point
!
! Function result:
!
Type(Cartesian)            :: &
   X  ! The Cartesian coordinates of the point
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: gd  ! Parallel curvature of ellipsoid
!----------------------------------------------------------


gd = Re/Sqrt(1.0_Double - flatfn*Sin(G%Phi*dtr)**2)

X%X(2) = Cos(G%Phi*dtr)*Sin(G%Lambda*dtr)*(gd + G%H)
X%X(1) = Cos(G%Phi*dtr)*Cos(G%Lambda*dtr)*(gd + G%H)
X%X(3) = Sin(G%Phi*dtr)*(gd*funsq + G%H)


End Function Cart_from_Geod



!==========================================================
Function Geod_from_Cart &
  (X)    & ! <-- Cartesion coordinates
Result(G)  ! --> Geodetic coordinates
!
! Conversion from Cartesian to geodetic coordinates
!----------------------------------------------------------
! Method:
!   Fast approximate solution.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 18 Sep 1998 | Function.
!   3.0   | 11 Jun 1999 | Fast approximation.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    rtd, dtr
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Cartesian), Intent(In)  :: &
   X  ! Cartesian coordinates of a point
!
! Function result:
!
Type (Geodetic)              :: &
   G  ! Geodetic coordinates of the point
!----------------------------------------------------------
! Local scalars:
!
Real(Double) :: Rxy
Real(Double) :: Theta, CosPhi, SinPhi
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

Rxy   = Sqrt(Sum(X%X(1:2)**2))
Theta = Atan2(X%X(3)*Re, Rxy*Rp)


!----------------------------------------------------------
! 2. FAST APPROXIMATE SOLUTION
!----------------------------------------------------------


!--- 2.1. Approximate calculation of longitude

G%Phi  = rtd*Atan2(X%X(3) + e*Rp*Sin(Theta)**3,   &
                   Rxy    - e*Re*Cos(Theta)**3)
CosPhi = Cos(dtr*G%Phi)
SinPhi = Sin(dtr*G%Phi)


!--- 2.2. Calculation of height from longitude

G%H = -(-CosPhi*Rp**2*Rxy - Re**2*SinPhi*X%X(3) +  &
        Re*Rp*Sqrt(CosPhi**2*Rp**2 +               &
           Re**2*SinPhi**2 - SinPhi**2*Rxy**2 +    &
           2*CosPhi*SinPhi*Rxy*X%X(3) -            &
           CosPhi**2*X%X(3)**2))/                  &
       (CosPhi**2*Rp**2 + Re**2*SinPhi**2)


!--- 2.3. Calculation of longitude

G%Lambda = rtd*Atan2(X%X(2),X%X(1))


End Function Geod_from_Cart



!==========================================================
Function GCLat_from_GDLat &
  (GDLat)    & ! <-- Geodetic latitude (degree)
Result(GCLat)  ! --> Geocentric latitude (radian)
!
! Conversion of a geodetic latitude to a geocedntric one.
!----------------------------------------------------------
! Method:
!   Geodetic -> Cartesian -> Spherical coordinate conversion
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Pi
!
Use Coordinates, only: &
! Imported Routines:
    Spher_from_Cart
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   GDLat  ! Geodetic latitude (degree)
!
! Function result:
!
Real(Double)             :: &
   GCLat  ! Geocentric latitude from Equator [radian]
!----------------------------------------------------------
! Local Variables:
!
Type (Geodetic)  :: G  ! Geodetic coordinates
Type (Cartesian) :: X  ! Cartesian coordinates
Type (Spherical) :: S  ! Spherical coordinates
!----------------------------------------------------------


G%H      = 0.0_Double
G%Phi    = GDLat
G%Lambda = 0.0_Double

X = Cart_from_Geod (G)
S = Spher_from_Cart(X)

GCLat = Pi/2 - S%Phi


End Function GCLat_from_GDLat



!==========================================================
Subroutine Curvature &
  (G,     & ! <-- Geodetic coordinates of a surface point
   Theta, & ! <-- Azimuth direction of the cross-section
   XLC,   & ! --> Local curvature center
   RC)      ! --> Curvature radius
!
! Calculation of the section curvature and its center
! for the Earth reference ellipsoid.
!----------------------------------------------------------
! Method:
!   The use of the theorems of Meusnier and Euler.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Oct 1997 | Original code.
!   1.1   | 20 Oct 1997 | The use of type Cartesian.
!   1.2   | 20 Sep 1998 | Projection of G to surface.
!----------------------------------------------------------
!
Use Defaults, only: &
! Imported Parameters:
    dtr
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type (Geodetic), Intent(In)    :: &
   G     ! Geodetic coordinates of a surface point
Real(Double), Intent(In)       :: &
   Theta ! Azimuts angle of cross-section [radian]
!
! Output arguments:
!
Type (Cartesian), Intent(Out)  :: &
   XLC   ! Cartesian coordinates of the local curvature center
Real(Double), Intent(Out)      :: &
   RC    ! Curvature radius
!----------------------------------------------------------
! Local Varibles:
!
Type (Cartesian) :: X     ! Cartesian coordinates
Real(Double)     :: Ak1   ! Main meridional curvature
Real(Double)     :: Ak2   ! Main parallel curvature
Type (Cartesian) :: N     ! Surface normal
!----------------------------------------------------------


!----------------------------------------------------------
! 1. PROJECTION OF POINT TO THE SURFACE
!----------------------------------------------------------

X   = Cart_from_Geod(Geodetic(0.0_Double, G%Phi, G%Lambda))


!----------------------------------------------------------
! 2. CALCULATION OF CURVATURE RADIUS
!----------------------------------------------------------

Ak1 = Rp*Re**4/(Sqrt(Re**4 + (Rp**2 - Re**2)*Sum(X%X(1:2)**2)))**3
Ak2 = Cos(G%Phi*dtr)/Sqrt(Sum(X%X(1:2)**2))
Rc  = 1.0_Double/Abs(Ak1*(Cos(Theta))**2 + Ak2*(Sin(Theta))**2)


!----------------------------------------------------------
! 3. CALCULATION OF CURVATURE CENTER
!----------------------------------------------------------


!--- 3.1. Calculation of normal to the surface

N%X(1) = Cos(dtr*G%Phi)*Cos(dtr*G%Lambda)
N%X(2) = Cos(dtr*G%Phi)*Sin(dtr*G%Lambda)
N%X(3) = Sin(dtr*G%Phi)


!--- 3.2. Calculation of curvature center

XLC%X(:) = X%X(:) - RC*N%X(:)


End Subroutine Curvature



!==========================================================
Function Jacobian_GC &
  (G)         & ! <-- Geodetic coordinates of a point
Result(JGC)     ! --> Jacobian d(Geodetic)/d(Cartesian)
!
! Calculation of Jacobian matrix of Geodetic coordinates
! with respect to Cartesian coordinates:
! d(z,Phi,Lambda)/d(x1,x2,x3).
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Nov 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    dtr, rtd
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G         ! Geodetic coordinates of a point
             ! G must not be polar.
!
! Function result:
!
Real(Double) :: &
   JGC(3,3)  ! Jacobian d(z,Phi,Lambda)/d(x1,x2,x3)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double)    :: CP    ! Cos(Phi)
Real(Double)    :: SP    ! Sin(Phi)
Real(Double)    :: CL    ! Cos(Lambda)
Real(Double)    :: SL    ! Sin(Lambda)
Type(Cartesian) :: VZ    ! Unit vector in Z-direction
Type(Cartesian) :: VPhi  ! Unit vector in Phi-direction
Type(Cartesian) :: VLam  ! Unit vector in Lambda-direction
Real(Double)    :: RPhi  ! Horizontal parllel curvature
Real(Double)    :: RLam  ! Meridional curvature
Real(Double)    :: gd    ! Basic parallel radius
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COORDINATE TANGENT VECTORS AND CURVATURES
!----------------------------------------------------------


!--- 1.0. Sines and cosines

CP = Cos(dtr*G%Phi)
SP = Sin(dtr*G%Phi)
CL = Cos(dtr*G%Lambda)
SL = Sin(dtr*G%Lambda)


!--- 1.1. Z-direction

VZ%X   = (/  CP*CL, CP*SL, SP  /)


!--- 1.2. Lambda-direction

VLam%X = (/ -SL,  CL,  0.0_Double  /)

gd     = Re/Sqrt(1.0_Double - flatfn*SP**2)
RLam   = CP*(gd + G%H)


!--- 1.3. Phi-direction

VPhi%X = (/ VZ%X(2)*VLam%X(3) - VZ%X(3)*VLam%X(2),  &
            VZ%X(3)*VLam%X(1) - VZ%X(1)*VLam%X(3),  &
            VZ%X(1)*VLam%X(2) - VZ%X(2)*VLam%X(1)   /)

RPhi   = G%H + (Sqrt(Re**4 + (Rp**2 - Re**2)*(CP*gd)**2))**3/(Rp*Re**4)


!----------------------------------------------------------
! 2. DERIVATIVES
!----------------------------------------------------------

JGC(1,:) = VZ%X(:)
JGC(2,:) = rtd*VPhi%X(:)/RPhi
JGC(3,:) = rtd*VLam%X(:)/RLam


End Function Jacobian_GC



!==========================================================
Subroutine Jacobian_GC_adj &
  (G,     & ! <-- Geodetic coordinates of a point
   JGC,   & ! --> Jacobian d(Geodetic(i))/d(Cartesian(j))
   JGC_G)   ! --> d(JGC(i,j))/d(Geodetic(k))
!
! Calculation of Jacobian matrix of Geodetic coordinates
! with respect to Cartesian coordinates
! JGC(i,j) = d(z,Phi,Lambda)/d(x(j)) and "mixed" hessian
! d(JGC(i,j))/d(z,Phi,Lambda).
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Apr 2000 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    dtr, rtd
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G             ! Geodetic coordinates of a point
                 ! G must not be polar.
!
! Output arguments:
!
Real(Double), Intent(Out)  :: &
   JGC(3,3)      ! Jacobian d(z,Phi,Lambda)/d(x1,x2,x3)
Real(Double), Intent(Out)  :: &
   JGC_G(3,3,3)  ! Mixed hessian
!----------------------------------------------------------
! Local Scalars:
!
Real(Double)    :: CP    ! Cos(Phi)
Real(Double)    :: SP    ! Sin(Phi)
Real(Double)    :: CL    ! Cos(Lambda)
Real(Double)    :: SL    ! Sin(Lambda)
Type(Cartesian) :: VZ    ! Unit vector in Z-direction
Type(Cartesian) :: VPhi  ! Unit vector in Phi-direction
Type(Cartesian) :: VLam  ! Unit vector in Lambda-direction
Real(Double)    :: RPhi  ! Horizontal parllel curvature
Real(Double)    :: RLam  ! Meridional curvature
Real(Double)    :: gd    ! Basic parallel radius
Integer         :: i     ! Dimension index
!
! Local Arrays:
!
Real(Double)    :: &
   CP_G(3),       & ! d(Cos(Phi))/d(Geodetic(i))
   SP_G(3),       & ! d(Sin(Phi))/d(Geodetic(i))
   CL_G(3),       & ! d(Cos(Lambda))/d(Geodetic(i))
   SL_G(3)          ! d(Sin(Lambda))/d(Geodetic(i))
Type(Cartesian) :: &
   VZ_G(3)          ! d(VZ)/d(Geodetic(i))
Type(Cartesian) :: &
   VLam_G(3)        ! d(VLam)/d(Geodetic(i))
Real(Double)    :: &
   gd_G(3),       & ! d(gd)/d(Geodetic(i))
   RLam_G(3)        ! d(RLam)/d(Geodetic(i))
Type(Cartesian) :: &
   VPhi_G(3)        ! d(VPhi)/d(Geodetic(i))
Real(Double)    :: &
   RPhi_G(3)        ! d(RPhi)/d(Geodetic(i))
!----------------------------------------------------------


!----------------------------------------------------------
! 1. COORDINATE TANGENT VECTORS AND CURVATURES
!----------------------------------------------------------


!--- 1.0. Sines and cosines

CP = Cos(dtr*G%Phi)
SP = Sin(dtr*G%Phi)
CL = Cos(dtr*G%Lambda)
SL = Sin(dtr*G%Lambda)

!------ 1.0.a. Adjoint calculations

CP_G(:) = 0
SP_G(:) = 0
CL_G(:) = 0
SL_G(:) = 0

CP_G(2) = -dtr*SP
SP_G(2) =  dtr*CP
CL_G(3) = -dtr*SL
SL_G(3) =  dtr*CL


!--- 1.1. Z-direction

VZ%X   = (/  CP*CL, CP*SL, SP  /)

!------ 1.1.a. Adjoint calculations

Do i=1,3
   VZ_G(i)%X   = (/  CP_G(i)*CL + CP*CL_G(i),  &
                     CP_G(i)*SL + CP*SL_G(i),  &
                     SP_G(i)  /)
End Do


!--- 1.2. Lambda-direction

VLam%X = (/ -SL,  CL,  0.0_Double  /)

gd     = Re/Sqrt(1.0_Double - flatfn*SP**2)
RLam   = CP*(gd + G%H)

!------ 1.2.a. Adjoint calculations

Do i=1,3
   VLam_G(i)%X = (/ -SL_G(i),  CL_G(i),  0.0_Double  /)
End Do

gd_G(:) = SP_G(:)*Re*flatfn*SP/              &
          (Sqrt(1.0_Double - flatfn*SP**2)**3)
RLam_G(:) = CP_G(:)*(gd + G%H) + CP*gd_G(:)
RLam_G(1) = RLam_G(1) + CP


!--- 1.3. Phi-direction

VPhi%X = (/ VZ%X(2)*VLam%X(3) - VZ%X(3)*VLam%X(2),  &
            VZ%X(3)*VLam%X(1) - VZ%X(1)*VLam%X(3),  &
            VZ%X(1)*VLam%X(2) - VZ%X(2)*VLam%X(1)   /)

RPhi   = G%H + (Sqrt(Re**4 + (Rp**2 - Re**2)*(CP*gd)**2))**3/(Rp*Re**4)

!------ 1.3.a. Adjoint calculations

Do i=1,3
   VPhi_G(i)%X = &
      (/ VZ_G(i)%X(2)*VLam%X(3) + VZ%X(2)*VLam_G(i)%X(3)     &
         - VZ_G(i)%X(3)*VLam%X(2) - VZ%X(3)*VLam_G(i)%X(2),  &
         VZ_G(i)%X(3)*VLam%X(1) + VZ%X(3)*VLam_G(i)%X(1)     &
         - VZ_G(i)%X(1)*VLam%X(3) - VZ%X(1)*VLam_G(i)%X(3),  &
         VZ_G(i)%X(1)*VLam%X(2) + VZ%X(1)*VLam_G(i)%X(2)     &
         - VZ_G(i)%X(2)*VLam%X(1) - VZ%X(2)*VLam_G(i)%X(1)   /)
End Do

RPhi_G(:) = (3*Sqrt(Re**4 + (Rp**2 - Re**2)*(CP*gd)**2)/(Rp*Re**4))* &
            (Rp**2 - Re**2)*(CP*gd)*   &
            (CP_G(:)*gd + CP*gd_G(:))
RPhi_G(1) = RPhi_G(1) + 1


!----------------------------------------------------------
! 2. DERIVATIVES
!----------------------------------------------------------


!--- 2.1. Non-adjoint calculations

JGC(1,:) = VZ%X(:)
JGC(2,:) = rtd*VPhi%X(:)/RPhi
JGC(3,:) = rtd*VLam%X(:)/RLam


!------ 2.1.a. Adjoint calculations

Do i=1,3
   JGC_G(1,:,i) = VZ_G(i)%X(:)
   JGC_G(2,:,i) = rtd*(VPhi_G(i)%X(:)*RPhi - &
                       VPhi%X(:)*RPhi_G(i))/RPhi**2
   JGC_G(3,:,i) = rtd*(VLam_G(i)%X(:)*RLam - &
                       VLam%X(:)*RLam_G(i))/RLam**2
End Do


End Subroutine Jacobian_GC_adj



!==========================================================
Function Hessian_GC &
  (G)         & ! <-- Geodetic coordinates of a point
Result(HGC)     ! --> Hessian d2(Geodetic)/d(Cartesian)2
!
! Calculation of Hessian matrix of Geodetic coordinates
! with respect to Cartesian coordinates:
! d2(z,Phi,Lambda)/d(x(j))d(x(k)).
!----------------------------------------------------------
! Method:
!   Transform of mixed coordinate derivatives to
!   cartesian coordinates.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Apr 2000 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G           ! Geodetic coordinates of a point
               ! G must not be polar.
!
! Function result:
!
Real(Double) :: &
   HGC(3,3,3)  ! Hessian d2(z,Phi,Lambda)/d(x(j))d(x(k))
!----------------------------------------------------------
! Local Scalars:
!
Integer   :: i, j, k  ! Dimension indexes
!
! Local Arrays:
!
Real(Double)   :: &
   JGC(3,3),        & ! Jacobian
   JGC_G(3,3,3)       ! Mixed hessian
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF MIXED HESSIAN
!----------------------------------------------------------

Call Jacobian_GC_adj &
  (G,     & ! <-- Geodetic coordinates of a point
   JGC,   & ! --> Jacobian d(Geodetic(i))/d(Cartesian(j))
   JGC_G)   ! --> d(JGC(i,j))/d(Geodetic(k))


!----------------------------------------------------------
! 2. TRANSFORM FROM MIXED TO CARTESIAN COORDINATES
!----------------------------------------------------------

Do i=1,3
   Do j=1,3
      Do k=1,3
         HGC(i,j,k) = Sum(JGC_G(i,j,:)*JGC(:,k))
      End Do
   End Do
End Do


End Function Hessian_GC



!==========================================================
Function Gravity &
  (Alt,   & ! <-- Altitude above reference ellipsoid [km]
   GCLat)   ! <-- Geocentric latitude [rad]
            ! --> Gravity acceleration [m/s^2]
!
! Calculation of the gravity acceleration at a given
! spatial location.
!----------------------------------------------------------
! Method:
!   Kurt Lambeck, Geophysical Geodesy, 1988
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Oct 1997 | Original code.
!   2.0   | 18 Sep 1998 | Function.
!   3.0   | 20 Sep 1998 | Argument order.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   GCLat    ! Geocentric latitude (rad)
Real(Double), Intent(In)  :: &
   Alt      ! Altitude above reference ellipsoid (km)
!
! Function result:
!
Real(Double) :: &
   Gravity  ! Gravity acceleration (m/s^2)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: g_sur  ! Surface gravity for given latitude
Real(Double) :: r_sur  ! g_sur/ge
Real(Double) :: R0     ! Effective Earth radius
Real(Double) :: s2_lat ! Sin(GCLat)**2
Real(Double) :: s2_2la ! Sin(2*GCLat)**2)/4
!----------------------------------------------------------


s2_lat  = Sin(GCLat)**2
s2_2la  = s2_lat * (1.0_Double - s2_lat)

!r_sur  = (1 + f2*(Sin(GCLat)**2) - f4/4*(Sin(2*GCLat)**2))
r_sur   = (1 + f2* s2_lat         - f4  * s2_2la          )
g_sur   = ge*r_sur

!R0     = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)*(Sin(GCLat)**2))
R0      = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)* s2_lat        )

Gravity = g_sur*(R0/(R0 + Alt))**2


End Function Gravity



!==========================================================
Function Geop_from_Alt &
   (H,     &  ! <-- Altitude above reference ellipsoid [km]
    GCLat) &  ! <-- Geocentric latitude [rad]  ]
Result(GPH)   ! --> Geopotential height [gpkm]
!
! Conversion of geometrical height to geopotential.
!----------------------------------------------------------
! Method:
!   Approximate formula from 'US Standard Atmosphere'
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Oct 1997 | Original code.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   H      ! Altitude above reference ellipsoid [km]
Real(Double), Intent(In)  :: &
   GCLat  ! Geodetic latitude [rad]
!
! Function result:
!
Real(Double)             :: &
   GPH    ! Geopotential height (gpkm)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: g_sur  ! Gravity acceleration on the surface
Real(Double) :: r_sur  ! g_sur/ge
Real(Double) :: R0     ! Effective Earth radius
Real(Double) :: s2_lat ! Sin(GCLat)**2
Real(Double) :: s2_2la ! Sin(2*GCLat)**2)/4
!----------------------------------------------------------


s2_lat  = Sin(GCLat)**2
s2_2la  = s2_lat * (1.0_Double - s2_lat)

!r_sur  = (1 + f2*(Sin(GCLat))**2 - f4/4*(Sin(2*GCLat))**2)
r_sur   = (1 + f2* s2_lat         - f4  * s2_2la          )
g_sur   = ge*r_sur

!R0     = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)*(Sin(GCLat))**2)
R0      = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)* s2_lat        )

GPH     = (g_sur/g_ave)*(R0*H/(R0 + H))


End Function Geop_from_Alt



!==========================================================
Function Alt_from_Geop &
  (GPH,   & ! <-- Geopotential height [gpkm]
   GCLat) & ! <-- Geocentric latitude [rad]
Result(H)   ! --> Altitude above reference ellipsoid [km]
!
! Conversion of geopotential height to geometrical
! altitude above reference ellipsoid.
!----------------------------------------------------------
! Method:
!   Approximate formula from 'US Standard Atmosphere'
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 10 Oct 1997 | Original code.
!   2.0   | 18 Sep 1997 | Function.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   GPH    ! Geopotential height [gpkm]
Real(Double), Intent(In)  :: &
   GCLat  ! Geodetic latitude [rad]
!
! Function result:
!
Real(Double)             :: &
   H      ! Altitude above reference ellipsoid [km]
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: g_sur  ! Gravity acceleration on the surface
Real(Double) :: r_sur  ! g_sur/ge
Real(Double) :: R0     ! Effective Earth radius
Real(Double) :: s2_lat ! Sin(GCLat)**2
Real(Double) :: s2_2la ! Sin(2*GCLat)**2)/4
!----------------------------------------------------------


s2_lat  = Sin(GCLat)**2
s2_2la  = s2_lat * (1.0_Double - s2_lat)

!r_sur  = (1 + f2*(Sin(GCLat))**2 - f4/4*(Sin(2*GCLat))**2)
r_sur   = (1 + f2* s2_lat         - f4  * s2_2la          )
g_sur   = ge*r_sur

!R0     = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)*(Sin(GCLat))**2)
R0      = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)* s2_lat        )

H       = R0*GPH/(g_sur*R0/g_ave - GPH)


End Function Alt_from_Geop



!==========================================================
Subroutine Alt_from_Geop_adj &
  (GPH,    & ! <-- Geopotential height [gpkm]
   GCLat,  & ! <-- Geocentric latitude [rad]
   H,      & ! --> Altitude above reference ellipsoid [km]
   H_GPH)    ! --> d(H)/d(GPH)
!
! Conversion of geopotential height to geometrical
! altitude above reference ellipsoid: adjoint version.
!----------------------------------------------------------
! Method:
!   Approximate formula from 'US Standard Atmosphere'
!----------------------------------------------------------
! (C) Copyright 1998-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   2.0   | 18 Sep 1997 | Basic non-adjoint version.
!   1.0   | 10 Apr 2000 | Adjoint version.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   GPH    ! Geopotential height [gpkm]
Real(Double), Intent(In)  :: &
   GCLat  ! Geodetic latitude [rad]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   H      ! Altitude above reference ellipsoid [km]
!
Real(Double), Intent(Out) :: &
   H_GPH  ! d(H)/d(GPH)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: g_sur  ! Gravity acceleration on the surface
Real(Double) :: r_sur  ! g_sur/ge
Real(Double) :: R0     ! Effective Earth radius
Real(Double) :: s2_lat ! Sin(GCLat)**2
Real(Double) :: s2_2la ! Sin(2*GCLat)**2)/4
!----------------------------------------------------------


s2_lat  = Sin(GCLat)**2
s2_2la  = s2_lat * (1.0_Double - s2_lat)

!r_sur  = (1 + f2*(Sin(GCLat))**2 - f4/4*(Sin(2*GCLat))**2)
r_sur   = (1 + f2* s2_lat         - f4  * s2_2la          )
g_sur   = ge*r_sur

!R0     = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)*(Sin(GCLat))**2)
R0      = r_sur*Re / (1 + f + m + (-3*f + 5*m/2)* s2_lat        )

H       = R0*GPH/(g_sur*R0/g_ave - GPH)

H_GPH   = (g_sur*R0**2/g_ave)/(g_sur*R0/g_ave - GPH)**2


End Subroutine Alt_from_Geop_adj



!==========================================================
Function GAST &
  (Year,   & ! <-- Year of occultation
   Month,  & ! <-- Month
   Day,    & ! <-- Day
   Hour,   & ! <-- Hour
   Minute, & ! <-- Minute
   Sec,    & ! <-- Second of Rpginning of occultation
   DSec)     ! <-- Second of occultation duration
             ! --> Frame rotation angle [rad]
!
! Calculation of Greenwich Apparent Sidereal Time Angle
! Rptween the Earth frame and absolute J2000 frame
! for rotation of the GPS/MET frame.
!----------------------------------------------------------
! Method:
!   Astronomical Alamanus, 1993
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 11 Oct 1997 | Original code.
!   2.0   | 18 Sep 1998 | Function.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Pi
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)        :: &
   Year   ! Year of occultation
Integer, Intent(In)        :: &
   Month  ! Month of occultation
Integer, Intent(In)        :: &
   Day    ! Day of occultation
Integer, Intent(In)        :: &
   Hour   ! Hour of occultation
Integer, Intent(In)        :: &
   Minute ! Minute of occultation
Real(Double), Intent(In)   :: &
   Sec    ! Second of Rpginning of occultation
Real(Double), Intent(In)   :: &
   DSec   ! Second of duration of occultation
!
! Function result:
!
Real(Double)              :: &
   GAST   ! GPS/MET frame rotation angle [rad]
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: JYear, JMonth, Cent, J
Real(Double) :: JDay
Real(Double) :: tu, gmst, utco, utc
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF JULIAN DAY NUMBER
!----------------------------------------------------------

JYear   = Year - (12-Month)/10
JMonth  = Month + 1 + 12*((12-Month)/10)
Cent    = JYear/100
J       = 2 - Cent + Cent/4 + Int(365.25_Double*JYear) + &
          Int(30.6001_Double*JMonth)
JDay    = Real(J + Day, Double) + 1720994.5_Double


!----------------------------------------------------------
! 2. CALCULATION OF GREENWICH MEAN SIDEREAL TIME
!----------------------------------------------------------

tu   = (JDay - 2451545.0_Double)/36525.0_Double
gmst = 24110.54841_Double + 8640184.812866_Double*tu + &
       0.093104_Double*tu**2 - 6.2e-6_Double*tu**3
utco = Real(Hour*3600, Double) + Real(Minute*60, Double) + Sec
utc  = (utco + DSec)*1.0027379093_Double
gmst = Modulo(gmst + utc, 86400.0_Double)


!----------------------------------------------------------
! 3. CALCULATION OF FRAME ROTATION ANGLE
!----------------------------------------------------------

GAST = gmst*2*Pi/86400.0_Double


End Function GAST



End Module Earth
