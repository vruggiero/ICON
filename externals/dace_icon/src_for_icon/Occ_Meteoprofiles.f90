!
!+ GNSS Radio occultation observation operator: Calculation of profiles
!
MODULE Occ_Meteoprofiles
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Calculation of refractivity, temperature, pressure,
!   and humidity profiles.
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
! V1_11        2010/06/16 Harald Anlauf
!  N_from_TPQ, T_from_NPQ, Q_from_TPN: use 3-term refractivity formula
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_43        2015-08-19 Harald Anlauf
!  implement refractivity expression from Aparicio & Laroche (2011)
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
! Module Occ_Meteoprofiles
!
! Calculation of refractivity, temperature, pressure,
! and humidity profiles.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Oct 1998 | Original code.
!   2.0   | 06 Jul 1999 | Rd, Rnu, C1, and C2 in Earth.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!
Use Earth, only: &
! Imported Parameters:
    Rd, Rnu, C1, C2, C3
!----------------------------------------------------------
Implicit None
Private
Public :: aq, bq, NQ_to_TP, T_from_NPQ, Rd, Eps
!----------------------------------------------------------
! Public Parameters:
!
      ! q  - specific humidity [kg/kg]
      ! T  - temperature [K]
      ! Tv - virtual temperature [K]
      ! P  - pressure [mb]
      ! Pw - water vapour pressure [mb]
      ! N  - refractivity [absolute]
!
Real(Double), Parameter :: &
   eps =       &           ! T = Tv/(1 + eps*q)
      Rnu/Rd - 1.0_Double
!
Real(Double), Parameter :: &
   aq = Rd/Rnu,   &        ! Pw = P*q/(aq + bq*q)
   bq = 1.0_Double - aq
!----------------------------------------------------------
! Private Scalars:
!
Real(Double), Private :: &
   GCLat                   ! Geocentric latitude [rad]
!
! Private Arrays:
!
Real(Double), Private, Pointer :: &
   PZ(:),           &
   PN(:),  P2N(:),  &      ! Pointers for access to arrays
   PQ(:),  P2Q(:)          ! N, Q, D2N, D2Q for interpolation.
!----------------------------------------------------------
! Public Procedures:
public :: N_from_TPQ            ! Refractivity from T,P,Q
public :: N_from_TPQ_AL2011     ! Refractivity after Aparicio & Laroche (2011)
public :: N_from_TPQ_AL2011_adj ! incl. adjoint version
!==========================================================
Contains
!==========================================================
Elemental &
Function N_from_TPQ &
  (T,     & ! <-- Temperature [K]
   P,     & ! <-- Pressure [mb]
   Q)     & ! <-- Specific humidity [kg/kg]
Result(N)   ! --> Refractivity [dimensionless]
!
! Calculation of refractivity from temperature, pressure,
! and humidity.
!----------------------------------------------------------
! Method:
!
!                 P            P Q                  P Q
!   N(T,P,Q) = C1--- + C2----------------- + C3--------------
!                 T       T**2 (aq + bq Q)      T (aq + bq Q)
!
!   (Bean and Datton, Radiometeorology)
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   T    ! Temperature [K]
!
Real(Double), Intent(In) :: &
   P    ! Pressure [mb]
!
Real(Double), Intent(In) :: &
   Q    ! Specific humdity [kg/kg]
!
! Function result:
!
Real(Double) :: &
   N    ! Refractivity
!----------------------------------------------------------

#ifdef _SMITH_WEINTRAUB_TWOTERM

N = (C1 + C2*Q/(T*(aq + bq*Q)))*P/T

#else

Real(Double) :: Pw_r  ! Pw/P

Pw_r = Q/(aq + bq*Q)

N    =  (C1 +   C2*Pw_r/T + C3*Pw_r)*P/T

#endif

End Function N_from_TPQ



!==========================================================
Function T_from_NPQ &
  (N,     & ! <-- Refractivity [dimensionless]
   P,     & ! <-- Pressure [mb]
   Q)     & ! <-- Specific humidity [kg/kg]
Result(T)   ! --> Temperature [K]
!
! Calculation of temperature from refractivity, pressure,
! and humidity.
!----------------------------------------------------------
! Method:
!   Analytic solution of equation N = N(T,P,Q).
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   N    ! Refractivity
!
Real(Double), Intent(In) :: &
   P    ! Pressure [mb]
!
Real(Double), Intent(In) :: &
   Q    ! Specific humdity [kg/kg]
!
! Function result:
!
Real(Double) :: &
   T    ! Temperature [K]
!----------------------------------------------------------

#ifdef _SMITH_WEINTRAUB_TWOTERM

T = (C1*P + Sqrt((C1*P)**2 + 4*C2*N*P*Q/(aq + bq*Q)))/(2*N)

#else

Real(Double) :: Pw_r  ! Pw/P

Pw_r = Q/(aq + bq*Q)

T = ((C1 + C3*Pw_r)*P + Sqrt(((C1 + C3*Pw_r)*P)**2 + 4*C2*N*P*Pw_r))/(2*N)

#endif

End Function T_from_NPQ



!==========================================================
Function Q_from_TPN &
  (T,     & ! <-- Temperature [K]
   P,     & ! <-- Pressure [mb]
   N)     & ! <-- Refractivity [dimensionless]
Result(Q)   ! --> Specific humidity [kg/kg]
!
! Calculation of humidity from temperature, pressure,
! and refractivity.
!----------------------------------------------------------
! Method:
!   Analytic solution of equation N = N(T,P,Q).
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Sep 1998 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   T    ! Temperature [K]
!
Real(Double), Intent(In) :: &
   P    ! Pressure [mb]
!
Real(Double), Intent(In) :: &
   N    ! Refractivity
!
! Function result:
!
Real(Double) :: &
   Q    ! Specific humdity [kg/kg]
!----------------------------------------------------------

#ifdef _SMITH_WEINTRAUB_TWOTERM

Q = aq*T*(C1*P - N*T)/(bq*N*T**2 - (C1*bq*T + C2)*P)

#else

Q = aq*T*(C1*P - N*T)/(bq*N*T**2 - (C1*bq*T + C2 + C3*T)*P)

#endif

End Function Q_from_TPN



!==========================================================
Elemental &
Function N_from_TPQ_AL2011 &
  (T,     & ! <-- Temperature       [K]
   P,     & ! <-- Pressure          [Pa]
   Q)     & ! <-- Specific humidity [kg/kg]
Result(N)   ! --> Refractivity [dimensionless]
!----------------------------------------------------------
! Calculation of refractivity from temperature, pressure,
! and humidity using the expression given in
! J.M. Aparicio and S. Laroche, JGR 116, D11104 (2011).
!----------------------------------------------------------
! Method:
!
!   N = (n-1) = 10^-6 * N0 * (1 + 10^-6/6 * N0),
!
!   N0(T,P,Q) = rho_d*(B1 + B2/T) + rho_w*(B3 + B4/T)
!
!   where rho_d = rho * (1-Q), rho_w = Q * rho,
!   and rho has been calculated using CIPM-2007.
!----------------------------------------------------------
Use CIPM_2007, only: Z           ! Compressibility factor
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: T    ! Temperature      [K]
Real(Double), Intent(In) :: P    ! Pressure         [Pa]
Real(Double), Intent(In) :: Q    ! Specific humdity [kg/kg]
!
! Function result:
!
Real(Double) ::             N    ! Refractivity [dimensionless]
!----------------------------------------------------------
  !====================
  ! Private parameters:
  !====================
  real(double), parameter :: &
       a1 =  222.682_double, &
       a2 =    0.069_double, &
       a3 = 6701.605_double, &
       a4 = 6385.886_double, &
       t0 =  273.15_double,  &
       b1 = a1 - a2,         &
       b2 = t0 * a2,         &
       b3 = a3 - a4,         &
       b4 = t0 * a4

  !----------------
  ! Local variables
  !----------------
  real(double) :: xv            ! Water vapor mixing ratio
  real(double) :: rho           ! Density
  real(double) :: N0            ! First-order refractivity [N-units]

  xv   = Q / (aq + bq*Q)
  rho  = P * (1 - bq*xv) / (Rd * T * Z(P,T,Q))

  N0   = rho * ((1-Q) * (b1 + b2 / T) + Q * (b3 + b4 / T))

  N    = 1.e-6_Double * N0 * (1 + (1.e-6_Double/6) * N0)

End Function N_from_TPQ_AL2011
!==========================================================
Elemental &
Subroutine N_from_TPQ_AL2011_adj &
  (T,     & ! <-- Temperature       [K]
   P,     & ! <-- Pressure          [Pa]
   Q,     & ! <-- Specific humidity [kg/kg]
   N,     & ! --> Refractivity [dimensionless]
   N_T,   & ! --> d(N)/d(T)
   N_P,   & ! --> d(N)/d(P)
   N_Q)     ! --> d(N)/d(Q)
!----------------------------------------------------------
! Calculation of refractivity from temperature, pressure,
! and humidity using the expression given in
! J.M. Aparicio and S. Laroche, JGR 116, D11104 (2011).
! Adjoint version of N_from_TPQ_AL2011.
!----------------------------------------------------------
Use CIPM_2007, only: Z_adj        ! Compressibility factor (adjoint version)
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: T    ! Temperature [K]
Real(Double), Intent(In)  :: P    ! Pressure [Pa]
Real(Double), Intent(In)  :: Q    ! Specific humdity [kg/kg]
!
! Output arguments:
!
Real(Double), Intent(Out) :: N    ! Refractivity [dimensionless]
Real(Double), Intent(Out) :: N_T  ! d(N)/d(T) [K^-1]
Real(Double), Intent(Out) :: N_P  ! d(N)/d(P) [Pa^-1]
Real(Double), Intent(Out) :: N_Q  ! d(N)/d(Q) [dimensionless]
!----------------------------------------------------------
  !====================
  ! Private parameters:
  !====================
  real(double), parameter :: &
       a1 =  222.682_double, &
       a2 =    0.069_double, &
       a3 = 6701.605_double, &
       a4 = 6385.886_double, &
       t0 =  273.15_double,  &
       b1 = a1 - a2,         &
       b2 = t0 * a2,         &
       b3 = a3 - a4,         &
       b4 = t0 * a4

  !----------------
  ! Local variables
  !----------------
  real(double) :: Z             ! Compressibility factor
  real(double) :: Z_P, Z_T, Z_Q ! and derivatives
  real(double) :: xv            ! Water vapor mixing ratio
  real(double) :: xv_Q          ! d(xv)/d(Q)
  real(double) :: rho           ! Density
  real(double) :: rho_T         ! d(rho)/d(T)
  real(double) :: rho_P         ! d(rho)/d(P)
  real(double) :: rho_Q         ! d(rho)/d(Q)
  real(double) :: N0            ! First-order refractivity [N-units]
  real(double) :: N0_T          ! d(N0)/d(T)
  real(double) :: N0_P          ! d(N0)/d(P)
  real(double) :: N0_Q          ! d(N0)/d(Q)

  call Z_adj (P, T, Q, Z, Z_P, Z_T, Z_Q)

  xv    = Q  / (aq + bq*Q)
  xv_Q  = aq / (aq + bq*Q)**2

  rho   = P   * (1 - bq*xv) / (Rd * T * Z)
  rho_T = rho * (-1 / T                   - Z_T / Z)
  rho_P = rho * ( 1 / P                   - Z_P / Z)
  rho_Q = rho * (-bq * xv_Q / (1 - bq*xv) - Z_Q / Z)

  N0    = rho   * ((1-Q) * (b1 + b2 / T) + Q * (b3 + b4 / T))

  N0_T  = rho_T * ((1-Q) * (b1 + b2 / T) + Q * (b3 + b4 / T))  + &
          rho   * ((1-Q) *     (-b2)     + Q *     (-b4))/T**2

  N0_P  = rho_P * ((1-Q) * (b1 + b2 / T) + Q * (b3 + b4 / T))

  N0_Q  = rho_Q * ((1-Q) * (b1 + b2 / T) + Q * (b3 + b4 / T))  + &
          rho   * (  -     (b1 + b2 / T) +     (b3 + b4 / T))

  N     = 1.e-6_Double * N0   * (1 + (1.e-6_Double/6) * N0)

  N_T   = 1.e-6_Double * N0_T * (1 + (1.e-6_Double/3) * N0)
  N_P   = 1.e-6_Double * N0_P * (1 + (1.e-6_Double/3) * N0)
  N_Q   = 1.e-6_Double * N0_Q * (1 + (1.e-6_Double/3) * N0)

End Subroutine N_from_TPQ_AL2011_adj
!==========================================================
Subroutine NQ_to_TP &
  (GDLat, & ! <-- Geodetic latitude [deg]
   Z,     & ! <-- Altitude above reference ellipsoid [km]
   N,     & ! <-- Profile of refractivity
   Q,     & ! <-- Profile of specific humidity
   T,     & ! --> Profile of temperature
   P)       ! ~~> Profile of pressure
!
! Calculation of profiles of temperature and pressure
! profiles of refactivity and humdity.
!----------------------------------------------------------
! Method:
!   Calculation of P(z) from numerical integration of
!   barometric formula:
!
!   d ln(P(z))                    g(z)
!   ---------- = - -------------------------------------
!      dz           Rd T(N(z),P(z),Q(z)) (1 + eps(Q(z))
!
!   with given N(z) and Q(z) and known dependence
!   T(N,P,Q). Calculation of T from profiles of
!   N(z), P(z), and Q(z):  T(z) = T(N(z),P(z),Q(z)).
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Sep 1998 | Original code
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Routines:
    GCLat_from_GDLat, Gravity
!
Use Interpolation, only: &
! Imported Routines:
    Init_Spline, Spline
!
Use Dif_equations, only: &
! Imported Routines:
    Runge_Kutta
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   GDLat    ! Geodetic latitude [deg]
!
Real(Double), Target, Intent(In) :: &
   Z(1:)    ! Altitudes above reference ellipsoid [km]
            ! in increasing order.
!
Real(Double), Intent(In) :: &
   N(1:)    ! Gridded refractivity
            ! Must be positive.
!
Real(Double), Target, Intent(In) :: &
   Q(1:)    ! Gridded humidity [kg/kg]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   T(1:)    ! Gridded temperature [K]
!
Real(Double), Optional, Intent(Out) :: &
   P(1:)    ! Gridded pressure [mb]
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   Z_max = 120.0_Double ! Beginning of integration
Integer, Parameter      :: &
   KZ    = 1200         ! Number of integeration steps
!
! Local Scalars:
!
Real(Double) :: DZ      ! Integration step
Real(Double) :: ZP      ! Integration variable
Real(Double) :: LnNZ    ! Interpolated Ln(N(Z))
Real(Double) :: DLnNZ   ! Derivative of Ln(N(Z))
Real(Double) :: QZ      ! Interpolated Q(Z)
Integer      :: i       ! Integration step index
!
! Local Arrays:
!
Real(Double) :: LnP(1)  ! Ln(P) for numerical integration.
Real(Double), Target :: &
   LnN(Size(Z)),  &     ! Ln(N(Z))
   D2N(Size(Z)),  &     ! Spline coefficients for ln(N(Z))
   D2Q(Size(Z))         ! Spline coefficients for Q(Z)
Real(Double) :: &
   ZI(0:KZ),      &     ! Integration grid of Z
   TZ(0:KZ),      &     ! Calculated T(Z)
   D2T(0:KZ),     &     ! Spline coefficients for T(Z)
   LnPZ(0:KZ),    &     ! Calculated Ln(P(Z))
   D2P(0:KZ)            ! Spline coefficients for Ln(P(Z))
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------


! 1.1. Calculation of geocentric latitude

GCLat  = GCLat_from_GDLat(GDLat)


! 1.2. Calculation of spline coefficients

LnN(:) = Log(N(:))

Call Init_Spline(Z, LnN, D2N)
Call Init_Spline(Z, Q,   D2Q)


! 1.3. Setting pointers for access to arrays

PZ  => Z
PN  => LnN
P2N => D2N
PQ  => Q
P2Q => D2Q


! 1.4. Setting initial conditions for integration

DZ      = -Z_max/KZ
ZP      = Z_max
ZI(KZ)  = ZP

Call Spline(Z, LnN, D2N, Z_max, LnNZ, DLnNZ)
TZ(KZ)   = -1d3*Gravity(ZP, GCLat)/(Rd*DLnNZ)
LnP(1)   = LnNZ + Log(TZ(KZ)/C1)
LnPZ(KZ) = LnP(1)


!----------------------------------------------------------
! 2. NUMERICAL INTEGRATION
!----------------------------------------------------------

Integration: Do i=KZ-1,0,-1
   Call Runge_Kutta(FTZP, DZ, ZP, LnP)
   Call Spline(PZ, PN, P2N, ZP, LnNZ)
   Call Spline(PZ, PQ, P2Q, ZP, QZ)
   ZI(i)   = ZP
   TZ(i)   = T_from_NPQ(Exp(LnNZ), Exp(LnP(1)), QZ)
   LnPZ(i) = LnP(1)
End Do Integration


!----------------------------------------------------------
! 3. INTERPOLATION OF CALCULATED T(Z) AND P(Z)
!----------------------------------------------------------

Call Init_Spline(ZI, TZ,   D2T)
Do i=1,Size(Z)
   Call Spline(ZI, TZ,   D2T, Z(i), T(i))
End Do

If (Present(P)) then
   Call Init_Spline(ZI, LnPZ, D2P)
   Do i=1,Size(Z)
      Call Spline(ZI, LnPZ, D2P, Z(i), P(i))
   End Do
   P(:) = Exp(P(:))
End If


!----------------------------------------------------------
! 4. NULLIFYING POINTERS
!----------------------------------------------------------

Nullify(PZ, PN, P2N, PQ, P2Q)


End Subroutine NQ_to_TP



!==========================================================
Function FTZP &
  (ZP,    & ! <-- Altitude abouve reference ellipsoid [km]
   LnP)     ! <-- Logarithm of pressure [mb]
            ! --> Right part of barometric formula

!
! Calculation of right part of barometric formula
! for calculation of T(Z) from N(Z) and Q(Z)
!----------------------------------------------------------
! Method:
!
!   d ln(P(z))                    g(z)
!   ---------- = - -------------------------------------
!      dz           Rd T(N(z),P(z),Q(z)) (1 + eps(Q(z))
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 01 Oct 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!
Use Interpolation, only: &
! Imported Routines:
    Spline
!
Use Earth, only: &
! Imported Routines:
    Gravity
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   ZP                 ! Altitude above reference ellipsoid
!
Real(Double), Intent(In) :: &
   LnP(1:)            ! Logarithm of pressure [mb]
!
! Function result:
!
Real(Double)             :: &
   FTZP(Size(LnP))    ! Right part of barometric formula
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: LnNZ  ! Interpolated ln(N(z))
Real(Double) :: QZ    ! Interpolated Q(z)
Real(Double) :: TZP   ! T(z,P(z))
Real(Double) :: GZ    ! Gravity acceleration g(z,lat)
!----------------------------------------------------------
! Global variables used:
!
!Real(Double), Private :: &
!   GCLat                   ! Geocentric latitude [rad]
!Real(Double), Private, Pointer :: &
!   PZ(:),           &      ! Pointers for access to arrays
!   PN(:),  P2N(:),  &      ! Z, LnN, Q, D2N, D2Q
!   PQ(:),  P2Q(:)          ! for interpolation.
!----------------------------------------------------------


Call Spline(PZ, PN, P2N, ZP, LnNZ)
Call Spline(PZ, PQ, P2Q, ZP, QZ)

TZP     = T_from_NPQ(Exp(LnNZ), Exp(LnP(1)), QZ)
GZ      = Gravity(ZP, GCLat)

FTZP(1) = -1d3*GZ/(Rd*TZP*(1 + eps*QZ))

End Function FTZP



End Module Occ_Meteoprofiles
