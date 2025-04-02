!
!+ GNSS Radio occultation observation operator: 1d Abel integral operator
!
#if __NEC_VERSION__ / 100 == 305    /* Work around issues with nfort 3.5.x */
!NEC$ options "-fno-reorder-logical-expression"
#endif

MODULE ECHAM_1Dvar
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   The 1D observational operator for the radiooccultation
!   measurements and its adjoint.
!   Calculation of Abel integral for one profile
!   and its derivatives wrt model variables.
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
!  optimisations, fix adjoint code
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  extended checks and diagnostics for invalid rays
! V1_15        2011/12/06 Harald Anlauf
!  modify status flag for 'borderline' rays
! V1_22        2013-02-13 Harald Anlauf
!  rename data_type -> grid_type, add vertical coordinate type
! V1_27        2013-11-08 Harald Anlauf
!  horizontal interpolation for ICON,
!  switch integration formula above model top.
! V1_31        2014-08-21 Harald Anlauf
!  improved, smoother integration level spacing for ICON
! V1_39        2015-01-07 Harald Anlauf
!  catch pathological profiles (GME only)
! V1_42        2015-06-08 Harald Anlauf
!  horint_mode
! V1_45        2015-12-15 Harald Anlauf
!  ECHAM_Refraction_1D: optimize calculation of adjoint variables
! V1_46        2016-02-05 Andreas Rhodin
!  mo_atm_state: base decisions on new flag 'vct', not 'ivctype'
! V1_51        2017-02-24 Harald Anlauf
!  ECHAM_1DVar: OpenMP optimization of adjoint calculation
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
! Diagnostics for vectorization
!
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION) && 0
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif
!==============================================================================

! Module ECHAM_1Dvar
!
! Routines for implementing 1DVar of
! bending angles
!----------------------------------------------------------
! (C) Copyright 2004, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Dec 2004 | Original code.
!----------------------------------------------------------
!         | 11 Oct 2006 | changes for 3dvar, A.Rhodin
!----------------------------------------------------------
! Modules used:
!
Use Defaults,    only: &
! Imported Parameters:
    Double, Pi

Use ICO_grid,    only: gf         ! global fields

Use mo_atm_grid, only: VCT_P_HYB  ! hybrid pressure vertical coordinate flag

Use mo_system,   only: flush
!----------------------------------------------------------
Implicit None
PRIVATE
PUBLIC :: echam_refraction_1d
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine ECHAM_Refraction_1D &
  (PRO,       & ! <-- RO impact parameter [km]
   GC,        & ! <-- Occultation point (geodetic)
   RC,        & ! <-- Curvature radius
   Vrb,       & ! <-- Verbose level
   Strict,    & ! <-- Strict checks on profile validity
   Z_oro_min, & ! <-- Min. height above surface [km]
   Dxdz_min,  & ! <-- Lower bound for dx/dz
   Dz_duct,   & ! <-- Min. distance above ducting layer [km]
   Adjoint,   & ! <-- Evaluate adjoint (derivatives)
   EM,        & ! --> Model refraction angle
   ZPRO,      & ! --> geometric height for impact parameter
   IDX,       & ! --> Subgrid indices [point, index]
   EM_T,      & ! --> d(EM(IP))/d(T(IGP,k))
   EM_Q,      & ! --> d(EM(IP))/d(Q(IGP,k))
   EM_P,      & ! --> d(EM(IP))/d(Psur(IGP))
   Stat)        ! --> Error status
!
! The 1D observational operator for the radiooccultation
! measurements and its adjoint.
!----------------------------------------------------------
! Method:
!   Calculation of Abel integral for one profile
!   and its derivatives wrt model variables.
!----------------------------------------------------------
! (C) Copyright 2004, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Dec 2004 | Original code.
!----------------------------------------------------------
!         | 11 Oct 2006 | changes for 3dvar, A.Rhodin
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
!   Cartesian,     &
! Imported Operators:
    Operator(*),   &
    Operator(+),   &
    Operator(-),   &
    Operator(/),   &
    Operator(.x.), &
    Assignment(=)
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use ECHAM_fields, only: &
! Imported Routines:
    Interpolate_Refractivity
!
Use ECHAM_fields_adj, only: &
! Imported Routines:
    Interpolate_Refractivity_adj
!
Use mo_wmo_tables, only: &
! Imported Parameters:
!   WMO6_LATLON,         &
!   WMO6_GAUSSIAN,       &
!   DWD6_ICON,           &
    DWD6_ICOSAHEDRON
!
#if defined (NAG) || defined(__PGI)
! For systems where the complementary error function is not defined
! we need an accurate fallback version.
! (On the NEC SX, we definitely want the vectorized intrinsic.)
!
use mo_algorithms, only: erfc => derfc
!
#endif
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)      :: &
   PRO(:)        ! Impact parameters [km]
!
Type(Geodetic), Intent(In)    :: &
   GC            ! Occultation point (geodetic)
!
Real(Double), Intent(In)      :: &
   RC            ! Curvature radius
!
Integer, Optional, Intent(In) :: &
   Vrb           ! Verbose level
!
Logical, Optional, Intent(in) :: &
   Strict        ! Apply strict checks on profile validity
!
Real(Double), Intent(In)      :: &
   Z_oro_min     ! Min. height above surface [km]
!
Real(Double), Intent(In)      :: &
   Dxdz_min      ! Lower bound for dx/dz

Real(Double), Optional, Intent(In) :: &
   Dz_duct       ! Min. distance above ducting layer

Logical, Optional, Intent(in) :: &
   Adjoint       ! Evaluate adjoint (derivatives)
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   EM(:)         ! Model refraction angles
!
Real(Double), Intent(Out) :: &
   ZPRO(:)       ! Geometric height for impact parameter
!
Integer, Intent(Out)      :: &
   IDX(:,:)      ! Subgrid indices [point, index]
!
Real(Double), Intent(Out) :: &
   EM_T(:,:,:)   ! d(EM(IP))/d(T(IGP,k))
!
Real(Double), Intent(Out) :: &
   EM_Q(:,:,:)   ! d(EM(IP))/d(Q(IGP,k))
!
Real(Double), Intent(Out) :: &
   EM_P(:,:)     ! d(EM(IP))/d(Psur(IGP))
!
Integer, Intent(Out)      :: &
   Stat(:)       ! Error status
                 ! 0 = OK
                 ! 1 = bending angle could not be evaluated
                 ! 2 = bending angle unreliable (numerically unstable)
                 ! 3 = bending angle extrapolated/below minimum height
!----------------------------------------------------------
#if defined (_CRAYFTN) || ((__GNUC__ * 100 + __GNUC_MINOR__) >= 409) \
  || (defined (__NEC__) && (__NEC_VERSION__ >= 30004))
! Assert compiler that these can be optimized:
contiguous :: IDX, EM, EM_T, EM_Q, EM_P
#endif
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   Zatm = 120.0_Double, & ! Maximum integration height
   DZ   = 0.05_Double     ! Integration step
!
! Split integration range, switching from constant integration
! step DZ (below Zsplit) to increasing linearly from DZ
! at Zsplit to DZtop at the upper integration limit.
! (Setting Zsplit = Zatm recovers the original constant step size).
!
Real(Double), Parameter :: &
   Zsplit = 30.0_Double, &! Height for change of integration steps
!  Zsplit = Zatm,        &! (Debug) Don't change integration steps
   DZtop  = 50 * DZ       ! Integration step size at top level (2.5 km)
!
! Parameter of exponential level spacing above splitting level (new scheme)
Real(Double), Parameter :: &
   DZ0    = 1.028_Double * DZ
!
Real(Double), Parameter :: &
   Zswitch = 80.0_Double  ! Height for switching integration formulas
                          ! (above model top or where dT/dz > 0)
Integer :: Nswitch        ! Index for switching
!
! Local Scalars:
!
! --- Dimensions
!
Integer        :: NLev    ! Number of levels
Integer        :: NZ      ! Number of integration steps
Integer        :: IZ      ! Index of vertical integration grid
Integer        :: KP      ! Number of impact parameters
Integer        :: IP      ! Impact parameter index
Integer        :: IZP     ! Integration start index
Integer        :: NGP     ! Subgrid dimension
Integer        :: IGP     ! Auxiliary subgrid index
Integer        :: k       ! Auxiliary level index
!
! --- Coordinates
!
Real(Double)   :: ZP      ! Altitude of point [km]
Real(Double)   :: Zmin    ! Minimum model Z for this lon/lat
Real(Double)   :: Zmax    ! Maximum model Z for this lon/lat
Real(Double)   :: PLon    ! Longitude
Real(Double)   :: PLat    ! Latitude
!
! --- Integrand variables
!
Real(Double)   :: NP       ! Interpolated N
Real(Double)   :: RX       ! R(X)
Real(Double)   :: NX       ! n(X)
Real(Double)   :: NR       ! dn/dR
Real(Double)   :: NRR      ! d2n/dR2
Real(Double)   :: F_N      ! dF/dn
Real(Double)   :: F_NR     ! dF/dNR
Real(Double)   :: X1, X2   ! Integration limits
Real(Double)   :: X1I, X2I ! Interpolation limits
Real(Double)   :: S1, S2   ! Sqrt(X1,2**2 - P**2)
Real(Double)   :: W1, W2   ! Weight for F1 and F2
Real(Double)   :: xm, sm   ! Averages    of X, S
Real(Double)   :: dx, ds   ! Differences of X, S
Real(Double)   :: xi, xi2  ! Expansion coefficient
Real(Double)   :: lmx      ! Temporary for expansion of logarithm
Real(Double)   :: fac, tmp ! Temporaries
!
! --- Work variables
!
Character(Len=80) :: Line  ! Terminal line
Integer        :: NZ1      ! Integration steps (lower part)
Integer        :: NZ2      ! Integration steps (upper part)
Integer        :: IZ0      ! Index of height where     Z < ZMIN + Z_oro_min
Integer        :: IG0      ! Index of height where dX/dZ < dXdZ_min
Integer        :: Verb     ! Verbosity level
Logical        :: Lstrict  ! Strictness check
Logical        :: Ladjoint ! Evaluate adjoint (derivatives)
Real(Double)   :: PROmin   ! Minimum impact parameter
Real(Double)   :: Z0       ! Reference height for interval generation
Real(Double)   :: alf, bet ! Parameters for interval generation
Real(Double)   :: lambda   ! Scale parameter of exponential decay
Real(Double)   :: e1,e2,r1,r2,t1,t2,u1,u2,p1,p2 ! Temporaries
Logical        :: lcompat  ! Backward compatibility mode for operational GME
integer        :: k_asr    ! Level index above superrefraction layer
real(Double)   :: x_asr    ! First level above superrefraction layer
real(Double)   :: dx_asr   ! Distance above superrefraction layer
!
!
! Local Arrays:
!
Real(Double), Allocatable :: &
   Z(:),                & ! Altitude grid
   X(:),                & ! Refractive radius grid
   X_Z(:),              & ! d(X(IZ))/d(Z(IZ))
   S(:)                   ! Sqrt(X**2 - P**2)
Real(Double) :: &
   NG(3),               & ! Interpolated dN/d(alt,lat,lon)
   NH(3,3)                ! Interpolated hessian matrix of N
Real(Double), Allocatable :: &
   NP_T(:,:),           & ! d(NP)/d(T(IGP,k))
   NP_Q(:,:),           & ! d(NP)/d(Q(IGP,k))
   NP_P(:),             & ! d(NP)/d(Psur(IGP))
   NG_T(:,:,:),         & ! d(NG(i))/d(T(IGP,k))
   NG_Q(:,:,:),         & ! d(NG(i))/d(Q(IGP,k))
   NG_P(:,:)              ! d(NG(i))/d(Psur(IGP))
Real(Double), Allocatable :: &
   F(:),                & ! NR/(n*(n+R*NR))
   F_T(:,:,:),          & ! d(F(IZ))/d(T(IGP,k))
   F_Q(:,:,:),          & ! d(F(IZ))/d(Q(IGP,k))
   F_P(:,:)               ! d(F(IZ))/d(Psur(IGP))
Real(Double), Allocatable :: &
   X_T(:,:,:),          & ! d(X(IZ))/d(T(IGP,k))
   X_Q(:,:,:),          & ! d(X(IZ))/d(Q(IGP,k))
   X_P(:,:)               ! d(X(IZ))/d(Psur(IGP))
Real(Double), Allocatable :: &
   W1_(:),              & ! Integration weights from lower bounds W1(IZ)
   W2_(:),              & ! Integration weights from upper bounds W2(IZ)
   W1_X(:,:),           & ! d(W1(IZ))/d(X(IZ)), d(W1(IZ))/d(X(IZ+1))
   W2_X(:,:)              ! d(W2(IZ))/d(X(IZ)), d(W2(IZ))/d(X(IZ+1))
Logical,      Allocatable :: &
   err(:)                 ! Error mask

Integer ::              & ! Internal work array:
   IDX_(gf% ngp,3)        ! Subgrid indices [point, index]
!
! Auxiliary arrays for code optimization (vectorization)
!
Real(Double), Allocatable :: &
   em_t_ (:,:),         & ! Section of EM_T etc., transposed
   em_q_ (:,:)
Real(Double), Allocatable :: &
   f_t_  (:,:,:),       & ! Same as F_T etc., but with
   f_q_  (:,:,:),       & ! permuted dimensions
   x_t_  (:,:,:),       & ! (NLev,NZ,NGP)
   x_q_  (:,:,:)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

Verb = 0
If (present (Vrb))     Verb = Vrb

Lstrict = .false.
If (present (Strict))  Lstrict = Strict

Ladjoint = .true.
If (present (Adjoint)) Ladjoint = Adjoint

dx_asr = -1._Double
if (present (Dz_duct)) dx_asr = Dz_duct

!--- 1.1. Dimensions

NLev = gf% nz

KP   = Size(PRO)
if (KP == 0) return

NGP = gf% ncol  ! Number of actually used column(s)

! Backward compatibility with operational GME
lcompat = (gf% Grid_Type == DWD6_ICOSAHEDRON .and. gf% vctype == VCT_P_HYB)


!--- 1.2. Computation of latitude and longitude

PLon = GC%Lambda
PLat = GC%Phi


!--- 1.3. Altitude grid size and dimension
FTRACE_BEGIN("ECHAM_Refraction_1D:1.3")

ZP = 0.0_Double

Call Interpolate_Refractivity  &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP)         ! --> Interpolated N


!-- Number of integration intervals with constant spacing
Zmin = floor (Zmin/DZ)*DZ     ! fix grid for consistency of nonlinear/adjoint
Zmin = min (Zmin, 0._Double)  ! Enable extrapolation for lowest ray(s)
NZ1  = 1 + Ceiling ((min (Zsplit, Zatm) - Zmin) / DZ)
NZ1  = max (NZ1, 1)
Z0   = Zmin + (NZ1-1)*DZ

if (lcompat) then
!-- Determine parameters for linearly increasing spacing (OLD)
   alf = DZ
   NZ2 = nint (2*(Zatm - Z0) / (alf + DZtop))
   NZ2 = max (NZ2, 0)
   bet = (Zatm - Z0 - NZ2*alf) / max (NZ2, 1)**2
else
!-- Determine parameters for exponentially increasing spacing (NEW)
!   Z(k) = Z(k0) + alf * (exp ((k-k0)*bet) - 1)
   alf = DZ0 * (Zatm - DZtop - Z0) / (DZtop - DZ0)
   bet = log (1 + DZ0 / alf)
   NZ2 = nint (log (1 + (Zatm-Z0) / alf) / bet)
   bet = log (1 + (Zatm-Z0) / alf) / max (NZ2, 1)
!  print *, "NZ2,alf,bet,dz0,dzt=", NZ2,alf,bet, &
!           alf*(exp(bet)-1),alf*exp(NZ2*bet)*(1-exp(-bet))
end if

NZ  = NZ1 + NZ2


FTRACE_END  ("ECHAM_Refraction_1D:1.3")

!--- 1.4. Integrand and its derivatives

Allocate(Z(NZ))
Allocate(X(NZ))
Allocate(F(NZ))
Allocate(X_Z(NZ))
Allocate(NP_T(NGP,NLev))
Allocate(NP_Q(NGP,NLev))
Allocate(NP_P(NGP))
Allocate(NG_T(3,NGP,NLev))
Allocate(NG_Q(3,NGP,NLev))
Allocate(NG_P(3,NGP))

if (Ladjoint) then
   Allocate(F_T(NZ,NGP,NLev))
   Allocate(F_Q(NZ,NGP,NLev))
   Allocate(F_P(NZ,NGP))
   Allocate(X_T(NZ,NGP,NLev))
   Allocate(X_Q(NZ,NGP,NLev))
   Allocate(X_P(NZ,NGP))
end if
!
! Auxiliary arrays with "rotated dimensions"
!
if (Ladjoint) then
   allocate(f_t_(NLev,NZ,NGP))
   allocate(f_q_(NLev,NZ,NGP))
   allocate(x_t_(NLev,NZ,NGP))
   allocate(x_q_(NLev,NZ,NGP))
end if


!-- Set up integration intervals: constant spacing
Do IZ = 1, NZ1
   Z(IZ) = Zmin + (IZ-1)*DZ
End Do

if (lcompat) then
!-- Linearly increasing spacing
   Do IZ = NZ1+1, NZ
      tmp   = IZ - nz1
      Z(IZ) = Z0 + alf * tmp + bet * tmp*tmp
   end do
else
!-- Grid with exponentially increasing spacing (NEW)
!   Z(k) = Z(k0) + alf * (exp ((k-k0)*bet) - 1)
   Do IZ = NZ1+1, NZ
      tmp   = IZ - nz1
      Z(IZ) = Z0 + alf * (exp (tmp*bet) - 1._Double)
   End Do
end if
!print *, "NZ2, dz(nz1..)=", NZ2,Z(NZ1)-Z(NZ1-1),Z(NZ1+1)-Z(NZ1),Z(NZ1+2)-Z(NZ1+1),&
!     "..",Z(NZ)-Z(NZ-1)

!write(0,*) "middle     =", z(nz1)
!if (nz1 < nz) &
!write(0,*) "middle(+1) =", z(nz1+1)
!write(0,*) "top   (-1) =", z(nz-1)
!write(0,*) "top        =", z(nz)
!write(0,*)
!NZ1 = NZ    ! Debug: enforce original integration weights below


PROmin = minval (PRO)
IZP    = 1
IZ0    = 0

Do IZ=NZ,1,-1

FTRACE_BEGIN("ECHAM_Refraction_1D:1.4int")
   Call Interpolate_Refractivity_adj  &
     (PLon,     & ! <-- Longitude of point [deg]
      PLat,     & ! <-- Latitude  of point [deg]
      Z(IZ),    & ! <-- Altitude  of point [km]
      Zmin,     & ! --> Minimum model Z for this lon/lat
      Zmax,     & ! --> Maximum model Z for this lon/lat
      NP,       & ! --> Interpolated N
      NG,       & ! --> Interpolated dN/d(alt,lat,lon)
      NH,       & ! --> Interpolated hessian matrix of N
      IDX_,     & ! --> Subgrid indices [point, index]
      NP_T,     & ! --> d(NP)/d(T(IGP,k))
      NP_Q,     & ! --> d(NP)/d(Q(IGP,k))
      NP_P,     & ! --> d(NP)/d(Psur(IGP))
      NG_T,     & ! --> d(NG(i))/d(T(IGP,k))
      NG_Q,     & ! --> d(NG(i))/d(Q(IGP,k))
      NG_P)       ! --> d(NG(i))/d(Psur(IGP))
FTRACE_END  ("ECHAM_Refraction_1D:1.4int")

FTRACE_BEGIN("ECHAM_Refraction_1D:1.4set")
   RX    = RC + Z(IZ)
   NX    = 1.0_Double + NP
   X(IZ) = RX*NX

   NR    = NG(1)
   NRR   = NH(1,1)

   X_Z(IZ) = (NX + RX*NR)

   F(IZ) = NR/(NX*(NX + RX*NR))
   F_N   = -NR*(2*NX + RX*NR)    &
            / (NX*(NX + RX*NR))**2
   F_NR  = 1.0_Double/(NX + RX*NR)**2

if (lcompat) then
if (F(IZ) >= 0 .and. Z(IZ) > 10) then
   write (0,'(2(A,f8.3))') " At: lat=", PLat, "  lon=", PLon
   write (0,*) "ECHAM_Refraction_1D: F(IZ)>=0:",IZ,F(IZ),Z(IZ),NR,NX,(NX+RX*NR)
end if
end if

if (Ladjoint) then

#if defined (__NEC__)
! Manual loop exchange to work around issues in nfort <= 3.0.4
!NEC$ loop_count(4)
   Do IGP = 1, NGP
!NEC$ ivdep
   Do k = 1, Nlev
      F_T(IZ,IGP,k) = F_N*NP_T(IGP,k) + F_NR*NG_T(1,IGP,k)
      F_Q(IZ,IGP,k) = F_N*NP_Q(IGP,k) + F_NR*NG_Q(1,IGP,k)
   End Do
   End Do
#else

   ! "Original" derivatives
!  F_T(IZ,:,:) = F_N*NP_T(:,:)*(NX/(NX + RX*NR)) + &
!                F_NR*(NG_T(1,:,:) - NP_T(:,:)*RX*NRR/(NX + RX*NR))
!  F_Q(IZ,:,:) = F_N*NP_Q(:,:)*(NX/(NX + RX*NR)) + &
!                F_NR*(NG_Q(1,:,:) - NP_Q(:,:)*RX*NRR/(NX + RX*NR))

   ! "Naive" derivatives
   F_T(IZ,:,:) = F_N*NP_T(:,:) + F_NR*NG_T(1,:,:)
   F_Q(IZ,:,:) = F_N*NP_Q(:,:) + F_NR*NG_Q(1,:,:)
#endif

   ! "Original" derivatives
!  F_P(IZ,:)   = F_N*NP_P(:)*(NX/(NX + RX*NR)) + &
!                F_NR*(NG_P(1,:) - NP_P(:)*RX*NRR/(NX + RX*NR))

   ! "Naive" derivatives
   F_P(IZ,:)   = F_N*NP_P(:)   + F_NR*NG_P(1,:)

   !-------------------------------------------------------------------
   ! Calculate derivatives of the integral boundaries w.r.t. the fields
   !-------------------------------------------------------------------
   X_T(IZ,:,:) = RX * NP_T(:,:)
   X_Q(IZ,:,:) = RX * NP_Q(:,:)
   X_P(IZ,:)   = RX * NP_P(:)

end if ! Ladjoint

FTRACE_END  ("ECHAM_Refraction_1D:1.4set")

   !--------------------------------------------------
   ! Debug cases with extrapolation for the lowest ray
   ! (Signs chosen such that good values are positive)
   !--------------------------------------------------
   if (Z(IZ) < ZMIN .and. IZ0 == 0) then
      IZ0 = IZ
      if (Verb > 1) then
         write(0,*)
         write(0,'(A,2F9.3)') &
              "ECHAM_Refraction_1D: Warning: Z < ZMIN   @ Lat,Lon=", Plat, Plon
         write(0,'(A,3F9.3,9ES12.3)') &
              "Z,ZMIN,X,-NR,NX+RX*NR,-F,F_N:", &
              Z(IZ),ZMIN,X(IZ)-RC,-NR,NX+RX*NR,-F(IZ),F_N
         write(0,'(A,F9.3)') "Z_oro_min (active) =", Z_oro_min
      end if
   end if

   if (z_oro_min > 0._Double) then
      if (Z(IZ) < (ZMIN + Z_oro_min) .and. IZ0 == 0)  IZ0 = IZ
   end if
   !---------------------------------------------------------------
   ! Skip interpolation for intervals below lowest impact parameter
   !---------------------------------------------------------------
   if (X(IZ) < PROmin) then
      IZP = IZ
      exit
   end if

End Do

if (Ladjoint) then
   do IGP = 1, NGP
      f_t_(:,:,IGP) = transpose (F_T(1:NZ,IGP,1:NLev))
      f_q_(:,:,IGP) = transpose (F_Q(1:NZ,IGP,1:NLev))
      x_t_(:,:,IGP) = transpose (X_T(1:NZ,IGP,1:NLev))
      x_q_(:,:,IGP) = transpose (X_Q(1:NZ,IGP,1:NLev))
   end do
end if


IDX(1:gf% ncol,:) = IDX_(1:gf% ncol,:)  ! Copy back used column indices only

X(:IZP-1) = 0._Double

!---------------------------------------------------
! Check for presence of ducting layer in model state
!---------------------------------------------------
k_asr = 0
x_asr = 0._Double
if (dx_asr >= 0._Double) then
   do IZ = NZ, IZP, -1
      if (X_Z(IZ) <= 0) then
         k_asr = IZ + 1         ! First level above ...
         x_asr = X(k_asr)       ! ... superrefraction layer
         exit
      end if
   end do
end if

if (k_asr > 0 .and. verb > 1) then
   IG0 = k_asr - 1
   write(0,*)
   write(0,'(A,2F9.3)') &
        "ECHAM_Refraction_1D: Ducting: (dX/dZ < 0) @ Lat,Lon=", Plat, Plon
   write(0,'(A,3F9.3)') "ECHAM_Refraction_1D: lowest impact heights:", &
        (PRO(IZ)-RC, iz = 1, min (size(pro), 3))
   do IZ = max (IG0,IZP), min (IG0+11,NZ)
      write(0,'(A,2F9.3,9ES12.3)')                       &
           "Z,X, dX/dZ,-F:", Z(IZ),X(IZ)-RC,X_Z(IZ),-F(IZ)
   end do
end if

!-------------------------------------------------------------------
! Check for non-monotonicity leading to divergent vertical gradients
!-------------------------------------------------------------------
IG0 = 0
tmp = max (Dxdz_min, 0._Double)
do IZ = NZ, IZP, -1
   if (X_Z(IZ) <= tmp) then
      IG0 = IZ
      exit
   end if
end do
if (IG0 > 0) then
   if (X_Z(IG0) <= 0._Double .and. k_asr == 0) then
      write(0,*)
      write(0,'(A,2F9.3)') &
           "ECHAM_Refraction_1D: Warning: dX/dZ < 0  @ Lat,Lon=", Plat, Plon
      if (Verb > 0) then
         write(0,'(A,3F9.3,9ES12.3)')                           &
              " Z, ZMIN, Z0, (NX+RX*NR), -F:",                  &
              Z(IG0),ZMIN,Z(1),X_Z(IG0),-F(IG0)
         write(0,'(A,3F9.3,9ES12.3)')                           &
              "Amin,X,X(+1),dX/dZ,dX/dZ(+1):",                  &
              PROmin-RC,X(IG0)-RC,X(IG0+1)-RC,X_Z(IG0),X_Z(IG0+1)
         do IZ = max (IG0-10,IZP), min (IG0+10,NZ)
            write(0,'(A,2F9.3,9ES12.3)')                        &
                 "Z,X, dX/dZ,-F:", Z(IZ),X(IZ)-RC,X_Z(IZ),-F(IZ)
         end do
         write(0,'(A,F9.3)') "Dxdz_min (active) =", Dxdz_min
      end if
   end if
end if

!-----------------------------------------------------------------------
! Set level for switching integration formulas.  (ICON: do not switch in
! the stratosphere where the vertical temperature gradient is positive.)
!-----------------------------------------------------------------------
Nswitch = NZ1
if (gf% vctype /= VCT_P_HYB) then
   Nswitch = NZ
   if (Z(NZ) > Zswitch) then
      Do IZ = NZ1, NZ-1
         if (Z(IZ) > Zswitch) then
            Nswitch = IZ
            exit
         end if
      End Do
   end if
end if

!----------------------------------------------------------
! 2. COMPUTATION OF BENDING ANGLE
!----------------------------------------------------------

Allocate(S(NZ))
Allocate(w1_(NZ-1))
Allocate(w2_(NZ-1))
Allocate(err(NZ-1))
Allocate(W1_X(2,NZ-1))
Allocate(W2_X(2,NZ-1))

if (Ladjoint) then
   allocate(em_t_(NLev,NGP))
   allocate(em_q_(NLev,NGP))
end if


!--- 2.0. Initialization

Stat(1:KP)     =  0
ZPRO(1:KP)     = -1.0_Double  ! invalid value: -1 km
EM  (1:KP)     =  0.0_Double

EM_T(1:KP,:,:) =  0.0_Double
EM_Q(1:KP,:,:) =  0.0_Double
EM_P(1:KP,:)   =  0.0_Double


Angles: Do IP = 1,KP

   !--- Explicitly exclude rays too close to superrefraction layer

   if (k_asr > 0 .and. PRO(IP) < x_asr + dx_asr) then
      if (Verb > 0) then
         write(0,'(A,2F9.3,I4,2F9.3)')                            &
              "ECHAM_Refraction_1D: Reject: Lat,Lon,IP,P,x_asr=", &
              Plat,Plon,IP,PRO(IP)-RC,x_asr-RC
      end if
      stat(IP) = 1
      cycle angles
   end if


   !--- 2.1. Computation of S(X)

   Where (X(:) > PRO(IP))
!     S(:) = Sqrt(X(:)**2 - PRO(IP)**2)
      !-----------------------------------------
      ! This variant is numerically more robust:
      !-----------------------------------------
      S(:) = Sqrt ((X(:)+PRO(IP)) * (X(:)-PRO(IP)))
   Elsewhere
      S(:) = 0.0_Double
   End Where


   !--- 2.2. Determination of starting point

   IZP = Sum(MaxLoc(Z(:), Mask = (X(:) <= PRO(IP))))

   If (IZP > NZ .or. IZP < 1) then
      Stat(IP) = 1
      If (Verb >= 3) then
         Write (Line,'(2X,A,I6,2X,A,F8.4,2X,A,F8.4,2X,A,F8.4,2X,A,I2)') &
              'IP = ', IP, 'P = ', PRO(IP)-RC,                          &
              'Xmin-A =  ', X(1)-PRO(IP), 'Zmin =', Zmin, 'Stat = ', Stat(IP)
         write (0,'(a)') trim (line)
      End If
   End If

   If (IZP <= IG0 .or. (IZP <= IZ0 .and. Lstrict)) then
      if (Verb > 1) then
         write(0,'(A,3F9.3,4I6)')                                       &
              "ECHAM_Refraction_1D: Reject : Lat,Lon,P,IP,IZ,IG0,IZ0=", &
              Plat,Plon,PRO(IP)-RC,IP,IZP,IG0,IZ0
      end if
      ! We evaluate the bending angle but flag it accordingly
      If (Stat(IP) == 0) then
         if (IZP <= IZ0) Stat(IP) = 3   ! Impact height below threshold
         if (IZP <= IG0) Stat(IP) = 2   ! Gradient(r*n) below threshold
      end If
   End If

   If (IZP < 1 .or. IZP >= NZ) Cycle Angles

   ZPRO(IP) = Z(IZP)

   !--- 2.3. Setting initial condition

!  X2I = X(IZP)
!  X2  = PRO(IP)
!  S2  = 0.0_Double
   err(IZP:) = .false.

   !--- 2.4. Numerical integration
FTRACE_BEGIN("ECHAM_Refraction_1D:2.4")

! Debug potentially pathological profiles that do not decay exponentially:
if (lcompat) then
iz = max (Nswitch, IZP)
if (maxval (f(iz:nz)) >= 0) then
   k = maxloc (f(iz:nz), 1) + iz-1
   write(0,'(2(A,f8.3))') " At: lat=", PLat, "  lon=", PLon
   write(0,*) "ECHAM_Refraction_1D: WARNING: f>=0:",f(k),f(k-1),f(k+1),k,iz,izp,nswitch,nz
   write(0,*) "ECHAM_Refraction_1D:              :",z(k),z(k-1),z(k+1)
end if
if (minval (f(iz:nz-1) / f(iz+1:nz), 1) <= 0) then
   k = minloc (f(iz:nz-1) / f(iz+1:nz), 1) + iz-1
   write(0,'(2(A,f8.3))') " At: lat=", PLat, "  lon=", PLon
   write(0,*) "ECHAM_Refraction_1D: WARNING: arg<=0:",k,f(k),f(k+1),iz,izp,nswitch,nz
   !----------------------------
   ! Scan for pathological range
   !----------------------------
   do k = nz-1, izp, -1
      if (f(k) / f(k+1) <= 0) then
         nswitch = k+1
         exit
      end if
   end do
   write(0,*) "ECHAM_Refraction_1D: WARNING: using nswitch=", nswitch
end if
end if

   !------------------------------------------------------------------
   ! Determine integration weights in region with exponential decrease
   !------------------------------------------------------------------
   Do IZ = max (Nswitch, IZP), NZ-1

      X1  = X(IZ)
      if (IZ == IZP) then
         X1 = PRO(IP)
      end if
      X2  = X(IZ+1)

      X1I = X(IZ)
      X2I = X(IZ+1)

      if (F(IZ)/F(IZ+1) > 1._Double) then
         lambda  = log (F(IZ)/F(IZ+1)) / (X2I - X1I)
      else
         lambda  = EPSILON (lambda)
         err(IZ) = .true.
      end if
      t1  =       lambda*(X1I - PRO(IP))
      t2  =       lambda*(X2I - PRO(IP))
      p1  = sqrt (lambda*(X1I + PRO(IP)))
      p2  = sqrt (lambda*(X2I + PRO(IP)))
      r1  = sqrt (lambda*(X1  - PRO(IP)))
      r2  = sqrt (lambda*(X2  - PRO(IP)))
      e1  = sqrt (Pi) * erfc (r1)
      e2  = sqrt (Pi) * erfc (r2)
      u1  = r1 * exp (- r1**2)
      u2  = r2 * exp (- r2**2)
      fac = -2*PRO(IP) / (lambda*(X2I - X1I))
      W1  =  fac * exp (t1) * ((0.5_Double-t2)*(e2-e1) + (u2-u1)) / p1
      W2  = -fac * exp (t2) * ((0.5_Double-t1)*(e2-e1) + (u2-u1)) / p2

      W1_(IZ) = W1
      W2_(IZ) = W2
      !----------------------------------------------------------------
      ! Derivatives of W1, W2 w.r.t. integration subinterval boundaries
      !----------------------------------------------------------------
!     if (r1 > 0._Double) then
      if (X1 == X1I) then
         tmp =     2*PRO(IP) * lambda / (p1 * r1)
      else
         tmp =     0
      end if
      W1_X(1,IZ) = (lambda + 1._Double  / (X2I - X1I)           &
                           - 0.5_Double / (X1I + PRO(IP))) * W1 &
                 + tmp
      W2_X(2,IZ) = (lambda - 1._Double  / (X2I - X1I)           &
                           - 0.5_Double / (X2I + PRO(IP))) * W2 &
                 - 2*PRO(IP) * lambda / (p2 * r2)

      W1_X(2,IZ) = - fac * exp (t1) * ((0.5_Double-t1)*(e2-e1) + (u2-u1)) &
                     / (p1 * (X2I - X1I))
      W2_X(1,IZ) = - fac * exp (t2) * ((0.5_Double-t2)*(e2-e1) + (u2-u1)) &
                     / (p2 * (X2I - X1I))
   End Do

   if (any (err(IZP:))) then
      if (Verb > 1) then
         write(0,*)
         write(0,*) "ECHAM_Refraction_1D: Integrand not decreasing exponentially!"
         write(0,*) "Problem occurred", count (err(IZP:)), "times"
         write(0,'(2(A,f8.3))') " At: lat=", PLat, "  lon=", PLon
         write(0,*) "    NZ1, NZ2, NZ, IZP, Nswitch =", NZ1, NZ2, NZ, IZP, Nswitch
         write(0,*) "    IZ     X(IZ)     Z(IZ)     S(IZ)   |F(IZ)|     |F(IZ+1)|"
         do IZ = max (NZ1, IZP), NZ-1
            if (err(IZ)) then
               write(0,'(i7,3F10.3,2ES12.3)') IZ,X(IZ),Z(IZ),S(IZ),-F(IZ),-F(IZ+1)
            end if
         end do
         write(0,*) "Switching to alternative integration formula."
      end if
      !-------------------------------------------------------
      ! Error occurred, switch to original integration formula
      ! for all potentially affected intervals.
      !-------------------------------------------------------
      k = NZ
      do IZ = max (NZ1, IZP), NZ-1
         if (err(IZ)) then
            k = IZ+1
         end if
      end do
      Nswitch = k
   end if

   !------------------------------------------------------------------
   ! Determine integration weights for arbitrary case, but small steps
   !------------------------------------------------------------------
   Do IZ = IZP, Nswitch-1

      X1  = X(IZ)
      if (IZ == IZP) then
         X1 = PRO(IP)
      end if
      X2  = X(IZ+1)

      X1I = X(IZ)
      X2I = X(IZ+1)

      S1  = S(IZ)
      S2  = S(IZ+1)

!--------------------------------------------------------------------------
! The original expression for the weights W1,W2 is numerically unstable
! for high vertical resolution, because the argument of the logarithm
! is then very close to 1:
!
!     W1  = -2*PRO(IP)*(S1 - S2 + X2I*Log((X2 + S2)/(X1 + S1)))/(X2I - X1I)
!     W2  = -2*PRO(IP)*(S2 - S1 - X1I*Log((X2 + S2)/(X1 + S1)))/(X2I - X1I)
!
! We replace it by a suitable expansion of the logarithm and perform
! the cancellations of leading terms explicitly.
!--------------------------------------------------------------------------

      fac = -2*PRO(IP) / (X2I - X1I)
      xm  = (X2 + X1) * 0.5_Double
      sm  = (S1 + S2) * 0.5_Double
      dx  =  X2 - X1
      ds  =  S2 - S1
      xi  = (dx + ds)  / (xm  + sm)
      xi2 = xi*xi
      if (xi2 < 1.e-3_Double) then
!        lmx = Log ((2._Double + xi)/(2._Double - xi)) - xi
         !----------------------------------------------------
         ! Use series expansion of [log ((2+xi)/(2-xi)) - xi]:
         !----------------------------------------------------
         lmx = ((((  1._Double/ 2304 )*xi2 &
                   + 1._Double/  448 )*xi2 &
                   + 1._Double/   80 )*xi2 &
                   + 1._Double/   12 )*xi2 * xi
         W1  = fac*(0.5_Double*dx*xi + X2*lmx + (X2I-X2)*(xi+lmx))
         W2  = fac*(0.5_Double*dx*xi - X1*lmx - (X1I-X1)*(xi+lmx))
      else
         ! Use original expression
         tmp = Log ((X2 + S2)/(X1 + S1))
         W1  = fac*(- ds + X2I*tmp)
         W2  = fac*(  ds - X1I*tmp)
      end if

      w1_(IZ) = W1
      w2_(IZ) = W2

      !----------------------------------------------------------------
      ! Derivatives of W1, W2 w.r.t. integration subinterval boundaries
      ! W1_X(i,IZ) = dW1(IZ)/dX(IZ+(i-1)),
      ! W2_X(i,IZ) = dW2(IZ)/dX(IZ+(i-1)), i=1,2
      !----------------------------------------------------------------
      r1 = W1 / (X2I - X1I)
      r2 = W2 / (X2I - X1I)
      if (S1 > 0._Double) then
         W1_X(1,IZ) =   r1 + 2*PRO(IP) / S1
      else
         W1_X(1,IZ) =   r1
      end if
      W1_X(2,IZ)    =   r2
      W2_X(1,IZ)    = - r1
      W2_X(2,IZ)    = - r2 - 2*PRO(IP) / S2
   End Do

FTRACE_BEGIN("ECHAM_Refraction_1D:2.4int")

L1:Do IZ = IZP, NZ-1
      W1 = w1_(IZ)
      W2 = w2_(IZ)

      EM(IP) = EM(IP) + W1*F(IZ) + W2*F(IZ+1)

      if (Ladjoint) then
#ifdef _CRAYFTN
!DIR$ NEXTSCALAR
#endif
!NEC$ loop_count(4)
      EM_P(IP,:)   = EM_P(IP,:)   + W1*F_P(IZ,:)   + W2*F_P(IZ+1,:) &
                   + ( F(IZ)  *W1_X(1,IZ)                           &
                     + F(IZ+1)*W2_X(1,IZ)) * X_P(IZ  ,:)            &
                   + ( F(IZ)  *W1_X(2,IZ)                           &
                     + F(IZ+1)*W2_X(2,IZ)) * X_P(IZ+1,:)

      end if ! Ladjoint
   End Do L1 ! IZ

   if (Ladjoint) then

#if defined (__NEC__)

L2:Do IZ = IZP, NZ-1
      W1 = w1_(IZ)
      W2 = w2_(IZ)

! Manual loop exchange to work around issues in nfort <= 3.0.4
!NEC$ loop_count(4)
L3:Do IGP = 1, NGP

!!NEC$ select_vector
!NEC$ ivdep
      Do k = 1, Nlev
!CDIR ARRAYCOMB
         EM_T(IP, IGP ,k) = EM_T(IP, IGP ,k) &
                          + W1*F_T(IZ, IGP ,k) + W2*F_T(IZ+1, IGP ,k) &
                          + ( F(IZ)  *W1_X(1,IZ)                      &
                            + F(IZ+1)*W2_X(1,IZ)) * X_T(IZ  , IGP ,k) &
                          + ( F(IZ)  *W1_X(2,IZ)                      &
                            + F(IZ+1)*W2_X(2,IZ)) * X_T(IZ+1, IGP ,k)
         EM_Q(IP, IGP ,k) = EM_Q(IP, IGP ,k) &
                          + W1*F_Q(IZ, IGP ,k) + W2*F_Q(IZ+1, IGP ,k) &
                          + ( F(IZ)  *W1_X(1,IZ)                      &
                            + F(IZ+1)*W2_X(1,IZ)) * X_Q(IZ  , IGP ,k) &
                          + ( F(IZ)  *W1_X(2,IZ)                      &
                            + F(IZ+1)*W2_X(2,IZ)) * X_Q(IZ+1, IGP ,k)
!CDIR END ARRAYCOMB
      End Do

   End Do L3 ! IGP

   End Do L2 ! IZ

#else

!$omp parallel do private(IGP,IZ,W1,W2)
#ifdef _CRAYFTN
!DIR$ NEXTSCALAR
!DIR$ BLOCKINGSIZE(3)
#endif
L3:Do IGP = 1, NGP

   em_t_(:,IGP) = 0._Double    ! Initialize work arrays
   em_q_(:,IGP) = 0._Double

L2:Do IZ = IZP, NZ-1
      W1 = w1_(IZ)
      W2 = w2_(IZ)

      em_t_(:,IGP) = em_t_(:,IGP) + W1*f_t_(:,IZ,IGP) + W2*f_t_(:,IZ+1,IGP)   &
                   + ( F(IZ)  *W1_X(1,IZ)                                     &
                     + F(IZ+1)*W2_X(1,IZ)) * x_t_(:,IZ  ,IGP)                 &
                   + ( F(IZ)  *W1_X(2,IZ)                                     &
                     + F(IZ+1)*W2_X(2,IZ)) * x_t_(:,IZ+1,IGP)
      em_q_(:,IGP) = em_q_(:,IGP) + W1*f_q_(:,IZ,IGP) + W2*f_q_(:,IZ+1,IGP)   &
                   + ( F(IZ)  *W1_X(1,IZ)                                     &
                     + F(IZ+1)*W2_X(1,IZ)) * x_q_(:,IZ  ,IGP)                 &
                   + ( F(IZ)  *W1_X(2,IZ)                                     &
                     + F(IZ+1)*W2_X(2,IZ)) * x_q_(:,IZ+1,IGP)

! Previous loop body (after simple loop interchange)
!     EM_T(IP,IGP,:) = EM_T(IP,IGP,:) + W1*F_T(IZ,IGP,:) + W2*F_T(IZ+1,IGP,:) &
!                  + ( F(IZ)  *W1_X(1,IZ)                                     &
!                    + F(IZ+1)*W2_X(1,IZ)) * X_T(IZ  ,IGP,:)                  &
!                  + ( F(IZ)  *W1_X(2,IZ)                                     &
!                    + F(IZ+1)*W2_X(2,IZ)) * X_T(IZ+1,IGP,:)
!     EM_Q(IP,IGP,:) = EM_Q(IP,IGP,:) + W1*F_Q(IZ,IGP,:) + W2*F_Q(IZ+1,IGP,:) &
!                  + ( F(IZ)  *W1_X(1,IZ)                                     &
!                    + F(IZ+1)*W2_X(1,IZ)) * X_Q(IZ  ,IGP,:)                  &
!                  + ( F(IZ)  *W1_X(2,IZ)                                     &
!                    + F(IZ+1)*W2_X(2,IZ)) * X_Q(IZ+1,IGP,:)

! Original loop body (before switching loop order):
!
!     EM_T(IP,:,:) = EM_T(IP,:,:) + W1*F_T(IZ,:,:) + W2*F_T(IZ+1,:,:) &
!                  + ( F(IZ)  *W1_X(1,IZ)                             &
!                    + F(IZ+1)*W2_X(1,IZ)) * X_T(IZ  ,:,:)            &
!                  + ( F(IZ)  *W1_X(2,IZ)                             &
!                    + F(IZ+1)*W2_X(2,IZ)) * X_T(IZ+1,:,:)
!     EM_Q(IP,:,:) = EM_Q(IP,:,:) + W1*F_Q(IZ,:,:) + W2*F_Q(IZ+1,:,:) &
!                  + ( F(IZ)  *W1_X(1,IZ)                             &
!                    + F(IZ+1)*W2_X(1,IZ)) * X_Q(IZ  ,:,:)            &
!                  + ( F(IZ)  *W1_X(2,IZ)                             &
!                    + F(IZ+1)*W2_X(2,IZ)) * X_Q(IZ+1,:,:)

   End Do L2 ! IZ

   End Do L3 ! IGP
!$omp end parallel do

   EM_T(IP,:,:) = transpose (em_t_)
   EM_Q(IP,:,:) = transpose (em_q_)

#endif

   end if   ! Ladjoint

FTRACE_END  ("ECHAM_Refraction_1D:2.4int")
FTRACE_END  ("ECHAM_Refraction_1D:2.4")


   !--- 2.5. Displaying progress message

!  If (Verb >= 3) then
   If (Verb >= 3 .or. EM(IP) <= 0._Double) then
      Write (Line,'(2X,A,I6,2X,A,F8.4,2X,A,ES13.6,2X,A,I2,A1)')   &
         'IP = ', IP, 'P = ', PRO(IP)-RC,                         &
         'EM = ', EM(IP), 'Stat = ', Stat(IP)
      write(0,'(a)') trim(line)
      call flush (0)
   End If

End Do Angles

If (Verb >= 3) then
   Write(0,'()')
End If


!----------------------------------------------------------
! 3. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(S)
Deallocate(w1_)
Deallocate(w2_)
Deallocate(Z)
Deallocate(X)
Deallocate(F)

if (Ladjoint) then
   Deallocate(NP_T)
   Deallocate(NP_Q)
   Deallocate(NP_P)
   Deallocate(NG_T)
   Deallocate(NG_Q)
   Deallocate(NG_P)
   Deallocate(F_T)
   Deallocate(F_Q)
   Deallocate(F_P)
end if

End Subroutine ECHAM_Refraction_1D

End Module ECHAM_1Dvar
