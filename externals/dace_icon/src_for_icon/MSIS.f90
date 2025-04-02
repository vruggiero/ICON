!
!+ GNSS Radio occultation observation operator: MSIS climate atmospheric model
!
MODULE MSIS
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   MSIS climate atmospheric model
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
! V1_9         2010/04/20 Harald Anlauf
!  optimize for SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Make FTRACE_REGIONs depending on macro DISABLE_FTRACE_REGION
! V1_14        2011/11/08 Harald Anlauf
!  Minor optimizations for NEC SX-9
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


!
! Module MSIS
!
! MSIS climate atmospheric model
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Nov 1996 | F77 version,
!         |             | S. Syndergaard (DMI).
!         |    Jun 1997 | adapted for inclusion in F90 module,
!         |             | G. Kirchengast, IMG/UoG.
!         |    Aug 1997 | vertical curvature of refractivity,
!         |             | G. Kirchengast, IMG/UoG.
!   1.0   | 28 Oct 1998 | Fortran 90 version.
!   2.0   | 20 Mar 1999 | Error check.
!   3.0   | 24 Jun 1999 | MSIS_Geop, MSIS_Pressure_Levels,
!         |             | MSIS_Num_Levels.
!         | 12 Oct 2006 | get path from mo_run_params, A.Rhodin
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Scalar Variables:
    Double
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Exit_Callee,   &
    Error
Use mo_run_params, only: &
    path_file,           &! concatenate path/filename
    Path => data          ! path for MSISxx.asc files
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
Private
Public :: Msis_Refractivity, Msis_Init, Msis_Num_Levels, &
          Msis_Pressure_Levels, Msis_Geop
!
! Public Parameters:
!
Integer, Parameter :: &
   err_Bad_Month  = 4001, &
   err_No_Coef    = 4002, &
   err_Read_Coef  = 4003
!
Real(Double), Parameter :: &
   LgScale = 1/3.0_Double  ! Lg pressure step
!
! Private Parameters:
!
Integer, Parameter, Private  :: &
   klim = 31,     &
   nlim = 7           ! Numbers of expansion coefficients
!
! Private Arrays:
!
Real(Double), Private, save ::  &
   Ac0(0:nlim,0:klim),  &
   Ac1(nlim,0:klim),    &
   Bc1(nlim,0:klim),    &
   Aa0(0:nlim),         &
   Aa1(nlim),           &
   Ba1(nlim),           &
   Ab0(0:nlim),         &
   Ab1(nlim),           &
   Bb1(nlim),           &
   Ad0(0:nlim),         &
   Ad1(nlim),           &
   Bd1(nlim)    ! Expansion coefficients.
Real(Double), Private, save ::  &
   n10(nlim),           &
   n11(nlim),           &
   n20(nlim),           &
   n21(nlim),           &
   n10x(nlim),          &
   n11x(nlim)   ! Constants for Clenshaw's recurrence formula.
Real(Double), Private, save :: &
   CosPhi,   SinPhi,  &
   CosTheta, SinTheta    ! Sin/Cos(Lat/Lon)
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine MSIS_Init &
  (Mon,     & ! <-- Month
   ErrStat)   ! --> Error status
!
! Initialization of MSIS atmospheric model.
!----------------------------------------------------------
! Method:
!   Reading coefficients from one of files
!   MSIS01.asc .. MSIS12.asc (a file for a month).
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Nov 1996 | F77 version,
!         |             | S. Syndergaard (DMI).
!   1.0   | 28 Oct 1998 | Fortran 90 version.
!   1.1   | 06 Jul 1999 | Accurate test for read errors.
!----------------------------------------------------------
! Modules used:
!
Use IO, only: &
! Imported Routines:
    FileOpen
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In)  :: &
   Mon      ! Month number
            ! Must be 1..12
!
! Output arguments:
!
Type(Error_Status), Pointer :: &
   ErrStat  ! Error status
!----------------------------------------------------------
! Local Parameters:
!
!Character(Len=*), Parameter :: &
!   Path =   &   ! MSISxx.asc file path
!      '../data/'                      ! DWD 3dvar
!      '/pf/m/mo/m214029/F90/MSIS/'    ! Regen
!      '/net/gardiken/scr6/scratch/m214089/F90/MSIS/'  ! Gardiken
!       'd:/F90/MSIS/'                 ! PC D:
!       'c:/F90/MSIS/'                 ! PC C:
!       'c:/Gorbunov/F90/MSIS/'        ! PC C: in Graz
!
! Local Scalars:
!
Character(Len=20) :: FileName  ! Coefficient file name
Integer           :: IU        ! File unit number
Integer           :: n         ! Coefficient index
Integer           :: Stat      ! Error code
!----------------------------------------------------------
! Saved value of previous month
Integer, save     :: cur_month = -1
!----------------------------------------------------------
! Global variables used:
!
! 1. Expansion coefficients.
! 2. Constants for Clenshaw's recurrence formula.
!----------------------------------------------------------


!----------------------------------------------------------
! 0. CHECK MONTH NUMBER
!----------------------------------------------------------

Call Enter_Callee &
  ('MSIS_Init',   & ! <-- User routine
   ErrStat)         ! <-> Pointer to callee status

If ((Mon < 1) .or. (Mon > 12)) then
   Call Exit_Callee &
     (ErrStat,           & ! <-> Pointer to callee status
      err_Bad_Month,     & ! <~~ User error code
      0,                 & ! <~~ System error code
      'Bad month number')  ! <~~ Error description
   Return
End If

!-- Return immediately if same month as in last call
if (Mon == cur_month) then
   Call Exit_Callee (ErrStat)
   return
end if
cur_month = Mon


!----------------------------------------------------------
! 1. READ COEFFICIENTS FOR SPECIFIED MONTH
!----------------------------------------------------------

Write(FileName,'(A,I2.2,A)') 'MSIS', Mon, '.asc'

Call FileOpen(  &
   path_file(Path,FileName),& ! <-- File pathname
   Status   = 'OLD',        & ! <~~ Open status
   Access   = 'SEQUENTIAL', & ! <~~ Access mode
   Form     = 'FORMATTED',  & ! <~~ File form
   RecL     = 5096,         & ! <~~ Record length
   Position = 'REWIND',     & ! <~~ Open position
   Action   = 'READ',       & ! <~~ Read/write permissions
   Unit     = IU,           & ! --> Unit number
   IOStat   = Stat)           ! --> Error code

If (Stat /= 0) then
   write(0,*) 'MSIS_Init: cannot open:',trim(path_file(Path,FileName))
   Call Exit_Callee &
     (ErrStat,           & ! <-> Pointer to callee status
      err_No_Coef,       & ! <~~ User error code
      0,                 & ! <~~ System error code
      'Open MSISxx.asc')   ! <~~ Error description
   Return
End If

ReadData: Do n=1,1
   Read (IU,*,IOStat=Stat) Ac0
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ac1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Bc1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Aa0
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Aa1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ba1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ab0
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ab1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Bb1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ad0
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Ad1
   If (Stat /= 0) Exit
   Read (IU,*,IOStat=Stat) Bd1
   If (Stat /= 0) Exit
End Do ReadData

Close (Unit = IU)

If (Stat /= 0) then
   Call Exit_Callee &
     (ErrStat,           & ! <-> Pointer to callee status
      err_Read_Coef,     & ! <~~ User error code
      0,                 & ! <~~ System error code
      'Read MSISxx.asc')   ! <~~ Error description
   Return
End If


!----------------------------------------------------------
! 2. INITIALIZATION OF CONSTANTS FOR CLENSHAW'S FORMULA
!----------------------------------------------------------

Do n=1,nlim
   n10(n)=(2*n+1)/dble(n+1)
   n11(n)=(2*n+1)/dble(n)
   n20(n)=(n+1)/dble(n+2)
   n21(n)=(n+2)/dble(n+1)
End Do


Call Exit_Callee(ErrStat)


End Subroutine MSIS_Init



!==========================================================
Subroutine MSIS_Refractivity &
  (G,      & ! <-- Geodetic coordinates
   N,      & ! --> Refractivity
   NG,     & ! ~~> Refractivity gradient
   NH)       ! ~~> Refractivity hessian
!
! Calulation of MSIS refractivity, its gradient and hessian.
!----------------------------------------------------------
! Method:
!   Invoke of MSIS_Expand, finite difference for gradient
!   and hessian.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 29 Oct 1998 | Original code
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G        ! Geodetic coordinates
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   N        ! Refractivity [absolute]
!
Real(Double), Intent(Out), Optional :: &
   NG(3)    ! Refractivity gradient
!
Real(Double), Intent(Out), Optional :: &
   NH(3,3)  ! Refractivity hessian
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   DH = 0.01_Double   ! Finite difference step
!
! Local Scalars:
!
Real(Double) :: &
   Nhi, Nlo           ! Refractivities for finite differences
!----------------------------------------------------------


N = MSIS_Expand(G)

If (Present(NG) .or. Present(NH)) then
   Nhi = MSIS_Expand(Geodetic(G%H + DH, G%Phi, G%Lambda))
   Nlo = MSIS_Expand(Geodetic(G%H - DH, G%Phi, G%Lambda))

   If (Present(NG)) then
      NG(1)   = (Nhi - Nlo)/(2*DH)
      NG(2:3) = 0
   End If

   If (Present(NH)) then
      NH(:,:) = 0
      NH(1,1) = (Nhi - 2*N + Nlo)/DH**2
   End If
End If


End Subroutine MSIS_Refractivity



!==========================================================
Real(Double) Function MSIS_Expand &
  (G)       ! <-- Geodetic coordinates
            ! --> MSIS refactivity
!
! Calculation of MSIS refractivity at given
! point
!----------------------------------------------------------
! Method:
!   Chebyshev polynomials and spherical harmonics
!   expansions.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Nov 1996 | F77 version,
!         |             | S. Syndergaard (DMI).
!   1.0   | 28 Oct 1998 | Fortran 90 version.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Scalar Variables:
    dtr
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G                  ! Geodetic coordinates
!
! Function result:
!
!Real(Double) :: &
!  MSIS_Expand  ! Refracitivity (absolute units)
!----------------------------------------------------------
! Local Parameters:
!
Real(Double) :: &
   zref = 100d0,   &
   h0   = 100d0
!
! Local Scalars:
!
Real(Double) :: Phi      ! Spheric longitude [rad]
Real(Double) :: Theta    ! Spheric latitude [rad]
Type(Geodetic), Save :: &
   G_old =    &          ! Last call coordinates
      Geodetic(0, -999, -999)
Integer      :: k        ! Expansion indices
Real(Double), Save   :: &
   Alfa, Beta, Refrac0   ! Work variables
Real(Double) :: x, z     ! Changed variable
Real(Double) :: &
   dk, dk1, dk2          ! Variable for recurrence
Real(Double) :: fx       ! Chebyshev's polynomial
Real(Double) :: Hs       ! Refractivity scale height
!
! Local Arrays:
!
Real(Double), Save :: &
   c(0:klim)
!
!----------------------------------------------------------
! Global variables used:
!
! 1. Expansion coefficients.
! 2. Constants for Clenshaw's recurrence formula.
! 3. Sin/Cos(Lat/Lon).
!----------------------------------------------------------


!----------------------------------------------------------
! 1. RECALCULATION OF LAT/LON-DEPENDENT PARAMETERS
!----------------------------------------------------------

If ((G_old%Phi /= G%Phi) .or. (G_old%Lambda /= G%Lambda)) then

FTRACE_BEGIN("MSIS_Expand:1")
   Phi      = G%Lambda*dtr
   CosPhi   = Cos(Phi)
   SinPhi   = Sin(Phi)
   Theta    = (90-G%Phi)*dtr
   CosTheta = Cos(Theta)
   SinTheta = Sin(Theta)

   n10x(:) = CosTheta*n10(:)
   n11x(:) = CosTheta*n11(:)

   Do k=0,klim
      c(k) = Spheric(Ac0(:,k), Ac1(:,k), Bc1(:,k))
   End Do

   Alfa    = Spheric(Aa0, Aa1, Ba1)
   Beta    = Spheric(Ab0, Ab1, Bb1)
   Refrac0 = Spheric(Ad0, Ad1, Bd1)

   G_old   = G
FTRACE_END  ("MSIS_Expand:1")

End If


!----------------------------------------------------------
! 2. CHEBYSHEV POLYNOMIALS EXPANSION
!----------------------------------------------------------

FTRACE_BEGIN("MSIS_Expand:2")

!--- 2.1. Change of variables

If (G%H < 0d0) then
   z = Tanh(G%H)/h0
Else
   z = G%H/h0
End If

x = 1 - 2*Exp(-z)


!--- 2.2. Clenshaw's recurrence formula (Chebychev polynomials).

dk1 = 0d0
dk2 = 0d0

!NEC$ unroll(klim)
Do k=klim,1,-1
   dk  = 2*x*dk1 - dk2 + c(k)
   dk2 = dk1
   dk1 = dk
End Do

fx = x*dk1 - dk2 + 0.5_Double*c(0)


!--- 2.3. Evaluate the refractivity scale height.

Hs = fx*Exp(-(z/zref)**2) + (alfa*z + beta)


!--- 2.4. Evaluate the refractivity.

MSIS_Expand = refrac0*exp(-G%H/Hs)

FTRACE_END  ("MSIS_Expand:2")

End Function MSIS_Expand



!==========================================================
Real(Double) Function Spheric &
  (A0, A1, B1)   ! <-- Coefficients
                 ! --> Spherical harmonic
!
! Calculation of spherical harmonics.
!----------------------------------------------------------
! Method:
!   Clenshaw's recurrence formula.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!         |    Nov 1996 | F77 version,
!         |             | S. Syndergaard (DMI).
!   1.0   | 28 Oct 1998 | Fortran 90 version.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   A0(0:nlim)
!
Real(Double), Intent(In) :: &
   A1(nlim)
!
Real(Double), Intent(In) :: &
   B1(nlim)
!
! Function result:
!
!Real(Double) :: &
!  Spheric  ! Spherical harmonic
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
    dn10, dn20, dn11, dn21, &
    dn0,  dn1              ! Variables for recurrence
Integer      :: n          ! Recurrence index
!----------------------------------------------------------
! Global variables used:
!
! 1. Constants for Clenshaw's recurrence formula.
! 2. Sin/Cos(Lat/Lon).
!----------------------------------------------------------


dn10 = 0d0
dn20 = 0d0
dn11 = 0d0
dn21 = 0d0

Do n=nlim,1,-1
   dn0  = n10x(n)*dn10 - n20(n)*dn20 + A0(n)
   dn1  = n11x(n)*dn11 - n21(n)*dn21 + A1(n)*CosPhi + B1(n)*SinPhi
   dn20 = dn10
   dn21 = dn11
   dn10 = dn0
   dn11 = dn1
End Do

Spheric = CosTheta*dn10 + SinTheta*dn11 - 0.5_Double*dn20 + A0(0)


End Function Spheric



!==========================================================
Subroutine MSIS_Geop &
  (G,    & ! <-- Geodetic coordinates
   P,    & ! <-- Pressure levels [mb]
   H,    & ! --> Geopotential heigths [gpkm]
   Z)      ! --> Altitudes [km]
!
! Calculation of geopotential heights and altitudes of
! given pressure levels for MSIS.
!----------------------------------------------------------
! Method:
!   Integration of static equation, inversion of
!   dependence p(z), calculation of geopotential heights.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Jun 1999 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
   Geodetic,          &
! Imported Routines:
   GCLat_from_GDLat,  &
   Gravity,           &
   Geop_from_Alt
!
Use Earth, only: &
! Imported Parameters:
   Rd,     &
   C1
!
Use Interpolation, only: &
! Imported Routines:
   Linear
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Geodetic), Intent(In) :: &
   G        ! Geodetic coordinates
            ! Only latitude and longitude are used
!
Real(Double), Intent(In) :: &
   P(:)     ! Pressure levels [mb]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   H(:)     ! Geopotential heights [gpkm]
!
Real(Double), Intent(Out) :: &
   Z(:)     ! Altitudes [km]
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   DZ   = 1.5_Double,   &    ! Integration step
   Zmax = 300.0_Double       ! Upper height limit
Integer, Parameter :: &
   KI = Zmax/DZ + 0.5_Double ! Number of integration steps
!
! Local Scalars:
!
Real(Double)   :: GCLat      ! Geocentric latitude
Type(Geodetic) :: Gi         ! Integration variable
Real(Double)   :: ZC         ! Characteristic scale of N [km]
Real(Double)   :: N0, N1, N2 ! Values of N in integration step
Integer        :: i          ! Integration step index
Integer        :: KP         ! Number of pressure levels
!
! Local Arrays:
!
!Real(Double), Allocatable, Dimension(:) :: &
Real(Double), Dimension(0:KI) :: &
   ZH,       &  ! Hires grid of heights
   PH           ! Hires grid of pressure levels
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------


!--- 1.1. Calculation of geocentric latitude

GCLat = GCLat_from_GDLat(G%Phi)


!--- 1.2. Array allocation

!KI = Ceiling(Zmax/DZ)
!Allocate(ZH(0:KI), PH(0:KI))


!----------------------------------------------------------
! 2. INTEGRATION OF PRESSURE
!----------------------------------------------------------


!--- 2.1. Initial conditions at upper boundary

FTRACE_BEGIN("MSIS_geop:2.1")

Gi    = G
GI%H  = Zmax + DZ
N0    = MSIS_Expand(Gi)*Gravity(Gi%H, GCLat)
ZH(0) = Zmax
GI%H  = ZH(0)
N2    = MSIS_Expand(Gi)*Gravity(Gi%H, GCLat)
ZC    = DZ/(Log(N2)- Log(N0))
PH(0) = N2*1d3*ZC/(Rd*C1)

FTRACE_END  ("MSIS_geop:2.1")

!--- 2.2. Numerical integration

FTRACE_BEGIN("MSIS_geop:2.2")

Do i=1,KI
   ZH(i) = ZH(i-1) - DZ
   N0    = N2
   Gi%H  = (ZH(i) + ZH(i-1))/2
   N1    = MSIS_Expand(Gi)*Gravity(Gi%H, GCLat)
   Gi%H  = ZH(i)
   N2    = MSIS_Expand(Gi)*Gravity(Gi%H, GCLat)
   PH(i) = PH(i-1) + (N0 + 4*N1 + N2)*1d3*DZ/(6*Rd*C1)
End Do

FTRACE_END  ("MSIS_geop:2.2")

!--- 2.3. Logarithmic scale of pressure

PH(:) = Log(PH(:))


!----------------------------------------------------------
! 3. ALTITUDES AND GEOPOTENTIAL HEIGHTS OF PRESSURE LEVELS
!----------------------------------------------------------

FTRACE_BEGIN("MSIS_geop:2.3")

KP = Size(P)

Do i=1,KP
   Call Linear  &
     (PH,        & ! <-- Argument grid
      ZH,        & ! <-- Gridded function
      Log(P(i)), & ! <-- Interpolation point
      Z(i))        ! --> Interpolated function value
   H(i) = Geop_from_Alt(Z(i), GCLat)
End Do

FTRACE_END  ("MSIS_geop:2.3")

!Deallocate(ZH, PH)

End Subroutine MSIS_Geop



!==========================================================
Subroutine MSIS_Pressure_Levels &
  (P0,   & ! <-- Basic pressure level
   P)      ! --> MSIS pressure levels
!
! Calculation of MSIS pressure level grid.
!----------------------------------------------------------
! Method:
!   Pressure levels in logarithmic scale.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 21 Jun 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   P0       ! Basic pressure levels
            ! Output pressure levels are: p(i)=p0*10**(-i/3), i=1,...,N
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   P(:)     ! MSIS pressure levels
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N  ! Number of pressure levels
Integer      :: i  ! Array index
Real(Double) :: S  ! Scaling factor
!----------------------------------------------------------


N    = Size(P)
S    = 10.0_Double**(-LgScale)

P(N) = P0*S

Do i=N-1,1,-1
   P(i) = P(i+1)*S
End Do


End Subroutine MSIS_Pressure_Levels



!==========================================================
Function MSIS_Num_Levels &
  (P0,     & ! <-- Basic pressure level [mb]
   PU)     & ! <-- Upper pressure level [mb]
Result(N)    ! --> Number of MSIS pressure levels
!
! Calculation of number of logarithmic-scaled pressure
! levels.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 21 Jun 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   P0       ! Basic pressure level
            ! Pressure levels are: p(i)=p0*10**(-i/3), i=1,...,N
!
Real(Double), Intent(In) :: &
   PU       ! Upper pressure level
            ! p(N)>=PU
!
! Function result:
!
Integer :: &
   N        ! Number of MSIS pressure levels
!----------------------------------------------------------


N = Ceiling((Log(P0/PU)/Log(10.0_Double))/LgScale)


End Function MSIS_Num_Levels



End Module MSIS
