!
!+ GNSS Radio occultation observation operator: Adjoint for ECHAM_fields
!
MODULE ECHAM_fields_adj
!
! Description:
!   GNSS Radio occultation bending angle observation operator.
!   Adjoint for ECHAM_fields.
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
!  optimize for SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_11        2010/06/16 Harald Anlauf
!  N_from_TPQ_adj: use 3-term refractivity formula
! V1_13        2011/11/01 Harald Anlauf
!  Make FTRACE_REGIONs depending on macro DISABLE_FTRACE_REGION
! V1_14        2011/11/08 Harald Anlauf
!  Minor cleanup
! V1_22        2013-02-13 Harald Anlauf
!  add vertical coordinate type, geopotential height
! V1_26        2013/06/27 Andreas Rhodin
!  remove "stat=" option in deallocate statements
! V1_27        2013-11-08 Harald Anlauf
!  for ICON: horizontal interpolation, remove ak,bk, waterload
! V1_42        2015-06-08 Harald Anlauf
!  horint_mode
! V1_43        2015-08-19 Harald Anlauf
!  implement refractivity expression from Aparicio & Laroche (2011)
! V1_46        2016-02-05 Andreas Rhodin
!  base decisions on new flag 'vct', not 'ivctype'
! V1_50        2017-01-09 Andreas Rhodin
!  adapt COSMO-MEC to ICON-LAM: return the 6 surrounding grid-points
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
! Module ECHAM_fields_adj
!
! Adjoint for ECHAM_fields.
!----------------------------------------------------------
! (C) Copyright 1999-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   8.0   | 21 Jun 1999 | Basic non-adjoint version.
!   1.0   | 14 Apr 2000 | Adjoint version.
!   2.0   | 07 Apr 2002 | Icosahedral grid included.
!   3.0   | 16 Apr 2002 | Field status.
!   3.1   | 20 Apr 2002 | Allocate_Profile_adj.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, WorkPr
!
Use ECHAM_fields, only: &
! Imported Type Definitions:
    Profile,   &
    Matrix,    &
! Imported Parameters:
    NG1,       & ! Number of points for lat/lon interpolation:
    fst_Null,        &
    fst_Allocated,   &
    fst_Initialized, &
! Imported Scalar Variables:
    NM,        & ! Number of MSIS half levels.
! Imported Array Variables:
!   A, B,      & ! Vertical coordinates.
    gf,        & ! global field data type
!   Hsur,      & ! Surface geopotential [gpm]
!   Psur,      & ! Surface pressure [Pa]
!   T,         & ! Temperature [K]
!   Q,         & ! Specific humidity [kg/kg]
!   XLon,      & ! Longitude grid [deg]
!   XLat,      & ! Latitude grid [deg]
!   GCLat,     & ! Geocentric latitudes [rad]
!   Zsur,      & ! Surface altitude [km]
    Z,         & ! Altitudes [km]
    LnN,       & ! Ln of refractive index
    D2N,       & ! Interpolation coefficients of Ln(N)
    PStat,     & ! Dynamical field status
! Imported Routines:
    Allocate_Profile
!
Use mo_wmo_tables, only: &
! Imported Parameters:
    WMO6_LATLON,         &
    WMO6_GAUSSIAN,       &
    DWD6_ICOSAHEDRON,    &
    DWD6_ICON
!
Use mo_exception, only: &
    finish

Use mo_atm_grid, only: &
    VCT_P_HYB  ! hybrid pressure vertical coordinate flag
!----------------------------------------------------------
Implicit None
Private
Public :: ECHAM_NGradN_adj, Interpolate_Refractivity_adj, &
          Allocate_Profile_adj, echam_init_adj, echam_cleanup_adj
!
! Public Arrays:
!
Integer, Allocatable, Public :: &
   PAStat(:,:,:)     ! Dynamical field status

!==============================================================================
#if   defined (_CRAYFTN) \
  || (defined (__GFORTRAN__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 406)  \
                             && ((__GNUC__ * 100 + __GNUC_MINOR__) <  800)) \
  || (defined (__INTEL_COMPILER) && (__INTEL_COMPILER >= 1500)) \
  || (defined (__NEC__) && (__NEC_VERSION__ >= 30004)) \
  || (defined (__PGI) && (__PGIC__ >= 15) && (__PGIC__ != 17))
#define _CONTIGUOUS         ,Contiguous
#define _POINTER    ,Pointer,Contiguous
#else
#define _CONTIGUOUS
#define _POINTER    ,Pointer
#endif
!==============================================================================
!
! Private Arrays: (originally public)
!
Type(Profile) _POINTER, Private :: &
   LnN_T(:,:,:),   & ! d(LnN(i))/d(T(i))
   LnN_Q(:,:,:),   & ! d(LnN(i))/d(Q(i))
   LnN_P(:,:,:),   & ! d(LnN(i))/d(Psur)
   LnNM_T(:,:,:),  & ! d(LnN_MSIS)/d(T(i))
   LnNM_Q(:,:,:)     ! d(LnN_MSIS)/d(Q(i))
Real(Double) _POINTER, Private  :: &
   LnNM_P(:,:,:)     ! d(LnN_MSIS)/d(Psur)
Type(Matrix) _POINTER, Private  :: &
   Z_T(:,:,:),     & ! d(Z(i))/d(T(j))
   Z_Q(:,:,:)        ! d(Z(i))/d(Q(j))
Type(Profile) _POINTER, Private :: &
   Z_P(:,:,:)        ! d(Z(i))/d(Psur)
Type(Matrix) _POINTER, Private  :: &
   D2N_Z(:,:,:),   & ! d(D2N(i))/d(Z(j))
   D2N_N(:,:,:)      ! d(D2N(i))/d(LnN(j))
!----------------------------------------------------------
!
! Private Scalars:
!
Logical, Private :: &
   Undefined_Pointers = .TRUE. ! Pointer status indicator.
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine ECHAM_cleanup_adj
!----------------------------------------------------------
! 2. CLEARING DATA
!----------------------------------------------------------

! This part was moved from ECHAM_init_adj to ECHAM_cleanup_adj

Integer   :: I1       ! 1st grid index
Integer   :: I2       ! 2nd grid index
Integer   :: ID       ! Diamond index

If (Undefined_Pointers) then
   Nullify(LnN_T)
   Nullify(LnN_Q)
   Nullify(LnN_P)
   Nullify(LnNM_T)
   Nullify(LnNM_Q)
   Nullify(LnNM_P)
   Nullify(Z_T)
   Nullify(Z_Q)
   Nullify(Z_P)
   Nullify(D2N_Z)
   Nullify(D2N_N)
   Undefined_Pointers = .FALSE.
Else
   If (allocated (PAStat)) then
FTRACE_BEGIN("ECHAM_cleanup_adj:Deallocate")
      Do ID = LBound(PAStat,3),UBound(PAStat,3)
         Do I2 = LBound(PAStat,2),UBound(PAStat,2)
            If (any (PAStat(:,I2,ID) /= fst_Null)) then
               Do I1 = LBound(PAStat,1),UBound(PAStat,1)
                  if (PAStat(I1, I2, ID) /= fst_Null) then
                     Deallocate(LnN_T (I1, I2, ID)%P)
                     Deallocate(LnN_Q (I1, I2, ID)%P)
                     Deallocate(LnN_P (I1, I2, ID)%P)
                     Deallocate(LnNM_T(I1, I2, ID)%P)
                     Deallocate(LnNM_Q(I1, I2, ID)%P)
                     Deallocate(Z_T   (I1, I2, ID)%M)
                     Deallocate(Z_Q   (I1, I2, ID)%M)
                     Deallocate(Z_P   (I1, I2, ID)%P)
                     Deallocate(D2N_Z (I1, I2, ID)%M)
                     Deallocate(D2N_N (I1, I2, ID)%M)
                     PAStat(I1, I2, ID) = fst_Null
                  end if
               End Do
            End If
         End Do
      End Do
      Deallocate (LnN_T )
      Deallocate (LnN_Q )
      Deallocate (LnN_P )
      Deallocate (LnNM_T)
      Deallocate (LnNM_Q)
      Deallocate (LnNM_P)
      Deallocate (Z_T   )
      Deallocate (Z_Q   )
      Deallocate (Z_P   )
      Deallocate (D2N_Z )
      Deallocate (D2N_N )
      Deallocate (PAStat)
FTRACE_END  ("ECHAM_cleanup_adj:Deallocate")
   End If
End If

end Subroutine ECHAM_cleanup_adj
!==========================================================
Subroutine ECHAM_Init_adj
!
! Initialization of ECHAM_fields_adj module.
!----------------------------------------------------------
! Method:
!   Clearing module global data.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Apr 2000 | Original code.
!   1.1   | 16 Apr 2002 | Field status.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Local Scalars:
!
Integer   :: N1s, N1e ! 1st grid index limits
Integer   :: N2       ! 2nd grid dimension
Integer   :: ND       ! Number of diamonds
!Integer   :: I1, I2   ! Grid indices
!Integer   :: ID       ! Diamond index
!----------------------------------------------------------
! Global variables used:
!
!  LnN_T(:,:)      ! d(LnN(i))/d(T(i))
!  LnN_Q(:,:)      ! d(LnN(i))/d(Q(i))
!  LnN_P(:,:)      ! d(LnN(i))/d(Psur)
!  LnNM_T(:,:)     ! d(LnN_MSIS)/d(T(i))
!  LnNM_Q(:,:)     ! d(LnN_MSIS)/d(Q(i))
!  LnNM_P(:,:)     ! d(LnN_MSIS)/d(Psur)
!  Z_T(:,:)        ! d(Z(i))/d(T(j))
!  Z_Q(:,:)        ! d(Z(i))/d(Q(j))
!  Z_P(:,:)        ! d(Z(i))/d(Psur)
!  D2N_Z(:,:)      ! d(D2N(i))/d(Z(j))
!  D2N_N(:,:)      ! d(D2N(i))/d(LnN(j))
!----------------------------------------------------------


!----------------------------------------------------------
! 1. DETERMINATION OF FIELD DIMENSIONS
!----------------------------------------------------------

N1S  = LBound(Z, 1)
N1E  = UBound(Z, 1)
N2   = Size(Z, 2)
ND   = Size(Z, 3)


!----------------------------------------------------------
! 2. CLEARING DATA
!----------------------------------------------------------

! This part was moved from ECHAM_init_adj to ECHAM_cleanup_adj

Call ECHAM_cleanup_adj

!----------------------------------------------------------
! 2. MEMORY ALLOCATION FOR SURFACE ARRAYS
!----------------------------------------------------------


!--- 2.1. Memory allocation

Allocate(PAStat(N1S:N1E, N2, ND))
PAStat(:,:,:) = fst_Null

FTRACE_BEGIN("ECHAM_Init_adj:Allocate")
Allocate(LnN_T (N1S:N1E, N2, ND),   &
         LnN_Q (N1S:N1E, N2, ND),   &
         LnN_P (N1S:N1E, N2, ND),   &
         LnNM_T(N1S:N1E, N2, ND),   &
         LnNM_Q(N1S:N1E, N2, ND),   &
         LnNM_P(N1S:N1E, N2, ND),   &
         Z_T   (N1S:N1E, N2, ND),   &
         Z_Q   (N1S:N1E, N2, ND),   &
         Z_P   (N1S:N1E, N2, ND),   &
         D2N_Z (N1S:N1E, N2, ND),   &
         D2N_N (N1S:N1E, N2, ND)    )
FTRACE_END  ("ECHAM_Init_adj:Allocate")


!--- 2.2. Nullifying vertical profiles and matrices

!Do I1=N1S,N1E
!   Do I2=1,N2
!      Do ID=1,ND
!         Nullify(LnN_T (I1, I2, ID)%P)
!         Nullify(LnN_Q (I1, I2, ID)%P)
!         Nullify(LnN_P (I1, I2, ID)%P)
!         Nullify(LnNM_T(I1, I2, ID)%P)
!         Nullify(LnNM_Q(I1, I2, ID)%P)
!         Nullify(Z_T   (I1, I2, ID)%M)
!         Nullify(Z_Q   (I1, I2, ID)%M)
!         Nullify(Z_P   (I1, I2, ID)%P)
!         Nullify(D2N_Z (I1, I2, ID)%M)
!         Nullify(D2N_N (I1, I2, ID)%M)
!      End Do
!   End Do
!End Do


!!! TEST
!Write (*,'(A)')    'ECHAM_field_adj status:'
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnN_T  ', Associated(LnN_T),  Shape(LnN_T)
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnN_Q  ', Associated(LnN_Q),  Shape(LnN_Q)
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnN_P  ', Associated(LnN_P),  Shape(LnN_P)
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnNM_T ', Associated(LnNM_T), Shape(LnNM_T)
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnNM_Q ', Associated(LnNM_Q), Shape(LnNM_Q)
!Write (*,'(A,L5,3(I4:'' x ''))') 'LnNM_P ', Associated(LnNM_P), Shape(LnNM_P)
!Write (*,'(A,L5,3(I4:'' x ''))') 'Z_T    ', Associated(Z_T),    Shape(Z_T)
!Write (*,'(A,L5,3(I4:'' x ''))') 'Z_Q    ', Associated(Z_Q),    Shape(Z_Q)
!Write (*,'(A,L5,3(I4:'' x ''))') 'Z_P    ', Associated(Z_P),    Shape(Z_P)
!Write (*,'(A,L5,3(I4:'' x ''))') 'D2N_Z  ', Associated(D2N_Z),  Shape(D2N_Z)
!Write (*,'(A,L5,3(I4:'' x ''))') 'D2N_N  ', Associated(D2N_N),  Shape(D2N_N)

End Subroutine ECHAM_Init_adj



!==========================================================
Elemental &
Subroutine N_from_TPQ_adj &
  (T,     & ! <-- Temperature [K]
   P,     & ! <-- Pressure [mb]
   Q,     & ! <-- Specific humidity [kg/kg]
   N,     & ! --> Refractivity [dimensionless]
   N_T,   & ! --> d(N)/d(T)
   N_P,   & ! --> d(N)/d(P)
   N_Q)     ! --> d(N)/d(Q)
!
! Calculation of refractivity from temperature, pressure,
! and humidity: Adjoint version.
!----------------------------------------------------------
! Method:
!
!                 P            P Q                  P Q
!   N(T,P,Q) = C1--- + C2----------------- + C3--------------
!                 T       T**2 (aq + bq Q)      T (aq + bq Q)
!
!   (Bean and Datton, Radiometeorology)
!----------------------------------------------------------
! (C) Copyright 1998-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 30 Sep 1998 | Basic non-adjoint version.
!   1.0   | 10 Apr 2000 | Adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Parameters:
    Rd, Rnu, C1, C2, C3
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
! Output arguments:
!
Real(Double), Intent(Out) :: &
   N    ! Refractivity [dimensionless]
!
Real(Double), Intent(Out) :: &
   N_T  ! d(N)/d(T) [K^-1]
!
Real(Double), Intent(Out) :: &
   N_P  ! d(N)/d(P) [mb^-1]
!
Real(Double), Intent(Out) :: &
   N_Q  ! d(N)/d(Q) [dimensionless]
!----------------------------------------------------------
! Local Parameters:
!
Real(Double), Parameter :: &
   aq = Rd/Rnu,   &        ! Pw = P*q/(aq + bq*q)
   bq = 1.0_Double - aq
!----------------------------------------------------------

#ifdef _SMITH_WEINTRAUB_TWOTERM

N = (C1 + C2*Q/(T*(aq + bq*Q)))*P/T

N_T = -(C1 + 2*C2*Q/(T*(aq + bq*Q)))*P/T**2

N_P = N/P

N_Q = C2*P*aq/(T*(aq + bq*Q))**2

#else

Real(Double) :: Pw_r  ! Pw/P

Pw_r = Q/(aq + bq*Q)

N    =  (C1 +   C2*Pw_r/T + C3*Pw_r)*P/T

N_T  = -(C1 + 2*C2*Pw_r/T + C3*Pw_r)*P/T**2

N_P  = N/P

N_Q  = P*(C2 + C3*T)*aq/(T*(aq + bq*Q))**2

#endif

End Subroutine N_from_TPQ_adj



!==========================================================
Subroutine Allocate_Profile_adj &
  (I1,       & ! <-- 1st grid index
   I2,       & ! <-- 2nd grid index
   ID)         ! <-- Diamond index

!
! Allocation of profiles for a grid point.
!----------------------------------------------------------
! (C) Copyright 2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 20 Apr 2002 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   I1            ! 1st grid index
!
Integer, Intent(In) :: &
   I2            ! 2nd grid index
!
Integer, Intent(In) :: &
   ID            ! Diamond index
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev    ! Number of levels
!----------------------------------------------------------


!----------------------------------------------------------
! 1. MEMORY ALLOCATION
!----------------------------------------------------------

NLev = gf% nz


Allocate(LnN_T (I1, I2, ID)%P(1:NLev),                &
         LnN_Q (I1, I2, ID)%P(1:NLev),                &
         LnN_P (I1, I2, ID)%P(1:NLev),                &
         LnNM_T(I1, I2, ID)%P(1:NLev),                &
         LnNM_Q(I1, I2, ID)%P(1:NLev),                &
         Z_T   (I1, I2, ID)%M(1:NLev,1:NLev),         &
         Z_Q   (I1, I2, ID)%M(1:NLev,1:NLev),         &
         Z_P   (I1, I2, ID)%P(1:NLev),                &
         D2N_Z (I1, I2, ID)%M(-NM+2:NLev,-NM+2:NLev), &
         D2N_N (I1, I2, ID)%M(-NM+2:NLev,-NM+2:NLev))

PAStat(I1, I2, ID) = fst_Allocated


End Subroutine Allocate_Profile_adj



!==========================================================
Subroutine Make_Refractivity_adj &
  (Hsur,    & ! <-- Surface geopotential [gpm]
   Gundu,   & ! <-- Geoid undulation     [m]
   Psur,    & ! <-- Surface pressure [Pa]
   Ph,      & ! <-- Pressure at half levels [Pa]
   T,       & ! <-- Temperature [K]
   Q,       & ! <-- Specific humidity [kg/kg]
   Qx,      & ! <-- Water load [kg/kg]
   G,       & ! <-- Geodetic coordinates
   GCLat,   & ! <-- Geocentric latitude [rad]
   LnN,     & ! --> Ln of refractive index
   LnN_T,   & ! --> d(LnN(i))/d(T(i))
   LnN_Q,   & ! --> d(LnN(i))/d(Q(i))
   LnN_P,   & ! --> d(LnN(i))/d(Psur)
   LnNM_T,  & ! --> d(LnN_MSIS)/d(T(i))
   LnNM_Q,  & ! --> d(LnN_MSIS)/d(Q(i))
   LnNM_P,  & ! --> d(LnN_MSIS)/d(Psur)
   Zsur,    & ! --> Surface altitude [km]
   Z,       & ! --> Altitudes of model levels [km]
   Z_T,     & ! --> d(Z(i))/d(T(j))
   Z_Q,     & ! --> d(Z(i))/d(Q(j))
   Z_P,     & ! --> d(Z(i))/d(Psur)
   H,       & ! ~~> Geopotential height at full levels
   P        ) ! ~~> Pressure at full levels
!
! Calculation of refractivity profile and geometrical
! altitudes of pressure levels from surface  geopotential,
! surface pressure, and profiles of temperature and humidity:
! Adjoint version.
!----------------------------------------------------------
! Method:
!   Described in:
!   1. Deutsches Klimarechenzentrum, Technical Report No.6,
!      The ECHAM3 Atmospheric General Circulation Model,
!      Revision 2, Hamburg, 1993.
!   2. M. E. Gorbunov, S. V. Sokolovsky, L. Bengtsson, Space
!      refractive tomography of the atmosphere: modeling of
!      the direct and inverse problems, Max-Planck Institut
!      fuer Meteorologie, Report No. 210, Hamburg, 1996.
!----------------------------------------------------------
! (C) Copyright 1999-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   3.0   | 01 Jun 1999 | Basic non-adjoint version.
!   1.0   | 11 Apr 2000 | Adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Parameters:
    g_ave,     &
! Imported Routines:
    Alt_from_Geop,      &
    Alt_from_Geop_adj
!
Use Occ_Meteoprofiles, only: &
! Imported Parameters:
    Rd, Eps,               &
! Imported Routines:
!   N_from_TPQ
    N_from_TPQ_AL2011_adj
!
Use MSIS, only: &
! Imported Routines:
    MSIS_Pressure_Levels,  &
    MSIS_Geop,             &
    MSIS_Refractivity
!
Use CIPM_2007, only: &
    Zeta_adj
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(WorkPr), Intent(In)    :: &
   Hsur    ! Surface geopotential [gpm]
!
Real(WorkPr), Intent(In)    :: &
   Gundu   ! Geoid undulation [m]
!
Real(WorkPr), Intent(In)    :: &
   Psur    ! Surface pressure [Pa]
!
Real(WorkPr), Intent(In)   :: &
   Ph(0:)  ! Pressure at half levels [Pa]
!
Real(WorkPr), Intent(In)    :: &
   T(1:)   ! Temperature profile [K]
!
Real(WorkPr), Intent(In)    :: &
   Q(1:)   ! Specific humidity [kg/kg]
!
Real(WorkPr), Intent(In)   :: &
   Qx(1:)  ! Water load [kg/kg]
!
Type(Geodetic), Intent(In)  :: &
   G       ! Geodetic coordinates
!
Real(Double), Intent(In)    :: &
   GCLat   ! Geocentric latitude [rad]
!
! Output arguments:
!
Real(Double), Intent(Out)   :: &
   LnN(-NM+2:) ! Ln of refractive index
!
Real(Double), Intent(Out)   :: &
   LnN_T(1:)   ! d(LnN(i))/d(T(i)) i=1..NLev
!
Real(Double), Intent(Out)   :: &
   LnN_Q(1:)   ! d(LnN(i))/d(Q(i)) i=1..NLev
!
Real(Double), Intent(Out)   :: &
   LnN_P(1:)   ! d(LnN(i))/d(Psur) i=1..NLev
!
Real(Double), Intent(Out)   :: &
   LnNM_T(1:)  ! d(LnN_MSIS)/d(T(i))
!
Real(Double), Intent(Out)   :: &
   LnNM_Q(1:)  ! d(LnN_MSIS)/d(Q(i))
!
Real(Double), Intent(Out)   :: &
   LnNM_P      ! d(LnN_MSIS)/d(Psur)
!
Real(Double), Intent(Out)   :: &
   Zsur        ! Surface altitude [km]
!
Real(Double), Intent(Out)   :: &
   Z(-NM+2:)   ! Altitudes of model levels [km]
!
Real(Double), Intent(Out)   :: &
   Z_T(1:,1:)  ! d(Z(i))/d(T(j)) i,j = 1..NLev
!
Real(Double), Intent(Out)   :: &
   Z_Q(1:,1:)  ! d(Z(i))/d(Q(j)) i,j = 1..NLev
!
Real(Double), Intent(Out)   :: &
   Z_P(1:)     ! d(Z(i))/d(Psur)
!
Real(Double), Intent(Out), Optional  :: &
   H(1:)       ! Geopotential height at full levels
!
Real(Double), Intent(Out), Optional  :: &
   P(1:)       ! Pressure at full levels
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: NLev    ! Number of levels
Integer         :: i       ! Level index
Real(Double)    :: alpha   ! Log(p(i+1)/p(i))
Real(Double)    :: alpha_P ! d(alpha)/d(Psur)
Type(Geodetic)  :: GM      ! Point in vertical profile
Real(Double)    :: DLn     ! Correction of Ln(P_MSIS)
Real(Double)    :: N       ! Refractivity [dimensionless]
Real(Double)    :: N_T     ! d(N)/d(T)
Real(Double)    :: N_P     ! d(N)/d(P)
Real(Double)    :: N_Q     ! d(N)/d(Q)
Real(Double)    :: LnScale ! Log of refractivity scaling factor
!
!
! Local Arrays:
!
! --- Half levels
!
Real(Double), Dimension(-NM+1:Size(T)) :: &
   Phalf,          & ! Half level pressure [Pa]
   Hhalf             ! Half level geopotential [gpm]
!
! --- Full levels
!
Real(Double), Dimension(-NM+2:Size(T)) :: &
   Pfull,          & ! Full level pressure [Pa]
   Hfull,          & ! Full level geopotential [gpm]
   Tvirt             ! Virtual temperature [K]
!
! --- Work arrays
!
Real(Double), Dimension(-NM+2:1) :: &
   PM,             & ! MSIS level pressure [Pa]
   HM,             & ! MSIS level geopotential [gpkm]
   ZM                ! MSIS level altitude [km]
Real(Double), Allocatable :: &
   Hhalf_T(:,:),   & ! d(Hhalf(i))/d(T(j)) [gpm/K]
   Hfull_T(:,:),   & ! d(Hfull(i))/d(T(j)) [gpm/K]
   Hhalf_Q(:,:),   & ! d(Hhalf(i))/d(Q(j)) [gpm/K]
   Hfull_Q(:,:),   & ! d(Hfull(i))/d(Q(j)) [gpm/K]
   Hhalf_P(:),     & ! d(Hfull(i))/d(Psur) [gpm/Pa]
   Pfull_P(:),     & ! d(Pfull(i))/d(Psur) [Pa/Pa]
   Hfull_P(:),     & ! d(Hfull(i))/d(Psur) [gpm/Pa]
   Z_GPH(:)          ! d(Z(i))/d(Hfull(i)) [km/gpkm]
Real(Double), Dimension(Size(T)) :: &
   Tv_T,           & ! d(Tvirt(i))/d(T(i))
   Tv_Q              ! d(Tvirt(i))/d(Q(i))
!
Real(Double), Dimension(Size(T),3) :: &
   Zeta,          & ! Compressibility correction function
   Zeta_P,        & ! and derivatives
   Zeta_T,        & !
   Zeta_Q
!
Real(Double), Dimension(size(T)) :: &
   DZhalf,         & ! Correction from compressibility
   DZfull,         & ! for half and full levels
   DZhalf_T,       & ! and derivatives
   DZfull_T,       & ! w.r.t.
   DZhalf_Q,       & ! full-level variables
   DZfull_Q          !
!----------------------------------------------------------
! Global variables used:
!
!   A(:), B(:)   ! Vertical coordinates (ECHAM/GME only).
!   NM           ! Number of MSIS half levels.
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

NLev = Size(T)

Allocate(Hhalf_T(0:NLev,NLev))
Allocate(Hfull_T(NLev,NLev))
Allocate(Hhalf_Q(0:NLev,NLev))
Allocate(Hfull_Q(NLev,NLev))
Allocate(Hhalf_P(0:NLev))
Allocate(Pfull_P(1:NLev))
Allocate(Hfull_P(1:NLev))
Allocate(Z_GPH(NLev))


!----------------------------------------------------------
! 1. CALCULATION OF ECHAM HALF AND FULL LEVELS
!----------------------------------------------------------

!--- 1.1. Calculation of ECHAM half levels

!Phalf(0:NLev) = gf% A(0:NLev) + gf% B(0:NLev)*Psur
!NEC$ SHORTLOOP
Phalf(0:NLev) = Ph(0:NLev)
!NEC$ SHORTLOOP
Tvirt(1:NLev) = T(1:NLev)*(1 + Eps*Q(1:NLev) - Qx(1:NLev))
!NEC$ SHORTLOOP
Tv_T (1:NLev) =           (1 + Eps*Q(1:NLev) - Qx(1:NLev))
!NEC$ SHORTLOOP
Tv_Q (1:NLev) = T(1:NLev)*     Eps

if (gf% vctype == VCT_P_HYB) then
!NEC$ SHORTLOOP
Pfull(1:NLev) =      (Phalf(0:NLev-1) + Phalf(1:NLev)) * 0.5_Double
else
!NEC$ SHORTLOOP
Pfull(1:NLev) = sqrt (Phalf(0:NLev-1) * Phalf(1:NLev))
end if

if (gf% vctype == VCT_P_HYB) then
   Phalf(0) = Phalf(1) / 4      ! ECHAM/GME
end if


Hhalf(NLev) = Hsur/g_ave + Gundu
Zsur        = Alt_from_Geop(1e-3_Double*Hhalf(NLev), GCLat)

!NEC$ SHORTLOOP
DZhalf   = 0._Double
!NEC$ SHORTLOOP
DZfull   = 0._Double
!NEC$ SHORTLOOP
DZhalf_T = 0._Double
!NEC$ SHORTLOOP
DZfull_T = 0._Double
!NEC$ SHORTLOOP
DZhalf_Q = 0._Double
!NEC$ SHORTLOOP
DZfull_Q = 0._Double
Zeta_P   = 0._Double
if (gf% ref_model == 2) then
   Zeta = 0._Double
!NEC$ SHORTLOOP
   do i = 1, NLev
      call Zeta_adj (Phalf(i-1), T(i), Q(i), Zeta(i,1),   Zeta_P(i,1), &
                                             Zeta_T(i,1), Zeta_Q(i,1))
      call Zeta_adj (Phalf(i  ), T(i), Q(i), Zeta(i,2),   Zeta_P(i,2), &
                                             Zeta_T(i,2), Zeta_Q(i,2))
      call Zeta_adj (Pfull(i  ), T(i), Q(i), Zeta(i,3),   Zeta_P(i,3), &
                                             Zeta_T(i,3), Zeta_Q(i,3))
      DZhalf  (i) = Zeta  (i,2) - Zeta  (i,1)
      DZhalf_T(i) = Zeta_T(i,2) - Zeta_T(i,1)
      DZhalf_Q(i) = Zeta_Q(i,2) - Zeta_Q(i,1)
      DZfull  (i) = Zeta  (i,2) - Zeta  (i,3)
      DZfull_T(i) = Zeta_T(i,2) - Zeta_T(i,3)
      DZfull_Q(i) = Zeta_Q(i,2) - Zeta_Q(i,3)
   end do
end if


!NEC$ SHORTLOOP
Do i=NLev-1,0,-1
   alpha    = Log(Phalf(i+1)/Phalf(i))
   Hhalf(i) = Hhalf(i+1) + Rd/g_ave*Tvirt(i+1) * (alpha + DZhalf(i+1))
End Do


!------ 1.1.1. Adjoint: T, Q, and P derivatives

!NEC$ SHORTLOOP
Hhalf_T(NLev,:) = 0
!NEC$ SHORTLOOP
Hhalf_Q(NLev,:) = 0
Hhalf_P(NLev)   = 0

Do i=NLev-1,0,-1

   alpha          = Log(Phalf(i+1)/Phalf(i))
   alpha_P        = gf% B_ad(i+1)/Phalf(i+1) - gf% B_ad(i)/Phalf(i)

!NEC$ SHORTLOOP
   Hhalf_T(i,:)   = Hhalf_T(i+1,:)
   Hhalf_T(i,i+1) = Hhalf_T(i,i+1) +  &
      Rd/g_ave*Tv_T (i+1) * (alpha + DZhalf  (i+1)) + &
      Rd/g_ave*Tvirt(i+1) *          DZhalf_T(i+1)

!NEC$ SHORTLOOP
   Hhalf_Q(i,:)   = Hhalf_Q(i+1,:)
   Hhalf_Q(i,i+1) = Hhalf_Q(i,i+1) +  &
      Rd/g_ave*Tv_Q (i+1) * (alpha + DZhalf  (i+1)) + &
      Rd/g_ave*Tvirt(i+1) *          DZhalf_Q(i+1)

   Hhalf_P(i)     = Hhalf_P(i+1)   +  &
      Rd/g_ave*Tvirt(i+1) *           &
      (alpha_P +   (gf% B_ad(i+1)*Zeta_P(i+1,2) - gf% B_ad(i)*Zeta_P(i+1,1)))

End Do

if (gf% vctype == VCT_P_HYB) then
   Phalf(0) = Ph(0)             ! ECHAM/GME: restore pressure at model top
end if


!NEC$ SHORTLOOP
Do i = 1,NLev
   if (gf% vctype == VCT_P_HYB) then
      Pfull_P(i) = (gf% B_ad(i-1) + gf% B_ad(i))/2
   else
!     Pfull_P(i) = (gf% B_ad(i-1) * sqrt (Phalf(i)/Phalf(i-1)) + &
!                   gf% B_ad(i)   * sqrt (Phalf(i-1)/Phalf(i)) )/2
      Pfull_P(i) = (gf% B_ad(i-1) * Phalf(i) + gf% B_ad(i) * Phalf(i-1)) / &
                   (2 * Pfull(i))
   end if
End Do


!--- 1.2. Calculation of ECHAM full levels

!NEC$ SHORTLOOP
Do i = 1,NLev

!  Pfull(i) = (Phalf(i-1) + Phalf(i))/2

   If (i == 1 .and. gf% vctype == VCT_P_HYB) then
      alpha = Log(2.0_Double)
   Else
      alpha = Log(Phalf(i)/Pfull(i))
   End If

   Hfull(i) = Hhalf(i) + Rd/g_ave*Tvirt(i) * (alpha + DZfull(i))

!------ 1.2.1. Adjoint: T, Q, and P derivatives

   If (i == 1 .and. gf% vctype == VCT_P_HYB) then
      alpha_P = 0
   Else
      alpha_P = gf% B_ad(i)/Phalf(i) - Pfull_P(i)/Pfull(i)
   End If
   Hfull_P(i)   = Hhalf_P(i)   + &
                  Rd/g_ave*Tvirt(i) * (alpha_P                 + &
                                       gf% B_ad(i)*Zeta_P(i,2) - &
                                       Pfull_P(i) *Zeta_P(i,3)   )

!NEC$ SHORTLOOP
   Hfull_T(i,:) = Hhalf_T(i,:)
   Hfull_T(i,i) = Hfull_T(i,i) + &
                  Rd/g_ave*Tv_T (i) * (alpha+DZfull  (i)) + &
                  Rd/g_ave*Tvirt(i) *        DZfull_T(i)

!NEC$ SHORTLOOP
   Hfull_Q(i,:) = Hhalf_Q(i,:)
   Hfull_Q(i,i) = Hfull_Q(i,i) + &
                  Rd/g_ave*Tv_Q (i) * (alpha+DZfull  (i)) + &
                  Rd/g_ave*Tvirt(i) *        DZfull_Q(i)

   Call Alt_from_Geop_adj &
     (1d-3*Hfull(i),   & ! <-- Geopotential height [gpkm]
      GCLat,           & ! <-- Geocentric latitude [rad]
      Z(i),            & ! --> Altitude above reference ellipsoid [km]
      Z_GPH(i))          ! --> d(H)/d(GPH)

   Z_T(i,:) = Z_GPH(i)*1d-3*Hfull_T(i,:)
   Z_Q(i,:) = Z_GPH(i)*1d-3*Hfull_Q(i,:)
   Z_P(i)   = Z_GPH(i)*1d-3*Hfull_P(i)

End Do

!----------------------------------------------------
! "Scaling" of refractivity: N = (n-1) -> ref_scale*N
! Note: the scaling factor cancels in derivatives, as
! d(log(ref_scale*N))/dX = (dN/dX)/N for ref_scale>0.
!----------------------------------------------------
LnScale = 0._Double
if (gf% ref_scale /= 1._Double) LnScale = log (gf% ref_scale)

select case (gf% ref_model)
case (1)

!NEC$ shortloop
!NEC$ ivdep
   Do i = 1,NLev
      Call N_from_TPQ_adj &
         (Real(T(i), Double),    & ! <-- Temperature [K]
          1e-2_Double*Pfull(i),  & ! <-- Pressure [mb]
          Real(Q(i), Double),    & ! <-- Specific humidity [kg/kg]
          N,                     & ! --> Refractivity [dimensionless]
          N_T,                   & ! --> d(N)/d(T) [K^-1]
          N_P,                   & ! --> d(N)/d(P) [mb^-1]
          N_Q)                     ! --> d(N)/d(Q) [dimensionless]

      LnN(i)   = Log(N) + LnScale
      LnN_T(i) = N_T/N
      LnN_Q(i) = N_Q/N
      LnN_P(i) = Pfull_P(i)*1e-2_Double*N_P/N
   End Do

case (2)

!NEC$ SHORTLOOP
   Do i = 1,NLev
      Call N_from_TPQ_AL2011_adj &
         (Real(T(i), Double),    & ! <-- Temperature [K]
          Pfull(i),              & ! <-- Pressure    [Pa]
          Real(Q(i), Double),    & ! <-- Specific humidity [kg/kg]
          N,                     & ! --> Refractivity [dimensionless]
          N_T,                   & ! --> d(N)/d(T) [K^-1]
          N_P,                   & ! --> d(N)/d(P) [Pa^-1]
          N_Q)                     ! --> d(N)/d(Q) [dimensionless]

      LnN(i)   = Log(N) + LnScale
      LnN_T(i) = N_T/N
      LnN_Q(i) = N_Q/N
      LnN_P(i) = Pfull_P(i)*N_P/N
   End Do

end select

!NEC$ SHORTLOOP
if (present(H)) H = Hfull(1:) - Gundu
!NEC$ SHORTLOOP
if (present(P)) P = Pfull(1:)

!----------------------------------------------------------
! 2. CALCULATION OF MSIS HALF AND FULL LEVELS
!----------------------------------------------------------

!--- 2.0. Unit conversion

!NEC$ SHORTLOOP
Phalf(0:NLev) = 1e-2_Double*Phalf(0:NLev) ! Pa  --> mb
!NEC$ SHORTLOOP
Hhalf(0:NLev) = 1e-3_Double*Hhalf(0:NLev) ! gpm --> gpkm
!NEC$ SHORTLOOP
Hfull(1:NLev) = 1e-3_Double*Hfull(1:NLev) ! gpm --> gpkm


!--- 2.1. Calculation of MSIS pressure levels

Call MSIS_Pressure_Levels &
! (1e-2_Double * gf% A(1) , & ! <-- Basic pressure level [mb]
  (1e-2_Double * Ph(1)    , & ! <-- Basic pressure level [mb]
   Phalf(-NM+1:0))            ! --> MSIS pressure levels [mb]

!NEC$ NOVECTOR
Do i=-NM+2,1
   Pfull(i) = (Phalf(i-1) + Phalf(i))/2
End Do


!--- 2.2. Calculation of MSIS full and half(0) level altitudes

!NEC$ NOVECTOR
PM(-NM+2:0) = Pfull(-NM+2:0)
PM(1)       = Phalf(0)

Call MSIS_Geop &
  (G,            & ! <-- Geodetic coordinates
   PM(-NM+2:1),  & ! <-- Pressure levels [mb]
   HM(-NM+2:1),  & ! --> Geopotential heights [gpkm]
   ZM(-NM+2:1))    ! --> Altitudes [km]

!NEC$ NOVECTOR
Hfull(-NM+2:0) = HM(-NM+2:0)
Hhalf(0)       = HM(1)
!NEC$ NOVECTOR
Z(-NM+2:0)     = ZM(-NM+2:0)


!--- 2.3. Calculation of MSIS full level refractivities

GM = G

Do i=-NM+2,1
   GM%H = Z(i)
   Call MSIS_Refractivity &
     (GM,     & ! <-- Geodetic coordinates
      LnN(i))   ! --> Refractivity
   LnN(i) = Log(LnN(i))
End Do


!--- 2.4. ECHAM-MSIS transfer area

Hfull(1) = Hhalf(1) + 1e-3_Double*(Rd/g_ave)*Tvirt(1)*Log(Phalf(1)/Pfull(1))

!------ 2.4.1. Adjoint: T, Q, and P derivatives

!NEC$ SHORTLOOP
Hfull_T(1,:) = Hhalf_T(1,:)
Hfull_T(1,1) = Hfull_T(1,1) + &
      (Rd/g_ave)*(1 + Eps*Q(1) - Qx(1))*Log(Phalf(1)/Pfull(1))
!NEC$ SHORTLOOP
Hfull_Q(1,:) = Hhalf_Q(1,:)
Hfull_Q(1,1) = Hfull_Q(1,1) + &
      (Rd/g_ave)*     Eps*T(1)         *Log(Phalf(1)/Pfull(1))
Hfull_P(1)   = Hhalf_P(1)

Call Alt_from_Geop_adj &
  (Hfull(1),   & ! <-- Geopotential height [gpkm]
   GCLat,      & ! <-- Geocentric latitude [rad]
   Z(1),       & ! --> Altitude above reference ellipsoid [km]
   Z_GPH(1))     ! --> d(H)/d(GPH)

!NEC$ SHORTLOOP
Z_T(1,:) = Z_GPH(1)*1d-3*Hfull_T(1,:)
!NEC$ SHORTLOOP
Z_Q(1,:) = Z_GPH(1)*1d-3*Hfull_Q(1,:)
Z_P(1)   = Z_GPH(1)*1d-3*Hfull_P(1)

Call N_from_TPQ_adj &
  (Real(T(1), Double),    & ! <-- Temperature [K]
   Pfull(1),              & ! <-- Pressure [mb]
   Real(Q(1), Double),    & ! <-- Specific humidity [kg/kg]
   N,                     & ! --> Refractivity [dimensionless]
   N_T,                   & ! --> d(N)/d(T)
   N_P,                   & ! --> d(N)/d(P)
   N_Q)                     ! --> d(N)/d(Q)

LnN(1)   = Log(N)
LnN_T(1) = N_T/N
LnN_Q(1) = N_Q/N


!--- 2.5. Correction of MSIS profile

DLn          = Log(Phalf(1)/Phalf(0)) - &
               1e3_Double*(g_ave/Rd)*(Hhalf(0) - Hhalf(1))/Tvirt(1)
!NEC$ NOVECTOR
LnN(-NM+2:0) = LnN(-NM+2:0) + DLn

!------ 2.5.1. Adjoint: T, Q, and P derivatives

!NEC$ SHORTLOOP
LnNM_T(:) = (g_ave/Rd)*Hhalf_T(1,:)/Tvirt(1)
LnNM_T(1) = LnNM_T(1) +                           &
   1e3_Double*(g_ave/Rd)*(Hhalf(0) - Hhalf(1))*   &
   (1 + Eps*Q(1) - Qx(1))/Tvirt(1)**2

!NEC$ SHORTLOOP
LnNM_Q(:) = (g_ave/Rd)*Hhalf_Q(1,:)/Tvirt(1)
LnNM_Q(1) = LnNM_Q(1) +                           &
   1e3_Double*(g_ave/Rd)*(Hhalf(0) - Hhalf(1))*   &
   Eps*     T(1)         /Tvirt(1)**2

LnNM_P    = (g_ave/Rd)*Hhalf_P(1)/Tvirt(1)


!----------------------------------------------------------
! 3. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(Hhalf_T)
Deallocate(Hfull_T)
Deallocate(Hhalf_Q)
Deallocate(Hfull_Q)
Deallocate(Hhalf_P)
Deallocate(Pfull_P)
Deallocate(Hfull_P)
Deallocate(Z_GPH)


End Subroutine Make_Refractivity_adj



!==========================================================
Subroutine Vertical_Interpolation_LnN_adj &
  (I1,       & ! <-- 1st grid index
   I2,       & ! <-- 2nd grid index
   ID,       & ! <-- Diamond index
   ZP,       & ! <-- Altitude [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   LnNP,     & ! --> Interpolated Ln(N)
   LnNP_T,   & ! --> d(LnNP)/d(T(i))
   LnNP_Q,   & ! --> d(LnNP)/d(Q(i))
   LnNP_P,   & ! --> d(LnNP)/d(Psur)
   LnNZ,     & ! --> Interpolated dLn(N)/dZ
   LnNZ_T,   & ! --> d(LnNZ)/d(T(i))
   LnNZ_Q,   & ! --> d(LnNZ)/d(Q(i))
   LnNZ_P,   & ! --> d(LnNZ)/d(Psur)
   LnNZ2,    & ! --> Interpolated d2Ln(N)/dZ2
   T,        & ! ~~> Interpolated temperature
   Q,        & ! ~~> Interpolated specific humidity
   Gp,       & ! ~~> Interpolated geopotential height
   P,        & ! ~~> Interpolated Pressure
   Ps,       & ! ~~> Surface Pressure at this Point
   Hs,       & ! ~~> Surface geopotential height
   Gu        ) ! ~~> Geoid undulation
!
! Vertical interpolation of refractivity profiles from ECHAM.
!----------------------------------------------------------
! Method:
!   Vertical spline interpolation to given altitude.
!----------------------------------------------------------
! (C) Copyright 1999-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   5.0   | 20 Feb 1999 | Basic version
!         |             | Interpolate_Refractivity.
!   1.0   | 19 Apr 2000 | Adjoint version of
!         |             | vertical interpolation.
!   2.0   | 08 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!   2.1   | 16 Apr 2002 | Field status.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation_adj, only: &
! Imported Routines:
    Init_Spline_adj,  &
    Spline_adj
!
Use Interpolation, only: &
! Imported Routines:
    Init_Spline, &
    Spline
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic, &
! Imported Parameters:
    g_ave

!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Integer, Intent(In) :: &
   I1            ! 1st grid index
!
Integer, Intent(In) :: &
   I2            ! 2nd grid index
!
Integer, Intent(In) :: &
   ID            ! Diamond index
!
Real(Double), Intent(In) :: &
   ZP            ! Altitude [km]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Zmin          ! Minimum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   Zmax          ! Maximum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   LnNP          ! Interpolated Ln(N)
!
Real(Double), Intent(Out) :: &
   LnNP_T(1:)    ! d(LnNP)/d(T(i))
!
Real(Double), Intent(Out) :: &
   LnNP_Q(1:)    ! d(LnNP)/d(Q(i))
!
Real(Double), Intent(Out) :: &
   LnNP_P        ! d(LnNP)/d(Psur)
!
Real(Double), Intent(Out) :: &
   LnNZ          ! Interpolated dLn(N)/dZ
!
Real(Double), Intent(Out) :: &
   LnNZ_T(1:)    ! d(LnNZ)/d(T(i))
!
Real(Double), Intent(Out) :: &
   LnNZ_Q(1:)    ! d(LnNZ)/d(Q(i))
!
Real(Double), Intent(Out) :: &
   LnNZ_P        ! d(LnNZ)/d(Psur)
!
Real(Double), Intent(Out) :: &
   LnNZ2         ! Interpolated d2Ln(N)/dZ2
!
Real(Double), Intent(Out), Optional :: &
T, Q, Gp, P, Ps, Hs, Gu

!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev  ! Number of levels
Type(Geodetic) :: G     ! Geodetic coordinates
Real(Double)   :: PGCL  ! Geocentric latitude of point
Integer        :: i     ! Array index
Integer        :: ix
#if defined (__NEC__)
Integer        :: k     ! Auxiliary level index
#endif
Real(Double)   :: sum_lnNP_0  ! Sum(lnNP_N(-NM+2:0))
Real(Double)   :: sum_lnNZ_0  ! Sum(lnNZ_N(-NM+2:0))
!
! Local Arrays:
!
Real(Double) _POINTER :: &
   ZLL(:),       &
   FLL(:),       &
   DLL(:)          ! Pointers for array cross sections
Real(Double), Dimension(-NM+2:gf% nz) :: &
   LnNP_Z,       & ! d(LnNP)/d(Z(i))
   LnNP_N,       & ! d(LnNP)/d(LnN(i))
   LnNZ_Z,       & ! d(LnNZ)/d(Z(i))
   LnNZ_N          ! d(LnNZ)/d(LnN(i))
Real(Double) :: d2 (gf% nz)
Real(Double) :: Hf (gf% nz)
Real(Double) :: Pf (gf% nz)
!----------------------------------------------------------

ix = gf% i(i1,i2,id)

!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE
!----------------------------------------------------------

NLev = gf% nz


!----------------------------------------------------------
! 2. VERTICAL INTERPOLATION
!----------------------------------------------------------


!--- 2.1. Vertical profile allocation

If (PStat(I1, I2, ID) == fst_Null) then

!FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.1a")
   Call Allocate_Profile &
     (I1,       & ! <-- 1st grid index
      I2,       & ! <-- 2nd grid index
      ID)         ! <-- Diamond index
!FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.1a")

End If


If (PAStat(I1, I2, ID) == fst_Null) then

!FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.1b")
   Call Allocate_Profile_adj &
     (I1,       & ! <-- 1st grid index
      I2,       & ! <-- 2nd grid index
      ID)         ! <-- Diamond index
!FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.1b")

End If


!--- 2.2. Vertical profile intialization

If (PAStat(I1, I2, ID) == fst_Allocated .or. present (Gp)) then

FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.2")
   G    = Geodetic(0.0_Double, gf% s(ix)% XLat, gf% s(ix)% XLon)
   PGCL = gf% s(ix)% GCLat

   Call Make_Refractivity_adj &
     (gf% s   (ix)% Hsur,        & ! <-- Surface geopotential [gpm]
      gf% s   (ix)% Gundu,       & !     Geoid undulation     [m]
      gf% s   (ix)% Psur,        & ! <-- Surface pressure [Pa]
      gf% Ph(:,ix),              & ! <-- Pressure at half levels [Pa]
      gf% T (:,ix),              & ! <-- Temperature [K]
      gf% Q (:,ix),              & ! <-- Specific humidity [kg/kg]
      gf% Qx(:,ix),              & ! <-- Water load [kg/kg]
      G,                         & ! <-- Geodetic coordinates
      PGCL,                      & ! <-- Geocentric latitude [rad]
      LnN   (I1, I2, ID)%P(:),   & ! --> Ln of refractive index
      LnN_T (I1, I2, ID)%P(:),   & ! --> d(LnN(i))/d(T(i))
      LnN_Q (I1, I2, ID)%P(:),   & ! --> d(LnN(i))/d(Q(i))
      LnN_P (I1, I2, ID)%P(:),   & ! --> d(LnN(i))/d(Psur)
      LnNM_T(I1, I2, ID)%P(:),   & ! --> d(LnN_MSIS)/d(T(i))
      LnNM_Q(I1, I2, ID)%P(:),   & ! --> d(LnN_MSIS)/d(Q(i))
      LnNM_P(I1, I2, ID),        & ! --> d(LnN_MSIS)/d(Psur)
      gf% s   (ix)% Zsur,        & ! --> Surface altitude [km]
      Z     (I1, I2, ID)%P(:),   & ! --> Altitudes of model levels [km]
      Z_T   (I1, I2, ID)%M(:,:), & ! --> d(Z(i))/d(T(j))
      Z_Q   (I1, I2, ID)%M(:,:), & ! --> d(Z(i))/d(Q(j))
      Z_P   (I1, I2, ID)%P(:),   & ! --> d(Z(i))/d(Psur)
      Hf,                        & ! ~~> Geopotential height, full model levels
      Pf                         ) ! ~~> Pressure at full model levels

   Call  Init_Spline_adj   &
     (Z  (I1, I2, ID)%P(:),      & ! <-- Argument grid
      LnN(I1, I2, ID)%P(:),      & ! <-- Gridded function
      D2N(I1, I2, ID)%P(:),      & ! --> 2nd derivative of spline
      D2N_Z(I1, I2, ID)%M(:,:),  & ! --> d(d2(i))/d(x(j))
      D2N_N(I1, I2, ID)%M(:,:))    ! --> d(d2(i))/d(f(j))

   PStat (I1, I2, ID) = fst_Initialized
   PAStat(I1, I2, ID) = fst_Initialized
FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.2")

End If


!--- 2.3. Setting pointers to vertical arrays

ZLL => Z  (I1, I2, ID)%P(:)
FLL => LnN(I1, I2, ID)%P(:)
DLL => D2N(I1, I2, ID)%P(:)


!--- 2.4. Vertical spline interpolation

!FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.4a")
!Allocate(LnNP_Z(-NM+2:NLev),    &
!         LnNP_N(-NM+2:NLev),    &
!         LnNZ_Z(-NM+2:NLev),    &
!         LnNZ_N(-NM+2:NLev))
!FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.4a")

FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.4")
Call Spline_adj &
  (ZLL,                  & ! <-- Argument grid
   FLL,                  & ! <-- Gridded function
   DLL,                  & ! <-- 2nd derivative of spline
   D2N_Z(I1, I2, ID)%M,  & ! <-- d(d2(i))/d(x(j))
   D2N_N(I1, I2, ID)%M,  & ! <-- d(d2(i))/d(f(j))
   ZP,                   & ! <-- Interpolation point
   LnNP,                 & ! --> Interpolated function value
   LnNP_Z,               & ! --> d(f_int)/d(x(i))
   LnNP_N,               & ! --> d(f_int)/d(f(i))
   LnNZ,                 & ! --> Interpolated 1st derivative
   LnNZ_Z,               & ! --> d(fd_int)/d(x(i))
   LnNZ_N,               & ! --> d(fd_int)/d(x(i))
   LnNZ2)                  ! --> Interpolated 2nd derivative
FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.4")


FTRACE_BEGIN("Vertical_Interpolation_LnN_adj:2.5_7")
!CDIRR ON_ADB(LnNP_N)
!CDIRR ON_ADB(LnNZ_N)
!CDIRR ON_ADB(LnNP_Z)
!CDIRR ON_ADB(LnNZ_Z)
!CDIRR ON_ADB(LnNP_T)
!CDIRR ON_ADB(LnNZ_T)
!CDIRR ON_ADB(LnNP_Q)
!CDIRR ON_ADB(LnNZ_Q)
!NEC$ NOVECTOR
sum_lnNP_0 = Sum (LnNP_N(-NM+2:0))
!NEC$ NOVECTOR
sum_lnNZ_0 = Sum (LnNZ_N(-NM+2:0))

#if defined (__NEC__)

!--- 2.5. Calculation of T derivatives
!--- 2.6. Calculation of Q derivatives

!NEC$ shortloop
!NEC$ ivdep
Do i = 1, NLev
   LnNP_T(i) = LnNP_N(i)  *  LnN_T(I1, I2, ID)%P(i) &
             + sum_lnNP_0 * LnNM_T(I1, I2, ID)%P(i)
   LnNZ_T(i) = LnNZ_N(i)  *  LnN_T(I1, I2, ID)%P(i) &
             + sum_lnNZ_0 * LnNM_T(I1, I2, ID)%P(i)
   LnNP_Q(i) = LnNP_N(i)  *  LnN_Q(I1, I2, ID)%P(i) &
             + sum_lnNP_0 * LnNM_Q(I1, I2, ID)%P(i)
   LnNZ_Q(i) = LnNZ_N(i)  *  LnN_Q(I1, I2, ID)%P(i) &
             + sum_lnNZ_0 * LnNM_Q(I1, I2, ID)%P(i)
End Do
!!!CDIR SHORTLOOP SELECT(VECTOR)
!NEC$ shortloop
Do i = 1, NLev
!NEC$ shortloop
   Do k = 1, Nlev
      LnNP_T(i) = LnNP_T(i) + LnNP_Z(k) * Z_T(I1, I2, ID)%M(k,i)
      LnNZ_T(i) = LnNZ_T(i) + LnNZ_Z(k) * Z_T(I1, I2, ID)%M(k,i)
      LnNP_Q(i) = LnNP_Q(i) + LnNP_Z(k) * Z_Q(I1, I2, ID)%M(k,i)
      LnNZ_Q(i) = LnNZ_Q(i) + LnNZ_Z(k) * Z_Q(I1, I2, ID)%M(k,i)
   End Do
End Do

!--- 2.7. Calculation of P derivatives

LnNP_P = sum_lnNP_0 * LnNM_P(I1, I2, ID)
LnNZ_P = sum_lnNZ_0 * LnNM_P(I1, I2, ID)
!NEC$ shortloop
Do k = 1, Nlev
   LnNP_P = LnNP_P &
          + LnNP_Z(k) *   Z_P(I1, I2, ID)%P(k) &
          + LnNP_N(k) * LnN_P(I1, I2, ID)%P(k)
   LnNZ_P = LnNZ_P &
          + LnNZ_Z(k) *   Z_P(I1, I2, ID)%P(k) &
          + LnNZ_N(k) * LnN_P(I1, I2, ID)%P(k)
end Do

#else /* !defined (__NEC__) */

!--- 2.5. Calculation of T derivatives

!$omp parallel private(i)
!$omp do schedule(static)
Do i = 1, NLev
   LnNP_T(i) = &
      Sum(LnNP_Z(1:NLev)*Z_T(I1, I2, ID)%M(1:NLev,i)) + &
      LnNP_N(i)*LnN_T(I1, I2, ID)%P(i)                + &
      sum_lnNP_0*LnNM_T(I1, I2, ID)%P(i)
   LnNZ_T(i) = &
      Sum(LnNZ_Z(1:NLev)*Z_T(I1, I2, ID)%M(1:NLev,i)) + &
      LnNZ_N(i)*LnN_T(I1, I2, ID)%P(i)                + &
      sum_lnNZ_0*LnNM_T(I1, I2, ID)%P(i)
End Do
!$omp end do nowait


!--- 2.6. Calculation of Q derivatives

!$omp do schedule(static)
Do i = 1, NLev
   LnNP_Q(i) = &
      Sum(LnNP_Z(1:NLev)*Z_Q(I1, I2, ID)%M(1:NLev,i)) + &
      LnNP_N(i)*LnN_Q(I1, I2, ID)%P(i)                + &
      sum_lnNP_0 * LnNM_Q(I1, I2, ID)%P(i)
   LnNZ_Q(i) = &
      Sum(LnNZ_Z(1:NLev)*Z_Q(I1, I2, ID)%M(1:NLev,i)) + &
      LnNZ_N(i)*LnN_Q(I1, I2, ID)%P(i)                + &
      sum_lnNZ_0 * LnNM_Q(I1, I2, ID)%P(i)
End Do
!$omp end do
!$omp end parallel


!--- 2.7. Calculation of P derivatives

LnNP_P = &
   Sum(LnNP_Z(1:NLev)*Z_P(I1, I2, ID)%P(1:NLev))    + &
   Sum(LnNP_N(1:NLev)*LnN_P(I1, I2, ID)%P(1:NLev))  + &
   sum_lnNP_0 * LnNM_P(I1, I2, ID)
LnNZ_P = &
   Sum(LnNZ_Z(1:NLev)*Z_P(I1, I2, ID)%P(1:NLev))    + &
   Sum(LnNZ_N(1:NLev)*LnN_P(I1, I2, ID)%P(1:NLev))  + &
   sum_lnNZ_0 * LnNM_P(I1, I2, ID)

#endif
FTRACE_END  ("Vertical_Interpolation_LnN_adj:2.5_7")


!--- 2.8. Calculation of limits of Z

Zmin = gf% s(ix)% Zsur
Zmax = Max(ZLL(1), ZLL(NLev))


!------------------------------
! 3. Process optional arguments
!------------------------------

if (present (T)) then

   Call  Init_Spline             &
     (Z  (I1, I2, ID)%P(1:),     & ! <-- Argument grid
      real(gf% T(:,ix),double),  & ! <-- Gridded function
      D2                         ) ! --> 2nd derivative of spline

   Call Spline &
    (Z  (I1, I2, ID)%P(1:),    & ! <-- Argument grid
     real(gf% T(:,ix),double), & ! <-- Gridded function
     d2,                       & ! <-- 2nd derivative of spline
     Real(ZP, Double),         & ! <-- Interpolation point
     T                         ) ! --> Interpolated function value

endif

if (present (Q)) then

   Call  Init_Spline             &
     (Z  (I1, I2, ID)%P(1:),     & ! <-- Argument grid
      real(gf% Q(:,ix),double),  & ! <-- Gridded function
      D2                         ) ! --> 2nd derivative of spline

   Call Spline &
    (Z  (I1, I2, ID)%P(1:),    & ! <-- Argument grid
     real(gf% Q(:,ix),double), & ! <-- Gridded function
     d2,                       & ! <-- 2nd derivative of spline
     Real(ZP, Double),         & ! <-- Interpolation point
     Q                         ) ! --> Interpolated function value

endif

if (present (Gp)) then

   Call  Init_Spline             &
     (Z  (I1, I2, ID)%P(1:),     & ! <-- Argument grid
      Hf,                        & ! <-- Gridded function
      D2                         ) ! --> 2nd derivative of spline

   Call Spline &
    (Z  (I1, I2, ID)%P(1:),    & ! <-- Argument grid
     Hf,                       & ! <-- Gridded function
     d2,                       & ! <-- 2nd derivative of spline
     Real(ZP, Double),         & ! <-- Interpolation point
     Gp                        ) ! --> Interpolated function value

endif

if (present (P)) then

   Call  Init_Spline             &
     (Z  (I1, I2, ID)%P(1:),     & ! <-- Argument grid
      Pf,                        & ! <-- Gridded function
      D2                         ) ! --> 2nd derivative of spline

   Call Spline &
    (Z  (I1, I2, ID)%P(1:),    & ! <-- Argument grid
     Pf,                       & ! <-- Gridded function
     d2,                       & ! <-- 2nd derivative of spline
     Real(ZP, Double),         & ! <-- Interpolation point
     P                         ) ! --> Interpolated function value

endif

if (present (Ps)) Ps = gf% s(ix)% Psur
if (present (Hs)) Hs = gf% s(ix)% Hsur / g_ave
if (present (Gu)) Gu = gf% s(ix)% Gundu

!----------------------------------------------------------
! 3. MEMORY DEALLOCATION
!----------------------------------------------------------

!Deallocate &
!  (LnNP_Z,    & ! d(LnNP)/d(Z(i))
!   LnNP_N,    & ! d(LnNP)/d(LnN(i))
!   LnNZ_Z,    & ! d(LnNZ)/d(Z(i))
!   LnNZ_N)      ! d(LnNZ)/d(LnN(i))


End Subroutine Vertical_Interpolation_LnN_adj



!==========================================================
Subroutine Interpolate_Refractivity_adj  &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NG,       & ! --> Interpolated dN/d(alt,lat,lon)
   NH,       & ! --> Interpolated hessian matrix of N
   IDX,      & ! --> Subgrid indices [point, index]
   NP_T,     & ! --> d(NP)/d(T(IGP,k))
   NP_Q,     & ! --> d(NP)/d(Q(IGP,k))
   NP_P,     & ! --> d(NP)/d(Psur(IGP))
   NG_T,     & ! --> d(NG(i))/d(T(IGP,k))
   NG_Q,     & ! --> d(NG(i))/d(Q(IGP,k))
   NG_P,     & ! --> d(NG(i))/d(Psur(IGP))
   T,        & ! ~~> Interpolated temperature
   Q,        & ! ~~> Interpolated specific humidity
   Gp,       & ! ~~> Interpolated geopotential height
   Pf,       & ! ~~> Interpolated pressure
   Ps,       & ! ~~> Interpolated surface pressure
   Hs,       & ! ~~> Interpolated surface geopotential height
   Gu        ) ! ~~> Interpolated geoid undulation
!
! Interpolation of global refractivity field from ECHAM:
! Adjoint version.
!----------------------------------------------------------
! Method:
!   1. Vertical spline interpolation to given altitude.
!   2. Summation of 2D Lon/Lat interpolation series:
!   f(z,Lat,Lon) = sum WG f(z,Lat ,Lon )
!                   i    i       i    i
!   Index i enumerates points of the interpolation
!   subgrid. For rectangular grid, weighting function WG
!   is calculated for polynomial interpolation with the
!   polynomial power defined by the interpolation subgrid
!   dimension NI. For icosahedral grid WG is defined as
!   symplectic weights of the vertices of surrounding
!   triangle.
!----------------------------------------------------------
! (C) Copyright 2000-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   5.1   | 13 Apr 2000 | Basic non-adjoint version.
!   1.0   | 19 Apr 2000 | Adjoint version.
!   2.0   | 08 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!   2.1   | 15 Apr 2002 | Multiple icosahedral
!         |             | indices taken into account.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_grid, only: &
! Imported Routines:
    Longitude_Interpolation,  &
    Latitude_Interpolation
!
Use ICO_grid, only: &
! Imported Array Variables:
    nspoke
! Imported Routines:
!   Multiple_Points
!
Use ICO_interpolation, only: &
! Imported Routines:
    Setup_Interpolation
!
Use mo_icon_grid,      only: &
    search_icon_global
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   PLon          ! Latitude [deg]
!
Real(Double), Intent(In) :: &
   PLat          ! Longitude [deg]
!
Real(Double), Intent(In) :: &
   ZP            ! Altitude [km]
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Zmin          ! Minimum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   Zmax          ! Maximum model Z for this lon/lat
!
Real(Double), Intent(Out) :: &
   NP            ! Interpolated N
!
Real(Double), Intent(Out) :: &
   NG(3)         ! Interpolated gradient dN/d(alt,lat,lon)
!
Real(Double), Intent(Out) :: &
   NH(3,3)       ! Interpolated hessian
                 ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))N
!
Integer, Intent(Out) :: &
   IDX(:,:)      ! Subgrid indices [point, index]
!
Real(Double), Intent(Out) :: &
   NP_T(:,:)     ! d(NP)/d(T(IGP,k))
!
Real(Double), Intent(Out) :: &
   NP_Q(:,:)     ! d(NP)/d(Q(IGP,k))
!
Real(Double), Intent(Out) :: &
   NP_P(:)       ! d(NP)/d(Psur(IGP))
!
Real(Double), Intent(Out) :: &
   NG_T(:,:,:)   ! d(NG(i))/d(T(IGP,k))
!
Real(Double), Intent(Out) :: &
   NG_Q(:,:,:)   ! d(NG(i))/d(Q(IGP,k))
!
Real(Double), Intent(Out) :: &
   NG_P(:,:)     ! d(NG(i))/d(Psur(KLon,KLat))
!
Real(Double), Intent(Out), Optional :: &
   T, Q, Gp, Pf, Ps, Hs, Gu
!----------------------------------------------------------
#if   defined (_CRAYFTN) \
  || (defined (__GFORTRAN__) && ((__GNUC__ * 100 + __GNUC_MINOR__) >= 406)  \
                             && ((__GNUC__ * 100 + __GNUC_MINOR__) <  800)) \
  || (defined (__INTEL_COMPILER) && (__INTEL_COMPILER >= 1500)) \
  || (defined (__NEC__) && (__NEC_VERSION__ >= 30004)) \
  || (defined (__PGI) && (__PGIC__ >= 15) && (__PGIC__ != 17))
contiguous :: NG_T, NG_Q, NG_P
contiguous :: NP_T, NP_Q, NP_P
#endif
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev    ! Number of levels
Integer        :: NGP     ! Subgrid dimension
Integer        :: IGP     ! Subgrid index
Integer        :: KLon    ! Longitude subgrid index
Integer        :: KLat    ! Latitude subgrid index
Integer        :: m1, m2  ! Interpolation triangle indices
!Integer       :: NMP     ! Number of multiple points
Real(Double)   :: FP      ! Interpolated ln(N)
Integer        :: i       ! Dimension index
#if defined (__NEC__)
Integer        :: k       ! Auxiliary level index
#endif
!
! Local Arrays:
!
! --- Interpolation on rectangular grid
!
Integer :: ILon(NG1)    ! Longitude interpolation subgrid
Integer :: ILat(NG1)    ! Latitude interpolation subgrid
Real(Double) :: &
   WLon(NG1),         & ! Weight for longitudinal interpolation
   WLon1(NG1),        & ! 1st derivative of WLon
   WLon2(NG1),        & ! 1st derivative of WLon
   WLat(NG1),         & ! Weight for latitudinal interpolation
   WLat1(NG1),        & ! 1st derivative of WLat
   WLat2(NG1)           ! 1st derivative of WLat
!
! --- Interpolation on icosahedral grid
!
Real(Double), save :: & ! Saved values of:
   PRLon(1)  = -999,  & ! Longitude of interpolation point
   PRLat(1)  = -999,  & ! Latitude  of interpolation point
   SWS(3,1)  =    0     ! Symplectic weights for triangle vertices
Integer, save :: &
   KNI       =   -1,  & ! GME resolution parameter
   KIDX(4,1) =    0     ! Interpolation triangle index
!Integer :: &
!  MIDX(5,3)            ! Multiple point grid indices
!
! --- Interpolation on ICON grid
!
Integer, save :: &
   jl(3,1)   =   -1,  & ! Line  indices
   jb(3,1)   =   -1     ! Block indices
Integer :: &
   npr (1)              ! number of neighbour points returned
!
! --- Interpolation weights
!
Real(Double), Dimension(gf% ngp,0:2,0:2) :: &
   W_                   ! Weighting function
Real(Double), Dimension(gf% ncol,0:2,0:2) :: &
   WG                   ! Weighting function
                        ! [Point, Lon.derivative, Lat.derivative]
!
! --- Interpolated functions
!
Real(Double), Dimension(gf% ncol) :: &
   ZNmin,             & ! Minimum model Z on subgrid
   ZNmax,             & ! Maximum model Z on subgrid
   FV,                & ! Vertically interpolated F on subgrid
   FZ,                & ! Vertical gradient on subgrid
   F2Z                  ! Second vertical derivative on subgrid
Real(Double) :: &
   FG(3),             & ! Interpolated dln(N)/d(alt,lat,lon)
   FH(3,3)              ! Interpolated hessian
                        ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))ln(N)
!
! --- T, Q, and P derivatives
!
Real(Double), Dimension(gf% ncol, gf% nz) :: &
   FV_T,              & ! d(LnNP)/d(T(i))
   FV_Q,              & ! d(LnNP)/d(Q(i))
   FZ_T,              & ! d(LnNZ)/d(T(i))
   FZ_Q                 ! d(LnNZ)/d(Q(i))
Real(Double), Dimension(gf% ncol) :: &
   FV_P,              & ! d(LnNP)/d(Psur)
   FZ_P                 ! d(LnNZ)/d(Psur)
!
! --- interpolated T, Q, Gp (diagnostics)
!
Real(Double), Allocatable :: TI(:),QI(:),GpI(:),PfI(:),PsI(:),HsI(:),GuI(:)


integer :: pe_i1_i2_id(4)
!----------------------------------------------------------
! Global variables used:
!
!   Hsur(:,:)    ! Surface geopotential [gpm]
!   Psur(:,:)    ! Surface pressure [Pa]
!   T(:,:,:)     ! Temperature [K]
!   Q(:,:,:)     ! Specific humidity [kg/kg]
!   XLon(:)      ! Longitude grid [deg]
!   XLat(:)      ! Latitude grid [deg]
!   GCLat(:)     ! Geocentric latitudes [rad]
!   Zsur(:,:)    ! Surface altitude [km]
!   Z(:,:,:)     ! Altitudes [km]
!   LnN(:,:,:)   ! Ln of refractive index
!   D2N(:,:,:)   ! Interpolation coefficients of Ln(N)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------

!--- 1.1. Subgrid size

!Select Case (gf% Grid_Type)
!   Case (WMO6_GAUSSIAN, WMO6_LATLON)
!      NGP = NG1**2
!   Case (DWD6_ICOSAHEDRON, DWD6_ICON)
!      NGP = 3
!   case default
!      call finish('Interpolate_Refractivity_adj','invalid grid type')
!End Select

!if (NGP /= gf% ngp) then
!   write(0,*) "gf% ngp /= NGP:", gf% ngp, NGP
!   call finish('Interpolate_Refractivity_adj','incorrect initialization of gf')
!end if

NGP = gf% ncol  ! Number of actually used column(s)


!--- 1.2. Number of vertical levels

NLev = gf% nz


!--- 1.3. Memory allocation

!FTRACE_BEGIN("Interpolate_Refractivity_adj:1.3")
!Allocate  &
!  (WG(NGP,0:2,0:2), & ! Weighting function
!   ZNmin(NGP),      & ! Minimum model Z on subgrid
!   ZNmax(NGP),      & ! Maximum model Z on subgrid
!   FV(NGP),         & ! Vertically interpolated F on subgrid
!   FZ(NGP),         & ! Vertical gradient on subgrid
!   F2Z(NGP),        & ! Second vertial derivative on subgrid
!   FV_T(NGP,NLev),  & ! d(LnNP)/d(T(i))
!   FV_Q(NGP,NLev),  & ! d(LnNP)/d(Q(i))
!   FV_P(NGP),       & ! d(LnNP)/d(Psur)
!   FZ_T(NGP,NLev),  & ! d(LnNZ)/d(T(i))
!   FZ_Q(NGP,NLev),  & ! d(LnNZ)/d(Q(i))
!   FZ_P(NGP))         ! d(LnNZ)/d(Psur)

if (present(T)) &
  allocate (TI(NGP),QI(NGP),GpI(NGP),PfI(NGP),PsI(NGP),HsI(NGP),GuI(NGP))
!FTRACE_END  ("Interpolate_Refractivity_adj:1.3")

!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION WEIGHTS
!----------------------------------------------------------


Select Case (gf% Grid_Type)

   !--- 2.1. Interpolation on rectangular grid

   Case (WMO6_GAUSSIAN, WMO6_LATLON)

      !--- 2.1.1. Longitudinal interpolation

      Call Longitude_Interpolation &
        (gf% glon(:),   & ! <-- Longitude grid [deg]
         PLon,          & ! <-- Longitude of point [deg]
         ILon,          & ! --> Interpolation subgrid
         WLon,          & ! --> Weights of subgrid points
         WLon1,         & ! ~~> Weight derivatives
         WLon2)           ! ~~> Weight 2nd derivatives

      !--- 2.1.2. Latitudinal interpolation

      Call Latitude_Interpolation &
        (gf% glat(:),   & ! <-- Latitude grid [deg]
         PLat,          & ! <-- Latitude of point [deg]
         ILat,          & ! --> Interpolation subgrid
         WLat,          & ! --> Weights of subgrid points
         WLat1,         & ! ~~> Weight derivatives
         WLat2)           ! ~~> Weight 2nd derivatives

      !--- 2.1.3. Interpolation indices and weights

      IGP = 0

      Do KLat = 1,NG1
         Do KLon = 1,NG1
            IGP         = IGP + 1
            IDX(IGP,1)  = ILon(KLon)
            IDX(IGP,2)  = ILat(KLat)
            IDX(IGP,3)  = 1
            W_(IGP,0,0) = WLon(KLon)*WLat(KLat)
            W_(IGP,1,0) = WLon1(KLon)*WLat(KLat)
            W_(IGP,0,1) = WLon(KLon)*WLat1(KLat)
            W_(IGP,1,1) = WLon1(KLon)*WLat1(KLat)
            W_(IGP,2,0) = WLon2(KLon)*WLat(KLat)
            W_(IGP,0,2) = WLon(KLon)*WLat2(KLat)
         End Do
      End Do

   !--- 2.2. Interpolation on icosahedral grid

   Case (DWD6_ICOSAHEDRON)

      !--- 2.2.1. Interpolation indices and weights

FTRACE_BEGIN("Interpolate_Refractivity_adj:2.2.1")
      if (PLon /= PRLon(1) .or. PLat /= PRLat(1) .or. gf% ni /= KNI) then
         PRLon(1) = PLon
         PRLat(1) = PLat
         KNI      = gf% ni

         Call Setup_Interpolation &
           (1,          & ! <-- Number of interpolation points
            PRLon,      & ! <-- Point longitudes [deg]
            PRLat,      & ! <-- Point latitudes [deg]
            KIDX,       & ! --> Index array for interpolation
            SWS)          ! --> Symplectic weights for triangle vertices
      end if

      W_(:,:,:) = 0.0_Double
      W_(:,0,0) = SWS(:,1)
      m1        = KIDX(4,1)
      m2        = Mod(m1,6) + 1
      IDX(1,:)  = KIDX(1:3,1)
      IDX(2,1)  = KIDX(1,1) + nspoke(1,m1)
      IDX(2,2)  = KIDX(2,1) + nspoke(2,m1)
      IDX(2,3)  = KIDX(3,1)
      IDX(3,1)  = KIDX(1,1) + nspoke(1,m2)
      IDX(3,2)  = KIDX(2,1) + nspoke(2,m2)
      IDX(3,3)  = KIDX(3,1)
FTRACE_END  ("Interpolate_Refractivity_adj:2.2.1")

      !--- 2.2.2. Computation of equivalent indices

!      Do IGP = 1, NGP
!         Call Multiple_Points &
!           (IDX(IGP,:),  & ! <-- Grid point indices
!            NMP,         & ! --> Number of multiple points
!            MIDX)          ! --> Multiple point indices
!         IDX(IGP,:) = MIDX(1,:)
!      End Do

      do igp = 1, 3
        pe_i1_i2_id  = gf% marr (:,IDX(igp,1),IDX(igp,2),IDX(igp,3))
        IDX(igp,1:3) = pe_i1_i2_id(2:4)
      end do

   !--- 2.?. Interpolation on ICON triangular grid

   Case (DWD6_ICON)
      if (PLon /= PRLon(1) .or. PLat /= PRLat(1) .or. gf% ni /= KNI) then
         PRLon(1) = PLon
         PRLat(1) = PLat
         KNI      = gf% ni
         call search_icon_global   &
              (gf% grid% icongrid, &! <-- Interpolation metadata
               PRLon,              &! <-- Point longitude [deg]
               PRLat,              &! <-- Point latitude  [deg]
               jl, jb,             &! --> Line and block indices
               SWS,                &! --> Interpolation weights
               NPR                 )! --> Number of neighbours
      end if
      W_(:,:,:)  = 0.0_Double
      W_(:,0,0)  = SWS(1:3,1)
      idx(1:3,1) = jl (1:3,1)
      idx(1:3,2) = jb (1:3,1)
      idx(1:3,3) = 1
      do igp = 1, 3
         pe_i1_i2_id  = gf% marr (:,IDX(igp,1),IDX(igp,2),IDX(igp,3))
         IDX(igp,1:3) = pe_i1_i2_id(2:4)
      end do

!   !--- 2.?. Single column: no horizontal interpolation.
!
!   Case (DWD6_GRID_NONE)
!      W_(:,:,:)  = 0.0_Double   ! Trivial weights
!      W_(1,0,0)  = 1.0_Double
!      IDX(1,1:3) = 1            ! Dummy indices

   case default
      call finish('Interpolate_Refractivity_adj','invalid grid type')
End Select

select case (gf% hint_mode)
case (0)
   !--------------------------
   ! Nearest model column only
   !--------------------------
   i = maxloc (W_(:,0,0), dim=1)
   IDX(1,1:3) = IDX(i,1:3)
   W_(:,:,:)  = 0.0_Double      ! Trivial weights
   W_(1,0,0)  = 1.0_Double
end select

WG(1:NGP,:,:) = W_(1:NGP,:,:)   ! Weights of actually used column(s)


!----------------------------------------------------------
! 3. VERTICAL INTERPOLATION
!----------------------------------------------------------


if(present(T).or.present(Q).or.present(Gp)) then
FTRACE_BEGIN("Interpolate_Refractivity_adj:3a")
 Do IGP = 1,NGP
   Call Vertical_Interpolation_LnN_adj  &
     (IDX(IGP,1),          & ! <-- 1st grid index
      IDX(IGP,2),          & ! <-- 2nd grid index
      IDX(IGP,3),          & ! <-- Diamond index
      ZP,                  & ! <-- Altitude [km]
      ZNmin(IGP),          & ! --> Minimum model Z for this lon/lat
      ZNmax(IGP),          & ! --> Maximum model Z for this lon/lat
      FV(IGP),             & ! --> Interpolated Ln(N)
      FV_T(IGP,:),         & ! --> d(LnNP)/d(T(i))
      FV_Q(IGP,:),         & ! --> d(LnNP)/d(Q(i))
      FV_P(IGP),           & ! --> d(LnNP)/d(Psur)
      FZ(IGP),             & ! --> Interpolated dLn(N)/dZ
      FZ_T(IGP,:),         & ! --> d(LnNZ)/d(T(i))
      FZ_Q(IGP,:),         & ! --> d(LnNZ)/d(Q(i))
      FZ_P(IGP),           & ! --> d(LnNZ)/d(Psur)
      F2Z(IGP),            & ! --> Interpolated d2Ln(N)/dZ2
      TI (IGP),            & ! ~~> Interpolated temperature
      QI (IGP),            & ! ~~> Interpolated humidity
      GpI(IGP),            & ! ~~> Interpolated geopotential height
      PfI(IGP),            & ! ~~> Interpolated pressure
      PsI(IGP),            & ! ~~> Surface pressure at this location
      HsI(IGP),            & ! ~~> Surface geopotential height
      GuI(IGP)             ) ! ~~> Geoid undulation
 End Do
FTRACE_END  ("Interpolate_Refractivity_adj:3a")
else
FTRACE_BEGIN("Interpolate_Refractivity_adj:3b")
 Do IGP = 1,NGP
   Call Vertical_Interpolation_LnN_adj  &
     (IDX(IGP,1),          & ! <-- 1st grid index
      IDX(IGP,2),          & ! <-- 2nd grid index
      IDX(IGP,3),          & ! <-- Diamond index
      ZP,                  & ! <-- Altitude [km]
      ZNmin(IGP),          & ! --> Minimum model Z for this lon/lat
      ZNmax(IGP),          & ! --> Maximum model Z for this lon/lat
      FV(IGP),             & ! --> Interpolated Ln(N)
      FV_T(IGP,:),         & ! --> d(LnNP)/d(T(i))
      FV_Q(IGP,:),         & ! --> d(LnNP)/d(Q(i))
      FV_P(IGP),           & ! --> d(LnNP)/d(Psur)
      FZ(IGP),             & ! --> Interpolated dLn(N)/dZ
      FZ_T(IGP,:),         & ! --> d(LnNZ)/d(T(i))
      FZ_Q(IGP,:),         & ! --> d(LnNZ)/d(Q(i))
      FZ_P(IGP),           & ! --> d(LnNZ)/d(Psur)
      F2Z(IGP))              ! --> Interpolated d2Ln(N)/dZ2
 End Do
FTRACE_END  ("Interpolate_Refractivity_adj:3b")
endif

!----------------------------------------------------------
! 4. HORIONTAL INTERPOLATION
!----------------------------------------------------------


!--- 4.1. Function and height limits
FTRACE_BEGIN("Interpolate_Refractivity_adj:4.1_2")

if(present(T))  T  = Sum(WG(:,0,0)*TI (:))
if(present(Q))  Q  = Sum(WG(:,0,0)*QI (:))
if(present(Gp)) Gp = Sum(WG(:,0,0)*GpI(:))
if(present(Pf)) Pf = Sum(WG(:,0,0)*PfI(:))
if(present(Ps)) Ps = Sum(WG(:,0,0)*PsI(:))
if(present(Hs)) Hs = Sum(WG(:,0,0)*HsI(:))
if(present(Gu)) Gu = Sum(WG(:,0,0)*GuI(:))

!!!CDIR BEGIN SHORTLOOP
!NEC$ NOVECTOR
FP    = Sum(WG(1:NGP,0,0)*FV(:))
!NEC$ NOVECTOR
Zmin  = Sum(WG(1:NGP,0,0)*ZNmin(:))
!NEC$ NOVECTOR
Zmax  = Sum(WG(1:NGP,0,0)*ZNmax(:))
!!!CDIR END

!--- 4.2. Gradient

!!!CDIR BEGIN SHORTLOOP
!NEC$ NOVECTOR
FG(1) = Sum(WG(1:NGP,0,0)*FZ(:))
!NEC$ NOVECTOR
FG(2) = Sum(WG(1:NGP,0,1)*FV(:))
!NEC$ NOVECTOR
FG(3) = Sum(WG(1:NGP,1,0)*FV(:))
!!!CDIR END

FTRACE_END  ("Interpolate_Refractivity_adj:4.1_2")

!--- 4.3. T, Q, and P derivatives
FTRACE_BEGIN("Interpolate_Refractivity_adj:4.3")

!NEC$ loop_count(4)
Do IGP = 1, NGP

   NP_P(IGP)     = WG(IGP,0,0)*FV_P(IGP)

!!!CDIR ARRAYCOMB
!NEC$ SHORTLOOP
   NP_T(IGP,:)   = WG(IGP,0,0)*FV_T(IGP,:)
!NEC$ SHORTLOOP
   NP_Q(IGP,:)   = WG(IGP,0,0)*FV_Q(IGP,:)

!NEC$ SHORTLOOP
   NG_T(1,IGP,:) = WG(IGP,0,0)*FZ_T(IGP,:)
!NEC$ SHORTLOOP
   NG_T(2,IGP,:) = WG(IGP,0,1)*FV_T(IGP,:)
!NEC$ SHORTLOOP
   NG_T(3,IGP,:) = WG(IGP,1,0)*FV_T(IGP,:)

!NEC$ SHORTLOOP
   NG_Q(1,IGP,:) = WG(IGP,0,0)*FZ_Q(IGP,:)
!NEC$ SHORTLOOP
   NG_Q(2,IGP,:) = WG(IGP,0,1)*FV_Q(IGP,:)
!NEC$ SHORTLOOP
   NG_Q(3,IGP,:) = WG(IGP,1,0)*FV_Q(IGP,:)
!!!CDIR END ARRAYCOMB

   NG_P(1,IGP)   = WG(IGP,0,0)*FZ_P(IGP)
   NG_P(2,IGP)   = WG(IGP,0,1)*FV_P(IGP)
   NG_P(3,IGP)   = WG(IGP,1,0)*FV_P(IGP)

End Do

FTRACE_END  ("Interpolate_Refractivity_adj:4.3")

!--- 4.4. Hessian

FTRACE_BEGIN("Interpolate_Refractivity_adj:4.4")
!!!CDIR BEGIN SHORTLOOP
!NEC$ NOVECTOR
FH(1,1) = Sum(WG(1:NGP,0,0)*F2Z(:))
!NEC$ NOVECTOR
FH(2,2) = Sum(WG(1:NGP,0,2)*FV(:))
!NEC$ NOVECTOR
FH(3,3) = Sum(WG(1:NGP,2,0)*FV(:))
!NEC$ NOVECTOR
FH(1,2) = Sum(WG(1:NGP,0,1)*FZ(:))
!NEC$ NOVECTOR
FH(1,3) = Sum(WG(1:NGP,1,0)*FZ(:))
!NEC$ NOVECTOR
FH(3,2) = Sum(WG(1:NGP,1,1)*FV(:))
FH(2,1) = FH(1,2)
FH(3,1) = FH(1,3)
FH(2,3) = FH(3,2)
!!!CDIR END
FTRACE_END  ("Interpolate_Refractivity_adj:4.4")


!----------------------------------------------------------
! 5. TRANSFORM FROM LN(N) TO N
!----------------------------------------------------------

!--- 5.1. Refractive index
FTRACE_BEGIN("Interpolate_Refractivity_adj:5.1")

NP = Exp(FP)

#if defined (__NEC__)
! Manual loop exchange to work around issues in nfort <= 3.0.4
!NEC$ loop_count(4)
do IGP = 1, NGP
!NEC$ shortloop
!NEC$ ivdep
  do k = 1, Nlev
     NP_T(IGP,k) = NP*NP_T(IGP,k)
     NP_Q(IGP,k) = NP*NP_Q(IGP,k)
  end do
end do
#else
NP_T(:,:)   = NP*NP_T(:,:)
NP_Q(:,:)   = NP*NP_Q(:,:)
#endif

NP_P(:)     = NP*NP_P(:)

FTRACE_END  ("Interpolate_Refractivity_adj:5.1")

!--- 5.2. Refractive gradient
FTRACE_BEGIN("Interpolate_Refractivity_adj:5.2")

NG(:) = NP*FG(:)

#if defined (__NEC__)
Do i=1,3
! Manual loop exchange to work around issues in nfort <= 3.0.4
!NEC$ loop_count(4)
   Do IGP = 1, NGP
!NEC$ shortloop
!NEC$ ivdep
     Do k = 1, Nlev
        NG_T(i,IGP,k) = NP_T(IGP,k)*FG(i) + NP*NG_T(i,IGP,k)
        NG_Q(i,IGP,k) = NP_Q(IGP,k)*FG(i) + NP*NG_Q(i,IGP,k)
     End Do
   End Do
   NG_P(i,:)   = NP_P(:)*FG(i)   + NP*NG_P(i,:)
End Do
#else
Do i=1,3
   NG_T(i,:,:) = NP_T(:,:)*FG(i) + NP*NG_T(i,:,:)
   NG_Q(i,:,:) = NP_Q(:,:)*FG(i) + NP*NG_Q(i,:,:)
   NG_P(i,:)   = NP_P(:)*FG(i)   + NP*NG_P(i,:)
End Do
#endif

FTRACE_END  ("Interpolate_Refractivity_adj:5.2")

!--- 5.3. Refractive hessian

Do i=1,3
   NH(i,:) = NP*(FH(i,:) + FG(i)*FG(:))
End Do


!----------------------------------------------------------
! 6. MEMORY DEALLOCATION
!----------------------------------------------------------

if (present(T)) deallocate (TI,QI,GpI,PfI,PsI,HsI,GuI)


End Subroutine Interpolate_Refractivity_adj



!==========================================================
Subroutine ECHAM_NGradN_adj &
  (X,        & ! <-- Cartesian coordinates of point
   NGN,      & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
   Stat,     & ! --> Error status
   IDX,      & ! --> Subgrid indices [point, index]
   NP_T,     & ! --> d(NP)/d(T(IGP,k))
   NP_Q,     & ! --> d(NP)/d(Q(IGP,k))
   NP_P,     & ! --> d(NP)/d(Psur(IGP))
   NGN_T,    & ! --> d(NGN(i))/d(T(IGP,k))
   NGN_Q,    & ! --> d(NGN(i))/d(Q(IGP,k))
   NGN_P)      ! --> d(NGN(i))/d(Psur(IGP))
!
! Calculation of interpolated (1 + N)*Grad(N) for
! ECHAM global fields: Adjoint version.
!----------------------------------------------------------
! Method:
!   Calculation of gradient in geodetic coordinates and
!   transform to Cartesian coordinates.
!----------------------------------------------------------
! (C) Copyright 1999-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   3.0   | 22 Jun 1999 | Basic non-adjoint version.
!   1.0   | 14 Apr 2000 | Adjoint version.
!   2.0   | 08 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!----------------------------------------------------------
! Modules used:
!
Use Coordinates, only: &
! Imported Type Definitions:
    Cartesian
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Routines:
    Geod_from_Cart,  &
    Jacobian_GC,     &
    Hessian_GC
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In)  :: &
   X               ! Cartesian coordinates of point
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   NGN             ! Interpolated (1 + N)*Grad(N)
!
Real(Double), Intent(Out)    :: &
   NP              ! Interpolated N
!
Real(Double), Intent(Out)    :: &
   NHN(3,3)        ! Interpolated Grad x (1+N)Grad(N)
!
Integer, Intent(Out)         :: &
   Stat            ! Error status:
                   !   0 - point above surface
                   !   1 - point under surface
!
Integer, Intent(Out) :: &
   IDX(:,:)        ! Subgrid indices [point, index]
!
Real(Double), Intent(Out) :: &
   NP_T(:,:)       ! d(NP)/d(T(IGP,i))
!
Real(Double), Intent(Out) :: &
   NP_Q(:,:)       ! d(NP)/d(Q(IGP,i))
!
Real(Double), Intent(Out) :: &
   NP_P(:)         ! d(NP)/d(Psur(IGP))
!
Real(Double), Intent(Out) :: &
   NGN_T(:,:,:)    ! d(NG(i))/d(T(IGP,i))
!
Real(Double), Intent(Out) :: &
   NGN_Q(:,:,:)    ! d(NG(i))/d(Q(IGP,i))
!
Real(Double), Intent(Out) :: &
   NGN_P(:,:)      ! d(NG(i))/d(Psur(IGP))
!----------------------------------------------------------
! Local Scalars:
!
Type(Geodetic) :: G     ! Geodetic coordinates of point
Real(Double)   :: PLon  ! Longitude of point
Real(Double)   :: PLat  ! Latitude of point
Real(Double)   :: ZP    ! Altitude of point
Real(Double)   :: Zmin  ! Minimum model Z
Real(Double)   :: Zmax  ! Maximum model Z
Integer        :: NGP   ! Number of subgrid points for interpolation
Integer        :: NLev  ! Number of model levels
Integer        :: i, j  ! Dimension indexes
!
! Local Arrays:
!
Real(Double)   :: &
   NG(3)        ! Interpolated Grad(N) in geodetic coordinates
Real(Double)   :: &
   NH(3,3)      ! Interpolated Hess(N) in geodetic coordinates
Real(Double)   :: &
   JGC(3,3),  & ! Jacobian d(Geodetic)/d(Cartesian)
   HGC(3,3,3)   ! Hessian d2(Geodetic)/d(Cartesian)2
Real(Double), Allocatable :: &
   NG_T(:,:,:),  & ! d(NG(i))/d(T(IGP,k))
   NG_Q(:,:,:),  & ! d(NG(i))/d(Q(IGP,k))
   NG_P(:,:)       ! d(NG(i))/d(Psur(IGP))
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INTERPOLATION OF ECHAM-MSIS FIELD
!----------------------------------------------------------


!--- 1.1. Coordinate calculation

G    = Geod_from_Cart(X)
PLon = G%Lambda
PLat = G%Phi
ZP   = G%H


!--- 1.2. Memory allocation

Select Case (gf% Grid_Type)
   Case (WMO6_GAUSSIAN, WMO6_LATLON)
      NGP = NG1**2
   Case (DWD6_ICOSAHEDRON, DWD6_ICON)
      NGP = 3
   case default
      call finish('ECHAM_NGradN_adj','invalid grid type')
End Select

NLev = gf% nz

Allocate &
  (NG_T(3,NGP,NLev),  & ! d(NG(i))/d(T(IGP,k))
   NG_Q(3,NGP,NLev),  & ! d(NG(i))/d(Q(IGP,k))
   NG_P(3,NGP))         ! d(NG(i))/d(Psur(IGP))


!--- 1.3. Interpolation and calculation
!---      of T, P, and Q derivatives

Call Interpolate_Refractivity_adj &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NG,       & ! --> Interpolated dN/d(alt,lat,lon)
   NH,       & ! --> Interpolated hessian matrix of N
   IDX,      & ! --> Subgrid indices [point, index]
   NP_T,     & ! --> d(NP)/d(T(IGP,k))
   NP_Q,     & ! --> d(NP)/d(Q(IGP,k))
   NP_P,     & ! --> d(NP)/d(Psur(IGP))
   NG_T,     & ! --> d(NG(i))/d(T(IGP,k))
   NG_Q,     & ! --> d(NG(i))/d(Q(IGP,k))
   NG_P)       ! --> d(NG(i))/d(Psur(IGP))


!----------------------------------------------------------
! 2. TRANSFORM TO CARTESIAN COORDINATES
!----------------------------------------------------------


!--- 2.1. Transform of gradient

JGC    = Jacobian_GC(G)

NGN%X(:) = (1+NP)*MatMul(NG(:), JGC(:,:))


!--- 2.2. Transform of gradient derivatives

Do i=1,3
   NGN_T(i,:,:) = NG_T(1,:,:)*JGC(1,i) + &
                  NG_T(2,:,:)*JGC(2,i) + &
                  NG_T(3,:,:)*JGC(3,i)
   NGN_Q(i,:,:) = NG_Q(1,:,:)*JGC(1,i) + &
                  NG_Q(2,:,:)*JGC(2,i) + &
                  NG_Q(3,:,:)*JGC(3,i)
   NGN_P(i,:)   = NG_P(1,:)*JGC(1,i) + &
                  NG_P(2,:)*JGC(2,i) + &
                  NG_P(3,:)*JGC(3,i)
End Do


!--- 2.3. Transform of hessian

HGC = Hessian_GC(G)

Do i=1,3
   Do j=1,3
      NHN(i,j) = &
         Sum(JGC(:,i)*NG(:))*Sum(JGC(:,j)*NG(:)) + &
         (1+NP)*Sum(HGC(:,i,j)*NG(:))            + &
         (1+NP)*Sum(JGC(:,i)*MatMul(NH(:,:),JGC(:,j)))
   End Do
End Do


!--- 2.3. Memory deallocation

Deallocate &
  (NG_T,    & ! d(NG(i))/d(T(IGP,k))
   NG_Q,    & ! d(NG(i))/d(Q(IGP,k))
   NG_P)      ! d(NG(i))/d(Psur(IGP))


!----------------------------------------------------------
! 3. STATUS DEFINITION
!----------------------------------------------------------

If (G%H >= Zmin) then
   Stat = 0
Else
   Stat = 1
End If


End Subroutine ECHAM_NGradN_adj



End Module ECHAM_fields_adj
