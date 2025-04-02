!
!+ GNSS Radio occultation observation operator: Interface to gridded fields
!
MODULE ECHAM_fields
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Interface to griddedfields of H, P, T, and Q,
!   and calculation of gridded field of N and interpolation
!   coefficients for N and Q
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
!  Improve performance of memory management for SX-9,
!  do not nullify array components of some derived types
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Harald Anlauf
!  Make FTRACE_REGIONs depending on macro DISABLE_FTRACE_REGION
! V1_22        2013-02-13 Harald Anlauf
!  rename data_type -> grid_type, add vertical coordinate type
! V1_27        2013-11-08 Harald Anlauf
!  adapt horizontal interpolation and vertical levels for ICON
! V1_29        2014/04/02 Harald Anlauf
!  NM should be SAVEd
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
! Module ECHAM_fields
!
! Reading ECHAM global gridded fields of H, P, T, and Q,
! and calculation of gridded field of N and interpolation
! coefficients for N and Q
!----------------------------------------------------------
! (C) Copyright 1998-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 08 Oct 1998 | Original code.
!   2.0   | 13 Nov 1998 | Surface altitude added.
!   3.0   | 16 Dec 1998 | Vertical_Coordinates in separate
!         |             | module.
!   4.0   | 23 Dec 1998 | Single precision global fields.
!   5.0   | 25 Dec 1998 | Type Profile.
!   6.0   | 17 Mar 1999 | Interpolate_Constituent.
!   7.0   | 20 Mar 1999 | Error check.
!   8.0   | 21 Jun 1999 | Combining phi(p) from ECHAM and MSIS.
!   8.1   | 12 Apr 2000 | GCLat and NI public.
!   9.0   | 03 May 2000 | ECHAM_NGHN.
!  10.0   | 14 May 2001 | Absorption included.
!  11.0   | 07 Apr 2002 | Icosahedral grid included.
!  12.0   | 16 Apr 2002 | Field status.
!  12.1   | 20 Apr 2002 | Allocate_Profile.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double, WorkPr
!
Use Errors, only: &
! Imported Type Definitions:
    Error_Status,  &
! Imported Routines:
    Enter_Callee,  &
    Exit_Callee,   &
    Error
!
Use mo_wmo_tables, only: &
! Imported Parameters:
    WMO6_LATLON,         &
    WMO6_GAUSSIAN,       &
    DWD6_ICOSAHEDRON,    &
    DWD6_ICON

Use mo_exception, only: finish, message

Use ICO_grid, only: &
    gf,       &! global fields
    t_global   ! global fields data type

Use mo_atm_grid, only: &
    VCT_P_HYB  ! hybrid pressure vertical coordinate flag
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! public entities
!
private
public :: gf       ! global fields
public :: t_global ! global fields data type
public :: matrix   ! matrix  data type
public :: profile  ! profile data type

public :: echam_cleanup ! deallocate atmospheric fields
public :: echam_init    ! initialize atmospheric fields
public :: interpolate_refractivity
                        ! Interpolation of global refractivity field from ECHAM

public :: NG1, NM, PStat
public :: FST_NULL, fst_ALLOCATED, FST_INITIALIZED
public :: allocate_profile
public :: echam_ngradn, echam_nghn

!==============================================================================
#ifndef __NULLIFY_POINTER__
! No initialization of pointers in derived types
#define __NULL__
#else
! Use default initialization of pointers with NULL (*VERY* slow on NEC SX-9)
#define __NULL__ => NULL()
#endif
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
! Public Parameters:
!
! --- Error codes
!
Integer, Parameter :: &
   err_Data_Lack     = 5201, &
   err_Size_Mismatch = 5202
!
! --- Field states
!
Integer, Parameter :: &
   fst_Null        = 0,  &
   fst_Allocated   = 1,  &
   fst_Initialized = 2
!
! --- Subgrid dimension
!
Integer, Parameter :: &
   NG1 = 2  ! Number of grid points for lat/lon interpolation:
            ! 2 - linear interpolation
            ! 4 - cubic polynomial interpolation
!
Real(Double), Parameter :: &
   P_up = 1e-5_Double   ! Pressure at uppermost half level [mbar]
!
! Public Type Definitions:
!
Type Profile
   Real(Double) _POINTER :: P(:)   __NULL__ ! Vertical profile at surface point
End Type Profile
!
Type Matrix
   Real(Double) _POINTER :: M(:,:) __NULL__ ! 2D matrix at surface point
End Type Matrix

!
! Public Scalars:
!
Integer :: NM = 0     ! Number of MSIS pressure levels
!
! Public Arrays:
!
! --- Grids and gridded GCM fields
!
!Real(WorkPr), Pointer, Public :: &
!   Hsur(:,:,:),     & ! Surface geopotential [gpm]
!   Psur(:,:,:),     & ! Surface pressure [Pa]
!   T(:,:,:,:),      & ! Temperature [K]
!   Q(:,:,:,:)         ! Specific humidity [kg/kg]
!   XLon(:,:,:),     & ! Longitude grid [deg]
!   XLat(:,:,:)        ! Latitude grid [deg]
!
!Real(Double), Allocatable, Public :: &
!   GCLat(:,:,:)       ! Geocentric latitudes [rad]
!
!Real(Double), Pointer, Public      :: &
!   A(:), B(:)         ! Vertical coordinates
!
! --- Altitudes, refractivity, gradient,
! --- interpolation coefficients
!
!Real(Double), Allocatable, Public  :: &
!   Zsur(:,:,:)        ! Surface altitude [km]
!
Type(Profile), Allocatable, Public :: &
   Z(:,:,:),        & ! Altitudes [km]
   LnN(:,:,:),      & ! Ln of refractive index
   D2N(:,:,:)         ! Interpolation coefficients of Ln(N)
!
Type(Matrix), Allocatable, Private :: &
   LnNI(:,:,:),     & ! Ln of dispersive imaginary part of refractive index
   D2NI(:,:,:)        ! Interpolation coefficients of Ln(N)
!
Real(Double), Allocatable, Private :: &
   Freq(:)            ! Frequency channels [Hz]
!
Integer, Allocatable :: &
   PStat(:,:,:)       ! Dynamical field status
!
! Index order of model fields on rectangular grid:
!   1. Longitude (0 --> 360)
!   2. Latitude  (North --> South)
!   3. Stub index = 1
!   4. Level     (highest --> lowest)
!
! Index order of model fields on icosahedral grid:
!   1. j1 inside diamond
!   2. j2 inside diamond
!   3. jd - diamond number
!   4. Level     (highest --> lowest)
!----------------------------------------------------------
!
Contains


!==========================================================
Subroutine deallocate_fields (g)
type (t_global) ,intent(inout) :: g
  type (t_global) :: default
  !------------------------------------------------------------
  ! deallocate components of global atmospheric field data type
  ! (also clears all metadata)
  !------------------------------------------------------------
  if (associated (g%i))      deallocate (g%i)
  if (associated (g%s))      deallocate (g%s)
  if (associated (g%t))      deallocate (g%t)
  if (associated (g%q))      deallocate (g%q)
  if (associated (g%qx))     deallocate (g%qx)
  if (associated (g%p))      deallocate (g%p)
  if (associated (g%ph))     deallocate (g%ph)
  if (associated (g%glon))   deallocate (g%glon)
  if (associated (g%glat))   deallocate (g%glat)
! if (associated (g%a))      deallocate (g%a)
! if (associated (g%b))      deallocate (g%b)
  if (associated (g%b_ad))   deallocate (g%b_ad)
  if (associated (g%z))      deallocate (g%z)
  if (associated (g%xnglob)) nullify    (g%xnglob)
  if (associated (g%marr))   nullify    (g%marr)
  if (associated (g%grid))   nullify    (g%grid)
  g = default
end Subroutine deallocate_fields
!------------------------------------------------------------------------------
!Subroutine allocate_fields (g, n, grid_type, lb, ub, nz, nd, ni, ni2, ni3)
!type (t_global) ,intent(inout) :: g         ! variable to allocate
!integer         ,intent(in)    :: n         ! number of columns to allocate
!integer         ,intent(in)    :: grid_type ! grid type
!integer         ,intent(in)    :: lb (3)    ! lower bounds of global field
!integer         ,intent(in)    :: ub (3)    ! upper bounds (x1,x2,diam)
!integer         ,intent(in)    :: nz        ! number of levels
!integer         ,intent(in)    :: nd
!integer         ,intent(in)    :: ni
!integer         ,intent(in)    :: ni2
!integer         ,intent(in)    :: ni3
!  !------------------------------------------------------------
!  ! allocate components of global atmospheric field data type
!  !------------------------------------------------------------
!  call deallocate_fields (g)
!  g% grid_type = grid_type
!  g% n         = n
!  g% lbg       = lb
!  g% ubg       = ub
!  g% nz        = nz
!  g% nd        = nd
!  g% ni        = ni
!  g% ni2       = ni2
!  g% ni3       = ni3
!  allocate (g% i    (lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
!  allocate (g% s    (n))
!  allocate (g% t    (nz,n))
!  allocate (g% q    (nz,n))
!  allocate (g% p    (nz,n))
!  allocate (g% a    (nz+1))
!  allocate (g% b    (nz+1))
!  allocate (g% b_ad (nz+1))
!  select case (grid_type)
!  case (WMO6_GAUSSIAN, WMO6_LATLON)
!    allocate (g% glon (lb(1):ub(1)))
!    allocate (g% glat (lb(2):ub(2)))
!  case (DWD6_ICOSAHEDRON)
!  end select
!End Subroutine allocate_fields
!------------------------------------------------------------------------------
Subroutine ECHAM_cleanup

!--- 0.2. Clearing old data

! this part was moved from Subroutine ECHAM_init to ECHAM_cleanup

integer :: i1, i2, id, ErrCode

call Deallocate_fields (gf)

!Deallocate(Hsur,  Stat = ErrCode)
!Deallocate(Psur,  Stat = ErrCode)
!Deallocate(T,     Stat = ErrCode)
!Deallocate(Q,     Stat = ErrCode)
!Deallocate(XLon,  Stat = ErrCode)
!Deallocate(XLat,  Stat = ErrCode)
!Deallocate(GCLat, Stat = ErrCode)
!Deallocate(Zsur,  Stat = ErrCode)

If (Allocated(Z)) then

FTRACE_BEGIN("ECHAM_cleanup:Deallocate")
   Do ID=LBound(Z,3),UBound(Z,3)
      Do I2=LBound(Z,2),UBound(Z,2)
         if (any (PStat(:,I2,ID) /= fst_Null)) then
            Do I1=LBound(Z,1),UBound(Z, 1)
               if (PStat(I1, I2, ID) /= fst_Null) then
                  Deallocate(Z  (I1, I2, ID)%P, Stat = ErrCode)
                  Deallocate(LnN(I1, I2, ID)%P, Stat = ErrCode)
                  Deallocate(D2N(I1, I2, ID)%P, Stat = ErrCode)
               end if
            End Do
         end if
      End Do
   End Do

   If (Allocated(Freq)) then
      Do ID=LBound(Z,3),UBound(Z,3)
         Do I2=LBound(Z,2),UBound(Z,2)
            if (any (PStat(:,I2,ID) /= fst_Null)) then
               Do I1=LBound(Z,1),UBound(Z, 1)
                  if (PStat(I1, I2, ID) /= fst_Null) then
                     Deallocate(LnNI(I1, I2, ID)%M, Stat = ErrCode)
                     Deallocate(D2NI(I1, I2, ID)%M, Stat = ErrCode)
                  end if
               End Do
            end if
         End Do
      End Do
   End If
FTRACE_END  ("ECHAM_cleanup:Deallocate")

End If

Deallocate(PStat, Stat = ErrCode)
Deallocate(Z,     Stat = ErrCode)
Deallocate(LnN,   Stat = ErrCode)
Deallocate(D2N,   Stat = ErrCode)
Deallocate(LnNI,  Stat = ErrCode)
Deallocate(D2NI,  Stat = ErrCode)
Deallocate(Freq,  Stat = ErrCode)

!Deallocate(A, B,  Stat = ErrCode)

!Nullify(Hsur, Psur, T,  Q,  XLon, XLat)

end Subroutine ECHAM_cleanup
!------------------------------------------------------------------------------
Subroutine ECHAM_Init ( &
!  Pathnames, & ! <-- Pathmnames of data files
   ErrStat,   & ! <-> Error status
   XFreq)       ! <~~ Frequency channels [Hz]
!
! Reading global fields of surface geopotential,
! surface pressure, temperature and humidity,
! calculation of refractive index and interpolation
! coefficients of refractivie index and humidity.
!----------------------------------------------------------
! Method:
!    Described in:
!    1. Deutsches Klimarechenzentrum, Technical Report No.6,
!       The ECHAM3 Atmospheric General Circulation Model,
!       Revision 2, Hamburg, 1993.
!    2. M. E. Gorbunov, S. V. Sokolovsky, L. Bengtsson, Space
!       refractive tomography of the atmosphere: modeling of
!       the direct and inverse problems, Max-Planck Institut
!       fuer Meteorologie, Report No. 210, Hamburg, 1996.
!----------------------------------------------------------
! (C) Copyright 1998-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 08 Oct 1998 | Original code.
!   2.0   | 13 Nov 1998 | Surface altitutde.
!   2.1   | 20 Feb 1999 | Relaxation to MSIS profile.
!   2.2   | 20 Mar 1999 | Check for Allocated(Z)
!         |             | when clearing old data.
!   3.0   | 14 May 2001 | Absorption included.
!   4.0   | 07 Apr 2002 | Reading data on icosahedral
!         |             | grid included.
!   4.1   | 16 Apr 2002 | Field status.
!----------------------------------------------------------
! Modules used:
!
!!Use ECHAM_IO, only: &
!!! Imported Routines:
!!    Get_GRIB
!
Use Earth, only: &
! Imported Routines:
    GCLat_from_GDLat
!
Use MSIS, only: &
! Imported Routines:
    MSIS_Init,       &
    MSIS_Num_Levels
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
!Character(Len=*), Intent(In)  :: &
!   Pathnames(1:)  ! Pathnames of data files
!
! InOut arguments:
!
Type(Error_Status), Pointer   :: &
   ErrStat        ! Error status
!
! Optional input arguments:
!
Real(Double), Intent(In), Optional :: &
   XFreq(:)       ! Frequency channels [Hz]
!----------------------------------------------------------
! Local Scalars:
!
!Integer  :: Date     ! Data set date
Integer   :: Mon      ! Data set month
!Integer  :: Time     ! Data set time
!Integer  :: IFile    ! Data file number
Integer   :: N1S, N1E ! 1st grid index limits
!Integer  :: N1       ! 1nd grid dimension
Integer   :: N2       ! 2nd grid dimension
Integer   :: ND       ! Number of diamonds
!Integer   :: NLev     ! Number of levels
!Integer   :: I1, I2   ! Grid indices
!Integer   :: ID       ! Diamond index
Integer   :: NC       ! Number of frequency channels
Integer   :: ErrCode  ! Error code
integer   :: ix
!----------------------------------------------------------
! Global variables used:
!
!   Hsur(:,:,:)    ! Surface geopotential [gpm]
!   Psur(:,:,:)    ! Surface pressure [Pa]
!   T(:,:,:,:)     ! Temperature [K]
!   Q(:,:,:,:)     ! Specific humidity [kg/kg]
!   XLon(:,:,:)    ! Longitude grid [deg]
!   XLat(:,:,:)    ! Latitude grid [deg]
!   GCLat(:,:,:)   ! Geocentric latitudes [rad]
!   Zsur(:,:,:)    ! Surface altitude [km]
!   Z(:,:,:,:)     ! Altitudes [km]
!   LnN(:,:,:,:)   ! Ln of refractive index
!   D2N(:,:,:,:)   ! Interpolation coefficients of Ln(N)
!   A(:), B(:)     ! Vertical coordinates
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------


!--- 0.1. Error status

Call Enter_Callee &
  ('ECHAM_Init',  & ! <-- User routine
   ErrStat)         ! <-> Pointer to callee status


!----------------------------------------------------------
! 1. READING DATA
!----------------------------------------------------------

!!Read: Do IFile = 1, Size(Pathnames)
!!   Call Get_GRIB &
!!     (PathNames(IFile),  & ! <-- Pathname of data file
!!      Grid_Type,         & ! --> Data representation type
!!      Hsur,              & ! --> Surface geopotential [gpm]
!!      Psur,              & ! --> Surface pressure [Pa]
!!      T,                 & ! --> Temperature [K]
!!      Q,                 & ! --> Specific humidity [kg/kg]
!!      XLon,              & ! --> Longitude grid [deg]
!!      XLat,              & ! --> Latitude grid [deg]
!!      A,                 & ! --> Vertical coordinate parameters A
!!      B,                 & ! --> Vertical coordinate parameters B
!!      Date,              & ! --> Field date
!!      Time,              & ! --> Field time
!!      ErrStat)             ! --> Error code
!!   If (Error(ErrStat)) then
!!      Call Exit_Callee(ErrStat)
!!      Return
!!   End If
!!End Do Read


!----------------------------------------------------------
! 2. CHECKING FOR PRESENSE OF ALL DATA ARRAYS
!----------------------------------------------------------

!If (.not. All((/               &
!            Associated(XLon),  &
!            Associated(XLat),  &
!            Associated(Hsur),  &
!            Associated(Psur),  &
!            Associated(T),     &
!            Associated(Q) /))) then
!   Call Exit_Callee &
!     (ErrStat,           & ! <-> Pointer to callee status
!      err_Data_Lack,     & ! <~~ User error code
!      0,                 & ! <~~ System error code
!      'Missing fields')    ! <~~ Error description
!   Return
!End If


!----------------------------------------------------------
! 3. CHECKING FOR SIZE FIT
!----------------------------------------------------------

N1S  = gf% lbg(1)
N1E  = gf% ubg(1)
!N1  = gf% ubg(1) - gf% lbg(1) + 1
N2   = gf% ubg(2) - gf% lbg(2) + 1
ND   = gf% ubg(3) - gf% lbg(3) + 1
!NLev = gf% nz

!If (.not. All((/                                     &
!        All(Shape(Hsur) == (/ N1, N2, ND /)),        &
!        All(Shape(Psur) == (/ N1, N2, ND /)),        &
!        All(Shape(T)    == (/ N1, N2, ND, NLev /)),  &
!        All(Shape(Q)    == (/ N1, N2, ND, NLev /)) /))) then
!   Call Exit_Callee &
!     (ErrStat,           & ! <-> Pointer to callee status
!      err_Size_Mismatch, & ! <~~ User error code
!      0,                 & ! <~~ System error code
!      'Size mismatch')     ! <~~ Error description
!   Return
!End If


!----------------------------------------------------------
! 4. CALCULATION OF VERTICAL COORDINATES
!----------------------------------------------------------

!If ((.not. Associated(A)) .and. (.not. Associated(B))) then
!
!   Allocate(A(0:NLev), B(0:NLev), Stat=ErrCode)
!   If (ErrCode /= 0) then
!      Call Exit_Callee &
!        (ErrStat,           & ! <-> Pointer to callee status
!         0,                 & ! <~~ User error code
!         ErrCode,           & ! <~~ System error code
!         'A/B allocation')    ! <~~ Error description
!      Return
!   End If
!
!   Call Vertical_Coordinates(A, B, ErrStat)
!   If (Error(ErrStat)) then
!      Call Exit_Callee(ErrStat)
!      Return
!   End If
!
!End If


!----------------------------------------------------------
! 5. CALCULATION OF GEOCENTRIC LATITUDES
!----------------------------------------------------------

!Select Case (Grid_Type)
!
!   Case (WMO6_GAUSSIAN, WMO6_LATLON)
!
!      Allocate(GCLat(1, N2, 1), Stat=ErrCode)
!
!   Case (DWD6_ICOSAHEDRON)
!
!      Allocate(GCLat(N1S:N1E, N2, ND), Stat=ErrCode)
!
!End Select

!If (ErrCode /= 0) then
!   Call Exit_Callee  &
!     (ErrStat,            & ! <-> Pointer to callee status
!      0,                  & ! <~~ User error code
!      ErrCode,            & ! <~~ System error code
!      'GCLAt allocation')   ! <~~ Error description
!   Return
!End If

do ix = 1, gf% n
  gf% s(ix)% GCLat = GCLat_from_GDLat(Real(gf% s(ix)% XLat, Double))
End Do


!----------------------------------------------------------
! 6. MEMORY ALLOCATION FOR SURFACE ARRAYS
!----------------------------------------------------------


!--- 6.1. Memory allocation

FTRACE_BEGIN("ECHAM_Init:Allocate")

Allocate(PStat(N1S:N1E, N2, ND))
PStat(:,:,:) = fst_Null

!Allocate(Zsur(N1S:N1E, N2, ND),  &
Allocate(   Z(N1S:N1E, N2, ND),  &
          LnN(N1S:N1E, N2, ND),  &
          D2N(N1S:N1E, N2, ND),  &
          Stat = ErrCode)

If (ErrCode /= 0) then
   Call Exit_Callee  &
     (ErrStat,                   & ! <-> Pointer to callee status
      0,                         & ! <~~ User error code
      ErrCode,                   & ! <~~ System error code
      'Surface field allocation')  ! <~~ Error description
   Return
End If

FTRACE_END  ("ECHAM_Init:Allocate")

!--- 6.2. Optional field allocation

If (Present(XFreq)) then

   NC = Size(XFreq)

   Allocate(Freq(NC))
   Freq(:) = XFreq(:)

   Allocate(LnNI(N1S:N1E, N2, ND),  &
            D2NI(N1S:N1E, N2, ND),  &
            Stat=ErrCode)

   If (ErrCode /= 0) then
      Call Exit_Callee  &
        (ErrStat,                   & ! <-> Pointer to callee status
         0,                         & ! <~~ User error code
         ErrCode,                   & ! <~~ System error code
         'Surface field allocation')  ! <~~ Error description
      Return
   End If

Else

   NC = 0

End If


!--- 6.3. Nullifying vertical profiles

!Do I1=N1S,N1E
!   Do I2=1,N2
!      Do ID=1,ND
!         if(allocated(z))   Nullify(Z  (I1, I2, ID)%P)
!         if(allocated(lnn)) Nullify(LnN(I1, I2, ID)%P)
!         if(allocated(d2n)) Nullify(D2N(I1, I2, ID)%P)
!      End Do
!   End Do
!End Do

!!If (NC <> 0) then
!If (NC /= 0) then
!   Do I1=N1S,N1E
!      Do I2=1,N2
!         Do ID=1,ND
!            Nullify(LnNI(I1, I2, ID)%M)
!            Nullify(D2NI(I1, I2, ID)%M)
!         End Do
!      End Do
!   End Do
!End If


!----------------------------------------------------------
! 7. MSIS INITIALIZATION
!----------------------------------------------------------

Mon = Mod (gf% yyyymmdd/100, 100)

Call MSIS_Init &
  (Mon,      & ! <-- Month
   ErrStat)    ! --> Error code

If (Error(ErrStat)) then
   Call Exit_Callee(ErrStat)
   call message ('echam_init','error in MSIS_Init')
   Return
End If


!----------------------------------------------------------
! 8. CALCULATION OF NUMBER OF MSIS PRESSURE LEVELS
!----------------------------------------------------------

!NM = MSIS_Num_Levels (1e-2_Double * gf% A(1), P_up)
NM = MSIS_Num_Levels (1e-2_Double * gf% p0_msis, P_up)


Call Exit_Callee(ErrStat)


End Subroutine ECHAM_Init



!==========================================================
Subroutine Allocate_Profile &
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

Allocate(  Z(I1, I2, ID)%P(-NM+2:NLev),     &
         LnN(I1, I2, ID)%P(-NM+2:NLev),     &
         D2N(I1, I2, ID)%P(-NM+2:NLev))

PStat(I1, I2, ID) = fst_Allocated


End Subroutine Allocate_Profile



!==========================================================
Subroutine Make_Refractivity &
  (Hsur,   & ! <-- Surface geopotential [gpm]
   Gundu,  & ! <-- Geoid undulation     [m]
   Psur,   & ! <-- Surface pressure [Pa]
   Ph,     & ! <-- Pressure at half levels [Pa]
   T,      & ! <-- Temperature [K]
   Q,      & ! <-- Specific humidity [kg/kg]
   Qx,     & ! <-- Water load [kg/kg]
   G,      & ! <-- Geodetic coordinates
   GCLat,  & ! <-- Geocentric latitude [rad]
   LnN,    & ! --> Ln of refractive index
   Zsur,   & ! --> Surface altitude [km]
   Z,      & ! --> Altitudes of model levels [km]
   LnNI)     ! ~~> Ln of dispersive imaginary part of refractive index
!
! Calculation of refractivity profile and geometrical
! altitudes of pressure levels from surface  geopotential,
! surface pressure, and profiles of temperature and humidity.
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
! (C) Copyright 1998-2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Oct 1998 | Original code.
!   2.0   | 13 Nov 1998 | Surface altitude.
!   3.0   | 01 Jun 1999 | Combining phi_ECHAM(P) and phi_MSIS(P).
!   4.0   | 14 May 2001 | Absorption included.
!----------------------------------------------------------
! Modules used:
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic,  &
! Imported Parameters:
    g_ave,     &
! Imported Routines:
    Alt_from_Geop
!
Use Occ_Meteoprofiles, only: &
! Imported Parameters:
    Rd, Eps,    &
! Imported Routines:
    N_from_TPQ,            &
    T_from_NPQ,            &
    N_from_TPQ_AL2011
!
Use MSIS, only: &
! Imported Routines:
    MSIS_Pressure_Levels,  &
    MSIS_Geop,             &
    MSIS_Refractivity
!
Use Occ_MPM93, only: &
! Imported Routines:
    N93AIR
!
Use CIPM_2007, only: &
    Zeta
!
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(WorkPr), Intent(In)   :: &
   Hsur           ! Surface geopotential [gpm]
!
Real(WorkPr), Intent(In)   :: &
   Gundu          ! Geoid undulation [m]
!
Real(WorkPr), Intent(In)   :: &
   Psur           ! Surface pressure [Pa]
!
Real(WorkPr), Intent(In)   :: &
   Ph(0:)         ! Pressure at half levels [Pa]
!
Real(WorkPr), Intent(In)   :: &
   T(1:)          ! Temperature profile [K]
!
Real(WorkPr), Intent(In)   :: &
   Q(1:)          ! Specific humidity [kg/kg]
!
Real(WorkPr), Intent(In)   :: &
   Qx(1:)         ! Water load [kg/kg]
!
Type(Geodetic), Intent(In) :: &
   G              ! Geodetic coordinates
!
Real(Double), Intent(In)   :: &
   GCLat          ! Geocentric latitude [rad]
!
! Output arguments:
!
Real(Double), Intent(Out)  :: &
   LnN(-NM+2:)    ! Ln of refractive index
!
Real(Double), Intent(Out)  :: &
   Zsur           ! Surface altitude [km]
!
Real(Double), Intent(Out)  :: &
   Z(-NM+2:)      ! Altitudes of model levels [km]
!
! Output optional arguments:
!
Real(Double), Intent(Out), Optional :: &
   LnNI(-NM+2:,:) ! Ln of dispersive imaginary part of refractive index
!----------------------------------------------------------
! Local Scalars:
!
Integer         :: NLev    ! Number of levels
Integer         :: i       ! Level index
Real(Double)    :: alpha   ! Log(p(i+1)/p(i))
Type(Geodetic)  :: GM      ! Point in vertical profile
Real(Double)    :: DLn     ! Correction of Ln(P_MSIS)
Integer         :: IC      ! Channel number
Complex(Double) :: ZN      ! Dispersive part of refractivity [N-units]
Real(Double)    :: NN      ! Non-dispersive part of refractivity [N-units]
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
Real(Double) :: &
   TM(-NM+2:0)       ! MSIS level temperature [K]
!
Real(Double), dimension(size(T)) :: &
   DZhalf,         & ! Correction from compressibility
   DZfull            ! for half and full levels
!
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


!----------------------------------------------------------
! 1. CALCULATION OF ECHAM HALF AND FULL LEVELS
!----------------------------------------------------------

!--- 1.1. Calculation of ECHAM half levels

!Phalf(0:NLev) = gf% A(0:NLev) + gf% B(0:NLev)*Psur
Phalf(0:NLev) = Ph(0:NLev)
Tvirt(1:NLev) = T(1:NLev)*(1 + Eps*Q(1:NLev) - Qx(1:NLev))

if (gf% vctype == VCT_P_HYB) then
Pfull(1:NLev) =      (Phalf(0:NLev-1) + Phalf(1:NLev)) * 0.5_Double
else
Pfull(1:NLev) = sqrt (Phalf(0:NLev-1) * Phalf(1:NLev))
end if

if (gf% vctype == VCT_P_HYB) then
   Phalf(0) = Phalf(1) / 4      ! ECHAM/GME
end if


Hhalf(NLev) = Hsur/g_ave + Gundu
Zsur        = Alt_from_Geop(1e-3_Double*Hhalf(NLev), GCLat)

DZhalf = 0._Double
DZfull = 0._Double
if (gf% ref_model == 2) then
   do i = 1, NLev
      DZhalf(i) = Zeta (Phalf(i),T(i),Q(i)) - Zeta (Phalf(i-1),T(i),Q(i))
      DZfull(i) = Zeta (Phalf(i),T(i),Q(i)) - Zeta (Pfull(i  ),T(i),Q(i))
   end do
end if


Do i=NLev-1,0,-1
   alpha    = Log(Phalf(i+1)/Phalf(i))
   Hhalf(i) = Hhalf(i+1) + Rd/g_ave*Tvirt(i+1) * (alpha + DZhalf(i+1))
End Do

if (gf% vctype == VCT_P_HYB) then
   Phalf(0) = Ph(0)             ! ECHAM/GME: restore pressure at model top
end if


!--- 1.2. Calculation of ECHAM full levels

Do i = 1,NLev

!  Pfull(i) = (Phalf(i-1) + Phalf(i))/2

   If (i == 1 .and. gf% vctype == VCT_P_HYB) then
      alpha = Log(2.0_Double)
   Else
      alpha = Log(Phalf(i)/Pfull(i))
   End If

   Hfull(i) = Hhalf(i) + Rd/g_ave*Tvirt(i) * (alpha + DZfull(i))
   Z(i)     = Alt_from_Geop(1d-3*Hfull(i), GCLat)

End Do

!----------------------------------------------------
! "Scaling" of refractivity: N = (n-1) -> ref_scale*N
!----------------------------------------------------
LnScale = 0._Double
if (gf% ref_scale /= 1._Double) LnScale = log (gf% ref_scale)

select case (gf% ref_model)
case (1)
   LnN(1:NLev) = Log (N_from_TPQ                     &
                        (Real(T(1:NLev), Double),    &
                         Pfull (1:NLev)*1e-2_Double, & ! <-- hPa
                         Real(Q(1:NLev), Double)))
case (2)
   LnN(1:NLev) = Log (N_from_TPQ_AL2011              &
                        (Real(T(1:NLev), Double),    &
                         Pfull (1:NLev),             & ! <-- Pa
                         Real(Q(1:NLev), Double)))
end select

LnN(1:NLev) = LnN(1:NLev) + LnScale

!----------------------------------------------------------
! 2. CALCULATION OF MSIS HALF AND FULL LEVELS
!----------------------------------------------------------

!--- 2.0. Unit conversion

Phalf(0:NLev) = 1e-2_Double*Phalf(0:NLev) ! Pa  --> mb
Pfull(1:NLev) = 1e-2_Double*Pfull(1:NLev) ! Pa  --> mb
Hhalf(0:NLev) = 1e-3_Double*Hhalf(0:NLev) ! gpm --> gpkm
Hfull(1:NLev) = 1e-3_Double*Hfull(1:NLev) ! gpm --> gpkm


!--- 2.1. Calculation of MSIS pressure levels

Call MSIS_Pressure_Levels &
! (1e-2_Double * gf% A(1) , & ! <-- Basic pressure level [mb]
  (1e-2_Double * Ph(1)    , & ! <-- Basic pressure level [mb]
   Phalf(-NM+1:0))            ! --> MSIS pressure levels [mb]

Do i=-NM+2,1
   Pfull(i) = (Phalf(i-1) + Phalf(i))/2
End Do


!--- 2.2. Calculation of MSIS full and half(0) level altitudes

PM(-NM+2:0) = Pfull(-NM+2:0)
PM(1)       = Phalf(0)

Call MSIS_Geop &
  (G,            & ! <-- Geodetic coordinates
   PM(-NM+2:1),  & ! <-- Pressure levels [mb]
   HM(-NM+2:1),  & ! --> Geopotential heights [gpkm]
   ZM(-NM+2:1))    ! --> Altitudes [km]

Hfull(-NM+2:0) = HM(-NM+2:0)
Hhalf(0)       = HM(1)
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
Z(1)     = Alt_from_Geop(Hfull(1), GCLat)
LnN(1)   = Log(N_from_TPQ               &
                 (Real(T(1), Double),   &
                  Pfull(1),             &
                  Real(Q(1), Double)))


!--- 2.5. Correction of MSIS profile

DLn          = Log(Phalf(1)/Phalf(0)) - &
               1e3_Double*(g_ave/Rd)*(Hhalf(0) - Hhalf(1))/Tvirt(1)
LnN(-NM+2:0) = LnN(-NM+2:0) + DLn
PM(-NM+2:0)  = PM(-NM+2:0)*Exp(DLn)


!----------------------------------------------------------
! 3. CALCULATION OF ABSORPTION PROFILE
!----------------------------------------------------------

If (Present(LnNI) .and. Allocated(Freq)) then


   !--- 3.1. Absorption computation at MSIS levels

   Do i=-NM+2,0
      TM(i) = T_from_NPQ &
        (Exp(LnN(i)),     & ! <-- Refractivity [dimensionless]
         PM(i),           & ! <-- Pressure [mb]
         0.0_Double)        ! <-- Specific humidity [kg/kg]
      Do IC=1,Size(Freq)
         Call N93AIR &
           (Freq(IC),     & ! <-- Frequency [Hz]
            PM(i),        & ! <-- Pressure [mbar]
            TM(i),        & ! <-- Temperature [K]
            0.0_Double,   & ! <-- Specific humidity [kg/kg]
            ZN,           & ! --> Dispersive part of refractivity [N-units]
            NN)             ! --> Non-dispersive part of refractivity [N-units]
         LnNI(i,IC) = Log(Aimag(ZN))
      End Do
   End Do


   !--- 3.2. Absorption computation at ECHAM levels

   Do i=1,NLev
      Do IC=1,Size(Freq)
         Call N93AIR &
           (Freq(IC),            & ! <-- Frequency [Hz]
            Pfull(i),            & ! <-- Pressure [mbar]
            Real(T(i), Double),  & ! <-- Temperature [K]
            Real(Q(i), Double),  & ! <-- Specific humidity [kg/kg]
            ZN,                  & ! --> Dispersive part of refractivity [N-units]
            NN)                    ! --> Non-dispersive part of refractivity [N-units]
         LnNI(i,IC) = Log(Aimag(ZN))
      End Do
   End Do

End If


End Subroutine Make_Refractivity



!==========================================================
Subroutine Interpolate_Refractivity &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NG,       & ! ~~> Interpolated dN/d(alt,lat,lon)
   NH,       & ! ~~> Interpolated hessian matrix of N
   NIP)        ! ~~> Interpolated Im(N)
!
! Interpolation of global refractivity field from ECHAM.
!----------------------------------------------------------
! Method:
!   1. Vertical spline interpolation to given altitude.
!   2. Summation of 2D Lon/Lat interpolation series:
!   f(z,Lat,Lon) = sum WG f(z,Lat ,Lon )
!                   i    i       i    i
!   Index i enumerates points of the interpolation
!   subgrid. For rectangular grid, weigting function WG
!   is calculated for polynomial interpolation with the
!   polynomial power defined by the interpolation subgrid
!   dimension NI. For icosahedral grid WG is defined as
!   symplectic weights of the vertices of surrounding
!   triangle.
!----------------------------------------------------------
! (C) Copyright 1998-2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 09 Oct 1998 | Original code.
!   2.0   | 13 Nov 1998 | Latitudinal and longitudinal
!         |             | interpolation separated,
!         |             | use of surface altitude.
!   3.0   | 17 Nov 1998 | d/d(alt,lat,lon)
!   4.0   | 25 Dec 1998 | Separate subroutines for
!         |             | refractivity and humidity.
!   5.0   | 20 Feb 1999 | Calculation of N instead of ln(N)
!   5.1   | 13 Apr 2000 | Vertical interpolation in a
!         |             | separate subroutine.
!   6.0   | 14 May 2001 | Optional interpolation of Im(N).
!   7.0   | 07 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!----------------------------------------------------------
! Modules used:
!
Use ECHAM_grid, only: &
! Imported Parameters:
    Longitude_Interpolation,  &
    Latitude_Interpolation
!
Use ICO_grid, only: &
! Imported Array Variables:
    nspoke
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
Real(Double), Optional, Intent(Out) :: &
   NG(3)         ! Interpolated gradient dN/d(alt,lat,lon)
!
Real(Double), Optional, Intent(Out) :: &
   NH(3,3)       ! Interpolated hessian
                 ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))N
!
Real(Double), Optional, Intent(Out) :: &
   NIP(:)        ! Interpolated Im(N)
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NGP     ! Subgrid dimension
Integer        :: IGP     ! Subgrid index
Integer        :: KLon    ! Longitude subgrid index
Integer        :: KLat    ! Latitude subgrid index
Integer        :: m1, m2  ! Interpolation triangle indices
Real(Double)   :: FP      ! Interpolated Ln(N)
Integer        :: i       ! Work index
Integer        :: NC      ! Number of frequency channels
Integer        :: IC      ! Channel number
Integer        :: Stat    ! Error status
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
   WLat1(NG1),        & ! 1st derivative of WLon
   WLat2(NG1)           ! 1st derivative of WLon
!
! --- Interpolation on icosahedral grid
!
Real(Double) :: &
   PRLon(1),          & ! Longitude of interpolation point
   PRLat(1),          & ! Latitude  of interpolation point
   SWS(3,1)             ! Symplectic weights for triangle vertices
Integer :: &
   KIDX(4,1)            ! Interpolation triangle index
!
! --- Interpolation on ICON grid
!
integer :: jl(3,1),   & ! Line  indices
           jb(3,1),   & ! Block indices
           npr (1)      ! number of neighbour points returned
!
! --- Interpolation indices and weights
!
Integer, Dimension(gf% ngp,3) :: &
   IDX                  ! Subgrid indices [point, index]
Real(Double), Dimension(gf% ngp,0:2,0:2)  :: &
   W_                   ! Weighting function
Real(Double), Dimension(gf% ncol,0:2,0:2) :: &
   WG                   ! Weighting function
                        ! [Point, Lon.derivative, Lat.derivative]
!
! --- Interpolated values
!
Real(Double), Allocatable :: &
   ZNmin(:),          & ! Minimum model Z on subgrid
   ZNmax(:),          & ! Maximum model Z on subgrid
   FV(:),             & ! Vertically interpolated Ln(N) on subgrid
   FZ(:),             & ! Vertical gradient on subgrid
   F2Z(:)               ! Second vertial derivative on subgrid
Real(Double), Allocatable :: &
   FVI(:,:),          & ! Vertically interpolated Ln(Im(N)) on subgrid
   FIP(:)               ! Interpolated Ln(Im(N))
Real(Double) :: &
   FG(3),             & ! Interpolated dln(N)/d(alt,lat,lon)
   FH(3,3)              ! Interpolated hessian
                        ! (d/d(alt,lat,lon))x(d/d(alt,lat,lon))ln(N)
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


!--- 1.1. Subgrid dimension

!Select Case (gf% Grid_Type)
!   Case (WMO6_GAUSSIAN, WMO6_LATLON)
!      NGP = NG1**2
!   Case (DWD6_ICOSAHEDRON, DWD6_ICON)
!      NGP = 3
!   case default
!      write (0,*) 'Interpolate_Refractivity: invalid grid type',gf% Grid_Type
!      call finish('Interpolate_Refractivity','invalid grid type')
!End Select

NGP = gf% ncol  ! Number of actually used column(s)


!--- 1.2. Memory allocation

Allocate  &
  (ZNmin(NGP), ZNmax(NGP),        &
   FV(NGP),    FZ(NGP),    F2Z(NGP))

If (Allocated(Freq) .and. Present(NIP)) then
   NC = Size(Freq)
   Allocate(FVI(NGP,NC), FIP(NC))
End If


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

      PRLon(1) = PLon
      PRLat(1) = PLat

      Call Setup_Interpolation &
        (1,          & ! <-- Number of interpolation points
         PRLon,      & ! <-- Point longitudes [deg]
         PRLat,      & ! <-- Point latitudes [deg]
         KIDX,       & ! --> Index array for interpolation
         SWS)          ! --> Symplectic weights for triangle vertices

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

      do igp = 1, 3
         pe_i1_i2_id  = gf% marr (:,IDX(igp,1),IDX(igp,2),IDX(igp,3))
         IDX(igp,1:3) = pe_i1_i2_id(2:4)
!        IDX(igp,4)   = pe_i1_i2_id(1)
      end do

   !--- 2.?. Interpolation on ICON triangular grid

   Case (DWD6_ICON)
      PRLon(1) = PLon
      PRLat(1) = PLat
      call search_icon_global (gf% grid% icongrid, &! <-- Interpolation metadata
                               PRLon,              &! <-- Point longitude [deg]
                               PRLat,              &! <-- Point latitude  [deg]
                               jl, jb,             &! --> Line and block indices
                               SWS,                &! --> Interpolation weights
                               NPR                 )! --> Number of neighbours
      W_(:,:,:)  = 0.0_Double
      W_(:,0,0)  = SWS(1:3,1)
      idx(1:3,1) = jl (1:3,1)
      idx(1:3,2) = jb (1:3,1)
      idx(1:3,3) = 1
      do igp = 1, 3
         pe_i1_i2_id  = gf% marr (:,IDX(igp,1),IDX(igp,2),IDX(igp,3))
         IDX(igp,1:3) = pe_i1_i2_id(2:4)
!        IDX(igp,4)   = pe_i1_i2_id(1)
      end do

!   !--- 2.?. Single column: no horizontal interpolation.
!
!   Case (DWD6_GRID_NONE)
!      W_(:,:,:)  = 0.0_Double   ! Trivial weights
!      W_(1,0,0)  = 1.0_Double
!      IDX(1,1:3) = 1            ! Dummy indices

   case default
      write (0,*) 'Interpolate_Refractivity: invalid grid type',gf% Grid_Type
      call finish('Interpolate_Refractivity','invalid grid type')
End Select

select case (gf% hint_mode)
case (0)
   !--------------------------
   ! Nearest model column only
   !--------------------------
   i = maxloc (W_(:,0,0), dim=1)
!print *, "#### i=", i, W_(i,0,0)
   IDX(1,1:3) = IDX(i,1:3)
   W_(:,:,:)  = 0.0_Double      ! Trivial weights
   W_(1,0,0)  = 1.0_Double
end select

WG(1:NGP,:,:) = W_(1:NGP,:,:)   ! Weights of actually used column(s)

!----------------------------------------------------------
! 3. VERTICAL INTERPOLATION
!----------------------------------------------------------

Do IGP = 1,NGP
   If (Allocated(Freq) .and. Present(NIP)) then
      Call Vertical_Interpolation_LnN &
        (IDX(IGP,1),         & ! <-- 1st grid index
         IDX(IGP,2),         & ! <-- 2nd grid index
         IDX(IGP,3),         & ! <-- Diamond index
         ZP,                 & ! <-- Altitude [km]
         ZNmin(IGP),         & ! --> Minimum model Z for this lon/lat
         ZNmax(IGP),         & ! --> Maximum model Z for this lon/lat
         FV(IGP),            & ! --> Interpolated Ln(N)
         FZ(IGP),            & ! --> Interpolated dLn(N)/dZ
         F2Z(IGP),           & ! --> Interpolated d2Ln(N)/dZ2
         FVI(IGP,:))           ! ~~> Interpolated Ln(NI)
   Else
      Call Vertical_Interpolation_LnN &
        (IDX(IGP,1),         & ! <-- 1st grid index
         IDX(IGP,2),         & ! <-- 2nd grid index
         IDX(IGP,3),         & ! <-- Diamond index
         ZP,                 & ! <-- Altitude [km]
         ZNmin(IGP),         & ! --> Minimum model Z for this lon/lat
         ZNmax(IGP),         & ! --> Maximum model Z for this lon/lat
         FV(IGP),            & ! --> Interpolated Ln(N)
         FZ(IGP),            & ! --> Interpolated dLn(N)/dZ
         F2Z(IGP))             ! --> Interpolated d2Ln(N)/dZ2
   End If
End Do


!----------------------------------------------------------
! 4. HORIONTAL INTERPOLATION
!----------------------------------------------------------


!--- 4.1. Function and height limits

FP   = Sum(WG(:,0,0)*FV(:))
Zmin = Sum(WG(:,0,0)*ZNmin(:))
Zmax = Sum(WG(:,0,0)*ZNmax(:))


!--- 4.2. Gradient

If (Present(NG)) then
   FG(1) = Sum(WG(:,0,0)*FZ(:))
   FG(2) = Sum(WG(:,0,1)*FV(:))
   FG(3) = Sum(WG(:,1,0)*FV(:))
End If


!--- 4.3. Hessian

If (Present(NH)) then
   FH(1,1) = Sum(WG(:,0,0)*F2Z(:))
   FH(2,2) = Sum(WG(:,0,2)*FV(:))
   FH(3,3) = Sum(WG(:,2,0)*FV(:))
   FH(1,2) = Sum(WG(:,0,1)*FZ(:))
   FH(2,1) = FH(1,2)
   FH(1,3) = Sum(WG(:,1,0)*FZ(:))
   FH(3,1) = FH(1,3)
   FH(3,2) = Sum(WG(:,1,1)*FV(:))
   FH(2,3) = FH(3,2)
End If


!--- 4.4. Imaginary part

If (Allocated(Freq) .and. Present(NIP)) then
   Do IC=1,NC
      FIP(IC) = Sum(WG(:,0,0)*FVI(:,IC))
   End Do
End If


!----------------------------------------------------------
! 5. TRANSFORM FROM LN(N) TO N
!----------------------------------------------------------

NP = Exp(FP)

If (Present(NG)) then
   NG(:) = NP*FG(:)
End If

If (Present(NH)) then
   Do i=1,3
      NH(i,:) = NP*(FH(i,:) + FG(i)*FG(:))
   End Do
End If

If (Present(NIP)) then
   If (Allocated(Freq)) then
      NIP(:) = Exp(FIP(:))
   Else
      NIP(:) = 0.0_Double
   End If
End If


!----------------------------------------------------------
! 6. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(ZNmin, ZNmax, FV, FZ, F2Z)
Deallocate(FVI, FIP, Stat=Stat)


End Subroutine Interpolate_Refractivity



!==========================================================
Subroutine Vertical_Interpolation_LnN &
  (I1,       & ! <-- 1st grid index
   I2,       & ! <-- 2nd grid index
   ID,       & ! <-- Diamond index
   ZP,       & ! <-- Altitude [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   LnNP,     & ! --> Interpolated Ln(N)
   LnNZ,     & ! --> Interpolated dLn(N)/dZ
   LnNZ2,    & ! --> Interpolated d2Ln(N)/dZ2
   LnNIP)      ! ~~> Interpolated Ln(NI)
!
! Vertical interpolation of refractivity profiles from ECHAM.
!----------------------------------------------------------
! Method:
!   Vertical spline interpolation to given altitude.
!----------------------------------------------------------
! (C) Copyright 1998-2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Apr 2000 | Extracted from
!         |             | Interpolate_Refractivity.
!   2.0   | 14 May 2001 | Optional argument LnNIP.
!   3.0   | 07 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!   3.1   | 16 Apr 2002 | Field status.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Init_Spline,  &
    Spline
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
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
   LnNZ          ! Interpolated dLn(N)/dZ
!
Real(Double), Intent(Out) :: &
   LnNZ2         ! Interpolated d2Ln(N)/dZ2
!
Real(Double), Optional, Intent(Out) :: &
   LnNIP(:)      ! Interpolated Ln(NI)[channel]
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NLev    ! Number of levels
Type(Geodetic) :: G       ! Geodetic coordinates of point
Real(Double)   :: PGCL    ! Geocentric latitude of point
Integer        :: NC      ! Number of channels
Integer        :: IC      ! Channel number
Real(Double)   :: LnNIZ   ! Interpolated 1st derivative of LnNI
Real(Double)   :: LnNIZ2  ! Interpolated 2nd derivative of LnNI
Integer        :: ix
!
! Local Arrays:
!
Real(Double) _POINTER :: &
   ZLL(:),       & ! Pointer for cross section of Z
   FLL(:),       & ! Pointer for cross section of LnN or LnNI
   DLL(:)          ! Pointer for cross section of D2N or D2NI
!----------------------------------------------------------

integer :: i

ix = gf% i(i1,i2,id)

if (ix==0) then
  write(0,*)
  write(0,*)   'Vertical_Interpolation_LnN:  invalid index:',i1,i2,id
  write(0,*)   'lbound:', lbound (gf% i)
  write(0,*)   'ubound:', ubound (gf% i)
  write(0,*)
  write(0,'(4x,(200i1))') (mod(i/100,10),i=lbound (gf% i,1),ubound (gf% i,1))
  write(0,'(4x,(200i1))') (mod(i/ 10,10),i=lbound (gf% i,1),ubound (gf% i,1))
  write(0,'(4x,(200i1))') (mod(i    ,10),i=lbound (gf% i,1),ubound (gf% i,1))
  where(gf% i(: ,: ,id) /= 0) gf% i(: ,: ,id) = 8
  where(gf% i(i1,: ,id) == 0) gf% i(i1,: ,id) = 1
  where(gf% i(: ,i2,id) == 0) gf% i(: ,i2,id) = 1
  gf%       i(i1,i2,id)                       = 99
  do i=lbound (gf% i,2), ubound (gf% i,2)
    write(0,'(i3,1x,(200i1))') i, gf% i(:,i,id)
  end do
  call finish ('Vertical_Interpolation_LnN','invalid index')
endif

!----------------------------------------------------------
! 1. CALCULATION OF GRID SIZE
!----------------------------------------------------------

NLev = gf% nz


!----------------------------------------------------------
! 2. VERTICAL INTERPOLATION
!----------------------------------------------------------


!--- 2.1. Vertical profile allocation

If (PStat(I1, I2, ID) == fst_Null) then

   Allocate(  Z(I1, I2, ID)%P(-NM+2:NLev),     &
            LnN(I1, I2, ID)%P(-NM+2:NLev),     &
            D2N(I1, I2, ID)%P(-NM+2:NLev))

   If (Allocated(Freq)) then

      NC = Size(Freq)

      Allocate(LnNI(I1, I2, ID)%M(-NM+2:NLev, NC),     &
               D2NI(I1, I2, ID)%M(-NM+2:NLev, NC))

   End If

   PStat(I1, I2, ID) = fst_Allocated

End If


!--- 2.2. Vertical profile intialization

If (PStat(I1, I2, ID) == fst_Allocated) then

   G    = Geodetic(0.0_Double, gf% s(ix)% XLat, gf% s(ix)% XLon)
   PGCL = gf% s(ix)% GCLat

   If (Allocated(Freq)) then

      NC = Size(Freq)

      Call Make_Refractivity &
        (gf% s   (ix)% Hsur,        & ! <-- Surface geopotential [gpm]
         gf% s   (ix)% Gundu,       & ! <-- Geoid undulation     [m]
         gf% s   (ix)% Psur,        & ! <-- Surface pressure [Pa]
         gf% Ph(:,ix),              & ! <-- Pressure at half levels [Pa]
         gf% T (:,ix),              & ! <-- Temperature [K]
         gf% Q (:,ix),              & ! <-- Specific humidity [kg/kg]
         gf% Qx(:,ix),              & ! <-- Water load [kg/kg]
         G,                         & ! <-- Geodetic coordinates
         PGCL,                      & ! <-- Geocentric latitude [rad]
         LnN   (I1, I2, ID)%P(:),   & ! --> Ln of refractive index
         gf% s   (ix)% Zsur,        & ! --> Surface altitude [km]
         Z     (I1, I2, ID)%P(:),   & ! --> Altitudes of model levels [km]
         LnNI  (I1, I2, ID)%M(:,:))   ! ~~> Ln of dispersive imaginary part of refractive index

      Do IC=1,NC
         Call  Init_Spline   &
           (Z   (I1, I2, ID)%P(:),      & ! <-- Argument grid
            LnNI(I1, I2, ID)%M(:,IC),   & ! <-- Gridded function
            D2NI(I1, I2, ID)%M(:,IC))     ! --> 2nd derivative of spline
      End Do

   Else

      Call Make_Refractivity &
        (gf% s   (ix)% Hsur,        & ! <-- Surface geopotential [gpm]
         gf% s   (ix)% Gundu,       & ! <-- Geoid undulation     [m]
         gf% s   (ix)% Psur,        & ! <-- Surface pressure [Pa]
         gf% Ph(:,ix),              & ! <-- Pressure at half levels [Pa]
         gf% T (:,ix),              & ! <-- Temperature [K]
         gf% Q (:,ix),              & ! <-- Specific humidity [kg/kg]
         gf% Qx(:,ix),              & ! <-- Water load [kg/kg]
         G,                         & ! <-- Geodetic coordinates
         PGCL,                      & ! <-- Geocentric latitude [rad]
         LnN   (I1, I2, ID)%P(:),   & ! --> Ln of refractive index
         gf% s   (ix)% Zsur,        & ! --> Surface altitude [km]
         Z     (I1, I2, ID)%P(:))     ! --> Altitudes of model levels [km]

   End If

   Call  Init_Spline   &
     (Z  (I1, I2, ID)%P(:),      & ! <-- Argument grid
      LnN(I1, I2, ID)%P(:),      & ! <-- Gridded function
      D2N(I1, I2, ID)%P(:))        ! --> 2nd derivative of spline

   PStat(I1, I2, ID) = fst_Initialized

End If


!--- 2.3. Vertical spline interpolation of LnN

ZLL => Z  (I1, I2, ID)%P(:)
FLL => LnN(I1, I2, ID)%P(:)
DLL => D2N(I1, I2, ID)%P(:)

Call Spline  &
  (ZLL,       & ! <-- Argument grid
   FLL,       & ! <-- Gridded function
   DLL,       & ! <-- 2nd derivative of spline
   ZP,        & ! <-- Interpolation point
   LnNP,      & ! --> Interpolated function value
   LnNZ,      & ! --> Interpolated 1st derivative
   LnNZ2)       ! --> Interpolated 2nd derivative


!--- 2.4. Calculation of limits of Z

Zmin = gf% s(ix)% Zsur
Zmax = Max(ZLL(1), ZLL(NLev))


!--- 2.5. Vertical spline interpolation of LnNI

If (Allocated(Freq) .and. Present(LnNIP)) then

   Do IC=1,NC

      FLL => LnNI(I1, I2, ID)%M(:, IC)
      DLL => D2NI(I1, I2, ID)%M(:, IC)

      Call Spline  &
        (ZLL,         & ! <-- Argument grid
         FLL,         & ! <-- Gridded function
         DLL,         & ! <-- 2nd derivative of spline
         ZP,          & ! <-- Interpolation point
         LnNIP(IC),   & ! --> Interpolated function value
         LnNIZ,       & ! --> Interpolated 1st derivative
         LnNIZ2)        ! --> Interpolated 2nd derivative

   End Do

End If



End Subroutine Vertical_Interpolation_LnN



!==========================================================
Subroutine Interpolate_Constituent &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   F,        & ! <-- Constitutent field
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   FP)         ! --> Interpolated constituent
!
! Interpolation of global refractivity field from ECHAM.
!----------------------------------------------------------
! Method:
!   1. Linear interpolation to given altitude.
!   2. Summation of 2D Lon/Lat interpolation series:
!   f(z,Lat,Lon) = sum WG f(z,Lat ,Lon )
!                   i    i       i    i
!   Index i enumerates points of the interpolation
!   subgrid. For rectangular grid, weigting function WG
!   is calculated for polynomial interpolation with the
!   polynomial power defined by the interpolation subgrid
!   dimension NI. For icosahedral grid WG is defined as
!   symplectic weights of the vertices of surrounding
!   triangle.
!----------------------------------------------------------
! (C) Copyright 1999-2002, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 17 Mar 1999 | Original code.
!   7.0   | 07 Apr 2002 | Interpolation of data on
!         |             | icosahedral grid included.
!   7.1   | 16 Apr 2002 | Field status.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Init_Spline,  &
    Linear
!
Use Earth, only: &
! Imported Type Definitions:
    Geodetic
!
Use ECHAM_grid, only: &
! Imported Parameters:
    Longitude_Interpolation,  &
    Latitude_Interpolation
!
Use ICO_grid, only: &
! Imported Array Variables:
    nspoke
!
Use ICO_interpolation, only: &
! Imported Routines:
   Setup_Interpolation
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
Real(WorkPr), Pointer :: &
   F(:,:,:,:)    ! Constitutent field
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
   FP            ! Interpolated N
!----------------------------------------------------------
! Local Scalars:
!
Integer        :: NGP     ! Number of subgrid points for interpolation
Integer        :: IGP     ! Subgrid point index
Integer        :: NLev    ! Number of levels
Integer        :: KLon    ! Longitude subgrid index
Integer        :: KLat    ! Latitude subgrid index
Integer        :: m1, m2  ! Interpolation triangle indices
Integer        :: I1      ! 1st grid index
Integer        :: I2      ! 2nd grid index
Integer        :: ID      ! Diamond index
Type(Geodetic) :: G       ! Geodetic coordinates
Real(Double)   :: PGCL    ! Geocentric latitude of point
integer        :: ix
!
! Local Arrays:
!
! --- Interpolation on rectangular grid
!
Integer :: ILon(NG1)    ! Longitude interpolation subgrid
Integer :: ILat(NG1)    ! Latitude cell indices
Real(Double) :: &
   WLon(NG1),         & ! Weight for longitudinal interpolation
   WLat(NG1)            ! Weight for latitudinal interpolation
!
! --- Interpolation on icosahedral grid
!
Real(Double) :: &
   PRLon(1),          & ! Longitude of interpolation point
   PRLat(1),          & ! Latitude  of interpolation point
   SWS(3,1)             ! Symplectic weights for triangle vertices
Integer :: &
   KIDX(4,1)            ! Interpolation triangle index
!
! --- Interpolation indices and weights
!
Integer, Allocatable :: &
   IDX(:,:)             ! Grid indices of grid points
                        ! [Point, index]
Real(Double), Allocatable :: &
   WG(:)                ! Weighting function
                        ! [Point, Lon.derivative, Lat.derivative]
!
! --- Interpolated values
!
Real(Double), Allocatable :: &
   ZNmin(:),          & ! Minimum model Z on subgrid
   ZNmax(:),          & ! Maximum model Z on subgrid
   FV(:)                ! Vertically interpolated F on subgrid
Real(Double), Pointer :: &
   ZLL(:)               ! Pointer for cross sections of Z
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

!--- 1.1. Number of vertical levels

NLev = gf% nz


!--- 1.2. Interpolation subgrid size

Select Case (gf% Grid_Type)
   Case (WMO6_GAUSSIAN, WMO6_LATLON)
      NGP = NG1**2
   Case (DWD6_ICOSAHEDRON)
      NGP = 3
   case default
      write (0,*) 'Interpolate_Constituent: invalid grid type',gf% Grid_Type
      stop
End Select


!--- 1.3. Memory allocation

Allocate  &
  (IDX(NGP,3), WG(NGP),              &
   ZNmin(NGP), ZNmax(NGP), FV(NGP))


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
         WLon)            ! --> Weights of subgrid points

      !--- 2.1.2. Latitudinal interpolation

      Call Latitude_Interpolation &
        (gf% glat(:),   & ! <-- Latitude grid [deg]
         PLat,          & ! <-- Latitude of point [deg]
         ILat,          & ! --> Interpolation subgrid
         WLat)            ! --> Weights of subgrid points

      !--- 2.1.3. Interpolation indices and weights

      IGP = 0

      Do KLat = 1,NG1
         Do KLon = 1,NG1
            IGP         = IGP + 1
            IDX(IGP,1)  = ILon(KLon)
            IDX(IGP,2)  = ILat(KLat)
            IDX(IGP,3)  = 1
            WG(IGP)     = WLon(KLon)*WLat(KLat)
         End Do
      End Do

   !--- 2.2. Interpolation on icosahedral grid

   Case (DWD6_ICOSAHEDRON)

      PRLon(1) = PLon
      PRLat(1) = PLat

      Call Setup_Interpolation &
        (1,          & ! <-- Number of interpolation points
         PRLon,      & ! <-- Point longitudes [deg]
         PRLat,      & ! <-- Point latitudes [deg]
         KIDX,       & ! --> Index array for interpolation
         SWS)          ! --> Symplectic weights for triangle vertices

      WG(:)     = SWS(:,1)
      m1        = KIDX(4,1)
      m2        = Mod(m1,6) + 1
      IDX(1,:)  = KIDX(1:3,1)
      IDX(2,1)  = KIDX(1,1) + nspoke(1,m1)
      IDX(2,2)  = KIDX(2,1) + nspoke(2,m1)
      IDX(2,3)  = KIDX(3,1)
      IDX(3,1)  = KIDX(1,1) + nspoke(1,m2)
      IDX(3,2)  = KIDX(2,1) + nspoke(2,m2)
      IDX(3,3)  = KIDX(3,1)

   case default
      write (0,*) 'Interpolate_Constituent: invalid grid type',gf% Grid_Type
      stop
End Select


!----------------------------------------------------------
! 3. VERTICAL INTERPOLATION
!----------------------------------------------------------

Do IGP = 1,NGP

   I1 = IDX(IGP, 1)
   I2 = IDX(IGP, 2)
   ID = IDX(IGP, 3)
   ix = gf% i(i1,i2,id)

   If (PStat(I1, I2, ID) == fst_Null) then
      Allocate(  Z(I1, I2, ID)%P(-NM+2:NLev),  &
               LnN(I1, I2, ID)%P(-NM+2:NLev),  &
               D2N(I1, I2, ID)%P(-NM+2:NLev))
      PStat(I1, I2, ID) = fst_Allocated
   End If
   If (PStat(I1, I2, ID) == fst_Allocated) then
      G    = Geodetic(0.0_Double, gf% s(ix)% XLat, gf% s(ix)% XLon)
      PGCL = gf% s(ix)% GCLat
      Call Make_Refractivity      &
        (gf% s   (ix)% Hsur,      &
         gf% s   (ix)% Gundu,     &
         gf% s   (ix)% Psur,      &
         gf% Ph(:,ix),            & ! <-- Pressure at half levels [Pa]
         gf% T (:,ix),            &
         gf% Q (:,ix),            &
         gf% Qx(:,ix),            & ! <-- Water load [kg/kg]
         G,                       &
         PGCL,                    &
         LnN  (I1, I2, ID)%P(:),  &
         gf% s   (ix)% Zsur,      &
         Z    (I1, I2, ID)%P(:))
      Call Init_Spline   &
        (Z  (I1, I2, ID)%P(:),    &
         LnN(I1, I2, ID)%P(:),    &
         D2N(I1, I2, ID)%P(:))
      PStat(I1, I2, ID) = fst_Initialized
   End If
   ZLL => Z(I1, I2, ID)%P(1:NLev)
   Call Linear(ZLL, Real(F(I1, I2, ID, :), Double), &
               ZP, FV(IGP), CExt =.True.)
   ZNmin(IGP) = gf% s(ix)% Zsur
   ZNmax(IGP) = Max(ZLL(1), ZLL(NLev))
End Do


!----------------------------------------------------------
! 4. HORIONTAL INTERPOLATION
!----------------------------------------------------------

FP   = Sum(WG(:)*FV(:))
Zmin = Sum(WG(:)*ZNmin(:))
Zmax = Sum(WG(:)*ZNmax(:))


!----------------------------------------------------------
! 6. MEMORY DEALLOCATION
!----------------------------------------------------------

Deallocate(IDX, WG, ZNmin, ZNmax, FV)


End Subroutine Interpolate_Constituent



!==========================================================
Subroutine ECHAM_NGradN &
  (X,        & ! <-- Cartesian coordinates of point
   NGradN,   & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   Stat,     & ! --> Error status
   NIP)        ! ~~> Interpolated Im(N)
!
! Calculation of interpolated (1 + N)*Grad(N) for
! ECHAM global fields.
!----------------------------------------------------------
! Method:
!   Calculation of gradient in geodetic coordinates and
!   transform to Cartesian coordinates.
!----------------------------------------------------------
! (C) Copyright 1999-2001, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 25 Jan 1999 | Original code.
!   2.0   | 18 Feb 1999 | Argument NP.
!   2.1   | 20 Feb 1999 | Relaxation to MSIS profile.
!   2.2   | 10 Mar 1999 | Stat=1 when G%H < Zmin.
!   3.0   | 22 Jun 1999 | No combination with N_MSIS.
!   4.0   | 15 May 2001 | Argument NIP.
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
    Jacobian_GC
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Type(Cartesian), Intent(In) :: &
   X        ! Cartesian coordinates of point
!
! Output arguments:
!
Type(Cartesian), Intent(Out) :: &
   NGradN   ! Interpolated (1 + N)*Grad(N)
!
Real(Double), Intent(Out) :: &
   NP       ! Interpolated N
!
Integer, Intent(Out)         :: &
   Stat     ! Error status:
            !   0 - point above surface
            !   1 - point under surface
!
Real(Double), Optional, Intent(Out) :: &
   NIP(:)   ! Interpolated Im(N)
!----------------------------------------------------------
! Local Scalars:
!
Type(Geodetic) :: G     ! Geodetic coordinates of point
Real(Double)   :: PLon  ! Longitude of point
Real(Double)   :: PLat  ! Latitude of point
Real(Double)   :: ZP    ! Altitude of point
Real(Double)   :: Zmin  ! Minimum model Z
Real(Double)   :: Zmax  ! Maximum model Z
!
! Local Arrays:
!
Real(Double)   :: &
   NGP(3)     ! Interpolated grad(N) in geodetic coordinates
Real(Double)   :: &
   JGC(3,3)   ! Jacoubian d(Geodetic)/d(Cartesian)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INTERPOLATION OF ECHAM-MSIS FIELD
!----------------------------------------------------------

G    = Geod_from_Cart(X)
PLon = G%Lambda
PLat = G%Phi
ZP   = G%H

Call Interpolate_Refractivity &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NGP,      & ! ~~> Interpolated dN/d(alt,lat,lon)
   NIP = NIP)  ! ~~> Interpolated Im(N)


!----------------------------------------------------------
! 2. TRANSFORM TO CARTESIAN COORDINATES
!----------------------------------------------------------

JGC    = Jacobian_GC(G)

NGradN%X(:) = (1+NP)*MatMul(NGP(:), JGC(:,:))


!----------------------------------------------------------
! 3. STATUS DEFINITION
!----------------------------------------------------------

If (G%H >= Zmin) then
   Stat = 0
Else
   Stat = 1
End If


End Subroutine ECHAM_NGradN



!==========================================================
Subroutine ECHAM_NGHN &
  (X,        & ! <-- Cartesian coordinates of point
   NGN,      & ! --> Interpolated (1 + N)*Grad(N)
   NP,       & ! --> Interpolated N
   NHN,      & ! --> Interpolated Grad x (1+N)Grad(N)
   Stat)       ! --> Error status
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
!   1.0   | 03 May 2000 | Version with hessian.
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
!----------------------------------------------------------
! Local Scalars:
!
Type(Geodetic) :: G     ! Geodetic coordinates of point
Real(Double)   :: PLon  ! Longitude of point
Real(Double)   :: PLat  ! Latitude of point
Real(Double)   :: ZP    ! Altitude of point
Real(Double)   :: Zmin  ! Minimum model Z
Real(Double)   :: Zmax  ! Maximum model Z
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
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INTERPOLATION OF ECHAM-MSIS FIELD
!----------------------------------------------------------


!--- 1.1. Coordinate calculation

G    = Geod_from_Cart(X)
PLon = G%Lambda
PLat = G%Phi
ZP   = G%H


!--- 1.2. Interpolation and calculation
!---      of T, P, and Q derivatives

Call Interpolate_Refractivity &
  (PLon,     & ! <-- Longitude of point [deg]
   PLat,     & ! <-- Latitude  of point [deg]
   ZP,       & ! <-- Altitude  of point [km]
   Zmin,     & ! --> Minimum model Z for this lon/lat
   Zmax,     & ! --> Maximum model Z for this lon/lat
   NP,       & ! --> Interpolated N
   NG,       & ! --> Interpolated dN/d(alt,lat,lon)
   NH)         ! --> Interpolated hessian matrix of N


!----------------------------------------------------------
! 2. TRANSFORM TO CARTESIAN COORDINATES
!----------------------------------------------------------


!--- 2.1. Transform of gradient

JGC    = Jacobian_GC(G)

NGN%X(:) = (1+NP)*MatMul(NG(:), JGC(:,:))


!--- 2.2. Transform of hessian

HGC = Hessian_GC(G)

Do i=1,3
   Do j=1,3
      NHN(i,j) = &
         Sum(JGC(:,i)*NG(:))*Sum(JGC(:,j)*NG(:)) + &
         (1+NP)*Sum(HGC(:,i,j)*NG(:))            + &
         (1+NP)*Sum(JGC(:,i)*MatMul(NH(:,:),JGC(:,j)))
   End Do
End Do


!----------------------------------------------------------
! 3. STATUS DEFINITION
!----------------------------------------------------------

If (G%H >= Zmin) then
   Stat = 0
Else
   Stat = 1
End If


End Subroutine ECHAM_NGHN



End Module ECHAM_fields
