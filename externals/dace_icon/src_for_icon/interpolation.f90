!
!+ GNSS Radio occultation operator: spline, linear, polynomial interpolation
!
MODULE Interpolation
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Module for spline, linear, and polynomial interpolation.
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
!  Seek_Index_Double: replace original algorithm by faster bisection
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  remove unused variables
! V1_14        2011/11/08 Harald Anlauf
!  Cleanup and minor optimization
! V1_19        2012-04-16 Harald Anlauf
!  interpolation, interpolation_adj: optimize for sxf90 rev.430+
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
! Module Interpolation
!
! Module for spline, linear, and polynomial interpolation.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 13 Oct 1997 | Original code.
!   1.1   | 14 Sep 1998 | Modified Linear.
!   2.0   | 29 Sep 1998 | Subroutine Polynomial.
!   3.0   | 23 Dec 1998 | Single and Double versions
!         |             | of spline interpolation.
!   4.0   | 15 Jan 1999 | Seek_Index.
!   5.0   | 18 Jan 1999 | Single and Double versions
!         |             | of linear interpolation.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Single, Double
!----------------------------------------------------------
Implicit None
Private
Public :: Init_Spline, Spline, Seek_Index, Polynomial, Linear
!----------------------------------------------------------
! Interfaces:
!
Interface Init_Spline
   Module Procedure Init_Spline_Double
   Module Procedure Init_Spline_Single
End Interface
!
Interface Spline
   Module Procedure Spline_Double
   Module Procedure Spline_Single
End Interface
!
Interface Linear
   Module Procedure Linear_Double
   Module Procedure Linear_Single
End Interface
!
Interface Seek_Index
   Module Procedure Seek_Index_Double
   Module Procedure Seek_Index_Single
End Interface
!----------------------------------------------------------
Contains

!==========================================================
Subroutine Init_Spline_Double &
  (x,      & ! <-- Argument grid
   f,      & ! <-- Gridded function
   d2,     & ! --> 2nd derivative of spline
   ErrCode)  ! ~~> Error code
!
! Calculation of 2nd derivative of spline.
!----------------------------------------------------------
! Method:
!   Drive-through solution.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   1.1   | 13 Oct 1997 | Optional argument ErrCode.
!   1.2   | 14 Sep 1998 | Improved interface.
!   2.0   | 06 Oct 1998 | Any monotonous argument possible.
!   3.0   | 23 Dec 1998 | Single and Double versions.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)  :: &
   x(1:)   ! Argument grid
           ! x(i) must be monotonous.
Real(Double), Intent(In)  :: &
   f(1:)   ! Gridded function
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   d2(1:)  ! 2nd derivative of spline
Integer, Intent(Out), Optional          :: &
   ErrCode ! Error code:
           !   0 - no error
           !   1 - x is not mononous
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   dfl, dfr, df, a, b, c
Integer      :: &
   i,        & ! Array index
   N           ! Number of array elements
!
! Local Arrays:
!
Real(Double) :: &
   d1(Size(x)) ! Work array for drive-through
!----------------------------------------------------------


!----------------------------------------------------------
! 0. CHECKING PARAMETERS
!----------------------------------------------------------

N = Size(x)

If (Present (ErrCode)) then
!CDIR NEIGHBORS
   If (All(x(1:N-1) < x(2:N)) .or. &
       All(x(1:N-1) > x(2:N))) then
      ErrCode = 0
   Else
      ErrCode = 1
      Return
   End If
End If


!----------------------------------------------------------
! 1. DRIVE-THROUGH CALCULATION OF SPLINE COEFFICIENTS
!----------------------------------------------------------

d2(1) = 0.0_Double
d2(N) = 0.0_Double
d1(1) = 0.0_Double

Do i = 2, N-1
   dfl   = (f(i) - f(i-1))/(x(i) - x(i-1))
   dfr   = (f(i+1) - f(i))/(x(i+1) - x(i))
   df    = (dfr - dfl)/(x(i+1) - x(i-1))
   a     = (x(i) - x(i-1))/(2*(x(i+1) - x(i-1)))
   b     = (x(i+1) - x(i))/(2*(x(i+1) - x(i-1)))
   c     = 1 + a*d1(i-1)
   d1(i) = -b/c
   d2(i) = (3*df - a*d2(i-1))/c
End Do

Do i = N-1, 2, -1
   d2(i) = d1(i)*d2(i+1) + d2(i)
End Do


End Subroutine Init_Spline_Double



!==========================================================
Subroutine Init_Spline_Single &
  (x,      & ! <-- Argument grid
   f,      & ! <-- Gridded function
   d2,     & ! --> 2nd derivative of spline
   ErrCode)  ! ~~> Error code
!
! Calculation of 2nd derivative of spline.
!----------------------------------------------------------
! Method:
!   Drive-through solution.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   1.1   | 13 Oct 1997 | Optional argument ErrCode.
!   1.2   | 14 Sep 1998 | Improved interface.
!   2.0   | 06 Oct 1998 | Any monotonous argument possible.
!   3.0   | 23 Dec 1998 | Single and Double versions.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Single), Intent(In)  :: &
   x(1:)   ! Argument grid
           ! x(i) must be monotonous.
Real(Single), Intent(In)  :: &
   f(1:)   ! Gridded function
!
! Output arguments:
!
Real(Single), Intent(Out) :: &
   d2(1:)  ! 2nd derivative of spline
Integer, Intent(Out), Optional          :: &
   ErrCode ! Error code:
           !   0 - no error
           !   1 - x is not mononous
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   dfl, dfr, df, a, b, c
Integer      :: &
   i,        & ! Array index
   N           ! Number of array elements
!
! Local Arrays:
!
Real(Double) :: &
   d1(Size(x)) ! Work array for drive-through
!----------------------------------------------------------


!----------------------------------------------------------
! 0. CHECKING PARAMETERS
!----------------------------------------------------------

N = Size(x)

If (Present (ErrCode)) then
!CDIR NEIGHBORS
   If (All(x(1:N-1) < x(2:N)) .or. &
       All(x(1:N-1) > x(2:N))) then
      ErrCode = 0
   Else
      ErrCode = 1
      Return
   End If
End If


!----------------------------------------------------------
! 1. DRIVE-THROUGH CALCULATION OF SPLINE COEFFICIENTS
!----------------------------------------------------------

d2(1) = 0.0_Single
d2(N) = 0.0_Single
d1(1) = 0.0_Double

Do i = 2, N-1
   dfl   = (f(i) - f(i-1))/(x(i) - x(i-1))
   dfr   = (f(i+1) - f(i))/(x(i+1) - x(i))
   df    = (dfr - dfl)/(x(i+1) - x(i-1))
   a     = (x(i) - x(i-1))/(2*(x(i+1) - x(i-1)))
   b     = (x(i+1) - x(i))/(2*(x(i+1) - x(i-1)))
   c     = 1 + a*d1(i-1)
   d1(i) = -b/c
   d2(i) = (3*df - a*d2(i-1))/c
End Do

Do i = N-1, 2, -1
   d2(i) = d1(i)*d2(i+1) + d2(i)
End Do


End Subroutine Init_Spline_Single



!==========================================================
Subroutine Spline_Double &
  (x,      & ! <-- Argument grid
   f,      & ! <-- Gridded function
   d2,     & ! <-- 2nd derivative of spline
   x_int,  & ! <-- Interpolation point
   f_int,  & ! --> Interpolated function value
   fd_int, & ! ~~> Interpolated 1st derivative
   fd2_int)  ! ~~> Interpolated 2nd derivative
!
! Spline interpolation of a gridded function with linear
! extrapolation outside grid extent.
!----------------------------------------------------------
! Method:
!   Searching for grid interval containing the
!   interpolation point and summation of polynomial
!   with given spline-coefficients.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   1.1   | 14 Sep 1998 | Improved interface.
!   2.0   | 06 Oct 1998 | Any monotonous argument possible.
!   3.0   | 23 Dec 1998 | Single and Double versions.
!   4.0   | 15 Jan 1999 | Use of Seek_Index.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   x(1:)    ! Argument grid
            ! x(i) must be monotonous.
Real(Double), Intent(In) :: &
   f(1:)    ! Gridded function
Real(Double), Intent(In) :: &
   d2(1:)   ! 2nd derivative of spline
Real(Double), Intent(In)               :: &
   x_int    ! Interpolation point
!
! Output arguments:
!
Real(Double), Intent(Out)              :: &
   f_int   ! Interpolated function value
Real(Double), Intent(Out), Optional    :: &
   fd_int  ! Interpolated 1st derivative
Real(Double), Intent(Out), Optional    :: &
   fd2_int ! Interpolated 2nd derivative
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   a1, a2, a3,  & ! Polynomial coefficients
   dx,          & ! Grid interval
   dx_t,        & ! Grid-point to interpolation-point distance
   x_t,         & ! Interpolation point projected to grid extent
   fd             ! Interpolated derivative
Integer  :: &
   i,           & ! Array index
   N              ! Number of data
!Integer :: Dir   ! Direction of argument change
Integer  :: i_int ! Interpolation interval index
!----------------------------------------------------------


!----------------------------------------------------------
! 1. LOCATION OF INTERPOLATION POINT INSIDE GRID
!----------------------------------------------------------

N     = Size(x)

x_t   = Min(Max(x_int, Min(x(1),x(N))), Max(x(1),x(N)))
!Dir   = Nint(Sign(1.0_Double, x(N)-x(1)))

! i_int = Sum(MaxLoc(Dir*x(1:N-1), Dir*x(1:N-1) <= Dir*x_t))

i_int = Seek_Index(x, x_t)
i     = Max(i_int, 1)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION COEFFICIENTS
!----------------------------------------------------------

dx    = x(i+1) - x(i)
a2    = d2(i)/2
a3    = (d2(i+1) - d2(i))/(6*dx)
a1    = (f(i+1) - f(i))/dx - dx*(a2 + dx*a3)


!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATED VALUE
!----------------------------------------------------------

dx_t  = x_t - x(i)
fd    = a1 + dx_t*(2*a2 + dx_t*3*a3)
f_int = f(i) + dx_t*(a1 + dx_t*(a2 + dx_t*a3)) + &
        fd*(x_int - x_t)

If (Present(fd_int)) then
   fd_int = fd
End if

If (Present(fd2_int)) then
   fd2_int  = 2*a2 + dx_t*6*a3
End if


End Subroutine Spline_Double



!==========================================================
Subroutine Spline_Single &
  (x,      & ! <-- Argument grid
   f,      & ! <-- Gridded function
   d2,     & ! <-- 2nd derivative of spline
   x_int,  & ! <-- Interpolation point
   f_int,  & ! --> Interpolated function value
   fd_int, & ! ~~> Interpolated 1st derivative
   fd2_int)  ! ~~> Interpolated 2nd derivative
!
! Spline interpolation of a gridded function with linear
! extrapolation outside grid extent.
!----------------------------------------------------------
! Method:
!   Searching for grid interval containing the
!   interpolation point and summation of polynomial
!   with given spline-coefficients.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   1.1   | 14 Sep 1998 | Improved interface.
!   2.0   | 06 Oct 1998 | Any monotonous argument possible.
!   3.0   | 23 Dec 1998 | Single and Double versions.
!   4.0   | 15 Jan 1999 | Use of Seek_Index.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Single), Intent(In) :: &
   x(1:)    ! Argument grid
            ! x(i) must be monotonous.
Real(Single), Intent(In) :: &
   f(1:)    ! Gridded function
Real(Single), Intent(In) :: &
   d2(1:)   ! 2nd derivative of spline
Real(Single), Intent(In)               :: &
   x_int    ! Interpolation point
!
! Output arguments:
!
Real(Single), Intent(Out)              :: &
   f_int   ! Interpolated function value
Real(Single), Intent(Out), Optional    :: &
   fd_int  ! Interpolated 1st derivative
Real(Single), Intent(Out), Optional    :: &
   fd2_int ! Interpolated 2nd derivative
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   a1, a2, a3,  & ! Polynomial coefficients
   dx,          & ! Grid interval
   dx_t,        & ! Grid-point to interpolation-point distance
   x_t,         & ! Interpolation point projected to grid extent
   fd             ! Interpolated derivative
Integer  :: &
   i,           & ! Array index
   N              ! Number of data
!Integer :: Dir   ! Direction of argument change
Integer  :: i_int ! Interpolation interval index
!----------------------------------------------------------


!----------------------------------------------------------
! 1. LOCATION OF INTERPOLATION POINT INSIDE GRID
!----------------------------------------------------------

N     = Size(x)

x_t   = Min(Max(x_int, Min(x(1),x(N))), Max(x(1),x(N)))
!Dir   = Nint(Sign(1.0_Single, x(N)-x(1)))

! i_int = Sum(MaxLoc(Dir*x(1:N-1), Dir*x(1:N-1) <= Dir*x_t))

i_int = Seek_Index(x, Real(x_t, Single))
i     = Max(i_int, 1)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION COEFFICIENTS
!----------------------------------------------------------

dx    = x(i+1) - x(i)
a2    = d2(i)/2
a3    = (d2(i+1) - d2(i))/(6*dx)
a1    = (f(i+1) - f(i))/dx - dx*(a2 + dx*a3)


!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATED VALUE
!----------------------------------------------------------

dx_t  = x_t - x(i)
fd    = a1 + dx_t*(2*a2 + dx_t*3*a3)
f_int = f(i) + dx_t*(a1 + dx_t*(a2 + dx_t*a3)) + &
        fd*(x_int - x_t)

If (Present(fd_int)) then
   fd_int = fd
End if

If (Present(fd2_int)) then
   fd2_int  = 2*a2 + dx_t*6*a3
End if


End Subroutine Spline_Single



!==========================================================
Subroutine Linear_Double &
  (x,       & ! <-- Argument grid
   f,       & ! <-- Gridded function
   x_int,   & ! <-- Interpolation point
   f_int,   & ! --> Interpolated function value
   ErrCode, & ! ~~> Error code
   CExt)      ! <~~ Constant/linear extrapolation
!
! Linear interpolation of a gridded function
!----------------------------------------------------------
! Method:
!   Searching for grid interval containing the
!   interpolation point and linear interpolation
!   within interval.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   2.0   | 14 Sep 1998 | Subroutine, parameter ErrCode.
!   3.0   | 06 Oct 1998 | Any monotonous argument possible.
!   4.0   | 04 Nov 1998 | Constant extrapolation possible.
!   5.0   | 15 Jan 1999 | Use of Seek_Index.
!   6.0   | 18 Jan 1999 | Single and Double versions.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   x(:)  ! Argument grid
         ! x(i) must increase monotonously.
Real(Double), Intent(In) :: &
   f(:)  ! Gridded function f(x)
Real(Double), Intent(In)               :: &
   x_int ! Interpolation point
!
! Output arguments:
!
Real(Double), Intent(Out)              :: &
   f_int   ! Interpolated function value
Integer, Intent(Out), Optional         :: &
   ErrCode ! Error code:
           !   0 - no error
           !   1 - x is not mononously increasing
!
! Input arguments:
!
Logical, Intent(In), Optional          :: &
   CExt    ! Constant/linear extrapolation:
           !   .True.  - constant extrapolation
           !   .False. - linear extrapolation (default)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: a1     ! Slope
Real(Double) :: x_t    ! Interpolation point projected to grid extent
Integer      :: i      ! Array index
Integer      :: N      ! Upper array index
!Integer     :: Dir    ! Direction of argument change
Logical      :: CX     ! Constant/linear extrapolation
Integer      :: i_int  ! Interpolation interval index
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION, CHECKING PARAMETERS
!----------------------------------------------------------

N = Size(x)

If (Present (ErrCode)) then
!CDIR NEIGHBORS
   If (All(x(1:N-1) < x(2:N)) .or. &
       All(x(1:N-1) > x(2:N))) then
      ErrCode = 0
   Else
      ErrCode = 1
      Return
   End If
End If

If (Present(CExt)) then
   CX = CExt
Else
   CX = .False.
End If


!----------------------------------------------------------
! 1. LOCATION OF INTERPOLATION POINT INSIDE GRID
!----------------------------------------------------------

x_t   = Min(Max(x_int, Min(x(1),x(N))), Max(x(1),x(N)))
!Dir   = Nint(Sign(1.0_Double, x(N)-x(1)))

! i_int = Sum(MaxLoc(Dir*x(1:N-1), Dir*x(1:N-1) <= Dir*x_t))

i_int = Seek_Index(x, x_t)
i     = Max(i_int, 1)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATED VALUE
!----------------------------------------------------------

a1 = (f(i+1) - f(i))/(x(i+1) - x(i))

If (CX) then
   f_int = f(i) + (x_t - x(i))*a1
Else
   f_int = f(i) + (x_int - x(i))*a1
End If


End Subroutine Linear_Double



!==========================================================
Subroutine Linear_Single &
  (x,       & ! <-- Argument grid
   f,       & ! <-- Gridded function
   x_int,   & ! <-- Interpolation point
   f_int,   & ! --> Interpolated function value
   ErrCode, & ! ~~> Error code
   CExt)      ! <~~ Constant/linear extrapolation
!
! Linear interpolation of a gridded function
!----------------------------------------------------------
! Method:
!   Searching for grid interval containing the
!   interpolation point and linear interpolation
!   within interval.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 12 Oct 1997 | Original code.
!   2.0   | 14 Sep 1998 | Subroutine, parameter ErrCode.
!   3.0   | 06 Oct 1998 | Any monotonous argument possible.
!   4.0   | 04 Nov 1998 | Constant extrapolation possible.
!   5.0   | 15 Jan 1999 | Use of Seek_Index.
!   6.0   | 18 Jan 1999 | Single and Double versions.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Single), Intent(In) :: &
   x(:)  ! Argument grid
         ! x(i) must increase monotonously.
Real(Single), Intent(In) :: &
   f(:)  ! Gridded function f(x)
Real(Single), Intent(In)               :: &
   x_int ! Interpolation point
!
! Output arguments:
!
Real(Single), Intent(Out)              :: &
   f_int   ! Interpolated function value
Integer, Intent(Out), Optional         :: &
   ErrCode ! Error code:
           !   0 - no error
           !   1 - x is not mononously increasing
!
! Input arguments:
!
Logical, Intent(In), Optional          :: &
   CExt    ! Constant/linear extrapolation:
           !   .True.  - constant extrapolation
           !   .False. - linear extrapolation (default)
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: a1     ! Slope
Real(Double) :: x_t    ! Interpolation point projected to grid extent
Integer      :: i      ! Array index
Integer      :: N      ! Upper array index
!Integer     :: Dir    ! Direction of argument change
Logical      :: CX     ! Constant/linear extrapolation
Integer      :: i_int  ! Interpolation interval index
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION, CHECKING PARAMETERS
!----------------------------------------------------------

N = Size(x)

If (Present (ErrCode)) then
!CDIR NEIGHBORS
   If (All(x(1:N-1) < x(2:N)) .or. &
       All(x(1:N-1) > x(2:N))) then
      ErrCode = 0
   Else
      ErrCode = 1
      Return
   End If
End If

If (Present(CExt)) then
   CX = CExt
Else
   CX = .False.
End If


!----------------------------------------------------------
! 1. LOCATION OF INTERPOLATION POINT INSIDE GRID
!----------------------------------------------------------

x_t   = Min(Max(x_int, Min(x(1),x(N))), Max(x(1),x(N)))
!Dir   = Nint(Sign(1.0_Single, x(N)-x(1)))

! i_int = Sum(MaxLoc(Dir*x(1:N-1), Dir*x(1:N-1) <= Dir*x_t))

i_int = Seek_Index(x, Real(x_t, Single))
i     = Max(i_int, 1)


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATED VALUE
!----------------------------------------------------------

a1 = (f(i+1) - f(i))/(x(i+1) - x(i))

If (CX) then
   f_int = f(i) + (x_t - x(i))*a1
Else
   f_int = f(i) + (x_int - x(i))*a1
End If


End Subroutine Linear_Single



!==========================================================
Function Seek_Index_Double &
  (x,         & ! <-- X-grid
   xp)        & ! <-- Point to locate
Result(ip)      ! --> Index of grid interval containing xp
!
! Locating grid interval containing given point.
!----------------------------------------------------------
! Method:
!   Combination of Newton-like and dichotomic
!   interative search.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Jan 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   x(1:)    ! X-grid
            ! Grid must be homogeneous, this is not checked.
!
Real(Double), Intent(In) :: &
   xp       ! Point to locate inside grid
!
! Function result:
!
Integer :: &
   ip       ! Index of grid interval containing xp
            ! >0 - ip: (x(ip) <= xp <= x(ip+1))
            !  0 - xp is not inside grid
            ! -1 - iterations do not converge
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N      ! Number of grid points
Integer      :: Imin   ! Upper estimate of index
Integer      :: Imax   ! Lower estimate of index
Real(Double) :: Dir    ! Direction of argument change
!Integer      :: It     ! Iteration count
!Integer      :: di     ! Index increment in iterations
!Integer      :: is     ! Step direction count
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------


!--- 1.1. Grid size and direction calculation

N    = Size(x)
Dir  = Sign(1.0_Double, x(N)-x(1))


!--- 1.2. Checking if point is inside grid

If ((Dir*xp < Dir*x(1)) .or. (Dir*xp > Dir*x(N))) then
   ip = 0
   Return
End If


#if 1

!------------------------------
! Find interval using bisection
!------------------------------
Imin = 1
Imax = N
if (Dir >= 0._Double) then
   do while (Imax - Imin > 1)
      ip = (Imax + Imin) / 2
      if (x(ip) > xp) then
         Imax = ip
      else
         Imin = ip
      endif
   end do
else
   do while (Imax - Imin > 1)
      ip = (Imax + Imin) / 2
      if (x(ip) < xp) then
         Imax = ip
      else
         Imin = ip
      endif
   end do
end if
ip = Imin

#else

!----------------------------------------------------------
! 2. INDEX SEARCH
!----------------------------------------------------------


!--- 2.1. Initial approximation

Imin = 1
Imax = N
ip   = Imin + Floor(Real(Imax - Imin, Double)*(xp - x(Imin))/(x(Imax) - x(Imin)))
ip   = Max(1, Min(ip, N-1))


!--- 2.2. Iterative index search

It = 0
is = 0

Search: Do
   If ((Dir*x(ip) <= Dir*xp) .and. (Dir*xp <= Dir*x(ip+1))) then
      Exit Search
   End If
   If (Abs(is) > 1) then
      ip = (Imax + Imin)/2
      is = 0
   End If
   If (Dir*x(ip+1) < Dir*xp) then
      Imin = ip + 1
      di   = Floor(Real(Imax - Imin, Double)*(xp - x(Imin))/(x(Imax) - x(Imin)))
      ip   = Imin + di
      is   = is + 1
   Else If (Dir*xp < Dir*x(ip)) then
      Imax = ip
      di   = Floor(Real(Imin - Imax, Double)*(xp - x(Imax))/(x(Imin) - x(Imax)))
      ip   = Imax + di
      is   = is - 1
   End If
   ip = Max(1, Min(ip, N-1))
   It = It + 1
   If (It > N) then
      ip = -1
      Exit Search
   End If
End Do Search

#endif


End Function Seek_Index_Double



!==========================================================
Function Seek_Index_Single &
  (x,         & ! <-- X-grid
   xp)        & ! <-- Point to locate
Result(ip)      ! --> Index of grid interval containing xp
!
! Locating grid interval containing given point.
!----------------------------------------------------------
! Method:
!   Combination of Newton-like and dichotomic
!   interative search.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 15 Jan 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Single), Intent(In) :: &
   x(1:)    ! X-grid
            ! Grid must be homogeneous, this is not checked.
!
Real(Single), Intent(In) :: &
   xp       ! Point to locate inside grid
!
! Function result:
!
Integer :: &
   ip       ! Index of grid interval containing xp
            ! >0 - ip: (x(ip) <= xp <= x(ip+1))
            !  0 - xp is not inside grid
            ! -1 - iterations do not converge
!----------------------------------------------------------
! Local Scalars:
!
Integer  :: N      ! Number of grid points
Integer  :: Imin   ! Upper estimate of index
Integer  :: Imax   ! Lower estimate of index
Integer  :: Dir    ! Direction of argument change
Integer  :: It     ! Iteration count
Integer  :: di     ! Index increment in iterations
Integer  :: is     ! Step direction count
!----------------------------------------------------------


!----------------------------------------------------------
! 1. INITIALIZATION
!----------------------------------------------------------


!--- 1.1. Grid size and direction calculation

N    = Size(x)
Dir  = Nint(Sign(1.0_Single, x(N)-x(1)))


!--- 1.2. Checking if point is inside grid

If ((Dir*xp < Dir*x(1)) .or. (Dir*xp > Dir*x(N))) then
   ip = 0
   Return
End If


!----------------------------------------------------------
! 2. INDEX SEARCH
!----------------------------------------------------------


!--- 2.1. Initial approximation

Imin = 1
Imax = N
ip   = Imin + Floor(Real(Imax - Imin, Single)*(xp - x(Imin))/(x(Imax) - x(Imin)))
ip   = Max(1, Min(ip, N-1))


!--- 2.2. Iterative index search

It = 0
is = 0

Search: Do
   If ((Dir*x(ip) <= Dir*xp) .and. (Dir*xp <= Dir*x(ip+1))) then
      Exit Search
   End If
   If (Abs(is) > 1) then
      ip = (Imax + Imin)/2
      is = 0
   End If
   If (Dir*x(ip+1) < Dir*xp) then
      Imin = ip + 1
      di   = Floor(Real(Imax - Imin, Single)*(xp - x(Imin))/(x(Imax) - x(Imin)))
      ip   = Imin + di
      is   = is + 1
   Else If (Dir*xp < Dir*x(ip)) then
      Imax = ip
      di   = Floor(Real(Imin - Imax, Single)*(xp - x(Imax))/(x(Imin) - x(Imax)))
      ip   = Imax + di
      is   = is - 1
   End If
   ip = Max(1, Min(ip, N-1))
   It = It + 1
   If (It > N) then
      ip = -1
      Exit Search
   End If
End Do Search


End Function Seek_Index_Single



!==========================================================
Subroutine Polynomial &
  (c,     & ! <-- Polynomial coefficients
   x,     & ! <-- Polynomial argument
   P,     & ! --> Polynomial value
   DP)      ! ~~> Polynomial derivative
!
! Calculation of polynomial and its derivative
!----------------------------------------------------------
! Method:
!   Horner scheme.
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
Real(Double), Intent(In) :: &
   c(0:)    ! Polynomial coefficients
!
Real(Double), Intent(In) :: &
   x        ! Polynomial argument
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   P        ! Polynomial value
            ! P = sum(c(i)*x**i, i=0..UBound(c,1))
!
Real(Double), Intent(Out), Optional :: &
   DP       ! Polynomial derivative
            ! DP = sum(i*c(i)*x**(i-1), i=0..UBound(c,1))
!----------------------------------------------------------
! Local Scalars:
!
Integer      :: N  ! Polynomial power
Integer      :: i  ! Power index
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF POLYNOMIAL VALUE
!----------------------------------------------------------

N = UBound(c,1)

P = c(N)
Do i = N-1, 0, -1
   P = c(i) + x*P
End Do


!----------------------------------------------------------
! 2. CALCULATION OF POLYNOMIAL DERIVATIVE
!----------------------------------------------------------

If (.not. Present(DP)) then
   Return
End If

DP = c(N)*N
Do i = N-1, 1, -1
   DP = c(i)*i + x*DP
End Do


End Subroutine Polynomial



!==========================================================
Subroutine Extremum &
  (X,     & ! <-- Arguments
   Y,     & ! <-- Function values
   Xext,  & ! --> Extremum location
   Yext)    ! --> Extremum value
!
! Calculation of extremum of a square polynomial
!----------------------------------------------------------
! Method:
!   Standard formulas.
!----------------------------------------------------------
! (C) Copyright 1999, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 26 Nov 1999 | Original code.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In) :: &
   X(0:2)   ! Arguments
!
Real(Double), Intent(In) :: &
   Y(0:2)   ! Function values
!
! Output arguments:
!
Real(Double), Intent(Out) :: &
   Xext     ! Extremum location
!
Real(Double), Intent(Out) :: &
   Yext     ! Extremum value
!----------------------------------------------------------
! Local Arrays:
!
Real(Double) :: A(0:2)  ! Interpolation coefficients.
!----------------------------------------------------------


!----------------------------------------------------------
! 1. CALCULATION OF INTERPOLATION COEFFICIENTS
!----------------------------------------------------------

A(0) = Y(0)*X(1)*X(2)/((X(0)-X(1))*(X(0)-X(2))) + &
       Y(1)*X(0)*X(2)/((X(1)-X(0))*(X(1)-X(2))) + &
       Y(2)*X(0)*X(1)/((X(2)-X(0))*(X(2)-X(1)))

A(1) = -Y(0)*(X(1)+X(2))/((X(0)-X(1))*(X(0)-X(2)))  &
       -Y(1)*(X(0)+X(2))/((X(1)-X(0))*(X(1)-X(2)))  &
       -Y(2)*(X(0)+X(1))/((X(2)-X(0))*(X(2)-X(1)))

A(2) = Y(0)/((X(0)-X(1))*(X(0)-X(2))) + &
       Y(1)/((X(1)-X(0))*(X(1)-X(2))) + &
       Y(2)/((X(2)-X(0))*(X(2)-X(1)))


!----------------------------------------------------------
! 2. CALCULATION OF EXTREMUM
!----------------------------------------------------------

Xext = -A(1)/(2*A(2))
Yext = A(0) + Xext*(A(1) + Xext*A(2))


End Subroutine Extremum



End Module Interpolation
