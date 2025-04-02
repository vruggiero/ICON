!
!+ GNSS Radio occultation operator: adjoints for spline interpolation
!
MODULE Interpolation_adj
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Linear adjoints for spline interpolation routines.
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
!  Spline_Double_adj: change allocatable work arrays to automatic
!  Optimizations for NEC SX-9
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_14        2011/11/08 Harald Anlauf
!  Cleanup and minor optimization
! V1_19        2012-04-16 Harald Anlauf
!  interpolation, interpolation_adj: optimize for sxf90 rev.430+
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
! Module Interpolation_adj
!
! Linear adjoints for spline interpolation routines.
!----------------------------------------------------------
! (C) Copyright 2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 06 Apr 2000 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double !, Single
!----------------------------------------------------------
Implicit None
    Private
    Public :: Init_Spline_adj, Spline_adj
!----------------------------------------------------------
! Interfaces:
!
Interface Init_Spline_adj
   Module Procedure Init_Spline_Double_adj
!   Module Procedure Init_Spline_Single
End Interface
!
Interface Spline_adj
   Module Procedure Spline_Double_adj
!   Module Procedure Spline_Single
End Interface
!----------------------------------------------------------
Contains


!==========================================================
Subroutine Init_Spline_Double_adj &
  (x,      & ! <-- Argument grid
   f,      & ! <-- Gridded function
   d2,     & ! --> 2nd derivative of spline
   d2_x,   & ! --> d(d2(i))/d(x(j))
   d2_f,   & ! --> d(d2(i))/d(f(j))
   ErrCode)  ! ~~> Error code
!
! Adjoint version of routine for calculation of
! 2nd derivative of spline.
!----------------------------------------------------------
! Method:
!   Adjoint for drive-through solution.
!----------------------------------------------------------
! (C) Copyright 1998-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   3.0   | 23 Dec 1998 | Basic non-adjoint version.
!   1.0   | 06 Apr 2000 | Adjoint version.
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)       :: &
   x(1:)       ! Argument grid
               ! x(i) must be monotonous.
Real(Double), Intent(In)       :: &
   f(1:)       ! Gridded function
!
! Output arguments:
!
Real(Double), Intent(Out)      :: &
   d2(1:)      ! 2nd derivative of spline
Real(Double), Intent(Out)      :: &
   d2_x(:,:)   ! d(d2(i))/d(x(j))
Real(Double), Intent(Out)      :: &
   d2_f(:,:)   ! d(d2(i))/d(f(j))
Integer, Intent(Out), Optional :: &
   ErrCode     ! Error code:
               !   0 - no error
               !   1 - x is not mononous
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   dfl, dfr, df,  &
   dxl, dxr, dx,  &
   a, b, c
Integer      :: &
   i,          & ! Array index
   N             ! Number of array elements
!
! Local Arrays:
!
Real(Double), Dimension(Size (x)) :: &
   d1,         & ! Work array for drive-through
   df_x,       & ! d(df)/d(x)
   a_x,        & ! d(a)/d(x)
   b_x,        & ! d(b)/d(x)
   c_x,        & ! d(c)/d(x)
   df_f          ! d(df)/d(x)
Real(Double), Dimension(Size (x),Size (x)) :: &
   d1_x          ! d(d1)/d(x)
!----------------------------------------------------------


!----------------------------------------------------------
! 0. INITIALIZATION
!----------------------------------------------------------

!--- 0.1. Checking parameters

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


!--- 0.2. Array allocation

!Allocate(d1(N))
!Allocate(d1_x(N,N), df_x(N), a_x(N), b_x(N), c_x(N))
!Allocate(df_f(N))


!----------------------------------------------------------
! 1. DRIVE-THROUGH CALCULATION OF SPLINE COEFFICIENTS
!----------------------------------------------------------

!--- 1.1. Boundary initialization

d2(1) = 0.0_Double
d2(N) = 0.0_Double
d1(1) = 0.0_Double

d2_x(1,:) = 0.0_Double
d2_x(N,:) = 0.0_Double
d1_x(1,:) = 0.0_Double
d2_f(1,:) = 0.0_Double
d2_f(N,:) = 0.0_Double


!--- 1.2. Drive forward

Do i = 2, N-1

   !--- 1.2.1. Non-adjoint calculations

   dxl   = x(i) - x(i-1)
   dfl   = (f(i) - f(i-1))/dxl
   dxr   = x(i+1) - x(i)
   dfr   = (f(i+1) - f(i))/dxr
   dx    = dxl + dxr
   df    = (dfr - dfl)/dx
   a     = dxl/(2*dx)
   b     = dxr/(2*dx)
   c     = 1 + a*d1(i-1)
   d1(i) = -b/c
   d2(i) = (3*df - a*d2(i-1))/c

   !--- 1.2.2. Adjoint calculations: x derivatives

   df_x(:)   = 0
   df_x(i-1) = (df - dfl/dxl)/dx
   df_x(i)   = (dfr/dxr + dfl/dxl)/dx
   df_x(i+1) = -(df + dfr/dxr)/dx

   a_x(:)    = 0
   a_x(i-1)  = (2*a - 1)/(2*dx)
   a_x(i)    = 1/(2*dx)
   a_x(i+1)  = -a/dx

   b_x(:)    = 0
   b_x(i-1)  = b/dx
   b_x(i)    = -1/(2*dx)
   b_x(i+1)  = (1 - 2*b)/(2*dx)

!CDIR ARRAYCOMB
   c_x(:)    = a_x(:)*d1(i-1) + a*d1_x(i-1,:)

   d1_x(i,:) = -b_x(:)/c + b*c_x(:)/c**2

   d2_x(i,:) = (3*df_x(:) - a_x(:)*d2(i-1) - a*d2_x(i-1,:))/c - &
               c_x(:)*(3*df - a*d2(i-1))/c**2
!CDIR END ARRAYCOMB


   !--- 1.2.3. Adjoint calculations: f derivatives

   df_f(:)   = 0
   df_f(i-1) = 1/(dxl*dx)
   df_f(i)   = (-1/dxr - 1/dxl)/dx
   df_f(i+1) = 1/(dxr*dx)

   d2_f(i,:) = (3*df_f(:) - a*d2_f(i-1,:))/c

End Do


!--- 1.3. Drive back

Do i = N-1, 2, -1

   !--- 1.3.1. Non-adjoint calculations

   d2(i) = d1(i)*d2(i+1) + d2(i)

   !--- 1.3.2. Adjoint calculations: x derivatives

   d2_x(i,:) = d1_x(i,:)*d2(i+1) + d1(i)*d2_x(i+1,:) + d2_x(i,:)

   !--- 1.3.3. Adjoint calculations: f derivatives

   d2_f(i,:) = d1(i)*d2_f(i+1,:) + d2_f(i,:)

End Do


!----------------------------------------------------------
! 2. MEMORY DEALLOCATION
!----------------------------------------------------------

!Deallocate(d1)
!Deallocate(d1_x, df_x, a_x, b_x, c_x)
!Deallocate(df_f)


End Subroutine Init_Spline_Double_adj


!==========================================================
Subroutine Spline_Double_adj &
  (x,          & ! <-- Argument grid
   f,          & ! <-- Gridded function
   d2,         & ! <-- 2nd derivative of spline
   d2_x,       & ! <-- d(d2(i))/d(x(j))
   d2_f,       & ! <-- d(d2(i))/d(f(j))
   x_int,      & ! <-- Interpolation point
   f_int,      & ! --> Interpolated function value
   f_int_x,    & ! --> d(f_int)/d(x(i))
   f_int_f,    & ! --> d(f_int)/d(f(i))
   fd_int,     & ! --> Interpolated 1st derivative
   fd_int_x,   & ! --> d(fd_int)/d(x(i))
   fd_int_f,   & ! --> d(fd_int)/d(x(i))
   fd2_int)      ! ~~> Interpolated 2nd derivative
!
! Adjoint version of spline interpolation of
! a gridded function.
!----------------------------------------------------------
! Method:
!   Searching for grid interval containing the
!   interpolation point and summation of polynomial
!   with given spline-coefficients.
!----------------------------------------------------------
! (C) Copyright 1998-2000, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   4.0   | 15 Jan 1999 | Basic non-adjoint version.
!   1.0   | 19 Apr 2000 | Adjoint version.
!----------------------------------------------------------
! Modules used:
!
Use Interpolation, only: &
! Imported Routines:
    Seek_Index
!----------------------------------------------------------
Implicit None
!----------------------------------------------------------
! Input arguments:
!
Real(Double), Intent(In)     :: &
   x(:)         ! Argument grid
                ! x(i) must be monotonous.
Real(Double), Intent(In)     :: &
   f(:)         ! Gridded function
Real(Double), Intent(In)     :: &
   d2(:)        ! 2nd derivative of spline
Real(Double), Intent(In)     :: &
   d2_x(:,:)    ! d(d2)/d(x)
Real(Double), Intent(In)     :: &
   d2_f(:,:)    ! d(d2)/d(f)
Real(Double), Intent(In)     :: &
   x_int        ! Interpolation point
!
! Output arguments:
!
Real(Double), Intent(Out)    :: &
   f_int        ! Interpolated function value
Real(Double), Intent(Out)    :: &
   f_int_x(:)   ! d(f_int)/d(x)
Real(Double), Intent(Out)    :: &
   f_int_f(:)   ! d(f_int)/d(f)
Real(Double), Intent(Out)    :: &
   fd_int       ! Interpolated 1st derivative
Real(Double), Intent(Out)    :: &
   fd_int_x(:)  ! d(fd_int)/d(x)
Real(Double), Intent(Out)    :: &
   fd_int_f(:)  ! d(fd_int)/d(x)
Real(Double), Intent(Out), Optional    :: &
   fd2_int      ! Interpolated 2nd derivative
!----------------------------------------------------------
! Local Scalars:
!
Real(Double) :: &
   a1, a2, a3,  & ! Polynomial coefficients
   dx,          & ! Grid interval
   dx_t,        & ! Grid-point to interpolation-point distance
   x_t,         & ! Interpolation point projected to grid extent
   Dir            ! Direction of argument change
Integer  :: &
   i,           & ! Array index
   N              ! Number of data
Integer  :: i_int ! Interpolation interval index
!
! Local Arrays:
!
Real(Double), Dimension(Size (x)) :: &
   a1_x,        & ! d(a1)/d(x)
   a2_x,        & ! d(a2)/d(x)
   a3_x,        & ! d(a3)/d(x)
   a1_f,        & ! d(a1)/d(f)
   a2_f,        & ! d(a2)/d(f)
   a3_f           ! d(a3)/d(f)
!----------------------------------------------------------


!----------------------------------------------------------
! 1. LOCATION OF INTERPOLATION POINT INSIDE GRID
!----------------------------------------------------------

!--- 1.1. Determination of data size

N     = Size(x)


!--- 1.2. Grid interval location

x_t   = Min(Max(x_int, Min(x(1),x(N))), Max(x(1),x(N)))
Dir   = Sign(1.0_Double, x(N)-x(1))

i_int = Seek_Index(x, x_t)
i     = Max(i_int, 1)


!--- 1.3. Memory allocation

!Allocate(a1_x(N), a2_x(N), a3_x(N))
!Allocate(a1_f(N), a2_f(N), a3_f(N))


!----------------------------------------------------------
! 2. CALCULATION OF INTERPOLATION COEFFICIENTS
!----------------------------------------------------------

!--- 2.1. Non-adjoint calculations

dx    = x(i+1) - x(i)
a2    = d2(i)/2
a3    = (d2(i+1) - d2(i))/(6*dx)
a1    = (f(i+1) - f(i))/dx - dx*(a2 + dx*a3)


!--- 2.2. Adjoint calculations: x derivatives

a2_x(:)   = d2_x(i,:)/2

a3_x(:)   = (d2_x(i+1,:) - d2_x(i,:))/(6*dx)
a3_x(i)   = a3_x(i) + a3/dx
a3_x(i+1) = a3_x(i+1) - a3/dx

a1_x(:)   = -dx*(a2_x(:) + dx*a3_x(:))
a1_x(i)   = a1_x(i) + &
            (f(i+1) - f(i))/dx**2 + a2 + 2*dx*a3
a1_x(i+1) = a1_x(i+1) - &
            (f(i+1) - f(i))/dx**2 - a2 - 2*dx*a3


!--- 2.3. Adjoint calculations: f derivatives

!CDIR ARRAYCOMB
a2_f(:)   = d2_f(i,:)/2

a3_f(:)   = (d2_f(i+1,:) - d2_f(i,:))/(6*dx)

a1_f(:)   = -dx*(a2_f(:) + dx*a3_f(:))
!CDIR END ARRAYCOMB

a1_f(i)   = a1_f(i) - 1/dx
a1_f(i+1) = a1_f(i+1) + 1/dx



!----------------------------------------------------------
! 3. CALCULATION OF INTERPOLATED VALUE
!----------------------------------------------------------

!--- 3.1. Non-adjoint calculations

dx_t   = x_t - x(i)
fd_int = a1 + dx_t*(2*a2 + 3*dx_t*a3)
f_int  = f(i) + dx_t*(a1 + dx_t*(a2 + dx_t*a3)) + &
         fd_int*(x_int - x_t)

If (Present(fd2_int)) then
   fd2_int  = 2*a2 + 6*dx_t*a3
End if


!--- 3.2. Adjoint calculations: x derivatives

!CDIR ARRAYCOMB
fd_int_x(:) = a1_x(:) + dx_t*(2*a2_x(:) + 3*dx_t*a3_x(:))

f_int_x(:)  = dx_t*(a1_x(:) + dx_t*(a2_x(:) + dx_t*a3_x(:))) + &
              fd_int_x(:)*(x_int - x_t)
!CDIR END ARRAYCOMB

fd_int_x(i) = fd_int_x(i) - 2*a2 - 6*dx_t*a3

f_int_x(i)  = f_int_x(i) - &
              a1 - dx_t*(2*a2 + 3*dx_t*a3)

If (Dir*x_int > Dir*x(N)) then
   fd_int_x(N) = fd_int_x(N) + 2*a2 + 6*dx_t*a3
   f_int_x(N)  = f_int_x(N) + &
                 a1 + dx_t*(2*a2 + 3*dx_t*a3) - &
                 fd_int
Else If (Dir*x_int < Dir*x(1)) then
   fd_int_x(1) = fd_int_x(1) + 2*a2 + 6*dx_t*a3
   f_int_x(1)  = f_int_x(1) + &
                 a1 + dx_t*(2*a2 + 3*dx_t*a3) - &
                 fd_int
End If


!--- 3.3. Adjoint calculations: f derivatives

!CDIR ARRAYCOMB
fd_int_f(:) = a1_f(:) + dx_t*(2*a2_f(:) + 3*dx_t*a3_f(:))

f_int_f(:)  = dx_t*(a1_f(:) + dx_t*(a2_f(:) + dx_t*a3_f(:))) + &
              fd_int_f(:)*(x_int - x_t)
!CDIR END ARRAYCOMB

f_int_f(i)  = f_int_f(i) + 1


!----------------------------------------------------------
! 4. MEMORY DEALLOCATION
!----------------------------------------------------------

!Deallocate(a1_x, a2_x, a3_x)
!Deallocate(a1_f, a2_f, a3_f)


End Subroutine Spline_Double_adj



End Module Interpolation_adj
