!
!+ GNSS Radio occultation operator: solution of normal differential equations
!
MODULE Dif_Equations
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Numerical solution of normal differential equations.
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
! Module Dif_Equations
!
! Numerical solution of normal differential equations.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 29 Sep 1998 | Original code.
!----------------------------------------------------------
! Modules used:
!
Use Defaults, only: &
! Imported Parameters:
    Double
!----------------------------------------------------------
Implicit None
Private
Public :: Runge_Kutta
!----------------------------------------------------------
!
Contains

!==========================================================
Subroutine Runge_Kutta &
  (F,      & ! <-- Right-part function
   Dt,     & ! <-- Integration step
   t,      & ! <-> Time variable
   X)        ! <-> Dynamic variable vector
!
! A step of numerical integration of a systme
! of differential equations X'=F(X,t).
!----------------------------------------------------------
! Method:
!   Runge - Kutta scheme of 4th order.
!   G. Korn and T. Korn, Handbook on mathematics.
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
Interface
   Function F(t, X)  ! Right-part function
      Use Defaults, only: &
         Double
      Real(Double), Intent(In) :: t
      Real(Double), Intent(In) :: X(1:)
      Real(Double)             :: F(Size(X))
   End Function F
End Interface
!
Real(Double), Intent(In) :: &
   Dt    ! Time integration step
!
! Inout arguments:
!
Real(Double), Intent(InOut) :: &
   t     ! Time variable
!
Real(Double), Intent(InOut) :: &
   X(1:) ! Dynamic variable vector
!----------------------------------------------------------
! Local Arrays:
!
Real(Double)  :: &
   K1(Size(X)), K2(Size(X)), K3(Size(X)), K4(Size(X))
!----------------------------------------------------------


K1(:) = F(t,        X(:))
K2(:) = F(t + Dt/2, X(:) + K1(:)*Dt/2)
K3(:) = F(t + Dt/2, X(:) + K2(:)*Dt/2)
K4(:) = F(t + Dt,   X(:) + K3(:)*Dt)

X(:)  = X(:) + (K1(:) + 2*K2(:) + 2*K3(:) + K4(:))*Dt/6
t     = t + Dt


End Subroutine Runge_Kutta




End Module Dif_Equations


