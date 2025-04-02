!
!+ GNSS Radio occultation operator: Default real type definition and constants
!
! $Id$
!
MODULE Defaults
!
! Description:
!   GNSS Radio occultation bending angle observation operator:
!   Default real type definition and useful constants.
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
!  Double: adjust arguments of Selected_Real_Kind for NEC SX;
!  set default precision for ECHAM fields to double
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
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


! Module Defaults
!
! Default real type definition and useful constants.
!----------------------------------------------------------
! (C) Copyright 1998, M. E. Gorbunov.
!
! Version |    Date     | Comment
!---------+-------------+--------------------
!   1.0   | 07 Oct 1997 | Original code.
!   2.0   | 08 Sep 1998 | Earth constants in a separate
!         |             |   module.
!   2.1   | 11 Sep 1998 | Parameter Single
!----------------------------------------------------------
Implicit None
Private
Public :: Single, Double, WorkPr, Pi, dtr, rtd
!----------------------------------------------------------
! Public Parameters:
!
Integer,      Parameter :: &
   Single =                &  ! Real kind of single precision
       Selected_Real_Kind(5,30)
!
Integer,      Parameter :: &
   Double =                &  ! Real kind of double precision
       Selected_Real_Kind(13,200)
!
Integer,      Parameter :: &
   WorkPr = Double            ! Double or Single (for ECHAM fields only)

!
Real(Double), Parameter :: &
   Pi     =                &  ! Value of Pi
       3.141592653589793238_Double
!
Real(Double), Parameter :: &
   dtr    = Pi/180.0_Double   ! Degree-to-radian conversion factor
!
Real(Double), Parameter :: &
   rtd    = 180.0_Double/Pi   ! Radian-to-degree conversion factor
!
!
end module Defaults
