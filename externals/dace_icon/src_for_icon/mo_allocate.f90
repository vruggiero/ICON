!
!+ keep track of the allocation level in a function call hierarchy
!
! $Id$
!
MODULE mo_allocate
!
! Description:
!   As Fortran 95 does not support allocatable array components of derived
!   types, pointer components must be used in many cases.  If a variable
!   of derived type goes out of scope these components must be deallocated
!   explicitly. The module variable call_level keeps track of the level in
!   a function call hierarchy so that a function body is able to determine
!   if arguments goe out of scope so that their components may be
!   deallocated properly.
!   The following example shows the usage of this module:
!
!  FUNCTION Y (X)
!  TYPE(T) :: X, Y
!    CALL ENTER_FUNCTION
!    X % ALLOCATION_LEVEL = CALL_LEVEL
!    ALLOCATE (X % POINTER)
!    ...
!    IF (X % ALLOCATION_LEVEL == CALL_LEVEL) DEALLOCATE (X % POINTER)
!    CALL LEAVE_FUNCTION
!  END FUNCTION Y
!
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
! Language: Fortran 95.
! Software Standards:
!
!------------------------------------------------------------------------------
  implicit none
  private
  !----------------
  ! Public entities
  !----------------
  public :: call_level
  public :: enter_function
  public :: leave_function

  !-----------------
  ! module variables
  !-----------------
  integer :: call_level = 0

  !----------------------
  ! contained subroutines
  !----------------------
contains

  subroutine enter_function
    call_level = call_level + 1
  end subroutine enter_function

  subroutine leave_function
    call_level = call_level - 1
  end subroutine leave_function

end module mo_allocate
