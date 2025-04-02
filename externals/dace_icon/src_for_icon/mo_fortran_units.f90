!
!+ keeps a list of unused Fortran unit numbers
!
MODULE mo_fortran_units
!
! Description:
!   This module keeps a list of unit numbers used in the application.
!   In order to avoid collisions the application should ask for
!   unit numbers by:
!       get_unit_number or
!       reserve_unit_number (unit) (if a specific unit number is required)
!   and release unit numbers by:
!       return_unit_number (unit)
!
!   Cray f90 Compiler:
!     0, 5, 6           : associated to stderr, stdin, stdout,
!                         but may be reassigned
!     0 - 99, 105 - 299 : available for user specification
!     100 - 104         : reserved for system use
!
!   NAG f90 Compiler:
!     0, 5, 6           : preconnected to stderr, stdin, stdout
!     63                : maximum unit number
!
!   This module:
!     1 - 4, 7 - 999    : available for user specification
!     5, 6              : reserved
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
! V1_31        2014-08-21 Harald Anlauf
!  Increase highest_unit_number to 999
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  MPIfM/DWD  2001-2007  original source
!------------------------------------------------------------------------------
implicit none
private

!----------------
! Public entities
!----------------
public :: get_unit_number,     &!  reserve a unit number
          return_unit_number,  &!  release a unit number
          reserve_unit_number, &!  reserve a specific unit number
          print_used_units,    &!  printout reserved unit numbers
          nin, nout, nerr

  integer, parameter :: highest_unit_number = 999       ! was: 63
  logical            :: initialized         = .false.

  type unit_number_list_type
    logical, dimension (highest_unit_number) :: free
  end type unit_number_list_type

  type (unit_number_list_type) :: unit_number_list = &
        unit_number_list_type (.false.)

  integer, save :: nin  = 5     ! standard input stream
  integer, save :: nout = 6     ! standard output stream
  integer, save :: nerr = 0     ! error output stream

!==============================================================================
contains
!==============================================================================
  subroutine init_unit_number_list
    !-------------------------------
    ! Initialize unit_number_list:
    ! all units besides 5,6 are free
    !-------------------------------
    unit_number_list% free    = .true.
    unit_number_list% free(5) = .false.
    unit_number_list% free(6) = .false.
    initialized               = .true.
  end subroutine init_unit_number_list
!------------------------------------------------------------------------------
  function reserve_unit_number (unit) result (ok)
  !-------------------------------------------------
  ! Reserves unit number 'unit' for the application
  ! Returns .true. if 'unit' is not in use
  !-------------------------------------------------
  integer :: unit
  logical :: ok
    if (.not. initialized) call init_unit_number_list
    select case (unit)
    case (1:highest_unit_number)
      ok = unit_number_list% free (unit)
      unit_number_list% free (unit) = .false.
    case default
      ok = .false.
    end select
  end function reserve_unit_number
!------------------------------------------------------------------------------
  subroutine return_unit_number (unit)
  !-------------------------------------------------
  ! Returns unit number 'unit' if not used further.
  ! aborts if 'unit' has never been required before.
  !-------------------------------------------------
  integer, intent (inout) :: unit
    if (.not. initialized) call init_unit_number_list
    select case (unit)
    case (1:highest_unit_number)
      if (unit_number_list% free (unit) ) then
        print *,'WARNING in return_unit_number'
        print *,'  unit number',unit,'was never used'
      endif
      unit_number_list% free (unit) = .true.
    case default
      print *,'WARNING in return_unit_number'
      print *,'  wrong unit number',unit
    end select
    unit = -1
  end subroutine return_unit_number
!------------------------------------------------------------------------------
  function get_unit_number () result (unit)
  !--------------------------------------------
  ! Returns a free unit number
  ! Returns '0' if all unit numbers are in use
  !--------------------------------------------
  integer :: unit
    integer :: i
    if (.not. initialized) call init_unit_number_list
    unit = -1
    do i = highest_unit_number, 1, -1
      if (unit_number_list% free (i)) then
        unit = i
        unit_number_list% free (i) = .false.
        exit
      endif
    end do
  end function get_unit_number
!------------------------------------------------------------------------------
  subroutine print_used_units
    integer :: i
    if (.not. initialized) call init_unit_number_list
    print *,'used unit numbers:'
    do i = 1, highest_unit_number
      if (.not. unit_number_list% free(i)) write(*,'(i4)',advance='NO') i
    end do
    print *
  end subroutine print_used_units
!------------------------------------------------------------------------------
end module mo_fortran_units

