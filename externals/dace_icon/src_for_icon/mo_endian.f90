!
!+ Utilities for little/big endian detection and conversion
!
! $Id$
!
MODULE mo_endian
!
! Description:
!   Utilities for little/big endian detection and conversion.
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
! Authors:
! Andreas Rhodin  DWD  2002-2006
!------------------------------------------------------------------------------

  !-------------
  ! used modules
  !-------------
  use mo_kind,          only: sp, i4, i8 ! precision kind parameters
  implicit none

  !----------------
  ! public entities
  !----------------
  private
  public :: little ! function:   returns .true. on little endian machine
  public :: flip   ! subroutine: converts little <-> big endian

  interface flip
    module procedure flip_r_sp   ! real*4
    module procedure flip_i_i4   ! integer*4
    module procedure flip_i_i8   ! integer*8
  end interface flip

contains
!------------------------------------------------------------------------------
  function little()
  logical :: little
    integer(i4)      ,parameter :: ii4 = ((49*256+50)*256+51)*256+52
    character(len=*) ,parameter :: cl  = '4321'
    little = (cl==transfer(ii4,cl))
  end function little
!------------------------------------------------------------------------------
  elemental subroutine flip_r_sp (x)
  real(sp) ,intent(inout) :: x
    character(len=4) :: c
    character(len=4) :: t
    c = transfer (x,c)
    t(1:1) = c(4:4)
    t(2:2) = c(3:3)
    t(3:3) = c(2:2)
    t(4:4) = c(1:1)
    x = transfer (t,x)
  end subroutine flip_r_sp
!------------------------------------------------------------------------------
  elemental subroutine flip_i_i4 (x)
  integer(i4) ,intent(inout) :: x
    character(len=4) :: c
    character(len=4) :: t
    c = transfer (x,c)
    t(1:1) = c(4:4)
    t(2:2) = c(3:3)
    t(3:3) = c(2:2)
    t(4:4) = c(1:1)
    x = transfer (t,x)
  end subroutine flip_i_i4
!------------------------------------------------------------------------------
  elemental subroutine flip_i_i8 (x)
  integer(i8) ,intent(inout) :: x
    character(len=8) :: c
    character(len=8) :: t
    c = transfer (x,c)
    t(1:1) = c(8:8)
    t(2:2) = c(7:7)
    t(3:3) = c(6:6)
    t(4:4) = c(5:5)
    t(5:5) = c(4:4)
    t(6:6) = c(3:3)
    t(7:7) = c(2:2)
    t(8:8) = c(1:1)
    x = transfer (t,x)
  end subroutine flip_i_i8
!------------------------------------------------------------------------------
end module mo_endian
