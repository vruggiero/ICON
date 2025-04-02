!
!+ Trace memory usage
!
MODULE mo_mem_usage
!
! Description:
!   Trace memory usage.
!   Bookkeeping of memory allocated and deallocated in various parts
!   of the program.
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
! V1_20        2012-06-18 Andreas Rhodin
!  changed comments and names of public routines in module mo_p_output
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2006
!------------------------------------------------------------------------------
use mo_kind,       only: i8                 ! integer*8 kind parameter
use mo_mpi_dace,   only: dace,             &! MPI group info
                         p_sum,p_max,p_min  ! MPI sum, MPI max
use mo_p_output,   only: oline,            &! output line buffer
                         iol,              &! index of next line to write
                         nextline,         &! routine to increment line number
                         flush_buf          ! routine to write buffer
use mo_dec_matrix, only: dec_matrix_mem     ! trace memory usage
use mo_t_bg_err_op,only: cov_mem            ! trace memory usage
implicit none

private
public :: trace_mem_usage !  report memory usage
public :: print_mem_usage !  print  memory usage
public :: tmu             !  memory usage report enable flag

integer     :: count     = 1  !  count calls to subroutine trace_mem_usage
integer     :: tmu       = 0  !  memory usage report enable flag
integer     :: walls     = 0  !  wall time start
integer     :: wallm     = 0  !  maxumum wall count value
integer     :: rate      = 0  !  wall count rate
integer(i8) :: bytes_all = 0

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
  subroutine trace_mem_usage (comment)
  character(len=*) ,intent(in) :: comment

    integer(i8)                 :: bytes

    if (tmu <= 0) return

    bytes     = 0

    call dec_matrix_mem (bytes, count, tmu, comment)
    call print_mem_usage (comment, 'mo_dec_matrix', bytes)

    call cov_mem (bytes)
    call print_mem_usage (comment, 't_bg_err_op', bytes)

    call print_mem_usage (comment, 'all', bytes_all)

    count     = count + 1
    bytes_all = 0

  end subroutine trace_mem_usage
!-----------------------------------------------------------------------------
  subroutine print_mem_usage (comment1, comment2, bytes)
  character(len=*) ,intent(in) :: comment1
  character(len=*) ,intent(in) :: comment2
  integer(i8)      ,intent(in) :: bytes

    integer                     :: wall
    integer(i8)                 :: bytes_sum, bytes_max, bytes_min
    character(len=*), parameter :: fm2 = '(2i8,1x,2(a16,2x),a4,1x,i21," #tmu")'
    character(len=*), parameter :: fmt = '(2i8,1x,2(a16,2x),i4,1x,i21," #tmu")'

    if (tmu <= 0) return

    if (rate==0) then
      call system_clock (count      = walls)
      call system_clock (count_rate = rate)
      call system_clock (count_max  = wallm)
      if(rate==0) rate = 1
    endif
    call system_clock (count=wall)
    wall = (wall-walls)
    if (wall < 0) wall = wall + wallm
    wall = wall / rate

    bytes_min = p_min (bytes)
    bytes_max = p_max (bytes)
    bytes_sum = p_sum (bytes)
    if (bytes_sum > 0) then
      if (tmu >= 2) then
        call nextline
        write (oline(iol),fmt) count, wall, comment1, comment2, dace% pe, bytes
        call flush_buf
      endif
      if (dace% lpio) then
        write (6,fm2) count, wall, comment1, comment2, 'min' , bytes_min
        write (6,fm2) count, wall, comment1, comment2, 'max' , bytes_max
        write (6,fm2) count, wall, comment1, comment2, 'mean', bytes_sum &
                                                             / dace% npe
        write (6,fm2) count, wall, comment1, comment2, 'all' , bytes_sum
      endif
    endif
    bytes_all = bytes_all + bytes

  end subroutine print_mem_usage
!-----------------------------------------------------------------------------
end module mo_mem_usage


