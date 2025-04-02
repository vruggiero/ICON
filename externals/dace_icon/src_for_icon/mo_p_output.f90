!
!+ Handle parallel job output (to unit 6, stdout)
!
MODULE mo_p_output
!
! Description:
!   Routines to handle parallel job output (to unit 6, stdout).
!   Usage:
!   Each line to be written to the output must be written to a buffer
!   (array of character variables) on each processor element by the
!   sequence:
!      CALL NEXTLINE                          ! increase line buffer index
!      WRITE (OLINE(IOL),"formatstring") line ! write line
!
!   Or by
!      CALL ADD_LINE (line)
!   Finally all lines are written to unit 6, PE by PE:
!      CALL FLUSH_BUF
!   In order to write lines only on one PE (e.g. header lines) routine
!   ADD_LINE_PIO may be used (equivalent to call CALL NEXTLINE; WRITE
!   (OLINE(IOL).. on the I/O processor element only.
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
! V1_5         2009/05/25 Harald Anlauf
!  write_p: optimize MPI communication
! V1_7         2009/08/24 Harald Anlauf
!  increase 'siz_oline' to 30000
! V1_9         2010/04/20 Harald Anlauf
!  write_pio_only: write info message when output buffer overflows
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Harald Anlauf
!  write_p: use p_gatherv
! V1_20        2012-06-18 Andreas Rhodin
!  changed comments and names of public routines
! V1_42        2015-06-08 Andreas Rhodin
!  mo_p_output: increase 'len_oline' to 256
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2002-2007
!==============================================================================
  !-------------
  ! Modules used
  !-------------
  use mo_mpi_dace, only: dace, p_gather, p_gatherv
  implicit none
  !----------------
  ! Public entities
  !----------------
  private
  public :: oline          ! flush buffer (array of strings)
  public :: iol            ! index of next line to write to buffer
  public :: nextline       ! routine to increment line number
  public :: flush_buf      ! routine to gather and write buffer to stdout
  public :: flush_buf_pio  ! as flush_buf, but only if dace% lpio
  public :: add_line       ! routine to add a line to the output buffer
  public :: add_line_pio   ! add a line to output buffer for dace% lpio only
  !-----------------
  ! Module variables
  !-----------------
  integer,parameter        :: len_oline =   256  ! length of one line
  integer,parameter        :: siz_oline = 30000  ! size of array of lines
  integer                  :: iol       =     0  ! index of last line used
  character(len=len_oline) :: oline(siz_oline+1) ! line buffer
  save                     :: oline
contains
!==============================================================================
#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif


  subroutine flush_buf
  !-------------------------------------------
  ! write and clear the output buffer.
  ! the routine must be called on all PEs.
  ! output is gathered and written by PE p_io.
  !-------------------------------------------
  integer                  :: pe, i, n, k
  integer                  :: lines(0:dace% npe-1)
  logical                  :: oflow(0:dace% npe-1)     ! Buffer overflow on pe?

    integer                  :: nl                      ! No.lines for gatherv
    character(len=len_oline), allocatable :: rbuf(:)    ! Receive buffer
    FTRACE_BEGIN('flush_buf')
    !-----------------------------------------
    ! Flush line buffer on I/O processor first
    !-----------------------------------------
    n = iol                     ! Save line count
    call flush_buf_pio ()
    !------------------------------------------
    ! Gather output lines from other processors
    !------------------------------------------
    call p_gather (iol, lines, root=dace% pio)
    if (dace% lpio) then
      nl = sum (lines)
    else
      nl = 1                   ! Senders need dummy buffer for p_gatherv
    end if
    allocate (rbuf(nl))
    call p_gatherv (sendbuf=oline(1:iol), recvbuf=rbuf, root=dace% pio)
    if (dace% lpio) then
      oflow       = lines > siz_oline
      oflow(dace% pio) = n     > siz_oline
      i = 0
      do pe = 0, dace% npe - 1
        do k = 1, lines(pe)
          i = i + 1
          write(6,'(a)') trim (rbuf(i))
        end do
        !---------------------------------------------
        ! Mark truncation of buffer in the right place
        !---------------------------------------------
        if (oflow(pe)) &
          write(6,'("!+!+!+! flush_buf: output buffer overflow on PE ",i5)') pe
      end do
    end if
    iol = 0
    FTRACE_END('flush_buf')
  end subroutine flush_buf
!------------------------------------------------------------------------------
  subroutine flush_buf_pio
  !----------------------------------------------------------------
  ! Write and clear the output buffer on PE p_io only.
  ! This routine does not gather the buffer content from other PEs.
  ! Recommended if p_io writes a lot to prevent buffer overflow.
  !----------------------------------------------------------------
    integer :: i
    if (dace% lpio) then
      do i = 1, iol
        write(6,'(a)') trim(oline(i))
      end do
      if (iol > siz_oline) &
        write(6,'("!+!+!+! flush_buf: output buffer overflow on PE ",i5)') dace% pe
      iol = 0
    endif
  end subroutine flush_buf_pio
!------------------------------------------------------------------------------
  subroutine nextline
  !-------------------------------------------------------
  ! increment the current line number of the output buffer
  ! to be used (written) by the actual processor element.
  !-------------------------------------------------------
    iol = min (iol+1, siz_oline+1)
    oline(iol) = ''
  end subroutine nextline
!------------------------------------------------------------------------------
  subroutine add_line (line)
  !--------------------------------
  ! Add a line to the output buffer
  !--------------------------------
  character(len=*) ,intent(in) :: line
    call nextline
    oline(iol) = line
  end subroutine add_line
!------------------------------------------------------------------------------
  subroutine add_line_pio (line)
  !----------------------------------------------
  ! Add a line to the output buffer on p_io only.
  !----------------------------------------------
  character(len=*) ,intent(in) :: line
    if (dace% lpio) then
      call nextline
      oline(iol) = line
    endif
  end subroutine add_line_pio
!==============================================================================
end module mo_p_output
