!
!+ bookkeeping of cpu-time and wall time
!
MODULE mo_cpu_time
!
! Description:
!   Stop cpu-time, wall-time, and memory used for program segments.
!   Print time and memory requirements of program segments.
!
! Current Maintainer: DWD, Harald Anlauf
!    phone: +49 69 8062 4941
!    fax:   +49 69 8062 3721
!    email: harald.anlauf@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_1         2008/11/05 Andreas Rhodin
!  First operational 3D-Var release
! V1_4         2009/03/26 Andreas Rhodin
!  Flush output buffer
! V1_5         2009/05/25 Harald Anlauf
!  Write PE id for cpu-time and rss info (experimental)
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  change printout: write seperator (-----) below timings
! V1_20        2012-06-18 Andreas Rhodin
!  changed comment lines
! V1_22        2013-02-13 Andreas Rhodin
!  prepare for emacs outline mode
! V1_23        2013-03-26 Andreas Rhodin
!  option to check synchronisation of PEs
! V1_28        2014/02/26 Harald Anlauf
!  RSS summary: use maxval of measured rss to improve experience on Linux
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Autors:
! Andreas Rhodin  MPIfM/DWD  2001-2008  original code
! Harald Anlauf   DWD        2008       trace maxrss, handle counter wraparound
!------------------------------------------------------------------------------
!-------------
! Modules used
!-------------
use mo_mpi_dace,   only: dace, p_barrier, p_gather, p_bcast
use mo_kind,       only: dp, i8        ! double precision kind parameters
use mo_exception,  only: finish        ! abort in case of error
use mo_system,     only: maxrss,      &! Max. resident set size (if supported)
                         flush         ! flush output buffer
use mo_run_params, only: barrier,     &! call mpibarrier for diagnostics
                         check_sync,  &! check syncronisation
                         debug_rss,   &! Memory usage debugging level
                         node_id       ! Node index for rank

#if defined (_CRAYFTN) || defined (__CRAYXC)
use mo_system,     only: get_hugepages ! Hugepage usage
#endif

implicit none
save

!----------------
! Public entities
!----------------
private
public :: stop_time   ! stop time and memory used by a program segment
public :: print_times ! print summary of time and memory requirementa
public :: iclk        ! integer kind type used for arguments of SYSTEM_CLOCK
public :: detail      ! Detailed output including pes

!-------------------------
! Private module variables
!-------------------------

!-------------------------------------------
! Kind type of arguments to system_clock ():
!-------------------------------------------
#if (defined (__SUNPRO_F90) || defined (__SUNPRO_F95) || \
     defined (NAGFOR)       || defined (__NEC__)      || \
     defined (__PGI)        || defined (__FLANG)      || \
     defined (__GFORTRAN__) || defined (__INTEL_COMPILER))
!--------------------------------------------------
! Sun f95 offers a high resolution, making multiple
! wrap-arounds likely with default integers.
! NAG and NEC recommend int64.
!--------------------------------------------------
integer, parameter         :: iclk = i8         ! Use long integers
#else
integer, parameter         :: iclk = kind (1)   ! Default integer kind
#endif

integer ,parameter         :: maxtimes = 1000
real                       :: cpus  (maxtimes)
real                       :: cpue  (maxtimes)
integer(iclk)              :: walls (maxtimes)
integer(iclk)              :: walle (maxtimes)
character(len=48)          :: names (maxtimes)
integer                    :: rss   (maxtimes) = 0
integer                    :: huge2m(maxtimes) = 0
integer                    :: smap2m(maxtimes) = 0
integer                    :: ntimes = 0
integer(iclk)              :: rate   = 0
character(len=*),parameter :: form1 = '(5f8.2,2i7,i8, 2x,a)'
character(len=*),parameter :: form0 = &
'(/" cpu(min)   (max)  (mean)   wall percent minrss maxrss  sumrss  task")'

character(len=*),parameter :: form2 = &
'(/" cpu(min)        (max)       (mean)   wall percent maxrss       task")'
character(len=*),parameter :: form3 = &
     '(2(f8.2," [",i2,"]"),3f8.2,i7," [",i2,"]",2x,a)'

logical                    :: detail = .false.        ! Detailed info including pes
logical                    :: ftrace_region = .false. ! ftrace_region open

!---------------------
! Contained subroutine
!---------------------
contains
!==============================================================================
  subroutine stop_time (name, outline)
  !---------------------------------------------------------------------
  ! Determine cpu-time, wall-time and memory requirements of sections of
  ! the program. Print out the requirements of the previous section.
  ! Print header line for emacs outline mode
  !---------------------------------------------------------------------
  character(len=*) ,intent(in)           :: name    ! description of segment
  character(len=*) ,intent(in) ,optional :: outline ! *s for emacs outline mode

    character(len=64) :: line
    integer(iclk)     :: wall
    real              :: cpu
    character(len=8)  :: c
    call cpu_time     (cpu)

!   call check_cpu_num ()       ! Validate process pinning
    if (barrier) call p_barrier
    call system_clock (count=wall)
    if (ntimes > 0 .and. ntimes <= maxtimes) then
      cpue  (ntimes) = cpu
      walle (ntimes) = wall
      rss   (ntimes) = maxrss ()
#if defined (_CRAYFTN) || defined (__CRAYXC)
      if (debug_rss > 0) &
           call get_hugepages (hugepages=huge2m(ntimes), smaps=smap2m(ntimes))
#endif
    endif
    if (ntimes < maxtimes) then
      ntimes = ntimes + 1
      names (ntimes) = name
      walls (ntimes) = wall
      call cpu_time (cpu)
      cpus  (ntimes) = cpu
    endif

#if defined (_FTRACE)
    if (ntimes == maxtimes ) call finish ('stop_time','ntimes reached maxtimes')
    if (ntimes > 1) then
      if (ftrace_region) call ftrace_region_end ( trim(names(ntimes-1)) )
    end if
    call ftrace_region_begin ( trim(name) )
    ftrace_region = .true.
#endif

    if (barrier .and. ntimes > 1) then
      if (dace% lpio) then
         if (.not. detail) then
            write (6,form0)
         else
            write (6,form2)
         end if
      end if
      !---------------------
      ! check syncronisation
      !---------------------
      if (check_sync) then
        line = name
        call p_bcast (line, dace% pio)
        if (name /= line) then
          write (0,*) dace% pe,'stop_time, out of sync: name='//trim(name)
          write (0,*) dace% pe,'stop_time, out of sync: pe0 ='//trim(line)
          call finish ('stop_time','out of sync: '//trim(line))
        endif
      endif
      !--------------
      ! print timings
      !--------------
      call print_line (ntimes-1, ntimes-1, names(ntimes-1))
      if (debug_rss>2) call print_rss_overview ()
      !------------------------------------
      ! print header for emacs outline mode
      !------------------------------------
      if (dace% lpio) then
        write (6,'(/a)') repeat('-',79)
        c = '**'
        if (present(outline)) c = outline
        if (name/='' .and. c/='') &
          write (6,'(/a," ",a/)') trim(c), trim(name)
      endif
    endif
  end subroutine stop_time
!------------------------------------------------------------------------------
  subroutine print_times
  !---------------------------------------------------------------------
  ! Print out summary of cpu-time, wall-time and memory requirements of
  ! the sections of the programm previously accumulated by calls to
  ! 'stop_time'
  !---------------------------------------------------------------------
    integer :: i
    if (dace% lpio) then
       if (.not. detail) then
          write (6,form0)
       else
          write (6,form2)
       end if
    end if
    do i=1,ntimes-1
      call print_line (i ,i ,names(i))
    end do
    if (ntimes > 1) call print_line (1, ntimes-1, 'Total')
    if (dace% lpio) call flush(6)
    if (debug_rss>0) call print_rss_overview (final=.true.)

#if defined (_FTRACE)
    if (ftrace_region) then
      ftrace_region = .false.
      call ftrace_region_end ( trim(names(ntimes)) )
    end if
#endif

  end subroutine print_times
!------------------------------------------------------------------------------
  subroutine print_line (i0 ,i1 ,name)
  integer          ,intent(in)           :: i0, i1
  character(len=*) ,intent(in)           :: name
    real          :: cpu ,wall ,percent, mi, ma, mn
    real(dp)      :: cpud(0:dace% npe-1)
    integer(iclk) :: count_max
    integer       :: sr, mir, mar, rssm, rssp(0:dace% npe-1)
    integer       :: pe_mi, pe_ma, pe_mar
    !---------------------------------
    ! determine count rate (wall time)
    !---------------------------------
    if (rate==0) then
      call system_clock (count_rate=rate)
      if(rate==0) rate = 1
    endif
    !---------------
    ! write one line
    !---------------
    cpu  =  cpue (i1) - cpus (i0)
    wall = (walle(i1) - walls(i0)) / real(rate)
    if (wall < 0) then
       !--------------------------------------------------------------------
       ! Try to handle simple wrap-around of integer counter of system clock
       !--------------------------------------------------------------------
       call system_clock (count_max=count_max)
       wall = (walle(i1) - walls(i0) + count_max + 1) / real(rate)
    end if
    rssm = maxval (rss(i0:i1))          ! max. rss measured on this pe
    call p_gather (real (cpu,dp), cpud, root=dace% pio)
    call p_gather (rssm         , rssp, root=dace% pio)
    if (dace% lpio) then
       mi     = minval (cpud)
       ma     = maxval (cpud)
       mn     =    sum (cpud) / dace% npe
       pe_mi  = minloc (cpud,1) - 1
       pe_ma  = maxloc (cpud,1) - 1
       mir    = minval (rssp)
       mar    = maxval (rssp)
       sr     =    sum (rssp)
       pe_mar = maxloc (rssp,1) - 1
       percent = 999.99
       if(wall/=0.) percent = 100. * mn/wall
       if (.not. detail) then
          write (6,form1) mi,ma,mn ,wall ,percent ,mir, mar, sr, trim(name)
       else
          write (6,form3) &
               mi,pe_mi, ma,pe_ma, mn ,wall ,percent, mar,pe_mar, trim(name)
       end if
    end if
  end subroutine print_line
!==============================================================================
  subroutine print_rss_overview (final)
    logical, optional, intent(in) :: final
    !------------------------------------------------
    ! Print memory usage (rss) overview for all nodes
    !------------------------------------------------
    integer :: i, nn, maxrss, srss, maxh2m, maxmap, smap, mh2m
    integer :: rss_(0:dace% npe-1)
    integer :: h2m_(0:dace% npe-1)
    integer :: map_(0:dace% npe-1)
    logical :: mask(0:dace% npe-1)
    logical :: lfinal

    if (ntimes == 0) return
    if (.not. allocated (node_id)) return

    lfinal = .false.; if (present (final)) lfinal = final

    if (dace% lpio) then
       write (6,'(/a)') " Memory usage overview at node level (id,total,pe)"
    endif
    maxrss = maxval (rss   (1:ntimes))          ! max. rss measured on this pe
    maxh2m = maxval (huge2m(1:ntimes))          ! max. hugepages(2M)
    maxmap = maxval (smap2m(1:ntimes))          ! max. smap(2M)     on this pe
    call p_gather (maxrss, rss_, root=dace% pio)
    call p_gather (maxh2m, h2m_, root=dace% pio)
    call p_gather (maxmap, map_, root=dace% pio)
    nn = maxval (node_id)
    if (dace% lpio) then
       do i = 1, nn
          mask = (node_id == i)
          srss = sum (rss_, mask)
          if (debug_rss > 1) then
             write (6,'(i4,i8," :",*(i6))') i, srss, pack (rss_, mask)
          else
             write (6,'(i4,i8)')            i, srss
          end if
       end do
       if (lfinal) write (6,'(a)') repeat('-',79)
    end if

#if defined (_CRAYFTN) || defined (__CRAYXC)
    if (dace% lpio) then
       if (.not. lfinal) write (6,'()')
       write (6,'(a)') " Hugepages usage overview at node level (id,total,pe)"
       do i = 1, nn
          mask = (node_id == i)
          smap = sum    (map_, mask)
          mh2m = maxval (h2m_, mask)
          if (debug_rss > 1) then
             write (6,'(i4,i8," :",*(i6))') i, smap, pack (map_, mask)
          else
             write (6,'(i4,i8)')            i, smap
          end if
          if (smap /= mh2m*2) then
             write (0,'(a,i4,i8," |",i8," :",*(i6))') &
                  "Hugepages mismatch!", i, mh2m, smap, pack (map_, mask)
          end if
       end do
       if (lfinal) write (6,'(a)') repeat('-',79)
    end if
#endif

    if (dace% lpio) call flush(6)
    call p_barrier ()
  end subroutine print_rss_overview
!==============================================================================
  subroutine get_cpu_num (processor)
    !-----------------------------------------------------
    ! Determine CPU on which the calling thread is running
    ! (This may be Linux/glibc specific, thus restricted.)
    !-----------------------------------------------------
    integer, intent(out) :: processor
#if (defined (__linux__) || defined (__GFORTRAN__)) && !defined (__NEC__)
    interface
       function cpu_num () bind(c,name="sched_getcpu")
         use iso_c_binding, only: c_int
         integer(c_int) :: cpu_num
       end function cpu_num
    end interface
    processor = cpu_num ()
#else
    processor = -1
#endif
  end subroutine get_cpu_num
!------------------------------------------------------------------------------
  subroutine check_cpu_num ()
    !-------------------------
    ! Validate process pinning
    ! (binding to fixed core)
    !-------------------------
    integer :: processor = -99
    integer :: i, j, mm, nn, cpumin, cpumax, newcpu
    integer :: cpu_(0:dace% npe-1)
    logical :: mask(0:dace% npe-1)
    integer :: rank(0:dace% npe-1)

    call get_cpu_num (newcpu)
    if (processor == -99) processor = newcpu
    if (processor ==  -1) return
    if (processor /= newcpu) then
       write(0,'("pe=",i4,"  cpu number change:",i4,"  >> ",i4)') &
            dace% pe, processor, newcpu
       processor = newcpu
    end if

    if (.not. allocated (node_id)) return

    call p_gather (processor, cpu_, root=dace% pio)
    if (dace% lpio) then
       rank = [(i,i=0,dace% npe-1)]
       nn   = maxval (node_id)
       do i = 1, nn
          mask   = (node_id == i)
          cpumin = minval (cpu_, mask)
          cpumax = maxval (cpu_, mask)
          do j = cpumin, cpumax
             mm = count (cpu_ == j .and. mask)
             if (mm > 1) then
                write(0,'("Node",i3, "  CPU",i3,"  Ranks:",*(i4))') &
                     i, j, pack (rank, cpu_ == j .and. mask)
             end if
          end do
       end do
    end if
  end subroutine check_cpu_num
!==============================================================================
end module mo_cpu_time
