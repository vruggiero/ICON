!
!+ Provide system call entries for different platforms
!
MODULE mo_system
!
! Description:
!   This module provides system call entries for different platforms.
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
! V1_7         2009/08/24 Andreas Rhodin
!  No changes
! V1_8         2009/12/09 Andreas Rhodin
!  replace _ECMWF_ by _EXTNAME; call DWD trap-handler only on IBM
! V1_10        2010/04/20 Michael Gertz
!  Adaptions for SVN
! V1_13        2011/11/01 Andreas Rhodin
!  changed 3dvar revision numbers (CVS->SVN)
! V1_19        2012-04-16 Harald Anlauf
!  add section for Cray
! V1_27        2013-11-08 Harald Anlauf
!  Workaround for Cray (doesn't support getlog())
! V1_28        2014/02/26 Harald Anlauf
!  Intel/Cray XC30: modify workaround for broken getlogin
!
! Code Description:
! Language: Fortran 95.
! Software Standards:
!
! Authors:
! Andreas Rhodin  DWD  2001-2008
! Harald Anlauf   DWD  2007-2008
!-------------------------------------------------------------------

#if defined (NAGFOR) || defined (__NEC__)

  !----------------------------------------------------------------------
  ! The NAG compiler provides system call routines in a number of modules
  ! and similarly the NEC Aurora compiler.
  !----------------------------------------------------------------------
  use f90_unix,      only: abort,     &! unix abort call
                           flush       ! flush output buffer
  use f90_unix_proc, only: system      ! pass command line to system/shell
  use f90_unix_env,  only: getlogin,  &! get user name
                           gethostname ! get host name
  implicit none

  public :: maxrss                     ! Max. resident set size (if supported)

contains


!=========================================================================
#elif defined(__SX__)

  private
  public :: abort
  public :: getlogin
  public :: gethostname
  public :: flush
  public :: maxrss
  public :: install_trap

  !------------------------------------------------------------------------
  ! The SX provides an intrinsic routine 'abort' without parameters
  ! We here define a module procedure with optional argument 'text'.
  !------------------------------------------------------------------------
  interface abort
    module procedure abort_text
  end interface abort

  interface flush
    module procedure flush_mp
  end interface flush

contains

  subroutine abort_text (text)
  character(len=*),intent(in),optional :: text
  external :: abort
  if(present(text)) then
    write (0,*) 'abort:',text
  endif
  call abort
  end subroutine abort_text

  subroutine getlogin (user)
  character(len=*) ,intent(out) :: user
    character(len=8) :: us
!   call GETLOG (us)
!   if (us == "") call getenv ("LOGNAME", us)   ! Try LOGNAME with MPI
    call get_environment_variable ("LOGNAME", us)
    user = us
  end subroutine getlogin

  subroutine gethostname (host)
  character(len=*) ,intent(out) :: host
    integer           :: ihost, HOSTNM
    character(len=32) :: chost
    ihost = HOSTNM (chost)
    if (ihost /= 0) chost = 'unknown'
    host = chost
  end subroutine gethostname

  subroutine flush_mp (unit)
  integer, intent(in) :: unit
    external :: flush
    call flush (unit)            ! use library function
  end subroutine flush_mp

!=========================================================================
#elif defined(__ibm__)

  use mo_kind,       only: i4
  implicit none
  !------------------------------------------------------------------------
  ! The IBM RS6000 provides an intrinsic routine 'abort' without parameters
  ! We here define a module procedure with optional argument 'text'.
  !------------------------------------------------------------------------
  interface abort
    module procedure abort_text
  end interface abort

  interface flush
    module procedure flush_mp
  end interface flush

#if defined (_EXTNAME)
#define HOSTNM hostnm
#define GETLOG getlog
#define FLUSH  flush
#else
#define HOSTNM hostnm_
#define GETLOG getlog_
#define FLUSH  flush_
#endif

  public :: abort
  public :: getlogin
  public :: gethostname
  public :: install_trap
  public :: flush
  public :: maxrss

contains

  subroutine abort_text (text)
  character(len=*),intent(in),optional :: text
  intrinsic :: abort
  if(present(text)) then
    write (0,*) 'abort:',text
  endif
  call abort
  end subroutine abort_text

  subroutine gethostname (host)
  character(len=*) ,intent(out) :: host
    integer           :: ihost, HOSTNM
    character(len=32) :: chost
    ihost = HOSTNM (chost)
    if (ihost /= 0) chost = 'unknown'
    host = chost
  end subroutine gethostname

  subroutine getlogin (user)
  character(len=*) ,intent(out) :: user
    character(len=8) :: us
    call GETLOG (us)
    user = us
  end subroutine getlogin

  subroutine flush_mp (unit)
    integer, intent(in) :: unit
    integer  :: ios
    external :: FLUSH
    if (unit == 6) then                 ! Work around xlf 9.1 bug with stdout
       call FLUSH (unit)               ! IBM extension for older XLF compilers
    else
       flush (unit=unit, iostat=ios)    ! Use the Fortran 2003 FLUSH statement
    end if
  end subroutine flush_mp

#elif defined(__INTEL_COMPILER)

  !------------------------------------------------------------------------
  ! The Intel compiler provides system call routines in a number of modules
  !------------------------------------------------------------------------
  use ifport,  only: abort, &!
                     hostnam
#if (__INTEL_COMPILER < 910)
  use ifport,  only: flush
#endif
  use ifposix, only: pxfgetlogin,     &! Get login name (when tty present)
                     pxfstructcreate, &! Create instance of a structure
                     pxfstructfree,   &! Delete instance of a structure
                     pxfgeteuid,      &! Get effective user ID (Linux)
                     pxfgetuid,       &! Get real      user ID (Linux)
                     pxfgetpwuid,     &! Get information from /etc/passwd
                     pxfstrget         ! Get string value from structure
  implicit none

  private
  public :: abort
  public :: getlogin
  public :: gethostname
  public :: install_trap
  public :: flush
  public :: maxrss
#if defined (_CRAYFTN) || defined (__CRAYXC)
  public :: get_hugepages
#endif

contains

  subroutine getlogin (s, lens)
  character(len=*) ,optional ,intent(out) :: s
  integer          ,optional ,intent(out) :: lens
    integer            :: ierror
    integer            :: ilen
    character (len=64) :: c
    if (present(s))    s    = ''
    if (present(lens)) lens = 0
    ilen = len (c)                      ! Buffer size (needed by ifort)
    c = ""

#if defined (__CRAYXE) || defined (CRAYPAT)
    call getenv ("LOGNAME", c)
    if (present(s))    s    = trim     (c)
    if (present(lens)) lens = len_trim (c)
#else
    call pxfgetlogin (c, ilen, ierror)
    if (ierror==0) then
       !print *, "pxfgetlogin: user = '", trim (c), "'", ilen, ierror
       if (present(s))    s    = c
       if (present(lens)) lens = ilen
    else
       !------------------------------------------------------
       ! When the current process has no controlling tty,
       ! (e.g. when running in batch), (pxf)getlogin may fail.
       ! Use the POSIX getpwuid (geteuid ()) workaround below:
       !------------------------------------------------------
       call posix_getlogin (c, ilen)
       !print *, "posix_getlogin: user = '", trim (c), "'", ilen, len_trim (c)
       if (present(s))    s    = trim (c)
       if (present(lens)) lens = ilen
    endif
#endif
  end subroutine getlogin
  !----------------------------------------------------------------------------
  subroutine posix_getlogin (user, ilen)
    !-------------------------------------------------------------
    ! Reliable getlogin (the POSIX way): use getpwuid (geteuid ())
    !-------------------------------------------------------------
    character(len=*), intent(out) :: user
    integer,          intent(out) :: ilen

    integer, parameter :: N = 64                    ! Max. length of user name
    integer            :: ierror, uid, jhandle, k
    character(len=N)   :: c

    call pxfgeteuid (uid, ierror)
    if (ierror /= 0) call abort ("posix_getlogin: pxfgeteuid")

    call pxfstructcreate ("passwd", jhandle, ierror)
    if (ierror /= 0) call abort ("posix_getlogin: pxfstructcreate")

    call pxfgetpwuid (uid, jhandle, ierror)
    if (ierror /= 0) call abort ("posix_getlogin: pxfgetpwuid")

    ilen = len (c)
    c = ""
    call pxfstrget (jhandle, "pw_name", c, ilen, ierror)
    if (ierror /= 0) call abort ("posix_getlogin: pxfstrget")
    !-----------------------------------------
    ! Remove garbage characters (runtime bug?)
    ! This is a very ugly hack, but it works!
    !-----------------------------------------
    do k = 1, N
       if (llt (c(k:k), " ") .or. lgt (c(k:k), achar (126))) c(k:k) = " "
    end do
    k = index (c, " ")
    if (k > 0) then
       user = c(1:k-1)  ! Truncate user name after first blank
    else
       user = trim (c)
    end if
    ilen = len_trim (user)

    call pxfstructfree (jhandle, ierror)
    if (ierror /= 0) call abort ("posix_getlogin: pxfstructfree")
  end subroutine posix_getlogin
  !----------------------------------------------------------------------------
  subroutine gethostname (name)
  character(len=*) :: name
    integer :: result
    result = hostnam (name)
  end subroutine gethostname

#if (__INTEL_COMPILER >= 910)
  subroutine flush (unit)
    integer, intent(in) :: unit
    if (unit /= 0) then
       ! ifort 10 does not like flushing unit 0 ...
       flush (unit=unit)        ! Use the Fortran 2003 FLUSH statement
    end if
  end subroutine flush
#endif

!==============================================================================
#elif defined(__PGI) || defined(__FLANG)

  implicit none

  private
  public :: abort
  public :: gethostname
  public :: getlogin
  public :: install_trap
  public :: flush
  public :: maxrss
  public :: get_hugepages

  interface abort
    module procedure abort_text
  end interface abort

contains

  subroutine abort_text (text)
  character(len=*),intent(in),optional :: text
  if(present(text)) then
    write (0,*) 'abort:',text
  endif
  call exit(1)
  end subroutine abort_text

  subroutine gethostname (name)
  character(len=*) :: name
    integer :: result
    integer :: hostnm
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
  character(len=*) ,intent(out) :: user
    call getlog (user)
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    flush (unit=unit)
  end subroutine flush

!==============================================================================
#elif defined(__sun) || defined(__SUNPRO_F90) || defined(__SUNPRO_F95)

  implicit none

  private
  public :: abort
  public :: gethostname
  public :: getlogin
  public :: install_trap
  public :: flush
  public :: maxrss

  interface abort
    module procedure abort_text
  end interface abort

contains

  subroutine abort_text (text)
  character(len=*),intent(in),optional :: text
  if(present(text)) then
    write (0,*) 'abort:',text
  endif
  call exit(1)
  end subroutine abort_text

  subroutine gethostname (name)
  character(len=*) :: name
    integer :: result
    integer :: hostnm
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
  character(len=*) ,intent(out) :: user
!   call getlog (user)              ! getlog broken on hpc(GFS)!?
    call getenv ("LOGNAME", user)
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    flush (unit=unit)
  end subroutine flush

!==============================================================================
#elif defined (_CRAYFTN) || defined (__CRAYXC)
  ! Untested

  implicit none

  private
  public :: abort
  public :: gethostname
  public :: getlogin
  public :: install_trap
  public :: flush
  public :: maxrss
  public :: get_hugepages

  !-------------------------------------------------
  ! build a wrapper around an external abort routine
  !-------------------------------------------------
  interface abort
    module procedure abort_text
  end interface abort

contains

  ! Untested, maybe someone with access to a Cray needs to check this...
  subroutine abort_text (text)
    character(len=*), intent(in) :: text
    intrinsic :: abort
    write (0,*) 'abort:',text
    call abort ()
  end subroutine abort_text

  subroutine gethostname (name)
    character(len=*) :: name
    integer :: result
    integer :: hostnm           ! The intrinsic must be declared on Cray!?
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
    character(len=*) ,intent(out) :: user
    call getenv ("LOGNAME", user)
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    flush (unit=unit)           ! Use the Fortran 2003 FLUSH statement
  end subroutine flush

!==============================================================================
#elif defined(__G95__)

  implicit none

  private
  public :: abort
  public :: gethostname
  public :: getlogin
  public :: install_trap
  public :: flush
  public :: maxrss

  interface abort
     module procedure abort_text
  end interface abort

contains

  subroutine abort_text (text)
    character(len=*),intent(in) :: text
    intrinsic :: abort

    write (0,*) 'abort:', text
    call abort ()
  end subroutine abort_text

  subroutine gethostname (name)
  character(len=*) :: name
    integer :: result
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
  character(len=*) ,intent(out) :: user
    call getlog (user)
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    flush (unit=unit)           ! Use the Fortran 2003 FLUSH statement
  end subroutine flush

!==============================================================================
#elif defined(__GFORTRAN__) || (defined(__GNUC__) && !defined(__G95__))
  ! Only recent versions of gfortran (>= 4.3) identify themselves
  ! via a specific preprocessor symbol...

  implicit none

  private
  public :: abort
  public :: gethostname
  public :: getlogin
  public :: install_trap
  public :: flush
  public :: maxrss
#if defined (_CRAYFTN) || defined (__CRAYXC)
  public :: get_hugepages
#endif

  !-------------------------------------------------
  ! build a wrapper around an external abort routine
  !-------------------------------------------------
  interface abort
    module procedure abort_text
  end interface abort

contains

  subroutine abort_text (text)
    character(len=*), intent(in) :: text
    intrinsic :: abort
    write (0,*) 'abort:',text
    call abort ()
  end subroutine abort_text

  subroutine gethostname (name)
    character(len=*) :: name
    integer :: result
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
    character(len=*) ,intent(out) :: user

#if defined (__CRAYXE)
    call getenv ("LOGNAME", user)
#else
    call getlog (user)
#endif
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    flush (unit=unit)           ! Use the Fortran 2003 FLUSH statement
  end subroutine flush

!==============================================================================
#else
  ! Vanilla Fortran compiler

  !-------------------------------------------------
  ! build a wrapper around an external abort routine
  !-------------------------------------------------
  interface abort
    module procedure abort_text
  end interface abort

contains

  subroutine abort_text (text)
    character(len=*),intent(in),optional :: text
    intrinsic :: abort
!   external  :: abort
    if (present (text)) then
       write (0,*) 'abort:',text
    endif
    call abort ()
!   call exit(1)        ! If abort() not available, try exit()
  end subroutine abort_text

  subroutine gethostname (name)
    character(len=*) :: name
    integer :: result
    result = hostnm (name)
  end subroutine gethostname

  subroutine getlogin (user)
    character(len=*) ,intent(out) :: user

#if defined (__CRAYXE)
    call getenv ("LOGNAME", user)
#else
    call getlog (user)
#endif
  end subroutine getlogin

  subroutine flush (unit)
    integer, intent(in) :: unit
    ! Dummy routine
  end subroutine flush

#endif

  !============================================================================

  subroutine install_trap (coredump, status)
    !------------------------------
    ! install signal handler on IBM
    !------------------------------
    logical, intent(in),  optional :: coredump      ! Request core dump?
    integer, intent(out), optional :: status        ! Quiet return of status

#if defined(__ibm__)
    integer(i4) :: signal_trap, sigs(1), ires, core
    external    :: signal_trap
    logical     :: dump

    dump = .false.; if (present (coredump)) dump = coredump
    if (dump) then
       core  = 1 ! >0: coredump
    else
       core  = 0 ! =0: no coredump
    end if
    sigs     = 0 ! list of signals, last element=0, all zero = default
    ires     = 0
    ires     = signal_trap (core, sigs)
    if (present (status)) then
       status = 1
       if (ires == 1) status = 0
    else
       if (ires /= 1) &                         ! Installation failed
       write(0,*) 'signal_trap (core, sigs) =', ires
    end if
#else
    if (present (status)) status = 0
#endif
  end subroutine install_trap

  !============================================================================

  function maxrss ()
    !---------------------------------------------------------------------
    ! Returns the maximum resident set size in MB on the current processor
    ! (works only on systems with working getrusage, i.e. no Linux)
    !---------------------------------------------------------------------
    integer :: maxrss

#if defined (HAVE_MAXRSS)

#if ( defined (_AIX) || defined (__ibm__) || \
      defined (__SX__) || defined (__NEC__) || \
    ((defined (__cygwin__) || defined (__CYGWIN__)) && defined (__G95__)) )
    !--------------------------------------------------------------------
    ! We have a platform that supports the Fortran 2003 ISO C bindings
    !--------------------------------------------------------------------
    interface
       function imaxrss () bind(C,NAME="imaxrss_c")
         use, intrinsic  :: iso_c_binding, only : C_LONG
         integer(C_LONG) :: imaxrss
       end function imaxrss
    end interface
#elif ( defined (__linux__) && \
        ( defined (__GFORTRAN__) || \
          defined (__G95__) || \
          defined (__PGI) || \
          defined (__FLANG) || \
          defined (_CRAYFTN) || \
          defined (__NEC__) || \
         (defined (__INTEL_COMPILER) && (__INTEL_COMPILER >= 1010)) || \
          defined (__SUNPRO_F90) || defined (__SUNPRO_F95) ))
    !--------------------------------------------------------------------
    ! We have a platform that supports the Fortran 2003 ISO C bindings,
    ! but we need to call a function that only reports the current RSS.
    !--------------------------------------------------------------------
    interface
       function imaxrss () bind(C,NAME="irss_c")
         use, intrinsic  :: iso_c_binding, only : C_LONG
         integer(C_LONG) :: imaxrss
       end function imaxrss
    end interface
#else
    !--------------------------------------------------------------------
    ! Hopefully cfortran.h did the right thing for us...
    !--------------------------------------------------------------------
    interface
       function imaxrss ()
#ifdef CP4
       USE mo_kind, ONLY: i4
       INTEGER(i4) :: imaxrss
#else
       USE mo_kind, ONLY: i8
       INTEGER(i8) :: imaxrss
#endif
       end function imaxrss
    end interface
#endif
    !--------------------------------------------------------------------

    maxrss = (imaxrss () + 512) / 1024

#else
    !--------------------------------------------------------------------
    ! Default: not provided or supported on this host/platform.  You lose..
    !--------------------------------------------------------------------
    maxrss = 0
#endif

  end function maxrss

  !============================================================================

  subroutine get_hugepages (hugepages, smaps)
    integer, intent(out), optional :: hugepages   ! Hugepages (2M), this node
    integer, intent(out), optional :: smaps       ! Megabytes, this process
    !------------------------------------------------------
    ! Inquire Linux kernel for total used 2MB hugepages
    ! on current node and hugepage usage on current process
    !
    ! Based on suggestions by D. Sternkopf (Cray).
    !
    ! Per node:
    ! cat /sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages
    !
    ! Per process:
    ! grep -B 11 'KernelPageSize:     2048 kB' /proc/$pid/smaps| \
    ! grep "^Size:" | awk 'BEGIN{sum=0}{sum+=$2}END{print sum/1024}'
    !------------------------------------------------------
    character(*), parameter :: smaps_file     = "/proc/self/smaps"
    character(*), parameter :: hugepages_file =                 &
         "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages"

    character(*), parameter :: keyword1       = "Size:"
    character(*), parameter :: keyword2       = "KernelPageSize:"
    integer,      parameter :: lenk1           = len (keyword1)
    integer,      parameter :: lenk2           = len (keyword2)

    integer        :: ios, unit = 99
    integer        :: n, cur_size, kpg_size, sum_4k, sum_2m
    character(512) :: buffer

    if (present (hugepages)) then
       n = -1
       open (unit,file=hugepages_file,status='old',action='read',iostat=ios)
       if (ios == 0) then
          read(unit,*,iostat=ios) n
          if (ios /= 0) then
             print *, "read failed on " // hugepages_file
          end if
       end if
       if (ios <= 0) close (unit)
       hugepages = n
    end if

    if (present (smaps)) then
       sum_4k = 0
       sum_2m = 0
       open (unit,file=smaps_file,status='old',action='read',iostat=ios)
       if (ios == 0) then
          cur_size = -1
          kpg_size = -1
          do
             buffer = ""
             read (unit,'(a)',iostat=ios,end=999) buffer
             if (ios /= 0) exit
             if      (buffer(:lenk1) == keyword1) then
                read (buffer(lenk1+1:),*) cur_size
             else if (buffer(:lenk2) == keyword2) then
                read (buffer(lenk2+1:),*) kpg_size
                select case (kpg_size)
                case (4)
                   sum_4k = sum_4k + cur_size   ! 4k pages
                case (2048)
                   sum_2m = sum_2m + cur_size   ! 2M pages
                case default
                   print *, "Unsupported KernelPageSize:", kpg_size
                end select
                !print *, "sum_4k, sum_2m=", sum_4k, sum_2m
             end if
          end do
999       continue
       end if
       if (ios <= 0) close (unit)
       smaps = sum_2m / 1024            ! Convert to MB
    end if
  end subroutine get_hugepages

end module mo_system

!--------------------------------------------------
! external routine as a substitute for the routine
! called by the COSMO observation operator library.
! (originally linked from dwdlib/grib1)
! provided by support/fsleep.c .
!--------------------------------------------------
! subroutine fsleep (isec)
! integer, intent(in) :: isec
!   ..
! end subroutine fsleep
